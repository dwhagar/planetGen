# planetgen/generation/steps.py

"""
Steps (PERF.51): the one way a generation or maintenance run reports the
progress of a piece of work.

A step has a name, a stats kind (what `generation_stats` records its speed
under) and a work count. At its start the helper predicts how long it will
take from that kind and count and the speeds this server recorded
(`GenerationStats.pool_rate`, PERF.32/PERF.33): a step predicted to take
longer than `tuning.PROGRESS_BAR_SECONDS` (15 s) draws its bar at once, one
predicted shorter draws nothing, and one that runs past 15 s anyway gets its
bar the moment it passes it. A step with nothing recorded for its kind has no
prediction and is treated as long, so a first run still shows its bar. On
leaving, a finished step records its measured speed back under its kind, so
the next prediction has history.

The first step shown on a display is its main bar, any step shown while
another is is the second bar under it (`detail`). A step draws on the
display it is given (`progress=`; `None` draws nothing), else the open one;
a maintenance command that has no display asks for one of its own for as long
as it draws (`own=True`).

A step can also run inside a worker: `worker_step` there sends its progress
down a channel (`planetgen.queue.redisqueue.Channel`, or any object with
`put`), and a `Relay` in the parent turns those reports into steps drawn
in the parent's display (PERF.50 is the first application).

```
with steps.step("Neighbours", "link", total, args=args) as bar:
    ...
    bar.update(advance=1)
```
"""

import queue as queue_module
import contextlib
import sys
import threading
import time

from planetgen import tuning
from planetgen.util import log

UNSET = object()
"""A parameter left to the helper to work out."""

REPORT_INTERVAL_SECONDS = 0.25
"""float: How often a worker's step reports its progress to the parent."""


def predict(stats, kind, total, workers=1):
    """
    Seconds a step of `kind` with `total` units is expected to take, from the
    pool rate recorded for it, or `None` when it can't be told (no history
    for the kind, or no known count).
    """
    if stats is None or kind is None or not total:
        return None
    rate = stats.pool_rate(kind, workers)
    return total / rate if rate else None


def _stats_for(args):
    if args is None:
        return None
    from planetgen.generation import run_common

    return run_common._generation_stats(args)


def _active_progress():
    from planetgen.generation import run_common

    return run_common.active_progress()


STEP_KINDS = {
    "sector": "Sector fill",
    "scatter": "Bright-star layer",
    "phenomena": "Phenomena layer",
    "phenomena-clear": "Phenomena scatter: clearing",
    "phenomena-special": "Phenomena scatter: special rows",
    "phenomena-insert": "Phenomena scatter: writing the special rows",
    "phenomena-stamp": "Phenomena scatter: epoch stamp",
    "link": "Neighbour linking",
    "paths": "Sector paths",
    "backfill": "Bright-star backfill",
    "topup": "Backfill top-up",
    "save": "Sector save",
    "population": "Population pass",
    "reset": "Reset",
    "orbits": "Orbit update",
    "migrate": "Migration",
    "stages": "Stages",
    "containment": "Containment refresh",
    "nearest": "Nearest-systems refresh",
    "warm-map": "Galaxy Map warm-up",
    "dedupe-sectors": "Sector name pass",
    "dedupe-systems": "System name pass",
}
"""dict: Every kind of step the program times (UX.84): the `generation_stats`
kind and what the Stats page calls it. A step built with a kind that is not
here fails `tests/test_step_registry.py`, so a new long step is registered
(and gets a bar) when it is written. A kind may carry a `:label` suffix for
a detail bar under its step."""


class Step:
    """
    One piece of work with its own progress bar, drawn only when it is
    long (see this module's docstring). Use as a context manager.

    Args:
        name (str): The bar's description.
        kind (str, optional): The `generation_stats` kind its speed is
            predicted from and recorded under; `None` neither predicts nor
            records.
        total (float, optional): Units of work expected, or `None` while
            unknown (`update(total=...)` later).
        args (argparse.Namespace, optional): The run's arguments, for its
            stored speeds (`stats` overrides).
        stats (GenerationStats, optional): The stored speeds.
        workers (int): Worker processes the step runs on (its rate is the
            pool's at that count).
        progress (Progress, optional): The display to draw on; `None` draws
            nothing; left out, the open display (`run_common.active_progress`), if any.
        own (bool): With no display to draw on, open one of its own while the step draws.
        predicted (float or None): Seconds expected, given by a caller that
            knows better than `kind` and `total` (`None`: unknown).
        prior (float or None): The recorded rate the bar's time left starts
            from (default: the pool rate for `kind`).
        tau (float, optional): The bar's decay time constant.
        percent (bool): A share-done bar (weighted work).
        record (bool): Record the speed on leaving (False when the caller
            records its own tasks).
        force (bool): Draw at once whatever the prediction.
        threshold (float): Seconds a step may take before it draws.
    """

    def __init__(self, name, kind=None, total=None, *, args=None, stats=None, workers=1, progress=UNSET,
                 predicted=UNSET, prior=UNSET, tau=None, percent=False, record=True, force=False,
                 threshold=None, own=False):
        self.name = name
        self.kind = kind
        self.total = total
        self.stats = stats if stats is not None else _stats_for(args)
        self.workers = max(1, int(workers or 1))
        self._display = progress
        self.own = own
        self.predicted = predict(self.stats, kind, total, self.workers) if predicted is UNSET else predicted
        if prior is UNSET:
            prior = self.stats.pool_rate(kind, self.workers) if self.stats is not None and kind else None
        self.prior, self.tau = prior, tau
        self.percent, self.record, self.force = percent, record, force
        self.threshold = tuning.PROGRESS_BAR_SECONDS if threshold is None else threshold
        self.done = 0.0
        self.task = None
        self.started = None
        self._lock = threading.RLock()
        self._timer = None
        self._own = None
        self._closed = False
        self._is_main = False

    # -- life ------------------------------------------------------------

    @property
    def shown(self):
        """Whether the step draws a bar now."""
        return self.task is not None

    def __enter__(self):
        self.started = time.monotonic()
        if self.force or self.predicted is None or self.predicted > self.threshold:
            self._show()
        else:
            # Predicted short: nothing is drawn, unless it runs past the limit anyway.
            self._timer = threading.Timer(self.threshold, self._show)
            self._timer.daemon = True
            self._timer.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close(success=exc_type is None)
        return False

    def close(self, success=True):
        """Ends the step: removes (or completes) its bar, and on success records its speed."""
        with self._lock:
            if self._closed:
                return
            self._closed = True
            if self._timer is not None:
                self._timer.cancel()
            seconds = time.monotonic() - self.started if self.started is not None else 0.0
            self._hide(success)
        if success and self.record and self.kind and self.stats is not None:
            self._record(seconds)

    def _record(self, seconds):
        units = self.total if self.total else self.done
        if not units or seconds <= 0:
            return
        try:
            self.stats.record(self.kind, 0.0, seconds, systems=units, stars=units, workers=self.workers)
            self.stats.flush()      # a step is rare and long: the next run's prediction should not wait for the 30 s flush
        except Exception as exc:  # noqa: BLE001 -- statistics never fail a run
            log.debug(f"Steps: could not record {self.kind}: {exc}")

    # -- the bar ---------------------------------------------------------

    def _show(self):
        with self._lock:
            if self._closed or self.task is not None:
                return
            progress = _active_progress() if self._display is UNSET else self._display
            if progress is None and self.own:
                from planetgen.generation import run_common

                progress = self._own = run_common._generation_progress(disable=not sys.stdout.isatty())
                progress.start()
                log.set_console(progress.console)
            if progress is None:
                return          # no display: nothing is drawn, the work is still tracked
            self._display = progress
            fields = {"percent": True} if self.percent else {}
            if self.prior:
                fields["prior"] = self.prior
            if self.tau:
                fields["tau"] = self.tau
            total = self.total
            self.task = progress.add_task(self.name, total=total, completed=min(self.done, total) if total else 0,
                                          **fields)
            main = getattr(progress, "main_step", None)
            if main is None or main is self:
                progress.main_step = self
                progress.main_task = self.task
                self._is_main = True
            else:
                progress.detail_task = self.task

    def _hide(self, success):
        progress = self._display
        if self.task is None or progress is None or progress is UNSET:
            self._stop_own()
            return
        if success and self.total:
            progress.update(self.task, completed=self.total)
        if self._is_main:
            progress.main_step = None
            if self.total is None:
                progress.remove_task(self.task)     # never measured: no bar left unmeasured
        else:
            progress.detail_task = None
            progress.remove_task(self.task)
        self._stop_own(keep_bar=self._is_main)
        self.task = None

    def _stop_own(self, keep_bar=False):
        if self._own is not None:
            self._own.stop()
            log.reset_console()
            self._own = None

    # -- progress --------------------------------------------------------

    def update(self, advance=None, completed=None, total=None, description=None):
        """Moves the bar like `Progress.update` (kept when no bar is
        drawn, so a bar that appears late starts where the work is)."""
        with self._lock:
            if total is not None:
                self.total = total
            if completed is not None:
                self.done = float(completed)
            if advance:
                self.done += advance
            if description is not None:
                self.name = description
            if self.task is None or self._display is None or self._display is UNSET:
                return
            kwargs = {}
            if total is not None:
                kwargs["total"] = total
            if completed is not None or advance:
                kwargs["completed"] = min(self.done, self.total) if self.total else self.done
            if description is not None:
                kwargs["description"] = description
            if kwargs:
                self._display.update(self.task, **kwargs)

    def advance(self, amount=1):
        self.update(advance=amount)


@contextlib.contextmanager
def hooked(on_progress, name, kind, total=None, **kwargs):
    """
    For a function that reports `on_progress(label, done, total)` to its
    caller: yields `on_progress` itself when the caller gave one (the caller
    draws it), and otherwise a callback of the same shape that moves a
    `Step(name, kind, total)` on the open display, so the work has a bar
    whoever calls it (UX.84). Nothing is drawn when there is no display.
    """
    if on_progress is not None:
        yield on_progress
        return
    if kwargs.get("progress", UNSET) is UNSET and not kwargs.get("own") and _active_progress() is None:
        yield lambda _label, _done, _size: None        # nobody is drawing: not even a timer
        return
    with Step(name, kind, total, **kwargs) as bar:
        def report(_label, done, size):
            bar.update(completed=done, total=size)

        yield report


def step(name, kind=None, total=None, **kwargs):
    """`Step(name, kind, total, ...)`."""
    return Step(name, kind, total, **kwargs)


# ---------------------------------------------------------------------------
# Steps inside a worker
# ---------------------------------------------------------------------------

class _NullStep:
    """What `worker_step` gives a task run without a channel."""

    def __enter__(self):
        return self

    def __exit__(self, *_exc):
        return False

    def update(self, **_kwargs):
        pass

    def advance(self, _amount=1):
        pass


class _WorkerStep(_NullStep):
    """A step in a worker: reports to the parent down `channel`, which draws it."""

    def __init__(self, channel, key, name, kind, total, workers):
        self.channel, self.key = channel, key
        self.start = (name, kind, total, workers)
        self.done = 0.0
        self._last = 0.0

    def _send(self, item):
        try:
            self.channel.put(item)
        except Exception as exc:  # noqa: BLE001 -- progress never fails the work
            log.debug(f"Steps: could not report {self.key}: {exc}")

    def __enter__(self):
        self._send(("start", self.key) + self.start)
        return self

    def __exit__(self, exc_type, _exc, _tb):
        self._send(("end", self.key, exc_type is None))
        return False

    def update(self, advance=None, completed=None, total=None, description=None):
        if completed is not None:
            self.done = float(completed)
        if advance:
            self.done += advance
        now = time.monotonic()
        if total is not None or description is not None or now - self._last >= REPORT_INTERVAL_SECONDS:
            self._last = now
            self._send(("progress", self.key, self.done, total, description))

    def advance(self, amount=1):
        self.update(advance=amount)


def worker_step(channel, key, name, kind=None, total=None, workers=1):
    """
    A step for work inside a worker process: its progress goes down
    `channel` (carried in the task's payload) to a `Relay` in the parent.
    With no channel (`None`) it does nothing.
    """
    if channel is None:
        return _NullStep()
    return _WorkerStep(channel, key, name, kind, total, workers)


class Relay:
    """
    The parent's end of `worker_step`: turns each worker's reports into a
    `Step` drawn on the run's display (the same 15 second rule), and ends
    any a worker left open. Give workers `Relay.put` directly with one
    worker, or a channel `drain` reads from.

    Args:
        args (argparse.Namespace, optional): The run's arguments (its stored speeds).
        progress (Progress, optional): The display to draw on.
    """

    def __init__(self, args=None, progress=UNSET):
        self.args, self.progress = args, progress
        self.steps = {}
        self._lock = threading.Lock()

    def put(self, item):
        kind = item[0]
        with self._lock:
            if kind == "start":
                _tag, key, name, stats_kind, total, workers = item
                if key not in self.steps:
                    self.steps[key] = Step(name, stats_kind, total, args=self.args, workers=workers,
                                           progress=self.progress).__enter__()
            elif kind == "progress":
                _tag, key, done, total, description = item
                step_ = self.steps.get(key)
                if step_ is not None:
                    step_.update(completed=done, total=total, description=description)
            elif kind == "end":
                _tag, key, success = item
                step_ = self.steps.pop(key, None)
                if step_ is not None:
                    step_.close(success=bool(success))

    def drain(self, channel, stop):
        """Reads `channel` into `put` until `stop` is set (run in a thread)."""
        while True:
            try:
                item = channel.get(timeout=0.25)
            except queue_module.Empty:
                if stop.is_set():
                    return
                continue
            except (EOFError, OSError):
                return
            self.put(item)

    def close(self):
        """Ends every step still open (a worker that died mid-step)."""
        with self._lock:
            for step_ in list(self.steps.values()):
                step_.close(success=False)
            self.steps.clear()
