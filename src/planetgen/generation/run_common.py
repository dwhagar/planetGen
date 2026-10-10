# planetgen/generation/run_common.py

"""
Run Common
==========

What every generation command shares: the console progress bar, the
`--strict` refuse-or-warn rule, the run's sector/system/phenomenon
counts, the work queue a run hands its sectors to, and the generation
statistics and estimate a bulk run checks before it starts.
"""

import json
import os
import sys
from collections import Counter

from rich.progress import (
    BarColumn, MofNCompleteColumn, Progress, ProgressColumn, TextColumn, TimeElapsedColumn,
)
from rich.text import Text

from planetgen.queue import progress_file, progress_rate, work as workQueue
from planetgen.db import store
from planetgen.generation import stats as generationStats
from planetgen import tuning as program_constants
from planetgen.util import log
from planetgen.galaxy.skeleton import expected_system_count_at_density_1



_ACTIVE_PROGRESS = []
"""list: The `_ReportingProgress` displays that are started, innermost last (`active_progress`)."""


def active_progress():
    """The open generation display (`_generation_progress` inside its `with`), or `None`; the steps
    (`planetgen.generation.steps`) draw on it."""
    return _ACTIVE_PROGRESS[-1] if _ACTIVE_PROGRESS else None


class _ReportingProgress(Progress):
    """
    `rich.progress.Progress` that keeps a `progress_rate.DecayingRate` per
    task (the rate behind `_DecayingRemainingColumn`'s ETA, PERF.7) and
    mirrors its most recently changed task, with that rate and ETA, to
    `planetgen.queue.progress_file` (a no-op unless the web interface
    started this run).
    """

    main_task = None
    """The task the progress file reports when `detail_task` changes."""

    main_step = None
    """The `steps.Step` drawing `main_task`, if one does (a step shown while it is set is the second bar)."""

    def start(self):
        super().start()
        _ACTIVE_PROGRESS.append(self)

    def stop(self):
        if self in _ACTIVE_PROGRESS:
            _ACTIVE_PROGRESS.remove(self)
        for task_id in list(self._tasks):
            self._report(task_id, force=True)
        super().stop()

    detail_task = None
    """A second bar under `main_task` (PERF.4's slow-layer bar), written
    to the progress file as its `detail`."""

    def add_task(self, description, *args, prior=None, tau=None, **kwargs):
        """Adds a bar. `prior` is the rate recorded for this kind of work in
        the bar's own units a second (PERF.33), which the live rate is
        blended into; `tau` the decay time constant (`progress_rate.time_constant`)."""
        rate = progress_rate.DecayingRate(prior=prior) if tau is None else progress_rate.DecayingRate(tau=tau,
                                                                                                    prior=prior)
        task_id = super().add_task(description, *args, rate=rate, **kwargs)
        self._report(task_id, force=True)
        return task_id

    def remove_task(self, task_id):
        super().remove_task(task_id)
        if self.main_task is not None:
            self._report(self.main_task, force=True)

    def _detail(self):
        task = self._tasks.get(self.detail_task) if self.detail_task is not None else None
        if task is None:
            return None
        rate = task.fields.get("rate")
        remaining = None if task.total is None else task.total - task.completed
        return {"description": task.description.strip(), "completed": task.completed, "total": task.total,
                "eta_s": rate.eta(remaining) if rate is not None else None}

    def _rate(self, task_id):
        task = self._tasks.get(task_id)
        return task.fields.get("rate") if task is not None else None

    def _report(self, task_id, force=False):
        if self.main_task is not None:
            # A run with a main bar always writes that one, with any
            # second bar as its `detail`.
            task_id = self.main_task
        task = self._tasks.get(task_id)
        if task is None:
            return
        rate = task.fields.get("rate")
        remaining = None if task.total is None else task.total - task.completed
        progress_file.report(task.completed, task.total, task.description, force=force,
                            rate=rate.rate if rate is not None else None,
                            eta_s=rate.eta(remaining) if rate is not None else None, detail=self._detail(),
                            percent=bool(task.fields.get("percent")))

    def update(self, task_id, **kwargs):
        before = self._completed(task_id)
        super().update(task_id, **kwargs)
        self._record(task_id, before)
        self._report(task_id, force=kwargs.get("total") is not None)

    def advance(self, task_id, advance=1):
        before = self._completed(task_id)
        super().advance(task_id, advance)
        self._record(task_id, before)
        self._report(task_id)

    def _completed(self, task_id):
        task = self._tasks.get(task_id)
        return task.completed if task is not None else 0

    def _record(self, task_id, before):
        rate = self._rate(task_id)
        if rate is not None:
            rate.add(self._completed(task_id) - before)



class _DecayingRemainingColumn(ProgressColumn):
    """
    Time left for a task at its decaying-average rate (PERF.7,
    `progress_rate.DecayingRate`), rather than rich's own estimate from
    its last few updates, which jumps about when many workers finish
    together. Blank until the first unit is done; 0:00:00 once finished.
    """

    max_refresh = 0.5

    def render(self, task):
        rate = task.fields.get("rate")
        remaining = None if task.total is None else task.total - task.completed
        seconds = rate.eta(remaining) if rate is not None else None
        if task.finished:
            seconds = 0
        if seconds is None:
            return Text("-:--:--", style="progress.remaining")
        hours, rest = divmod(int(round(seconds)), 3600)
        minutes, secs = divmod(rest, 60)
        return Text(f"{hours}:{minutes:02d}:{secs:02d}", style="progress.remaining")


class _CountColumn(MofNCompleteColumn):
    """`MofNCompleteColumn`, or the share done for a task added with
    `percent=True` (PERF.9's bright-star bar, whose units are weighted
    work rather than anything worth counting)."""

    def render(self, task):
        if task.fields.get("percent"):
            share = min(task.completed / task.total, 1.0) if task.total else 0.0
            return Text(f"{100 * share:.0f}%", style="progress.download")
        return super().render(task)


def _generation_progress(disable=False):
    """
    Builds the shared `rich.progress.Progress` used by `run_galaxy`'s three
    modes -- one "Sectors" task per run tracking how many sectors have been
    generated so far. Deliberately not used by `run_sector`'s own
    `--num-sectors` loop: that command's own run is normally short enough
    (and its sector count small enough) that a bar added more visual noise
    than it was worth, whereas a `galaxy` run (a whole ring, a
    neighborhood, or a random start's default 12 pc one) can mean
    thousands of sectors and legitimately benefit from a progress display.
    There also used to be a second, nested bar for "systems in the current
    sector," added/removed once per sector; that one is gone for good --
    most sectors hold anywhere from zero to a handful of systems, and
    system generation itself is fast, so a bar that flashed on and off
    again within a single frame for nearly every sector was pure noise,
    not something worth reviving alongside this one.

    Every task shows both elapsed time and an estimated time remaining,
    for as long as it runs (`TimeElapsedColumn`/`_DecayingRemainingColumn`)
    -- the estimate comes from a decaying average of units finished per
    second (`progress_rate.DecayingRate`, a 60 s time constant), counted
    in the bar's own unit (sectors, or layers of bright stars), so it
    stays steady while several workers report at once.

    Callers use this as a context manager (`with _generation_progress() as
    progress:`); `rich.progress.Progress` is a `Live` display under the
    hood, so **every status line logged while it's active must go through
    `progress.console`, never a raw stdout write** -- printing directly to
    stdout fights with the `Live` region's own redraws (each plain `print`
    call forces the bar to erase itself, scroll up with the new text, and
    get redrawn at the bottom again), which is exactly what caused the
    flicker/scrolling a real terminal used to show once before. `run_galaxy`
    calls `log.set_console(progress.console)` right after opening this
    context manager (and `log.reset_console()` once it's closed), so every
    `log.normal(...)`/`log.debug(...)` call made anywhere during a galaxy
    run -- including deep inside `planetgen` modules -- is
    automatically routed through `progress.console.print(...)` instead of a
    raw stdout write for as long as the bar is live. Routed that way, rich
    prints each line safely *above* the live region and leaves the bar
    itself pinned at the bottom, redrawn in place with no flicker. This only
    matters when stdout is a real interactive terminal in the first place --
    `Console` auto-detects that (`Console.is_terminal`) and falls back to
    plain, periodic line-by-line bar output otherwise (piped to a file, a CI
    log, etc.), so no separate handling is needed for that case.

    The same counts also go to `$PLANETGEN_PROGRESS_FILE` when it is set
    (`planetgen.queue.progress_file`), so the web interface's Generate page
    can show a run it started in the background.

    Args:
        disable (bool): No bars at all (`--estimate-only`).

    Returns:
        Progress: Not yet started.
    """
    return _ReportingProgress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        _CountColumn(),
        TextColumn("[dim]elapsed"),
        TimeElapsedColumn(),
        TextColumn("[dim]remaining"),
        _DecayingRemainingColumn(),
        disable=disable,
    )


def _strict(args):
    """Whether this run stops on what the console otherwise only warns
    about (`--strict`, GEN.81). The console never tells the user no: it
    warns, then does what was asked."""
    return bool(getattr(args, "strict", False))


def _refuse_or_warn(args, message, code=1):
    """GEN.81: under `--strict`, logs `message` as an error and exits with
    `code`; otherwise logs it as a warning and returns, so the run goes
    ahead."""
    if _strict(args):
        log.error(message)
        raise SystemExit(code)
    log.normal(f"WARNING: {message} Going ahead (--strict stops here instead).")


RUN_COUNTS = Counter()
"""Counter: What this run saved (`sectors`, `systems`, `phenomena`), for
the activity log's `generate.finish` line."""


def _count_sector(sector):
    RUN_COUNTS["sectors"] += 1
    RUN_COUNTS["systems"] += len(sector.entries)
    RUN_COUNTS["phenomena"] += len(sector.phenomena)


def _edge_pc():
    """
    The sector edge length used for every grid computation, in parsecs --
    the one standard, `program_constants.DEFAULT_SECTOR_EDGE_PC` (not a
    CLI option for `galaxy` or `plan`, so the grid and the skeleton always
    agree).
    """
    return float(program_constants.DEFAULT_SECTOR_EDGE_PC)


def _log_level(args):
    """The console severity `main` configured from `--quiet`/`--debug`,
    for the work queue's workers."""
    if getattr(args, "quiet", False):
        return log.SILENT
    if getattr(args, "debug", None) is not None:
        return log.DEBUG
    return log.NORMAL


def _work_queue(args, title):
    """
    The `workQueue.WorkQueue` a run hands its sectors to (PERF.8):
    `--workers` (or `PLANETGEN_WORKERS`) worker processes, by default 80%
    of the cores less one when MySQL runs on this machine, with the
    control database's lease so only one run's workers use the machine
    at a time. One worker generates every sector right here, in order
    (still recorded in the job tree, ADM.12, without the lease).
    """
    mysql_config = store.mysql_config_from_args(args)
    workers = _worker_count(args)
    log.debug(f"{title}: {workers} worker process(es) ({workQueue.cpu_count()} cores).")
    return workQueue.WorkQueue(
        title, workers=workers,
        control_config=store.control_mysql_config(mysql_config),
        log_level=_log_level(args), debug_file=getattr(args, "debug", None) or None,
    )


def _worker_count(args):
    """The worker processes a run uses (`workQueue.worker_count`)."""
    return workQueue.worker_count(getattr(args, "workers", None), store.mysql_config_from_args(args).host)


_RUN_STATS = {}
"""dict: The run's `generationStats.GenerationStats`, per control
database, loaded on first use (`_generation_stats`)."""


ESTIMATE_PREFIX = "ESTIMATE "
"""str: Starts the one JSON line `--estimate-only` prints (the Generate
page reads it)."""


class _EstimateOnly(Exception):
    """Raised by `_check_estimate` under `--estimate-only`: the run stops
    there, before anything is written, and `run_galaxy`/`run_sector`
    print the estimate."""

    def __init__(self, what, result):
        super().__init__(what)
        self.what = what
        self.result = result


def _generation_stats(args):
    """The stored speeds and sizes (`generationStats.GenerationStats`), loaded once per run; see
    `_stats_for_config`."""
    return _stats_for_config(store.mysql_config_from_args(args))


def _stats_for_config(mysql_config):
    """The stored speeds and sizes of the server `mysql_config` names, loaded once per run; none (defaults
    only, nothing recorded) when `PLANETGEN_GENERATION_STATS` is `0` (the test suite's setting)."""
    if os.environ.get(generationStats.STATS_ENV_VAR, "1").strip().lower() in ("0", "off", "no", "false"):
        return _RUN_STATS.setdefault(None, generationStats.GenerationStats())
    control = store.control_mysql_config(mysql_config)
    key = (control.host, control.port, control.database)
    if key not in _RUN_STATS:
        _RUN_STATS[key] = generationStats.GenerationStats(control)
    return _RUN_STATS[key]


def _e_value():
    """Systems per standard sector at density 1."""
    return expected_system_count_at_density_1(program_constants.DEFAULT_SECTOR_EDGE_LY)


def _sector_density(sector_args):
    """A sector's density (`relative_density`), from `--density` or
    `--num-systems`."""
    if sector_args.density is not None:
        return float(sector_args.density)
    return float(sector_args.num_systems or 0) / _e_value()


def _expected_systems(sector_args):
    if sector_args.num_systems is not None:
        return float(sector_args.num_systems)
    return float(sector_args.density or 0.0) * _e_value()


def _record_sector(args, result, seconds):
    """Adds one filled sector to its density's speed bucket (PERF.53: not one that made no systems, which
    would pull the per-sector and per-system averages toward nothing)."""
    if not result.get("systems"):
        return
    _generation_stats(args).record("sector", result.get("density"), seconds,
                                   systems=result.get("systems", 0), stars=result.get("stars", 0),
                                   workers=_worker_count(args))


def _finish_stats(args):
    """After a run: measures the galaxy database's size per system and
    writes the run's speeds back."""
    for stats in _RUN_STATS.values():
        if not stats.available:
            continue
        config = store.mysql_config_from_args(args)
        try:
            conn = store.get_connection(config)
        except Exception as exc:  # noqa: BLE001 -- stats never fail a run
            log.debug(f"Generation stats: can't measure the database ({exc}).")
        else:
            try:
                stats.measure_size(conn, config.database)
            finally:
                conn.close()
        stats.flush()
    # The next run (or web request) reads them fresh.
    _RUN_STATS.clear()


def _interactive():
    return sys.stdin.isatty() and sys.stdout.isatty()


def _ask(progress, question):
    """Asks `question` on the terminal with the progress bars paused."""
    if progress is not None:
        progress.stop()
    try:
        return input(question).strip().lower()
    except EOFError:
        return ""
    finally:
        if progress is not None:
            progress.start()


def _estimate_sectors(args, sector_args_list):
    """The `generationStats.Estimate` for filling these sectors, with the
    database disk checked."""
    config = store.mysql_config_from_args(args)
    stats = _generation_stats(args)
    result = generationStats.estimate(
        [(_sector_density(sector_args), _expected_systems(sector_args)) for sector_args in sector_args_list],
        stats, config.database, workers=_worker_count(args),
    )
    try:
        conn = store.get_connection(config)
        try:
            disk = generationStats.database_disk(conn, config.host, config.database)
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- no disk reading: nothing refused for space
        log.debug(f"Generation estimate: no disk reading ({exc}).")
        disk = None
    generationStats.check_disk(result, disk)
    return result


def _check_estimate(args, sector_args_list, what, progress=None):
    """
    PERF.3: before a bulk run writes anything, estimates its size (+10%)
    and time and shows them; refuses it when it would take more than a
    quarter of the database disk or leave less than 5 GB free; and, on a
    terminal, asks before filling more than one sector (`--yes` skips the
    question, never the refusal). Off a terminal (the Generate page's
    jobs, scripts) it shows the estimate and goes on: the page asked
    first. Once per run (`args._estimate_checked`).

    Raises:
        _EstimateOnly: Under `--estimate-only`.
        SystemExit: Refused (1), or the answer was no (0).
    """
    if getattr(args, "_estimate_checked", False):
        return
    result = _estimate_sectors(args, sector_args_list)
    if getattr(args, "estimate_only", False):
        raise _EstimateOnly(what, result)
    args._estimate_checked = True
    args._sector_prior = _sector_prior(result)
    args._sector_estimate = result
    log.normal(f"Estimate for {what}: {result.summary()}")
    if result.refusal:
        _refuse_or_warn(args, result.refusal)
    if result.sectors > 1 and not getattr(args, "yes", False) and _interactive():
        if _ask(progress, f"Generate these {result.sectors:,} sectors? [y/N] ") not in ("y", "yes"):
            log.normal("Nothing was generated.")
            raise SystemExit(0)


def _sector_prior(result):
    """`(sectors a second, decay time constant)` the stored speeds predict for the run `result`
    estimated, or `None` when nothing was recorded for it yet (PERF.33)."""
    if not result.measured or not result.sectors or not result.seconds or result.seconds <= 0:
        return None
    mean_task_seconds = result.seconds * result.workers / result.sectors
    return (result.sectors / result.seconds, progress_rate.time_constant(mean_task_seconds, result.workers))


def sector_bar(progress, args, description, total):
    """The step that counts the sectors of a run (PERF.51), predicted and started from the rate this
    server recorded for such sectors (`_check_estimate` worked it out): drawn at once when the run is
    expected to take over 15 seconds or has no history, `record=False` because each sector is recorded
    as it is saved. Use as a context manager."""
    from planetgen.generation import steps

    prior = getattr(args, "_sector_prior", None)
    estimate = getattr(args, "_sector_estimate", None)
    predicted = estimate.seconds if estimate is not None and estimate.measured else None
    return steps.Step(description, "sector", total, progress=progress, stats=_stats_or_none(args),
                      predicted=predicted, prior=prior[0] if prior else None, tau=prior[1] if prior else None,
                      record=False)


def _stats_or_none(args):
    try:
        return _generation_stats(args)
    except Exception:  # noqa: BLE001 -- args without a database (a test's stand-in)
        return None


def _print_estimate(exc):
    """`--estimate-only`'s output: the summary, then one JSON line."""
    log.normal(f"Estimate for {exc.what}: {exc.result.summary()}")
    if exc.result.refusal:
        log.normal(exc.result.refusal)
    body = exc.result.as_dict()
    body["what"] = exc.what
    sys.stdout.write(ESTIMATE_PREFIX + json.dumps(body) + "\n")
    sys.stdout.flush()
