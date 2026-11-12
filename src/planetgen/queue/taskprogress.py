# planetgen/queue/taskprogress.py

"""
A task's own progress bar and ETA (PERF.147): one pure function from a task's counter (`done` of `total`, `started`), a clock and the rate the server has recorded for the same
kind of work, to whether the task shows a bar, how far along it is and how long it has left.

* **When there is a bar.** A task that is expected to take longer than `SHOW_SECONDS` (the recorded rate says `total / rate` over it) shows its bar from the start, with the recorded
  rate's ETA. A task with no recorded rate, or one expected to be quick, gains its bar once it has run longer than `SHOW_SECONDS`. A task that ends in less shows no bar.
* **The ETA.** The rate to estimate with moves from the recorded one to the running task's own as it proceeds: `rate = w * live + (1 - w) * stored`, `w = done / total / LIVE_FULL_AT`
  (at most 1), so the running task counts for more with every item and fully from `LIVE_FULL_AT` of the way. The live side is held back until the task has run `HOLD_SECONDS`
  (`progress_rate.HOLD_SECONDS`: the first seconds of a unit say little about its speed). With no recorded rate the ETA is the live rate alone, with no live rate yet the recorded one.
  The same blend as `progress_rate.blended_rate` (a stored rate moving to the live one), with the weight taken from the task's own progress, because a unit is one task, not many.

The record is the heartbeat's: a unit's `unitbeat.Beat` carries `done`, `estimate` (the total) and `started`, so the Queue page, the API (`/units`) and the log read one set of numbers
(`unitwatch.view`). The recorded rate is items a second per command (`unitrates`).
"""

import math
from dataclasses import dataclass
from typing import Optional

from planetgen.queue import progress_rate

SHOW_SECONDS = 15.0
"""float: A task expected to take, or taking, longer than this shows its own bar (Boss, 2026-10-21)."""

LIVE_FULL_AT = 0.25
"""float: The share of a task done at which its ETA trusts the running rate alone; before it, the running rate's weight is the share done over this."""

HOLD_SECONDS = progress_rate.HOLD_SECONDS
"""float: The running rate is not used until the task has run this long."""

STORED, BLENDED, LIVE = "stored", "blended", "live"
SOURCES = (STORED, BLENDED, LIVE)
"""tuple: Where an ETA came from: the recorded rate, the blend of the two, the running rate alone."""


@dataclass(frozen=True)
class TaskProgress:
    """
    One task's bar.

    Attributes:
        show (bool): Whether the task shows a bar (expected or actually longer than `SHOW_SECONDS`).
        fraction (float, optional): `done / total` (0 to 1), `None` where the task has no total.
        eta_seconds (float, optional): Seconds left, `None` where no rate is known yet.
        source (str, optional): One of `SOURCES`, `None` with no ETA.
        rate (float, optional): The items a second the ETA used.
        expected_seconds (float, optional): How long the whole task takes at the recorded rate, `None` where none is recorded.
        weight (float): The running rate's share of the blended rate (0 to 1).
    """
    show: bool
    fraction: Optional[float] = None
    eta_seconds: Optional[float] = None
    source: Optional[str] = None
    rate: Optional[float] = None
    expected_seconds: Optional[float] = None
    weight: float = 0.0


def _good(value):
    return isinstance(value, (int, float)) and math.isfinite(value) and value > 0


def live_weight(fraction):
    """The running rate's share of the blended rate at `fraction` of the task done: it grows with the progress and is 1 from `LIVE_FULL_AT`."""
    if fraction is None or not math.isfinite(fraction) or fraction <= 0:
        return 0.0
    return min(1.0, fraction / LIVE_FULL_AT)


def progress(done, total, started, now, stored_rate=None):
    """
    The bar of a task.

    Args:
        done (float): The work finished, in the task's own items.
        total (float): The work expected (0 or `None` for a task with no counter: no bar and no ETA).
        started (float): When the task began, on the clock `now` is read from.
        now (float): The time.
        stored_rate (float, optional): Items a second the server recorded for this kind of work (`unitrates`).

    Returns:
        TaskProgress: See there.
    """
    ran = max(now - started, 0.0) if _good(now) or now == 0 else 0.0
    if not _good(total):
        return TaskProgress(show=ran > SHOW_SECONDS)
    done = min(max(done, 0.0), total) if isinstance(done, (int, float)) and math.isfinite(done) else 0.0
    fraction = done / total
    stored = stored_rate if _good(stored_rate) else None
    expected = total / stored if stored else None
    show = (expected is not None and expected > SHOW_SECONDS) or ran > SHOW_SECONDS
    if done >= total:
        return TaskProgress(show=show, fraction=1.0, eta_seconds=0.0, source=LIVE if stored is None else BLENDED, expected_seconds=expected, weight=1.0)
    live = done / ran if done > 0 and ran > 0 else None
    weight = live_weight(fraction) if live is not None and ran >= HOLD_SECONDS else 0.0
    if live is None or weight == 0.0:
        rate, source = (stored, STORED) if stored is not None else (live, LIVE if live is not None else None)
    elif stored is None or weight >= 1.0:
        rate, source = live, LIVE
    else:
        rate, source = weight * live + (1.0 - weight) * stored, BLENDED
    eta = (total - done) / rate if _good(rate) else None
    return TaskProgress(show=show, fraction=fraction, eta_seconds=eta, source=source if eta is not None else None, rate=rate if eta is not None else None,
                        expected_seconds=expected, weight=weight if source in (BLENDED, LIVE) else 0.0)


def as_dict(task):
    """A `TaskProgress` as the API and the pages show it (plain numbers; `fraction` as `percent`, 0 to 100)."""
    return {"show": task.show, "percent": None if task.fraction is None else round(100.0 * task.fraction, 1), "eta_seconds": task.eta_seconds, "source": task.source,
            "rate": task.rate, "expected_seconds": task.expected_seconds}
