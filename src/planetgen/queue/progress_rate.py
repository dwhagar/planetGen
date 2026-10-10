# planetgen/queue/progress_rate.py

"""
The rate behind every generation progress bar's ETA (PERF.7, PERF.37): work
finished per second as an exponentially weighted average over time, so
the estimate follows the run's current speed (a dense stretch of the
galaxy, more workers, a busy database) instead of rich's own short
window, and stays steady while many workers report at once.

`planetgen`'s progress bars keep one `DecayingRate` per bar and show
`eta()`; `progressFile` writes the same rate and ETA for the Generate
page, so the terminal and the page agree.

PERF.33: a bar may start from the rate this server recorded for the same
kind of work (`prior`, PERF.32) and blend the live rate into it as units
finish: `rate = w * live + (1 - w) * prior` with `w = n / (n + 15)`, the
live side held back until `HOLD_UNITS` units and `HOLD_SECONDS` seconds.
With a prior the ETA exists before the first unit is done; without one the
bar behaves as before (live only, no ETA until a unit finishes). One pure
function, `blended_rate`, serves the terminal bar, the Generate and Queue
pages and the stall rule.
"""

import math
import time

TIME_CONSTANT_SECONDS = 60.0
"""float: How fast the average forgets: work finished this long ago
counts about a third (1/e) as much as work finished just now."""

PRIOR_WEIGHT_UNITS = 15
"""int: Finished units at which the live rate and the recorded one weigh
the same (`w = n / (n + 15)`)."""

HOLD_UNITS = 5
HOLD_SECONDS = 20.0
"""The live rate is not used with a prior until this many units have
finished and this long has passed (a first unit finishing says little)."""

ETA_LOW_FACTOR = 0.88
ETA_HIGH_FACTOR = 1.15
"""Where the middle 80% of estimates lie, in times the true time left
(docs/design/performance-eta-queue-and-caching.md 2.3)."""

STALL_FLOOR_SECONDS = 60.0
"""float: No unit finishing for at least this long (and 3 times the longest
unit so far) is a stall: the ETA stops counting down."""


def time_constant(mean_task_seconds=None, workers=1):
    """The decay time constant for a bar of tasks that take `mean_task_seconds`
    each on `workers` workers: `max(60 s, 20 x mean task seconds / workers)`,
    so the average spans many tasks however slow they are."""
    if mean_task_seconds is None or not _finite(mean_task_seconds) or mean_task_seconds <= 0:
        return TIME_CONSTANT_SECONDS
    return max(TIME_CONSTANT_SECONDS, 20.0 * mean_task_seconds / max(int(workers or 1), 1))


def blended_rate(live, n, prior, elapsed=None):
    """
    The rate to estimate with: the live rate `live` after `n` finished units
    blended with the recorded `prior` by `w = n / (n + 15)`, the live side
    held back until `HOLD_UNITS` units are done and `elapsed` seconds have
    passed (`HOLD_SECONDS`; `None` skips that test). Either may be `None`:
    no prior gives the live rate, no live rate gives the prior, neither
    gives `None`. A pure function.
    """
    prior = prior if prior is not None and _finite(prior) and prior > 0 else None
    live = live if live is not None and _finite(live) and live > 0 else None
    if prior is None:
        return live
    if live is None or n < HOLD_UNITS or (elapsed is not None and elapsed < HOLD_SECONDS):
        return prior
    weight = n / (n + PRIOR_WEIGHT_UNITS)
    return weight * live + (1.0 - weight) * prior


def eta_range(eta):
    """`(low, high)` seconds around `eta` (`ETA_LOW_FACTOR`, `ETA_HIGH_FACTOR`), or `None`."""
    if eta is None or not _finite(eta):
        return None
    return eta * ETA_LOW_FACTOR, eta * ETA_HIGH_FACTOR


class DecayingRate:
    """
    Units finished per second, as two decayed sums (PERF.37): `N`, the
    units finished with weight `exp(-age / tau)`, and `D`, the time that
    has passed with the same weight; the rate is `N / D`.

    Each `add(amount)` spreads `amount` evenly over the interval since the
    previous one, so a burst of tasks finishing together (several workers
    at once) counts exactly as much as the same tasks finishing evenly over
    the same time. Because the average is a ratio of sums rather than a
    running average that started at the first single completion, the early
    rate is not up to twice too low and the ETA not up to twice too long.

    Args:
        tau (float): The time constant, seconds (`time_constant`).
        clock (callable): Returns the time in seconds (tests pass their own).
        prior (float, optional): The rate recorded for this kind of work
            (units a second, PERF.32), blended in by `blended_rate`.
    """

    def __init__(self, tau=TIME_CONSTANT_SECONDS, clock=time.monotonic, prior=None):
        self.tau = float(tau)
        self.clock = clock
        self.prior = prior if prior is not None and _finite(prior) and prior > 0 else None
        self.live = None
        self.n = 0
        self.started = clock()
        self.longest_gap = 0.0
        self.last = self.started
        self._units = 0.0
        self._time = 0.0
        # The latest interval's average weight (decayed time over time), so
        # more units finishing at that same instant join that interval.
        self._share = None
        # Units finished before any time has passed carry into the first interval.
        self._pending = 0.0

    def add(self, amount, now=None):
        """Records `amount` more units finished (at `now`). A NaN or
        infinite amount or time is ignored, so one bad report can't turn
        the rate (and every later ETA) into NaN; a time earlier than the
        last one (a clock stepped back) joins the latest interval."""
        if not _finite(amount) or amount <= 0:
            return
        now = self.clock() if now is None else now
        if not _finite(now):
            return
        self.n += 1
        if now > self.last:
            interval = now - self.last
            self.longest_gap = max(self.longest_gap, interval)
            decay = math.exp(-interval / self.tau)
            weight = self.tau * (1.0 - decay)
            self._share = weight / interval
            self._time = self._time * decay + weight
            self._units = self._units * decay + (self._pending + amount) * self._share
            self._pending = 0.0
            self.last = now
        elif self._share is None:
            self._pending += amount
            return
        else:
            self._units += amount * self._share
        self.live = self._units / self._time

    @property
    def rate(self):
        """Units a second to estimate with: the live rate, blended with
        the recorded one when there is a prior (`blended_rate`); `None`
        while neither is known."""
        if self.prior is None:
            return self.live
        return blended_rate(self.live, self.n, self.prior, self.clock() - self.started)

    @property
    def source(self):
        """Where the rate comes from: `"live"`, `"recorded"` (the prior alone, or
        blended in), or `None`."""
        if self.rate is None:
            return None
        return "recorded" if self.prior is not None else "live"

    def stalled(self, now=None):
        """Whether no unit has finished for `max(60 s, 3 x the longest
        gap so far)`: the ETA then stops counting down."""
        now = self.clock() if now is None else now
        return _finite(now) and now - self.last > max(STALL_FLOOR_SECONDS, 3.0 * self.longest_gap)

    def eta(self, remaining, now=None):
        """
        Seconds until `remaining` more units are done at the current
        rate, counting down between updates; `None` until the first
        unit is done (or at once, with a prior), 0 when nothing remains.
        """
        if remaining is None or not _finite(remaining):
            return None
        if remaining <= 0:
            return 0.0
        rate = self.rate
        if not rate:
            return None
        now = self.clock() if now is None else now
        left = remaining
        # A clock stepped back (or a NaN one) doesn't add time to the
        # estimate: it just hasn't started counting down yet.
        elapsed = now - self.last if _finite(now) and now > self.last else 0.0
        if self.stalled(now):
            elapsed = 0.0  # a stall: hold the estimate rather than count it down
        # Counts down between updates, but never below the time the rest
        # takes once the next unit is overdue, so a slow unit holds the
        # estimate rather than letting it reach zero early.
        return max(0.0, left / rate - elapsed, (left - 1) / rate)


def _finite(value):
    try:
        return math.isfinite(value)
    except TypeError:
        return False
