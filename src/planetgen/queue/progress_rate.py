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
"""

import math
import time

TIME_CONSTANT_SECONDS = 60.0
"""float: How fast the average forgets: work finished this long ago
counts about a third (1/e) as much as work finished just now."""


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
        tau (float): The time constant, seconds.
        clock (callable): Returns the time in seconds (tests pass their own).
    """

    def __init__(self, tau=TIME_CONSTANT_SECONDS, clock=time.monotonic):
        self.tau = float(tau)
        self.clock = clock
        self.rate = None
        self.last = clock()
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
        if now > self.last:
            interval = now - self.last
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
        self.rate = self._units / self._time

    def eta(self, remaining, now=None):
        """
        Seconds until `remaining` more units are done at the current
        rate, counting down between updates; `None` until the first
        unit is done, 0 when nothing remains.
        """
        if remaining is None or not _finite(remaining):
            return None
        if remaining <= 0:
            return 0.0
        if not self.rate:
            return None
        now = self.clock() if now is None else now
        left = remaining
        # A clock stepped back (or a NaN one) doesn't add time to the
        # estimate: it just hasn't started counting down yet.
        elapsed = now - self.last if _finite(now) and now > self.last else 0.0
        # Counts down between updates, but never below the time the rest
        # takes once the next unit is overdue, so a slow unit holds the
        # estimate rather than letting it reach zero early.
        return max(0.0, left / self.rate - elapsed, (left - 1) / self.rate)


def _finite(value):
    try:
        return math.isfinite(value)
    except TypeError:
        return False
