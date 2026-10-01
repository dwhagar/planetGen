# stellarObjects/progressRate.py

"""
The rate behind every generation progress bar's ETA (PERF.7): work
finished per second as an exponentially weighted average over time, so
the estimate follows the run's current speed (a dense stretch of the
galaxy, more workers, a busy database) instead of rich's own short
window, and stays steady while many workers report at once.

`generate.py`'s progress bars keep one `DecayingRate` per bar and show
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
    Units finished per second, averaged with weight `exp(-age / tau)`.

    Each `add(amount)` folds in the speed since the previous one
    (`amount / interval`), weighted by how long that interval was, so a
    burst of tasks finishing together moves the average no more than the
    same tasks finishing evenly over the same time.

    Args:
        tau (float): The time constant, seconds.
        clock (callable): Returns the time in seconds (tests pass their own).
    """

    def __init__(self, tau=TIME_CONSTANT_SECONDS, clock=time.monotonic):
        self.tau = float(tau)
        self.clock = clock
        self.rate = None
        self.last = clock()
        # The average before the latest interval was folded in, so more
        # units finishing at that same instant join that interval.
        self._before = (None, self.last)
        self._amount = 0.0

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
            # Units that finished at the bar's very first instant (a
            # coarse clock) had no interval to fold into yet: they carry
            # over into this one instead of being dropped.
            if self.rate is not None or not self._amount:
                self._before = (self.rate, self.last)
                self._amount = 0.0
            self.last = now
        self._amount += amount
        rate, since = self._before
        interval = self.last - since
        if interval <= 0:
            return
        speed = self._amount / interval
        if rate is None:
            self.rate = speed
        else:
            self.rate = rate + (1.0 - math.exp(-interval / self.tau)) * (speed - rate)

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
