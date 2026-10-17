# planetgen/queue/hourdisplay.py

"""
The time left shown as whole hours (UX.3), by the rule of section 2.4 of docs/design/performance-eta-queue-and-caching.md.

Boss asked for "an ETA until it will be finished to the nearest hour rounding up". Rounding each estimate up on
its own makes the figure jump: on a 10 hour run `ceil(eta / 1 h)` rose 9 to 18 times and showed fewer hours than
were left for 15% to 25% of the time. So the figure shown keeps a memory of what it showed last:

* a fresh figure is `ceil(ETA_HIGH_FACTOR x eta / 1 h)`, never below 1: the top of the likely range, so the run
  seldom outlasts it;
* it **falls freely** as the estimate falls;
* it **rises** only when even the bottom of the range, `ETA_LOW_FACTOR x eta`, is longer than the figure shown, so
  an estimate that merely wobbles never raises it.

Pure: the caller keeps the last figure (the banner keeps it in Redis, so every visitor sees the same one).
"""

import math

from planetgen.queue.progress_rate import ETA_HIGH_FACTOR, ETA_LOW_FACTOR

HOUR_SECONDS = 3600.0
"""float: Seconds in an hour."""


def hours_left(eta_seconds, shown=None):
    """
    The whole hours to show for a run with `eta_seconds` left.

    Args:
        eta_seconds (float): The run's time left (`runeta.estimate`), at least 0.
        shown (int, optional): The figure shown last time for this run, or `None` for the first.

    Returns:
        int: At least 1.
    """
    fresh = max(1, math.ceil(ETA_HIGH_FACTOR * eta_seconds / HOUR_SECONDS))
    if shown is None or fresh <= shown:
        return fresh
    if ETA_LOW_FACTOR * eta_seconds / HOUR_SECONDS > shown:
        return fresh
    return shown
