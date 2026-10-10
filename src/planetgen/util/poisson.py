# planetgen/util/poisson.py

"""
Poisson counts for any mean in constant time (PERF.60).

Below `tuning.POISSON_REJECTION_MEAN` the count is Knuth's exact
multiplicative method (a handful of uniform draws); from there up it is
Hormann's transformed rejection with squeeze ("PTRS", Insurance: Mathematics
and Economics 12, 1993), which is exact too and takes about the same few
draws whether the mean is 10 or 10 million. The object-first scatter
(`generation/bright_stars.py`) draws one count per layer, with means in the
thousands to millions, which Knuth's O(mean) method cannot do.

Every draw comes from the `rng` it is given (a `draw.Stream`, or the `draw`
module's run stream), so the same seed gives the same counts (GEN.56).
"""

import math

from planetgen import tuning
from planetgen.util import draw


def poisson_count(mean, rng=draw):
    """
    A Poisson-distributed non-negative integer with the given mean.

    Args:
        mean (float): The Poisson mean (lambda); at or below 0 always gives 0.
        rng (draw.Stream or the draw module): The generator to draw from.

    Returns:
        int: A Poisson-distributed sample.

    Raises:
        ValueError: If `mean` is not finite.
    """
    if not math.isfinite(mean):
        raise ValueError(f"poisson_count: mean must be finite, got {mean!r}")
    if mean <= 0.0:
        return 0
    if mean < tuning.POISSON_REJECTION_MEAN:
        return _knuth(mean, rng)
    return _transformed_rejection(mean, rng)


def _knuth(mean, rng):
    threshold = math.exp(-mean)
    count = 0
    product = 1.0
    while True:
        product *= rng.random()
        if product <= threshold:
            return count
        count += 1


def _transformed_rejection(mean, rng):
    """Hormann's PTRS, valid for a mean of 10 or more."""
    slam = math.sqrt(mean)
    loglam = math.log(mean)
    b = 0.931 + 2.53 * slam
    a = -0.059 + 0.02483 * b
    inv_alpha = 1.1239 + 1.1328 / (b - 3.4)
    vr = 0.9277 - 3.6224 / (b - 2.0)
    log_inv_alpha = math.log(inv_alpha)
    while True:
        u = rng.random() - 0.5
        v = rng.random()
        us = 0.5 - abs(u)
        k = math.floor((2.0 * a / us + b) * u + mean + 0.43)
        if us >= 0.07 and v <= vr:
            return k
        if k < 0 or (us < 0.013 and v > us):
            continue
        if v > 0.0 and math.log(v) + log_inv_alpha - math.log(a / (us * us) + b) <= -mean + k * loglam - math.lgamma(k + 1.0):
            return k
