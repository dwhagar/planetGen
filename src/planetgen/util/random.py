# planetgen/util/random.py

"""
Random
======

Seeded sampling helpers: power-law and bounded-bell draws, and the one
log-uniform draw. Each uses the `random` module's global generator (or the
`rng` it is given), so a seeded galaxy draws the same values.
"""

import math
import random

from planetgen.util.checks import finite_domain


def log_uniform(low, high, rng=random):
    """A value drawn uniformly in log space between `low` and `high`: one
    `rng.uniform` draw (the `random` module's global generator unless an
    `rng` is given)."""
    return math.exp(rng.uniform(math.log(low), math.log(high)))


def power_law_share(low, high, slope):
    """The unnormalized weight of [low, high] under dN/dlogM ~ M^-slope
    (`slope` > 0): `low^-slope - high^-slope`."""
    return low ** -slope - high ** -slope


def sample_power_law(low, high, slope):
    """A value in [low, high] drawn from dN/dlogM ~ M^-slope (`slope` > 0)
    by inverting its cumulative distribution: one `random.random()`."""
    low_term, high_term = low ** -slope, high ** -slope
    return (low_term + random.random() * (high_term - low_term)) ** (-1 / slope)


@finite_domain(clamped=("mode_fraction",))
def sample_bounded_bell(min_val, max_val, mode_fraction, spread_divisor=3.0, max_attempts=1000):
    """
    Draws a random value in [min_val, max_val] from a bounded bell-curve
    (Gaussian) distribution peaking at `min_val + mode_fraction * (max_val -
    min_val)`, instead of a flat uniform draw across the whole range.

    `mode_fraction` (0.0-1.0) is "what fraction through the available range
    is the statistically most common (modal) value" -- e.g. a class whose
    real-world single-body analog sits 27% of the way from its declared
    radius minimum to its maximum uses `mode_fraction=0.27` so generated
    instances cluster around that real value instead of being spread flatly
    across the whole declared range (see `tuning.PLANET_CLASSES`'
    own `size_mode` values and the real-world analogs their docstrings
    cite).

    The standard deviation self-adjusts to whichever bound is nearer the
    mode (reaching it in `spread_divisor` steps, ~3-sigma by default) so a
    mode pinned near one edge of the range still produces a legible bell
    shape instead of wasting most of its probability mass outside the
    range entirely. Uses rejection sampling (redraw until the result falls
    in [min_val, max_val]) rather than clamping a raw Gaussian draw to the
    bounds, which would pile spillover probability mass up at the boundary
    and destroy the bell shape there; `max_attempts` is a termination
    safety net (falls back to the clamped mean if never satisfied, which
    the self-adjusting spread should make effectively unreachable in
    practice for any reasonable `mode_fraction`).

    Args:
        min_val (float): Lower bound of the available range for this draw
                         (not necessarily a class's full declared range --
                         e.g. a moon's available radius may be further
                         capped by its parent's Hill sphere).
        max_val (float): Upper bound of the available range for this draw.
        mode_fraction (float): 0.0-1.0, how far through [min_val, max_val]
                               the distribution's peak sits.
        spread_divisor (float, optional): How many standard deviations
                                          reach the nearer bound from the
                                          mode. Defaults to 3.0.
        max_attempts (int, optional): Rejection-sampling attempt cap.

    Returns:
        float: The sampled value, guaranteed within [min_val, max_val].
    """
    if max_val <= min_val:
        return min_val
    mode_fraction = min(1.0, max(0.0, mode_fraction))
    mean = min_val + mode_fraction * (max_val - min_val)
    span = max_val - min_val
    # Floored so a mode pinned exactly at an edge (fraction 0.0 or 1.0)
    # doesn't collapse the spread to zero.
    nearer_bound_distance = max(min(mean - min_val, max_val - mean), span * 0.02)
    stdev = nearer_bound_distance / spread_divisor
    for _ in range(max_attempts):
        value = random.gauss(mean, stdev)
        if min_val <= value <= max_val:
            return value
    return min(max(mean, min_val), max_val)
