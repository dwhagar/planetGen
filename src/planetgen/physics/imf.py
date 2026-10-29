# planetgen/physics/imf.py

"""
The Initial Mass Function (GEN.926)
===================================

The one place the generator's initial mass function lives: `xi(M) dM`, the number of stars born in `dM` solar masses, a
continuous piecewise power law `k M^-alpha` between the edges in `tuning.IMF_BREAKS_SOL` with the slopes in
`tuning.IMF_SLOPES` (the rows the mass function is set by).

Everything that needs the function asks here: the draw of a star's mass (`sample`), the number or mass of stars in a range
(`number`, `integral`), the segments themselves for a table built over a mass grid (`segments`, `segment_at`), and the
primitives of a single power-law segment (`power_law_integral`, `sample_power_law`). The same function used to be written
out in six modules (the draw and the count in `stellar_evolution`, the moments in `system_mass`, the stars per solar mass
in `hii`, the mass return in `galaxy.disk_history`, the grid in `star_population`, the cluster check in `mathcheck`), each
joining the segments its own way; a change of the curve is now a change of the rows and of this module.

`mathcheck` keeps its own written-out copy of the Kroupa (2001) function on purpose: it is the check that this one still
is what the literature says.
"""

import functools
import math

import numpy as np

from planetgen import tuning
from planetgen.util import draw


def lowest_mass_sol():
    """The lightest initial mass, in solar masses, the function reaches (the floor the generator makes a star at)."""
    return tuning.IMF_BREAKS_SOL[0]


def highest_mass_sol():
    """The heaviest initial mass, in solar masses, the function reaches."""
    return tuning.IMF_BREAKS_SOL[-1]


@functools.lru_cache(maxsize=None)
def _segments(breaks, slopes):
    norms = [1.0]
    for i in range(1, len(slopes)):
        norms.append(norms[-1] * breaks[i] ** (slopes[i] - slopes[i - 1]))
    return tuple((breaks[i], breaks[i + 1], slopes[i], norms[i]) for i in range(len(slopes)))


def segments():
    """`((low, high, alpha, k), ...)`: the function's segments, each `k M^-alpha` between `low` and `high` solar masses,
    joined continuously (`k` of the first is 1)."""
    return _segments(tuple(tuning.IMF_BREAKS_SOL), tuple(tuning.IMF_SLOPES))


def segment_at(mass_sol):
    """`(low, high, alpha, k)` of the segment that holds `mass_sol` (a mass on an edge belongs to the segment above it,
    the top one for the top edge)."""
    found = segments()
    for segment in found:
        if mass_sol < segment[1]:
            return segment
    return found[-1]


def density(mass, breaks=None, slopes=None):
    """`xi(M)` at each of the masses in the array `mass` (solar masses; 0 outside the edges), for the function of the rows
    in `tuning` or, given them, of other `breaks` and `slopes` (a check against another function, as `galaxy.disk_history`
    does for the mass a population returns)."""
    breaks = tuning.IMF_BREAKS_SOL if breaks is None else breaks
    slopes = tuning.IMF_SLOPES if slopes is None else slopes
    out = np.zeros_like(mass)
    for low, high, alpha, k in _segments(tuple(breaks), tuple(slopes)):
        inside = (mass >= low) & (mass <= high)
        out[inside] = k * mass[inside] ** -alpha
    return out


def power_law_integral(a, b, alpha):
    """The integral of `M^-alpha` from `a` to `b`."""
    if alpha == 1.0:
        return math.log(b / a)
    return (b ** (1 - alpha) - a ** (1 - alpha)) / (1 - alpha)


def sample_power_law(a, b, alpha, rng=draw):
    """A mass in `[a, b]` drawn from `M^-alpha` by inverting its cumulative distribution (one `rng.random()`)."""
    u = rng.random()
    if alpha == 1.0:
        return a * (b / a) ** u
    p = 1 - alpha
    return (a ** p + u * (b ** p - a ** p)) ** (1 / p)


def integral(power, low, high):
    """The integral of `M^power xi(M) dM` from `low` to `high` solar masses: the number of stars for power 0, their mass
    for power 1. Only the ratio of two integrals means anything (the function has no overall scale)."""
    total = 0.0
    for seg_low, seg_high, alpha, k in segments():
        a, b = max(low, seg_low), min(high, seg_high)
        if a >= b:
            continue
        exponent = power + 1.0 - alpha
        total += k * (math.log(b / a) if exponent == 0.0 else (b ** exponent - a ** exponent) / exponent)
    return total


def number(low_sol, high_sol):
    """The number of stars born in `[low_sol, high_sol]` solar masses per unit of the function; 0 for an empty range."""
    return integral(0.0, low_sol, high_sol)


def sample(min_mass_sol=None, max_mass_sol=None, rng=draw):
    """
    Draws an initial mass (solar masses) from the function, optionally truncated to `[min_mass_sol, max_mass_sol]`.

    Each segment's number weight is its power-law integral (the segments joined continuously at the breaks); a segment is
    picked by weight and the mass drawn from it by inverse transform. Two draws from `rng`.
    """
    lo = lowest_mass_sol() if min_mass_sol is None else max(min_mass_sol, lowest_mass_sol())
    hi = highest_mass_sol() if max_mass_sol is None else min(max_mass_sol, highest_mass_sol())
    if not lo < hi:
        raise ValueError(f"empty IMF range [{lo}, {hi}] Msun")

    weighted = []
    for seg_low, seg_high, alpha, k in segments():
        a, b = max(seg_low, lo), min(seg_high, hi)
        if a < b:
            weighted.append((a, b, alpha, k * power_law_integral(a, b, alpha)))

    total = sum(weight for *_, weight in weighted)
    pick = rng.random() * total
    for a, b, alpha, weight in weighted:
        if pick < weight:
            break
        pick -= weight
    return sample_power_law(a, b, alpha, rng)
