# planetgen/physics/imf.py

"""
The Initial Mass Function (GEN.926, GEN.927, GEN.928)
=====================================================

The one place the generator's initial mass function lives: `xi(M) dM`, the number of stars born in `dM` solar masses, a
continuous piecewise power law `k M^-alpha` between the edges in `tuning.IMF_BREAKS_SOL` with the slopes in
`tuning.IMF_SLOPES`. Those rows are the [R0705] function for individual stars of the Galactic disk (GEN.927): the
log-normal below one solar mass (`tuning.IMF_BELL_*`) cut into segments 0.05 dex wide, each a power law through the log-normal
at both ends, and a `M^-2.3` power law above one solar mass up to 150.

There are two curves (`SINGLE`, `PRIMARY`). `SINGLE` is `xi`: every star born, companions included. `PRIMARY` is the mass function
of the stars that are the primary of a system, derived from `xi` and the multiplicity rows (`generation.primary_mass_function`,
GEN.928); the stars of the systems together are `xi`. A star drawn as the primary of a system, or counted as a placed system,
is on `PRIMARY`; a rate per star born (the ionizing stars, the black holes, the mass returned) is on `SINGLE`.

Everything that needs the function asks here: the draw of a star's mass (`sample`), the number or mass of stars in a range
(`number`, `integral`), the segments themselves for a table built over a mass grid (`segments`, `segment_at`), and the
primitives of a single power-law segment (`power_law_integral`, `sample_power_law`). The same function used to be written
out in six modules (the draw and the count in `stellar_evolution`, the moments in `system_mass`, the stars per solar mass
in `hii`, the mass return in `galaxy.disk_history`, the grid in `star_population`, the cluster check in `mathcheck`), each
joining the segments its own way; a change of the curve is now a change of the rows and of this module.

`mathcheck` keeps its own written-out copy of the [R0705] function on purpose: it is the check that this one still
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


SINGLE = "single"
"""str: The curve of every star born, companions included (`xi(M)`, the rows of `tuning`)."""

PRIMARY = "primary"
"""str: The curve of the stars that are the primary of a system, derived from `SINGLE` and the companions the heavier primaries
bring (`generation.primary_mass_function`, GEN.928). Same scale as `SINGLE`, so the ratio of their numbers is the systems per star."""


@functools.lru_cache(maxsize=None)
def _segments(breaks, slopes):
    norms = [1.0]
    for i in range(1, len(slopes)):
        norms.append(norms[-1] * breaks[i] ** (slopes[i] - slopes[i - 1]))
    return tuple((breaks[i], breaks[i + 1], slopes[i], norms[i]) for i in range(len(slopes)))


@functools.lru_cache(maxsize=None)
def _primary_segments():
    # Deferred import: the derivation reads the companion tables of the generator, which import this module.
    from planetgen.generation import primary_mass_function
    return primary_mass_function.segments()


def segments(curve=SINGLE):
    """`((low, high, alpha, k), ...)`: the function's segments, each `k M^-alpha` between `low` and `high` solar masses,
    joined continuously (`k` of the first of `SINGLE` is 1)."""
    if curve == PRIMARY:
        return _primary_segments()
    return _segments(tuple(tuning.IMF_BREAKS_SOL), tuple(tuning.IMF_SLOPES))


def segment_at(mass_sol, curve=SINGLE):
    """`(low, high, alpha, k)` of the segment that holds `mass_sol` (a mass on an edge belongs to the segment above it,
    the top one for the top edge)."""
    found = segments(curve)
    for segment in found:
        if mass_sol < segment[1]:
            return segment
    return found[-1]


def density(mass, breaks=None, slopes=None, curve=SINGLE):
    """`xi(M)` at each of the masses in the array `mass` (solar masses; 0 outside the edges), for the function of the rows
    in `tuning` or, given them, of other `breaks` and `slopes` (a check against another function, as `galaxy.disk_history`
    does for the mass a population returns), or of the `curve` named."""
    found = segments(curve) if breaks is None and slopes is None else _segments(tuple(breaks), tuple(slopes))
    out = np.zeros_like(mass)
    for low, high, alpha, k in found:
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


def integral(power, low, high, curve=SINGLE):
    """The integral of `M^power xi(M) dM` from `low` to `high` solar masses: the number of stars for power 0, their mass
    for power 1. Only the ratio of two integrals means anything (the function has no overall scale), within one `curve`
    or between the two (`PRIMARY` over `SINGLE` is the systems per star)."""
    total = 0.0
    for seg_low, seg_high, alpha, k in segments(curve):
        a, b = max(low, seg_low), min(high, seg_high)
        if a >= b:
            continue
        exponent = power + 1.0 - alpha
        total += k * (math.log(b / a) if exponent == 0.0 else (b ** exponent - a ** exponent) / exponent)
    return total


def number(low_sol, high_sol, curve=SINGLE):
    """The number of stars born in `[low_sol, high_sol]` solar masses per unit of the function; 0 for an empty range."""
    return integral(0.0, low_sol, high_sol, curve)


def sample(min_mass_sol=None, max_mass_sol=None, rng=draw, curve=SINGLE):
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
    for seg_low, seg_high, alpha, k in segments(curve):
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
