# planetgen/physics/remnant_progenitors.py

"""
Which massive stars leave a black hole and which a neutron star (GEN.554; docs/design/systems/explosion-islands.md)
===========================================================================================================

Every star born with 8 to 150 solar masses (`tuning.CORE_COLLAPSE_MASS_RANGE_SOL`) ends as one or the other, and which
depends on the star's initial mass in windows, not at one cut: "islands of 'explodability' in a sea of black hole
formation" [R4251]. The share of the massive stars that fail to explode and leave a black hole is the
galaxy's failed share (`black_hole_formation.failed_share`); the windows are the initial masses that make up that share:

* the two verified islands of `tuning.EXPLOSION_FAILURE_ISLANDS_SOL` (23 to 25 and 35 to 45 solar masses), shrunk about
  their centres when the share is smaller than theirs (about 7% of the massive stars);
* when it is larger, the islands and the top of the mass range down to the mass that makes up the share.

The neutron stars are the massive stars outside the windows. The progenitor of a remnant is drawn from the IMF over its
kind's windows (`sample_mass_sol`), and a black hole whose progenitor lies in a window failed outright (a direct collapse).
The windows follow the share alone; a draw with no galaxy around it takes the share of the middle draw at solar metallicity
(`reference_failed_share`).

The age of a remnant follows from its progenitor (`collapse_age_gy`).
"""

import functools

from planetgen import tuning
from planetgen.physics import imf, stellar_evolution
from planetgen.util import draw

NEUTRON_STAR = "neutron-star"
BLACK_HOLE = "black-hole"
KINDS = (NEUTRON_STAR, BLACK_HOLE)
"""tuple: The kinds of remnant a massive star leaves (`remnant_mass.KINDS`)."""

_SHARE_DIGITS = 6
_BISECTIONS = 60


def reference_failed_share():
    """The failed share of a draw with no galaxy around it: the solar share of the reference structure, from its draws
    (`black_hole_formation.reference_failed_shares`)."""
    from planetgen.galaxy import structure   # a late import: the structure imports the physics of the shares
    return structure.REFERENCE_STRUCTURE.black_hole_solar_share


def _number(low, high):
    return imf.number(low, high)


def _union(intervals):
    """`intervals` sorted, with those that overlap joined."""
    merged = []
    for low, high in sorted(intervals):
        if merged and low <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], high))
        else:
            merged.append((low, high))
    return tuple(merged)


def _count(intervals):
    return sum(_number(low, high) for low, high in intervals)


def _scaled_islands(factor):
    """The islands with their log widths times `factor` about their geometric centres."""
    scaled = []
    for low, high in tuning.EXPLOSION_FAILURE_ISLANDS_SOL:
        centre = (low * high) ** 0.5
        scaled.append((centre * (low / centre) ** factor, centre * (high / centre) ** factor))
    return scaled


@functools.lru_cache(maxsize=512)
def _failure_windows(share):
    bottom, top = tuning.CORE_COLLAPSE_MASS_RANGE_SOL
    target = share * _number(bottom, top)
    if target <= 0.0:
        return ()
    if target >= _number(bottom, top):
        return ((bottom, top),)
    if target <= _count(tuning.EXPLOSION_FAILURE_ISLANDS_SOL) * (1.0 + 1.0e-4):
        low, high = 0.0, 1.0
        for _ in range(_BISECTIONS):
            middle = 0.5 * (low + high)
            if _count(_scaled_islands(middle)) < target:
                low = middle
            else:
                high = middle
        return _union(_scaled_islands(0.5 * (low + high)))
    low, high = bottom, top
    for _ in range(_BISECTIONS):
        middle = 0.5 * (low + high)
        if _count(_union(list(tuning.EXPLOSION_FAILURE_ISLANDS_SOL) + [(middle, top)])) > target:
            low = middle
        else:
            high = middle
    return _union(list(tuning.EXPLOSION_FAILURE_ISLANDS_SOL) + [(0.5 * (low + high), top)])


def failure_windows(failed_share=None):
    """
    The windows `((low, high), ...)` of initial mass, solar masses, in which a massive star fails to explode and leaves a
    black hole, for a population that fails `failed_share` of its massive stars (the reference share when omitted):
    their IMF share of the 8 to 150 solar mass stars is `failed_share`.
    """
    share = reference_failed_share() if failed_share is None else failed_share
    return _failure_windows(round(min(max(share, 0.0), 1.0), _SHARE_DIGITS))


def kind_windows(kind, failed_share=None):
    """The windows of initial mass that end as `kind`: the failure windows for a black hole, what lies between them in 8 to
    150 solar masses for a neutron star."""
    if kind not in KINDS:
        raise ValueError(f"{kind!r} is not a remnant of a massive star ({', '.join(KINDS)})")
    windows = failure_windows(failed_share)
    if kind == BLACK_HOLE:
        return windows
    bottom, top = tuning.CORE_COLLAPSE_MASS_RANGE_SOL
    complement, edge = [], bottom
    for low, high in windows:
        if low > edge:
            complement.append((edge, low))
        edge = high
    if edge < top:
        complement.append((edge, top))
    return tuple(complement)


def fails(initial_mass_sol, failed_share=None):
    """Whether a star of this initial mass fails to explode (a direct collapse to a black hole)."""
    return any(low <= initial_mass_sol <= high for low, high in failure_windows(failed_share))


def number(kind, lowest_sol=0.0, failed_share=None):
    """The stars per unit of the IMF whose progenitors end as `kind`, no lighter than `lowest_sol`."""
    return sum(_number(max(low, lowest_sol), high) for low, high in kind_windows(kind, failed_share))


def sample_mass_sol(kind, lowest_sol=0.0, rng=draw, failed_share=None):
    """
    The initial mass, solar masses, of a progenitor of a `kind` that has had time to form: the IMF over the kind's windows,
    no lighter than `lowest_sol`. `None` when no window reaches that high.
    """
    pieces = [(max(low, lowest_sol), high) for low, high in kind_windows(kind, failed_share)]
    pieces = [(low, high) for low, high in pieces if low < high]
    if not pieces:
        return None
    weights = [_number(low, high) for low, high in pieces]
    pick = rng.random() * sum(weights)
    for (low, high), weight in zip(pieces, weights):
        if pick < weight:
            break
        pick -= weight
    return imf.sample(low, high, rng)


def collapse_age_gy(kind, age_bias=None, rng=draw, population=None, age_law=None, failed_share=None):
    """
    The age, Gy, of a neutron star or stellar black hole in a population (GEN.553): the time since the star that made it
    collapsed. The progenitor is drawn from the IMF over the kind's windows (`sample_mass_sol`), no lighter than the lightest
    star that can have collapsed inside the population's window; the star's birth age is drawn as any star's is
    (`sample_star_age_gy`: the population's window and age law) but no younger than the progenitor took to collapse,
    and the remnant's age is the birth age less that time. A window where no progenitor can have collapsed gives the
    youngest age a remnant has (`tuning.COMPACT_REMNANT_AGE_RANGE_GY`).

    Args:
        kind (str): "neutron-star" or "black-hole" (`remnant_mass.KINDS`).
        age_bias (str, optional): "young", "old" or None, as `sample_star_age_gy`.
        rng: The random source.
        population (str, optional): The stellar population (`population_age_range_gy`); None: the whole disk's history.
        age_law (AgeLaw, optional): The density of birth ages inside the window (GEN.546); None: uniform.
        failed_share (float, optional): The failed share of the place (`failure_windows`).
    """
    window_low, window_high = stellar_evolution.age_window_gy(*stellar_evolution.population_age_range_gy(population), age_bias)
    floor = tuning.COMPACT_REMNANT_AGE_RANGE_GY[0]
    progenitor = sample_mass_sol(kind, stellar_evolution.youngest_progenitor_sol(window_high * 1.0e9), rng, failed_share)
    if progenitor is None:
        return floor
    formed_gy = stellar_evolution.remnant_formation_time_gy(progenitor)
    window_low = max(window_low, formed_gy)
    born_gy = rng.uniform(window_low, window_high) if age_law is None else age_law.sample(window_low, window_high, rng)
    return max(born_gy - formed_gy, floor)
