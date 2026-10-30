# planetgen/physics/substellar_ceiling.py

"""
The hottest a brown dwarf can be at its age, from the Baraffe et al. 2015 (BHAC15) isochrones (GEN.494, part 1).

A body below the hydrogen-burning limit cools for ever, so the hottest brown dwarf of an age is the one at the limit, and its
temperature falls as the population ages. Baraffe et al. 2015 (A&A 577, A42, section 3.3) place the limit at 0.07 solar masses and
give its temperature as about 2900 K at 1 Myr and 2950 K at 10 Myr; the table has it at 2897 K and 2946 K and then falling, 2794 K at
100 Myr, 2047 K at 1 Gyr and 1626 K at 10 Gyr. So a "late-M-like" brown dwarf (the row g1 of the atmosphere recipes, 2238 K to the
top of the map) exists only while it is young, and an old brown dwarf cannot be hotter than the mid-L dwarfs.

`ceiling_k(age_gyr)` is the temperature of the limit mass at that age, `HOTTEST_K` the highest it ever reaches (the hot edge of the
young late-M-like row), and `youngest_age_at_gyr(temperature_k)` the age at which a body of a temperature can no longer be a brown dwarf
because the limit mass has cooled below it. The table runs 0.5 Myr to 10 Gyr; outside it the functions answer None.

Limits. The limit mass is the paper's round 0.07 (the grid has 0.07 and 0.072, and the temperature of 0.072 is 0 to 8 percent higher:
2954 K at 10 Myr and 1790 K at 10 Gyr, so the true limit at old ages is a little above 0.07 and the ceiling is a lower bound there). Solar
metallicity only; the boundary moves with [Fe/H] (GEN.936). `cap_k(mass, age)` is what generation reads (GEN.494's wiring, Sources 2): `brown_dwarf` and `rogue_surface` cap the temperature of a body below the limit mass at `hottest_k(age)` and of one at or above it at `HOTTEST_K`,
in place of the 2800 K clip every age shared, and the hot edge of the atmosphere row g1 is `HOTTEST_K`. The cooling tracks (`baraffe_isochrones.cooling_track`, GEN.948) are themselves BHAC15, so the cap acts only at the peak
and where a sub-limit track runs a few kelvin over the limit mass's.
"""

import math

from planetgen.datasets import baraffe_isochrones, fetch
from planetgen.util.checks import finite_domain

LIMIT_MASS_SOLAR = 0.07
"""float: The hydrogen-burning limit of the table's solar-metallicity models (Baraffe et al. 2015, section 3.3: M = 0.07 solar masses,
Teff about 2900 K at 1 Myr and 2950 K at 10 Myr). [Stand-in SC1 of GEN.494, tier T1, tweak L0.]"""

HOTTEST_K = 2950.0
"""float: The highest temperature the limit mass reaches at any age: 2946 K at 10 Myr on the table, 2950 K as the paper rounds it. The
young late-M-like brown dwarf row ends here instead of the 2800 K of the rogue-body clip. [Stand-in SC2 of GEN.494, tier T1, tweak L0.]"""


@finite_domain()
def ceiling_k(age_gyr, grid=None):
    """The effective temperature, K, of a body at the hydrogen-burning limit at this age (the hottest brown dwarf there is); None when
    the age is outside the table's 0.5 Myr to 10 Gyr."""
    state = baraffe_isochrones.interpolate(LIMIT_MASS_SOLAR, age_gyr, grid)
    return None if state is None else state["teff_k"]


@finite_domain()
def hottest_k(age_gyr, grid=None):
    """The hottest a brown dwarf of this age can be, K: `ceiling_k` with the age held to the table's range (0.5 Myr to 10 Gyr). Past 10 Gyr the limit mass is still cooling, so the
    10 Gyr answer is an upper bound there. A run told to disregard missing data, with the table not downloaded, answers `HOTTEST_K` at every age (the fallback, recorded in
    `fetch.stand_ins_used`); any other run raises `fetch.DataMissing`."""
    if grid is None and fetch.locate(baraffe_isochrones.DATASET).use_stand_in("brown dwarf cap at the 2,950 K peak at every age"):
        return HOTTEST_K
    grid = baraffe_isochrones.isochrones() if grid is None else grid
    ages = [iso.age_gyr for iso in grid]
    return ceiling_k(min(max(age_gyr, min(ages)), max(ages)), grid)


@finite_domain()
def cap_k(mass_solar, age_gyr, grid=None):
    """The hottest the effective temperature of a body of this mass and age is allowed to be (GEN.494's wiring): `hottest_k` of the age below the hydrogen-burning limit; for a body at or
    above the limit mass (a very low mass star, which settles on the plateau and does not cool away) only the peak of the curve, `HOTTEST_K`."""
    return hottest_k(age_gyr, grid) if mass_solar < LIMIT_MASS_SOLAR else HOTTEST_K


@finite_domain()
def can_be_this_hot(temperature_k, age_gyr, grid=None):
    """Whether a brown dwarf of this age can be this hot: True when it is at or below the ceiling of the age; None outside the table."""
    ceiling = ceiling_k(age_gyr, grid)
    return None if ceiling is None else temperature_k <= ceiling


@finite_domain()
def oldest_age_gyr(temperature_k, grid=None):
    """The age, Gyr, past which no brown dwarf is as hot as this: the age at which the limit mass has cooled to this temperature (log-age
    interpolation on the table's ages after the peak). Zero when the temperature is above `HOTTEST_K` (never); the table's last age when
    the limit mass is still hotter at 10 Gyr."""
    grid = baraffe_isochrones.isochrones() if grid is None else grid
    if temperature_k > max(ceiling_k(iso.age_gyr, grid) for iso in grid):
        return 0.0
    ages = [iso.age_gyr for iso in grid]
    peak = max(range(len(ages)), key=lambda i: ceiling_k(ages[i], grid))
    previous = None
    for age in ages[peak:]:
        value = ceiling_k(age, grid)
        if value < temperature_k:
            if previous is None:
                return age
            (age_hot, hot), (age_cold, cold) = previous, (age, value)
            weight = (hot - temperature_k) / (hot - cold)
            return math.exp(math.log(age_hot) + weight * (math.log(age_cold) - math.log(age_hot)))
        previous = (age, value)
    return ages[-1]
