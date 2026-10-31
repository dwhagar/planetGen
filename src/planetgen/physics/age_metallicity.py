# planetgen/physics/age_metallicity.py

"""
The age-metallicity relation of the disc, by branch (GEN.744 part 2, GEN.944 part 1).

A star's [Fe/H] was drawn from its population alone (`metallicity.sample_fe_h`), independent of its age inside the population's window.
In the solar neighbourhood the two are tied, and differently on the two branches (Haywood et al. 2013, A&A 560, A109, Section 3.3 and
Fig. 9):

- **The thin disc** (`stellar_abundances.BRANCH_THIN`) gains metals slowly, 0.025 dex per Gyr, with a wide scatter.
- **The thick disc** (`BRANCH_HIGH_ALPHA`, older than `tuning.STELLAR_THICK_START_GY`, 8 Gyr) shows a tight relation that rises 0.15
  dex per Gyr over the 4 to 5 Gyr it took to form, "also implying a decrease by a factor of 5-6 in the production of iron after 8 Gyr".

The relation is the mean `mean_fe_h(age, branch)`: a line through the branch's anchor, at the slope the paper gives. The anchors (the thin
disc's at birth, the thick disc's at 8 Gyr), the scatters and the thick disc's share of the old stars are stand-ins (AR1 to AR4) chosen to
match the populations the generator has today (`tuning.METALLICITY_BY_POPULATION`: 0 for the young and intermediate disc, -0.325 for the old, whose ages reach the thick disc's 13 Gyr)
and the thin disc's lower limit near -0.7 (Bensby et al. 2014).

**How a star takes it.** The [Fe/H] is still drawn first, from the place and the population (the scatter's tables, the places' own
distributions and the bins of GEN.651 all stand on that order), and the age follows it: `tilted_age` draws `AGE_CANDIDATES` ages the way the
caller draws one and keeps one in proportion to the density of the relation at the star's [Fe/H] (`log_density`: the thin disc, and the
thick disc beside it for ages past 8 Gyr, mixed by `thick_share`). So a metal-poor star of the old disc is older than a metal-rich one,
and the population's [Fe/H] and the age window stay what they were. Only the disc's three populations take it (`takes`: each one's [Fe/H] distribution is about the relation's mix over its window, so the tilt moves the ages with the metals and leaves the population's mean age where it was, to 0.05 Gyr); the halo, the bulge, the
clusters and the satellites have their own histories. The relation is that of the solar neighbourhood, and the radial gradient is not taken
out of the star's [Fe/H], so a metal-poor place tilts its stars older, as an inside-out disc does.
"""

import math

from planetgen import tuning
from planetgen.physics import stellar_abundances

THIN_SLOPE_DEX_PER_GY = -0.025
"""float: d[Fe/H]/d(age) of the thin disc, dex per Gyr (Haywood et al. 2013: 0.025 dex per Gyr of metal gain, so older is poorer)."""
THICK_SLOPE_DEX_PER_GY = -0.15
"""float: d[Fe/H]/d(age) of the thick disc (Haywood et al. 2013: about 0.15 dex per Gyr)."""
THIN_BIRTH_FE_H = 0.05
"""float: The thin disc's mean [Fe/H] at age 0. [Stand-in AR1, tier T6, tweak L1: the young and intermediate disc of today average 0.]"""
THICK_ANCHOR_FE_H = -0.10
"""float: The thick disc's mean [Fe/H] at `tuning.STELLAR_THICK_START_GY`, where it ends. [Stand-in AR2, tier T6, tweak L1: Haywood's Fig. 9
runs from about -0.8 at 12 Gyr to about 0 at 8 Gyr.]"""
SIGMA = {stellar_abundances.BRANCH_THIN: 0.15, stellar_abundances.BRANCH_HIGH_ALPHA: 0.10}
"""dict: The scatter about the mean, dex, by branch. [Stand-in AR3, tier T6, tweak L1: the thin disc's is the generator's own 0.15; the
thick disc's is tighter, "a much tighter correlation".]"""
THICK_OLDEST_GY = 13.0
"""float: The oldest thick-disc age; the relation is held at its end beyond it."""


def mean_fe_h(age_gy, branch):
    """The mean [Fe/H] of the disc at this age on this branch."""
    if branch == stellar_abundances.BRANCH_THIN:
        return THIN_BIRTH_FE_H + THIN_SLOPE_DEX_PER_GY * age_gy
    if branch == stellar_abundances.BRANCH_HIGH_ALPHA:
        age = min(max(age_gy, tuning.STELLAR_THICK_START_GY), THICK_OLDEST_GY)
        return THICK_ANCHOR_FE_H + THICK_SLOPE_DEX_PER_GY * (age - tuning.STELLAR_THICK_START_GY)
    raise ValueError(f"unknown branch {branch!r}; expected one of {stellar_abundances.BRANCHES}")


THICK_SHARE = 0.8
"""float: The share of the disc's stars older than `tuning.STELLAR_THICK_START_GY` that are thick-disc, the rest being the thin disc's oldest.
[Stand-in AR4, tier T6, tweak L1: Haywood et al. 2013 have the thin disc start forming 8 to 9 Gyr ago, so little of it is older.]"""
AGE_CANDIDATES = 8
"""int: How many ages `tilted_age` draws to pick one from. A numerical resolution of the tilt, not a physical value: more candidates follow
the relation more closely."""


def takes(population):
    """Whether a star of this population follows the relation: one of the disc's (`tuning.DISK_AGE_LAW_POPULATIONS`). A star of no population
    (the whole disc, ages from 0 to 13 Gyr with the intermediate disc's [Fe/H]) does not: its [Fe/H] is not the mix of the relation over
    those ages, so a tilt would crowd it toward the young."""
    return population in tuning.DISK_AGE_LAW_POPULATIONS


def thick_share(age_gy):
    """The share of the disc's stars of this age that are thick-disc: none before `tuning.STELLAR_THICK_START_GY`, `THICK_SHARE` after."""
    return 0.0 if age_gy < tuning.STELLAR_THICK_START_GY else THICK_SHARE


def _log_normal(value, mean, sigma):
    return -0.5 * ((value - mean) / sigma) ** 2 - math.log(sigma)


def log_density(fe_h, age_gy):
    """The log of the density of [Fe/H] `fe_h` among the disc's stars of this age: the thin disc's normal about its mean, mixed with the thick
    disc's by `thick_share`. Up to a constant that does not depend on the age (the sqrt(2 pi) of the normals)."""
    thick = thick_share(age_gy)
    terms = [math.log(1.0 - thick) + _log_normal(fe_h, mean_fe_h(age_gy, stellar_abundances.BRANCH_THIN),
                                                  SIGMA[stellar_abundances.BRANCH_THIN])]
    if thick > 0.0:
        terms.append(math.log(thick) + _log_normal(fe_h, mean_fe_h(age_gy, stellar_abundances.BRANCH_HIGH_ALPHA),
                                                    SIGMA[stellar_abundances.BRANCH_HIGH_ALPHA]))
    peak = max(terms)
    return peak + math.log(sum(math.exp(term - peak) for term in terms))


def tilted_age(fe_h, population, draw_age, rng):
    """
    An age for a star of birth [Fe/H] `fe_h`, drawn as `draw_age()` draws one but kept in proportion to the relation's density at `fe_h`:
    `AGE_CANDIDATES` ages, one picked by a uniform number of `rng` with weights `exp(log_density(fe_h, age))`. A population that does not
    follow the relation (`takes`) gets one `draw_age()` and no more draws.
    """
    if not takes(population):
        return draw_age()
    ages = [draw_age() for _ in range(AGE_CANDIDATES)]
    logs = [log_density(fe_h, age) for age in ages]
    peak = max(logs)
    weights = [math.exp(value - peak) for value in logs]
    pick = rng.random() * sum(weights)
    for age, weight in zip(ages, weights):
        if pick < weight:
            return age
        pick -= weight
    return ages[-1]
