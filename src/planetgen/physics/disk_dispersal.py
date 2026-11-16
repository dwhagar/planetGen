# planetgen/physics/disk_dispersal.py

"""
How long a star holds its disk (GEN.1551, slice 0: the law, not yet wired)
==========================================================================

A young star holds a class 0, I or II disk only until its own dispersal time; after that it is class III and bodies may stand where the disk was (GEN.1552). Today every star
loses its disk at one age, 3 Myr (`tuning.DISK_CLASS_AGES_YR`), which holds disks too long under 3 Myr and none between 3 and 10 Myr, where the surveys still see them.

**The law.** The share of young stars that still hold a primordial disk falls about exponentially with age, with a characteristic time of 2.5 Myr (half-life 1.7 Myr) [R5298], 2 to 3 Myr
at the short infrared wavelengths of 22 clusters from 1 to 100 Myr [R5296], and a half-life of 1.3 to 2 Myr (3.5 Myr with magnetic tracks) in 69 clusters of 5 Myr and under [R5297]. So a star's dispersal time is
`t = -tau ln(u)`, `u` a draw in (0, 1] from the star's own stream, and `tau` the characteristic time of its mass: [R5298] has it fall from about 3 Myr for brown dwarfs (under
0.08 solar masses) to about 1 Myr for stars over 1.3, so `tau` runs log-linear in mass through `DISPERSAL_TAU_MYR_BY_MASS` (the 0.5 solar mass knot is the survey's 2.5 Myr, where most of
the young stars of the surveys lie), flat outside.

**A close companion** (inside `disk_truncation.CLOSE_COMPANION_AXIS_AU`) disperses two thirds of the disks within 1 Myr [R5225]: that is an exponential of characteristic time `1 / ln 3` Myr, and it is the
same draw, `u`: the dispersal time is the shorter of the two (GEN.1475's `gas_survival_factor` reads this time, it is not a second draw).

**The draw is keyed to the star** (`draw_dispersal_time_yr`), never to its age, so the time does not re-roll when the age moves.

Not the whole of the survey curve: the shares a lane remembered (0.8 at 1 Myr, half at 3, 0.1 at 5 to 6) are not an exponential (this law gives 0.67, 0.30 and 0.11 for a star of 0.5 solar masses); the exponential is what the papers
say, and a lane that reads the cluster tables of [R5296] can replace it by a Weibull or a table without touching the callers. Disk fractions of the surveys also differ by half between wavelengths and
by the model of the pre-main-sequence tracks [R5297], so the 2.5 Myr is a middle, not a measurement of any star.
"""

import math

from planetgen.util import draw

DISPERSAL_TAU_MYR_BY_MASS = ((0.08, 3.0), (0.5, 2.5), (1.3, 1.0))
"""tuple: `(mass in solar masses, characteristic time in Myr)` of the exponential fall of the disk fraction: 3 Myr for brown dwarfs and 1 Myr over 1.3 solar masses, the ends [R5298] gives, and the
2.5 Myr of the survey as a whole at 0.5 solar masses, where the surveyed stars lie. The line between is a placement of the knots, nudged inside the spread of the three papers. [Stand-in 6.264 DZ1, tier T2, tweak L1.]"""

CLOSE_COMPANION_SHARE_WITHIN_1_MYR = 2.0 / 3.0
"""float: The share of the stars with a companion inside `disk_truncation.CLOSE_COMPANION_AXIS_AU` whose disk is gone within 1 Myr ([R5225], read in `disk_truncation`). [Stand-in 6.264 DZ3, tier T2, tweak L0.]"""

CLOSE_COMPANION_TAU_MYR = 1.0 / math.log(1.0 / (1.0 - CLOSE_COMPANION_SHARE_WITHIN_1_MYR))
"""float: The characteristic time of the exponential that loses `CLOSE_COMPANION_SHARE_WITHIN_1_MYR` of its disks within 1 Myr, `1 / ln 3`, 0.91 Myr. [Stand-in 6.264 DZ3.]"""

YEAR_PER_MYR = 1.0e6


def characteristic_time_yr(mass_sol, close_companion=False):
    """The characteristic time `tau` of the disk fraction of a star of this mass, years: `DISPERSAL_TAU_MYR_BY_MASS` (log-linear in mass, flat outside the knots); with a close companion the shorter of it and `CLOSE_COMPANION_TAU_MYR`."""
    if mass_sol <= 0.0:
        raise ValueError(f"a disk time needs a mass above zero, not {mass_sol}")
    knots = DISPERSAL_TAU_MYR_BY_MASS
    if mass_sol <= knots[0][0]:
        tau = knots[0][1]
    elif mass_sol >= knots[-1][0]:
        tau = knots[-1][1]
    else:
        tau = knots[-1][1]
        for (m0, t0), (m1, t1) in zip(knots, knots[1:]):
            if mass_sol <= m1:
                share = math.log(mass_sol / m0) / math.log(m1 / m0)
                tau = t0 + share * (t1 - t0)
                break
    if close_companion:
        tau = min(tau, CLOSE_COMPANION_TAU_MYR)
    return tau * YEAR_PER_MYR


def dispersal_time_yr(mass_sol, u, close_companion=False):
    """
    The age at which a star of this mass loses its disk, years, from the draw `u` in (0, 1]: `-tau ln u`. A larger `u` is an earlier loss: `u = 1` is 0 and `u` near 0 a long-lived disk.

    A close companion changes `tau`, not the draw, so the same star has the shorter time (never the longer) with a companion than without.
    """
    if not 0.0 < u <= 1.0:
        raise ValueError(f"the draw must be in (0, 1], not {u}")
    return -characteristic_time_yr(mass_sol, close_companion) * math.log(u)


def disk_fraction(mass_sol, age_yr, close_companion=False):
    """The share of stars of this mass still holding a disk at this age: `exp(-age / tau)` (the share of draws whose `dispersal_time_yr` is later than the age)."""
    return math.exp(-max(age_yr, 0.0) / characteristic_time_yr(mass_sol, close_companion))


def draw_dispersal_time_yr(star_key, mass_sol, close_companion=False):
    """
    A star's own dispersal time, years: the draw `u` is read from the stream of `star_key` (anything that names the star and does not move with its age: its stored seed or its place) so the same
    star always has the same time, with the age, the mass and the companion moving only what is done with `u`.
    """
    u = 1.0 - draw.Stream(f"disk-dispersal:{star_key}").random()        # (0, 1]
    return dispersal_time_yr(mass_sol, u, close_companion)
