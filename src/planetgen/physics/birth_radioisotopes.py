# planetgen/physics/birth_radioisotopes.py

"""
Aluminium-26 and iron-60 at a star's birth, from the gamma-ray masses, the star formation of its region and the granularity of the
sources (GEN.476, part 1; docs/design/composition/radioisotope-tuning.md section 5).

The radiogenic module drew the aluminium-26 ratio log-uniform from 0.1 to 10 times the Solar System's 5.23e-5 (`radiogenic`, a stand-in
taken from the Solar System's own value and range). Here the ratio is built from what makes it:

1. **The yield per mass of stars formed.** The gamma-ray lines of 26Al (1809 keV) and 60Fe (1173 and 1332 keV) measure the mass of each
   in the Galaxy: 1.7 to 3 solar masses of 26Al, best 2.0 (Diehl et al. 2021, review, arXiv:2109.08558, section 2.4; Limongi and Chieffi
   2006 give 2.0), and 60Fe from the line ratio 0.15 (flux, 60Fe over 26Al; the same review). The stars of the Galaxy form at 1.65 +- 0.19
   solar masses a year (Licquia and Newman 2015, arXiv:1407.1078), and each isotope is in steady state: its mass is the yield times the
   rate times its mean life. So the yield per solar mass of stars formed is `mass / (rate * mean_life)`, about 1.2e-6 for 26Al and
   4e-7 for 60Fe. That is a property of massive-star nucleosynthesis, not of the Milky Way as a place; the Milky Way is only where it
   was measured, and the metallicity dependence of the yield (stronger Wolf-Rayet winds at high Z) is left out.
2. **The mean ratio in a region.** The 26Al mass fraction of the gas is `yield * mean_life / depletion_time` (the gas mass over the star
   formation rate of the region), and the ratio to stable 27Al is that over the aluminium mass fraction (27/26 for number). For the
   Galaxy's 0.7e10 solar masses of cold gas (Bland-Hawthorn and Gerhard 2016, arXiv:1602.07702, the value they adopt) the depletion
   time is 4.2 Gyr and the mean 26Al/27Al comes to about 5e-6, ten times under the Solar System's 5.23e-5 (Jacobsen et al. 2008, via
   Lugaro et al. 2018, arXiv:1808.00233, Table 1); the Solar System is a test point and never an input.
3. **Granularity.** The sources are discrete: Lugaro et al. 2018 (section 3.4, Eq. 12, READ) show that when the recurrence time delta
   between additions is longer than the mean life tau the abundance at a random moment is dominated by the last addition and its decay,
   not by the steady-state mean. `sample_ratio` draws the events backward from the star's birth (Poisson in time with recurrence
   `delta`, 5 to 10 Myr for core-collapse events, from Lugaro et al.), gives each a yield that keeps the
   mean (`mean * delta / tau`) times a log-normal factor for the richness of the cluster that made it, and sums the decays. 26Al has
   tau 1.035 Myr, so it is spiky (the relative scatter is about 1.6 to 2.2 before the cluster factor); 60Fe has 3.78 Myr and
   averages over more events.

Limits. The yield does not vary with metallicity or with the shape of the initial mass function; the recurrence is one number for a
whole region (it should follow the region's cluster formation rate); the cluster factor `EVENT_YIELD_SIGMA_LN` is tuned (tweak L2) so that
the Solar System ratio (about ten times the mean) lies inside the distribution's upper tail, which the equal-events model cannot reach
(its maximum is about `delta/tau` times the mean, 5 to 10 times). Nothing in generation calls this yet: the wiring step is GEN.476's second
step (`radiogenic.draw_content` taking the ratio from here, with the host's region and formation time) on a merge train.
"""

import math

from planetgen.util import draw
from planetgen.util.checks import finite_domain

MEAN_LIFE_YR = {"Al26": 1.035e6, "Fe60": 3.78e6}
"""dict: The mean life of each isotope, years (Lugaro et al. 2018 Table 1: 26Al half-life 0.717 Myr, mean life 1.035; 60Fe 2.62 and 3.78).
[Stand-in BR1, tier T1, tweak L0.]"""

ATOMIC_MASS = {"Al26": 26.0, "Al27": 27.0, "Fe60": 60.0, "Fe56": 56.0}
"""dict: Mass numbers used to turn a mass ratio into a number ratio."""

GALAXY_AL26_MASS_SOLAR = 2.0
"""float: The mass of 26Al in the Galaxy, solar masses (Diehl et al. 2021 section 2.4: 1.7 to 3 from the 1809 keV line, with 30 percent from the
assumed source distribution; Limongi and Chieffi 2006: 2.0). [Stand-in BR2, tier T1, tweak L0.]"""

FE60_AL26_FLUX_RATIO = 0.15
"""float: The gamma-ray line flux ratio 60Fe over 26Al of the Galaxy (Diehl et al. 2021, 60Fe section; the ratio of the yields by number once
converted). [Stand-in BR2, tier T1, tweak L0.]"""

GALAXY_STAR_FORMATION_SOLAR_MASSES_PER_YR = 1.65
"""float: The Galaxy's star formation rate, solar masses a year (Licquia and Newman 2015: 1.65 +- 0.19). [Stand-in BR3, tier T1, tweak L0.]"""

RECURRENCE_RANGE_YR = (5.0e6, 1.0e7)
"""tuple: The recurrence time between core-collapse additions to a region's gas, years (Lugaro et al. 2018, section 3.4: 5 to 10 Myr; 2ND-HAND). A star draws one value log-uniformly. [Stand-in BR4, tier T2, tweak L0.]"""

EVENT_YIELD_SIGMA_LN = 1.0
"""float: The log-normal spread of the yield of one addition about the mean yield, the richness of the cluster that made it (tuned so the
Solar System ratio lies in the upper tail of the distribution; the mean of the factor is held at 1). [Stand-in BR5, tier T6, tweak L2.]"""

WINDOW_MEAN_LIVES = 12.0
"""float: How far back the events are summed, in mean lives (a contribution beyond is under 1e-5 of its first value)."""


@finite_domain()
def yield_per_solar_mass_formed(isotope):
    """The mass of an isotope in steady state per solar mass of stars formed per unit mean life: the Galaxy's mass of it over its star
    formation rate and its mean life (dimensionless, solar masses of isotope per solar mass of stars formed)."""
    if isotope == "Al26":
        mass = GALAXY_AL26_MASS_SOLAR
    elif isotope == "Fe60":
        # the flux ratio is the ratio of the decay rates; the masses follow with the mean lives and the mass numbers
        mass = (GALAXY_AL26_MASS_SOLAR * FE60_AL26_FLUX_RATIO * MEAN_LIFE_YR["Fe60"] / MEAN_LIFE_YR["Al26"]
                * ATOMIC_MASS["Fe60"] / ATOMIC_MASS["Al26"])
    else:
        raise ValueError(f"no gamma-ray mass for {isotope!r}")
    return mass / (GALAXY_STAR_FORMATION_SOLAR_MASSES_PER_YR * MEAN_LIFE_YR[isotope])


@finite_domain()
def mean_ratio(isotope, depletion_time_yr, stable_mass_fraction):
    """
    The steady-state number ratio of the radioactive isotope to its stable partner (26Al/27Al, 60Fe/56Fe) in the gas of a region.

    Args:
        isotope (str): "Al26" or "Fe60".
        depletion_time_yr (float): The region's gas mass over its star formation rate, years.
        stable_mass_fraction (float): The mass fraction of the stable partner in the gas (27Al or 56Fe).
    """
    if not (depletion_time_yr > 0.0 and stable_mass_fraction > 0.0):
        raise ValueError("the depletion time and the stable mass fraction must be positive")
    mass_ratio = yield_per_solar_mass_formed(isotope) * MEAN_LIFE_YR[isotope] / depletion_time_yr / stable_mass_fraction
    partner = "Al27" if isotope == "Al26" else "Fe56"
    return mass_ratio * ATOMIC_MASS[partner] / ATOMIC_MASS[isotope]


def sample_ratio(isotope, mean, stream, recurrence_yr=None):
    """
    A star's ratio at birth: the sum of the decayed additions of a Poisson sequence of events with the given recurrence, each carrying the
    mean yield `mean * recurrence / mean_life` times a log-normal cluster factor of mean 1. The recurrence is drawn log-uniformly from
    `RECURRENCE_RANGE_YR` when not given. Deterministic for one stream.
    """
    if not mean >= 0.0:
        raise ValueError("the mean ratio must be zero or more")
    tau = MEAN_LIFE_YR[isotope]
    if recurrence_yr is None:
        low, high = RECURRENCE_RANGE_YR
        recurrence_yr = math.exp(stream.uniform(math.log(low), math.log(high)))
    each = mean * recurrence_yr / tau
    sigma = EVENT_YIELD_SIGMA_LN
    window = WINDOW_MEAN_LIVES * tau
    total, age = 0.0, 0.0
    while True:
        age += -recurrence_yr * math.log(1.0 - stream.random())
        if age > window:
            return total
        factor = math.exp(sigma * stream.gauss() - 0.5 * sigma * sigma)
        total += each * factor * math.exp(-age / tau)


def fraction_above(isotope, mean, threshold, trials=4000, seed="birth-radioisotopes"):
    """The share of stars whose birth ratio exceeds `threshold`, by `trials` draws on a fixed stream (a test and a note use it, not the generator)."""
    stream = draw.Stream(seed)
    return sum(sample_ratio(isotope, mean, stream) > threshold for _ in range(trials)) / trials
