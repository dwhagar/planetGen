# planetgen/physics/debris_disk.py

"""
The debris disk of a main-sequence star (GEN.787): the cold belt of planetesimals and dust that 13 to 24 percent of the
Sun-like and A stars show in the infrared and 2 percent of the M stars do.

Whether a star has one is a draw at the observed fraction of its kind (Sibthorpe et al. 2018, MNRAS 475, 3046, the Herschel
DEBRIS survey: 47 discs round 275 FGK stars, F 23.9, G 14.3 and K 13.0 percent, M 2.2 percent of 89, 24 percent of the A
stars; Surveys find no clear dependence of the dust on the spectral type from late B to M, so a more massive
star takes the A figure). Where it is, the disk is a ring: its dust temperature is drawn round 58 K, the blackbody radius
the dust sits at is `r = (278.3 / T)^2 sqrt(L)` (Sibthorpe 2018, section 2.5, for a star of `L` solar luminosities), the
belt lies about 1.7 times further out than a blackbody grain would (their 7 to 40 AU blackbody radii "scale to about 40 AU
for realistic dust properties") and spans 0.7 to 1.5 of that radius. The descriptor holds no bodies and no mass: the mass
of a disk is not measured, only the light of its dust, so `disk_mass_solar` stays empty (the sheet rule says why).
"""

import math

from planetgen import tuning
from planetgen.util import draw

DEBRIS_KIND = "debris"
"""str: The `disk_kind` of a main-sequence star's cold belt."""

DEBRIS_CLASSES = ("V", "IV")
"""tuple: The Yerkes classes whose stars can hold a debris disk (dwarfs and subgiants; a giant has lost its belt to its wind and
its heat, and a white dwarf's dust is a different disk)."""


def incidence(mass_sol):
    """The share of stars of `mass_sol` solar masses that hold a debris disk (`tuning.DEBRIS_DISK_INCIDENCE`, one share for each band
    that ends at a mass of `tuning.DEBRIS_DISK_INCIDENCE_MASS_LIMITS_SOL`, the last band open)."""
    for upper, share in zip(tuning.DEBRIS_DISK_INCIDENCE_MASS_LIMITS_SOL, tuning.DEBRIS_DISK_INCIDENCE):
        if mass_sol < upper:
            return share
    return tuning.DEBRIS_DISK_INCIDENCE[-1]


def blackbody_radius_au(temperature_k, luminosity_sol):
    """The distance in AU at which a blackbody grain round a star of `luminosity_sol` is at `temperature_k`: 278.3 K at 1 AU
    falling as the square root of the distance and rising as the fourth root of the luminosity."""
    return (tuning.BLACKBODY_AT_ONE_AU_K / temperature_k) ** 2 * math.sqrt(luminosity_sol)


def draw_disk(mass_sol, luminosity_sol, rng=draw):
    """
    The debris disk of a main-sequence star, or None for one without (most of them): `disk_kind` "debris", `disk_inner_au`,
    `disk_outer_au`, and `disk_mass_solar` None. Takes two draws, one for whether the star has a disk and one for its dust
    temperature, in that order.
    """
    if rng.random() >= incidence(mass_sol):
        return None
    sigma, median = tuning.DEBRIS_DISK_TEMPERATURE_LOG_SIGMA, tuning.DEBRIS_DISK_TEMPERATURE_MEDIAN_K
    temperature = median * math.exp(rng.gauss(0.0, sigma))
    belt_au = tuning.DEBRIS_RADIUS_OVER_BLACKBODY * blackbody_radius_au(temperature, luminosity_sol)
    inner, outer = tuning.DEBRIS_DISK_EDGES_OVER_RADIUS
    return {"disk_kind": DEBRIS_KIND, "disk_mass_solar": None, "disk_inner_au": inner * belt_au, "disk_outer_au": outer * belt_au}
