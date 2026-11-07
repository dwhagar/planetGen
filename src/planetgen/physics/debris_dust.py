# planetgen/physics/debris_dust.py

"""
The dust a cold belt makes, and whether a survey would find it (GEN.1071)
==========================================================================

A debris disk is the dust that the collisions of a belt's planetesimals grind off, and a survey finds the belt by that dust's light. The outer disk of a system is a belt like the others
(`belt_origin.outer_disk`), so whether the star "has a debris disk" is not a second roll: it is whether the belt's present dust, `f` times the star's light, is above what a survey can tell
from the star.

The law is that of [R4785]. A belt that has ground down as far as the collisions allow (`M = A / t`, `belt_mass.collisional_constant`) gives, in the size law `q = 11/6` and with the
planetesimals' eccentricity equal to their inclination, the largest fractional luminosity `f_max = 0.004 r^1.5 (dr/r) Dc^0.5 L^-0.5 / (t G(Xc))` (Eq. 18 there, `r` in AU, `Dc` in km, `t` in Myr,
`L` in solar luminosities, the same `G` and `Xc` as the mass). A belt that has not yet ground that far has the dust in proportion to its mass (`K`, the cross-section and so the light, goes with
`M_tot`, their Eq. 14), so `f = f_max M / M_max`. A belt whose largest bodies are too strong for the collisions to start the cascade makes no dust (not reached once `Dc` is capped, below, except by a belt so wide and eccentric that `Xc` passes 1).

The cascade's largest body
--------------------------
The law takes `Dc` as the largest planetesimal of the cascade, 2000 km at most in the paper's nominal belt. A belt's own largest body (`belt_ledger`) can be thousands of kilometres across in a heavy
outer disk: a body that size is a protoplanet that grows alone and does not sit in the cascade, and it would put the cascade past `Xc = 1` where the law gives no dust at all. The dust is
therefore reckoned for `min(Dc, CASCADE_MAX_KM)`; a belt heavier than that cascade's ceiling `M_max` is in its steady state and shines at `f_max`.

Detection
---------
`DETECTABLE_FRACTIONAL_LUMINOSITY` is the dust a survey finds against the star's light, a stand-in. [R4081] find most of their disks near a fractional luminosity of 1e-5, but the belts the
generator builds are all heavy and old and nearly all at their steady state, so about half the K, G and F stars would shine above 1e-5 where the survey finds 13 to 24 percent. The cut is
tuned to the combined share of the K, G and F dwarfs, 17 percent [R4081]: it is the number the generator's disks must pass to be as rare as the surveys' (L3, hand-tuned against a sweep of
the generator). It leaves the A stars (young, bright, their dust diluted by `L^-0.5`) with fewer disks than the survey's 24 percent and the F stars with about half theirs; both are gaps of
the outer disk's mass and width, not of this law. [Stand-in DU1, tier Marginal educated guess from needed output, tweak L3.]
"""

import math

from planetgen.physics import belt_mass
from planetgen.util.checks import finite_domain

DEBRIS_KIND = "debris"
"""str: The `disk_kind` of a main-sequence star's cold belt."""

DEBRIS_CLASSES = ("V", "IV")
"""tuple: The Yerkes classes whose stars can hold a debris disk (dwarfs and subgiants; a giant has lost its belt to its wind and its heat, and a white dwarf's dust is a different disk)."""

F_MAX_COEFFICIENT = 0.004
"""float: The coefficient of the largest fractional luminosity of a belt ([R4785: Eq. 18], `r` in AU, `Dc` in km, `t` in Myr). [Stand-in DU2, tier T1, tweak L0.]"""

CASCADE_MAX_KM = 2000.0
"""float: The largest body the dust cascade is reckoned for, km: the nominal `Dc` of [R4785] (their Fig. 2 and Eq. 21). [Stand-in DU3, tier T1, tweak L1.]"""

DETECTABLE_FRACTIONAL_LUMINOSITY = 3.0e-5
"""float: The fractional luminosity, dust over star, at which the generator's disks are as rare as the surveys' (the real survey limit is near 1e-5). [Stand-in DU1, tier Marginal educated guess from needed output, tweak L3.]"""


@finite_domain()
def fractional_luminosity(lower_au, upper_au, luminosity_sol, age_myr, largest_km, mass_earth, star_mass_solar):
    """
    The dust of a belt between two distances, as a share of its star's light: `f_max M / M_max` of the module. Zero for a belt whose largest bodies the collisions cannot disrupt
    (`belt_mass.collisional_constant` is not positive there) or one with no mass.

    Args:
        lower_au, upper_au (float): The belt's edges.
        luminosity_sol (float): The star's luminosity, solar units (the age it has now).
        age_myr (float): The belt's age, Myr.
        largest_km (float): The belt's largest planetesimal (`Dc`).
        mass_earth (float): The belt's mass now, Earth masses.
        star_mass_solar (float): The star's mass at the belt's formation, solar masses.
    """
    if not (0.0 < lower_au < upper_au and luminosity_sol > 0.0 and age_myr > 0.0 and largest_km > 0.0):
        raise ValueError("the belt, the light, the age and the largest body must all be positive")
    if not mass_earth > 0.0:
        return 0.0
    largest_km = min(largest_km, CASCADE_MAX_KM)
    constant = belt_mass.collisional_constant(lower_au, upper_au, star_mass_solar, largest_km)
    if not constant > 0.0:
        return 0.0
    r = 0.5 * (lower_au + upper_au)
    e = (upper_au - lower_au) / (upper_au + lower_au)
    xc = 1.0e-3 * (r * belt_mass.disruption_threshold_j_per_kg(largest_km) / e ** 2) ** (1.0 / 3.0)
    steady = F_MAX_COEFFICIENT * r ** 1.5 * ((upper_au - lower_au) / r) * math.sqrt(largest_km) / math.sqrt(luminosity_sol) / (age_myr * belt_mass._g(xc))
    return steady * min(1.0, mass_earth * age_myr / constant)


def is_detectable(fractional):
    """Whether a disk of fractional luminosity `fractional` is one a survey finds."""
    return fractional >= DETECTABLE_FRACTIONAL_LUMINOSITY
