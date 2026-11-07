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

The belt's width (GEN.1278)
---------------------------
The law above is for a ring. A belt that spans the outer disk (11 to 400 AU around an A star) is not one: its dust is the sum of the rings in it, and the rings do not weigh the same. The solids fall as the
minimum-mass nebula does (`planet_mass.solids_surface_density_gcm2`), so the inner rings hold the mass and, with the dust going as the mass over the radius squared (`f ~ M r^-2 L^-1/2 Dc^-1/2` when the ring has
not ground down), give most of the light. `belt_fractional_luminosity` therefore cuts the belt into `ANNULI` rings of equal width in log radius, shares the belt's present mass among them by the solids they held, and adds their
dust, each ring keeping the belt's eccentricity (the stirring is the belt's: a thin ring has no eccentricity of its own, and `Xc` would pass 1 and the cascade never start). The old rule put all the mass at the middle radius of the belt
(200 AU for an A star) and gave far too little dust. [Stand-in DU4, tier Confident educated guess from physics, tweak L1: 32 rings.]

Detection
---------
A survey finds the dust as an excess over the star at a long wavelength, and what it can tell from the star is a fixed share of the star's own flux there. The star's flux in the Rayleigh-Jeans tail goes as
`R^2 T`, its luminosity as `R^2 T^4`, so the dust (a share `f` of the luminosity) stands out of the photosphere as `f T^3`: the fractional luminosity a survey can find falls as the cube of the star's temperature
(`detectable_fractional_luminosity`). A hot A star shows a disk of a fraction the Sun's would hide. `DETECTABLE_FRACTIONAL_LUMINOSITY` is the limit at the Sun's temperature, a stand-in: [R4081] find most of their disks
near a fractional luminosity of 1e-5 but the belts the generator builds are all heavy and old, so the limit is tuned to the surveys' shares by mass, 14 percent of G and 13 percent of K dwarfs, 24 percent of F and A stars
[R4081] (L3, hand-tuned against a sweep of 100 stars each of A1V, F5V, G2V and K3V: 17, 20, 17 and 18 percent). [Stand-in DU1, tier Marginal educated guess from needed output, tweak L3.]
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

DETECTABLE_FRACTIONAL_LUMINOSITY = 5.0e-5
"""float: The fractional luminosity, dust over star, a survey finds around a star of `REFERENCE_TEMPERATURE_K`, at which the generator's disks are as common as the surveys' (the real survey limit is near 1e-5).
[Stand-in DU1, tier Marginal educated guess from needed output, tweak L3.]"""

REFERENCE_TEMPERATURE_K = 5772.0
"""float: The Sun's effective temperature, K: the star `DETECTABLE_FRACTIONAL_LUMINOSITY` is the limit for. [Stand-in DU5, tier T1, tweak L0.]"""

DETECTION_TEMPERATURE_EXPONENT = 3.0
"""float: How the detectable fractional luminosity falls with the star's temperature: the dust against the photosphere at a long wavelength goes as `f T^3` (Rayleigh-Jeans flux `R^2 T`, luminosity `R^2 T^4`).
[Stand-in DU5, tier Confident educated guess from physics, tweak L1.]"""

ANNULI = 32
"""int: The rings of equal log width a belt is cut into for its dust (`belt_fractional_luminosity`). [Stand-in DU4, tier Confident educated guess from physics, tweak L1.]"""


@finite_domain()
def fractional_luminosity(lower_au, upper_au, luminosity_sol, age_myr, largest_km, mass_earth, star_mass_solar, eccentricity=None):
    """
    The dust of a ring between two distances, as a share of its star's light: `f_max M / M_max` of the module. Zero for a ring whose largest bodies the collisions cannot disrupt
    (`belt_mass.collisional_constant` is not positive there) or one with no mass. The planetesimals' eccentricity is the ring's half-width over its radius unless `eccentricity` is given.

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
    constant = belt_mass.collisional_constant(lower_au, upper_au, star_mass_solar, largest_km, eccentricity=eccentricity)
    if not constant > 0.0:
        return 0.0
    r = 0.5 * (lower_au + upper_au)
    e = eccentricity if eccentricity is not None else (upper_au - lower_au) / (upper_au + lower_au)
    xc = 1.0e-3 * (r * belt_mass.disruption_threshold_j_per_kg(largest_km) / e ** 2) ** (1.0 / 3.0)
    steady = F_MAX_COEFFICIENT * r ** 1.5 * ((upper_au - lower_au) / r) * math.sqrt(largest_km) / math.sqrt(luminosity_sol) / (age_myr * belt_mass._g(xc))
    return steady * min(1.0, mass_earth * age_myr / constant)


@finite_domain()
def belt_fractional_luminosity(lower_au, upper_au, luminosity_sol, age_myr, largest_km, mass_earth, star_mass_solar, surface_density):
    """
    The dust of a belt as the sum of the dust of the rings in it (GEN.1278, "The belt's width" in the module note): `ANNULI` rings of equal width in log radius share `mass_earth` by
    `surface_density(r) * area` and each is `fractional_luminosity` with the belt's own eccentricity.

    Args:
        lower_au, upper_au (float): The belt's edges.
        luminosity_sol, age_myr, largest_km, mass_earth, star_mass_solar: As `fractional_luminosity` for the whole belt.
        surface_density (callable): The solids' surface density at a distance in AU (any unit: only the ratios count).
    """
    if not (0.0 < lower_au < upper_au):
        raise ValueError("the belt needs 0 < lower < upper")
    if not mass_earth > 0.0:
        return 0.0
    eccentricity = (upper_au - lower_au) / (upper_au + lower_au)
    step = (math.log(upper_au) - math.log(lower_au)) / ANNULI
    edges = [lower_au * math.exp(step * i) for i in range(ANNULI + 1)]
    weights = [surface_density(math.sqrt(a * b)) * (b * b - a * a) for a, b in zip(edges, edges[1:])]
    total = sum(weights)
    if not total > 0.0:
        return 0.0
    return sum(fractional_luminosity(a, b, luminosity_sol, age_myr, largest_km, mass_earth * weight / total, star_mass_solar, eccentricity)
               for (a, b), weight in zip(zip(edges, edges[1:]), weights) if weight > 0.0)


def detectable_fractional_luminosity(temperature_k):
    """The fractional luminosity a survey can find around a star of `temperature_k`: `DETECTABLE_FRACTIONAL_LUMINOSITY` at the Sun's temperature, falling as the cube of the temperature."""
    if not temperature_k > 0.0:
        raise ValueError("the star needs a positive temperature")
    return DETECTABLE_FRACTIONAL_LUMINOSITY * (REFERENCE_TEMPERATURE_K / temperature_k) ** DETECTION_TEMPERATURE_EXPONENT


def is_detectable(fractional, temperature_k):
    """Whether a disk of fractional luminosity `fractional` around a star of `temperature_k` is one a survey finds."""
    return fractional >= detectable_fractional_luminosity(temperature_k)
