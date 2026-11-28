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
(200 AU for an A star) and gave far too little dust. [Stand-in 6.53 DU4, tier T6, tweak L1: 32 rings.]

Detection (GEN.1884, GEN.1892)
------------------------------
A survey finds the dust by its flux, not by its share of the star's light. DEBRIS observed every star at 100 and 160 micrometres to a flux its instrument could tell from the star's noise, and its sample is volume limited
by spectral type: the nearest 90 or so stars of each, the largest distance of the F, G and K stars being 24, 21 and 16 pc ([R4081] section 2.1 and 2.2). A star that shines more is therefore seen at a larger distance, and
the dust a survey can find around it is fainter as a share of its light: the dust's luminosity `f L` over the square of the sample's distance must pass one flux, `DETECTABLE_DUST_FLUX`. The distance of the sample goes with
the star's luminosity as the three published distances do (`survey_distance_pc`, fitted to [R4081] section 2.2: the power 0.164 through the F and K stars, the G star 3 percent off). The first cut of this module (GEN.1278)
used the star's own light as the yardstick, the share of its flux at a long wavelength, `f T^3`; it gave A, F, G and K stars shares that needed a cut rising 4 times from A to K, which a flux limit at a distance that falls with the
star's light explains (GEN.1892: dust temperature and wavelength, tried first, did not). `DETECTABLE_DUST_FLUX` is the one dial, tuned to the surveys' shares by type [R4081] (L3, GEN.1884: the belts since GEN.1797
sit nearer and are heavier, so the cut of GEN.1278 held 40 percent of them). [Stand-in 6.53 DU1, tier T9, tweak L3; the distance law DU5, tier T3, tweak L1.]
"""

import math

from planetgen.physics import belt_mass, constants
from planetgen.util.checks import finite_domain

DEBRIS_KIND = "debris"
"""str: The `disk_kind` of a main-sequence star's cold belt."""

DEBRIS_CLASSES = ("V", "IV")
"""tuple: The Yerkes classes whose stars can hold a debris disk (dwarfs and subgiants; a giant has lost its belt to its wind and its heat, and a white dwarf's dust is a different disk)."""

F_MAX_COEFFICIENT = 0.004
"""float: The coefficient of the largest fractional luminosity of a belt ([R4785: Eq. 18], `r` in AU, `Dc` in km, `t` in Myr). [Stand-in 6.53 DU2, tier T1, tweak L0.]"""

CASCADE_MAX_KM = 2000.0
"""float: The largest body the dust cascade is reckoned for, km: the nominal `Dc` of [R4785] (their Fig. 2 and Eq. 21). [Stand-in 6.53 DU3, tier T1, tweak L1.]"""

DETECTABLE_DUST_FLUX = 1.9e-7
"""float: The least dust a survey finds, as the dust's luminosity in solar units over the square of the distance of the sample in parsecs (`f L / d^2`): the flux its instrument tells from the star, one for every spectral type.
Tuned on 100 stars each of A1V, F5V, G2V and K3V so that the F and G shares sit in the survey's bands and the K and A shares come nearest to theirs (A 13, F 28, G 17, K 7 percent against 20 to 33, 19 to 29, 10 to 19 and 9 to 17, [R4081]); a cut on the mass of the generator's heavy old outer disks, not the instrument's own limit (which is
not given as a fractional luminosity in the paper, only as completeness contours; a fractional luminosity near 1e-5 is where its detected disks lie). [Stand-in 6.53 DU1, tier T9, tweak L3.]"""

SURVEY_DISTANCE_PC = 21.0
"""float: The distance of the survey's sample of a star of one solar luminosity, pc: the farthest G star of the Herschel DEBRIS sample is 21 pc away ([R4081] section 2.2: F, G and K stars at 24, 21 and 16 pc).
[Stand-in 6.53 DU5, tier T3, tweak L1.]"""

SURVEY_DISTANCE_LUMINOSITY_EXPONENT = 0.164
"""float: How the distance of the sample goes with the star's luminosity, the power through the F and K stars of the survey (24 and 16 pc, at the generator's F5V and K3V luminosities of 3.3 and 0.28 solar), which puts the G
star 3 percent from its 21 pc; the volume-limited sample reaches farther for stars that are rarer and shine more. For an A star (about 30 solar luminosities) it gives 36 pc; the survey says only that its A stars are the farthest.
[Stand-in 6.53 DU5, tier T3, tweak L1.]"""

ANNULI = 32
"""int: The rings of equal log width a belt is cut into for its dust (`belt_fractional_luminosity`). [Stand-in 6.53 DU4, tier T6, tweak L1.]"""


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


def survey_distance_pc(luminosity_sol):
    """The distance, pc, to which a survey's volume-limited sample of stars of `luminosity_sol` solar luminosities reaches (`SURVEY_DISTANCE_PC`, `SURVEY_DISTANCE_LUMINOSITY_EXPONENT`)."""
    if not luminosity_sol > 0.0:
        raise ValueError("the star needs a positive luminosity")
    return SURVEY_DISTANCE_PC * luminosity_sol ** SURVEY_DISTANCE_LUMINOSITY_EXPONENT


def detectable_fractional_luminosity(luminosity_sol):
    """The fractional luminosity a survey can find around a star of `luminosity_sol` solar luminosities: the flux `DETECTABLE_DUST_FLUX` at the distance of its sample, over the star's light."""
    return DETECTABLE_DUST_FLUX * survey_distance_pc(luminosity_sol) ** 2 / luminosity_sol


def is_detectable(fractional, luminosity_sol):
    """Whether a disk of fractional luminosity `fractional` around a star of `luminosity_sol` solar luminosities is one a survey finds."""
    return fractional >= detectable_fractional_luminosity(luminosity_sol)
