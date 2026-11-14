# planetgen/physics/influx.py

"""
The dust and comet influx a planet takes from every belt and comet of its system (GEN.1440, GEN.1453 slice A)
==============================================================================================================

Until now the only reservoir a planet's impacts knew was the late-accretion number of one Solar-System-sized leftover (`impact_history`, GEN.739). A planet's real supply is what the belts and comets of its
own system send in, less what a giant between them takes. This module works that out as a pure function of stored data: the reservoirs (belts with their mass and ring, comets with their mass and
perihelion), the giants (mass, radius, distance), the star's age and the planet's radius and distance. Nothing is stored and nothing is drawn; a sheet works it out when it is made (the dust is a pure
function of the stored belts and planets, so a column would only be a second copy of them that a rescan could leave stale).

**The law (a stand-in throughout; rows IF1 to IF6 of the gap plan).**

    influx = sum over reservoirs j of  supply_j  x  reach_j(planet)  x  cross_section_share(planet)  x  shield_j(planet)

* *Supply.* A belt loses mass to collisions as `M / t` at its age `t` (the collisional cascade `M = A / t` of `belt_mass`, so `|dM/dt| = M / t`); the share `BELT_DUST_SHARE` of that loss is dust or
  small bodies fine enough to drift inward past the belt's edge (Poynting-Robertson drag and resonance leaks). A comet reservoir sends its nucleus in over `COMET_LIFETIME_YR`.
* *Reach.* Only what drifts inward arrives: a planet inside the belt's outer edge `upper` is reached, one outside it is not (nothing is sent outward). A planet inside the ring itself gets the share of
  the ring that lies outside it (the ring is shared by the part of its width the planet's orbit leaves outside). A comet reaches a planet whose orbit its perihelion passes inside.
* *Cross-section.* An isotropic stream through the sphere of the planet's orbit hits the planet with the share `R_eff^2 / (4 a^2)` of it, `R_eff` the planet's radius (gravitational focusing is not counted: the
  escape speed of a planet against the speed of the stream is a second-order change for the dust and not known for the comets).
* *Crossings (GEN.1543).* A grain of belt dust does not cross the planet's orbit once. Poynting-Robertson drag takes a grain inward at `da/dt = 2 beta G M / (c a)` (`beta` the light-to-gravity ratio of
  the grains, `GRAIN_BETA` for the Sun), and a grain of eccentricity `GRAIN_ECCENTRICITY` crosses the orbit radius twice a period for as long as its semi-major axis is within that share of the planet's: `2 e a / (da/dt)`.
  The number of crossings is `N = e c sqrt(a) / (pi beta sqrt(G M))` (about 3e4 at 1 AU), and the planet takes `1 - exp(-N R_eff^2 / (4 a^2))` of the grain's passes. Without it the Earth got 45 kg a year
  from a Solar-System belt against the 1e7 measured (Plane 2012, Chem. Soc. Rev. 41, 6492, recalled); with the asteroid belt alone about 8e5 (the Jupiter-family comets are the larger share of the real dust). A comet is a nucleus that crosses once and is not multiplied.
* *Shield.* Each giant that lies between the reservoir and the planet ejects a share of what drifts past it, `1 - exp(-Theta / THETA_EJECT)` with `Theta` its Safronov number (`scatter`, the rule the
  swept disk of the giants follows); what is left is the product over the giants in the way.

What it leaves out (stated; rows IF5 and IF6): the giants' scattering of comets inward (a giant that throws a comet in as well as out), a companion star's stirring of a reservoir (GEN.1492), the
outward flux of a belt inside the planet's orbit, gravitational focusing, and the dust a planet's own moons and rings make.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, scatter

BELT_DUST_SHARE = 0.1
"""float: The share of a belt's collisional mass loss that reaches the inner system as dust and small bodies. A stand-in, in the range of the share of a cascade's mass that is in grains small enough for
radiation forces (the rest is ground to dust that leaves or stays). [Stand-in IF1, tier Marginal educated guess from physics, tweak L4.]"""

COMET_LIFETIME_YR = 1.0e7
"""float: The time over which a comet nucleus of a reservoir delivers its mass to the planets (its dynamical lifetime before ejection or sublimation, in the tens of millions of years for the
short-period and Halley kinds). [Stand-in IF2, tier Marginal educated guess from physics, tweak L4.]"""

GRAIN_BETA = 0.01
"""float: The ratio of radiation force to gravity for the grains that carry a belt's dust to the inner system, round the Sun (the zodiacal dust is mostly grains of tens to hundreds of microns, beta about 0.002 to
0.05); for another star it scales as its luminosity over its mass. [Stand-in IP1, tier Marginal educated guess from physics, tweak L2.]"""

GRAIN_ECCENTRICITY = 0.1
"""float: The orbital eccentricity of the dust grains, that of the parent belt's bodies stirred by its planets (the asteroid belt's proper eccentricities are 0.1 to 0.2). It sets how many times a grain crosses
a planet's orbit while it drifts. [Stand-in IP2, tier Marginal educated guess from physics, tweak L2.]"""

SOLAR_BETA_RATIO = 1.0
"""float: `GRAIN_BETA` is for a star with a luminosity-to-mass ratio of the Sun's; `grain_beta` multiplies it by this ratio of the star's."""


@dataclass(frozen=True)
class Reservoir:
    """A belt or a comet of the system: `kind` is `belt` or `comet`; `mass_kg`; the ring `lower_au` to `upper_au` (a comet: both its perihelion); `age_yr` is the belt's age."""

    kind: str
    mass_kg: float
    lower_au: float
    upper_au: float
    age_yr: float = 0.0


@dataclass(frozen=True)
class Giant:
    """A giant planet in the way: its mass (Earth masses), radius (Earth radii) and distance (AU)."""

    mass_earth: float
    radius_earth: float
    distance_au: float


def supply_kg_per_yr(reservoir):
    """What the reservoir sends inward each year, kg: a belt's collisional loss `M / t` times `BELT_DUST_SHARE`; a comet's mass over `COMET_LIFETIME_YR`."""
    if reservoir.kind == "belt":
        return BELT_DUST_SHARE * reservoir.mass_kg / reservoir.age_yr if reservoir.age_yr > 0.0 else 0.0
    return reservoir.mass_kg / COMET_LIFETIME_YR


def reach(reservoir, planet_au):
    """The share of the reservoir's inward supply that is sent inside the planet's orbit: 1 for a planet inside the ring, the part of the ring's width that lies outside the planet's orbit for one inside the ring, 0 beyond it."""
    if reservoir.kind == "comet":
        return 1.0 if planet_au >= reservoir.lower_au else 0.0
    if planet_au <= reservoir.lower_au:
        return 1.0
    if planet_au >= reservoir.upper_au:
        return 0.0
    return (reservoir.upper_au - planet_au) / (reservoir.upper_au - reservoir.lower_au)


def cross_section_share(radius_km, planet_au):
    """The share of an isotropic stream through the sphere of the planet's orbit that the planet intercepts: `R^2 / (4 a^2)`."""
    radius_au = radius_km / constants.AU_TO_KM
    return (radius_au / planet_au) ** 2 / 4.0 if planet_au > 0.0 else 0.0


def grain_beta(star_mass_solar, luminosity_solar=None):
    """The radiation-to-gravity ratio of the belt's grains round a star: `GRAIN_BETA` times the star's luminosity over its mass in solar units (`mass^2.5` for a main-sequence star when `luminosity_solar` is left out)."""
    luminosity = star_mass_solar ** 3.5 if luminosity_solar is None else luminosity_solar
    return GRAIN_BETA * SOLAR_BETA_RATIO * luminosity / star_mass_solar


def crossings(planet_au, star_mass_solar, luminosity_solar=None):
    """How many times a grain of belt dust crosses the orbit radius `planet_au` while Poynting-Robertson drag takes it past: `e c sqrt(a) / (pi beta sqrt(G M))` (twice a period for `2 e a / (da/dt)` of time)."""
    if planet_au <= 0.0 or star_mass_solar <= 0.0:
        return 0.0
    gm = constants.G * star_mass_solar * constants.SOLAR_MASS_TO_KG
    return (GRAIN_ECCENTRICITY * constants.SPEED_OF_LIGHT_M_S * math.sqrt(planet_au * constants.AU_TO_M)
            / (math.pi * grain_beta(star_mass_solar, luminosity_solar) * math.sqrt(gm)))


def shield(giants, star_mass_solar, planet_au, source_au):
    """What is left of the inward flux after every giant between the planet (`planet_au`) and the source (`source_au`) has ejected its share: the product of `exp(-Theta / THETA_EJECT)`."""
    left = 1.0
    for giant in giants:
        if planet_au < giant.distance_au < source_au:
            theta = scatter.safronov_number(giant.mass_earth, giant.radius_earth, giant.distance_au, star_mass_solar)
            left *= math.exp(-theta / scatter.THETA_EJECT)
    return left


def planet_influx_kg_per_yr(radius_km, planet_au, star_mass_solar, reservoirs, giants, luminosity_solar=None):
    """
    The mass, kg each year, of belt dust and comets a planet of `radius_km` at `planet_au` takes in from `reservoirs` (`Reservoir`) past `giants` (`Giant`, the planet itself left out), and the
    share of it that each reservoir gives, as `(total, [(reservoir, kg_per_yr)])`. A belt's dust crosses the planet's orbit `crossings` times (the planet takes `1 - exp(-N share)` of its passes; `luminosity_solar`
    sets the grains' drift), a comet's nucleus once.
    """
    parts = []
    for reservoir in reservoirs:
        source_au = math.sqrt(reservoir.lower_au * reservoir.upper_au) if reservoir.lower_au > 0.0 else reservoir.upper_au
        share = cross_section_share(radius_km, planet_au)
        if reservoir.kind == "belt":
            share = -math.expm1(-crossings(planet_au, star_mass_solar, luminosity_solar) * share)
        flux = supply_kg_per_yr(reservoir) * reach(reservoir, planet_au) * share * shield(giants, star_mass_solar, planet_au, source_au)
        parts.append((reservoir, flux))
    return sum(flux for _, flux in parts), parts
