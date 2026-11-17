# planetgen/physics/influx.py

"""
The dust and comet influx a planet takes from every belt and comet of its system (GEN.1440, GEN.1453 slice A)
==============================================================================================================

Until now the only reservoir a planet's impacts knew was the late-accretion number of one Solar-System-sized leftover (`impact_history`, GEN.739). A planet's real supply is what the belts and comets of its
own system send in, less what a giant between them takes. This module works that out as a pure function of stored data: the reservoirs (belts with their mass and ring, comets with their mass and
perihelion), the giants (mass, radius, distance), the star's age and the planet's radius and distance. Nothing is stored and nothing is drawn; a sheet works it out when it is made (the dust is a pure
function of the stored belts and planets, so a column would only be a second copy of them that a rescan could leave stale).

**The law (a stand-in throughout; rows IF1 to IF10 of the gap plan).**

    influx = sum over reservoirs j of  supply_j  x  reach_j(planet)  x  capture_share(planet)  x  shield_j(planet)

* *Supply.* A belt loses mass to collisions as `M / t` at its age `t` (the collisional cascade `M = A / t` of `belt_mass`, so `|dM/dt| = M / t`); the share `BELT_DUST_SHARE` of that loss is dust or
  small bodies fine enough to drift inward past the belt's edge (Poynting-Robertson drag and resonance leaks). A comet reservoir sends its nucleus in over `COMET_LIFETIME_YR`.
* *Reach.* Only what drifts inward arrives: a planet inside the belt's outer edge `upper` is reached, one outside it is not (nothing is sent outward). A planet inside the ring itself gets the share of
  the ring that lies outside it (the ring is shared by the part of its width the planet's orbit leaves outside). A comet reaches a planet whose orbit its perihelion passes inside.
* *Capture* (GEN.1543). The dust does not fill the sphere of the planet's orbit. It is a thin disk: grains of radius `GRAIN_RADIUS_M` leave the belt on orbits of eccentricity `DUST_ECCENTRICITY` and inclination
  `DUST_INCLINATION_RAD` and spiral inward under Poynting-Robertson drag, at the speed `drift_speed_m_s` the star's light gives them (`radiation_ratio`, the ratio of the light's push to the star's pull; the drag law
  `da/dt = -2 beta G M / (c a)`). The steady flow `supply` through the annulus at the planet's distance, of thickness `2 a sin(i)`, gives a mass density `supply / (2 pi a H v_drift)`, and the planet sweeps it
  at the speed `v_rel = v_orbit sqrt(e^2 + i^2)` of the grains past it (the encounter speed of two near-circular orbits) over its focused cross-section `pi R^2 (1 + v_esc^2 / v_rel^2)`. The share of the stream a planet takes
  is that rate over the supply: `capture_share`. The old law, a sphere's `R^2 / (4 a^2)`, spread the dust over every direction and put the Earth at 45 kg a year against the 1e7 kg a year measured.
* *Shield.* Each giant that lies between the reservoir and the planet ejects a share of what drifts past it, `1 - exp(-Theta / THETA_EJECT)` with `Theta` its Safronov number (`scatter`, the rule the
  swept disk of the giants follows); what is left is the product over the giants in the way.

What it leaves out (stated; rows IF5 and IF6): the giants' scattering of comets inward (a giant that throws a comet in as well as out), a companion star's stirring of a reservoir (GEN.1492), the
outward flux of a belt inside the planet's orbit, the grains' loss to collisions on the way in and their spread of sizes (one size carries the supply), and the dust a planet's own moons and rings make.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, scatter

BELT_DUST_SHARE = 0.03
"""float: The share of a belt's collisional mass loss that reaches the inner system as dust and small bodies. A stand-in, in the range of the share of a cascade's mass that is in grains small enough for
radiation forces (the rest is ground to dust that leaves or stays). 0.03 is the share that puts the Sun's own main belt, at the grain size and geometry below, at the 1e7 kg a year the Earth takes ([R1397], recalled; 28 tonnes a day). [Stand-in 6.247 IF1, tier T6, tweak L2.]"""

COMET_LIFETIME_YR = 1.0e7
"""float: The time over which a comet nucleus of a reservoir delivers its mass to the planets (its dynamical lifetime before ejection or sublimation, in the tens of millions of years for the
short-period and Halley kinds). [Stand-in 6.247 IF2, tier T6, tweak L4.]"""

GRAIN_RADIUS_M = 3.0e-5
"""float: The radius, m, of the grain that carries a reservoir's dust inward (30 micrometers: the size the zodiacal cloud's light and its Poynting-Robertson drift are set by). One size stands for the spread of sizes the cascade gives. [Stand-in 6.247 IF7, tier T6, tweak L1.]"""

GRAIN_DENSITY_KG_M3 = 2500.0
"""float: The density of a grain, kg per cubic meter: stony dust. [Stand-in 6.247 IF8, tier T6, tweak L0.]"""

DUST_ECCENTRICITY = 0.15
"""float: The eccentricity the dust has when it crosses a planet's orbit: the belt's own scatter, kept by the grains (Poynting-Robertson drag circularizes them over a time much longer than the crossing). [Stand-in 6.247 IF9, tier T6, tweak L1.]"""

DUST_INCLINATION_RAD = 0.25
"""float: The inclination, radians, of the dust to the plane (14 degrees): the disk's thickness is `2 a sin(i)`. The comet-fed dust of the Solar System stands at 10 to 30 degrees. [Stand-in 6.247 IF10, tier T6, tweak L1.]"""


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


def radiation_ratio(star_mass_solar, star_luminosity_solar):
    """`beta`, the ratio of the light's push on a grain of `GRAIN_RADIUS_M` to the star's pull: `3 L / (16 pi G M c rho s)` (a perfectly absorbing grain)."""
    mass_kg = star_mass_solar * constants.SOLAR_MASS_TO_KG
    return 3.0 * star_luminosity_solar * constants.SOLAR_LUMINOSITY / (16.0 * math.pi * constants.G * mass_kg * constants.SPEED_OF_LIGHT_M_S * GRAIN_DENSITY_KG_M3 * GRAIN_RADIUS_M)


def drift_speed_m_s(planet_au, star_mass_solar, star_luminosity_solar):
    """How fast a grain spirals in at `planet_au`, m/s: Poynting-Robertson drag, `2 beta G M / (c a)`."""
    mass_kg = star_mass_solar * constants.SOLAR_MASS_TO_KG
    return 2.0 * radiation_ratio(star_mass_solar, star_luminosity_solar) * constants.G * mass_kg / (constants.SPEED_OF_LIGHT_M_S * planet_au * constants.AU_TO_M)


def capture_share(radius_km, mass_kg, planet_au, star_mass_solar, star_luminosity_solar):
    """
    The share of the dust that drifts in past `planet_au` that the planet takes (GEN.1543): its focused cross-section times the grains' speed past it, over the annulus the stream fills times the speed it drifts at,
    `pi R^2 (1 + v_esc^2 / v_rel^2) v_rel / (2 pi a H v_drift)` with `H = 2 a sin(i)` and `v_rel = v_orbit sqrt(e^2 + i^2)`; at most 1.

    Args:
        radius_km, mass_kg (float): The planet's radius and mass.
        planet_au (float): Its distance.
        star_mass_solar, star_luminosity_solar (float): The star's mass and luminosity.
    """
    if planet_au <= 0.0 or radius_km <= 0.0 or mass_kg <= 0.0 or star_luminosity_solar <= 0.0:
        return 0.0
    a_m = planet_au * constants.AU_TO_M
    star_kg = star_mass_solar * constants.SOLAR_MASS_TO_KG
    v_orbit = math.sqrt(constants.G * star_kg / a_m)
    v_rel = v_orbit * math.hypot(DUST_ECCENTRICITY, DUST_INCLINATION_RAD)
    radius_m = radius_km * 1000.0
    v_escape_sq = 2.0 * constants.G * mass_kg / radius_m
    cross_section = math.pi * radius_m ** 2 * (1.0 + v_escape_sq / v_rel ** 2)
    thickness = 2.0 * a_m * math.sin(DUST_INCLINATION_RAD)
    return min(1.0, cross_section * v_rel / (2.0 * math.pi * a_m * thickness * drift_speed_m_s(planet_au, star_mass_solar, star_luminosity_solar)))


def shield(giants, star_mass_solar, planet_au, source_au):
    """What is left of the inward flux after every giant between the planet (`planet_au`) and the source (`source_au`) has ejected its share: the product of `exp(-Theta / THETA_EJECT)`."""
    left = 1.0
    for giant in giants:
        if planet_au < giant.distance_au < source_au:
            theta = scatter.safronov_number(giant.mass_earth, giant.radius_earth, giant.distance_au, star_mass_solar)
            left *= math.exp(-theta / scatter.THETA_EJECT)
    return left


def planet_influx_kg_per_yr(radius_km, mass_kg, planet_au, star_mass_solar, star_luminosity_solar, reservoirs, giants):
    """
    The mass, kg each year, of belt dust and comets a planet of `radius_km` and `mass_kg` at `planet_au` takes in from `reservoirs` (`Reservoir`) past `giants` (`Giant`, the planet itself left out), round a star of
    `star_mass_solar` and `star_luminosity_solar`, and the share of it that each reservoir gives, as `(total, [(reservoir, kg_per_yr)])`.
    """
    capture = capture_share(radius_km, mass_kg, planet_au, star_mass_solar, star_luminosity_solar)
    parts = []
    for reservoir in reservoirs:
        source_au = math.sqrt(reservoir.lower_au * reservoir.upper_au) if reservoir.lower_au > 0.0 else reservoir.upper_au
        flux = (supply_kg_per_yr(reservoir) * reach(reservoir, planet_au) * capture
                * shield(giants, star_mass_solar, planet_au, source_au))
        parts.append((reservoir, flux))
    return sum(flux for _, flux in parts), parts
