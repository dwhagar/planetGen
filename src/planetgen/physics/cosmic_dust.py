# planetgen/physics/cosmic_dust.py

"""
The cosmic dust that falls on a planet or moon, kg per year (GEN.1542).

V8 of the Venus sulfur-cloud set (iron-bearing polycyclic aromatic hydrocarbons) is carbon and iron that the interplanetary dust brings to the cloud top: the one paper that
proposes it needs the cloud-top carbon and iron to be 600 and 3 thousand years of the cosmic influx at the least ([R5170], abstract only). So the supply that gates it is the
influx on the body, not the dust its own wind lifts (`aerosols.draw_dust` makes that, and it carries no carbon). The influx follows from the body's size, its mass and its place:

    influx = `EARTH_INFLUX_KG_YR` x (cross-section x focusing) / (the Earth's) x (distance_au)^-`DISTANCE_EXPONENT` x sqrt(M_star / Msun)

* The Earth anchors the scale: about 28 tonnes a day (`EARTH_INFLUX_KG_YR`, [R5220], recalled; the estimates of the different methods run from 5 to 300 tonnes a day).
* The dust a body meets thins with the star's distance: the density of the zodiacal cloud falls as the distance to the power 1.3 (recalled) and the speed
  at which a body meets it rises as the circular speed, which falls as the distance to the power one half, so the mass per area and time falls as `DISTANCE_EXPONENT` 1.8.
* A heavier star gives the grains a higher speed, so the flux rises as the square root of its mass; how much dust a star's belts and comets make is not modelled (every star is the Sun's).
* The gravitational focusing of the body is `ring_origin.focusing` (the same encounter speed as the comets and the ring-erosion meteoroids), and the whole of the infall is taken to reach the air.

The ring's own meteoroid law (`ring_origin.meteoroid_flux_kg_m2_s`, Saturn-anchored with exponent 1.5) gives the Earth 0.7 tonnes a day, 40 times less than the measured 28: it describes the outer
system, so the inner system is anchored here at the Earth.
"""

import math

from planetgen.physics import constants, ring_origin

EARTH = dict(mass_kg=5.972e24, radius_m=6.371e6, distance_au=1.0, star_mass_kg=constants.SOLAR_MASS_TO_KG)
"""dict: The Earth, the body `EARTH_INFLUX_KG_YR` is measured on. Source tier T1, tweak level L0."""

VENUS = dict(mass_kg=4.867e24, radius_m=6.0518e6, distance_au=0.72333, star_mass_kg=constants.SOLAR_MASS_TO_KG)
"""dict: Venus, the one planet whose cloud holds the species V8 is proposed for: the reference of the supply floor (`sulfur_cloud_set` number `dust_min_share_of_reference`). Source tier T1, tweak level L0."""

EARTH_INFLUX_KG_YR = 1.0e7
"""float: The mass of cosmic dust that reaches the Earth, kg per year (28 tonnes a day, [R5220], recalled; other methods give 5 to 300 tonnes a day). Stand-in CV11, tier T3, tweak L0."""

DISTANCE_EXPONENT = 1.8
"""float: The mass flux per area falls with the star's distance to this power: the density of the zodiacal dust (1.3, recalled) times the encounter speed (0.5). Stand-in CV12, tier T3, tweak L1."""

STAR_MASS_EXPONENT = 0.5
"""float: The flux rises with the star's mass to this power, through the encounter speed (the circular speed goes as the square root of the star's mass); the dust the star's belts and comets make is taken as the Sun's. Stand-in CV13, tier T6, tweak L1."""


def influx_kg_yr(mass_kg, radius_m, distance_au, star_mass_kg):
    """The cosmic dust that reaches a body of `mass_kg` and `radius_m` at `distance_au` from a star of `star_mass_kg`, kg per year. 0 for a body with no size or place."""
    if not (mass_kg and radius_m and distance_au and star_mass_kg) or min(mass_kg, radius_m, distance_au, star_mass_kg) <= 0.0:
        return 0.0
    cross_section = (radius_m / EARTH["radius_m"]) ** 2
    focus = ring_origin.focusing(mass_kg, radius_m, distance_au, star_mass_kg) \
        / ring_origin.focusing(EARTH["mass_kg"], EARTH["radius_m"], EARTH["distance_au"], EARTH["star_mass_kg"])
    place = distance_au ** -DISTANCE_EXPONENT * (star_mass_kg / EARTH["star_mass_kg"]) ** STAR_MASS_EXPONENT
    return EARTH_INFLUX_KG_YR * cross_section * focus * place


def influx_per_area_kg_m2_yr(mass_kg, radius_m, distance_au, star_mass_kg):
    """`influx_kg_yr` over the body's whole surface (4 pi r^2), kg per m2 per year."""
    if not radius_m or radius_m <= 0.0:
        return 0.0
    return influx_kg_yr(mass_kg, radius_m, distance_au, star_mass_kg) / (4.0 * math.pi * radius_m ** 2)


def reference_per_area_kg_m2_yr():
    """The per-area influx on Venus, kg per m2 per year: what the supply floor is a share of."""
    return influx_per_area_kg_m2_yr(**VENUS)


def share_of_reference(mass_kg, radius_m, distance_au, star_mass_kg):
    """The body's per-area influx over Venus's."""
    return influx_per_area_kg_m2_yr(mass_kg, radius_m, distance_au, star_mass_kg) / reference_per_area_kg_m2_yr()
