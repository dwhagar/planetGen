# planetgen/physics/giant_draw.py

"""
A giant's heavy elements, radius and interior, drawn together (GEN.460, GEN.642; docs/design/composition/giant-interiors.md,
giant-heavy-element-law.md)
=================================================================================================================

Called from `planets.generate_planet_properties` once a giant's mass is settled (`planets._apply_giant_mass`). The order is
fixed and is part of the generator epoch:

1. **The envelope's metallicity**: a multiple of the star's own Z, log-uniform over the class's range (`ENVELOPE_Z`).
2. **The heavy-element mass M_Z.** A class J giant takes Chachan et al. 2025's law with its scatter
   (`giant_radius.heavy_mass_earth`); an ice giant or a gas dwarf, whose mass is mostly or partly heavy elements and lies
   below that law's 20 Earth masses, takes a share of its mass from `HEAVY_FRACTION`.
3. **The radius**, unless the caller gave one: where the Fortney, Marley and Barnes 2007 tables hold the body (`giant_radius.in_range`
   and `covers`), their radius for the mass, core, starlight and age plus the Sestovic et al. 2018 inflation; elsewhere the old
   anchored mass-radius relation, which the caller passes as `fallback_radius_km` (a class T, a class I, a giant under 17 Earth
   masses, a moon, a body whose heavy elements outweigh what the tables hold).
4. **The interior** (`giant_interior.draw_interior`): the layers, whose whole-body totals are the body's bulk composition.

The one-bar temperature that decides helium rain is the effective temperature (starlight and the planet's own heat) times
`ONE_BAR_OVER_EFFECTIVE`.
"""

import math

from planetgen import tuning
from planetgen.physics import bulk_composition, bulk_draw, climate_draws, constants, giant_interior, giant_radius
from planetgen.util import draw
from planetgen.util.random import log_uniform

SOURCES = {
    "chachan2025": "Chachan et al. 2025, arXiv:2509.20428 (the heavy-element mass of giants)",
    "fortney2007": "Fortney, Marley and Barnes 2007, ApJ 659:1661 (giant radii, Tables 2 to 4)",
    "sestovic2018": "Sestovic, Demory and Queloz 2018, A&A 612:A74 (hot-giant inflation)",
}

HEAVY_FRACTION = {"I": (0.75, 0.95), "T": (0.10, 0.50)}
"""dict: A class's heavy-element share of its mass, uniform, for the classes below the heavy-element law's 20 Earth masses
(Chachan et al. 2025 give M_Z = M below 14.7); the complement of the stand-in envelope shares that stood here before. Ice
giants are 5 to 25 percent hydrogen and helium (Uranus and Neptune about 10 to 20 percent, Helled and Fortney 2020, READ as a
range, not a fit). Class J takes the law."""

ENVELOPE_Z = {"I": (5.0, 30.0), "J": (1.0, 8.0), "T": (1.0, 12.0)}
"""dict: The envelope's heavy-element mass fraction as a multiple of the star's Z, log-uniform; a stand-in range (Jupiter's
atmosphere is 3 to 4 times solar, Saturn's about 10, Uranus and Neptune 20 to 100 in carbon) until a source gives the law."""

MAX_ENVELOPE_Z = 0.5
"""float: The most heavy elements an envelope holds, as a fraction of its mass."""

WATER_SHARE_OF_HEAVY = {"I": (0.30, 0.60), "J": (0.10, 0.30), "T": (0.20, 0.40)}
"""dict: For a giant born beyond the snow line, the water's share of its heavy elements, uniform; a stand-in range. Ammonia and
methane follow at `AMMONIA_OVER_WATER` and `METHANE_OVER_WATER` of the water."""
AMMONIA_OVER_WATER = 0.10
METHANE_OVER_WATER = 0.15

BOND_ALBEDO = 0.503
"""float: Jupiter's Bond albedo (Li et al. 2018, Nature Communications 9:3709, the newer of two measurements; Hanel et al. 1981 had
0.343), taken for every giant until the albedo is drawn by class."""

INTERNAL_TEMPERATURE_K = 100.0
"""float: The temperature of Jupiter's own heat (Li et al. 2018: 100 K), scaled by `(M / Jupiter)^0.5` in `effective_temperature_k`,
a stand-in until the cooling tracks of the radius tables give the heat of a body of this mass and age."""

ONE_BAR_OVER_EFFECTIVE = 1.33
"""float: The temperature at one bar over the effective temperature (Jupiter: 166 K over 125 K; Saturn: 134 K over 98.8 K, which
is 1.36), a stand-in."""

MOON_ONE_BAR_K = 130.0
"""float: The one-bar temperature taken for a giant moon, whose distance from the star the body does not carry (a stand-in)."""

REFERENCE_FRACTIONS = (0.08, 0.07)
"""tuple: The mantle's FeO share and the core's light-element share used to read the iron share of the rock-and-iron from the
star's ratios (`bulk_draw.dry_partition`), so that no further draw is made."""


def starlight_w_m2(planet):
    """The starlight at a planet (not a moon), W/m2; `None` for a moon, which does not carry its parent's distance."""
    if planet.is_moon:
        return None
    return climate_draws.starlight_w_m2(planet.star.luminosity, planet.distance * constants.AU_TO_KM)


def effective_temperature_k(flux_w_m2, mass_earth):
    """The temperature a giant radiates at: its starlight (`BOND_ALBEDO`) and its own heat added in the fourth power."""
    own_k = INTERNAL_TEMPERATURE_K * (mass_earth / giant_interior.EARTH_MASS_JUPITER) ** 0.5
    starlight_k = climate_draws.equilibrium_temperature_k(BOND_ALBEDO, flux_w_m2)
    return (starlight_k ** 4 + own_k ** 4) ** 0.25


def _heavy_mass_earth(planet_class, mass_earth):
    if planet_class == "J":
        return giant_radius.heavy_mass_earth(mass_earth)
    low, high = HEAVY_FRACTION[planet_class]
    return draw.uniform(low, high) * mass_earth


def _radius_km(planet, mass_earth, heavy_earth, envelope_z, flux_w_m2, fallback_radius_km):
    if flux_w_m2 is not None and giant_radius.in_range(mass_earth):
        core_earth = giant_radius.table_core_mass_earth(mass_earth, heavy_earth, envelope_z)
        if giant_radius.covers(mass_earth, core_earth):
            flux_earth = flux_w_m2 / giant_radius.flux_w_m2(1.0)
            low, high = tuning.PLANET_CLASSES[planet.planet_class]["radius_range"]
            radius_km = giant_radius.radius_km(mass_earth, core_earth, flux_earth, planet.star.age)
            return min(high, max(low, radius_km))
    return fallback_radius_km(planet.planet_class, planet.mass)


def _ice_shares(planet):
    if planet.is_moon or planet.distance < planet.star.snow_line_au:
        return {}
    low, high = WATER_SHARE_OF_HEAVY[planet.planet_class]
    water = draw.uniform(low, high)
    return {"water_liquid": water, "nh3_ice": water * AMMONIA_OVER_WATER, "ch4_ice": water * METHANE_OVER_WATER}


def draw_giant(planet, radius_given, fallback_radius_km):
    """
    Draws a giant's heavy elements, radius (unless given) and interior, and stores its bulk composition and density.

    Args:
        planet (Planet): A giant class (`planets.uses_giant_mass_radius`) with its mass set (kg), its class and `star`.
        radius_given (bool): The caller gave the radius (km); it is kept.
        fallback_radius_km (callable): `(planet_class, mass_kg) -> km`, the anchored mass-radius relation (with its scatter,
            kept inside the class's range) for a body the tables do not hold.

    Sets:
        `radius` (unless given), `density` (g/cm3) and the attributes of `bulk_composition.BODY_FIELDS`.
    """
    mass_earth = planet.mass / constants.EARTH_MASS_TO_KG
    ratios = bulk_draw.host_ratios(planet.star)
    z_low, z_high = ENVELOPE_Z[planet.planet_class]
    envelope_z = min(MAX_ENVELOPE_Z, ratios.z * log_uniform(z_low, z_high))
    heavy_earth = _heavy_mass_earth(planet.planet_class, mass_earth)
    flux_w_m2 = starlight_w_m2(planet)
    if not radius_given:
        planet.radius = _radius_km(planet, mass_earth, heavy_earth, envelope_z, flux_w_m2, fallback_radius_km)
    volume_m3 = (4 / 3) * math.pi * (planet.radius * constants.KM_TO_M_FACTOR) ** 3
    planet.density = planet.mass / volume_m3 / 1000

    part = bulk_draw.dry_partition(ratios, *REFERENCE_FRACTIONS, False)
    temperature_k = (MOON_ONE_BAR_K if flux_w_m2 is None
                     else ONE_BAR_OVER_EFFECTIVE * effective_temperature_k(flux_w_m2, mass_earth))
    subject = giant_interior.GiantSubject(
        mass_earth=mass_earth, radius_km=planet.radius, heavy_mass_earth=heavy_earth, envelope_z=envelope_z,
        helium_share=ratios.y_envelope, temperature_1bar_k=temperature_k, iron_share_of_refractory=part.metal,
        ice_shares=_ice_shares(planet))
    interior = giant_interior.draw_interior(subject)
    shares = interior.whole_body_shares()
    bulk = bulk_composition.BulkFractions(
        metal=shares["metal"], rock=shares["rock"], water=shares["water"], ice=shares["ice"], envelope=shares["envelope"],
        envelope_helium=shares["envelope_helium"], envelope_z=shares["envelope_z"], origin="drawn")
    for attribute, value in bulk.to_attributes().items():
        setattr(planet, attribute, value)
