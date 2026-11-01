# planetgen/physics/giant_draw.py

"""
A giant's heavy elements, radius and interior, drawn together (GEN.460, GEN.642; docs/design/composition/giant-interiors.md,
giant-heavy-element-law.md)
=================================================================================================================

Called from `planets.generate_planet_properties` once a giant's mass is settled (`planets._apply_giant_mass`). A giant has one of two
regimes by its mass (`regime_of`): an ice giant (at most `ICE_GIANT_MAX_EARTH`, a core that found no gas) or a gas giant. The order is
fixed and is part of the generator epoch:

1. **The envelope's metallicity**: a multiple of the star's own Z, log-uniform over the regime's range (`ENVELOPE_Z`).
2. **The heavy-element mass M_Z.** A gas giant takes [R0710]'s law with its scatter
   (`giant_radius.heavy_mass_earth`); an ice giant, whose mass is mostly heavy elements and lies below that law's 20 Earth masses,
   takes a share of its mass from `HEAVY_FRACTION`.
3. **The radius**, unless the caller gave one: where the [R1429] tables hold the body (`giant_radius.in_range`
   and `covers`), their radius for the mass, core, starlight and age plus the [R4029] inflation; elsewhere the old
   anchored mass-radius relation, which the caller passes as `fallback_radius_km` (an ice giant, a giant under 17 Earth
   masses, a moon, a body whose heavy elements outweigh what the tables hold).
4. **The interior** (`giant_interior.draw_interior`): the layers, whose whole-body totals are the body's bulk composition.

The star's age is final only after the planets are (`Star.adjust_age_for_planets`), so a table radius is drawn at a provisional age and
`settle` takes it again once the age and the orbit are final (GEN.961); the interior, which is made to that radius, is drawn again
with it.

The one-bar temperature that decides helium rain is the effective temperature (starlight and the planet's own heat) times
`ONE_BAR_OVER_EFFECTIVE`.
"""

import math

from planetgen import tuning
from planetgen.physics import (bulk_composition, bulk_draw, climate_draws, composition_text, constants, giant_interior, giant_radius,
                               orbit_elements, rogue_surface)
from planetgen.util import draw
from planetgen.util.random import log_uniform
from planetgen.physics import unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")

SOURCES = {
    "chachan2025": "[R0710] (the heavy-element mass of giants)",
    "fortney2007": "[R1429]:1661 (giant radii, Tables 2 to 4)",
    "sestovic2018": "[R4029]:A74 (hot-giant inflation)",
}

ICE, GAS = "ice", "gas"
"""str: The two regimes of a giant (`regime_of`)."""
ICE_GIANT_MAX_EARTH = tuning.GIANT_NO_GAS_CORE_MAX_EARTH
"""float: The heaviest ice giant, Earth masses: a core that found no gas stops growing here (`planet_mass`); heavier is a gas giant."""

HEAVY_FRACTION = (0.75, 0.95)
"""tuple: An ice giant's heavy-element share of its mass, uniform, below the heavy-element law's 20 Earth masses
([R0710] give M_Z = M below 14.7). Ice giants are 5 to 25 percent hydrogen and helium (Uranus and Neptune about 10 to 20
percent, [R1892], READ as a range, not a fit). A gas giant takes the law."""

ENVELOPE_Z = {ICE: (5.0, 30.0), GAS: (1.0, 8.0)}
"""dict: The envelope's heavy-element mass fraction as a multiple of the star's Z, log-uniform, by regime; a stand-in range (Jupiter's
atmosphere is 3 to 4 times solar, Saturn's about 10, Uranus and Neptune 20 to 100 in carbon) until a source gives the law."""

MAX_ENVELOPE_Z = 0.5
"""float: The most heavy elements an envelope holds, as a fraction of its mass."""

WATER_SHARE_OF_HEAVY = {ICE: (0.30, 0.60), GAS: (0.10, 0.30)}
"""dict: For a giant born beyond the snow line, the water's share of its heavy elements, uniform; a stand-in range. Ammonia and
methane follow at `AMMONIA_OVER_WATER` and `METHANE_OVER_WATER` of the water."""
AMMONIA_OVER_WATER = 0.10
METHANE_OVER_WATER = 0.15

BOND_ALBEDO = 0.503
"""float: Jupiter's Bond albedo (Li et al. 2018, Nature Communications 9:3709, the newer of two measurements), taken for every giant's effective temperature."""

ONE_BAR_OVER_EFFECTIVE = 1.33
"""float: The temperature at one bar over the effective temperature (Jupiter: 166 K over 125 K; Saturn: 134 K over 98.8 K, which
is 1.36), a stand-in."""

MOON_ONE_BAR_K = 130.0
"""float: The one-bar temperature taken for a giant moon, whose distance from the star the body does not carry (a stand-in)."""

SETTLE_TOLERANCE = 0.001
"""float: A table giant whose radius at the final age and the orbit's mean light is within this share of the one it was drawn with
keeps its radius and interior (`settle`); the tables themselves are good to a few percent. A stand-in, registry row OE13."""

SETTLE_ATTRIBUTE = "_giant_draw"
"""str: The attribute of a giant that `draw_giant` leaves for `settle` (`DrawnGiant`); not stored, it lives for the generation."""

REFERENCE_FRACTIONS = (0.08, 0.07)
"""tuple: The mantle's FeO share and the core's light-element share used to read the iron share of the rock-and-iron from the
star's ratios (`bulk_draw.dry_partition`), so that no further draw is made."""


def starlight_w_m2(planet):
    """The starlight at a planet (not a moon), W/m2; `None` for a moon, which does not carry its parent's distance."""
    if planet.is_moon:
        return None
    return climate_draws.starlight_w_m2(planet.star.luminosity, planet.distance * constants.AU_TO_KM)


def own_flux_w_m2(planet):
    """A giant's own heat flux at its 1 bar level, W/m2, at its star's age: the cooling law through Jupiter and Neptune of
    `rogue_surface.giant_internal_flux_w_m2` (L ~ M^1.65 below a Jupiter mass, t^-1.3 in age; GEN.816). It stands in for the heat a
    body of this mass and age has left until the cooling tracks of the radius tables give it, and it has no helium-rain share yet
    (GEN.810), which is why a Saturn-like giant comes out cooler than Saturn."""
    age_gy = getattr(getattr(planet, "star", None), "age", None) or tuning.ROGUE_EARTH_REFERENCE_AGE_GY
    return rogue_surface.giant_internal_flux_w_m2(planet.mass, planet.radius, age_gy)


def effective_temperature_k(flux_w_m2, own_w_m2):
    """The temperature a giant radiates at, K: the starlight it absorbs (`BOND_ALBEDO`) and its own heat flux `own_w_m2` added as
    fluxes."""
    starlight_k = climate_draws.equilibrium_temperature_k(BOND_ALBEDO, flux_w_m2)
    return (starlight_k ** 4 + own_w_m2 / constants.STEFAN_BOLTZMANN_CONSTANT) ** 0.25


def regime_of(mass_earth):
    """`ICE` for a giant of at most `ICE_GIANT_MAX_EARTH` Earth masses, else `GAS`."""
    return ICE if mass_earth <= ICE_GIANT_MAX_EARTH else GAS


def _heavy_mass_earth(regime, mass_earth):
    if regime == GAS:
        return giant_radius.heavy_mass_earth(mass_earth)
    low, high = HEAVY_FRACTION
    return draw.uniform(low, high) * mass_earth


def _table_core_earth(mass_earth, heavy_earth, envelope_z):
    """The core mass (Earth masses) the radius tables are read at, or `None` where they do not hold the body."""
    if not giant_radius.in_range(mass_earth):
        return None
    core_earth = giant_radius.table_core_mass_earth(mass_earth, heavy_earth, envelope_z)
    return core_earth if giant_radius.covers(mass_earth, core_earth) else None


def _table_radius_km(mass_earth, core_earth, flux_w_m2, age_gy):
    """The tables' radius (no draw) at a starlight in W/m2 and an age."""
    return giant_radius.radius_km(mass_earth, core_earth, flux_w_m2 / giant_radius.flux_w_m2(1.0), age_gy)


def _ice_shares(planet, regime):
    if planet.is_moon or planet.distance < planet.star.snow_line_au:
        return {}
    low, high = WATER_SHARE_OF_HEAVY[regime]
    water = draw.uniform(low, high)
    return {"water_liquid": water, "nh3_ice": water * AMMONIA_OVER_WATER, "ch4_ice": water * METHANE_OVER_WATER}


class DrawnGiant:
    """What `settle` needs of a giant whose radius came from the tables: the star's ratios, the envelope's Z and the heavy-element mass it
    was drawn with (the ratios are read from the star's composition at the age of the draw, and a later age reads another draw of
    them, so they are kept) and the radius it was drawn with."""

    __slots__ = ("ratios", "envelope_z", "heavy_earth", "radius_km")

    def __init__(self, ratios, envelope_z, heavy_earth, radius_km):
        self.ratios, self.envelope_z, self.heavy_earth, self.radius_km = ratios, envelope_z, heavy_earth, radius_km


def draw_giant(planet, radius_given, fallback_radius_km):
    """
    Draws a giant's heavy elements, radius (unless given) and interior, and stores its bulk composition and density.

    Args:
        planet (Planet): A giant (`body_type` "g") with its mass set (kg) and `star`.
        radius_given (bool): The caller gave the radius (km); it is kept.
        fallback_radius_km (callable): `(mass_kg) -> km`, the anchored mass-radius relation (with its scatter) for a body the
            tables do not hold.

    Sets:
        `radius` (unless given), `density` (g/cm3), the attributes of `bulk_composition.BODY_FIELDS` and the interior's
        (`giant_interior.body_values`, DB.89), and for a giant sized from the tables the `DrawnGiant` that `settle` reads.
    """
    mass_earth = _kg_to_mearth(planet.mass)
    ratios = bulk_draw.host_ratios(planet.star)
    regime = regime_of(mass_earth)
    z_low, z_high = ENVELOPE_Z[regime]
    envelope_z = min(MAX_ENVELOPE_Z, ratios.z * log_uniform(z_low, z_high))
    heavy_earth = _heavy_mass_earth(regime, mass_earth)
    flux_w_m2 = starlight_w_m2(planet)
    setattr(planet, SETTLE_ATTRIBUTE, None)
    if not radius_given:
        core_earth = None if flux_w_m2 is None else _table_core_earth(mass_earth, heavy_earth, envelope_z)
        if core_earth is None:
            planet.radius = fallback_radius_km(planet.mass)
        else:
            planet.radius = _table_radius_km(mass_earth, core_earth, flux_w_m2, planet.star.age)
            setattr(planet, SETTLE_ATTRIBUTE, DrawnGiant(ratios, envelope_z, heavy_earth, planet.radius))
    _store_body(planet, mass_earth, ratios, regime, heavy_earth, envelope_z, flux_w_m2)


def _store_body(planet, mass_earth, ratios, regime, heavy_earth, envelope_z, flux_w_m2):
    """Stores the density of the radius the giant has, draws its interior to that radius and stores the bulk composition, the text
    and the interior's values."""
    volume_m3 = (4 / 3) * math.pi * (planet.radius * constants.KM_TO_M_FACTOR) ** 3
    planet.density = planet.mass / volume_m3 / 1000

    part = bulk_draw.dry_partition(ratios, *REFERENCE_FRACTIONS, False)
    temperature_k = (MOON_ONE_BAR_K if flux_w_m2 is None
                     else ONE_BAR_OVER_EFFECTIVE * effective_temperature_k(flux_w_m2, own_flux_w_m2(planet)))
    subject = giant_interior.GiantSubject(
        mass_earth=mass_earth, radius_km=planet.radius, heavy_mass_earth=heavy_earth, envelope_z=envelope_z,
        helium_share=ratios.y_envelope, temperature_1bar_k=temperature_k, iron_share_of_refractory=part.metal,
        ice_shares=_ice_shares(planet, regime))
    interior = giant_interior.draw_interior(subject)
    shares = interior.whole_body_shares()
    bulk = bulk_composition.BulkFractions(
        metal=shares["metal"], rock=shares["rock"], water=shares["water"], ice=shares["ice"], envelope=shares["envelope"],
        envelope_helium=shares["envelope_helium"], envelope_z=shares["envelope_z"], origin="drawn")
    for attribute, value in bulk.to_attributes().items():
        setattr(planet, attribute, value)
    planet.composition = composition_text.text(planet)  # GEN.430: the text states the stored fractions
    for attribute, value in giant_interior.body_values(interior).items():
        setattr(planet, attribute, value)


def settle_stream(planet, drawn):
    """The stream a settled giant's interior is drawn from: its own, keyed by what it is once the star's age is final (the run's stream
    is not touched, so nothing drawn after it moves). The key is the giant's mass, distance, star and the age and radius it settles
    from, which no other giant shares."""
    return draw.Stream(f"settle-giant:{float(planet.mass)!r}:{float(planet.distance)!r}:{float(planet.star.mass)!r}:"
                       f"{float(planet.star.age)!r}:{drawn.radius_km!r}")


def settle(planet):
    """
    Takes a table-sized giant's radius again at the star's final age and the mean light of the orbit, and draws its interior again to it
    (GEN.961, GEN.976).

    A giant is drawn when the star's age is provisional, and the eccentricity it is lit by comes after the radius; the age rises when the
    planets are in (`Star.adjust_age_for_planets`) and the table radius falls with it, by up to 9 percent. The radius is the tables' (no
    draw) at the age now and the light of `orbit_elements.flux_distance_au`. A giant that the tables moved by under `SETTLE_TOLERANCE`, one
    sized any other way and a moon are left as drawn. Otherwise the radius, the volume, the density and the gravity are set again, and the
    interior (the layers, the bulk composition and the text) is drawn from `settle_stream`. The star's ratios, the envelope's Z and the
    heavy-element mass are the original draws (`DrawnGiant`).

    Returns:
        bool: Whether the giant's radius was moved. The caller makes the air again (it reads the gravity).
    """
    drawn = getattr(planet, SETTLE_ATTRIBUTE, None)
    if drawn is None or planet.is_moon:
        return False
    mass_earth = _kg_to_mearth(planet.mass)
    regime = regime_of(mass_earth)
    core_earth = _table_core_earth(mass_earth, drawn.heavy_earth, drawn.envelope_z)
    if core_earth is None:
        return False
    flux_w_m2 = climate_draws.starlight_w_m2(
        planet.star.luminosity, orbit_elements.flux_distance_au(planet.distance, planet.eccentricity) * constants.AU_TO_KM)
    radius_km = _table_radius_km(mass_earth, core_earth, flux_w_m2, planet.star.age)
    if abs(radius_km - drawn.radius_km) <= SETTLE_TOLERANCE * drawn.radius_km:
        return False
    stream = settle_stream(planet, drawn)
    planet.radius = radius_km
    with draw.bound(stream):
        _store_body(planet, mass_earth, drawn.ratios, regime, drawn.heavy_earth, drawn.envelope_z, flux_w_m2)
    setattr(planet, SETTLE_ATTRIBUTE, DrawnGiant(drawn.ratios, drawn.envelope_z, drawn.heavy_earth, radius_km))
    return True
