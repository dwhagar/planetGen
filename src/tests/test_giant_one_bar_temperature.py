"""
GEN.816: the stored temperature of a giant, at its 1 bar level, is the greenhouse of the starlight it absorbs and its own heat radiated
together (`surface_air.Setting.internal_flux_w_m2`), not of the starlight alone.
"""

import copy

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import constants, giant_draw, planets
from tests import class_support

# (mass in Earth masses, radius km, distance AU, measured 1 bar K, tolerance, measured internal heat W/m2 or None for the generator's law): three of the
# four giants round a Sun of 4.5 Gy. Uranus is given its measured heat, 0.042 W/m2 [R3466]: the cooling law through Jupiter and Neptune
# gives a Uranus-mass giant 0.33 W/m2, eight times more, an anomaly of Uranus the law does not hold (registry row CC14, GEN.1033), which this test of the
# greenhouse should not carry.
GIANTS = ((317.8, 69911.0, 5.2, 166.0, 0.10, None), (14.5, 25362.0, 19.2, 76.0, 0.10, 0.042), (17.1, 24622.0, 30.07, 72.0, 0.10, None))


def _star():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    star = copy.copy(Star(config))
    star.age = 4.5
    star.luminosity = constants.SOLAR_LUMINOSITY
    return config, star


def _stored_temperature_k(mass_earth, radius_km, distance_au):
    config, star = _star()
    planet = class_support.draw_planet(star, distance_au, 1, config, mass=mass_earth * constants.EARTH_MASS_TO_KG, zone_override="c",
                                       moon_count=0)
    planet.mass = mass_earth * constants.EARTH_MASS_TO_KG
    planet.radius = radius_km
    planet.distance = distance_au
    planets.calculate_surface_gravity(planet)
    planets.calculate_atmospheric_conditions(planet)
    return planet.surface_temperature


@pytest.mark.parametrize("mass_earth, radius_km, distance_au, measured_k, tolerance, internal_w_m2", GIANTS)
def test_the_stored_temperature_lands_near_the_measured_one(monkeypatch, mass_earth, radius_km, distance_au, measured_k, tolerance, internal_w_m2):
    if internal_w_m2 is not None:
        monkeypatch.setattr(giant_draw, "own_flux_w_m2", lambda planet: internal_w_m2)
    assert _stored_temperature_k(mass_earth, radius_km, distance_au) == pytest.approx(measured_k, rel=tolerance)


def test_the_cooling_law_gives_a_uranus_mass_giant_far_more_heat_than_uranus_has(monkeypatch):
    """CC14: the anomaly the Uranus row is spared; when GEN.1033 brings the heat from the interior of each giant, this test goes."""
    config, star = _star()
    planet = class_support.draw_planet(star, 19.2, 1, config, mass=14.5 * constants.EARTH_MASS_TO_KG, zone_override="c", moon_count=0)
    planet.mass, planet.radius = 14.5 * constants.EARTH_MASS_TO_KG, 25362.0
    assert giant_draw.own_flux_w_m2(planet) > 3.0 * 0.042


def test_without_its_own_heat_a_jupiter_would_be_colder(monkeypatch):
    with_heat = _stored_temperature_k(317.8, 69911.0, 5.2)
    monkeypatch.setattr(giant_draw, "own_flux_w_m2", lambda planet: 0.0)
    assert _stored_temperature_k(317.8, 69911.0, 5.2) < with_heat - 10.0
