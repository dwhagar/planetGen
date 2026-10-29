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

# (mass in Earth masses, radius km, distance AU, measured 1 bar K, tolerance): the four giants round a Sun of 4.5 Gy.
GIANTS = ((317.8, 69911.0, 5.2, 166.0, 0.10), (14.5, 25362.0, 19.2, 76.0, 0.12), (17.1, 24622.0, 30.07, 72.0, 0.10))


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


@pytest.mark.parametrize("mass_earth, radius_km, distance_au, measured_k, tolerance", GIANTS)
def test_the_stored_temperature_lands_near_the_measured_one(mass_earth, radius_km, distance_au, measured_k, tolerance):
    assert _stored_temperature_k(mass_earth, radius_km, distance_au) == pytest.approx(measured_k, rel=tolerance)


def test_without_its_own_heat_a_jupiter_would_be_colder(monkeypatch):
    with_heat = _stored_temperature_k(317.8, 69911.0, 5.2)
    monkeypatch.setattr(giant_draw, "own_flux_w_m2", lambda planet: 0.0)
    assert _stored_temperature_k(317.8, 69911.0, 5.2) < with_heat - 10.0
