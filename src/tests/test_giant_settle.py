# tests/test_giant_settle.py

"""
GEN.961 and GEN.976: a giant sized from the radius tables is drawn when the star's age is provisional, and sized again when the age and the
orbit are final (`planets.settle_giant`, the "giant settle pass" of `StarSystem`). The unit rules of `giant_draw.settle` are in
test_giant_draw.py; here the planet it is wired to and the finished system.
"""

import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, giant_draw, giant_radius, orbit_elements, planets
from planetgen.util import draw
from tests import class_support

GIANT_KG = 100.0 * constants.EARTH_MASS_TO_KG
SYSTEMS = 40
SEED = 7


@pytest.fixture()
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("giant-settle:star"):
        return Star(cfg)


def _table_giant(star, seeds=40):
    for seed in range(seeds):
        planet = class_support.draw_planet(star, 0.4, seed, SystemConfig(), mass=GIANT_KG, zone_override="h", moon_count=0)
        if getattr(planet, giant_draw.SETTLE_ATTRIBUTE, None) is not None:
            return planet
    pytest.fail("no giant sized from the tables in the draws")


def _table_radius_km(planet):
    drawn = getattr(planet, giant_draw.SETTLE_ATTRIBUTE)
    mass_earth = planet.mass / constants.EARTH_MASS_TO_KG
    core = giant_radius.table_core_mass_earth(mass_earth, drawn.heavy_earth, drawn.envelope_z)
    distance_km = orbit_elements.flux_distance_au(planet.distance, planet.eccentricity) * constants.AU_TO_KM
    light = planet.star.luminosity / (4.0 * math.pi * distance_km ** 2) / 1e6
    return giant_radius.radius_km(mass_earth, core, light / giant_radius.flux_w_m2(1.0), planet.star.age)


def test_settling_a_giant_sets_its_radius_volume_gravity_and_air_together(host_star):
    planet = _table_giant(host_star)
    young = planet.radius
    host_star.age = max(host_star.age, 1.0) + 6.0
    assert planets.settle_giant(planet) is True
    assert planet.radius < young
    assert planet.radius == pytest.approx(_table_radius_km(planet), rel=1e-9)
    assert planet.volume == pytest.approx(4 / 3 * math.pi * planet.radius ** 3, rel=1e-12)
    assert planet.gravity == pytest.approx(constants.G * planet.mass / (planet.radius * 1000) ** 2 / constants.EARTH_GRAVITY, rel=1e-6)
    assert planet.density == pytest.approx(planet.mass / (4 / 3 * math.pi * (planet.radius * 1000) ** 3) / 1000, rel=1e-9)
    assert planet.surface_temperature is not None and planet.scale_height is not None
    assert planet.planet_class is not None


def test_a_settled_giant_is_settled_for_good(host_star):
    planet = _table_giant(host_star)
    host_star.age = max(host_star.age, 1.0) + 6.0
    planets.settle_giant(planet)
    stored = (planet.radius, planet.gravity, planet.scale_height, planet.composition)
    assert planets.settle_giant(planet) is False
    assert (planet.radius, planet.gravity, planet.scale_height, planet.composition) == stored


def test_settling_draws_nothing_from_the_run_stream(host_star):
    planet = _table_giant(host_star)
    host_star.age = max(host_star.age, 1.0) + 6.0
    draw.set_run_seed(5)
    planets.settle_giant(planet)
    after = draw.random()
    draw.set_run_seed(5)
    assert after == draw.random()


@pytest.mark.slow
def test_every_table_giant_of_a_finished_system_has_the_radius_of_its_final_age_and_orbit():
    draw.set_run_seed(SEED)
    checked = 0
    for _ in range(SYSTEMS):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        for planet in StarSystem(system_config=cfg).all_planets:
            if getattr(planet, giant_draw.SETTLE_ATTRIBUTE, None) is None:
                continue
            checked += 1
            assert planet.radius == pytest.approx(_table_radius_km(planet), rel=giant_draw.SETTLE_TOLERANCE * 1.01)
    assert checked, "the sweep made no giant sized from the tables"
