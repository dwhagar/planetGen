# tests/test_circumbinary_light.py

"""GEN.1484: a planet of a close pair keeps the share of its light the two stars' eclipses leave, stored as `light_factor`; every other body keeps all of it."""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import climate_draws, orbit_elements, radiation_field
from planetgen.util import draw

from tests.fuzz_support import first_seed_with

PAIR_SEEDS = range(100000, 100200)
"""range: The run seeds searched for a close pair with planets (TEST.714: no seed is pinned, since a change to what the seeds draw moves the case)."""


def _is_a_close_pair_with_planets(system):
    return any(hasattr(p, "planet_class") and p.light_factor is not None for p in system.primary.orbiters)


@pytest.fixture(scope="module")
def pair_system():
    _seed, system = first_seed_with(lambda seed: StarSystem(system_config=SystemConfig()), _is_a_close_pair_with_planets,
                                    seeds=PAIR_SEEDS, what="a close pair with planets")
    return system


def test_the_planets_of_a_close_pair_store_a_factor_below_one(pair_system):
    planets = [p for p in pair_system.primary.orbiters if hasattr(p, "planet_class")]
    assert planets and all(p.light_factor is not None and 0.99 < p.light_factor < 1.0 for p in planets)
    assert all(moon.light_factor is None for p in planets for moon in p.moons)


def test_the_temperature_of_a_circumbinary_planet_follows_the_light_it_keeps(pair_system):
    checked = 0
    for planet in pair_system.primary.orbiters:
        if getattr(planet, "body_type", None) != "t" or planet.bond_albedo is None:
            continue
        flux = radiation_field.starlight_w_m2(planet.star, orbit_elements.flux_distance_au(float(planet.distance), planet.eccentricity)) * planet.light_factor
        assert planet.surface_temperature == pytest.approx(climate_draws.surface_temperature_k(planet.bond_albedo, planet.greenhouse_factor, flux), rel=1e-6)
        checked += 1
    assert checked > 0


def test_a_lone_star_stores_no_factor():
    draw.set_run_seed(100000)
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    system = StarSystem(system_config=config)
    for planet in system.primary.orbiters:
        if hasattr(planet, "planet_class"):
            assert planet.light_factor is None
