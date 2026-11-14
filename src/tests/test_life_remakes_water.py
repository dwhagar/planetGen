# tests/test_life_remakes_water.py

"""
GEN.1590: the water is made again when the life pass moves the temperature or the pressure of a world. A liquid ocean stored at 432 K under 63 kPa (water boils at 360 K there) was left as it was.
"""

import pytest

from planetgen.generation import life as planetLife
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import atmosphere, hydrosphere
from tests import class_support

LIQUID = ("surface ocean", "brine seep", "hycean", "ice-covered ocean")


@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


def _ocean_world(host_star):
    inner, outer = host_star.habitable_zone
    for seed in range(60):
        planet = class_support.draw_planet(host_star, (inner + outer) / 2.0, seed)
        if planet.hydrosphere == "surface ocean" and getattr(planet, "_surface_w", None):
            return planet
    pytest.fail("(TEST.737) no ocean world among the first 60 seeds")


def _heat(planet, kelvin, pascal):
    planet.surface_temperature, planet.atmospheric_pressure = kelvin, pascal


def test_a_world_the_life_heated_past_the_boiling_point_loses_its_ocean(host_star):
    planet = _ocean_world(host_star)
    _heat(planet, 432.0, 63000.0)
    assert hydrosphere.water_boiling_k(63000.0) < 432.0
    planetLife._rebuild(planet, remake_water=False)
    assert planet.hydrosphere in LIQUID          # the old rebuild kept the ocean it had drawn
    planetLife._rebuild(planet, remake_water=True)
    assert planet.hydrosphere not in LIQUID, planet.hydrosphere
    assert planet.ocean_fraction == 0.0 and planet.land_fraction == 1.0
    assert planet.planet_class != "M"


def test_a_world_the_life_did_not_move_keeps_the_water_it_had(host_star):
    planet = _ocean_world(host_star)
    state, cover = planet.hydrosphere, planet.ocean_fraction
    planetLife._rebuild(planet, remake_water=False)
    assert (planet.hydrosphere, planet.ocean_fraction) == (state, cover)
