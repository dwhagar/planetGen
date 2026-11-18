# tests/test_life_remakes_water.py

"""
GEN.1590: the water is made again when the life pass moves the temperature or the pressure of a world. A liquid ocean stored at 432 K under 63 kPa (water boils at 360 K there) was left as it was.
"""

import pytest

from planetgen.generation import life as planetLife
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import atmosphere, hydrosphere, planets
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
        if planet.hydrosphere == "surface ocean" and planet.build.surface_w:
            return planet
    pytest.fail("no ocean world among the first 60 seeds")


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
    planets.apply_class(planet)                   # the class is read once, from the finished body (GEN.1745), not by `_rebuild`
    assert planet.planet_class != "M"


def test_a_world_the_life_did_not_move_keeps_the_water_it_had(host_star):
    planet = _ocean_world(host_star)
    state, cover = planet.hydrosphere, planet.ocean_fraction
    planetLife._rebuild(planet, remake_water=False)
    assert (planet.hydrosphere, planet.ocean_fraction) == (state, cover)


def test_the_class_is_read_once_after_the_life_pass_and_never_during_it(host_star, monkeypatch):
    """GEN.1745 (Boss, 11:25 AM EST 2026-11-15): no class is assigned until the evolution is finished. `apply_life_data` makes the air, the water and the score, then reads the class once from the finished body; the rebuilds in between leave it unread."""
    planet = _ocean_world(host_star)
    calls = []
    real = planets.apply_class
    monkeypatch.setattr(planets, "apply_class", lambda body: (calls.append((body.life_stage, body.planet_class)), real(body))[1])
    planetLife.apply_life_data(planet)
    assert len(calls) == 1
    assert planet.planet_class is not None
