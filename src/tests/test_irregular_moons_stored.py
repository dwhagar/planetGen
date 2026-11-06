# tests/test_irregular_moons_stored.py

"""
GEN.1084: a planet keeps its irregular moons in the database and reads them back as generated, and the moons that no longer stay in its zone are dropped when it moves.
"""

import pytest

from planetgen.db import store
from planetgen.generation import irregular_moons
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

from tests.test_irregular_moons import JUPITER_KG, _planet, _source


def _jupiter_set():
    return irregular_moons.for_planet(_planet(5.2, JUPITER_KG, 69911.0), _source())


def _approx_dict(value):
    if isinstance(value, dict):
        return {key: _approx_dict(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_approx_dict(item) for item in value]
    return pytest.approx(value) if isinstance(value, float) else value


def test_a_set_round_trips_through_a_dict():
    held = _jupiter_set()
    assert held.moons
    again = irregular_moons.IrregularSet.from_dict(held.to_dict())
    assert again == held
    assert irregular_moons.IrregularSet.from_dict(None) is irregular_moons.NO_IRREGULARS


def _config():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.MOONS = False
    cfg.BINARY_SYSTEM = False
    return cfg


def test_a_saved_planet_reads_back_its_irregular_moons(mysql_config):
    config = _config()
    held = _jupiter_set()
    for seed in range(1, 30):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=config)
        planets = [body for body in system.planets if body.body_type != 'a' and not body.is_moon]
        if not planets:
            continue
        planets[0].irregulars = held
        saved = store.save_system(system, config, config=mysql_config)
        with store.get_connection(mysql_config) as conn:
            loaded = store.load_star_system(conn, saved)
        read = [body for body in loaded.planets if body.body_type != 'a' and not body.is_moon]
        assert _approx_dict(held.to_dict()) == read[0].irregulars.to_dict()
        assert all(body.irregulars is irregular_moons.NO_IRREGULARS for body in read[1:])
        return
    pytest.fail("no system with a planet in 29 tries")


def test_moving_a_planet_in_drops_the_irregular_moons_beyond_its_new_zone():
    star = type("Star", (), {"mass": 1.989e30})()
    planet = _planet(5.2, JUPITER_KG, 69911.0)
    planet.star = star
    planet.irregulars = _jupiter_set()
    before = len(planet.irregulars.moons)
    assert irregular_moons.keep_stable(planet) == 0 and len(planet.irregulars.moons) == before
    planet.distance = 0.3          # a much smaller Hill sphere
    dropped = irregular_moons.keep_stable(planet)
    assert dropped == before and not planet.irregulars.moons
    assert planet.irregulars.lost >= dropped

