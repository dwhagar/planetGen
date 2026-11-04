# tests/test_belt_origin_stored.py

"""
GEN.1069: a belt keeps its origin, its history and its largest body in the database, and reads them back as generated.
"""

import pytest

from planetgen.db import store
from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

FIELDS = ("origin", "co_causes", "formation_snow_line_au", "primordial_mass_earth", "cleared_mass_earth", "scattered_mass_earth", "stirred", "origin_history",
          "origin_sentence", "scatter_fates", "largest_body_km")


def _config():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.MOONS = False
    cfg.BINARY_SYSTEM = False
    return cfg


def test_a_belt_round_trips_its_origin_and_largest_body_through_a_dict():
    config = _config()
    with deterministic_entropy(1069):
        for seed in range(1, 40):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=config)
            belt = next((b for b in system.planets if b.body_type == 'a' and b.origin is not None), None)
            if belt is None:
                continue
            again = AsteroidBelt.from_dict(belt.to_dict(), config)
            assert belt.largest_body_km > 0.0
            for name in FIELDS:
                assert getattr(again, name) == getattr(belt, name), name
            return
    pytest.fail("no belt with an origin in 39 systems")


def test_a_saved_belt_reads_back_what_it_was_made_with(mysql_config):
    config = _config()
    seen = 0
    for seed in range(1, 25):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=config)
        belts = [b for b in system.planets if b.body_type == 'a' and b.origin is not None]
        if not belts:
            continue
        saved = store.save_system(system, config, config=mysql_config)
        with store.get_connection(mysql_config) as conn:
            loaded = store.load_star_system(conn, saved)
        read = [b for b in loaded.planets if b.body_type == 'a' and b.origin is not None]
        assert len(read) == len(belts)
        for before, after in zip(belts, read):
            for name in FIELDS:
                value, got = getattr(before, name), getattr(after, name)
                if isinstance(value, float):
                    assert got == pytest.approx(value), name
                elif value:
                    assert got == value, name
                else:
                    assert not got, name
            assert after.origin == before.origin and list(after.origin_history) == list(before.origin_history)
            seen += 1
        if seen >= 3:
            return
    assert seen, "no belt with an origin was saved"
