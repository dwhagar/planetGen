# tests/test_migration_stage_epochs.py

"""
GEN.764: migration 0164 gives a galaxy made at v163 the epochs it can be given: each kind of scatter that has rows gets the
plan's epoch, and a species gets the epoch of the sector its homeworld's system is in. `epoch_record` then counts the sectors and
species by epoch.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy import epoch_record
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.test_bright_star_scatter import _seed_galaxy, light_steps  # noqa: F401 -- a fixture


pytestmark = pytest.mark.fresh_database


def _system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if system.planets:
            return system, cfg
    pytest.fail("no system with a planet in 60 seeds")


def _make_species(conn, system_id, name="Tester"):
    planet_id = conn.execute("SELECT id FROM planets WHERE star_system_id = ? ORDER BY id LIMIT 1", (system_id,)).fetchone()["id"]
    conn.execute("INSERT INTO species (name, homeworld_planet_id, star_system_id, life_stage, build, climate, size) "
                 "VALUES (?, ?, ?, 'technological_civilization', 'stocky', 'temperate', 'medium')", (name, planet_id, system_id))


def test_a_stale_stage_is_found_by_its_own_epochs(mysql_config):
    """Fills made at 94 are behind epoch 96 (it changed the fill), but a scatter made at 94 is not (it did not)."""
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
            sector_id = conn.execute("INSERT INTO sectors (name, edge_mpc, generator_epoch) VALUES ('Old sector', 1000, 94)").lastrowid
            conn.execute("UPDATE star_systems SET sector_id = ? WHERE id = ?", (sector_id, system_id))
            _make_species(conn, system_id)
            conn.execute("UPDATE species SET generator_epoch = 94")
            conn.execute("INSERT INTO scatter_epochs (kind, generator_epoch) VALUES ('mass', 94)")
        found = epoch_record.behind(conn, current=96)
        assert found["stale"] == {"scatter": {}, "fill": {94: 1}, "population": {}}
        assert found["stages"] == ("fill", "population")       # a refilled sector's species are made again
        conn.execute("UPDATE scatter_epochs SET generator_epoch = 58 WHERE kind = 'mass'")
        conn.commit()
        found = epoch_record.behind(conn, current=96)
        assert found["stale"]["scatter"] == {"mass": 58} and found["stages"] == ("scatter", "fill", "population")
    finally:
        conn.rollback()
        conn.close()
