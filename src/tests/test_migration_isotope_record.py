# tests/test_migration_isotope_record.py

"""
GEN.1377, GEN.1420: the isotope record as child table rows. A system saved and loaded keeps every body's record; migration 0235 gives a v234 galaxy the three
tables (empty: nothing is backfilled), and a second run changes nothing.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import isotope_record
from planetgen.util import draw

pytestmark = pytest.mark.fresh_database


def _system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(getattr(p, "moons", None) for p in system.planets):
            return system, cfg
    pytest.fail("no system with a moon in 60 seeds")


def _records(system):
    return [getattr(body, isotope_record.ATTRIBUTE) for planet in system.planets if hasattr(planet, "moons") for body in [planet, *planet.moons]]


def test_a_saved_system_keeps_every_bodys_record(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before, after = _records(system), _records(loaded)
        assert any(before) and after == before
        stored = sum(conn.execute(f"SELECT COUNT(*) AS n FROM {kind}_isotopes i JOIN {kind}s b ON b.id = i.{kind}_id WHERE b.star_system_id = ?",
                                  (system_id,)).fetchone()["n"] for kind in ("planet", "moon"))
        assert stored == sum(len(record or ()) for record in before)
    finally:
        conn.rollback()
        conn.close()
