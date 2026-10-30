# tests/test_migration_temperature_swings.py

"""
GEN.848: every place with a temperature has a day and night swing of its own. A system saved and loaded keeps them; migration 0178 gives
a v177 galaxy the seven columns, empty.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import regional_temperature
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

SWINGS = tuple(column for column in regional_temperature.SWING_OF.values() if column != "temp_swing_k")


def _bodies(system):
    return [body for planet in system.planets if hasattr(planet, "moons") for body in [planet, *planet.moons]]


def _system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(200):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        bodies = _bodies(system)
        if (any(getattr(b, "temp_polar_swing_k", None) for b in bodies)
                and any(getattr(b, "temp_day_side_swing_k", None) == 0.0 for b in bodies)):
            return system, cfg
    pytest.fail("no system with a turning world and a locked one in 200 seeds")


def test_a_saved_system_keeps_the_swing_of_every_place(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before = [{column: getattr(body, column, None) for column in SWINGS} for body in _bodies(system)]
        after = [{column: getattr(body, column, None) for column in SWINGS} for body in _bodies(loaded)]
        assert after == before
        assert any(item["temp_polar_swing_k"] for item in before) and any(item["temp_day_side_swing_k"] == 0.0 for item in before)
    finally:
        conn.rollback()
        conn.close()


def test_a_v177_galaxy_gets_the_columns_empty(mysql_config):
    load_old_schema(mysql_config, 177)
    assert store.schema_status(mysql_config)[0] == 177
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert not conn.execute("SHOW COLUMNS FROM " + table + " LIKE 'temp_polar_swing_k'").fetchall()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            for column in SWINGS:
                assert conn.execute("SHOW COLUMNS FROM " + table + " LIKE ?", (column,)).fetchall(), (table, column)
    finally:
        conn.close()
