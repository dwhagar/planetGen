# tests/test_migration_temperature_coverage.py

"""
GEN.1265: what is read off the temperature of every place of a body is stored. A system saved and loaded keeps the columns; the migration gives a v204 galaxy the
columns, empty.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import temperature_coverage
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


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
        if any(getattr(body, "temp_livable_share", None) for body in _bodies(system)):
            return system, cfg
    pytest.fail("no system with a body that has a livable share in 200 seeds")


def test_a_saved_system_keeps_the_coverage_of_every_body(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before = [{column: getattr(body, column, None) for column in temperature_coverage.FIELDS} for body in _bodies(system)]
        after = [{column: getattr(body, column, None) for column in temperature_coverage.FIELDS} for body in _bodies(loaded)]
        assert after == before
        assert any(item["temp_livable_share"] for item in before)
    finally:
        conn.rollback()
        conn.close()


def test_a_v204_galaxy_gets_the_columns_empty(mysql_config):
    load_old_schema(mysql_config, 204)
    assert store.schema_status(mysql_config)[0] == 204
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert not conn.execute("SHOW COLUMNS FROM " + table + " LIKE 'temp_livable_share'").fetchall()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            for column in temperature_coverage.FIELDS:
                assert conn.execute("SHOW COLUMNS FROM " + table + " LIKE ?", (column,)).fetchall(), (table, column)
    finally:
        conn.close()
