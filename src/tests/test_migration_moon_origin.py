# tests/test_migration_moon_origin.py

"""
DB.70: a moon's origin and its giant parent's circumplanetary disk are stored. A system saved and loaded keeps both; migration gives the galaxy of the schema before it
galaxy the five columns, empty, and a moon without them reads back as having no origin.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import moon_composition
from planetgen.util import draw
from tests.db_schema_support import load_old_schema, migration_number

pytestmark = pytest.mark.fresh_database


def _system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(120):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(p.moons and p.moon_disk_lifetime_myr is not None for p in _planets(system)):
            return system, cfg
    pytest.fail("no system with a giant that has moons in 120 seeds")


def _planets(system):
    return [planet for planet in system.planets if hasattr(planet, "moons")]


def _origins(system):
    return [{field: getattr(body, field) for field in moon_composition.BODY_FIELDS}
            for planet in _planets(system) for body in [planet, *planet.moons]]


def test_a_saved_system_keeps_every_moons_origin_and_its_giants_disk(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before, after = _origins(system), _origins(loaded)
        assert any(item["moon_origin"] in moon_composition.ORIGINS for item in before)
        assert any(item["moon_disk_ice_line_km"] and item["moon_disk_ice_line_km"] > 0.0 for item in before)
        assert after == before
        for item in before:
            if item["moon_origin"] is not None:      # only a moon has an origin, and only a giant has a disk
                assert item["moon_disk_lifetime_myr"] is None
    finally:
        conn.rollback()
        conn.close()


def test_the_previous_schemas_galaxy_gets_the_columns_empty(mysql_config):
    load_old_schema(mysql_config, migration_number("moon_origin_and_disk") - 1)
    assert store.schema_status(mysql_config)[0] == migration_number("moon_origin_and_disk") - 1
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert not conn.execute("SHOW COLUMNS FROM " + table + " LIKE 'moon_origin'").fetchall()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            for column in moon_composition.BODY_FIELDS:
                assert conn.execute("SHOW COLUMNS FROM " + table + " LIKE ?", (column,)).fetchall(), (table, column)
    finally:
        conn.close()
