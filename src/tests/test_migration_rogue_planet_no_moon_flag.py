# tests/test_migration_rogue_planet_no_moon_flag.py

"""
GEN.558: migration 0244 drops `rogue_planets.has_moons` from a v243 galaxy and keeps its rows; a second run changes nothing, and a fresh schema has no such column.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import insert_filled_row, load_old_schema

pytestmark = pytest.mark.fresh_database


def _columns(conn):
    return [row["Field"] for row in conn.execute("SHOW COLUMNS FROM rogue_planets").fetchall()]


def _v243_galaxy(mysql_config):
    load_old_schema(mysql_config, 243)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "has_moons" in _columns(conn)
        insert_filled_row(conn, "rogue_planets", name="Old dwarf", planet_type="g", mass_kg=3.0e28, has_internal_heat=1, has_moons=1)
        conn.commit()
    finally:
        conn.close()


def test_a_v243_galaxy_loses_the_flag_and_keeps_its_rogue_planets(mysql_config):
    _v243_galaxy(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "has_moons" not in _columns(conn)
        row = conn.execute("SELECT name, mass_kg FROM rogue_planets").fetchone()
        assert row["name"] == "Old dwarf" and row["mass_kg"] == pytest.approx(3.0e28)
    finally:
        conn.close()


def test_a_second_run_changes_nothing(mysql_config):
    _v243_galaxy(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
