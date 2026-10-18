# tests/test_migration_life_structure.py

"""DB.50: migration 0110 adds the structure columns for life and flavour to `planets` and `moons`, touches no row,
and can be run again."""

import importlib

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema
from tests.test_migration_point_columns import _columns, _insert

pytestmark = pytest.mark.fresh_database

MIGRATION = importlib.import_module("planetgen.db.migrations.versions.0110_life_structure_columns")
NAMES = [name for name, _kind in MIGRATION.COLUMNS]


def test_a_v109_database_gains_the_columns_and_keeps_its_rows(mysql_config):
    load_old_schema(mysql_config, 109)
    assert store.schema_status(mysql_config)[0] == 109
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert not set(NAMES) & _columns(conn, table)
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        _insert(conn, "planets", body_type="t", star_system_id=7, star_id=4, flavor_text="a quiet shore",
                evolutionary_speed="Standard")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert set(NAMES) <= _columns(conn, table)
        row = conn.execute("SELECT evolutionary_speed, " + ", ".join(NAMES) + " FROM planets").fetchone()
        assert store.detect_schema_version(conn) == store.SCHEMA_VERSION
    finally:
        conn.close()
    assert row["evolutionary_speed"] == "Standard"
    assert all(row[name] is None for name in NAMES)  # nothing reads or writes them yet


def test_a_migration_run_again_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 109)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0109')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION


def test_the_columns_round_trip_a_value(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        _insert(conn, "planets", body_type="t", star_system_id=7, star_id=4, life_stage="complex_cells",
                evolutionary_pace="normal", life_milestone_age_gyr=2.5, life_system_age_gyr=3.1, flavor_key="planet:4")
        conn.commit()
        row = conn.execute("SELECT " + ", ".join(NAMES) + " FROM planets").fetchone()
    finally:
        conn.close()
    assert tuple(row[name] for name in NAMES) == ("complex_cells", "normal", 2.5, 3.1, "planet:4")
