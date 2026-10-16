# tests/test_migration_drop_galactic_speed.py

"""GEN.336: migration 0103 drops the scalar galactic speed (the movement vector of 0102 replaces it) and can be run again."""

import importlib

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

MIGRATION = importlib.import_module("planetgen.db.migrations.versions.0103_drop_galactic_speed")


def _has(conn, table, column):
    return bool(conn.execute(
        "SELECT COUNT(*) AS n FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = ? AND COLUMN_NAME = ?", (table, column)).fetchone()["n"])


def test_every_scalar_speed_column_is_gone_and_the_vector_stays(mysql_config):
    load_old_schema(mysql_config, 102)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert all(_has(conn, table, column) for table, column in MIGRATION.DROPS.items())
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table, column in MIGRATION.DROPS.items():
            assert not _has(conn, table, column), (table, column)
        for table in ("black_holes", "neutron_stars", "rogue_planets", "interstellar_comets", "asteroid_fields",
                      "nebulae", "supernova_remnants"):
            assert _has(conn, table, "velocity_x_kms"), table
        assert store.detect_schema_version(conn) == store.SCHEMA_VERSION
    finally:
        conn.close()


def test_a_migration_run_again_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 102)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0103')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
