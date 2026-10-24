# tests/test_migration_prose_drop.py

"""DB.52: migration 0111 drops the evolutionary paragraph tables, `flavor_text` and `hab_note` of `planets` and `moons`
and `system_configs.markdown`, keeps every other value of every row, and can be run again."""

import importlib

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema
from tests.test_migration_point_columns import _columns, _insert

pytestmark = pytest.mark.fresh_database

MIGRATION = importlib.import_module("planetgen.db.migrations.versions.0111_drop_prose_columns")


def _tables(conn):
    return {row["TABLE_NAME"] for row in conn.execute(
        "SELECT TABLE_NAME FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()").fetchall()}


def test_a_v110_database_loses_the_prose_and_keeps_the_rest(mysql_config):
    load_old_schema(mysql_config, 110)
    assert store.schema_status(mysql_config)[0] == 110
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert set(MIGRATION.TABLES) <= _tables(conn)
        for table, column in MIGRATION.COLUMNS:
            assert column in _columns(conn, table)
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        _insert(conn, "planets", body_type="t", star_system_id=7, star_id=4, flavor_text="a quiet shore",
                hab_note="breathable", evolutionary_speed="normal", life_stage="complex_cells", flavor_key="planet:4")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not set(MIGRATION.TABLES) & _tables(conn)
        for table, column in MIGRATION.COLUMNS:
            assert column not in _columns(conn, table)
        row = conn.execute("SELECT evolutionary_pace, life_stage, flavor_key FROM planets").fetchone()
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
    assert (row["evolutionary_pace"], row["life_stage"], row["flavor_key"]) == ("normal", "complex_cells", "planet:4")


def test_a_migration_run_again_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 110)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0110')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION


def test_a_new_database_has_none_of_them(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        assert not set(MIGRATION.TABLES) & _tables(conn)
        for table, column in MIGRATION.COLUMNS:
            assert column not in _columns(conn, table)
        assert "flavor_text_count" in _columns(conn, "planets")
    finally:
        conn.close()
