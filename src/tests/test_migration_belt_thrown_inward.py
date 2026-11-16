# tests/test_migration_belt_thrown_inward.py

"""
GEN.1231: migration 0265 adds `asteroid_belts.thrown_inward_kg` (DOUBLE, NULL), the mass of the belt the star's mass loss threw inward. The rows made before keep their figures and get NULL, and the
migration can be run over a table that already has the column (a stopped run is started again).
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _has_column(conn):
    return bool(conn.execute("SHOW COLUMNS FROM asteroid_belts LIKE 'thrown_inward_kg'").fetchall())


def test_the_column_is_added_empty_and_the_rows_keep_their_figures(mysql_config):
    load_old_schema(mysql_config, 264)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not _has_column(conn)
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        conn.execute("INSERT INTO asteroid_belts (id, star_system_id, orbital_index, distance_km, lower_limit_km, upper_limit_km, density, composition_summary, parent_type, parent_id) "
                     "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (1, 1, 1, 1.0, 0.5, 1.5, "typical", "rock", "star", 1))
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _has_column(conn)
        row = conn.execute("SELECT distance_km, thrown_inward_kg FROM asteroid_belts WHERE id = 1").fetchone()
        assert row["distance_km"] == 1.0 and row["thrown_inward_kg"] is None
    finally:
        conn.close()


def test_the_migration_can_be_run_over_a_table_that_has_the_column(mysql_config):
    load_old_schema(mysql_config, 264)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("ALTER TABLE asteroid_belts ADD COLUMN thrown_inward_kg DOUBLE NULL")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _has_column(conn)
    finally:
        conn.close()
