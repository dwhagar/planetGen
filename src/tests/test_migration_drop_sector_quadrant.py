# tests/test_migration_drop_sector_quadrant.py

"""
DB.112: migration 0264 drops `sector_system_counts.quadrant` and `idx_sector_system_counts_quadrant`. The rows keep their count and density, and the migration can be run over a table that has
already lost them (a stopped run is started again).
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _has_column(conn):
    return bool(conn.execute("SHOW COLUMNS FROM sector_system_counts LIKE 'quadrant'").fetchall())


def _has_index(conn):
    return any(row["Key_name"] == "idx_sector_system_counts_quadrant" for row in conn.execute("SHOW INDEX FROM sector_system_counts").fetchall())


def test_the_column_and_its_index_go_and_the_rows_keep_their_figures(mysql_config):
    load_old_schema(mysql_config, 263)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _has_column(conn)
        assert _has_index(conn)
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        conn.execute("INSERT INTO sector_system_counts (sector_id, system_count, system_density, quadrant) VALUES (1, 3, 0.5, 'II')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not _has_column(conn)
        assert not _has_index(conn)
        row = conn.execute("SELECT system_count, system_density FROM sector_system_counts WHERE sector_id = 1").fetchone()
        assert (row["system_count"], row["system_density"]) == (3, 0.5)
    finally:
        conn.close()


def test_the_migration_can_be_run_over_a_table_that_has_lost_them(mysql_config):
    load_old_schema(mysql_config, 263)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("ALTER TABLE sector_system_counts DROP INDEX idx_sector_system_counts_quadrant")
        conn.execute("ALTER TABLE sector_system_counts DROP COLUMN quadrant")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
