# tests/test_migration_cloud_decks.py

"""
GEN.1577: migration 0269 adds `planets.cloud_decks` (TEXT, NULL), a giant's cloud decks and trace-gas ceilings as JSON. The rows made before keep their figures and get NULL, and the migration can be run
over a table that already has the column (a stopped run is started again).
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _has_column(conn, table):
    return bool(conn.execute(f"SHOW COLUMNS FROM {table} LIKE 'cloud_decks'").fetchall())


def test_the_column_is_added_empty_to_planets_only(mysql_config):
    load_old_schema(mysql_config, 268)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not _has_column(conn, "planets")
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _has_column(conn, "planets")
        assert not _has_column(conn, "moons") and not _has_column(conn, "rogue_planets")
        assert conn.execute("SELECT COUNT(*) AS n FROM planets WHERE cloud_decks IS NOT NULL").fetchone()["n"] == 0
    finally:
        conn.close()
