# tests/test_migration_ring_age.py

"""
GEN.1057: migration 0192 adds `ring_age_yr` to `planets`, `moons` and `rogue_planets`, after `ring_ice_fraction`. A body from before v191 keeps its rows with NULL there.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

TABLES = ("planets", "moons", "rogue_planets")


def _position(conn, table, column):
    row = conn.execute("SELECT ORDINAL_POSITION AS p FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ? AND COLUMN_NAME = ?",
                       (table, column)).fetchone()
    return row["p"] if row else None


def test_a_v191_galaxy_gets_the_ring_age_column_empty_after_the_ice_share(mysql_config):
    load_old_schema(mysql_config, 190)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert all(_position(conn, table, "ring_age_yr") is None for table in TABLES)
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in TABLES:
            assert _position(conn, table, "ring_age_yr") == _position(conn, table, "ring_ice_fraction") + 1, table
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
