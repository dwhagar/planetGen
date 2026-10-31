# tests/test_migration_ocean_tides.py

"""
GEN.1034: migration 0191 adds `tide_range_spring_m`, `tide_range_neap_m` and `ocean_liquid` to `planets` and `moons`. A body from before v190 keeps its
rows with NULL there.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

COLUMNS = {"tide_range_spring_m": "double", "tide_range_neap_m": "double", "ocean_liquid": "varchar"}


def _type(conn, table, column):
    row = conn.execute("SELECT DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ? AND COLUMN_NAME = ?",
                       (table, column)).fetchone()
    return row["DATA_TYPE"] if row else None


def test_a_v189_galaxy_gets_the_tide_columns_empty(mysql_config):
    load_old_schema(mysql_config, 189)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert all(_type(conn, "planets", column) is None for column in COLUMNS)
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            for column, kind in COLUMNS.items():
                assert _type(conn, table, column) == kind, (table, column)
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
