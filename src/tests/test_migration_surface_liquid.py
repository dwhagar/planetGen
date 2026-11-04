# tests/test_migration_surface_liquid.py

"""
GEN.1035: migration 0188 gives a v187 galaxy the three surface liquid columns, empty, on both bodies tables.
"""

import pytest

from planetgen.db import store
from planetgen.physics import surface_liquid_world
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def test_a_v187_galaxy_gets_the_columns_empty(mysql_config):
    load_old_schema(mysql_config, 187)
    assert store.schema_status(mysql_config)[0] == 187
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert not conn.execute("SHOW COLUMNS FROM " + table + " LIKE 'surface_liquid'").fetchall()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            for column in surface_liquid_world.COLUMN_FIELDS:
                assert conn.execute("SHOW COLUMNS FROM " + table + " LIKE ?", (column,)).fetchall(), (table, column)
    finally:
        conn.close()
