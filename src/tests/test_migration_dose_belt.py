# tests/test_migration_dose_belt.py

"""
GEN.983: migration 0195 gives a v186 galaxy the belt dose column, empty, on both bodies tables.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def test_a_v186_galaxy_gets_the_column_empty(mysql_config):
    load_old_schema(mysql_config, 186)
    assert store.schema_status(mysql_config)[0] == 186
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert not conn.execute("SHOW COLUMNS FROM " + table + " LIKE 'dose_belt_msv_yr'").fetchall()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert conn.execute("SHOW COLUMNS FROM " + table + " LIKE 'dose_belt_msv_yr'").fetchall(), table
            assert conn.execute("SELECT COUNT(*) AS n FROM " + table + " WHERE dose_belt_msv_yr IS NOT NULL").fetchone()["n"] == 0
    finally:
        conn.close()
