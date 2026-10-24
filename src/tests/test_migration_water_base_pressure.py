# tests/test_migration_water_base_pressure.py

"""
GEN.461: migration 0151 adds `water_base_pressure_pa` to `planets` and `moons`. A body from before v151 keeps its rows with NULL
there.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def test_a_v150_galaxy_keeps_its_planets_with_no_base_pressure(mysql_config):
    load_old_schema(mysql_config, 150)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not conn.execute("SELECT COUNT(*) AS n FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                                " AND TABLE_NAME = 'planets' AND COLUMN_NAME = 'water_base_pressure_pa'").fetchone()["n"]
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            row = conn.execute(
                "SELECT DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?"
                " AND COLUMN_NAME = 'water_base_pressure_pa'", (table,)).fetchone()
            assert row["DATA_TYPE"] == "double", table
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
