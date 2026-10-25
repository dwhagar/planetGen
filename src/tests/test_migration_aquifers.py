# tests/test_migration_aquifers.py

"""
GEN.654, DB.92: migration 0157 adds `aquifer_depth_km` and `aquifer_thickness_km` to `planets` and `moons`; a body saved before keeps NULL.
"""

import pytest

from planetgen.db import store
from planetgen.physics import aquifers

pytestmark = pytest.mark.fresh_database


def test_a_v156_galaxy_keeps_its_planets_with_no_aquifer(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        for table in ("planets", "moons"):
            drops = ", ".join(f"DROP COLUMN {column}" for column in aquifers.AQUIFER_FIELDS)
            conn.execute(f"ALTER TABLE {table} {drops}, ALGORITHM=COPY")
        conn.execute("UPDATE alembic_version SET version_num = '0156'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (156)")
        conn.commit()
        assert store._stored_version(conn, "schema_migrations") == 156
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            for column in aquifers.AQUIFER_FIELDS:
                row = conn.execute(
                    "SELECT DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?"
                    " AND COLUMN_NAME = ?", (table, column)).fetchone()
                assert row["DATA_TYPE"] == "double", (table, column)
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
