# tests/test_migration_radiogenic_content.py

"""
GEN.459: migration 0141 adds `c_k0_ppm`, `c_th0_ppm`, `c_u0_ppm`, `al26_ratio` and `urey_ratio` to `planets` and `moons`. A body from
before v141 keeps its rows with NULL there; the version is detected from `planets.c_k0_ppm`.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

COLUMNS = ("c_k0_ppm", "c_th0_ppm", "c_u0_ppm", "al26_ratio", "urey_ratio")


def test_a_v140_galaxy_keeps_its_planets_with_no_content(mysql_config):
    load_old_schema(mysql_config, 140)
    assert store.schema_status(mysql_config)[0] == 140
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not conn.execute("SELECT COUNT(*) AS n FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                                " AND TABLE_NAME = 'planets' AND COLUMN_NAME = 'c_k0_ppm'").fetchone()["n"]
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            found = {row["COLUMN_NAME"]: row["DATA_TYPE"] for row in conn.execute(
                "SELECT COLUMN_NAME, DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?",
                (table,)).fetchall()}
            assert all(found[name] == "double" for name in COLUMNS), table
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
