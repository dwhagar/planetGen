# tests/test_migration_sector_density_used.py

"""
DB.60: migration 0125 adds `sectors.density_multiplier_used` (DECIMAL(10,4), NULL) and `sectors.override_id` (BIGINT UNSIGNED,
NULL) in one `ALGORITHM=COPY` rebuild. A sector saved before the revision has neither, a migration run again leaves the columns
as they are, and the columns hold what is written to them.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _column(conn, name):
    return conn.execute(
        "SELECT column_type AS t, is_nullable AS n FROM information_schema.COLUMNS"
        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'sectors' AND COLUMN_NAME = ?", (name,)).fetchone()


def test_a_v124_sector_has_no_density_after_the_migration(mysql_config):
    load_old_schema(mysql_config, 124)
    assert store.schema_status(mysql_config)[0] == 124
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _column(conn, "override_id") is None
        conn.execute("INSERT INTO sectors (name, edge_mpc) VALUES ('old', 4.0)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT density_multiplier_used, override_id FROM sectors WHERE name = 'old'").fetchone()
        assert row == {"density_multiplier_used": None, "override_id": None}
        assert _column(conn, "density_multiplier_used")["t"] == "decimal(10,4)"
        assert _column(conn, "override_id")["t"] == "bigint unsigned"
        assert {_column(conn, c)["n"] for c in ("density_multiplier_used", "override_id")} == {"YES"}
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_columns_as_they_are(mysql_config):
    load_old_schema(mysql_config, 124)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        override = conn.execute("INSERT INTO region_density_overrides (shape, params, multiplier, set_by_name)"
                                " VALUES ('sphere', '{}', 2.5, 'Boss')").lastrowid
        conn.execute("INSERT INTO sectors (name, edge_mpc, density_multiplier_used, override_id) VALUES ('kept', 4.0, 2.5, ?)",
                     (override,))
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0124')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT density_multiplier_used AS d, override_id AS o FROM sectors WHERE name = 'kept'").fetchone()
        assert float(row["d"]) == 2.5
        assert row["o"] == override
    finally:
        conn.close()
