# tests/test_migration_cluster_parent_link.py

"""
DB.61: migration 0127 makes `galaxy_components.name` nullable and adds `cluster_id` (BIGINT UNSIGNED, NULL, a key and a foreign
key to `galaxy_components.id`) to `stars` and `bright_stars`. Existing clusters keep their names, a star has no link until one is
written, the unique key stays over the named rows, a removed cluster clears the link, and a migration run again changes nothing.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

COMPONENT = ("INSERT INTO galaxy_components (kind, name, source, center_x_pc, center_y_pc, center_z_pc, velocity_x_kms,"
             " velocity_y_kms, velocity_z_kms, epoch_unix, mass_solar, half_mass_radius_pc, profile, age_gy, metallicity_fe_h,"
             " seed) VALUES ('cluster', ?, 'generated', 0, 0, 0, 0, 0, 0, 0, 1e5, 5, 'king', 10, -1, 1)")


def _column(conn, table, name):
    return conn.execute(
        "SELECT column_type AS t, is_nullable AS n FROM information_schema.COLUMNS"
        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ? AND COLUMN_NAME = ?", (table, name)).fetchone()


def test_a_v126_galaxy_keeps_its_cluster_names_and_stars_have_no_link(mysql_config):
    load_old_schema(mysql_config, 126)
    assert store.schema_status(mysql_config)[0] == 126
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _column(conn, "stars", "cluster_id") is None
        assert _column(conn, "galaxy_components", "name")["n"] == "NO"
        conn.execute(COMPONENT, ("NGC 104",))
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT name FROM galaxy_components").fetchall() == [{"name": "NGC 104"}]
        assert _column(conn, "galaxy_components", "name")["n"] == "YES"
        for table in ("stars", "bright_stars"):
            assert _column(conn, table, "cluster_id") == {"t": "bigint unsigned", "n": "YES"}
    finally:
        conn.close()


def test_unnamed_clusters_do_not_collide_but_a_name_stays_unique(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(COMPONENT, (None,))
        conn.execute(COMPONENT, (None,))
        conn.execute(COMPONENT, ("M 13",))
        with pytest.raises(Exception, match="uq_galaxy_components_name"):
            conn.execute(COMPONENT, ("M 13",))
    finally:
        conn.close()


def test_a_removed_cluster_clears_the_star_link(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        keys = {row["k"] for row in conn.execute(
            "SELECT CONSTRAINT_NAME AS k FROM information_schema.TABLE_CONSTRAINTS WHERE TABLE_SCHEMA = DATABASE()"
            " AND CONSTRAINT_TYPE = 'FOREIGN KEY'").fetchall()}
        assert {"fk_stars_cluster", "fk_bright_stars_cluster"} <= keys
        rules = {row["t"]: row["r"] for row in conn.execute(
            "SELECT TABLE_NAME AS t, DELETE_RULE AS r FROM information_schema.REFERENTIAL_CONSTRAINTS"
            " WHERE CONSTRAINT_SCHEMA = DATABASE() AND REFERENCED_TABLE_NAME = 'galaxy_components'"
            " AND TABLE_NAME IN ('stars', 'bright_stars')").fetchall()}
        assert rules == {"stars": "SET NULL", "bright_stars": "SET NULL"}
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_columns_as_they_are(mysql_config):
    load_old_schema(mysql_config, 126)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(COMPONENT, (None,))
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0126')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT name FROM galaxy_components").fetchall() == [{"name": None}]
    finally:
        conn.close()
