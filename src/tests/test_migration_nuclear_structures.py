# tests/test_migration_nuclear_structures.py

"""
GEN.280: migration 0099 renames `galaxy_shape.core_amplitude` to `nuclear_amplitude` (every galaxy planned before it
had no nuclear structures, so 0) and widens `bright_stars.population` to the two nuclear populations.
"""

import importlib

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema
from tests.test_migration_halo_populations import _insert_star

pytestmark = pytest.mark.fresh_database

migration = importlib.import_module("planetgen.db.migrations.versions.0099_nuclear_structures")


def _columns(conn):
    return {row["COLUMN_NAME"] for row in conn.execute(
        "SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'galaxy_shape'").fetchall()}


def test_a_v98_galaxy_with_a_core_gets_no_nuclear_structures_and_takes_nuclear_stars(mysql_config):
    load_old_schema(mysql_config, 98)
    assert store.schema_status(mysql_config)[0] == 98
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(
            "INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc, bulge_scale_radius_pc,"
            " bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude, arm_level, core_amplitude,"
            " spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
            " expected_system_count_at_density_1, outer_shell_index)"
            " VALUES (1, 2600, 300, 1580, 3.11, 2, 0.26, 0.4, 1, 15, 2600, 0, 36.5, 4, 6.3, 100)")
        conn.commit()
        with pytest.raises(Exception, match="chk_bright_stars_population"):
            _insert_star(conn, "nuclear_disc")
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        columns = _columns(conn)
        assert "nuclear_amplitude" in columns and "core_amplitude" not in columns
        assert conn.execute("SELECT nuclear_amplitude FROM galaxy_shape").fetchone()["nuclear_amplitude"] == 0
        for population in migration.POPULATIONS:
            _insert_star(conn, population)
        with pytest.raises(Exception, match="chk_bright_stars_population"):
            _insert_star(conn, "disk")
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_column_and_the_check(mysql_config):
    load_old_schema(mysql_config, 98)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0098')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "nuclear_amplitude" in _columns(conn)
        clause = conn.execute(
            "SELECT check_clause AS c FROM information_schema.check_constraints WHERE constraint_schema = DATABASE()"
            " AND constraint_name = 'chk_bright_stars_population'").fetchone()["c"]
        assert clause.count("_utf8mb4") == len(migration.POPULATIONS)
    finally:
        conn.close()
