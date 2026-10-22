# tests/test_migration_galaxy_family_columns.py

"""
DB.66: migration 0117 adds `galaxy_shape.family` and the component columns of docs/design/galaxy-types.md section 10. A
galaxy planned before the revision becomes a `spiral` with no spheroid and no clump share and keeps every other number,
the checks refuse a family or a number the design does not allow, and a migration run again leaves the columns as they are.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

COLUMNS = (
    "family", "bar_present", "bar_axis_ratio_b", "bar_axis_ratio_c", "bar_angle_rad", "thick_disk_amplitude",
    "spheroid_amplitude", "spheroid_sersic_n", "spheroid_effective_radius_pc", "spheroid_axis_ratio_b",
    "spheroid_axis_ratio_c", "spheroid_orientation_1_rad", "spheroid_orientation_2_rad", "spheroid_orientation_3_rad",
    "body_clump_share")

SHAPE = ("INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc, bulge_scale_radius_pc, bulge_amplitude,"
         " arm_count, pitch_angle_rad, arm_amplitude, spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
         " expected_system_count_at_density_1, outer_shell_index) VALUES (1, 2600, 300, 500, 1, 4, 0.2, 0.3, 3000, 0, 1, 100, 50, 7)")


def _columns(conn):
    return {row["c"]: row for row in conn.execute(
        "SELECT column_name AS c, column_type AS t, is_nullable AS n FROM information_schema.COLUMNS"
        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'galaxy_shape'").fetchall()}


def test_a_v116_galaxy_becomes_a_spiral_with_no_spheroid_and_keeps_its_numbers(mysql_config):
    load_old_schema(mysql_config, 116)
    assert store.schema_status(mysql_config)[0] == 116
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "family" not in _columns(conn)
        conn.execute(SHAPE)
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT * FROM galaxy_shape WHERE id = 1").fetchone()
        assert row["family"] == "spiral"
        for column in ("bar_present", "bar_axis_ratio_b", "bar_axis_ratio_c", "bar_angle_rad", "thick_disk_amplitude"):
            assert row[column] is None, column
        for column in ("spheroid_amplitude", "spheroid_sersic_n", "spheroid_effective_radius_pc", "spheroid_orientation_1_rad",
                       "spheroid_orientation_2_rad", "spheroid_orientation_3_rad", "body_clump_share"):
            assert row[column] == 0, column
        assert row["spheroid_axis_ratio_b"] == 1 and row["spheroid_axis_ratio_c"] == 1
        assert row["disk_scale_length_pc"] == 2600 and row["arm_count"] == 4 and row["k_norm"] == 1
        columns = _columns(conn)
        assert set(COLUMNS) <= set(columns)
        assert columns["family"]["n"] == "NO" and columns["bar_present"]["n"] == "YES"
        assert columns["spheroid_amplitude"]["n"] == "NO"
    finally:
        conn.close()


@pytest.mark.parametrize("assignment", [
    "family = 'ring'", "bar_axis_ratio_b = 0", "bar_axis_ratio_c = 1.5", "thick_disk_amplitude = -1",
    "spheroid_amplitude = -1", "spheroid_sersic_n = -1", "spheroid_effective_radius_pc = -1",
    "spheroid_axis_ratio_b = 0", "spheroid_axis_ratio_c = 2", "body_clump_share = 1.5", "body_clump_share = -0.1"])
def test_the_checks_refuse_a_family_or_a_number_the_design_does_not_allow(mysql_config, assignment):
    load_old_schema(mysql_config, 116)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(SHAPE)
        conn.commit()
        with pytest.raises(Exception, match="(?i)check"):
            conn.execute(f"UPDATE galaxy_shape SET {assignment} WHERE id = 1")
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_columns_as_they_are(mysql_config):
    load_old_schema(mysql_config, 116)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(SHAPE)
        conn.execute("UPDATE galaxy_shape SET family = 'elliptical', spheroid_amplitude = 2.5, spheroid_sersic_n = 4,"
                     " bar_present = 0, thick_disk_amplitude = 0 WHERE id = 1")
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0116')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT family, spheroid_amplitude, spheroid_sersic_n, bar_present, thick_disk_amplitude"
                           " FROM galaxy_shape WHERE id = 1").fetchone()
        assert (row["family"], row["spheroid_amplitude"], row["spheroid_sersic_n"]) == ("elliptical", 2.5, 4)
        assert row["bar_present"] == 0 and row["thick_disk_amplitude"] == 0
    finally:
        conn.close()
