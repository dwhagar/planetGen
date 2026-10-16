# tests/test_migration_cell_star_ceiling.py

"""
GEN.281: migration 0100 gives `galaxy_shape.cell_star_ceiling` the default 5,000 (from 50,000), moves a galaxy still at
the old default to it, and adds the CHECK that the ceiling is at least 1.
"""

import importlib

import pytest

from planetgen import tuning
from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

migration = importlib.import_module("planetgen.db.migrations.versions.0100_cell_star_ceiling_5000")

_INSERT = (
    "INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc, bulge_scale_radius_pc,"
    " bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude, arm_level, nuclear_amplitude,"
    " spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
    " expected_system_count_at_density_1, outer_shell_index, cell_star_ceiling)"
    " VALUES (1, 2600, 300, 1580, 3.11, 2, 0.26, 0.4, 1, 1, 2600, 0, 36.5, 4, 6.3, 100, ?)")


def _ceiling(conn):
    return conn.execute("SELECT cell_star_ceiling AS c FROM galaxy_shape").fetchone()["c"]


def test_the_default_is_the_tuning_and_5000(mysql_config):
    assert tuning.CELL_STAR_CEILING == 5000.0
    store.get_connection(mysql_config).close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        default = conn.execute(
            "SELECT COLUMN_DEFAULT AS d FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
            " AND TABLE_NAME = 'galaxy_shape' AND COLUMN_NAME = 'cell_star_ceiling'").fetchone()["d"]
        assert float(default) == 5000.0
        with pytest.raises(Exception, match="chk_galaxy_shape_cell_star_ceiling"):
            conn.execute(_INSERT, (0.5,))
    finally:
        conn.close()


@pytest.mark.parametrize("stored, moved", [(50000.0, 5000.0), (20000.0, 20000.0)])
def test_a_v99_galaxy_at_the_old_default_moves_and_one_with_its_own_ceiling_keeps_it(mysql_config, stored, moved):
    load_old_schema(mysql_config, 99)
    assert store.schema_status(mysql_config)[0] == 99
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(_INSERT, (stored,))
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _ceiling(conn) == moved
        with pytest.raises(Exception, match="chk_galaxy_shape_cell_star_ceiling"):
            conn.execute("UPDATE galaxy_shape SET cell_star_ceiling = 0")
    finally:
        conn.close()


def test_a_migration_run_again_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 99)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0099')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
