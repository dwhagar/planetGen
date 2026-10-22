# tests/test_migration_nebula_block_stats.py

"""
MAP.172: migration 0116 makes `nebula_block_stats`, `nebula_block_classes` and `nebula_block_state` and fills them from
the `nebula_field` rows and the stored nebulae the galaxy already holds. A migration run again leaves what it made.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

TABLES = ("nebula_block_stats", "nebula_block_classes", "nebula_block_state")

FIELD = [(1, 0, 0, 0, "M", 40.0, 810.0, 200.0, 50.0), (1, 0, 0, 1, "N", 12.0, 830.0, 190.0, 45.0)]


def _exists(conn, table):
    return conn.execute("SELECT COUNT(*) AS n FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
                        " AND TABLE_NAME = ?", (table,)).fetchone()["n"]


def _old_galaxy(mysql_config):
    load_old_schema(mysql_config, 115)
    assert store.schema_status(mysql_config)[0] == 115
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not any(_exists(conn, table) for table in TABLES)
        for row in FIELD:
            conn.execute("INSERT INTO nebula_field (cell_i, cell_j, cell_k, n, nebula_class, radius_ly, center_x_pc,"
                         " center_y_pc, center_z_pc) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", row)
        conn.execute(
            "INSERT INTO nebulae (name, nebula_class, nebula_type, radius_ly, composition, formation_cause,"
            " dominant_species, density_cm3, temperature_k, extinction_av, galactic_orbital_period_gy,"
            " galactic_orbital_phase_deg, galactic_min_update_interval_years, center_x_pc, center_y_pc, center_z_pc,"
            " galactic_radius_pc) VALUES ('Old', 'C', 'emission', 2.0, 'gas', 'cloud', 'H', 1.0, 9000.0, 0.5, 0.2,"
            " 10.0, 1000.0, 805.0, 195.0, 52.0, 830.0)")
        conn.commit()
    finally:
        conn.close()


def test_the_migration_fills_the_totals_from_the_field_and_the_stored_nebulae(mysql_config):
    _old_galaxy(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert all(_exists(conn, table) for table in TABLES)
        level_three = conn.execute("SELECT nebula_count, dominant_class, cover FROM nebula_block_stats WHERE m = 3"
                                   ).fetchall()
        assert sum(row["nebula_count"] for row in level_three) == 3
        assert all(0 < row["cover"] <= 1 for row in level_three)
        assert {row["dominant_class"] for row in level_three} <= {"M", "N", "C", None}
        top = conn.execute("SELECT nebula_count FROM nebula_block_stats WHERE m = 243").fetchall()
        assert sum(row["nebula_count"] for row in top) == 3
        assert conn.execute("SELECT nebulae_through_id AS n FROM nebula_block_state").fetchone()["n"] == 1
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_totals_as_they_are(mysql_config):
    _old_galaxy(mysql_config)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        before = conn.execute("SELECT * FROM nebula_block_stats ORDER BY m, shell_index, band_index, wedge_index"
                              ).fetchall()
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0115')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT * FROM nebula_block_stats ORDER BY m, shell_index, band_index, wedge_index"
                            ).fetchall() == before
    finally:
        conn.close()
