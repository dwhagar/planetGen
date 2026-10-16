# tests/test_migration_star_metallicity.py

"""
GEN.158: migration 0095 adds `metallicity_fe_h` to `stars` and `bright_stars`, NULL on the rows already there.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _columns(conn, table):
    return {row["COLUMN_NAME"]: row for row in conn.execute(
        "SELECT COLUMN_NAME, DATA_TYPE, IS_NULLABLE FROM information_schema.COLUMNS"
        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,)).fetchall()}


def test_a_v94_galaxy_gets_a_nullable_metallicity_column_on_both_tables(mysql_config):
    load_old_schema(mysql_config, 94)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(
            "INSERT INTO bright_stars (shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc,"
            " position_z_mpc, population, star_type, yerkes_class, mass_kg, radius_km, temperature_k, luminosity_w,"
            " age_gy, initial_mass_sol, seed) VALUES (0, 0, 0, 0, 0, 0, 'old', 'K0III', 'III', 1, 1, 4000, 1, 11, 1, 1)")
        conn.commit()
        for table in ("stars", "bright_stars"):
            assert "metallicity_fe_h" not in _columns(conn, table)
    finally:
        conn.close()
    assert store.schema_status(mysql_config)[0] == 94
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("stars", "bright_stars"):
            column = _columns(conn, table)["metallicity_fe_h"]
            assert (column["DATA_TYPE"], column["IS_NULLABLE"]) == ("double", "YES")
        assert conn.execute("SELECT metallicity_fe_h AS m FROM bright_stars").fetchone()["m"] is None
    finally:
        conn.close()


def test_a_stopped_migration_starts_again(mysql_config):
    """The first table done and the second not: the run goes on from there."""
    load_old_schema(mysql_config, 94)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("ALTER TABLE stars ADD COLUMN metallicity_fe_h DOUBLE NULL")
        conn.commit()
    finally:
        conn.close()
    assert store.schema_status(mysql_config)[0] == 94  # the marker is the last step's
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("stars", "bright_stars"):
            assert "metallicity_fe_h" in _columns(conn, table)
    finally:
        conn.close()
