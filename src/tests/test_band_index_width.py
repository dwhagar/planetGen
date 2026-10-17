# tests/test_band_index_width.py

"""
DB.72: a band number is INT wherever the scatter keeps it. A shell's bands run from `-(h + 1)` to `h + 1`, which passes
SMALLINT's 32,767 at shell 20,861 and reaches 102,942 at the outermost shell of the 262 kpc reach, so the whole-sphere
scatter stopped at 97 percent with "Out of range value for column 'band_index'" while the tests, which stay in the first
few hundred shells, never saw it.
"""

import importlib

import pytest

from planetgen.db import store
from planetgen.galaxy import spherical
from tests.db_schema_support import load_old_schema
from tests.test_migration_halo_populations import _insert_star

OUTERMOST_SHELL = 65535
"""int: The last shell of the 262 kpc reach (65,536 shells of 4 pc)."""

migration = importlib.import_module("planetgen.db.migrations.versions.0108_scatter_band_int")


def _largest_band(shell):
    return spherical.band_range(shell)[1]


def _store_one_row_in_each_table(conn, shell, band):
    conn.execute(
        "INSERT INTO bright_stars (shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc, position_z_mpc,"
        " population, star_type, yerkes_class, mass_kg, radius_km, temperature_k, luminosity_w, age_gy, initial_mass_sol,"
        " seed) VALUES (?, ?, 0, 0, 0, 0, 'young', 'O5V', 'V', 1, 1, 40000, 1, 0.001, 40, 1)", (shell, band))
    conn.execute("INSERT INTO sector_stats (shell_index, band_index, wedge_index) VALUES (?, ?, 0)", (shell, band))
    conn.execute(
        "INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, position_x_mpc, position_y_mpc,"
        " position_z_mpc, seed) VALUES (?, ?, 0, 'black-hole', 0, 0, 0, 1)", (shell, band))
    conn.commit()


def _read_back(conn, shell, band):
    return {table: conn.execute(f"SELECT COUNT(*) AS n FROM {table} WHERE shell_index = ? AND band_index = ?",
                                (shell, band)).fetchone()["n"]
            for table in migration.TABLES}


def test_the_outermost_shell_keeps_its_polar_cap_band_in_every_scatter_table(mysql_config):
    band = _largest_band(OUTERMOST_SHELL)
    assert band > 2 ** 16, "the outermost shell's cap band is past SMALLINT, even unsigned"
    conn = store.get_connection(mysql_config)
    try:
        for signed in (band, -band):
            _store_one_row_in_each_table(conn, OUTERMOST_SHELL, signed)
            assert _read_back(conn, OUTERMOST_SHELL, signed) == {table: 1 for table in migration.TABLES}
    finally:
        conn.rollback()
        conn.execute("DELETE FROM bright_stars WHERE shell_index = ?", (OUTERMOST_SHELL,))
        conn.execute("DELETE FROM sector_stats WHERE shell_index = ?", (OUTERMOST_SHELL,))
        conn.execute("DELETE FROM phenomenon_scatter WHERE shell_index = ?", (OUTERMOST_SHELL,))
        conn.commit()
        conn.close()


def test_the_first_shell_past_smallint_is_shell_20861():
    assert _largest_band(20860) == 32766
    assert _largest_band(20861) == 32768


def test_no_column_named_for_a_band_is_narrower_than_int(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        narrow = conn.execute(
            "SELECT table_name AS t, column_name AS c, column_type AS y FROM information_schema.columns"
            " WHERE table_schema = DATABASE() AND LOCATE('band', column_name) > 0"
            " AND data_type IN ('tinyint', 'smallint', 'mediumint')").fetchall()
        assert [(row["t"], row["c"], row["y"]) for row in narrow] == []
    finally:
        conn.close()


@pytest.mark.fresh_database
def test_a_v107_galaxy_takes_the_outermost_band_after_the_migration(mysql_config):
    band = _largest_band(OUTERMOST_SHELL)
    load_old_schema(mysql_config, 107)
    assert store.schema_status(mysql_config)[0] == 107
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert_star(conn, "young")        # a row at band 0 that the migration must keep as it is
        conn.commit()
        with pytest.raises(Exception, match="band_index"):
            _store_one_row_in_each_table(conn, OUTERMOST_SHELL, band)
        conn.rollback()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 1
        _store_one_row_in_each_table(conn, OUTERMOST_SHELL, band)
        assert _read_back(conn, OUTERMOST_SHELL, band) == {table: 1 for table in migration.TABLES}
    finally:
        conn.close()


@pytest.mark.fresh_database
def test_a_migration_run_again_leaves_the_columns_as_they_are(mysql_config):
    load_old_schema(mysql_config, 107)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0107')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in migration.TABLES:
            assert conn.execute(
                "SELECT column_type AS y FROM information_schema.columns WHERE table_schema = DATABASE()"
                " AND table_name = ? AND column_name = 'band_index'", (table,)).fetchone()["y"] == "int"
    finally:
        conn.close()
