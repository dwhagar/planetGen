# tests/test_migration_halo_populations.py

"""
GEN.234: migration 0094 widens `bright_stars.population` to the halo and cluster populations, so a pre-placed bright
star of the stellar halo can be stored. Nothing else changes.
"""

import importlib

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

migration = importlib.import_module("planetgen.db.migrations.versions.0094_halo_populations")


def _insert_star(conn, population):
    conn.execute(
        "INSERT INTO bright_stars (shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc, position_z_mpc,"
        " population, star_type, yerkes_class, mass_kg, radius_km, temperature_k, luminosity_w, age_gy, initial_mass_sol,"
        " seed) VALUES (0, 0, 0, 0, 0, 0, ?, 'K0III', 'III', 1, 1, 4000, 1, 11, 1, 1)", (population,))


def test_a_v93_galaxy_refuses_a_halo_star_until_it_is_migrated(mysql_config):
    load_old_schema(mysql_config, 93)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        with pytest.raises(Exception, match="chk_bright_stars_population"):
            _insert_star(conn, "halo")
    finally:
        conn.close()
    assert store.schema_status(mysql_config)[0] == 93
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for population in migration.POPULATIONS:
            _insert_star(conn, population)
        with pytest.raises(Exception, match="chk_bright_stars_population"):
            _insert_star(conn, "disk")
    finally:
        conn.close()


def test_the_check_names_exactly_the_six_populations(mysql_config):
    load_old_schema(mysql_config, 93)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        clause = conn.execute(
            "SELECT check_clause AS c FROM information_schema.check_constraints WHERE constraint_schema = DATABASE()"
            " AND constraint_name = 'chk_bright_stars_population'").fetchone()["c"]
    finally:
        conn.close()
    assert all(name in clause for name in migration.POPULATIONS)
    assert clause.count("_utf8mb4") == len(migration.POPULATIONS)
