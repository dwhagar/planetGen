"""
DB.108: migration 0168 adds the mass columns and fills them from the stored bodies: a sector's `total_mass_kg` (NULL where a body has
no stored mass), and the block and region rows above it.
"""

import pytest

from planetgen.db import region_stats, store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

HOLDING = (40, 20, 0)
KNOWN_NOTHING = (40, 20, 1)
EMPTY = (41, 21, 0)
SECTORS = (HOLDING, KNOWN_NOTHING, EMPTY)
COMET_KG = (3.0e12, 4.5e12)
UNKNOWN_KG = None


def _comet(conn, sector_id, mass_kg):
    conn.execute(
        "INSERT INTO interstellar_comets (sector_id, name, nucleus_diameter_km, velocity_kms, is_active, composition_summary,"
        " galactic_orbital_period_gy, galactic_orbital_phase_deg, galactic_min_update_interval_years, mass_kg)"
        " VALUES (?, 'Test', 1.0, 26.0, 0, 'water ice', 0.2, 10.0, 1000.0, ?)", (sector_id, mass_kg))


def _old_galaxy(mysql_config):
    load_old_schema(mysql_config, 167)
    assert store.schema_status(mysql_config)[0] == 167
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        ids = {}
        for n, (shell, band, wedge) in enumerate(SECTORS):
            conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES (?, 4000, ?, ?, ?)",
                         (f"Sector {n}", shell, band, wedge))
            ids[(shell, band, wedge)] = conn.execute(
                "SELECT id FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
                (shell, band, wedge)).fetchone()["id"]
            conn.execute(
                "INSERT INTO sector_stats (shell_index, band_index, wedge_index, expected_systems, actual_systems, actual_stars,"
                " mean_age_gy, total_luminosity_sol, bright_level_sol) VALUES (?, ?, ?, 1.5, 1, 2, 1.0, 10.0, 0)",
                (shell, band, wedge))
        _comet(conn, ids[HOLDING], COMET_KG[0])
        _comet(conn, ids[HOLDING], COMET_KG[1])
        _comet(conn, ids[KNOWN_NOTHING], UNKNOWN_KG)
        conn.commit()
    finally:
        conn.close()


def test_the_migration_fills_every_sector_block_and_region(mysql_config):
    _old_galaxy(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        def total(address):
            return conn.execute("SELECT total_mass_kg FROM sector_stats WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
                                address).fetchone()["total_mass_kg"]
        assert total(HOLDING) == sum(COMET_KG)
        assert total(KNOWN_NOTHING) is None                # a comet with no stored mass: not summed, not guessed
        assert total(EMPTY) == 0.0
        assert conn.execute("SELECT MAX(mass_sum_kg) AS m FROM galaxy_block_stats").fetchone()["m"] == sum(COMET_KG)
        regions = {tuple(r["block"])[1:]: r["mass_sum_kg"] for r in region_stats.regions(conn)}
        assert regions[tuple(region_stats.region_cell(HOLDING))[1:]] == sum(COMET_KG)
        assert store.schema_status(mysql_config)[0] == store.SCHEMA_VERSION
    finally:
        conn.close()
