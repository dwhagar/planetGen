"""
OPS.180: a placed sector that holds a star system or a recorded fill is charted whether or not the chart ledger holds its cell. The panels that count charted sectors
(the Generate page's "Galaxy at a glance", `shell_specs`; the Galaxy page's, `galaxy_placed_sectors`) read the ledger, so a sector with no row vanished from them. `chart_ledger.repair`
(migration 0253) writes the row and says what it found.
"""

import math

import pytest

from planetgen import tuning
from planetgen.db import chart_ledger, query, shell_specs, store
from planetgen.galaxy import spherical
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from tests.db_schema_support import insert_filled_row, load_old_schema

pytestmark = pytest.mark.db

EDGE_PC = float(tuning.DEFAULT_SECTOR_EDGE_PC)
CHARTED = (40, 1, 0)
FILLED = (40, 1, 1)           # a fill was recorded, the sector holds no system
ROGUE_ONLY = (41, 0, 2)       # a placeholder a rogue planet made: uncharted by design


def _position(address):
    center = spherical.sector_position_pc(*address, EDGE_PC)
    return {"center_x_pc": float(center[0]), "center_y_pc": float(center[1]), "center_z_pc": float(center[2]),
            "galactic_radius_pc": math.dist(center, (0.0, 0.0, 0.0)),
            "shell_index": address[0], "band_index": address[1], "wedge_index": address[2]}


def _save(config, address):
    sector = SpaceSector(f"Repair {address}", edge_ly=20.0)
    sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    return store.save_sector(sector, config=config, galaxy_position=_position(address))


@pytest.fixture
def galaxy(mysql_config):
    _save(mysql_config, CHARTED)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("INSERT INTO galaxy_shell (shell_index, high_band_index, cell_count) VALUES (40, 3, 50), (41, 3, 50)")
            store.uncharted_sector_id(conn, ROGUE_ONLY, EDGE_PC)
    finally:
        conn.close()
    return mysql_config


def _system(conn, sector_id):
    conn.execute("SET FOREIGN_KEY_CHECKS = 0")      # the system points at a configuration this test does not build
    insert_filled_row(conn, "star_systems", sector_id=sector_id, name="Held")
    conn.execute("SET FOREIGN_KEY_CHECKS = 1")
    conn.commit()


def _lose_ledger(conn, address):
    conn.execute("DELETE FROM chart_ledger WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", address)
    conn.commit()


def test_a_saved_sector_is_counted_as_charted_by_the_panels(galaxy):
    conn = store.get_connection(galaxy)
    try:
        specs = shell_specs.galaxy_shell_specs(conn, EDGE_PC)
        assert specs["charted"] == 1
        assert [row["id"] for row in query.galaxy_placed_sectors(conn)] == [conn.execute("SELECT id FROM sectors WHERE shell_index = 40").fetchone()["id"]]
    finally:
        conn.close()


def test_a_sector_with_content_and_no_ledger_row_is_found_and_repaired(galaxy):
    conn = store.get_connection(galaxy)
    try:
        sector_id = conn.execute("SELECT id FROM sectors WHERE shell_index = 40 AND band_index = 1 AND wedge_index = 0").fetchone()["id"]
        _system(conn, sector_id)
        _lose_ledger(conn, CHARTED)
        assert shell_specs.galaxy_shell_specs(conn, EDGE_PC)["charted"] == 0 and query.galaxy_placed_sectors(conn) == []      # the bug: the panels read none
        found = chart_ledger.diagnosis(conn)
        assert (found["placed"], found["charted"], found["missing"]) == (2, 0, 1)
        done = chart_ledger.repair(conn)
        assert done["repaired"] == 1 and done["missing"] == 1
        assert shell_specs.galaxy_shell_specs(conn, EDGE_PC)["charted"] == 1 and len(query.galaxy_placed_sectors(conn)) == 1
        assert chart_ledger.state(conn, CHARTED, 1) == chart_ledger.CHARTED
        again = chart_ledger.repair(conn)
        assert again["repaired"] == 0 and again["missing"] == 0
    finally:
        conn.close()


def test_a_recorded_fill_is_charted_and_a_rogue_planet_alone_is_not(galaxy):
    conn = store.get_connection(galaxy)
    try:
        store.uncharted_sector_id(conn, FILLED, EDGE_PC)
        conn.execute("INSERT INTO sector_stats (shell_index, band_index, wedge_index, bright_level_sol, level_before_fill_sol, relative_density, expected_systems)"
                     " VALUES (?, ?, ?, 0, -1, 1.0, 1.0)", FILLED)
        conn.commit()
        sector_id = conn.execute("SELECT id FROM sectors WHERE shell_index = 41").fetchone()["id"]
        insert_filled_row(conn, "rogue_planets", sector_id=sector_id, name="Wanderer")
        conn.commit()
        done = chart_ledger.repair(conn)
        assert done["repaired"] == 1
        assert chart_ledger.state(conn, FILLED, 1) == chart_ledger.CHARTED and chart_ledger.state(conn, ROGUE_ONLY, 1) == chart_ledger.UNCHARTED
    finally:
        conn.close()


@pytest.mark.fresh_database
def test_the_migration_writes_the_missing_rows_and_says_what_it_found(mysql_config, caplog):
    import logging

    load_old_schema(mysql_config, 246)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("INSERT INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc, shell_index, band_index, wedge_index)"
                     " VALUES ('Old', 4000, 1.0, 2.0, 3.0, 3.7, 40, 1, 0)")
        sector_id = conn.execute("SELECT id FROM sectors WHERE shell_index = 40").fetchone()["id"]
        conn.commit()
        _system(conn, sector_id)
    finally:
        conn.close()
    with caplog.at_level(logging.INFO):
        assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert chart_ledger.state(conn, (40, 1, 0), 1) == chart_ledger.CHARTED
    finally:
        conn.close()
    assert "1 placed sectors" in caplog.text and "1 held content with no row" in caplog.text
