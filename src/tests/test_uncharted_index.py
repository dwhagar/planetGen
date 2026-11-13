# tests/test_uncharted_index.py

"""
OPS.154: a cell a scatter left a waiting bright star or an unbuilt scattered phenomenon in has no sector row of its own: the block light stage makes none, the
Sectors list reads the cell from the parts (`db.query.virtual_cells`), and migration 0247 deletes the placeholder rows MAP.402 made (`db.uncharted_index.clear_unoccupied`).
"""

import pytest

from planetgen.api.config import Config
from planetgen.db import chart_ledger, query, store, uncharted_index
from planetgen.generation import run_plan
from planetgen.web.app import create_app
from tests.db_schema_support import insert_filled_row, load_old_schema

pytestmark = pytest.mark.db

EDGE_PC = 4.0
STAR_CELL = (9, 2, 7)         # a waiting bright star
BUILT_STAR_CELL = (9, 2, 8)   # its star became a system
SCATTER_CELL = (10, 3, 4)     # an unbuilt scattered phenomenon
BUILT_SCATTER_CELL = (11, 1, 2)
BOTH_CELL = (12, 1, 3)        # a star and a phenomenon: one row
CHARTED = (9, 2, 5)           # charted, with a waiting star nobody built (it lies in a sector that exists)
ROGUE_CELL = (13, 1, 1)       # a placeholder a rogue planet made earlier, with a waiting star too


def _star(conn, address, system_id=None):
    conn.execute(
        "INSERT INTO bright_stars (shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc, position_z_mpc, population, star_type, yerkes_class,"
        " mass_kg, radius_km, temperature_k, luminosity_w, age_gy, initial_mass_sol, seed, star_system_id) VALUES (?, ?, ?, 1.0, 2.0, 3.0, 'young', 'B2', 'V',"
        " 3e31, 5e6, 2e4, 1e30, 0.1, 15.0, 9, ?)", (*address, system_id))


def _scatter(conn, address, built):
    conn.execute(
        "INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, position_x_mpc, position_y_mpc, position_z_mpc, seed, built_at)"
        f" VALUES (?, ?, ?, 'black-hole', 1.0, 2.0, 3.0, 7, {'CURRENT_TIMESTAMP' if built else 'NULL'})", address)


@pytest.fixture
def galaxy(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = conn.execute("SELECT COALESCE(MAX(id), 0) + 1 AS n FROM star_systems").fetchone()["n"]
            conn.execute("SET FOREIGN_KEY_CHECKS = 0")      # the star points at a system this test does not build
            _star(conn, STAR_CELL)
            _star(conn, BUILT_STAR_CELL, system_id=system_id)
            conn.execute("SET FOREIGN_KEY_CHECKS = 1")
            _scatter(conn, SCATTER_CELL, built=False)
            _scatter(conn, BUILT_SCATTER_CELL, built=True)
            _star(conn, BOTH_CELL)
            _scatter(conn, BOTH_CELL, built=False)
            _star(conn, CHARTED)
            charted_id = store.uncharted_sector_id(conn, CHARTED, EDGE_PC)
            chart_ledger.record(conn, CHARTED, store.sector_run_id(conn), 1)
            conn.execute("UPDATE sectors SET name = 'Charted' WHERE id = ?", (charted_id,))
            store.uncharted_sector_id(conn, ROGUE_CELL, EDGE_PC)
            _star(conn, ROGUE_CELL)
    finally:
        conn.close()
    return mysql_config


def _addresses(conn):
    return {(r["shell_index"], r["band_index"], r["wedge_index"]) for r in conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors").fetchall()}


def _placeholders(conn):
    """The rows an earlier run made (MAP.402) for the three waiting cells."""
    for cell in (STAR_CELL, SCATTER_CELL, BOTH_CELL):
        store.uncharted_sector_id(conn, cell, EDGE_PC)
    conn.commit()


def _listed(conn, **options):
    return {(r["shell_index"], r["band_index"], r["wedge_index"]) for r in query.list_sectors(conn, coverage="uncharted", show_empty=False, limit=50, **options)}


def test_the_block_light_stage_makes_no_sector_row(galaxy):
    conn = store.get_connection(galaxy)
    try:
        before = _addresses(conn)
        run_plan.make_block_lights(conn, lambda done, total: None, lambda text: None, rebuild=True)
        conn.commit()
        assert _addresses(conn) == before
        assert not hasattr(uncharted_index, "index")
    finally:
        conn.close()


def test_the_cells_are_listed_and_counted_from_the_parts_without_a_row(galaxy):
    conn = store.get_connection(galaxy)
    try:
        run_plan.make_block_lights(conn, lambda done, total: None, lambda text: None, rebuild=True)
        conn.commit()
        assert _listed(conn) == {STAR_CELL, SCATTER_CELL, BOTH_CELL, ROGUE_CELL}      # the rogue's row is the cell's own; the charted cell is not uncharted
        assert BUILT_STAR_CELL not in _listed(conn) and BUILT_SCATTER_CELL not in _listed(conn)
        assert query.count_uncharted_sectors(conn) == 4
    finally:
        conn.close()


def test_the_rows_an_earlier_run_made_give_the_same_list_as_the_parts_alone(galaxy):
    conn = store.get_connection(galaxy)
    try:
        run_plan.make_block_lights(conn, lambda done, total: None, lambda text: None, rebuild=True)
        conn.commit()
        without = _listed(conn)
        _placeholders(conn)
        assert _listed(conn) == without
        assert query.count_uncharted_sectors(conn) == 4
    finally:
        conn.close()


def test_a_row_a_rogue_made_is_one_row_per_cell(galaxy):
    conn = store.get_connection(galaxy)
    try:
        first = store.uncharted_sector_id(conn, STAR_CELL, EDGE_PC)
        assert store.uncharted_sector_id(conn, STAR_CELL, EDGE_PC) == first
        row = conn.execute("SELECT * FROM sectors WHERE id = ?", (first,)).fetchone()
        assert row["uid"] is not None and row["center_x_pc"] is not None and row["edge_mpc"] == EDGE_PC * 1000.0 and query.sector_is_uncharted(conn, row)
    finally:
        conn.close()


def test_the_api_lists_the_cells_only_when_asked(galaxy):
    conn = store.get_connection(galaxy)
    try:
        run_plan.make_block_lights(conn, lambda done, total: None, lambda text: None, rebuild=True)
        conn.commit()
    finally:
        conn.close()

    class TestConfig(Config):
        MYSQL_CONFIG = galaxy
        WRITE_MYSQL_CONFIG = galaxy
        CONTROL_MYSQL_CONFIG = galaxy

    app = create_app(TestConfig)
    app.testing = True
    client = app.test_client()
    assert client.get("/api/sectors?limit=50").get_json()["total"] == 1
    body = client.get("/api/sectors?limit=50&coverage=uncharted").get_json()
    assert body["total"] == 4 and body["uncharted_total"] == 4 and body["empty_total"] == 0


@pytest.mark.fresh_database
def test_the_migration_deletes_the_rows_of_a_galaxy_that_has_them_and_keeps_the_rows_something_refers_to(mysql_config):
    load_old_schema(mysql_config, 245)
    assert store.schema_status(mysql_config)[0] == 245
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _star(conn, STAR_CELL)
        _scatter(conn, SCATTER_CELL, built=False)
        for cell in (STAR_CELL, SCATTER_CELL):
            store.uncharted_sector_id(conn, cell, EDGE_PC)
        kept = store.uncharted_sector_id(conn, ROGUE_CELL, EDGE_PC)
        insert_filled_row(conn, "rogue_planets", sector_id=kept, name="Wanderer")
        conn.commit()
        assert {STAR_CELL, SCATTER_CELL, ROGUE_CELL} <= _addresses(conn)
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not ({STAR_CELL, SCATTER_CELL} & _addresses(conn))
        assert ROGUE_CELL in _addresses(conn), "the rogue planet's row is its address"
    finally:
        conn.close()


def test_a_placeholder_does_not_count_as_a_charted_sector_for_the_seed_guard(galaxy):
    conn = store.get_connection(galaxy)
    try:
        _placeholders(conn)
        assert store.has_galaxy_sectors(conn), "the galaxy fixture holds one charted sector"
        conn.execute("DELETE FROM chart_ledger")
        conn.commit()
        assert not store.has_galaxy_sectors(conn)
    finally:
        conn.close()


def test_a_new_plan_clears_the_rows_nothing_lives_in(galaxy):
    conn = store.get_connection(galaxy)
    try:
        _placeholders(conn)
        before = _addresses(conn)
        deleted = uncharted_index.clear_unoccupied(conn)
        after = _addresses(conn)
        assert deleted == len(before - after) and {STAR_CELL, SCATTER_CELL, BOTH_CELL} <= before
        assert CHARTED in after, "a charted sector stays"
        assert not ({STAR_CELL, SCATTER_CELL, BOTH_CELL} & after)
    finally:
        conn.close()


def test_placeholder_rows_add_nothing_to_the_block_region_and_shell_totals(galaxy):
    from planetgen.db import block_stats

    conn = store.get_connection(galaxy)
    try:
        before_blocks = block_stats.rebuild(conn)
        before_regions = conn.execute("SELECT COUNT(*) AS n FROM region_stats").fetchone()["n"]
        _placeholders(conn)
        assert block_stats.rebuild(conn) == before_blocks
        assert conn.execute("SELECT COUNT(*) AS n FROM region_stats").fetchone()["n"] == before_regions
        assert conn.execute("SELECT COUNT(*) AS n FROM galaxy_block_stats WHERE m = 3 AND shell_index IN (10, 12)").fetchone()["n"] == 0
        box = ((-1e6, -1e6, -1e6), (1e6, 1e6, 1e6))
        assert query.galaxy_sector_count_in_box(conn, *box) == 1      # only the charted cell: the placeholders are not generated sectors
    finally:
        conn.close()
