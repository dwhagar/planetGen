# tests/test_uncharted_index.py

"""
MAP.402: every cell a scatter left a waiting bright star or an unbuilt scattered phenomenon in has a placeholder sector row (`db.uncharted_index`), made when
a scatter ends (`run_plan.refresh_block_lights_unit`), by `planetgen fill-block-totals` and by migration 0226, so the Sectors list and the search reach it.
"""

import pytest

from planetgen.api.config import Config
from planetgen.db import chart_ledger, query, store, uncharted_index
from planetgen.web.app import create_app
from tests.db_schema_support import load_old_schema

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


def test_a_row_is_made_for_each_cell_that_holds_a_waiting_object_and_for_no_other(galaxy):
    conn = store.get_connection(galaxy)
    try:
        before = _addresses(conn)
        made = uncharted_index.index(conn, edge_pc=EDGE_PC)
        assert made == 3
        assert _addresses(conn) - before == {STAR_CELL, SCATTER_CELL, BOTH_CELL}
        assert BUILT_STAR_CELL not in _addresses(conn) and BUILT_SCATTER_CELL not in _addresses(conn)
        assert query.count_uncharted_sectors(conn) == 4          # the rogue's placeholder and the three made here, none of them charted
    finally:
        conn.close()


def test_a_second_run_adds_nothing_and_a_made_row_is_a_complete_placeholder(galaxy):
    conn = store.get_connection(galaxy)
    try:
        uncharted_index.index(conn, edge_pc=EDGE_PC)
        assert uncharted_index.index(conn, edge_pc=EDGE_PC) == 0
        row = conn.execute("SELECT * FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", STAR_CELL).fetchone()
        assert row["uid"] is not None and row["center_x_pc"] is not None and row["edge_mpc"] == EDGE_PC * 1000.0
        assert row["name"] and query.sector_is_uncharted(conn, row)
        again = store.uncharted_sector_id(conn, STAR_CELL, EDGE_PC)
        assert again == row["id"], "the rogue landing and the index share one row per cell"
    finally:
        conn.close()


def test_the_cells_it_made_are_listed_as_uncharted_and_not_empty(galaxy):
    conn = store.get_connection(galaxy)
    try:
        uncharted_index.index(conn, edge_pc=EDGE_PC)
        listed = query.list_sectors(conn, coverage="uncharted", show_empty=False)
        assert {(r["shell_index"], r["band_index"], r["wedge_index"]) for r in listed} == {STAR_CELL, SCATTER_CELL, BOTH_CELL, ROGUE_CELL}
        assert all(r["uncharted"] for r in listed)
    finally:
        conn.close()


def test_a_cell_whose_object_is_built_later_keeps_its_row_as_an_empty_one(galaxy):
    conn = store.get_connection(galaxy)
    try:
        uncharted_index.index(conn, edge_pc=EDGE_PC)
        conn.execute("UPDATE phenomenon_scatter SET built_at = CURRENT_TIMESTAMP WHERE shell_index = ?", (SCATTER_CELL[0],))
        conn.commit()
        assert uncharted_index.index(conn, edge_pc=EDGE_PC) == 0
        shown = {(r["shell_index"], r["band_index"], r["wedge_index"]) for r in query.list_sectors(conn, coverage="uncharted", show_empty=False)}
        assert SCATTER_CELL not in shown
        every = {(r["shell_index"], r["band_index"], r["wedge_index"]) for r in query.list_sectors(conn, coverage="uncharted", show_empty=True)}
        assert SCATTER_CELL in every
    finally:
        conn.close()


def test_a_run_with_nothing_waiting_makes_no_row_and_no_run(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        runs = conn.execute("SELECT COUNT(*) AS n FROM generation_runs").fetchone()["n"]
        assert uncharted_index.index(conn) == 0
        assert conn.execute("SELECT COUNT(*) AS n FROM generation_runs").fetchone()["n"] == runs
    finally:
        conn.close()


def test_the_api_lists_the_made_cells_only_when_asked(galaxy):
    conn = store.get_connection(galaxy)
    try:
        uncharted_index.index(conn, edge_pc=EDGE_PC)
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
def test_the_migration_makes_the_rows_of_a_galaxy_scattered_before_it(mysql_config):
    load_old_schema(mysql_config, 225)
    assert store.schema_status(mysql_config)[0] == 225
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _star(conn, STAR_CELL)
        _scatter(conn, SCATTER_CELL, built=False)
        _scatter(conn, BUILT_SCATTER_CELL, built=True)
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _addresses(conn) == {STAR_CELL, SCATTER_CELL}
    finally:
        conn.close()


def test_a_placeholder_does_not_count_as_a_charted_sector_for_the_seed_guard(galaxy):
    conn = store.get_connection(galaxy)
    try:
        uncharted_index.index(conn, edge_pc=EDGE_PC)
        assert store.has_galaxy_sectors(conn), "the galaxy fixture holds one charted sector"
        conn.execute("DELETE FROM chart_ledger")
        conn.commit()
        assert not store.has_galaxy_sectors(conn)
    finally:
        conn.close()


def test_a_new_plan_clears_the_rows_nothing_lives_in(galaxy):
    conn = store.get_connection(galaxy)
    try:
        uncharted_index.index(conn, edge_pc=EDGE_PC)
        before = _addresses(conn)
        deleted = uncharted_index.clear_unoccupied(conn)
        after = _addresses(conn)
        assert deleted == len(before - after) and {STAR_CELL, SCATTER_CELL, BOTH_CELL, ROGUE_CELL} <= before
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
        uncharted_index.index(conn, edge_pc=EDGE_PC)
        assert block_stats.rebuild(conn) == before_blocks
        assert conn.execute("SELECT COUNT(*) AS n FROM region_stats").fetchone()["n"] == before_regions
        assert conn.execute("SELECT COUNT(*) AS n FROM galaxy_block_stats WHERE m = 3 AND shell_index IN (10, 12)").fetchone()["n"] == 0
        box = ((-1e6, -1e6, -1e6), (1e6, 1e6, 1e6))
        assert query.galaxy_sector_count_in_box(conn, *box) == 1      # only the charted cell: the placeholders are not generated sectors
    finally:
        conn.close()
