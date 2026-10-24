# tests/test_migration_block_ladder_root.py

"""
MAP.258: migration 0146 gives the block ladder its root of three blocks. `galaxy_block_stats.m` becomes INT UNSIGNED and
the table is filled again with a row for every level, the old rows (made for blocks that moved when the levels above 243
were put in) are gone, and the nebula block tables are made again when they were filled. (A finished migration run again
changes nothing: `test_db_migration_rerun.py`.)
"""

import pytest

from planetgen.db import store
from planetgen.galaxy.drill import DRILL_LEVELS, DRILL_TOP, drill_chain_of
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

SECTORS = [(30, 3, 5), (30, 3, 6), (31, 3, 5), (400, 40, 11)]
"""Generated sectors `(shell, band, wedge)`, two of them in one level-3 block and all four in the root's wedge 0."""


def _old_galaxy(mysql_config):
    load_old_schema(mysql_config, 134)
    assert store.schema_status(mysql_config)[0] == 134
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for n, (shell, band, wedge) in enumerate(SECTORS):
            conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES (?, 4000, ?, ?, ?)",
                         (f"Sector {n}", shell, band, wedge))
            conn.execute(
                "INSERT INTO sector_stats (shell_index, band_index, wedge_index, expected_systems, actual_systems, actual_stars,"
                " mean_age_gy, total_luminosity_sol, bright_level_sol) VALUES (?, ?, ?, 2.5, 2, 3, 4.0, 12.5, 0)",
                (shell, band, wedge))
        # A row of the old ladder: made for a block that is another block now.
        conn.execute("INSERT INTO galaxy_block_stats (m, shell_index, band_index, wedge_index, sector_count, systems,"
                     " expected_systems, stars, age_sum, luminosity_sol) VALUES (243, 7, 0, 0, 99, 99, 99, 99, 99, 99)")
        conn.commit()
    finally:
        conn.close()


def _column_type(conn):
    return conn.execute("SELECT COLUMN_TYPE AS t FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                        " AND TABLE_NAME = 'galaxy_block_stats' AND COLUMN_NAME = 'm'").fetchone()["t"].lower()


def test_the_migration_fills_every_level_and_drops_the_rows_of_the_old_ladder(mysql_config):
    _old_galaxy(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _column_type(conn) == "int unsigned"
        rows = conn.execute("SELECT m, shell_index, band_index, wedge_index, sector_count, systems FROM galaxy_block_stats"
                            ).fetchall()
        assert {r["m"] for r in rows} == set(DRILL_LEVELS[:-1])
        assert not [r for r in rows if r["sector_count"] == 99]
        stored = {(r["m"], r["shell_index"], r["band_index"], r["wedge_index"]): r["sector_count"] for r in rows}
        expected = {}
        for address in SECTORS:
            for block in drill_chain_of(*address)[:-1]:
                expected[tuple(block)] = expected.get(tuple(block), 0) + 1
        assert stored == expected
        roots = [count for key, count in stored.items() if key[0] == DRILL_TOP]
        assert sum(roots) == len(SECTORS)
        assert store.schema_status(mysql_config)[0] == store.SCHEMA_VERSION
    finally:
        conn.close()


def test_the_nebula_blocks_are_made_again_only_when_they_were_filled(mysql_config):
    _old_galaxy(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("INSERT INTO nebula_field (cell_i, cell_j, cell_k, n, nebula_class, radius_ly, center_x_pc, center_y_pc,"
                     " center_z_pc) VALUES (1, 0, 0, 0, 'M', 40.0, 810.0, 200.0, 50.0)")
        conn.execute("INSERT INTO nebula_block_state (id, nebulae_through_id) VALUES (1, 0)")
        conn.execute("INSERT INTO nebula_block_stats (m, shell_index, band_index, wedge_index, nebula_count, volume_pc3,"
                     " volume_diffuse_pc3, volume_emission_pc3, volume_reflection_pc3, volume_planetary_pc3, volume_dark_pc3,"
                     " volume_remnant_pc3, cover_volume_pc3, cover) VALUES (243, 7, 0, 0, 5, 1, 0, 0, 0, 0, 0, 0, 1, 0.5)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not conn.execute("SELECT 1 FROM nebula_block_stats WHERE nebula_count = 5").fetchall()
        assert conn.execute("SELECT nebula_count FROM nebula_block_stats WHERE m = 3").fetchall()
        assert {r["m"] for r in conn.execute("SELECT DISTINCT m FROM nebula_block_stats").fetchall()} == {243, 81, 27, 9, 3}
    finally:
        conn.close()
