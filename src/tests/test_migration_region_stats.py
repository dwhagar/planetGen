"""
DB.107: migration 0161 creates `region_stats` and fills it from the stored sectors: a row for every 9-block that holds a
generated sector, the same rows a refresh of each sector keeps.
"""

import pytest

from planetgen.db import region_stats, store
from planetgen.galaxy import spherical
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

SECTORS = [(40, 20, 0), (40, 20, 1), (41, 21, 0), (300, 150, 5), (91, 45, 0)]


def _old_galaxy(mysql_config):
    load_old_schema(mysql_config, 160)
    assert store.schema_status(mysql_config)[0] == 160
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for n, (shell, band, wedge) in enumerate(SECTORS):
            conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES (?, 4000, ?, ?, ?)",
                         (f"Sector {n}", shell, band, wedge))
        conn.commit()
    finally:
        conn.close()


def test_the_migration_fills_a_row_for_every_region_that_holds_a_sector(mysql_config):
    _old_galaxy(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        stored = {tuple(r["block"])[1:]: r for r in region_stats.regions(conn)}
        expected = {}
        for address in SECTORS:
            expected.setdefault(tuple(region_stats.region_cell(address))[1:], []).append(address)
        assert set(stored) == set(expected)
        assert {key: row["sector_count"] for key, row in stored.items()} == {key: len(v) for key, v in expected.items()}
        assert all(0.0 < row["charted_percent"] <= 100.0 and row["extent_radius_pc"] > 0.0 for row in stored.values())
        kept = dict(stored)
        for address in SECTORS:
            region_stats.refresh_chain(conn, address)
        assert {tuple(r["block"])[1:]: r for r in region_stats.regions(conn)} == kept
        assert store.schema_status(mysql_config)[0] == store.SCHEMA_VERSION
    finally:
        conn.close()
