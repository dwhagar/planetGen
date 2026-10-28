"""
DB.107: the stored regions (`region_stats`). A row is a 9-block holding at least one generated sector, with its sector count,
capacity, clump count and extent (`named_regions.region_extent`); `block_stats`' hooks keep it, a refresh of every chain
in any order equals a rebuild, one generated sector makes a region at a small share charted, and deleting the last sector
of a region takes the region away.
"""

import random

import numpy as np
import pytest

from planetgen.db import block_stats, edits, query, region_stats, store
from planetgen.galaxy import drill, named_regions, spherical

pytestmark = pytest.mark.db


def _cluster(start, size):
    """`size` sectors grown by face-neighbours from `start`, in a fixed order."""
    found = [start]
    frontier = [start]
    while frontier and len(found) < size:
        for neighbour in sorted(spherical.neighbor_addresses(*frontier.pop(0))):
            if neighbour not in found and len(found) < size:
                found.append(neighbour)
                frontier.append(neighbour)
    return found


def _addresses():
    """A cluster of 18 sectors, a lone sector in a middle shell and a lone far sector, as sorted addresses."""
    low, high = spherical.band_range(40)
    cluster = _cluster((40, (low + high) // 2, 0), 18)
    low, high = spherical.band_range(300)
    return sorted(set(cluster + [(300, (low + high) // 2, 5), (91, sum(spherical.band_range(91)) // 2, 0)]))


def _by_region(addresses):
    cells = {}
    for address in addresses:
        cells.setdefault(region_stats.region_cell(address), []).append(address)
    return cells


@pytest.fixture
def galaxy(mysql_config):
    """A database holding the sectors of `_addresses`, the block totals and the regions made."""
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            for n, (shell, band, wedge) in enumerate(_addresses()):
                conn.execute(
                    "INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES (?, 4000, ?, ?, ?)",
                    (f"Sector {n}", shell, band, wedge))
                conn.execute(
                    "INSERT INTO sector_stats (shell_index, band_index, wedge_index, expected_systems, actual_systems,"
                    " actual_stars, mean_age_gy, total_luminosity_sol, bright_level_sol) VALUES (?, ?, ?, 1.5, 1, 2, 1.0, 10.0, 0)",
                    (shell, band, wedge))
            block_stats.rebuild(conn)
    finally:
        conn.close()
    return mysql_config


def _stored(conn):
    return {tuple(region["block"])[1:]: region for region in region_stats.regions(conn)}


def test_a_rebuild_stores_one_row_for_every_region_that_holds_a_sector(galaxy):
    conn = store.get_connection(galaxy)
    try:
        stored = _stored(conn)
        expected = _by_region(_addresses())
        assert set(stored) == {tuple(cell)[1:] for cell in expected}
        edge = query.galaxy_edge_pc(conn)
        for cell, addresses in expected.items():
            row = stored[tuple(cell)[1:]]
            extent = named_regions.region_extent(cell, addresses, edge)
            assert row["block"] == cell
            assert row["sector_count"] == len(addresses)
            assert row["capacity"] == sum(1 for _ in drill.block_addresses(cell))
            assert row["clump_count"] == extent.clumps
            assert row["extent_radius_pc"] == pytest.approx(extent.radius_pc)
            assert row["extent_centre_pc"] == pytest.approx(tuple(extent.centre_pc))
            assert row["extent_radius_pc"] >= named_regions.framing_radius_pc(cell, edge) - 1e-9
    finally:
        conn.close()


def test_the_share_charted_is_the_sectors_over_the_capacity(galaxy):
    conn = store.get_connection(galaxy)
    try:
        for cell, addresses in _by_region(_addresses()).items():
            capacity = sum(1 for _ in drill.block_addresses(cell))
            assert region_stats.charted_percent(conn, cell) == pytest.approx(100.0 * len(addresses) / capacity)
        empty = drill.DrillBlock(9, 5000, 0, 0)
        assert region_stats.charted_percent(conn, empty) == 0.0
        assert region_stats.region_row(conn, empty) is None
    finally:
        conn.close()


def test_one_sector_makes_a_region_at_a_small_share(galaxy):
    conn = store.get_connection(galaxy)
    try:
        lone = (300, sum(spherical.band_range(300)) // 2, 5)
        cell = region_stats.region_cell(lone)
        row = region_stats.region_row(conn, cell)
        assert row["sector_count"] == 1 and row["clump_count"] == 1
        assert 0.0 < row["charted_percent"] < 5.0
        edge = query.galaxy_edge_pc(conn)
        found = region_stats.region_of(conn, spherical.sector_position_pc(*lone, edge), edge)
        assert found["block"] == cell
        assert np.allclose(found["extent_centre_pc"], spherical.sector_position_pc(*lone, edge))
        assert region_stats.region_of(conn, (123456.0, 0.0, 0.0), edge) is None
    finally:
        conn.close()


def test_a_refresh_of_every_chain_in_any_order_equals_a_rebuild(galaxy):
    conn = store.get_connection(galaxy)
    try:
        kept = _stored(conn)
        for order in (list(_addresses()), list(reversed(_addresses())), random.Random(3).sample(_addresses(), len(_addresses()))):
            conn.execute("DELETE FROM region_stats")
            for address in order:
                region_stats.refresh_chain(conn, address)
            assert _stored(conn) == kept
        conn.execute("DELETE FROM region_stats")
        assert region_stats.rebuild(conn) == len(kept)
        assert _stored(conn) == kept
        assert region_stats.refresh_chain(conn, (None, 0, 0)) == 0
    finally:
        conn.close()


def test_the_block_hooks_keep_the_regions(galaxy):
    conn = store.get_connection(galaxy)
    try:
        kept = _stored(conn)
        conn.execute("DELETE FROM region_stats")
        for address in _addresses():
            block_stats.refresh_chain(conn, address)
        assert _stored(conn) == kept
        conn.execute("DELETE FROM region_stats")
        block_stats.refresh_sectors(conn, _addresses())
        assert _stored(conn) == kept
        conn.execute("DELETE FROM region_stats")
        block_stats.rebuild(conn)
        assert _stored(conn) == kept
    finally:
        conn.close()


def test_deleting_sectors_shrinks_a_region_and_the_last_one_takes_it_away(galaxy):
    conn = store.get_connection(galaxy)
    try:
        regions = _by_region(_addresses())
        big = max(regions, key=lambda cell: len(regions[cell]))
        lone = region_stats.region_cell((300, sum(spherical.band_range(300)) // 2, 5))
        for victim in regions[big][:2]:
            sector = conn.execute("SELECT id FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
                                  victim).fetchone()["id"]
            edits.delete_sector_with_contents(conn, sector)
        conn.commit()
        assert region_stats.region_row(conn, big)["sector_count"] == len(regions[big]) - 2
        sector = conn.execute("SELECT id FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
                              regions[lone][0]).fetchone()["id"]
        edits.delete_sector_with_contents(conn, sector)
        conn.commit()
        assert region_stats.region_row(conn, lone) is None
        assert region_stats.charted_percent(conn, lone) == 0.0
        kept = _stored(conn)
        conn.execute("DELETE FROM region_stats")
        region_stats.rebuild(conn)
        assert _stored(conn) == kept
    finally:
        conn.close()


def test_the_regions_come_most_charted_first_and_can_be_asked_for_by_cell(galaxy):
    conn = store.get_connection(galaxy)
    try:
        everything = region_stats.regions(conn)
        shares = [region["charted_percent"] for region in everything]
        assert shares == sorted(shares, reverse=True)
        wanted = [everything[-1]["block"], drill.DrillBlock(9, 5000, 0, 0)]
        assert [region["block"] for region in region_stats.regions(conn, wanted)] == [everything[-1]["block"]]
    finally:
        conn.close()
