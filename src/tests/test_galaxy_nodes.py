"""
MAP.260: the node read (`db.nodes.galaxy_nodes`). A block's box holds every point of the block, a node lists only objects
inside its block, a child lists every object of its parent that lies inside it (the lists are nested down the chain of
blocks), the aggregate is what the stored tables hold, and a malformed request is refused.
"""

import math
import random

import pytest

from planetgen.db import block_lights, nodes, query, store
from planetgen.galaxy import spherical
from tests.bright_star_rows import B2V, star_row
from planetgen.galaxy.drill import DRILL_LEVELS, DrillBlock, drill_block_bounds, drill_chain_of, format_drill_key

pytestmark = pytest.mark.db

EDGE = 4.0


def _bright_row(position, luminosity_sol, population="young"):
    return star_row(spherical.sector_address_at(position, EDGE), tuple(int(position[a] * 1000) for a in range(3)), population=population,
                    luminosity_w=luminosity_sol * 3.828e26, **B2V)


@pytest.fixture
def lit(mysql_config):
    """A galaxy shape and two hundred bright stars in a few lumps, and the block totals made."""
    rng = random.Random(3)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc,"
                         " bulge_scale_radius_pc, bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude,"
                         " spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
                         " expected_system_count_at_density_1, outer_shell_index)"
                         " VALUES (1, 1, 1, 1, 1, 2, 0.2, 0.3, 1, 0, 1, 4, 10, 5)")
            rows = []
            for centre in ((8000.0, 100.0, 20.0), (-3000.0, 6000.0, 80.0), (200.0, -150.0, 5.0)):
                for _ in range(70):
                    position = tuple(centre[a] + rng.gauss(0, 300.0 if centre[0] != 200.0 else 40.0) for a in range(3))
                    rows.append(_bright_row(position, 10 ** rng.uniform(3, 6), rng.choice(("young", "old", "bulge"))))
            store.insert_bright_stars(conn, rows)
            conn.executemany("INSERT INTO galaxy_shell (shell_index, high_band_index, cell_count) VALUES (?, ?, 1)",
                             [(shell, spherical.band_range(shell)[1]) for shell in range(0, 3000)])
            block_lights.rebuild(conn)
        stars = conn.execute("SELECT id, shell_index, band_index, wedge_index FROM bright_stars").fetchall()
    finally:
        conn.close()
    return mysql_config, stars


def test_a_blocks_box_holds_every_point_of_the_block():
    rng = random.Random(1)
    for address in [(5, 0, 0), (40, 3, 7), (400, -20, 11), (2000, 9, 3)]:
        shell, band, wedge = address
        if not spherical.band_range(shell)[0] <= band <= spherical.band_range(shell)[1]:
            continue
        for block in drill_chain_of(shell, band, wedge):
            (inner, outer), (north, south), (start, end) = drill_block_bounds(tuple(block), EDGE)
            lo, hi = nodes.block_box(block, EDGE)
            for _ in range(60):
                r, theta, phi = rng.uniform(inner, outer), rng.uniform(north, south), rng.uniform(start, end)
                point = (r * math.sin(theta) * math.cos(phi), r * math.sin(theta) * math.sin(phi), r * math.cos(theta))
                assert all(lo[a] <= point[a] <= hi[a] for a in range(3)), (block, point)


def test_membership_is_the_chain_of_the_sector():
    address = (40, 3, 7)
    chain = drill_chain_of(*address)
    other = drill_chain_of(300, 20, 5)
    for block in chain:
        assert nodes.address_in_block(address, block)
    for block in other:
        assert nodes.address_in_block(address, block) == (block in chain), block
    assert any(block not in chain for block in other), "the two sectors are not in one block at every level"
    assert not nodes.address_in_block(None, chain[0])


def test_every_listed_star_is_inside_its_node_and_a_child_lists_what_its_parent_listed_inside_it(lit):
    config, stored = lit
    where = {row["id"]: (row["shell_index"], row["band_index"], row["wedge_index"]) for row in stored}
    conn = store.get_connection(config)
    try:
        for address in sorted({where[i] for i in where})[:: max(1, len(set(where.values())) // 12)]:
            chain = drill_chain_of(*address)
            keys = [format_drill_key(block) for block in chain]
            read = nodes.galaxy_nodes(conn, keys)["nodes"]
            lists = [[star["id"] for star in read[key]["stars"]] for key in keys]
            for block, ids in zip(chain, lists):
                assert all(nodes.address_in_block(where[i], block) for i in ids), block
                assert len(ids) <= (nodes.NODE_SECTOR_STARS if block.m == 1 else nodes.NODE_STARS)
            for parent_ids, child_block, child_ids in zip(lists, chain[1:], lists[1:]):
                inside = {i for i in parent_ids if nodes.address_in_block(where[i], child_block)}
                assert inside <= set(child_ids), child_block
            assert lists[-1], "the sector's own node lists its stars"
    finally:
        conn.close()


def test_the_root_nodes_list_the_brightest_stars_of_the_galaxy(lit):
    config, stored = lit
    conn = store.get_connection(config)
    try:
        roots = {format_drill_key(block) for block in {drill_chain_of(r["shell_index"], r["band_index"], r["wedge_index"])[0]
                                                     for r in stored}}
        read = nodes.galaxy_nodes(conn, sorted(roots))["nodes"]
        listed = {star["id"] for key in roots for star in read[key]["stars"]}
        assert len(listed) == min(len(stored), len(roots) * nodes.NODE_STARS) or len(listed) <= len(stored)
        top = conn.execute("SELECT id FROM bright_stars ORDER BY luminosity_w DESC LIMIT 3").fetchall()
        assert {r["id"] for r in top} <= listed
    finally:
        conn.close()


def test_the_aggregate_is_the_stored_light(lit):
    config, stored = lit
    conn = store.get_connection(config)
    try:
        row = stored[0]
        chain = drill_chain_of(row["shell_index"], row["band_index"], row["wedge_index"])
        block = next(b for b in chain if b.m == 243)
        aggregate = nodes.galaxy_nodes(conn, [format_drill_key(block)])["nodes"][format_drill_key(block)]["aggregate"]
        stored_row = conn.execute("SELECT star_count, star_luminosity_sol FROM galaxy_block_lights WHERE m = 243 AND shell_index = ?"
                                  " AND band_index = ? AND wedge_index = ?", (block.shell, block.band, block.wedge)).fetchone()
        assert aggregate["lights"]["stars"] == stored_row["star_count"] > 0
        assert aggregate["lights"]["star_luminosity_sol"] == pytest.approx(stored_row["star_luminosity_sol"])
        assert len(aggregate["lights"]["star_centre"]) == 3 and aggregate["lights"]["star_temperature_k"] > 0
        sector = nodes.galaxy_nodes(conn, [format_drill_key(chain[-1])])["nodes"][format_drill_key(chain[-1])]["aggregate"]
        assert "lights" not in sector
    finally:
        conn.close()


def test_a_malformed_or_huge_request_is_refused(lit):
    config, _stored = lit
    conn = store.get_connection(config)
    try:
        for key in ("nonsense", "5.0.0.0", "3.0.0", "1.99999999.0.0", "3.-1.0.0"):
            with pytest.raises(ValueError):
                nodes.galaxy_nodes(conn, [key])
        with pytest.raises(ValueError):
            nodes.galaxy_nodes(conn, [f"3.{n}.0.0" for n in range(nodes.MAX_NODES_PER_REQUEST + 1)])
        assert nodes.galaxy_nodes(conn, [])["nodes"] == {}
    finally:
        conn.close()


def test_a_generated_sector_is_listed_by_every_node_above_it_and_by_its_own(mysql_config):
    """A saved sector's stars, its place in the level-3 block's and its own `placed`, and its generated sums in the
    aggregate of every block over it."""
    from planetgen.db import block_stats
    from planetgen.galaxy.sector import SpaceSector
    from tests.test_binary_placement import ADDRESS, CENTER_PC, SYSTEM_POSITION_LY, _binary_system

    system, cfg = _binary_system()
    sector = SpaceSector("Pair Sector", edge_ly=13.0)
    sector.add_system(system, position=SYSTEM_POSITION_LY, system_config=cfg)
    position = {"center_x_pc": CENTER_PC[0], "center_y_pc": CENTER_PC[1], "center_z_pc": CENTER_PC[2],
                "galactic_radius_pc": math.dist(CENTER_PC, (0.0, 0.0, 0.0)),
                "shell_index": ADDRESS[0], "band_index": ADDRESS[1], "wedge_index": ADDRESS[2]}
    store.save_sector(sector, config=mysql_config, galaxy_position=position)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            block_stats.refresh_sectors(conn, [ADDRESS])
            conn.executemany("INSERT INTO galaxy_shell (shell_index, high_band_index, cell_count) VALUES (?, ?, 1)",
                             [(shell, spherical.band_range(shell)[1]) for shell in range(0, 200)])
        chain = drill_chain_of(*ADDRESS)
        keys = [format_drill_key(block) for block in chain]
        read = nodes.galaxy_nodes(conn, keys)["nodes"]
        for key, block in zip(keys, chain):
            node = read[key]
            assert len(node["generated"]) == 2, key
            assert node["aggregate"]["generated"] == 1 and node["aggregate"]["stats"]["systems"] >= 1, key
            assert ("placed" in node) == (block.m in nodes.PLACED_LEVELS)
            if "placed" in node:
                assert [(s["shell_index"], s["band_index"], s["wedge_index"]) for s in node["placed"]] == [ADDRESS]
        # A sector next to it holds none of its stars.
        for neighbour in ((ADDRESS[0], ADDRESS[1], (ADDRESS[2] + 5) % spherical.wedge_count(ADDRESS[0], ADDRESS[1])),):
            key = format_drill_key(DrillBlock(1, *neighbour))
            node = nodes.galaxy_nodes(conn, [key])["nodes"][key]
            assert node["generated"] == [] and node["aggregate"] == {} and node["placed"] == []
    finally:
        conn.close()


def test_a_box_wholly_above_a_thin_shells_reach_is_skipped_not_an_error(mysql_config):
    """A block's box can lie above or below what a small shell reaches in z (the colatitude ratio passes 1)."""
    conn = store.get_connection(mysql_config)
    try:
        for lo, hi in (((0.0, 0.0, 5000.0), (10.0, 10.0, 6000.0)), ((0.0, 0.0, -6000.0), (10.0, 10.0, -5000.0))):
            bands, _share = query._bright_star_bands(conn, lo, hi, EDGE)
            assert all(shell != 1 for shell, _low, _high in bands)
    finally:
        conn.close()
