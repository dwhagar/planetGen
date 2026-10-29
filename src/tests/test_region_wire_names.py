"""
API.35: the names of charted regions on the wire. A node of a 9-block that holds a generated sector carries `name`
("Fardoo Beez Region") and `parent` (its 81-block family, "Fardoo Reach"), the stage that lists such blocks as children names
each of them (and a 9-block's own stage names it), the packed node keeps both, and nothing else carries them: not a block of
another size, not a 9-block with no generated sector, and nothing at all while no naming key is in force.
"""

import pytest

from planetgen.db import block_stats, nodes, query, region_stats, store
from planetgen.galaxy import drill, named_regions
from planetgen.names import naming_key
from planetgen.web.lib import tilecache, tilepack
from tests.test_region_stats import _addresses, _by_region

pytestmark = pytest.mark.db

KEY = "00AB12CD"
EDGE = 4.0


@pytest.fixture
def galaxy(mysql_config):
    """A database holding the sectors of `test_region_stats._addresses` and the block totals and regions made."""
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


@pytest.fixture
def keyed():
    naming_key.set_resolver(lambda: KEY)
    yield KEY
    naming_key.set_resolver(None)


def _cells(config):
    conn = store.get_connection(config)
    try:
        return [region["block"] for region in region_stats.regions(conn)]
    finally:
        conn.close()


def _node(config, block):
    conn = store.get_connection(config)
    try:
        return nodes.galaxy_nodes(conn, [drill.format_drill_key(block)])["nodes"][drill.format_drill_key(block)]
    finally:
        conn.close()


def test_the_node_of_a_charted_region_carries_its_name_and_its_family(galaxy, keyed):
    for cell in _cells(galaxy):
        node = _node(galaxy, cell)
        assert node["name"] == named_regions.region_name(cell, KEY, EDGE)
        assert node["parent"] == named_regions.family_name(named_regions.family_of(cell), KEY, EDGE)
        assert node["name"].endswith(" Region") and node["parent"].endswith(" Reach")
        assert node["name"].split(" ")[0] == node["parent"].split(" ")[0]


def test_no_other_node_carries_a_name(galaxy, keyed):
    cell = _cells(galaxy)[0]
    for level in (1, 3, 27, 81, 243):
        block = drill.drill_chain_of(*(cell.shell, cell.band, cell.wedge))[drill.DRILL_LEVELS.index(level)]
        node = _node(galaxy, block)
        assert "name" not in node and "parent" not in node, level


def test_a_nine_block_with_no_generated_sector_has_no_name(galaxy, keyed):
    charted = set(_cells(galaxy))
    sibling = next(child for _shell, group in drill.drill_children(drill.drill_parent(next(iter(charted))))
                   for child in group if child not in charted and child.m == named_regions.REGION_LEVEL)
    node = _node(galaxy, sibling)
    assert "name" not in node and "parent" not in node


def test_without_a_naming_key_nothing_is_named(galaxy):
    assert naming_key.active_key() is None
    node = _node(galaxy, _cells(galaxy)[0])
    assert "name" not in node and "parent" not in node


def test_the_node_of_a_charted_region_carries_its_extent_with_or_without_a_key(galaxy):
    """MAP.319: the bounding sphere the map frames the region with rides on the node, named or not."""
    assert naming_key.active_key() is None
    for cell in _cells(galaxy):
        conn = store.get_connection(galaxy)
        try:
            region = region_stats.region_row(conn, cell)
        finally:
            conn.close()
        node = _node(galaxy, cell)
        assert node["extent"] == {"centre": list(region["extent_centre_pc"]), "radius": region["extent_radius_pc"]}
        assert node["extent"]["radius"] >= named_regions.framing_radius_pc(cell, EDGE)
        assert "name" not in node


def test_no_other_node_carries_an_extent(galaxy, keyed):
    cell = _cells(galaxy)[0]
    for level in (1, 3, 27, 81, 243):
        block = drill.drill_chain_of(*(cell.shell, cell.band, cell.wedge))[drill.DRILL_LEVELS.index(level)]
        assert "extent" not in _node(galaxy, block), level
    charted = set(_cells(galaxy))
    sibling = next(child for _shell, group in drill.drill_children(drill.drill_parent(cell))
                   for child in group if child not in charted and child.m == named_regions.REGION_LEVEL)
    assert "extent" not in _node(galaxy, sibling)


def test_the_stage_names_each_charted_nine_block_it_lists_and_a_nine_block_names_itself(galaxy, keyed):
    cell = _cells(galaxy)[0]
    conn = store.get_connection(galaxy)
    try:
        above = query.galaxy_stage(conn, drill.format_drill_key(drill.drill_parent(cell)))
        assert above["child_m"] == named_regions.REGION_LEVEL
        named = {(c["shell"], c["band"], c["wedge"]): c for c in above["children"]}
        assert named[(cell.shell, cell.band, cell.wedge)]["name"] == named_regions.region_name(cell, KEY, EDGE)
        assert all("name" in child and "parent" in child for child in above["children"])
        own = query.galaxy_stage(conn, drill.format_drill_key(cell))
        assert own["name"] == named_regions.region_name(cell, KEY, EDGE) and "name" not in own["children"][0]
        top = query.galaxy_stage(conn, None)
        assert "name" not in top and all("name" not in child for child in top["children"])
    finally:
        conn.close()


def test_the_packed_node_keeps_the_name_and_the_family(galaxy, keyed):
    cell = _cells(galaxy)[0]
    node = _node(galaxy, cell)
    unpacked = tilepack.unpack_node(tilecache._node_bytes(node))
    assert unpacked["name"] == node["name"] and unpacked["parent"] == node["parent"]
    assert unpacked["extent"]["centre"] == pytest.approx(node["extent"]["centre"], abs=1e-3)
    assert unpacked["extent"]["radius"] == pytest.approx(node["extent"]["radius"], abs=1e-3)
    plain = tilepack.unpack_node(tilecache._node_bytes({"aggregate": {}, "stars": [], "generated": []}))
    assert "name" not in plain and "parent" not in plain and "extent" not in plain
