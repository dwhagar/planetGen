"""
MAP.314: what a node of 27, 81 or 243 sectors lists on the wire so the Galaxy Map can name the charted regions inside it from a
distance. Such a node carries `regions`, one entry per charted 9-block it holds (`key`, `name`, `parent`, `centre`, `radius`, `sectors`,
`percent`), most sectors first; no other node has the list, nor has any node while no naming key is in force or when the block
holds no generated sector.
"""

import pytest

from planetgen.db import nodes, region_stats, store
from planetgen.galaxy import drill, named_regions
from tests.test_region_wire_names import EDGE, KEY, _cells, _node, galaxy, keyed  # noqa: F401 -- fixtures

pytestmark = pytest.mark.db


def _block_over(cell, level):
    """The block of `level` sectors that holds the 9-block `cell` (the 9-block itself at 9)."""
    block = cell
    while block.m < level:
        block = drill.drill_parent(block)
    return block


def _stored(config):
    conn = store.get_connection(config)
    try:
        return {region["block"]: region for region in region_stats.regions(conn)}
    finally:
        conn.close()


@pytest.mark.parametrize("level", region_stats.LABEL_LEVELS)
def test_a_node_lists_every_charted_region_it_holds_and_only_those(galaxy, keyed, level):
    stored = _stored(galaxy)
    for cell in _cells(galaxy):
        block = _block_over(cell, level)
        inside = {other for other in stored if _block_over(other, level) == block}
        listed = _node(galaxy, block)["regions"]
        assert {entry["key"] for entry in listed} == {drill.format_drill_key(other) for other in inside}, (level, block)
        for entry in listed:
            region = stored[drill.parse_node_key(entry["key"])]
            assert entry["name"] == named_regions.region_name(region["block"], KEY, EDGE)
            assert entry["parent"] == named_regions.family_name(named_regions.family_of(region["block"]), KEY, EDGE)
            assert entry["centre"] == list(region["extent_centre_pc"]) and entry["radius"] == region["extent_radius_pc"]
            assert entry["sectors"] == region["sector_count"]
            assert entry["percent"] == pytest.approx(region["charted_percent"], abs=1e-3)
        assert [entry["sectors"] for entry in listed] == sorted((entry["sectors"] for entry in listed), reverse=True)


def test_a_charted_regions_own_node_says_how_much_of_it_is_charted(galaxy, keyed):
    stored = _stored(galaxy)
    for cell, region in stored.items():
        node = _node(galaxy, cell)
        assert node["charted"] == {"sectors": region["sector_count"], "percent": pytest.approx(region["charted_percent"], abs=1e-3)}
    other = next(entry for entry in _node(galaxy, _block_over(next(iter(stored)), 27))["regions"])
    own = _node(galaxy, drill.parse_node_key(other["key"]))["charted"]
    assert (other["sectors"], other["percent"]) == (own["sectors"], own["percent"])      # the list and the node agree
    assert "charted" not in _node(galaxy, _block_over(next(iter(stored)), 27))


def test_the_list_is_longest_first_and_a_node_holding_several_regions_lists_each(galaxy, keyed):
    stored = _stored(galaxy)
    assert len(stored) >= 3
    top = max(region_stats.LABEL_LEVELS)
    listed = _node(galaxy, _block_over(next(iter(stored)), top))["regions"]
    assert len(listed) >= 1
    assert all(left["sectors"] >= right["sectors"] for left, right in zip(listed, listed[1:]))


def test_no_other_node_has_the_list(galaxy, keyed):
    cell = _cells(galaxy)[0]
    for level in (9, 729, 2187):
        assert "regions" not in _node(galaxy, _block_over(cell, level)), level
    inner = drill.drill_children(cell)[0][1][0]                                   # a 3-block of the region
    assert inner.m == 3 and "regions" not in _node(galaxy, inner)


def test_a_block_with_no_generated_sector_has_no_list(galaxy, keyed):
    charted = set(_cells(galaxy))
    cell = next(iter(charted))
    parent = _block_over(cell, 243)
    empty = next(child for _shell, group in drill.drill_children(parent)
                 for child in group if child.m == 81 and not any(_block_over(other, 81) == child for other in charted))
    assert "regions" not in _node(galaxy, empty)


def test_without_a_naming_key_there_is_no_list(galaxy):
    cell = _cells(galaxy)[0]
    for level in region_stats.LABEL_LEVELS:
        assert "regions" not in _node(galaxy, _block_over(cell, level)), level


def test_cells_in_a_block_are_its_nine_blocks_and_a_smaller_block_is_refused(galaxy):
    cell = _cells(galaxy)[0]
    family = _block_over(cell, 81)
    cells = region_stats.cells_in(family)
    assert cell in cells and len(cells) == len(set(cells)) and all(other.m == 9 for other in cells)
    assert all(drill.drill_parent(drill.drill_parent(other)) == family for other in cells)
    assert region_stats.cells_in(cell) == [cell]
    with pytest.raises(ValueError):
        region_stats.cells_in(drill.drill_children(cell)[0][1][0])


def test_the_packed_node_keeps_the_list(galaxy, keyed):
    from planetgen.web.lib import tilecache, tilepack
    node = _node(galaxy, _block_over(_cells(galaxy)[0], 81))
    unpacked = tilepack.unpack_node(tilecache._node_bytes(node))
    assert [entry["name"] for entry in unpacked["regions"]] == [entry["name"] for entry in node["regions"]]
    assert unpacked["regions"][0]["centre"] == pytest.approx(node["regions"][0]["centre"], abs=1e-3)
    assert "regions" not in tilepack.unpack_node(tilecache._node_bytes({"aggregate": {}, "stars": [], "generated": []}))
