"""
MAP.173: a node lists the groups of the nebulae inside it (`nebula_regions`): its child blocks that hold nebulae, the
largest by volume first, from the stored aggregates (`db/nebula_blocks.py`).
"""

import pytest

from planetgen import tuning
from planetgen.db import nebula_blocks, nodes, store
from planetgen.galaxy.drill import DrillBlock, drill_children, format_drill_key

from tests.test_nebula_blocks_db import _clouds, _plant_field, _store_nebula

pytestmark = pytest.mark.db

EDGE_PC = float(tuning.DEFAULT_SECTOR_EDGE_PC)


@pytest.fixture
def built(mysql_config):
    """Sixty field clouds and three stored nebulae among them, and the block aggregates made."""
    conn = store.get_connection(mysql_config)
    try:
        _plant_field(conn, _clouds(5, 60))
        for letter, radius_ly, center in (("C", 3.0, (805.0, 195.0, 52.0)), ("F", 8.0, (830.0, 210.0, 40.0)),
                                          ("H", 2.5, (770.0, 190.0, 60.0))):
            _store_nebula(conn, letter, radius_ly, center)
        nebula_blocks.rebuild(conn, EDGE_PC)
        conn.commit()
    finally:
        conn.close()
    return mysql_config


def _blocks_at(conn, m):
    return [DrillBlock(m, r["shell_index"], r["band_index"], r["wedge_index"]) for r in conn.execute(
        "SELECT shell_index, band_index, wedge_index FROM nebula_block_stats WHERE m = ? AND nebula_count > 0",
        (m,)).fetchall()]


def _regions(conn, block):
    return nodes.galaxy_nodes(conn, [format_drill_key(block)])["nodes"][format_drill_key(block)].get("nebula_regions", [])


def test_a_nodes_regions_hold_every_nebula_of_the_node_and_the_largest_come_first(built):
    conn = store.get_connection(built)
    try:
        for m in (243, 81, 27, 9):
            checked = 0
            for block in _blocks_at(conn, m)[:6]:
                regions = _regions(conn, block)
                own = nebula_blocks.read(conn, m, [(block.shell, block.band, block.wedge)])
                assert regions, (m, block)
                assert sum(r["count"] for r in regions) == next(iter(own.values()))["count"]
                volumes = [r["volume_pc3"] for r in regions]
                assert volumes == sorted(volumes, reverse=True)
                keys = {format_drill_key(child) for _shell, group in drill_children(block) for child in group}
                assert {r["key"] for r in regions} <= keys
                for region in regions:
                    assert region["count"] > 0 and len(region["centre"]) == 3 and region["radius_pc"] > 0
                    assert region["class"] and sum(region["shares"].values()) == pytest.approx(1.0)
                checked += 1
            assert checked
    finally:
        conn.close()


def test_a_regions_children_are_regions_of_its_own_node(built):
    conn = store.get_connection(built)
    try:
        block = max(_blocks_at(conn, 81), key=lambda b: 0)
        regions = _regions(conn, block)
        child = DrillBlock(27, *next((c.shell, c.band, c.wedge) for _s, group in drill_children(block) for c in group
                                     if format_drill_key(c) == regions[0]["key"]))
        below = _regions(conn, child)
        assert sum(r["count"] for r in below) == regions[0]["count"]
    finally:
        conn.close()


def test_the_regions_of_a_node_are_cut_to_the_budget(built, monkeypatch):
    conn = store.get_connection(built)
    try:
        block = max(_blocks_at(conn, 243), key=lambda b: len(_regions(conn, b)))
        everything = _regions(conn, block)
        assert len(everything) > 2
        monkeypatch.setattr(nodes, "NODE_NEBULA_REGIONS", 2)
        assert [r["key"] for r in nodes.node_nebula_regions(conn, block, 2)] == [r["key"] for r in everything[:2]]
        assert nodes.NODE_NEBULA_REGIONS == 2
    finally:
        conn.close()


def test_a_block_with_no_nebulae_or_no_finer_level_has_no_regions(built):
    conn = store.get_connection(built)
    try:
        assert nodes.node_nebula_regions(conn, DrillBlock(3, 0, 0, 0)) == []
        assert nodes.node_nebula_regions(conn, DrillBlock(1, 0, 0, 0)) == []
        empty = DrillBlock(243, 40, 0, 0)
        assert empty not in _blocks_at(conn, 243)
        assert nodes.node_nebula_regions(conn, empty) == []
        key = format_drill_key(empty)
        assert "nebula_regions" not in nodes.galaxy_nodes(conn, [key])["nodes"][key]
    finally:
        conn.close()
