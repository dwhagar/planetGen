"""
Tests for the Galaxy Map drill-down's block ladder:
`planetgen/galaxy/drill.py` and its page twin, the `drill*` functions in
`html/static/galaxyprisms.js` (run under node; those parity tests are
skipped where node isn't installed). The ladder must nest -- every block
inside exactly one parent, a parent's sectors exactly its children's --
and both sides must agree on every count and address.
"""

import json
import os
import random
import shutil
import subprocess

import pytest

from planetgen.galaxy.drill import (
    DRILL_LEVELS, DrillBlock, drill_block_sectors, drill_chain_of, drill_children, drill_parent, drill_slabs,
    drill_wedge_count, format_drill_key, parse_drill_key,
)
from planetgen.galaxy.geometry import ring_sector_count

NODE = shutil.which("node")
MODULE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "html", "static", "galaxyprisms.js")


def _run(script):
    """Runs `script` as an ES module with galaxyprisms.js imported as `P`;
    returns its JSON output."""
    source = f"import * as P from {json.dumps('file://' + os.path.abspath(MODULE))};\n{script}\n"
    result = subprocess.run([NODE, "--input-type=module", "-e", source], capture_output=True, text=True, check=True)
    return json.loads(result.stdout)

EDGE_RING = 3855
"""The default galaxy's edge ring (docs/design/galaxy-drilldown-navigation.md, section 3.2)."""
EDGE_LAYER = 317


def _rings(m):
    return range(0, EDGE_RING // m + 1)


def _sample_sectors(count=400, seed=11):
    rng = random.Random(seed)
    rings = [0, 1, 2, 3, 8, 9, 26, 27, 80, 242, 243, 728, 729, EDGE_RING] + [rng.randrange(EDGE_RING + 1) for _ in range(count)]
    sectors = []
    for i in rings:
        n = ring_sector_count(i)
        for k in {0, n - 1, rng.randrange(n)}:
            sectors.append((i, rng.randint(-EDGE_LAYER, EDGE_LAYER), k))
    sectors += [(0, 0, 0), (1705, -20, 3225), (5, 1, 2), (5, -1, 2), (5, -2, 2), (5, 2, 2)]
    return sectors


def test_the_design_docs_example_chain():
    """Section 5.4's breadcrumb: Sector 1,705·-20·3,225."""
    assert drill_chain_of(1705, -20, 3225) == [
        DrillBlock(243, 7, 14, 0), DrillBlock(27, 63, 115, -1), DrillBlock(3, 568, 1036, -7), DrillBlock(1, 1705, 3225, -20),
    ]


@pytest.mark.parametrize("m", [27, 3])
def test_child_wedge_counts_are_whole_multiples_of_their_parents(m):
    for i in _rings(m):
        child, parent = drill_wedge_count(m, i), drill_wedge_count(m * 9, i // 9)
        assert child % parent == 0 and child >= parent


@pytest.mark.parametrize("m, low, high", [(243, 0.85, 1.18), (27, 0.93, 1.07), (3, 0.94, 1.07)])
def test_wedge_arcs_stay_near_one_block_edge(m, low, high):
    """Section 3.4: a wedge's centerline arc against a block edge."""
    import math

    for i in _rings(m):
        ratio = 2 * math.pi * (i + 0.5) / drill_wedge_count(m, i)
        assert low <= ratio <= high, (m, i, ratio)


def test_slabs_nest_and_center_on_the_plane():
    assert drill_slabs(None, EDGE_LAYER) == [-1, 0, 1]
    assert drill_slabs(None, 121) == [0] and drill_slabs(None, 122) == [-1, 0, 1]
    assert drill_slabs(None) == [0]
    assert drill_slabs(DrillBlock(243, 0, 0, -1)) == list(range(-13, -4))
    assert drill_slabs(DrillBlock(3, 0, 0, 2)) == [5, 6, 7]
    for m in (243, 27, 3):
        f = 3 if m == 3 else 9
        for slab in range(-5, 6):
            for child in drill_slabs(DrillBlock(m, 0, 0, slab)):
                # The parent slab of every child slab is the slab it came from.
                assert drill_parent(DrillBlock(m // f, 0, 0, child)).slab == slab


@pytest.mark.parametrize("sector", _sample_sectors(40, seed=3))
def test_every_sector_is_in_exactly_its_chains_blocks(sector):
    """Walking down from the top through drill_children reaches the sector
    along exactly drill_chain_of's blocks."""
    ring, layer, slot = sector
    chain = drill_chain_of(ring, layer, slot)
    for parent, child in zip(chain, chain[1:]):
        members = [b for _slab, blocks in drill_children(parent) for b in blocks]
        assert members.count(child) == 1
        assert drill_parent(child) == parent


@pytest.mark.parametrize("block", [DrillBlock(27, 0, 0, 0), DrillBlock(27, 1, 2, -1), DrillBlock(27, 63, 115, -1),
                                   DrillBlock(27, 40, 200, 3)])
def test_a_blocks_sectors_are_exactly_its_childrens(block):
    """Every sector in a level-27 block's ring and layer range has that
    block in its chain exactly when it is a sector of one of its level-3
    children."""
    via_children = set()
    for _slab, children in drill_children(block):
        for child in children:
            for layer in drill_slabs(child):
                via_children.update(drill_block_sectors(child, layer))
    # A sector's level-27 block depends on its layer only through the
    # slab, so test each column once and the layers separately.
    layers = range(block.slab * 27 - 13, block.slab * 27 + 14)
    assert all(drill_chain_of(0, j, 0)[1].slab == block.slab for j in layers)
    assert drill_chain_of(0, layers[0] - 1, 0)[1].slab != block.slab
    assert drill_chain_of(0, layers[-1] + 1, 0)[1].slab != block.slab
    direct = set()
    for i in range(block.ring * 27, block.ring * 27 + 27):
        for k in range(ring_sector_count(i)):
            b27 = drill_chain_of(i, block.slab * 27, k)[1]
            if (b27.ring, b27.wedge) == (block.ring, block.wedge):
                direct.update(DrillBlock(1, i, k, j) for j in layers)
    assert direct and via_children == direct


def test_level_three_blocks_hold_about_nine_sectors_a_layer():
    sizes = [len(drill_block_sectors(DrillBlock(3, ring, wedge, 0), 0))
             for ring in range(1, 1286, 37) for wedge in range(0, drill_wedge_count(3, ring), 97)]
    assert min(sizes) >= 1 and max(sizes) <= 15
    assert 8.0 <= sum(sizes) / len(sizes) <= 10.5


def test_galaxy_children_are_every_top_block():
    groups = drill_children(None, max_ring=EDGE_RING, max_layer=EDGE_LAYER)
    assert [slab for slab, _blocks in groups] == [-1, 0, 1]
    assert len(groups[1][1]) == sum(drill_wedge_count(243, i) for i in range(16))
    with pytest.raises(ValueError):
        drill_children(None)


def test_keys_round_trip_and_reject_bad_blocks():
    block = DrillBlock(27, 63, 115, -1)
    assert format_drill_key(block) == "27.63.115.-1"
    assert parse_drill_key("27.63.115.-1") == block
    for bad in ("", None, "27.63.115", "1.0.0.0", "81.0.0.0", "27.-1.0.0", "27.0.99.0", "a.b.c.d", "27.0.0.0.0"):
        with pytest.raises(ValueError):
            parse_drill_key(bad)


def test_levels():
    assert DRILL_LEVELS == (243, 27, 3, 1)
    with pytest.raises(ValueError):
        drill_wedge_count(81, 0)


# --- The page's copy agrees --------------------------------------------------

node_only = pytest.mark.skipif(NODE is None, reason="node is not installed")


@node_only
def test_page_wedge_counts_match():
    rings = {m: list(_rings(m)) for m in (243, 27, 3)}
    rings[1] = list(range(0, EDGE_RING + 1, 7))
    got = _run(f"""
const rings = {json.dumps({str(m): r for m, r in rings.items()})};
const out = {{}};
for (const m of Object.keys(rings)) out[m] = rings[m].map(i => P.drillWedgeCount(Number(m), i));
console.log(JSON.stringify(out));
""")
    for m, ring_list in rings.items():
        assert got[str(m)] == [drill_wedge_count(m, i) for i in ring_list], m


@node_only
def test_page_chains_match():
    sectors = _sample_sectors()
    got = _run(f"""
console.log(JSON.stringify({json.dumps(sectors)}.map(([i, j, k]) =>
  P.drillChainOf(i, j, k).map(b => [b.m, b.ring, b.wedge, b.slab]))));
""")
    assert got == [[list(b) for b in drill_chain_of(*sector)] for sector in sectors]


@node_only
def test_page_children_and_slabs_match():
    chains = [drill_chain_of(*sector) for sector in _sample_sectors(30, seed=5)]
    blocks = sorted({b for chain in chains for b in chain[:3]})
    got = _run(f"""
const blocks = {json.dumps([list(b) for b in blocks])};
console.log(JSON.stringify({{
  children: blocks.map(([m, ring, wedge, slab]) =>
    P.drillChildren({{m, ring, wedge, slab}}).map(g => [g.slab, g.blocks.map(b => [b.m, b.ring, b.wedge, b.slab])])),
  galaxy: P.drillChildren(null, {EDGE_RING}, {EDGE_LAYER}).map(g => [g.slab, g.blocks.length]),
  galaxySlabs: [P.drillSlabs(null), P.drillSlabs(null, 121), P.drillSlabs(null, 122)],
  keys: ["27.63.115.-1", "1.0.0.0", "27.0.99.0", "x", "243.15.0.1"].map(k => P.parseDrillKey(k)),
  formatted: P.formatDrillKey({{m: 3, ring: 568, wedge: 1036, slab: -7}}),
}}));
""")
    assert got["children"] == [[[slab, [list(b) for b in group]] for slab, group in drill_children(b)] for b in blocks]
    assert got["galaxy"] == [[slab, len(group)] for slab, group in drill_children(None, EDGE_RING, EDGE_LAYER)]
    assert got["galaxySlabs"] == [[0], [0], [-1, 0, 1]]
    assert got["keys"][0] == {"m": 27, "ring": 63, "wedge": 115, "slab": -1}
    assert got["keys"][1:4] == [None, None, None]
    assert got["keys"][4] == {"m": 243, "ring": 15, "wedge": 0, "slab": 1}
    assert got["formatted"] == "3.568.1036.-7"
