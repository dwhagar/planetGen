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
    """Section 5.4's breadcrumb: Sector 1,705·-20·3,225, now down the 3-ary ladder (MAP.151b)."""
    assert drill_chain_of(1705, -20, 3225) == [
        DrillBlock(243, 7, 14, 0), DrillBlock(81, 21, 43, 0), DrillBlock(27, 63, 129, -1), DrillBlock(9, 189, 388, -2),
        DrillBlock(3, 568, 1166, -7), DrillBlock(1, 1705, 3225, -20),
    ]


@pytest.mark.parametrize("m", [81, 27, 9, 3])
def test_child_wedge_counts_are_whole_multiples_of_their_parents(m):
    for i in _rings(m):
        child, parent = drill_wedge_count(m, i), drill_wedge_count(m * 3, i // 3)
        assert child % parent == 0 and child >= parent


@pytest.mark.parametrize("m, low, high", [(243, 0.85, 1.18), (81, 0.82, 1.23), (27, 0.82, 1.25), (9, 0.83, 1.25),
                                          (3, 0.83, 1.25)])
def test_wedge_arcs_stay_near_one_block_edge(m, low, high):
    """Section 3.4: a wedge's centerline arc against a block edge (0.83 to 1.25 on the 3-ary ladder)."""
    import math

    for i in _rings(m):
        ratio = 2 * math.pi * (i + 0.5) / drill_wedge_count(m, i)
        assert low <= ratio <= high, (m, i, ratio)


def test_slabs_nest_and_center_on_the_plane():
    assert drill_slabs(DrillBlock(243, 0, 0, -1)) == [-4, -3, -2]
    assert drill_slabs(DrillBlock(3, 0, 0, 2)) == [5, 6, 7]
    for m in (243, 81, 27, 9, 3):
        for slab in range(-5, 6):
            for child in drill_slabs(DrillBlock(m, 0, 0, slab)):
                # The parent slab of every child slab is the slab it came from.
                assert drill_parent(DrillBlock(m // 3, 0, 0, child)).slab == slab


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


def _sectors_through_children(block):
    """A block's sectors by walking down its children to the level-3 blocks."""
    found = set()
    for _slab, children in drill_children(block):
        for child in children:
            if child.m == 3:
                for layer in drill_slabs(child):
                    found.update(drill_block_sectors(child, layer))
            elif child.m == 1:
                found.add(child)
            else:
                found.update(_sectors_through_children(child))
    return found


@pytest.mark.parametrize("block", [DrillBlock(27, 0, 0, 0), DrillBlock(27, 1, 2, -1), DrillBlock(27, 63, 129, -1),
                                   DrillBlock(27, 40, 200, 3), DrillBlock(9, 5, 3, 0), DrillBlock(81, 2, 7, -1)])
def test_a_blocks_sectors_are_exactly_its_childrens(block):
    """Every sector in a block's ring and layer range has that block in its
    chain exactly when it is a sector of one of its children, down to the
    level-3 blocks."""
    m = block.m
    level = DRILL_LEVELS.index(m)
    via_children = _sectors_through_children(block)
    # A sector's block depends on its layer only through the slab, so test
    # each column once and the layers separately.
    half = (m - 1) // 2
    layers = range(block.slab * m - half, block.slab * m + half + 1)
    assert all(drill_chain_of(0, j, 0)[level].slab == block.slab for j in layers)
    assert drill_chain_of(0, layers[0] - 1, 0)[level].slab != block.slab
    assert drill_chain_of(0, layers[-1] + 1, 0)[level].slab != block.slab
    direct = set()
    for i in range(block.ring * m, block.ring * m + m):
        for k in range(ring_sector_count(i)):
            b = drill_chain_of(i, block.slab * m, k)[level]
            if (b.ring, b.wedge) == (block.ring, block.wedge):
                direct.update(DrillBlock(1, i, k, j) for j in layers)
    assert direct and via_children == direct


def test_level_three_blocks_hold_about_nine_sectors_a_layer():
    sizes = [len(drill_block_sectors(DrillBlock(3, ring, wedge, 0), 0))
             for ring in range(1, 1286, 37) for wedge in range(0, drill_wedge_count(3, ring), 97)]
    assert min(sizes) >= 1 and max(sizes) <= 15
    assert 8.0 <= sum(sizes) / len(sizes) <= 10.5


def test_keys_round_trip_and_reject_bad_blocks():
    block = DrillBlock(27, 63, 115, -1)
    assert format_drill_key(block) == "27.63.115.-1"
    assert parse_drill_key("27.63.115.-1") == block
    for bad in ("", None, "27.63.115", "1.0.0.0", "729.0.0.0", "27.-1.0.0", "27.0.99.0", "a.b.c.d", "27.0.0.0.0"):
        with pytest.raises(ValueError):
            parse_drill_key(bad)


def test_levels():
    assert DRILL_LEVELS == (243, 81, 27, 9, 3, 1)
    with pytest.raises(ValueError):
        drill_wedge_count(729, 0)


# --- The page's copy agrees --------------------------------------------------

node_only = pytest.mark.skipif(NODE is None, reason="node is not installed")


@node_only
def test_page_wedge_counts_match():
    rings = {m: list(_rings(m)) for m in (243, 81, 27, 9, 3)}
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
def test_page_galaxy_block_groups_match():
    got = _run(f"""
console.log(JSON.stringify({{
  galaxy: P.galaxyBlockGroups({EDGE_RING}, {EDGE_LAYER}).map(g => [g.slab, g.blocks.length, [g.blocks[0].m, g.blocks[0].ring, g.blocks[0].wedge, g.blocks[0].slab], [g.blocks.at(-1).m, g.blocks.at(-1).ring, g.blocks.at(-1).wedge, g.blocks.at(-1).slab]]),
  slabs: [P.galaxyBlockGroups(10), P.galaxyBlockGroups(10, 121), P.galaxyBlockGroups(10, 122)].map(gs => gs.map(g => g.slab)),
  formatted: P.formatDrillKey({{m: 3, ring: 568, wedge: 1036, slab: -7}}),
}}));
""")
    count = sum(drill_wedge_count(243, i) for i in range(EDGE_RING // 243 + 1))
    assert got["galaxy"] == [[slab, count, [243, 0, 0, slab], [243, EDGE_RING // 243, drill_wedge_count(243, EDGE_RING // 243) - 1, slab]]
                             for slab in (-1, 0, 1)]
    assert got["slabs"] == [[0], [0], [-1, 0, 1]]
    assert got["formatted"] == "3.568.1036.-7"
