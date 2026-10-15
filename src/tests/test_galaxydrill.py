"""
Tests for the Galaxy Map's block ladder on the spherical grid:
`planetgen/galaxy/drill.py`. The ladder must nest -- every block inside exactly
one parent, a parent's sectors exactly its children's, every parent row edge
and wedge line a child's -- and every sector must belong to exactly one block
of each level.
"""

import math
import random
from collections import Counter

import pytest

from planetgen.galaxy import spherical
from planetgen.galaxy.drill import (
    DRILL_LEVELS, DRILL_STEP, DrillBlock, block_addresses, drill_band_count, drill_block_bounds, drill_block_sectors, drill_chain_of,
    drill_children, drill_parent, drill_wedge_count, format_drill_key, parse_drill_key,
)

EDGE = 4.0
OUTLINE_SHELL = 3764
"""The default galaxy's outermost shell (docs/design/spherical-sector-coordinates.md)."""


def _block_shells(m):
    return range(0, OUTLINE_SHELL // m + 1)


def _sample_sectors(count=40, seed=11):
    rng = random.Random(seed)
    shells = [1, 2, 3, 4, 5, 6, 9, 10, 26, 27, 28, 80, 242, 243, 244, 728, 729, OUTLINE_SHELL]
    shells += [rng.randrange(1, OUTLINE_SHELL + 1) for _ in range(count)]
    sectors = [(0, 0, 0)]
    for shell in shells:
        low, high = spherical.band_range(shell)
        for band in {low, high, 0, rng.randint(low, high), rng.randint(low, high)}:
            n = spherical.wedge_count(shell, band)
            for wedge in {0, n - 1, rng.randrange(n)}:
                sectors.append((shell, band, wedge))
    return sectors


def test_the_ladder_is_243_81_27_9_3_then_a_sector():
    assert DRILL_LEVELS == (243, 81, 27, 9, 3, 1)
    assert DRILL_STEP == 3


def test_a_block_shell_has_an_odd_number_of_rows_so_the_plane_runs_through_the_middle():
    for m in DRILL_LEVELS[:-1]:
        shells = _block_shells(m)
        for shell in shells[:: max(1, len(shells) // 200)]:
            assert drill_band_count(m, shell) % 2 == 1, (m, shell)


def test_the_top_counts_aim_at_m_edges_of_arc():
    for shell in (5, 40, 200, 1000, 3000):
        radius = shell * 243 + 121
        rows = drill_band_count(243, shell)
        assert rows == pytest.approx(math.pi * radius / 243, abs=1)
        equator = drill_wedge_count(243, shell, 0)
        assert equator == pytest.approx(2 * math.pi * radius / 243, rel=0.05, abs=1)


@pytest.mark.parametrize("m", DRILL_LEVELS[1:-1])
def test_rows_and_wedges_nest_exactly(m):
    """Every child row lies in one parent row and every child row's wedge count is a whole multiple of its
    parent row's, so no parent line is lost."""
    parent_m = m * DRILL_STEP
    shells = _block_shells(m)
    for shell in list(shells[:: max(1, len(shells) // 120)]) + [0, 1, 2, 3, 4]:
        parent_shell = shell // DRILL_STEP
        rows, parent_rows = drill_band_count(m, shell), drill_band_count(parent_m, parent_shell)
        assert rows % parent_rows == 0 and (rows // parent_rows) % 2 == 1, (m, shell)
        for band in range(-(rows // 2), rows // 2 + 1, max(1, rows // 40)):
            parent = drill_parent(DrillBlock(m, shell, band, 0))
            assert parent.band in range(-(parent_rows // 2), parent_rows // 2 + 1)
            assert drill_wedge_count(m, shell, band) % drill_wedge_count(parent_m, parent_shell, parent.band) == 0, (m, shell, band)


def test_a_sectors_chain_runs_from_the_top_to_the_sector():
    chain = drill_chain_of(1705, 120, 8000)
    assert [b.m for b in chain] == list(DRILL_LEVELS)
    assert chain[-1] == DrillBlock(1, 1705, 120, 8000)
    assert [b.shell for b in chain[:-1]] == [7, 21, 63, 189, 568]
    for parent, child in zip(chain, chain[1:]):
        assert drill_parent(child) == parent


def test_every_sector_chain_agrees_with_the_child_lists():
    for shell, band, wedge in _sample_sectors():
        chain = drill_chain_of(shell, band, wedge)
        assert chain[0].m == 243 and chain[-1] == DrillBlock(1, shell, band, wedge)
        for parent, child in zip(chain, chain[1:]):
            assert drill_parent(child) == parent
            assert child in [c for _s, group in drill_children(parent) for c in group], (parent, child)


def test_the_core_and_the_caps_join_wedge_zero_of_the_polar_and_equatorial_rows():
    core = drill_chain_of(0, 0, 0)
    assert all(b.shell == 0 and b.band == 0 and b.wedge == 0 for b in core[:-1])
    for shell in (4, 50, 600):
        high = spherical.band_range(shell)[1]
        north, south = drill_parent((1, shell, high, 0)), drill_parent((1, shell, -high, 0))
        assert north.wedge == south.wedge == 0
        rows = drill_band_count(3, shell // 3)
        assert north.band == rows // 2 and south.band == -(rows // 2)


@pytest.mark.parametrize("shells", [(3, 4, 5), (12, 13, 14), (30, 31, 32)])
def test_each_sector_is_in_exactly_one_level_three_block(shells):
    """The blocks' sector lists partition the shells: nothing missing, nothing twice, and the lists agree with the
    chain."""
    seen = Counter()
    blocks = set()
    for shell in shells:
        for cell in spherical.shell_cells(shell):
            blocks.add(drill_parent((1, *cell)))
    for block in blocks:
        for sector in drill_block_sectors(block):
            seen[sector] += 1
            assert drill_parent(sector) == block
    assert set(seen.values()) == {1}
    assert sum(seen.values()) == sum(spherical.cell_count(k) for k in shells)


def test_block_shell_zero_holds_the_core_and_the_first_shells():
    sectors = []
    rows = drill_band_count(3, 0)
    for row in range(-(rows // 2), rows // 2 + 1):
        for wedge in range(drill_wedge_count(3, 0, row)):
            sectors += drill_block_sectors(DrillBlock(3, 0, row, wedge))
    assert DrillBlock(1, 0, 0, 0) in sectors
    assert len(sectors) == 1 + spherical.cell_count(1) + spherical.cell_count(2)


def test_a_level_three_block_holds_about_a_cube_of_sectors():
    blocks = sorted({drill_parent((1, *cell)) for cell in spherical.shell_cells(120)})
    counts = sorted(len(drill_block_sectors(b)) for b in blocks[:: max(1, len(blocks) // 400)])
    assert counts[len(counts) // 2] in range(24, 34)
    assert counts[0] >= 1 and counts[-1] <= 60


def test_children_partition_their_parent_and_nest_in_its_bounds():
    rng = random.Random(3)
    for m in DRILL_LEVELS[:-1]:
        for _ in range(12):
            shell = rng.randrange(0, OUTLINE_SHELL // m + 1)
            rows = drill_band_count(m, shell)
            band = rng.randint(-(rows // 2), rows // 2)
            block = DrillBlock(m, shell, band, rng.randrange(drill_wedge_count(m, shell, band)))
            (r0, r1), (c0, c1), (l0, l1) = drill_block_bounds(block, EDGE)
            children = [c for _s, group in drill_children(block) for c in group]
            assert len(set(children)) == len(children)
            if m == 3:  # a polar or thin block may hold no sector of its own
                continue
            assert children
            for child in children:
                assert drill_parent(child) == block
                (cr0, cr1), (cc0, cc1), (cl0, cl1) = drill_block_bounds(child, EDGE)
                assert r0 - 1e-9 <= cr0 and cr1 <= r1 + 1e-9
                assert c0 - 1e-12 <= cc0 and cc1 <= c1 + 1e-12
                assert l0 - 1e-12 <= cl0 and cl1 <= l1 + 1e-12


def test_a_block_cuts_the_sphere_in_whole_rows_and_wedges():
    """The rows of a block shell tile colatitude 0 to pi, and the wedges of a row tile the circle."""
    for m, shell in ((243, 4), (81, 40), (27, 300), (9, 90), (3, 1000)):
        rows = drill_band_count(m, shell)
        edges = [drill_block_bounds(DrillBlock(m, shell, b, 0), EDGE)[1] for b in range(rows // 2, -(rows // 2) - 1, -1)]
        assert edges[0][0] == 0 and edges[-1][1] == pytest.approx(math.pi)
        assert all(a[1] == pytest.approx(b[0]) for a, b in zip(edges, edges[1:]))
        band = rows // 3
        wedges = drill_wedge_count(m, shell, band)
        assert drill_block_bounds(DrillBlock(m, shell, band, wedges - 1), EDGE)[2][1] == pytest.approx(2 * math.pi)


def test_keys_round_trip_and_bad_ones_are_refused():
    block = DrillBlock(81, 12, -3, 40)
    assert format_drill_key(block) == "81.12.-3.40"
    assert parse_drill_key("81.12.-3.40") == block
    for bad in ("", "81.12.-3", "x.1.2.3", "729.0.0.0", "1.0.0.0", "81.-1.0.0", "81.12.999.0", "81.12.0.999999"):
        with pytest.raises(ValueError):
            parse_drill_key(bad)


@pytest.mark.parametrize("m", [3, 9, 27])
def test_the_addresses_of_a_block_are_its_sectors_and_each_names_the_block(m):
    block = drill_chain_of(40, 3, 21)[DRILL_LEVELS.index(m)]
    addresses = list(block_addresses(block))
    assert addresses and len(addresses) == len(set(addresses))
    assert (40, 3, 21) in addresses
    assert all(drill_chain_of(*address)[DRILL_LEVELS.index(m)] == block for address in addresses)
    assert [a for a in block_addresses(DrillBlock(1, 40, 3, 21))] == [(40, 3, 21)]
    if m > DRILL_STEP:
        assert sorted(addresses) == sorted(
            a for _s, group in drill_children(block) for child in group for a in block_addresses(child))
