"""
Tests for `spherical.sectors_along_segment` (NAV.38), every sector a
straight segment passes through: along an axis, diagonally, through the core,
over a pole, out into the halo, in a face and through an edge, against dense
sampling, and each consecutive pair of a walk against the cells' neighbors.

The browser twin of the old cylinder walk went with the cylinder grid
(MAP.206).
"""

import math
import random

import pytest

from planetgen.galaxy import spherical
from planetgen.galaxy.spherical import (
    band_range, neighbor_addresses, sector_address_at, sector_position_pc, sectors_along_segment, shell_bounds_pc,
    wedge_count,
)

EDGE = 4.0


def _sampled(start, end, steps=20000):
    """The cells dense sampling of the segment finds."""
    cells = set()
    for k in range(steps + 1):
        t = k / steps
        point = tuple(a + t * (b - a) for a, b in zip(start, end))
        cells.add(sector_address_at(point, EDGE))
    return cells


def _adjacent(a, b):
    """Whether two cells of one walk touch: the same cell's neighbor across a face, or one
    that shares only an edge or a corner (a walk can pass through those), so the cells' shells
    are at most one apart and their bands and wedges overlap or touch."""
    if abs(a[0] - b[0]) > 1:
        return False
    if b in neighbor_addresses(*a):
        return True
    # An edge or corner: some neighbor of `a` is a neighbor of `b` (the walk went through the line between).
    return bool(set(neighbor_addresses(*a)) & set(neighbor_addresses(*b)))


def test_a_point_segment_is_its_own_cell():
    point = (8000.0, 13.0, -7.0)
    assert sectors_along_segment(point, point, EDGE) == [sector_address_at(point, EDGE)]


def test_along_the_x_axis_every_shell_once_in_order():
    cells = sectors_along_segment((0.5, 0.0, 0.0), (100.0, 0.0, 0.0), EDGE)
    assert cells == [(shell, 0, 0) for shell in range(26)]


def test_along_the_x_axis_backwards_is_the_same_walk_reversed():
    cells = sectors_along_segment((100.0, 0.0, 0.0), (0.5, 0.0, 0.0), EDGE)
    assert cells == [(shell, 0, 0) for shell in reversed(range(26))]


def test_straight_up_the_galactic_axis_goes_core_then_the_north_cap_of_every_shell():
    cells = sectors_along_segment((0.0, 0.0, 0.5), (0.0, 0.0, 100.0), EDGE)
    assert cells == [(0, 0, 0)] + [(shell, band_range(shell)[1], 0) for shell in range(1, 26)]


def test_a_steep_line_through_the_core_comes_out_of_a_south_cap_and_into_a_north_one():
    """The segment leaves the galaxy's south pole region and enters the north one: it lists cells of both caps."""
    cells = sectors_along_segment((20.0, 0.0, -30.0), (-20.0, 0.0, 30.0), EDGE)
    assert set(cells) == _sampled((20.0, 0.0, -30.0), (-20.0, 0.0, 30.0))
    assert any(spherical.is_cap(cell[0], cell[1]) and cell[1] < 0 for cell in cells)
    assert any(spherical.is_cap(cell[0], cell[1]) and cell[1] > 0 for cell in cells)


def test_ending_exactly_on_a_shell_face_includes_the_next_shell():
    _inner, outer = shell_bounds_pc(10, EDGE)
    cells = sectors_along_segment((10 * EDGE, 0.0, 0.0), (outer, 0.0, 0.0), EDGE)
    assert cells == [(10, 0, 0), (11, 0, 0)]
    cells = sectors_along_segment((10 * EDGE, 0.0, 0.0), (math.nextafter(outer, 0.0), 0.0, 0.0), EDGE)
    assert cells == [(10, 0, 0)]


def test_around_a_shell_centerline_every_wedge_crossed():
    """A chord inside one band, from one wedge's middle to another's, crosses
    each wedge between them once, in turning order."""
    shell = 40
    n = wedge_count(shell, 0)
    for first, last, wedges in ((0, 3, [0, 1, 2, 3]), (5, 2, [5, 4, 3, 2]), (n - 2, 1, [n - 2, n - 1, 0, 1])):
        start = sector_position_pc(shell, 0, first, EDGE)
        end = sector_position_pc(shell, 0, last, EDGE)
        assert sectors_along_segment(start, end, EDGE) == [(shell, 0, wedge) for wedge in wedges]


def test_diagonal_through_the_core_matches_dense_sampling():
    start, end = (-50.0, -3.0, -9.0), (50.0, 3.0, 9.0)
    cells = sectors_along_segment(start, end, EDGE)
    assert _sampled(start, end) <= set(cells)
    assert (0, 0, 0) in cells
    assert cells[0] == sector_address_at(start, EDGE)
    assert cells[-1] == sector_address_at(end, EDGE)
    assert len(cells) == len(set(cells))


def test_through_the_core_the_shells_run_down_then_up():
    cells = sectors_along_segment((-20.0, 0.0, 1.0), (20.0, 0.0, 1.0), EDGE)
    shells = [cell[0] for cell in cells]
    assert shells[0] == shells[-1] == 5 and min(shells) == 0
    low = shells.index(0)
    assert shells[:low + 1] == sorted(shells[:low + 1], reverse=True)
    assert shells[low:] == sorted(shells[low:])


def test_through_a_cell_edge_lists_the_cells_meeting_there():
    """A segment through the line where a shell face meets a band face: the
    cells on both sides of both faces it actually touches."""
    shell = 30
    inner, _outer = shell_bounds_pc(shell, EDGE)
    layout = spherical.shell_layout(shell)
    colatitude = layout.edge(layout.band_count // 2)    # the face between two bands of the shell
    point = (inner * math.sin(colatitude), 0.0, inner * math.cos(colatitude))
    along = (math.sin(colatitude), 0.0, math.cos(colatitude))
    start = tuple(p - 2.0 * u for p, u in zip(point, along))
    end = tuple(p + 2.0 * u for p, u in zip(point, along))
    cells = sectors_along_segment(start, end, EDGE)
    assert cells[0] == sector_address_at(start, EDGE)
    assert cells[-1] == sector_address_at(end, EDGE)
    assert sector_address_at(point, EDGE) in cells
    for a, b in zip(cells, cells[1:]):
        assert _adjacent(a, b)


@pytest.mark.parametrize("seed", range(6))
def test_random_segments_match_dense_sampling(seed):
    rng = random.Random(seed)
    for _ in range(15):
        center = (rng.uniform(-300, 300), rng.uniform(-300, 300), rng.uniform(-300, 300))
        start = tuple(c + rng.uniform(-40, 40) for c in center)
        end = tuple(c + rng.uniform(-40, 40) for c in center)
        cells = sectors_along_segment(start, end, EDGE)
        assert len(cells) == len(set(cells))
        assert _sampled(start, end, 6000) <= set(cells)
        assert cells[0] == sector_address_at(start, EDGE)
        for a, b in zip(cells, cells[1:]):
            assert _adjacent(a, b), (start, end, a, b)


def test_out_into_the_halo_and_far_shells():
    start, end = (8000.0, 10.0, 0.0), (8600.0, -400.0, 900.0)
    cells = sectors_along_segment(start, end, EDGE)
    assert _sampled(start, end, steps=100000) <= set(cells)
    for a, b in zip(cells, cells[1:]):
        assert _adjacent(a, b)
