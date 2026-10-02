"""
Tests for `galaxyGeometry.sectors_along_segment` (NAV.38), every sector a
straight segment passes through, and its browser twin
`galaxyprisms.sectorsAlongSegment` (run under node, skipped where node
isn't installed): along an axis, diagonally, through the core, out into
the halo, in a face and through an edge, against dense sampling, and the
two languages against each other.
"""

import json
import math
import os
import random
import shutil
import subprocess

import pytest

from stellarObjects.galaxyGeometry import (
    layer_bounds_pc, ring_bounds_pc, ring_radius_pc, ring_sector_count, sector_address_at,
    sector_position_pc, sectors_along_segment, slot_angle_bounds,
)

EDGE = 4.0
NODE = shutil.which("node")
MODULE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "html", "static", "galaxyprisms.js")


def _sampled(start, end, steps=20000):
    """The cells dense sampling of the segment finds."""
    cells = set()
    for k in range(steps + 1):
        t = k / steps
        point = tuple(a + t * (b - a) for a, b in zip(start, end))
        cells.add(sector_address_at(point, EDGE))
    return cells


def _adjacent(a, b):
    """Whether two cells of one walk share a face, edge or corner, as far
    as their bounds go (each consecutive pair of a walk must)."""
    (ra, la, sa), (rb, lb, sb) = a, b
    if abs(ra - rb) > 1 or abs(la - lb) > 1:
        return False
    a0, a1 = slot_angle_bounds(ra, sa)
    b0, b1 = slot_angle_bounds(rb, sb)
    tau = 2 * math.pi
    eps = 1e-9

    def overlap(p0, p1, q0, q1):
        return p0 <= q1 + eps and q0 <= p1 + eps

    return any(overlap(a0, a1, b0 + k * tau, b1 + k * tau) for k in (-1, 0, 1))


def test_a_point_segment_is_its_own_cell():
    point = (8000.0, 13.0, -7.0)
    assert sectors_along_segment(point, point, EDGE) == [sector_address_at(point, EDGE)]


def test_along_the_x_axis_every_ring_once_in_order():
    cells = sectors_along_segment((0.5, 0.0, 0.0), (100.0, 0.0, 0.0), EDGE)
    assert cells == [(ring, 0, 0) for ring in range(26)]


def test_along_the_x_axis_backwards_is_the_same_walk_reversed():
    cells = sectors_along_segment((100.0, 0.0, 0.0), (0.5, 0.0, 0.0), EDGE)
    assert cells == [(ring, 0, 0) for ring in reversed(range(26))]


def test_straight_up_into_the_halo_every_layer_once():
    x, y, _z = sector_position_pc(2000, 0, 7, EDGE)
    cells = sectors_along_segment((x, y, 0.0), (x, y, 3000.0), EDGE)
    assert cells == [(2000, layer, 7) for layer in range(0, 751)]


def test_segment_lying_in_a_layer_face_is_in_the_layer_above():
    _bottom, top = layer_bounds_pc(0, EDGE)
    cells = sectors_along_segment((10.0, 0.0, top), (10.0, 0.0, top), EDGE)
    assert {c[1] for c in sectors_along_segment((10.0, 1.0, top), (60.0, 1.0, top), EDGE)} == {1}
    assert cells[0][1] == 1
    under = math.nextafter(top, -math.inf)
    assert {c[1] for c in sectors_along_segment((10.0, 1.0, under), (60.0, 1.0, under), EDGE)} == {0}


def test_ending_exactly_on_a_ring_face_includes_the_next_ring():
    _inner, outer = ring_bounds_pc(10, EDGE)
    cells = sectors_along_segment((ring_radius_pc(10, EDGE), 0.0, 0.0), (outer, 0.0, 0.0), EDGE)
    assert cells == [(10, 0, 0), (11, 0, 0)]
    cells = sectors_along_segment((ring_radius_pc(10, EDGE), 0.0, 0.0), (math.nextafter(outer, 0.0), 0.0, 0.0), EDGE)
    assert cells == [(10, 0, 0)]


def test_around_a_ring_centerline_every_slot_crossed():
    """A chord inside one ring, from one slot's middle to another's, crosses
    each slot between them once, in turning order."""
    ring = 40
    n = ring_sector_count(ring)
    # A chord over a few slots at r = 162 pc sags well under the 2 pc to
    # the ring's inner face, so it stays in the ring.
    for first, last, slots in ((0, 3, [0, 1, 2, 3]), (5, 2, [5, 4, 3, 2]), (n - 2, 1, [n - 2, n - 1, 0, 1])):
        start = sector_position_pc(ring, 0, first, EDGE)
        end = sector_position_pc(ring, 0, last, EDGE)
        assert sectors_along_segment(start, end, EDGE) == [(ring, 0, slot) for slot in slots]


def test_diagonal_through_the_core_matches_dense_sampling():
    start, end = (-50.0, -3.0, -9.0), (50.0, 3.0, 9.0)
    cells = sectors_along_segment(start, end, EDGE)
    assert set(cells) == _sampled(start, end)
    assert any(c[0] == 0 for c in cells)
    assert cells[0] == sector_address_at(start, EDGE)
    assert cells[-1] == sector_address_at(end, EDGE)
    assert len(cells) == len(set(cells))


def test_through_the_axis_itself():
    """A line through the axis turns half a turn there without crossing a
    slot face: both ring-0 wedges it runs through are listed."""
    cells = sectors_along_segment((-20.0, 0.0, 1.0), (20.0, 0.0, 1.0), EDGE)
    assert set(cells) == _sampled((-20.0, 0.0, 1.0), (20.0, 0.0, 1.0))
    assert [c[0] for c in cells] == [5, 4, 3, 2, 1, 0, 0, 1, 2, 3, 4, 5]


def test_through_a_cell_edge_lists_the_cells_meeting_there():
    """A segment through the line where a ring face meets a layer face:
    the cells on both sides of both faces it actually touches."""
    inner, _outer = ring_bounds_pc(30, EDGE)
    _bottom, top = layer_bounds_pc(0, EDGE)
    start = (inner - 2.0, 0.5, top - 2.0)
    end = (inner + 2.0, 0.5, top + 2.0)
    mid = (inner, 0.5, top)
    cells = sectors_along_segment(start, end, EDGE)
    assert cells[0] == sector_address_at(start, EDGE)
    assert sector_address_at(mid, EDGE) in cells
    assert cells[-1] == sector_address_at(end, EDGE)
    for a, b in zip(cells, cells[1:]):
        assert _adjacent(a, b)


@pytest.mark.parametrize("seed", range(6))
def test_random_segments_match_dense_sampling(seed):
    rng = random.Random(seed)
    for _ in range(15):
        center = (rng.uniform(-300, 300), rng.uniform(-300, 300), rng.uniform(-20, 20))
        start = tuple(c + rng.uniform(-40, 40) for c in center)
        end = tuple(c + rng.uniform(-40, 40) for c in center)
        cells = sectors_along_segment(start, end, EDGE)
        assert len(cells) == len(set(cells))
        assert _sampled(start, end) <= set(cells)
        assert cells[0] == sector_address_at(start, EDGE)
        for a, b in zip(cells, cells[1:]):
            assert _adjacent(a, b), (start, end, a, b)


def test_out_into_the_halo_and_far_rings():
    start, end = (8000.0, 10.0, 0.0), (8600.0, -400.0, 900.0)
    cells = sectors_along_segment(start, end, EDGE)
    assert _sampled(start, end, steps=200000) <= set(cells)
    for a, b in zip(cells, cells[1:]):
        assert _adjacent(a, b)


@pytest.mark.skipif(NODE is None, reason="node is not installed")
def test_javascript_twin_agrees():
    rng = random.Random(38)
    segments = [
        ((0.5, 0.0, 0.0), (100.0, 0.0, 0.0)),
        ((-50.0, -3.0, -9.0), (50.0, 3.0, 9.0)),
        ((-20.0, 0.0, 1.0), (20.0, 0.0, 1.0)),
        ((8000.0, 10.0, 0.0), (8600.0, -400.0, 900.0)),
        ((10.0, 1.0, 2.0), (60.0, 1.0, 2.0)),
        ((3.0, 4.0, 5.0), (3.0, 4.0, 5.0)),
    ]
    for _ in range(40):
        center = (rng.uniform(-3000, 3000), rng.uniform(-3000, 3000), rng.uniform(-50, 50))
        segments.append(tuple(tuple(c + rng.uniform(-60, 60) for c in center) for _ in range(2)))
    source = (
        f"import * as P from {json.dumps('file://' + os.path.abspath(MODULE))};\n"
        f"const segs = {json.dumps(segments)};\n"
        f"console.log(JSON.stringify(segs.map(s => P.sectorsAlongSegment(s[0], s[1], {EDGE}))));\n"
    )
    result = subprocess.run([NODE, "--input-type=module", "-e", source], capture_output=True, text=True, check=True)
    got = json.loads(result.stdout)
    for (start, end), js in zip(segments, got):
        py = sectors_along_segment(start, end, EDGE)
        assert [(c["ring"], c["layer"], c["slot"]) for c in js] == py, (start, end)
