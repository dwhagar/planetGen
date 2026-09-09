# tests/test_sector_geometry.py

"""
Tests for `stellarObjects.sectorGeometry` -- cube orientation, naive
8-vertex corners, and the neighbor-relaxation pass that nudges corners
toward matching ones on nearby sectors (see that module's own docstring
for what "relaxed" does and does not guarantee).
"""

import math

import pytest

from stellarObjects.galaxyGeometry import galactic_radius_pc, sector_position_pc
from stellarObjects.sectorGeometry import (
    cube_orientation,
    naive_cube_vertices,
    relax_vertices,
)

EDGE_PC = 3.526  # ~DEFAULT_SECTOR_EDGE_LY (11.5 ly) converted to parsecs.


# ---------------------------------------------------------------------------
# cube_orientation
# ---------------------------------------------------------------------------

def test_cube_orientation_is_orthonormal_and_right_handed():
    position = sector_position_pc(50, 12000, EDGE_PC)
    local_x, local_y, local_z = cube_orientation(position)

    for axis in (local_x, local_y, local_z):
        length = math.sqrt(sum(c * c for c in axis))
        assert length == pytest.approx(1.0)

    def dot(a, b):
        return sum(ac * bc for ac, bc in zip(a, b))

    assert dot(local_x, local_y) == pytest.approx(0.0, abs=1e-9)
    assert dot(local_y, local_z) == pytest.approx(0.0, abs=1e-9)
    assert dot(local_x, local_z) == pytest.approx(0.0, abs=1e-9)

    # Right-handed: local_x cross local_y == local_z.
    cross = (
        local_x[1] * local_y[2] - local_x[2] * local_y[1],
        local_x[2] * local_y[0] - local_x[0] * local_y[2],
        local_x[0] * local_y[1] - local_x[1] * local_y[0],
    )
    for c, z in zip(cross, local_z):
        assert c == pytest.approx(z, abs=1e-9)


def test_cube_orientation_local_z_points_radially_outward():
    position = sector_position_pc(50, 12000, EDGE_PC)
    _, _, local_z = cube_orientation(position)
    r = galactic_radius_pc(position)
    expected = tuple(c / r for c in position)
    for actual, exp in zip(local_z, expected):
        assert actual == pytest.approx(exp)


def test_cube_orientation_handles_near_pole_degeneracy():
    # A position whose radial direction is (numerically) galactic north --
    # the coordinate doc's own flagged edge case. Must not raise or return
    # NaNs; falls back to projecting galaxy +X instead.
    position = (0.0, 0.0, 1000.0)
    local_x, local_y, local_z = cube_orientation(position)
    for axis in (local_x, local_y, local_z):
        for c in axis:
            assert not math.isnan(c)
    assert local_z == pytest.approx((0.0, 0.0, 1.0))
    assert local_x == pytest.approx((1.0, 0.0, 0.0))


# ---------------------------------------------------------------------------
# naive_cube_vertices
# ---------------------------------------------------------------------------

def test_naive_cube_vertices_returns_eight_points():
    position = sector_position_pc(50, 12000, EDGE_PC)
    vertices = naive_cube_vertices(position, EDGE_PC)
    assert len(vertices) == 8
    assert len(set(vertices)) == 8  # all distinct


def test_naive_cube_vertices_centroid_is_the_center():
    position = sector_position_pc(50, 12000, EDGE_PC)
    vertices = naive_cube_vertices(position, EDGE_PC)
    centroid = tuple(sum(v[axis] for v in vertices) / 8 for axis in range(3))
    for actual, expected in zip(centroid, position):
        assert actual == pytest.approx(expected)


def test_naive_cube_vertices_has_correct_edge_length_and_volume():
    position = sector_position_pc(50, 12000, EDGE_PC)
    vertices = naive_cube_vertices(position, EDGE_PC)

    def dist(a, b):
        return math.sqrt(sum((ac - bc) ** 2 for ac, bc in zip(a, b)))

    # Corner 0 (---) to corner 1 (+--) differs only in local x -- one edge.
    assert dist(vertices[0], vertices[1]) == pytest.approx(EDGE_PC)
    # Corner 0 (---) to corner 7 (+++) is the cube's space diagonal.
    assert dist(vertices[0], vertices[7]) == pytest.approx(EDGE_PC * math.sqrt(3))


def test_naive_cube_vertices_at_origin_axis_aligned():
    # A degenerate position still produces a valid, correctly-sized cube
    # (cube_orientation's own fallback axes), just not physically meaningful
    # as a real sector (never reached by sector_position_pc in practice).
    vertices = naive_cube_vertices((0.0, 0.0, 0.0), EDGE_PC)
    assert len(vertices) == 8
    half = EDGE_PC / 2.0
    expected_signs = [(-1, -1, -1), (1, -1, -1), (-1, 1, -1), (1, 1, -1),
                       (-1, -1, 1), (1, -1, 1), (-1, 1, 1), (1, 1, 1)]
    for vertex, signs in zip(vertices, expected_signs):
        for c, s in zip(vertex, signs):
            assert c == pytest.approx(s * half)


# ---------------------------------------------------------------------------
# relax_vertices
# ---------------------------------------------------------------------------

def test_relax_vertices_returns_eight_points_and_stays_near_naive():
    position = sector_position_pc(50, 12000, EDGE_PC)
    naive = naive_cube_vertices(position, EDGE_PC)
    relaxed = relax_vertices(50, 12000, EDGE_PC)

    assert len(relaxed) == 8
    # Relaxation only ever averages with corners within CORNER_MATCH_RADIUS_FACTOR
    # * edge_pc of the naive one, so it can never move a vertex further than
    # that from its naive position.
    for n, r in zip(naive, relaxed):
        d = math.sqrt(sum((nc - rc) ** 2 for nc, rc in zip(n, r)))
        assert d <= 0.5 * EDGE_PC + 1e-9


def test_relax_vertices_is_deterministic():
    a = relax_vertices(50, 12000, EDGE_PC)
    b = relax_vertices(50, 12000, EDGE_PC)
    assert a == b


def test_relax_vertices_shrinks_gap_between_two_actual_neighbors():
    # Find a real pair of nearest-neighbor addresses via a small brute-force
    # scan of nearby shells, then confirm relaxation pulls their nearest
    # naive corners closer together than they started (or leaves an exact
    # match untouched) -- the concrete, sector-scale version of the module
    # docstring's "shrinks, does not eliminate" claim.
    from stellarObjects.galaxyGeometry import shell_sector_count, shell_radius_pc

    shell_index = 50
    n_k = shell_sector_count(shell_index)
    base_slot = 12000
    base_position = sector_position_pc(shell_index, base_slot, EDGE_PC)

    def dist(a, b):
        return math.sqrt(sum((ac - bc) ** 2 for ac, bc in zip(a, b)))

    nearest_slot = None
    nearest_dist = None
    for slot in range(n_k):
        if slot == base_slot:
            continue
        candidate = sector_position_pc(shell_index, slot, EDGE_PC)
        d = dist(base_position, candidate)
        if nearest_dist is None or d < nearest_dist:
            nearest_dist = d
            nearest_slot = slot

    naive_a = naive_cube_vertices(base_position, EDGE_PC)
    naive_b = naive_cube_vertices(sector_position_pc(shell_index, nearest_slot, EDGE_PC), EDGE_PC)

    def min_corner_gap(vertices_a, vertices_b):
        return min(dist(va, vb) for va in vertices_a for vb in vertices_b)

    naive_gap = min_corner_gap(naive_a, naive_b)

    relaxed_a = relax_vertices(shell_index, base_slot, EDGE_PC)
    relaxed_b = relax_vertices(shell_index, nearest_slot, EDGE_PC)
    relaxed_gap = min_corner_gap(relaxed_a, relaxed_b)

    assert relaxed_gap <= naive_gap + 1e-9
