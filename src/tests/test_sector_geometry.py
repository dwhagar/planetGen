# tests/test_sector_geometry.py

"""
Tests for `stellarObjects.sectorGeometry` -- cube orientation, and the
exact local same-shell Voronoi tessellation (`local_lateral_cell`) each
sector's vertices are built from, extruded radially (`prism_vertices`).

The headline property under test is exactness: two real same-shell
neighbors must compute *identical* shared vertices independently (see the
module's own docstring for why this follows from computing genuine 3D
circumcenters rather than nudging two guesses toward each other), and the
fast same-shell neighbor search (`_same_shell_neighbors`) must agree with
a brute-force, unconditionally-correct reference search across a wide
sample of shells and polar angles -- this is what was actually validated
by hand during development (840 cases, zero mismatches) before this
module was written; the parametrized test below re-runs a slice of that
same cross-check so a future change can't silently regress it.
"""

import math

import pytest

from stellarObjects.galaxyGeometry import (
    galactic_radius_pc, sector_position_pc, shell_sector_count,
    enumerate_sectors_within_radius,
)
from stellarObjects.sectorGeometry import (
    cube_orientation,
    local_lateral_cell,
    prism_vertices,
    _same_shell_neighbors,
    _circumcenter_3d,
)

EDGE_PC = 3.526  # ~DEFAULT_SECTOR_EDGE_LY (11.5 ly) converted to parsecs.


def _dist(a, b):
    return math.sqrt(sum((ac - bc) ** 2 for ac, bc in zip(a, b)))


def _brute_force_same_shell_neighbors(shell_index, shell_slot_index, edge_pc, radius_pc):
    position = sector_position_pc(shell_index, shell_slot_index, edge_pc)
    return {
        slot for k2, slot, x, y, z, dist in
        enumerate_sectors_within_radius(position, radius_pc, edge_pc)
        if k2 == shell_index and slot != shell_slot_index
    }


# ---------------------------------------------------------------------------
# cube_orientation
# ---------------------------------------------------------------------------

def test_cube_orientation_is_orthonormal_and_right_handed():
    position = sector_position_pc(50, 12000, EDGE_PC)
    local_x, local_y, local_z = cube_orientation(position)

    for axis in (local_x, local_y, local_z):
        assert math.sqrt(sum(c * c for c in axis)) == pytest.approx(1.0)

    def dot(a, b):
        return sum(ac * bc for ac, bc in zip(a, b))

    assert dot(local_x, local_y) == pytest.approx(0.0, abs=1e-9)
    assert dot(local_y, local_z) == pytest.approx(0.0, abs=1e-9)
    assert dot(local_x, local_z) == pytest.approx(0.0, abs=1e-9)

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
    for actual, expected in zip(local_z, tuple(c / r for c in position)):
        assert actual == pytest.approx(expected)


def test_cube_orientation_handles_near_pole_degeneracy():
    position = (0.0, 0.0, 1000.0)
    local_x, local_y, local_z = cube_orientation(position)
    for axis in (local_x, local_y, local_z):
        for c in axis:
            assert not math.isnan(c)
    assert local_z == pytest.approx((0.0, 0.0, 1.0))
    assert local_x == pytest.approx((1.0, 0.0, 0.0))


# ---------------------------------------------------------------------------
# _circumcenter_3d
# ---------------------------------------------------------------------------

def test_circumcenter_3d_matches_hand_computed_right_triangle():
    # A right triangle's circumcenter is the midpoint of its hypotenuse.
    p, q, r = (0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0)
    cc = _circumcenter_3d(p, q, r)
    assert cc == pytest.approx((1.0, 1.0, 0.0))
    # Equidistant from all three, by construction.
    assert _dist(cc, p) == pytest.approx(_dist(cc, q))
    assert _dist(cc, q) == pytest.approx(_dist(cc, r))


def test_circumcenter_3d_returns_none_for_collinear_points():
    p, q, r = (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)
    assert _circumcenter_3d(p, q, r) is None


# ---------------------------------------------------------------------------
# _same_shell_neighbors: fast path vs. brute-force ground truth
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("shell_index,phi_frac", [
    (0, 0.5), (1, 0.5), (5, 0.1), (5, 0.9),
    (50, 0.05), (50, 0.5), (50, 0.95),
    (500, 0.02), (500, 0.5), (500, 0.98),
    (3000, 0.01), (3000, 0.5), (3000, 0.99),
    (3600, 0.3), (3600, 0.5), (3600, 0.7),
])
def test_same_shell_neighbors_matches_brute_force(shell_index, phi_frac):
    n_k = shell_sector_count(shell_index)
    slot = min(n_k - 1, max(0, int(phi_frac * n_k)))
    radius_pc = 3.0 * EDGE_PC

    fast = {s for s, _pos in _same_shell_neighbors(shell_index, slot, EDGE_PC)
            if _dist(_pos, sector_position_pc(shell_index, slot, EDGE_PC)) <= radius_pc}
    expected = _brute_force_same_shell_neighbors(shell_index, slot, EDGE_PC, radius_pc)
    assert fast == expected


# ---------------------------------------------------------------------------
# local_lateral_cell / prism_vertices
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("shell_index,shell_slot_index", [
    (0, 0), (1, 5), (5, 100), (50, 12000), (500, 500000), (3000, 50000000),
])
def test_local_lateral_cell_has_plausible_vertex_count(shell_index, shell_slot_index):
    cell = local_lateral_cell(shell_index, shell_slot_index, EDGE_PC)
    # Real Voronoi cells on a near-uniform point set concentrate on 4-8
    # sides (mean 6, per Euler's formula); shell 0's degenerate 3-sector
    # case can differ.
    assert 3 <= len(cell) <= 10


@pytest.mark.parametrize("shell_index,slot_a", [
    (1, 0), (1, 12), (5, 100), (50, 12000),
])
def test_local_lateral_cell_shares_exact_vertices_with_a_real_neighbor(shell_index, slot_a):
    cell_a = local_lateral_cell(shell_index, slot_a, EDGE_PC)

    # Find slot_a's nearest same-shell neighbor and confirm cell_b (computed
    # entirely independently, from slot_b's own perspective) shares at
    # least one vertex with cell_a to within floating-point noise. Small,
    # sparse shells (n_k small enough to use the brute-force candidate
    # path) are included deliberately: a past bug in the half-plane
    # bisector formula used after gnomonic projection was only wrong for
    # widely-separated candidates, so it passed unnoticed with only
    # large/dense shells (shell 50+) under test.
    neighbors = _same_shell_neighbors(shell_index, slot_a, EDGE_PC)
    position_a = sector_position_pc(shell_index, slot_a, EDGE_PC)
    slot_b, _pos_b = min(neighbors, key=lambda item: _dist(item[1], position_a))
    cell_b = local_lateral_cell(shell_index, slot_b, EDGE_PC)

    best = min(_dist(va, vb) for va in cell_a for vb in cell_b)
    assert best < 1e-6  # floating-point noise, not an approximation residual


def test_local_lateral_cell_fully_tiles_a_small_shell_with_no_orphan_vertices():
    # Regression test for a bug where the half-plane test used after
    # gnomonic projection was the *flat*-plane bisector of (0, 0) and
    # (u, v) -- a good approximation only when every candidate is close
    # (large/dense shells), but wrong enough on a small, sparse, fully
    # populated shell (every sector present, so every vertex must be
    # shared with at least one other sector's cell) that it silently
    # produced extra, unshared "orphan" vertices: two non-adjacent
    # sectors' cells met at a point a third, genuinely closer sector
    # should have cut off first.
    shell_index = 1
    n_k = shell_sector_count(shell_index)
    cells = {i: local_lateral_cell(shell_index, i, EDGE_PC) for i in range(n_k)}

    for i, cell in cells.items():
        for v in cell:
            shared = any(
                any(_dist(v, ov) < 1e-6 for ov in other_cell)
                for j, other_cell in cells.items() if j != i
            )
            assert shared, f"sector {i}'s vertex {v} is not shared with any other same-shell sector"


def test_prism_vertices_inner_and_outer_sit_on_the_right_spheres():
    shell_index, slot = 50, 12000
    result = prism_vertices(shell_index, slot, EDGE_PC)
    assert len(result["inner"]) == len(result["outer"])
    assert len(result["inner"]) >= 3

    inner_radius = shell_index * EDGE_PC
    outer_radius = (shell_index + 1) * EDGE_PC
    for v in result["inner"]:
        assert galactic_radius_pc(v) == pytest.approx(inner_radius)
    for v in result["outer"]:
        assert galactic_radius_pc(v) == pytest.approx(outer_radius)


def test_prism_vertices_shell_zero_inner_collapses_to_origin():
    result = prism_vertices(0, 0, EDGE_PC)
    for v in result["inner"]:
        assert v == (0.0, 0.0, 0.0)
    for v in result["outer"]:
        assert galactic_radius_pc(v) == pytest.approx(EDGE_PC)


def test_prism_vertices_is_deterministic():
    a = prism_vertices(50, 12000, EDGE_PC)
    b = prism_vertices(50, 12000, EDGE_PC)
    assert a == b


def test_local_lateral_cell_radially_adjacent_shells_share_outer_radius():
    # Shell k's outer bound and shell k+1's inner bound are the same
    # sphere -- confirms the "area matches, vertices needn't" radial
    # design directly on real output.
    shell_index, slot = 50, 12000
    outer_of_k = prism_vertices(shell_index, slot, EDGE_PC)["outer"]
    inner_of_k_plus_1 = prism_vertices(shell_index + 1, slot % shell_sector_count(shell_index + 1), EDGE_PC)["inner"]
    for v in outer_of_k:
        assert galactic_radius_pc(v) == pytest.approx((shell_index + 1) * EDGE_PC)
    for v in inner_of_k_plus_1:
        assert galactic_radius_pc(v) == pytest.approx((shell_index + 1) * EDGE_PC)
