# tests/test_galaxy_geometry.py

"""
Tests for `stellarObjects.galaxyGeometry` -- the shell/Fibonacci-sphere
radial tiling scheme (`docs/design/galaxy-coordinate-system.md` section 3)
and the radius-based neighborhood-enumeration primitive built on top of it
(section 8, `enumerate_sectors_within_radius`).

The enumeration tests are deliberately deterministic and checkable, not
statistical: `test_enumerate_matches_brute_force_over_several_shells`
cross-checks the pruned algorithm's output against an independent
brute-force scan of every slot in the same shells (same underlying
`sector_position_pc`, different traversal logic -- the "same ground truth
function, independently reached two ways" cross-check confirms the
pruning logic includes/excludes exactly the same slots brute force would,
not that the position formula itself is correct). The other tests use
small, hand-computable shells (`shell_index=0`, 3 slots) with explicitly
worked-out expected values.
"""

import math

import pytest

from stellarObjects.galaxyGeometry import (
    GOLDEN_RATIO,
    _candidate_shell_range,
    slot_index_bounds_for_phi_range,
    enumerate_sectors_within_radius,
    galactic_radius_pc,
    sector_position_pc,
    sector_wedge_vertices_pc,
    shell_radius_pc,
    shell_sector_count,
)

EDGE_PC = 3.526  # ~DEFAULT_SECTOR_EDGE_LY (11.5 ly) converted to parsecs.


# ---------------------------------------------------------------------------
# Shell math -- cross-checked against docs/design/galaxy-coordinate-system.md
# section 3's own worked table.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("shell_index,expected_n_k", [
    (0, 3), (1, 28), (2, 79), (3, 154), (5, 380), (10, 1385),
])
def test_shell_sector_count_matches_design_doc_table(shell_index, expected_n_k):
    assert shell_sector_count(shell_index) == expected_n_k


def test_shell_radius_pc_is_nominal_midpoint_of_shell_thickness():
    # r_k = (k + 0.5) * edge_pc -- shell 0 spans [0, edge_pc), nominal
    # radius is exactly half the edge length.
    assert shell_radius_pc(0, EDGE_PC) == pytest.approx(0.5 * EDGE_PC)
    assert shell_radius_pc(49, EDGE_PC) == pytest.approx(49.5 * EDGE_PC)


def test_sector_position_pc_matches_design_doc_worked_example():
    # docs/design/galaxy-coordinate-system.md section 6's worked example:
    # shell k=50, slot i=12000 (of N_50=32047), edge_pc=3.526.
    n_50 = shell_sector_count(50)
    assert n_50 == 32047

    x, y, z = sector_position_pc(50, 12000, 3.526)
    assert galactic_radius_pc((x, y, z)) == pytest.approx(178.06, abs=0.05)
    assert x == pytest.approx(-144.4, abs=0.5)
    assert y == pytest.approx(94.2, abs=0.5)
    assert z == pytest.approx(44.6, abs=0.5)


def test_sector_position_pc_rejects_out_of_range_slot_index():
    n_0 = shell_sector_count(0)
    with pytest.raises(ValueError):
        sector_position_pc(0, n_0, EDGE_PC)  # exactly one past the end
    with pytest.raises(ValueError):
        sector_position_pc(0, -1, EDGE_PC)


# ---------------------------------------------------------------------------
# sector_wedge_vertices_pc -- the sector's approximate on-shell cell shape,
# for the "proper 3D shape" sector-map outline (html/lib/starmap.py).
# ---------------------------------------------------------------------------

def test_sector_wedge_vertices_pc_returns_8_points_at_the_two_radial_bounds():
    shell_index, slot_index = 50, 12000
    r_k = shell_radius_pc(shell_index, EDGE_PC)
    r_min, r_max = r_k - EDGE_PC / 2, r_k + EDGE_PC / 2

    vertices = sector_wedge_vertices_pc(shell_index, slot_index, EDGE_PC)
    assert len(vertices) == 8

    radii = sorted(galactic_radius_pc(v) for v in vertices)
    for radius in radii[:4]:
        assert radius == pytest.approx(r_min, abs=1e-6)
    for radius in radii[4:]:
        assert radius == pytest.approx(r_max, abs=1e-6)


def test_sector_wedge_vertices_pc_is_centered_near_the_sector_position():
    # The wedge's own angular half-widths are small for a shell this
    # populated (n_50 = 32047), so its 8 vertices should average out
    # close to the exact center point `sector_position_pc` returns --
    # loosely (this is an approximate patch, not an exact one), but not
    # off by anything close to a full sector edge.
    shell_index, slot_index = 50, 12000
    center = sector_position_pc(shell_index, slot_index, EDGE_PC)
    vertices = sector_wedge_vertices_pc(shell_index, slot_index, EDGE_PC)

    mean = tuple(sum(v[axis] for v in vertices) / 8 for axis in range(3))
    for axis in range(3):
        assert mean[axis] == pytest.approx(center[axis], abs=EDGE_PC / 2)


def test_sector_wedge_vertices_pc_rejects_out_of_range_slot_index():
    n_0 = shell_sector_count(0)
    with pytest.raises(ValueError):
        sector_wedge_vertices_pc(0, n_0, EDGE_PC)
    with pytest.raises(ValueError):
        sector_wedge_vertices_pc(0, -1, EDGE_PC)


# ---------------------------------------------------------------------------
# Hand-computed small case: shell 0 has exactly 3 slots. Position formula,
# worked by hand:
#   phi_i   = acos(1 - 2*(i+0.5)/3)
#   theta_i = (2*pi*i / golden_ratio) mod 2*pi
# ---------------------------------------------------------------------------

def _hand_computed_shell_0_positions(edge_pc):
    r_0 = 0.5 * edge_pc
    positions = []
    for i in range(3):
        phi = math.acos(1 - 2 * (i + 0.5) / 3)
        theta = (2 * math.pi * i / GOLDEN_RATIO) % (2 * math.pi)
        x = r_0 * math.sin(phi) * math.cos(theta)
        y = r_0 * math.sin(phi) * math.sin(theta)
        z = r_0 * math.cos(phi)
        positions.append((x, y, z))
    return positions


def test_sector_position_pc_matches_hand_computed_shell_0():
    expected = _hand_computed_shell_0_positions(EDGE_PC)
    for i, expected_position in enumerate(expected):
        actual = sector_position_pc(0, i, EDGE_PC)
        assert actual == pytest.approx(expected_position)


def test_enumerate_from_origin_over_shell_0_matches_hand_computed_radius_filter():
    """A deterministic, hand-checkable case: with P at the origin and R
    exactly shell 0's own radius, every one of shell 0's 3 known,
    hand-computed slots must be returned (all sit at exactly r_0 from the
    origin); with R set just under that radius, none should be."""
    r_0 = shell_radius_pc(0, EDGE_PC)
    expected_positions = _hand_computed_shell_0_positions(EDGE_PC)

    included = list(enumerate_sectors_within_radius((0.0, 0.0, 0.0), r_0, EDGE_PC))
    included_addresses = {(shell, slot) for shell, slot, *_ in included}
    assert included_addresses == {(0, 0), (0, 1), (0, 2)}
    for shell, slot, x, y, z, dist in included:
        assert (x, y, z) == pytest.approx(expected_positions[slot])
        assert dist == pytest.approx(r_0)

    excluded = list(enumerate_sectors_within_radius((0.0, 0.0, 0.0), r_0 * 0.99, EDGE_PC))
    assert excluded == []


# ---------------------------------------------------------------------------
# Cross-check against independent brute force, for a P that is NOT the
# origin (the harder case per the design note) -- both the pruned
# algorithm and the brute-force scan build on the same
# `sector_position_pc`, but reach their answer via genuinely different
# traversal logic, so agreement confirms the pruning includes/excludes
# exactly the slots brute force would.
# ---------------------------------------------------------------------------

def _brute_force_within_radius(center, radius_pc, edge_pc, max_shell):
    """Iterates every slot of every shell 0..max_shell (inclusive) and
    filters by exact distance -- the "obviously correct, too slow to use
    for real" reference implementation `enumerate_sectors_within_radius`
    is checked against."""
    results = []
    for shell_index in range(max_shell + 1):
        n_k = shell_sector_count(shell_index)
        for slot_index in range(n_k):
            position = sector_position_pc(shell_index, slot_index, edge_pc)
            dist = math.dist(position, center)
            if dist <= radius_pc:
                results.append((shell_index, slot_index))
    return set(results)


@pytest.mark.parametrize("center,radius_pc", [
    ((5.0, 0.0, 0.0), 4.0),
    ((0.0, 6.0, 2.0), 5.0),
    ((-3.0, 3.0, 3.0), 6.0),
    ((10.0, -10.0, 5.0), 3.0),
])
def test_enumerate_matches_brute_force_over_several_shells(center, radius_pc):
    max_shell = 6  # generous -- shells 0-6 span radius up to ~22.9 pc at EDGE_PC
    expected = _brute_force_within_radius(center, radius_pc, EDGE_PC, max_shell)

    actual = {
        (shell, slot)
        for shell, slot, *_ in enumerate_sectors_within_radius(center, radius_pc, EDGE_PC)
    }

    assert actual == expected
    assert expected, "test parameters must produce at least one match to be a meaningful check"


def test_enumerate_distances_are_exact():
    center = (2.0, -1.0, 4.0)
    radius_pc = 5.0
    for shell, slot, x, y, z, dist in enumerate_sectors_within_radius(center, radius_pc, EDGE_PC):
        assert dist == pytest.approx(math.dist((x, y, z), center))
        assert dist <= radius_pc
        # And the returned position matches recomputing it directly.
        assert (x, y, z) == pytest.approx(sector_position_pc(shell, slot, EDGE_PC))


def test_enumerate_rejects_invalid_radius_or_edge():
    with pytest.raises(ValueError):
        list(enumerate_sectors_within_radius((0, 0, 0), -1.0, EDGE_PC))
    with pytest.raises(ValueError):
        list(enumerate_sectors_within_radius((0, 0, 0), 1.0, 0.0))


def test_enumerate_zero_radius_returns_at_most_the_exact_point():
    # A radius of exactly 0 should never match a slot unless P sits
    # exactly on that slot's own computed position (astronomically never
    # true for arbitrary P, but must not crash or over-match).
    results = list(enumerate_sectors_within_radius((123.456, -78.9, 0.1), 0.0, EDGE_PC))
    assert results == []


# ---------------------------------------------------------------------------
# Internal helper correctness (the two pruning steps in isolation).
# ---------------------------------------------------------------------------

def test_candidate_shell_range_is_exact_for_a_simple_case():
    # P at radius 10 pc, R = 2 pc, edge_pc = 1 pc -- shells are 1 pc thick,
    # nominal radius r_k = k + 0.5. Only k in [8, 11] can have
    # |r_k - 10| <= 2 (r_k in [8, 12]).
    k_min, k_max = _candidate_shell_range(p_norm=10.0, radius_pc=2.0, edge_pc=1.0)
    for k in range(k_min, k_max + 1):
        r_k = shell_radius_pc(k, 1.0)
        assert abs(r_k - 10.0) <= 2.0 + 1e-9
    # And the immediate neighbors just outside the range must fail the test.
    assert abs(shell_radius_pc(k_min - 1, 1.0) - 10.0) > 2.0
    assert abs(shell_radius_pc(k_max + 1, 1.0) - 10.0) > 2.0


def test_slot_index_bounds_round_trip_through_phi_for_index():
    n_k = shell_sector_count(20)
    for i in (0, 1, n_k // 2, n_k - 2, n_k - 1):
        phi = math.acos(1 - 2 * (i + 0.5) / n_k)
        i_min, i_max = slot_index_bounds_for_phi_range(phi, phi, n_k)
        assert i_min <= i <= i_max
