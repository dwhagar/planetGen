# tests/test_galaxy_geometry.py

"""
Tests for `stellarObjects.galaxyGeometry` -- the cylindrical
ring/layer/slot sector grid (`docs/design/galaxy-coordinate-system.md`,
"Cylindrical sector grid") and the radius-based neighborhood enumeration
built on it.

The enumeration tests cross-check the pruned algorithm against an
independent brute-force scan of every slot in the same rings and layers,
so they confirm the pruning includes/excludes exactly what brute force
would. The rest use small, hand-computable rings.
"""

import math
import random

import pytest

from stellarObjects.galaxyGeometry import (
    RING_SLOT_MULTIPLE,
    SectorCell,
    enumerate_sectors_within_radius,
    galactic_radius_pc,
    layer_bounds_pc,
    layer_center_z_pc,
    local_to_galaxy_pc,
    neighbor_addresses,
    parse_sector_designation,
    provisional_sector_designation,
    ring_bounds_pc,
    ring_radius_pc,
    ring_sector_count,
    sector_address_at,
    sector_cell_vertices_pc,
    sector_orientation,
    sector_position_pc,
    sector_quadrant,
    sector_zone,
    slot_angle_bounds,
)

EDGE_PC = 11.5 / 3.26156


@pytest.mark.parametrize("ring_index, expected", [
    (0, 4), (1, 8), (2, 16), (3, 20), (10, 64), (100, 632),
])
def test_ring_sector_count_rounds_centerline_circumference_to_a_multiple_of_4(ring_index, expected):
    assert ring_sector_count(ring_index) == expected


def test_ring_sector_count_is_always_a_positive_multiple_and_keeps_arcs_near_one_edge():
    for ring_index in range(0, 5000, 7):
        n = ring_sector_count(ring_index)
        assert n >= RING_SLOT_MULTIPLE and n % RING_SLOT_MULTIPLE == 0
        if ring_index >= 5:
            arc_edges = 2 * math.pi * (ring_index + 0.5) / n
            assert abs(arc_edges - 1.0) < 0.1


def test_ring_sector_count_rejects_negative_ring():
    with pytest.raises(ValueError):
        ring_sector_count(-1)


def test_ring_and_layer_bounds():
    assert ring_bounds_pc(3, 2.0) == (6.0, 8.0)
    assert ring_radius_pc(3, 2.0) == 7.0
    assert layer_bounds_pc(0, 2.0) == (-1.0, 1.0)
    assert layer_bounds_pc(-2, 2.0) == (-5.0, -3.0)
    assert layer_center_z_pc(-2, 2.0) == -4.0


def test_slot_angle_bounds_tile_the_full_circle():
    n = ring_sector_count(7)
    assert slot_angle_bounds(7, 0)[0] == 0.0
    assert slot_angle_bounds(7, n - 1)[1] == pytest.approx(2 * math.pi)
    for k in range(n - 1):
        assert slot_angle_bounds(7, k)[1] == pytest.approx(slot_angle_bounds(7, k + 1)[0])


def test_sector_position_pc_is_ring_centerline_slot_middle_layer_midplane():
    x, y, z = sector_position_pc(0, 2, 1, 2.0)
    # Ring 0 has 4 slots; slot 1's middle is at 135 degrees, radius 1.
    assert (x, y, z) == pytest.approx((-math.sqrt(0.5), math.sqrt(0.5), 4.0))


def test_sector_position_pc_rejects_out_of_range_slot():
    with pytest.raises(ValueError):
        sector_position_pc(0, 0, 4, EDGE_PC)
    with pytest.raises(ValueError):
        sector_position_pc(0, 0, -1, EDGE_PC)


def test_sector_address_at_round_trips_every_sector_center():
    for ring_index in range(0, 40):
        for slot_index in range(ring_sector_count(ring_index)):
            for layer_index in (-3, 0, 5):
                center = sector_position_pc(ring_index, layer_index, slot_index, EDGE_PC)
                assert sector_address_at(center, EDGE_PC) == (ring_index, layer_index, slot_index)


def test_sector_address_at_matches_bounds_for_random_points():
    rng = random.Random(7)
    for _ in range(500):
        point = (rng.uniform(-300, 300), rng.uniform(-300, 300), rng.uniform(-50, 50))
        ring, layer, slot = sector_address_at(point, EDGE_PC)
        r_lo, r_hi = ring_bounds_pc(ring, EDGE_PC)
        z_lo, z_hi = layer_bounds_pc(layer, EDGE_PC)
        t_lo, t_hi = slot_angle_bounds(ring, slot)
        theta = math.atan2(point[1], point[0]) % (2 * math.pi)
        assert r_lo <= math.hypot(point[0], point[1]) < r_hi
        assert z_lo <= point[2] < z_hi
        assert t_lo - 1e-12 <= theta < t_hi + 1e-12


def test_sector_orientation_is_radial_tangential_north():
    radial, tangential, north = sector_orientation((0.0, 5.0, 3.0))
    assert radial == pytest.approx((0.0, 1.0, 0.0))
    assert tangential == pytest.approx((-1.0, 0.0, 0.0))
    assert north == (0.0, 0.0, 1.0)
    # On the axis (only a hand-placed sector) it falls back to the
    # galaxy's own axes.
    assert sector_orientation((0.0, 0.0, 1.0)) == ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))


def test_local_to_galaxy_pc_applies_the_sector_axes():
    center = (0.0, 10.0, 2.0)
    assert local_to_galaxy_pc(center, (1.0, 2.0, 3.0)) == pytest.approx((-2.0, 11.0, 5.0))


def test_sector_cell_vertices_pc_sit_on_the_cell_bounds():
    vertices = sector_cell_vertices_pc(5, -1, 3, EDGE_PC)
    assert len(vertices) == 8
    r_bounds = ring_bounds_pc(5, EDGE_PC)
    z_bounds = layer_bounds_pc(-1, EDGE_PC)
    for index, (x, y, z) in enumerate(vertices):
        assert math.hypot(x, y) == pytest.approx(r_bounds[(index >> 2) & 1])
        assert z == pytest.approx(z_bounds[(index >> 1) & 1])
    with pytest.raises(ValueError):
        sector_cell_vertices_pc(0, 0, 9, EDGE_PC)


@pytest.mark.parametrize("ring_index", [0, 1, 12, 400])
def test_sector_cell_samples_land_inside_the_true_cell(ring_index):
    edge_ly = 11.5
    cell = SectorCell.for_ring(ring_index, edge_ly)
    rng = random.Random(ring_index)
    slot = ring_sector_count(ring_index) // 2
    center = sector_position_pc(ring_index, 2, slot, edge_ly)
    for _ in range(300):
        local = cell.sample(rng)
        assert cell.contains(local)
        point = local_to_galaxy_pc(center, local)
        assert sector_address_at(point, edge_ly) == (ring_index, 2, slot)


def test_sector_cell_volume_matches_an_edge_cubed_away_from_the_core():
    cell = SectorCell.for_ring(1000, 11.5)
    assert cell.volume == pytest.approx(11.5 ** 3, rel=0.01)
    # Every ring's cells together make the full annulus.
    for ring_index in (0, 3, 50):
        cell = SectorCell.for_ring(ring_index, 2.0)
        annulus = math.pi * (((ring_index + 1) * 2.0) ** 2 - (ring_index * 2.0) ** 2) * 2.0
        assert cell.volume * ring_sector_count(ring_index) == pytest.approx(annulus)


def test_sector_cell_contains_rejects_points_outside():
    cell = SectorCell.for_ring(10, 1.0)
    assert not cell.contains((0.0, 0.0, 0.6))
    assert not cell.contains((0.6, 0.0, 0.0))
    assert not cell.contains((0.0, 0.6, 0.0))
    assert cell.contains((0.0, 0.0, 0.0))


def _faces_touch(a, b):
    """Brute-force face adjacency between two cells: they share a face of
    nonzero area."""
    (ra, la, sa), (rb, lb, sb) = a, b
    if a == b:
        return False
    if ra == rb and la == lb:
        n = ring_sector_count(ra)
        return (sa - sb) % n in (1, n - 1)
    if ra == rb and sa == sb:
        return abs(la - lb) == 1
    if la == lb and abs(ra - rb) == 1:
        a0, a1 = slot_angle_bounds(ra, sa)
        b0, b1 = slot_angle_bounds(rb, sb)
        return min(a1, b1) - max(a0, b0) > 1e-12
    return False


@pytest.mark.parametrize("address", [(0, 0, 0), (0, 3, 2), (1, 0, 5), (2, -1, 15), (37, 4, 100), (38, 0, 0)])
def test_neighbor_addresses_match_brute_force_face_adjacency(address):
    ring, layer, _slot = address
    expected = set()
    for r in range(max(0, ring - 2), ring + 3):
        for l in range(layer - 2, layer + 3):
            for s in range(ring_sector_count(r)):
                if _faces_touch(address, (r, l, s)):
                    expected.add((r, l, s))
    got = neighbor_addresses(*address)
    assert len(got) == len(set(got))
    assert set(got) == expected


def test_sector_zone_groups_rings_into_about_100_ly_bands():
    assert sector_zone(0, 11.5) == 0
    assert sector_zone(8, 11.5) == 0
    assert sector_zone(9, 11.5) == 1
    assert sector_zone(5, 200.0) == 5


@pytest.mark.parametrize("x_pc, y_pc, expected", [(1, 1, 1), (-1, 1, 2), (-1, -1, 3), (1, -1, 4)])
def test_sector_quadrant_counts_counterclockwise_from_plus_x(x_pc, y_pc, expected):
    assert sector_quadrant(x_pc, y_pc) == expected


@pytest.mark.parametrize("address", [(0, 0, 0), (0, -4096, 3), (4073, 340, 25000), (12, -7, 50), (1, 4095, 7)])
def test_provisional_sector_designation_round_trips(address):
    code = provisional_sector_designation(*address)
    assert code == code.upper()
    assert parse_sector_designation(code) == address


def test_provisional_sector_designation_is_unique_over_a_block():
    seen = set()
    for ring in range(0, 6):
        for layer in range(-3, 4):
            for slot in range(ring_sector_count(ring)):
                seen.add(provisional_sector_designation(ring, layer, slot))
    assert len(seen) == sum(ring_sector_count(r) for r in range(6)) * 7


def test_provisional_sector_designation_rejects_out_of_range_layer():
    with pytest.raises(ValueError):
        provisional_sector_designation(0, 4096, 0)
    with pytest.raises(ValueError):
        parse_sector_designation(f"{(0 << 33) | (4096 << 20) | 9:X}")


def _brute_force_within_radius(center, radius_pc, edge_pc):
    r_c = math.hypot(center[0], center[1])
    max_ring = int((r_c + radius_pc) / edge_pc) + 1
    max_layer = int((abs(center[2]) + radius_pc) / edge_pc) + 1
    result = set()
    for ring in range(max_ring + 1):
        for layer in range(-max_layer, max_layer + 1):
            for slot in range(ring_sector_count(ring)):
                pos = sector_position_pc(ring, layer, slot, edge_pc)
                if math.dist(pos, center) <= radius_pc:
                    result.add((ring, layer, slot))
    return result


@pytest.mark.parametrize("center, radius_pc", [
    ((0.0, 0.0, 0.0), 20.0),
    ((30.0, -12.0, 4.0), 15.0),
    ((-70.0, 55.0, -20.0), 25.0),
    ((100.0, 0.0, 0.0), 3.0),
    ((5.0, 5.0, 0.0), 60.0),
])
def test_enumerate_matches_brute_force(center, radius_pc):
    got = {(r, l, s) for r, l, s, *_ in enumerate_sectors_within_radius(center, radius_pc, EDGE_PC)}
    assert got == _brute_force_within_radius(center, radius_pc, EDGE_PC)


def test_enumerate_randomized_against_brute_force():
    rng = random.Random(3)
    for _ in range(25):
        center = (rng.uniform(-120, 120), rng.uniform(-120, 120), rng.uniform(-30, 30))
        radius = rng.uniform(0, 30)
        got = {(r, l, s) for r, l, s, *_ in enumerate_sectors_within_radius(center, radius, EDGE_PC)}
        assert got == _brute_force_within_radius(center, radius, EDGE_PC)


def test_enumerate_reports_exact_centers_and_distances():
    center = (40.0, 10.0, 3.0)
    for ring, layer, slot, x, y, z, dist in enumerate_sectors_within_radius(center, 20.0, EDGE_PC):
        assert (x, y, z) == pytest.approx(sector_position_pc(ring, layer, slot, EDGE_PC))
        assert dist == pytest.approx(math.dist((x, y, z), center))


def test_enumerate_rejects_invalid_radius_or_edge():
    with pytest.raises(ValueError):
        list(enumerate_sectors_within_radius((0, 0, 0), -1.0, EDGE_PC))
    with pytest.raises(ValueError):
        list(enumerate_sectors_within_radius((0, 0, 0), 1.0, 0.0))


def test_enumerate_zero_radius_on_a_center_returns_just_that_sector():
    center = sector_position_pc(9, -2, 11, EDGE_PC)
    got = [(r, l, s) for r, l, s, *_ in enumerate_sectors_within_radius(center, 0.0, EDGE_PC)]
    assert got == [(9, -2, 11)]


def test_galactic_radius_pc_is_3d_distance():
    assert galactic_radius_pc((3.0, 4.0, 12.0)) == pytest.approx(13.0)
