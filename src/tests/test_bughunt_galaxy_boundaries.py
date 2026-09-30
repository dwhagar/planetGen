# tests/test_bughunt_galaxy_boundaries.py

"""
Tier 1/Tier 2 bug-hunt coverage: boundary conditions in the cylindrical
sector grid (`galaxyGeometry.py`) -- ring 0 (whose cells are pie wedges
meeting on the galactic axis), the last slot in a ring (where the angle
wraps back to 0), cells far above the plane, and points exactly on a
cell boundary.

Tier 1 here is crash-freedom and basic finiteness/sanity; Tier 2
cross-checks every sector center against the ring radius and layer
height it claims, since a systematic drift there wouldn't crash anything
but would be a real placement bug.
"""

import math

import pytest

from stellarObjects.galaxyGeometry import (
    cylindrical_radius_pc, layer_center_z_pc, ring_bounds_pc, ring_radius_pc, ring_sector_count,
    sector_address_at, sector_cell_vertices_pc, sector_position_pc,
)
from tests.bughunt_support import tier2

EDGE_PC = 10.0


@pytest.mark.parametrize("ring_index", [0, 1, 2])
def test_every_slot_in_ring_produces_finite_position(ring_index):
    for slot_index in range(ring_sector_count(ring_index)):
        for layer_index in (-4095, 0, 4095):
            x, y, z = sector_position_pc(ring_index, layer_index, slot_index, EDGE_PC)
            assert all(math.isfinite(v) for v in (x, y, z))


@pytest.mark.parametrize("ring_index", [0, 1, 2, 50, 500])
def test_first_and_last_slot_cell_vertices_are_finite_and_nondegenerate(ring_index):
    n = ring_sector_count(ring_index)
    for slot_index in (0, n - 1):
        vertices = sector_cell_vertices_pc(ring_index, 0, slot_index, EDGE_PC)
        assert len(vertices) == 8
        for vertex in vertices:
            assert all(math.isfinite(v) for v in vertex)
        # The outer face's two angular corners are distinct points.
        assert math.dist(vertices[4], vertices[5]) > 0


def test_ring_0_inner_corners_sit_on_the_axis():
    for slot_index in range(ring_sector_count(0)):
        vertices = sector_cell_vertices_pc(0, 0, slot_index, EDGE_PC)
        for vertex in vertices[:4]:
            assert math.hypot(vertex[0], vertex[1]) == pytest.approx(0.0, abs=1e-12)


def test_sector_address_at_handles_exact_boundaries():
    # On the axis, on a ring boundary, on the +X axis (angle 0) and just
    # below 2*pi -- each must map to a valid address.
    assert sector_address_at((0.0, 0.0, 0.0), EDGE_PC) == (0, 0, 0)
    assert sector_address_at((EDGE_PC, 0.0, 0.0), EDGE_PC) == (1, 0, 0)
    ring, _layer, slot = sector_address_at((50.0, -1e-12, 0.0), EDGE_PC)
    assert 0 <= slot < ring_sector_count(ring)
    assert sector_address_at((0.0, 0.0, EDGE_PC / 2), EDGE_PC)[1] == 1
    assert sector_address_at((0.0, 0.0, -EDGE_PC / 2), EDGE_PC)[1] == 0


@tier2
def test_sector_centers_sit_on_their_ring_and_layer():
    for ring_index in (0, 1, 10, 200):
        nominal = ring_radius_pc(ring_index, EDGE_PC)
        lo, hi = ring_bounds_pc(ring_index, EDGE_PC)
        n = ring_sector_count(ring_index)
        for slot_index in range(0, n, max(1, n // 8)):
            for layer_index in (-7, 0, 3):
                position = sector_position_pc(ring_index, layer_index, slot_index, EDGE_PC)
                assert cylindrical_radius_pc(position) == pytest.approx(nominal, rel=1e-9)
                assert lo < cylindrical_radius_pc(position) < hi
                assert position[2] == pytest.approx(layer_center_z_pc(layer_index, EDGE_PC))
