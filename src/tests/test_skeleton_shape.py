# tests/test_skeleton_shape.py

"""
Checks that the unfilled-sector skeleton -- the real, not-yet-generated
sector slots the density model says qualify (`predicted_star_count >= 1`,
the Galaxy Map's "planned" tier) -- traces a flattened disk, not a ball
around the galactic center.

Slots are sampled straight from the address scheme (random points in a
ball, snapped to their cell's center with `sector_address_at`) rather
than through
`galaxyViewport.planned_slots_in_view`, whose own view-radius cap
(`PLANNED_RADIUS_CAP_PC`) only ever returns a small sphere around the view
center; this test is about where qualifying slot centers really are.
Uses `generate.py plan`'s own default shape parameters.
"""

import math
import random

import pytest

from stellarObjects.galaxyDensity import build_galaxy_shape, predicted_star_count
from stellarObjects.galaxyGeometry import (
    cylindrical_radius_pc, layer_center_z_pc, ring_radius_pc, sector_address_at, sector_position_pc,
)
from stellarObjects.galaxySkeleton import expected_system_count_at_density_1

EDGE_LY = 11.5
EDGE_PC = EDGE_LY / 3.26156

# generate.py's own `plan` defaults (--disk-scale-length-pc etc.).
DEFAULT_SHAPE = build_galaxy_shape(
    disk_scale_length_pc=2800.0,
    disk_scale_height_pc=350.0,
    bulge_scale_radius_pc=200.0,
    bulge_amplitude=1.0,
    arm_count=2,
    pitch_angle_rad=math.radians(15.0),
    arm_amplitude=0.4,
)

SAMPLE_COUNT = 40000
MAX_SAMPLE_RADIUS_PC = 20000.0


def _sampled_slots():
    """(address, position) for SAMPLE_COUNT random slots: points drawn
    uniformly by radius (in a random direction) out past the galaxy's
    edge, each snapped to its own cell's center -- so nothing about the
    sampling itself favors the plane."""
    rng = random.Random(20260924)
    slots = []
    for _ in range(SAMPLE_COUNT):
        r = rng.uniform(0.0, MAX_SAMPLE_RADIUS_PC)
        cos_polar = rng.uniform(-1.0, 1.0)
        sin_polar = math.sqrt(1.0 - cos_polar * cos_polar)
        theta = rng.uniform(0.0, 2 * math.pi)
        point = (r * sin_polar * math.cos(theta), r * sin_polar * math.sin(theta), r * cos_polar)
        address = sector_address_at(point, EDGE_PC)
        slots.append((address, sector_position_pc(*address, EDGE_PC)))
    return slots


def _qualifying_positions(slots):
    expected = expected_system_count_at_density_1(EDGE_LY)
    return [position for _address, position in slots if predicted_star_count(position, DEFAULT_SHAPE, expected) >= 1.0]


def _approx(value):
    return pytest.approx(value, rel=1e-9, abs=1e-6)


def test_slot_centers_are_galaxy_frame_parsecs_on_their_own_ring_and_layer():
    for (ring_index, layer_index, _slot), position in _sampled_slots()[:2000]:
        assert cylindrical_radius_pc(position) == _approx(ring_radius_pc(ring_index, EDGE_PC))
        assert position[2] == _approx(layer_center_z_pc(layer_index, EDGE_PC))


def test_qualifying_slots_trace_a_flattened_disk_not_a_ball():
    positions = _qualifying_positions(_sampled_slots())
    assert len(positions) > 1000

    in_plane = sorted(math.hypot(x, y) for x, y, _z in positions)
    heights = sorted(abs(z) for _x, _y, z in positions)
    max_r = in_plane[-1]
    p95_height = heights[int(0.95 * len(heights))]

    # The disk reaches out to many scale lengths...
    assert max_r > 10000.0
    # ...but stays thin: a ball would have its heights comparable to its
    # in-plane extent.
    assert heights[-1] < 0.15 * max_r
    assert p95_height < 0.1 * max_r


def test_disk_thins_toward_its_outer_edge():
    positions = _qualifying_positions(_sampled_slots())

    def max_height(r_lo, r_hi):
        return max(abs(z) for x, y, z in positions if r_lo <= math.hypot(x, y) < r_hi)

    # Qualification needs a denser spot further out, so the qualifying
    # layer narrows with galactic radius instead of staying a fixed slab
    # (or bulging out like a sphere).
    assert max_height(8000.0, 12000.0) < max_height(0.0, 2000.0)
