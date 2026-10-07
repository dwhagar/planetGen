# tests/test_galaxy_viewport.py

"""
Tests for `planetgen.galaxy.viewport` -- the data layer behind the
interactive 3D Galaxy Map's live viewport queries (`queryDb.galaxy_view`).
See that module's own docstring for the tiers (placed/planned) this
covers the pure, database-free "planned" half of, plus the cube tiles.
"""

import math
import random
import sys

import pytest

from planetgen.galaxy.density import build_galaxy_shape, predicted_star_count
from planetgen.galaxy.geometry import enumerate_sectors_within_radius
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from planetgen.galaxy.viewport import (
    PLANNED_RADIUS_CAP_PC,
    PLANNED_TILE_MAX_EDGE_PC,
    TILE_MAX_LEVEL,
    TILE_ROOT_EDGE_PC,
    parse_tile_key,
    planned_slots_in_tile,
    tile_bounds_pc,
    tile_edge_pc,
    tile_key,
    tile_keys_containing,
    tile_level_for_view_radius,
    tiles_intersecting_sphere,
    planned_slots_in_view,
    qualifying_threshold_star_count,
)

EDGE_LY = 11.5
EDGE_PC = 3.526

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0,
    disk_scale_height_pc=12.0,
    bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0,
    arm_count=2,
    pitch_angle_rad=math.radians(15.0),
    arm_amplitude=0.4,
)
E = expected_system_count_at_density_1(edge_ly=EDGE_LY)

CENTER = (0.0, 0.0, 0.0)
RADIUS_PC = 15.0


def test_qualifying_threshold_is_one_star():
    assert qualifying_threshold_star_count() == 1.0


def test_unfiltered_without_shape_includes_every_enumerated_slot():
    expected = {
        (ring_index, layer_index, ring_slot_index)
        for ring_index, layer_index, ring_slot_index, *_ in enumerate_sectors_within_radius(CENTER, RADIUS_PC, EDGE_PC)
    }
    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(),
    )
    got = {(entry["ring_index"], entry["layer_index"], entry["ring_slot_index"]) for entry in results}
    assert got == expected
    assert all(entry["predicted_star_count"] is None and entry["relative_density"] is None for entry in results)


def test_filtered_with_shape_only_keeps_qualifying_slots():
    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, shape=SHAPE,
        expected_system_count_at_density_1=E, exclude_addresses=set(),
    )
    assert results  # this toy shape/radius combination has real qualifying content
    for entry in results:
        position = (entry["x"], entry["y"], entry["z"])
        assert predicted_star_count(position, SHAPE, E) >= 1.0
        assert entry["predicted_star_count"] >= 1.0


def test_excluded_addresses_are_skipped():
    baseline = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(),
    )
    assert len(baseline) >= 2
    excluded = {(baseline[0]["ring_index"], baseline[0]["layer_index"], baseline[0]["ring_slot_index"])}

    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=excluded,
    )
    got = {(entry["ring_index"], entry["layer_index"], entry["ring_slot_index"]) for entry in results}
    assert got.isdisjoint(excluded)
    assert len(results) == len(baseline) - 1


def test_results_are_sorted_closest_first():
    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(),
    )
    distances = [entry["distance_pc"] for entry in results]
    assert distances == sorted(distances)


def test_cap_truncates_to_the_closest_entries():
    full = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(), cap=1000,
    )
    assert len(full) > 3
    capped = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(), cap=3,
    )
    assert len(capped) == 3
    assert capped == full[:3]


def test_radius_is_clamped_to_the_planned_radius_cap():
    # A radius far beyond PLANNED_RADIUS_CAP_PC must not make this scan the
    # galaxy's outer rings (thousands of layers of slots each) -- if it
    # weren't clamped, this call alone would never return.
    huge_radius = PLANNED_RADIUS_CAP_PC * 500
    clamped = planned_slots_in_view(
        CENTER, PLANNED_RADIUS_CAP_PC, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(), cap=100000,
    )
    unclamped_request = planned_slots_in_view(
        CENTER, huge_radius, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(), cap=100000,
    )
    got_clamped = {(e["ring_index"], e["layer_index"], e["ring_slot_index"]) for e in clamped}
    got_unclamped_request = {(e["ring_index"], e["layer_index"], e["ring_slot_index"]) for e in unclamped_request}
    assert got_clamped == got_unclamped_request


def test_designation_is_present_and_stable():
    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(),
    )
    for entry in results:
        assert isinstance(entry["designation"], str)
        assert entry["designation"]


def test_planned_radius_cap_bounds_the_work_far_from_the_core():
    # A huge requested radius 8 kpc out -- the view that used to take
    # ~140 s and ~400 MB -- is clamped and returns promptly, capped.
    import time
    started = time.monotonic()
    result = planned_slots_in_view(
        (8000.0, 0.0, 0.0), 15000.0, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(),
    )
    assert len(result) == 4000
    assert max(entry["distance_pc"] for entry in result) <= PLANNED_RADIUS_CAP_PC
    assert time.monotonic() - started < 10


# --- Cube tiles --------------------------------------------------------------

def test_tile_key_round_trips_and_validates():
    assert parse_tile_key(tile_key(3, 1, 2, 7)) == (3, 1, 2, 7)
    for bad in ("", "1/2/3", "a/0/0/0", "-1/0/0/0", f"{TILE_MAX_LEVEL + 1}/0/0/0", "2/4/0/0", "2/0/-1/0"):
        with pytest.raises(ValueError):
            parse_tile_key(bad)


def test_tile_levels_halve_the_edge_and_tile_the_root_cube():
    assert tile_edge_pc(0) == TILE_ROOT_EDGE_PC
    assert tile_edge_pc(TILE_MAX_LEVEL) == PLANNED_TILE_MAX_EDGE_PC
    lo, hi = tile_bounds_pc(1, 1, 0, 1)
    assert lo == (0.0, -TILE_ROOT_EDGE_PC / 2, 0.0)
    assert hi == (TILE_ROOT_EDGE_PC / 2, 0.0, TILE_ROOT_EDGE_PC / 2)


@pytest.mark.parametrize("radius", [5.0, 16.0, 100.0, 777.0, 24000.0, 80000.0])
def test_view_level_tiles_are_at_least_the_view_radius_and_few(radius):
    level = tile_level_for_view_radius(radius)
    if level > 0:
        assert tile_edge_pc(level) >= radius or level == TILE_MAX_LEVEL
    assert level == TILE_MAX_LEVEL or tile_edge_pc(level + 1) < radius
    center = (1234.5, -987.6, 12.3)
    keys = tiles_intersecting_sphere(level, center, radius)
    assert 1 <= len(keys) <= 27


def test_tile_keys_containing_gives_one_holding_tile_per_level():
    point = (5.0, -1234.5, 16.0)
    keys = tile_keys_containing(point)
    assert len(keys) == TILE_MAX_LEVEL + 1
    for level, key in enumerate(keys):
        parsed = parse_tile_key(key)
        assert parsed[0] == level
        lo, hi = tile_bounds_pc(*parsed)
        assert all(lo[a] <= point[a] < hi[a] for a in range(3))
    # A point exactly on a boundary belongs to the tile above it, the same
    # half-open rule the sector-in-box query uses.
    assert tile_keys_containing((0.0, 0.0, 0.0))[1] == "1/1/1/1"
    assert tile_keys_containing((TILE_ROOT_EDGE_PC, 0.0, 0.0)) == []


def test_tiles_intersecting_sphere_covers_the_sphere():
    rng = random.Random(3)
    center = (-40.0, 25.0, 3.0)
    radius = 30.0
    keys = set(tiles_intersecting_sphere(TILE_MAX_LEVEL, center, radius))
    edge = tile_edge_pc(TILE_MAX_LEVEL)
    for _ in range(500):
        point = [center[axis] + rng.uniform(-radius, radius) for axis in range(3)]
        if math.dist(point, center) > radius:
            continue
        index = [math.floor((point[axis] + TILE_ROOT_EDGE_PC / 2) / edge) for axis in range(3)]
        assert tile_key(TILE_MAX_LEVEL, *index) in keys


def test_planned_slots_partition_into_tiles():
    # Every slot within a region shows up in exactly one small tile.
    center = (20.0, -10.0, 5.0)
    radius = 12.0
    by_tile = []
    for key in tiles_intersecting_sphere(TILE_MAX_LEVEL, center, radius):
        level, ix, iy, iz = parse_tile_key(key)
        by_tile.extend(
            (entry["ring_index"], entry["layer_index"], entry["ring_slot_index"])
            for entry in planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, None, None, set())
        )
    assert len(by_tile) == len(set(by_tile))
    expected = {(k, j, i) for k, j, i, *_ in enumerate_sectors_within_radius(center, radius, EDGE_PC)}
    assert expected <= set(by_tile)


def test_planned_slots_in_tile_respects_shape_and_exclusions():
    key = tiles_intersecting_sphere(TILE_MAX_LEVEL, (0.0, 0.0, 0.0), 0.0)[0]
    level, ix, iy, iz = parse_tile_key(key)
    unfiltered = planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, None, None, set())
    assert unfiltered and all(entry["predicted_star_count"] is None for entry in unfiltered)
    excluded = {(unfiltered[0]["ring_index"], unfiltered[0]["layer_index"], unfiltered[0]["ring_slot_index"])}
    remaining = planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, None, None, excluded)
    assert len(remaining) == len(unfiltered) - 1
    filtered = planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, SHAPE, E, set())
    assert all(entry["predicted_star_count"] >= qualifying_threshold_star_count() for entry in filtered)


def test_planned_slots_in_tile_stay_inside_the_galaxy_outline():
    # MAP.118: a slot outside the stored outline (a layer above the top,
    # or a ring past a layer's edge) is never offered as planned, since
    # generation refuses it.
    from planetgen.galaxy.skeleton import GalaxyBounds

    key = tiles_intersecting_sphere(TILE_MAX_LEVEL, (0.0, 0.0, 0.0), 0.0)[0]
    level, ix, iy, iz = parse_tile_key(key)
    unfiltered = planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, None, None, set())
    layers = sorted({entry["layer_index"] for entry in unfiltered})
    assert len(layers) > 1
    bounds = GalaxyBounds([(layer, 0) for layer in layers[:-1]], EDGE_PC)
    inside = planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, None, None, set(), bounds)
    assert inside
    assert all(bounds.contains(e["ring_index"], e["layer_index"]) for e in inside)
    assert all(e["layer_index"] != layers[-1] and e["ring_index"] == 0 for e in inside)
    # An empty outline (none built yet) filters nothing.
    assert planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, None, None, set(), GalaxyBounds([], EDGE_PC)) \
        == unfiltered


def test_planned_slots_skip_big_tiles_and_tiny_sectors():
    assert planned_slots_in_tile(TILE_MAX_LEVEL - 1, 2047, 2047, 2047, EDGE_PC, None, None, set()) == []
    # A 0.5 pc sector edge would put thousands of slots in a 16 pc tile.
    assert planned_slots_in_tile(TILE_MAX_LEVEL, 2048, 2048, 2048, 0.5, None, None, set()) == []


@pytest.mark.parametrize("radius", [5e-324, 2.2e-311, 1e-308, sys.float_info.min])
def test_tile_level_for_a_tiny_positive_radius_is_the_finest(radius):
    """MAP.90: a subnormal radius used to overflow log2; any tiny positive
    radius gets the finest level, here and in the page module's copy."""
    from planetgen.web.maps.galaxymap3d import _tile_level_for_view_radius

    assert tile_level_for_view_radius(radius) == TILE_MAX_LEVEL
    assert _tile_level_for_view_radius(radius) == TILE_MAX_LEVEL
