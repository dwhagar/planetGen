# tests/test_galaxy_viewport.py

"""
Tests for `stellarObjects.galaxyViewport` -- the data layer behind the
interactive 3D Galaxy Map's live viewport queries (`queryDb.galaxy_view`).
See that module's own docstring for the three tiers (placed/planned/
density) this covers the pure, database-free "planned"/"density" halves
of.
"""

import math
import random

import pytest

from stellarObjects.galaxyDensity import build_galaxy_shape, predicted_star_count, relative_density
from stellarObjects.galaxyGeometry import enumerate_sectors_within_radius
from stellarObjects.galaxySkeleton import expected_system_count_at_density_1
from stellarObjects.galaxyViewport import (
    DENSITY_SAMPLE_COUNT,
    DENSITY_TILE_SAMPLE_COUNT,
    PLANNED_RADIUS_CAP_PC,
    PLANNED_TILE_MAX_EDGE_PC,
    TILE_MAX_LEVEL,
    TILE_ROOT_EDGE_PC,
    density_points_for_tile,
    parse_tile_key,
    planned_slots_in_tile,
    tile_bounds_pc,
    tile_edge_pc,
    tile_key,
    tile_keys_containing,
    tile_level_for_view_radius,
    tiles_intersecting_sphere,
    _sample_bulge_point_pc,
    _sample_disk_point_pc,
    density_sample_points,
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


def test_density_sample_points_empty_without_shape():
    assert density_sample_points(CENTER, RADIUS_PC, shape=None) == []


def test_density_sample_points_empty_for_nonpositive_radius():
    assert density_sample_points(CENTER, 0.0, shape=SHAPE) == []
    assert density_sample_points(CENTER, -5.0, shape=SHAPE) == []


def test_density_sample_points_count_and_bounds():
    points = density_sample_points(CENTER, RADIUS_PC, shape=SHAPE, count=250)
    assert len(points) == 250
    for point in points:
        dx = point["x"] - CENTER[0]
        dy = point["y"] - CENTER[1]
        dz = point["z"] - CENTER[2]
        distance = math.sqrt(dx * dx + dy * dy + dz * dz)
        assert distance <= RADIUS_PC + 1e-9
        assert point["relative_density"] >= 0.0


def test_density_sample_points_default_count():
    points = density_sample_points(CENTER, RADIUS_PC, shape=SHAPE)
    assert len(points) == DENSITY_SAMPLE_COUNT


def test_density_sample_points_deterministic_for_same_view():
    first = density_sample_points((10.0, -5.0, 2.0), 20.0, shape=SHAPE, count=50)
    second = density_sample_points((10.0, -5.0, 2.0), 20.0, shape=SHAPE, count=50)
    assert first == second


def test_density_sample_points_varies_with_view():
    at_origin = density_sample_points(CENTER, RADIUS_PC, shape=SHAPE, count=50)
    elsewhere = density_sample_points((500.0, 0.0, 0.0), RADIUS_PC, shape=SHAPE, count=50)
    assert at_origin != elsewhere


def test_density_sample_points_all_within_radius_even_at_wide_zoom():
    # A view radius far past this toy shape's own disk/bulge scale
    # (40pc/10pc) -- the mixture proposal is centered on the galactic
    # origin, not the view center, so this also confirms points still get
    # filtered down to the requested sphere rather than leaking outside it.
    wide_radius = 2000.0
    points = density_sample_points(CENTER, wide_radius, shape=SHAPE, count=200)
    assert len(points) == 200
    for point in points:
        distance = math.sqrt(point["x"] ** 2 + point["y"] ** 2 + point["z"] ** 2)
        assert distance <= wide_radius + 1e-6


def test_density_sample_points_biased_toward_higher_density_than_uniform():
    """
    The whole point of switching from a uniform-in-sphere draw to
    bulge/disk importance sampling: at a view radius much wider than the
    galaxy's own real scale, a uniform draw wastes nearly its entire
    budget in near-empty outskirts (most of a big sphere's *volume* is far
    from the core), while sampling from where the real mass concentrates
    should land a plotted cloud with a much higher average density instead
    -- otherwise the map reads as a sparse, shapeless scatter rather than
    a recognizable galaxy shape.
    """
    from stellarObjects.galaxyViewport import _sample_point_in_sphere

    wide_radius = 2000.0
    importance_points = density_sample_points(CENTER, wide_radius, shape=SHAPE, count=300)
    mean_importance_density = sum(p["relative_density"] for p in importance_points) / len(importance_points)

    rng = random.Random("uniform-baseline-for-comparison")
    uniform_densities = []
    for _ in range(300):
        x, y, z = _sample_point_in_sphere(rng, CENTER, wide_radius)
        uniform_densities.append(relative_density((x, y, z), SHAPE))
    mean_uniform_density = sum(uniform_densities) / len(uniform_densities)

    assert mean_importance_density > mean_uniform_density * 5


def test_sample_disk_point_pc_matches_the_disk_envelope_scale():
    rng = random.Random("disk-envelope-check")
    n = 4000
    points = [_sample_disk_point_pc(rng, SHAPE) for _ in range(n)]
    mean_r_cyl = sum(math.hypot(x, y) for x, y, _z in points) / n
    mean_abs_z = sum(abs(z) for _x, _y, z in points) / n
    # Exponential/Laplace means equal their own scale parameter exactly,
    # in the limit -- a generous tolerance covers finite-sample noise.
    assert mean_r_cyl == pytest.approx(SHAPE.disk_scale_length_pc, rel=0.15)
    assert mean_abs_z == pytest.approx(SHAPE.disk_scale_height_pc, rel=0.15)


def test_sample_bulge_point_pc_matches_the_bulge_envelope_scale():
    rng = random.Random("bulge-envelope-check")
    n = 4000
    points = [_sample_bulge_point_pc(rng, SHAPE) for _ in range(n)]
    mean_r = sum(math.sqrt(x * x + y * y + z * z) for x, y, z in points) / n
    assert mean_r == pytest.approx(SHAPE.bulge_scale_radius_pc, rel=0.15)



# --- Bounded planned search --------------------------------------------------

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


def test_planned_slots_skip_big_tiles_and_tiny_sectors():
    assert planned_slots_in_tile(TILE_MAX_LEVEL - 1, 2047, 2047, 2047, EDGE_PC, None, None, set()) == []
    # A 0.5 pc sector edge would put thousands of slots in a 16 pc tile.
    assert planned_slots_in_tile(TILE_MAX_LEVEL, 2048, 2048, 2048, 0.5, None, None, set()) == []


def test_density_points_for_tile_is_deterministic_and_sized():
    first = density_points_for_tile(8, 128, 128, 128, SHAPE)
    assert len(first) == DENSITY_TILE_SAMPLE_COUNT
    assert first == density_points_for_tile(8, 128, 128, 128, SHAPE)
    lo, hi = tile_bounds_pc(8, 128, 128, 128)
    center = tuple((lo[axis] + hi[axis]) / 2 for axis in range(3))
    reach = 2 * tile_edge_pc(8)
    assert all(math.dist((p["x"], p["y"], p["z"]), center) <= reach + 1e-6 for p in first)
    assert density_points_for_tile(8, 128, 128, 128, None) == []


def test_zoomed_in_density_cloud_follows_the_disk_not_a_ball():
    """
    A view a few hundred parsecs across, out in a Milky-Way-sized disk,
    is far too small for the galaxy-wide mixture proposal: it used to keep
    only a handful of candidates and top up the rest uniformly, drawing a
    ball of points around the view. Sampled locally against the real
    density instead, the cloud stays a slab that thins away from the
    plane.
    """
    milky_way = build_galaxy_shape(
        disk_scale_length_pc=2800.0,
        disk_scale_height_pc=350.0,
        bulge_scale_radius_pc=200.0,
        bulge_amplitude=1.0,
        arm_count=2,
        pitch_angle_rad=math.radians(15.0),
        arm_amplitude=0.4,
    )
    center = (8000.0, 0.0, 0.0)
    radius = 1000.0
    points = density_sample_points(center, radius, shape=milky_way, count=1600)
    assert len(points) == 1600
    assert all(math.dist((p["x"], p["y"], p["z"]), center) <= radius + 1e-6 for p in points)

    # A uniform ball puts ~21% of its volume more than 600 pc off the plane
    # (two caps of height 400 pc); the disk puts far less there.
    far_off_plane = sum(1 for p in points if abs(p["z"]) > 600.0) / len(points)
    assert far_off_plane < 0.1
    near_plane = sum(1 for p in points if abs(p["z"]) < 200.0) / len(points)
    assert near_plane > 0.4
