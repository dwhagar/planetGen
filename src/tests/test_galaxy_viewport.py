# tests/test_galaxy_viewport.py

"""
Tests for `planetgen.galaxy.viewport` -- the data layer behind the
interactive 3D Galaxy Map's live viewport queries (`queryDb.galaxy_view`).
See that module's own docstring for the tiers (placed/planned) this
covers the pure, database-free "planned" half of, plus the cube tiles.
"""

import itertools
import math
import random

import pytest

from planetgen.galaxy.cellkey import REACH_PC
from planetgen.galaxy.density import build_galaxy_shape, predicted_star_count
from planetgen.galaxy.spherical import enumerate_sectors_within_radius
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from planetgen.galaxy.viewport import (
    PLANNED_RADIUS_CAP_PC,
    PLANNED_TILE_MAX_EDGE_PC,
    TILE_MAX_LEVEL,
    TILE_ROOT_EDGE_PC,
    TILE_ROOT_EXTRA_LEVELS,
    parse_tile_key,
    planned_slots_in_tile,
    tile_bounds_pc,
    tile_edge_pc,
    tile_key,
    tile_keys_containing,
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
        (shell_index, band_index, wedge_index)
        for shell_index, band_index, wedge_index, *_ in enumerate_sectors_within_radius(CENTER, RADIUS_PC, EDGE_PC)
    }
    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(),
    )
    got = {(entry["shell_index"], entry["band_index"], entry["wedge_index"]) for entry in results}
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
    excluded = {(baseline[0]["shell_index"], baseline[0]["band_index"], baseline[0]["wedge_index"])}

    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=excluded,
    )
    got = {(entry["shell_index"], entry["band_index"], entry["wedge_index"]) for entry in results}
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
    got_clamped = {(e["shell_index"], e["band_index"], e["wedge_index"]) for e in clamped}
    got_unclamped_request = {(e["shell_index"], e["band_index"], e["wedge_index"]) for e in unclamped_request}
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


@pytest.mark.parametrize("level", [0, 3, 8, 12])
def test_a_ball_no_wider_than_a_tile_touches_at_most_27_tiles(level):
    radius = tile_edge_pc(level)
    keys = tiles_intersecting_sphere(level, (1234.5, -987.6, 12.3), radius)
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


def test_tiles_cover_the_reach():
    # MAP.239: the root is the grid's whole reach, so every address in it, to its very edge, is in a tile of every level.
    reach = REACH_PC
    assert TILE_ROOT_EDGE_PC == 2 * reach
    below = math.nextafter(reach, 0.0)
    for corner in itertools.product((-reach, below), repeat=3):
        keys = tile_keys_containing(corner)
        assert len(keys) == TILE_MAX_LEVEL + 1, corner
    assert tile_keys_containing((reach, 0.0, 0.0)) == []
    assert tile_keys_containing((0.0, -reach - 1.0, 0.0)) == []


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
            (entry["shell_index"], entry["band_index"], entry["wedge_index"])
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
    excluded = {(unfiltered[0]["shell_index"], unfiltered[0]["band_index"], unfiltered[0]["wedge_index"])}
    remaining = planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, None, None, excluded)
    assert len(remaining) == len(unfiltered) - 1
    filtered = planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, SHAPE, E, set())
    assert all(entry["predicted_star_count"] >= qualifying_threshold_star_count() for entry in filtered)


def test_planned_slots_in_tile_stay_inside_the_galaxy_outline():
    # MAP.118: a slot outside the stored outline (a shell past the last, or
    # a band beyond a shell's highest) is never offered as planned, since
    # generation refuses it.
    from planetgen.galaxy import spherical
    from planetgen.galaxy.skeleton import GalaxyBounds

    key = tiles_intersecting_sphere(TILE_MAX_LEVEL, (0.0, 0.0, 0.0), 0.0)[0]
    level, ix, iy, iz = parse_tile_key(key)
    unfiltered = planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, None, None, set())
    shells = sorted({entry["shell_index"] for entry in unfiltered})
    assert len(shells) > 2
    bounds = GalaxyBounds([(shell, min(1, spherical.band_range(shell)[1])) for shell in shells[:-1]], EDGE_PC)
    inside = planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, None, None, set(), bounds)
    assert inside and len(inside) < len(unfiltered)
    assert all(bounds.qualifies(e["shell_index"], e["band_index"]) for e in inside)
    assert all(e["shell_index"] != shells[-1] and abs(e["band_index"]) <= 1 for e in inside)
    # An empty outline (none built yet) filters nothing.
    assert planned_slots_in_tile(level, ix, iy, iz, EDGE_PC, None, None, set(), GalaxyBounds([], EDGE_PC)) \
        == unfiltered


def test_planned_slots_skip_big_tiles_and_tiny_sectors():
    half = 2 ** (TILE_MAX_LEVEL - 1)
    assert planned_slots_in_tile(TILE_MAX_LEVEL - 1, half - 1, half - 1, half - 1, EDGE_PC, None, None, set()) == []
    # A 0.5 pc sector edge would put thousands of slots in a 16 pc tile.
    assert planned_slots_in_tile(TILE_MAX_LEVEL, half, half, half, 0.5, None, None, set()) == []


# --- The distance cut (MAP.152) ---------------------------------------------

CUT_BOUNDS = ((-12600.0, -12600.0, -960.0), (12600.0, 12600.0, 960.0))


def _box_distance(key, point):
    lo, hi = tile_bounds_pc(*parse_tile_key(key))
    return math.sqrt(sum((point[a] - min(max(point[a], lo[a]), hi[a])) ** 2 for a in range(3)))


def test_the_distance_cut_covers_the_bounds_once_with_smaller_tiles_nearer():
    from planetgen.galaxy.viewport import DISTANCE_CUT_EDGES, DISTANCE_CUT_FINEST_LEVEL, distance_cut_tiles

    camera = (8000.0, 0.0, 0.0)
    keys = distance_cut_tiles(camera, camera, CUT_BOUNDS, half_angle_deg=180.0, max_tiles=10 ** 6)
    parsed = [parse_tile_key(key) for key in keys]
    # No tile lies inside another: the cut covers its space once.
    boxes = [tile_bounds_pc(*p) for p in parsed]
    for i, (lo_a, hi_a) in enumerate(boxes):
        for lo_b, hi_b in boxes[i + 1:]:
            assert not all(lo_a[a] < hi_b[a] and lo_b[a] < hi_a[a] for a in range(3))
    for key, (level, *_index) in zip(keys, parsed):
        edge = tile_edge_pc(level)
        if level < DISTANCE_CUT_FINEST_LEVEL:
            assert edge <= DISTANCE_CUT_EDGES * _box_distance(key, camera) + 1e-9, key
        else:
            assert level == DISTANCE_CUT_FINEST_LEVEL
    levels = [p[0] for p in parsed]
    assert max(levels) == DISTANCE_CUT_FINEST_LEVEL and min(levels) <= 3 + TILE_ROOT_EXTRA_LEVELS


def test_the_distance_cut_is_nearest_first_and_capped():
    from planetgen.galaxy.viewport import distance_cut_tiles

    camera = (8000.0, 0.0, 0.0)
    keys = distance_cut_tiles(camera, camera, CUT_BOUNDS, half_angle_deg=180.0, max_tiles=40)
    assert len(keys) == 40
    distances = [_box_distance(key, camera) for key in keys]
    assert distances == sorted(distances)


def test_the_cone_leaves_out_tiles_behind_the_camera():
    from planetgen.galaxy.viewport import distance_cut_tiles

    camera = (8000.0, 0.0, 0.0)
    ahead = distance_cut_tiles(camera, (8001.0, 0.0, 0.0), CUT_BOUNDS)
    everywhere = distance_cut_tiles(camera, (8001.0, 0.0, 0.0), CUT_BOUNDS, half_angle_deg=180.0)
    assert set(ahead) < set(everywhere)
    assert len(ahead) < 0.6 * len(everywhere)
    # Nothing behind: every kept tile touches the half space in front (its far corner is ahead).
    for key in ahead:
        assert tile_bounds_pc(*parse_tile_key(key))[1][0] > camera[0] - 1e-9


def test_the_fetch_set_holds_the_cut_and_the_cuts_either_side():
    from planetgen.galaxy.viewport import (
        DISTANCE_CUT_EDGES, DISTANCE_CUT_FADE_BAND, distance_cut_tiles, fetch_cut_tiles,
    )

    camera, look = (8012.0, 21.0, 9.0), (8000.0, 4.0, 4.0)
    fetch = fetch_cut_tiles(camera, look, CUT_BOUNDS)
    assert len(fetch) == len(set(fetch))
    for factor in (DISTANCE_CUT_FADE_BAND, 1.0 / DISTANCE_CUT_FADE_BAND):
        assert set(distance_cut_tiles(camera, look, CUT_BOUNDS, k=DISTANCE_CUT_EDGES * factor)) <= set(fetch)
    distances = [_box_distance(key, camera) for key in fetch]
    assert distances == sorted(distances)


def test_a_camera_far_out_needs_a_handful_of_coarse_tiles():
    from planetgen.galaxy.viewport import distance_cut_tiles

    keys = distance_cut_tiles((34000.0, -20000.0, 28000.0), (0.0, 0.0, 0.0), CUT_BOUNDS)
    assert 1 <= len(keys) <= 20
    assert max(parse_tile_key(key)[0] for key in keys) <= 4 + TILE_ROOT_EXTRA_LEVELS
