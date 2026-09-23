# tests/test_galaxy_viewport.py

"""
Tests for `stellarObjects.galaxyViewport` -- the data layer behind the
interactive 3D Galaxy Map's live viewport queries (`queryDb.galaxy_view`).
See that module's own docstring for the three tiers (placed/planned/
density) this covers the pure, database-free "planned"/"density" halves
of.
"""

import math

from stellarObjects.galaxyDensity import build_galaxy_shape, predicted_star_count
from stellarObjects.galaxyGeometry import enumerate_sectors_within_radius
from stellarObjects.galaxySkeleton import expected_system_count_at_density_1
from stellarObjects.galaxyViewport import (
    DENSITY_SAMPLE_COUNT,
    PLANNED_RADIUS_CAP_PC,
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
        (shell_index, shell_slot_index)
        for shell_index, shell_slot_index, *_ in enumerate_sectors_within_radius(CENTER, RADIUS_PC, EDGE_PC)
    }
    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, EDGE_LY, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(),
    )
    got = {(entry["shell_index"], entry["shell_slot_index"]) for entry in results}
    assert got == expected
    assert all(entry["predicted_star_count"] is None and entry["relative_density"] is None for entry in results)


def test_filtered_with_shape_only_keeps_qualifying_slots():
    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, EDGE_LY, shape=SHAPE,
        expected_system_count_at_density_1=E, exclude_addresses=set(),
    )
    assert results  # this toy shape/radius combination has real qualifying content
    for entry in results:
        position = (entry["x"], entry["y"], entry["z"])
        assert predicted_star_count(position, SHAPE, E) >= 1.0
        assert entry["predicted_star_count"] >= 1.0


def test_excluded_addresses_are_skipped():
    baseline = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, EDGE_LY, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(),
    )
    assert len(baseline) >= 2
    excluded = {(baseline[0]["shell_index"], baseline[0]["shell_slot_index"])}

    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, EDGE_LY, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=excluded,
    )
    got = {(entry["shell_index"], entry["shell_slot_index"]) for entry in results}
    assert got.isdisjoint(excluded)
    assert len(results) == len(baseline) - 1


def test_results_are_sorted_closest_first():
    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, EDGE_LY, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(),
    )
    distances = [entry["distance_pc"] for entry in results]
    assert distances == sorted(distances)


def test_cap_truncates_to_the_closest_entries():
    full = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, EDGE_LY, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(), cap=1000,
    )
    assert len(full) > 3
    capped = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, EDGE_LY, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(), cap=3,
    )
    assert len(capped) == 3
    assert capped == full[:3]


def test_radius_is_clamped_to_the_planned_radius_cap():
    # A radius far beyond PLANNED_RADIUS_CAP_PC must not make this scan the
    # galaxy's outer shells (hundreds of millions of slots each) -- if it
    # weren't clamped, this call alone would never return.
    huge_radius = PLANNED_RADIUS_CAP_PC * 500
    clamped = planned_slots_in_view(
        CENTER, PLANNED_RADIUS_CAP_PC, EDGE_PC, EDGE_LY, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(), cap=100000,
    )
    unclamped_request = planned_slots_in_view(
        CENTER, huge_radius, EDGE_PC, EDGE_LY, shape=None,
        expected_system_count_at_density_1=None, exclude_addresses=set(), cap=100000,
    )
    got_clamped = {(e["shell_index"], e["shell_slot_index"]) for e in clamped}
    got_unclamped_request = {(e["shell_index"], e["shell_slot_index"]) for e in unclamped_request}
    assert got_clamped == got_unclamped_request


def test_designation_is_present_and_stable():
    results = planned_slots_in_view(
        CENTER, RADIUS_PC, EDGE_PC, EDGE_LY, shape=None,
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
