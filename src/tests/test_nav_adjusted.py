# tests/test_nav_adjusted.py
"""NAV.6: a course bent around the keep-out spheres it crosses (`nav_between(adjusted=True)`, `GET /api/nav?adjusted=1`)."""

import math

import pytest

from planetgen.db import obstacles, query, store
from planetgen.galaxy import unresolved_crossings
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.physics.units import ly_to_pc
from tests.publicids import pid
from tests.test_api import client  # noqa: F401  (the fixture)
from tests.test_corridor import _centers
from tests.test_route_edge_cases import _save_grid_sector


def _field_between(config, radius_ly):
    """An asteroid field on the line between the two sectors, 0.1 pc off it; returns its id."""
    centers = _centers(config)
    west, east = centers["West"], centers["East"]
    middle = [(a + b) / 2.0 for a, b in zip(west, east)]
    field = AsteroidField(SystemConfig())
    field.radius_ly = radius_ly
    placement = {"center_x_pc": middle[0], "center_y_pc": middle[1] + 0.1, "center_z_pc": middle[2],
                 "galactic_radius_pc": math.hypot(middle[0], middle[1])}
    conn = store.get_connection(config)
    try:
        field_id = store.insert_asteroid_field(conn, field, placement=placement)
        store.assign_uids(conn, phenomenon=("asteroid_fields", field_id))     # its name is its Object ID, which it takes here (GEN.1256)
        conn.commit()
    finally:
        conn.close()
    return field_id


def _systems(config):
    west = _save_grid_sector(config, "West", 0, [(0.0, 0.0, 0.0)])
    east = _save_grid_sector(config, "East", 4, [(0.0, 0.0, 0.0)])
    return west[0], east[0]


def _nav(config, west, east, **kwargs):
    conn = store.get_connection(config)
    try:
        return query.nav_between(conn, west, east, **kwargs)
    finally:
        conn.close()


def test_a_field_on_the_line_bends_the_direct_course_and_the_hop(mysql_config):
    west, east = _systems(mysql_config)
    field_id = _field_between(mysql_config, 20.0)
    result = _nav(mysql_config, west, east, adjusted=True)
    adjusted = result["adjusted"]
    assert adjusted["frame"] == "galactic" and adjusted["unit"] == "pc"
    assert adjusted["avoided"] == [f"asteroid_field:{field_id}"]
    assert [w["kind"] for w in adjusted["waypoints"]] == ["start", "tangent_in", "tangent_out", "end"]
    assert adjusted["length"] > adjusted["straight_length"]
    assert adjusted["length_ly"] == pytest.approx(adjusted["length"] / ly_to_pc(1.0))
    assert adjusted["extra_ly"] == pytest.approx(adjusted["length_ly"] - adjusted["straight_length_ly"])
    assert adjusted["warp_times"] and adjusted["fold_times"]
    assert adjusted["straight_length_ly"] == pytest.approx(result["direct"].distance_ly, rel=1e-6)
    route = result["route"]
    hops = route["hops"]
    assert all("adjusted" in hop for hop in hops)
    assert route["adjusted_distance_ly"] == pytest.approx(sum(hop["adjusted"]["length_ly"] for hop in hops))
    assert route["adjusted_distance_ly"] >= route["distance_ly"] - 1e-9
    assert len(route["adjusted_warp_times"]) == len(route["warp_times"])


def test_without_a_sphere_in_the_way_the_adjusted_course_is_the_straight_line(mysql_config):
    west, east = _systems(mysql_config)
    result = _nav(mysql_config, west, east, adjusted=True)
    adjusted = result["adjusted"]
    assert adjusted["avoided"] == [] and adjusted["extra"] == 0.0 and adjusted["warnings"] == []


def test_nothing_is_planned_unless_asked_for(mysql_config):
    west, east = _systems(mysql_config)
    _field_between(mysql_config, 20.0)
    result = _nav(mysql_config, west, east)
    assert "adjusted" not in result
    assert all("adjusted" not in hop for hop in result["route"]["hops"])
    assert "adjusted_distance_ly" not in result["route"]


def test_the_api_carries_the_adjusted_course_only_when_asked(mysql_config, client):  # noqa: F811
    west, east = _systems(mysql_config)
    field_id = _field_between(mysql_config, 20.0)
    base = f"/api/nav?from={pid('system', west)}&to={pid('system', east)}"
    plain = client.get(base).get_json()
    assert "adjusted" not in plain and "adjusted_distance_ly" not in plain["route"]
    shown = client.get(base + "&adjusted=1").get_json()
    assert shown["adjusted"]["avoided"] == [f"asteroid_field:{field_id}"]
    assert shown["route"]["adjusted_distance_ly"] >= shown["route"]["distance_ly"]
    assert len(shown["route"]["adjusted_warp_times"]) == len(shown["route"]["warp_times"])
    assert all("adjusted" in hop for hop in shown["route"]["hops"])


def _nebula_between(config, radius_ly, density):
    from planetgen.generation.phenomena.nebula import Nebula
    centers = _centers(config)
    west, east = centers["West"], centers["East"]
    middle = [(a + b) / 2.0 for a, b in zip(west, east)]
    nebula = Nebula(SystemConfig())
    nebula.radius_ly = radius_ly
    nebula.density_cm3 = density
    conn = store.get_connection(config)
    try:
        nebula_id = store.insert_nebula(conn, nebula, placement={
            "center_x_pc": middle[0], "center_y_pc": middle[1] + 0.05, "center_z_pc": middle[2],
            "galactic_radius_pc": math.hypot(middle[0], middle[1])})
        conn.commit()
    finally:
        conn.close()
    return nebula_id


def test_a_cloud_on_the_line_is_passed_through_with_its_length_and_column(mysql_config):  # NAV.28
    west, east = _systems(mysql_config)
    nebula_id = _nebula_between(mysql_config, 10.0, 100.0)
    result = _nav(mysql_config, west, east, adjusted=True)
    adjusted = result["adjusted"]
    assert adjusted["avoided"] == [] and adjusted["extra"] == 0.0           # a cloud is no obstacle
    [crossed] = adjusted["passed_through"]
    assert crossed["ref"] == f"nebula:{nebula_id}" and crossed["kind"] == "nebula"
    assert crossed["note"] and 0.0 < crossed["inside_length_ly"] <= 20.0 + 1e-6
    assert crossed["inside_length_ly"] == pytest.approx(20.0, rel=0.01)       # the line passes 0.05 pc from a 10 ly radius center
    assert crossed["column_cm2"] == pytest.approx(100.0 * crossed["inside_length"] * 3.0857e18, rel=1e-3)


def test_a_cloud_off_the_line_is_not_listed_and_the_names_of_bent_objects_come_with_the_course(mysql_config):
    west, east = _systems(mysql_config)
    field_id = _field_between(mysql_config, 20.0)
    adjusted = _nav(mysql_config, west, east, adjusted=True)["adjusted"]
    assert adjusted["passed_through"] == []
    assert list(adjusted["names"]) == [f"asteroid_field:{field_id}"] and adjusted["names"][f"asteroid_field:{field_id}"]


def test_a_field_between_two_systems_of_one_sector_bends_the_sector_scale_course_too(mysql_config):  # NAV.51
    first, second = _save_grid_sector(mysql_config, "Same", 0, [(-3.0, 0.0, 0.0), (3.0, 0.0, 0.0)])
    straight = _nav(mysql_config, first, second, adjusted=True)
    assert straight["scope"] == "sector" and straight["adjusted"]["avoided"] == []
    start, end = straight["adjusted"]["waypoints"][0]["p"], straight["adjusted"]["waypoints"][-1]["p"]
    middle = [(a + b) / 2.0 for a, b in zip(start, end)]
    field = AsteroidField(SystemConfig())
    field.radius_ly = 1.0
    placement = {"center_x_pc": middle[0], "center_y_pc": middle[1] + 0.01, "center_z_pc": middle[2],
                 "galactic_radius_pc": math.hypot(middle[0], middle[1])}
    conn = store.get_connection(mysql_config)
    try:
        field_id = store.insert_asteroid_field(conn, field, placement=placement)
        store.assign_uids(conn, phenomenon=("asteroid_fields", field_id))     # its name is its Object ID, which it takes here (GEN.1256)
        conn.commit()
    finally:
        conn.close()
    result = _nav(mysql_config, first, second, adjusted=True)
    assert result["scope"] == "sector"
    assert result["adjusted"]["avoided"] == [f"asteroid_field:{field_id}"]
    assert result["adjusted"]["length"] > result["adjusted"]["straight_length"]


# NAV.75: the unresolved stars a course passes through

def _core_field(ceiling=50.0):
    """An `UnresolvedField` of a toy galaxy whose core is capped at a low ceiling, so the cells are few and quick."""
    from tests.test_unresolved_crossings import SEED, SHAPE, EDGE
    skeleton = store.GalaxySkeletonInfo(
        shape=SHAPE, edge_pc=EDGE, outer_shell_index=200, expected_system_count_at_density_1=6.0, galaxy_seed=SEED,
        reach_shell_index=500, cell_star_ceiling=ceiling)
    return unresolved_crossings.UnresolvedField(skeleton)


def test_a_course_through_the_dense_core_says_how_many_unresolved_stars_it_crosses(mysql_config, monkeypatch):
    field = _core_field()
    monkeypatch.setattr(obstacles, "unresolved_field", lambda conn: field)
    conn = store.get_connection(mysql_config)
    try:
        course = obstacles.adjusted_course(conn, (-40.0, 3.0, 0.0), (40.0, 3.0, 0.0))
    finally:
        conn.close()
    shown = course["passed_through_unresolved"]
    assert set(shown) == {"crossings", "length_pc", "length_ly", "unresolved_systems", "cells"}
    assert shown["crossings"] > 0.0 and shown["cells"] > 5 and shown["length_pc"] == pytest.approx(80.0, rel=0.05)
    assert shown["length_ly"] == pytest.approx(shown["length_pc"] / ly_to_pc(1.0))
    assert shown["crossings"] == pytest.approx(field.crossings((-40.0, 3.0, 0.0), (40.0, 3.0, 0.0)).crossings, rel=0.02)


def test_a_course_outside_the_core_or_in_a_galaxy_without_a_shape_says_nothing_about_it(mysql_config, monkeypatch):
    field = _core_field()
    conn = store.get_connection(mysql_config)
    try:
        monkeypatch.setattr(obstacles, "unresolved_field", lambda c: field)
        far = field.reach_pc + 300.0
        assert "passed_through_unresolved" not in obstacles.adjusted_course(conn, (far, -50.0, 0.0), (far, 50.0, 0.0))
        monkeypatch.setattr(obstacles, "unresolved_field", lambda c: None)
        assert "passed_through_unresolved" not in obstacles.adjusted_course(conn, (-40.0, 3.0, 0.0), (40.0, 3.0, 0.0))
        monkeypatch.undo()
        assert obstacles.unresolved_field(conn) is None   # these tests save sectors but no galaxy shape
    finally:
        conn.close()


def test_the_count_rides_with_the_course_into_the_nav_result_and_the_api(mysql_config, client, monkeypatch):  # noqa: F811
    west, east = _systems(mysql_config)
    fake = {"crossings": 740.0, "length_pc": 200.0, "length_ly": 652.3, "unresolved_systems": 1.5e7, "cells": 51}
    monkeypatch.setattr(obstacles, "passed_through_unresolved", lambda conn, polyline: dict(fake))
    result = _nav(mysql_config, west, east, adjusted=True)
    assert result["adjusted"]["passed_through_unresolved"] == fake
    assert all(hop["adjusted"]["passed_through_unresolved"] == fake for hop in result["route"]["hops"])
    shown = client.get(f"/api/nav?from={pid('system', west)}&to={pid('system', east)}&adjusted=1").get_json()
    assert shown["adjusted"]["passed_through_unresolved"] == fake
    plain = client.get(f"/api/nav?from={pid('system', west)}&to={pid('system', east)}").get_json()
    assert "adjusted" not in plain
