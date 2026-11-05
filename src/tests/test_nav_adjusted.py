# tests/test_nav_adjusted.py
"""NAV.6: a course bent around the keep-out spheres it crosses (`nav_between(adjusted=True)`, `GET /api/nav?adjusted=1`)."""

import math

import pytest

from planetgen.db import query, store
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
