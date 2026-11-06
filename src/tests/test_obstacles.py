"""NAV.25: the keep-out spheres a straight segment touches, in order along it."""

import pytest

from planetgen.db import obstacles, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.physics.units import ly_to_pc
from tests.test_corridor import _centers
from tests.test_route_edge_cases import _save_grid_sector


def _field(conn, x, y, radius_ly):
    field = AsteroidField(SystemConfig())
    field.radius_ly = radius_ly
    placement = {"center_x_pc": x, "center_y_pc": y, "center_z_pc": 0.0, "galactic_radius_pc": (x * x + y * y) ** 0.5}
    return store.insert_asteroid_field(conn, field, placement=placement)


def test_a_sphere_the_segment_touches_is_listed_and_one_it_misses_is_not(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        near = _field(conn, 50.0, 0.2, 1.0)        # 1 ly is 0.307 pc: the line passes through it
        far = _field(conn, 60.0, 5.0, 1.0)         # 5 pc off the line
        beyond = _field(conn, 150.0, 0.0, 1.0)     # on the line's extension, past the end
        conn.commit()
        found = obstacles.obstacles_along_segment(conn, (0.0, 0.0, 0.0), (100.0, 0.0, 0.0), kinds=["asteroid_field"])
    finally:
        conn.close()
    assert [row["id"] for row in found] == [near]
    row = found[0]
    assert row["kind"] == "asteroid_field" and row["basis"] == "radius" and row["note"] is None
    assert row["radius_pc"] == pytest.approx(ly_to_pc(1.0)) and row["distance_pc"] == pytest.approx(0.2)
    assert row["along_pc"] == pytest.approx(50.0) and row["center_pc"] == (50.0, 0.2, 0.0)
    assert not row["contains_start"] and not row["contains_end"]
    assert far not in [r["id"] for r in found] and beyond not in [r["id"] for r in found]


def test_a_sphere_that_holds_an_end_is_flagged_and_the_list_is_in_order_along_the_segment(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        late = _field(conn, 80.0, 0.0, 1.0)
        early = _field(conn, 20.0, 0.1, 1.0)
        around_start = _field(conn, 0.1, 0.0, 1.0)
        conn.commit()
        found = obstacles.obstacles_along_segment(conn, (0.0, 0.0, 0.0), (100.0, 0.0, 0.0), kinds=["asteroid_field"])
    finally:
        conn.close()
    assert [row["id"] for row in found] == [around_start, early, late]
    assert [row["contains_start"] for row in found] == [True, False, False]
    assert all(not row["contains_end"] for row in found)


def test_compact_objects_and_systems_are_obstacles_by_their_own_radius(mysql_config):
    west = _save_grid_sector(mysql_config, "West", 0, [(0.0, 0.0, 0.0)])
    east = _save_grid_sector(mysql_config, "East", 4, [(0.0, 0.0, 0.0)])
    centers = _centers(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        hole = BlackHole(SystemConfig())
        placement = {"center_x_pc": 8000.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 8000.0}
        hole_id = store.insert_black_hole(conn, hole, placement=placement)
        rogue_id = store.insert_rogue_planet(conn, RoguePlanet(SystemConfig()), placement=placement)
        conn.commit()
        systems = obstacles.obstacles_along_segment(conn, centers["West"], centers["East"], kinds=["system"])
        near_hole = obstacles.obstacles_along_segment(conn, (7000.0, 0.0, 0.0), (9000.0, 0.0, 0.0))
    finally:
        conn.close()
    assert west[0] in [row["id"] for row in systems] or east[0] in [row["id"] for row in systems]
    assert all(row["kind"] == "system" and row["basis"] == "perimeter" and row["radius_pc"] > 0 for row in systems)
    listed = {(row["kind"], row["id"]): row for row in near_hole}
    assert ("black_hole", hole_id) in listed and ("rogue_planet", rogue_id) in listed
    assert listed[("black_hole", hole_id)]["radius_km"] > 0 and listed[("rogue_planet", rogue_id)]["radius_km"] > 0
    assert listed[("black_hole", hole_id)]["basis"] in ("galactic_hill", "radius")


def test_a_kind_with_no_sphere_is_refused_and_an_empty_galaxy_has_no_obstacles(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        assert obstacles.obstacles_along_segment(conn, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)) == []
        with pytest.raises(ValueError):
            obstacles.obstacles_along_segment(conn, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), kinds=["nebula"])
    finally:
        conn.close()
