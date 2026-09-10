# tests/test_navigation.py

"""
NAV feature tests: the pure math modules (`stellarObjects.navigation`,
`stellarObjects.navGraph`) need no database and run unconditionally;
`queryDb.nav_between` (the DB query/orchestration layer both `/api/nav`
and `html/nav.py` build on) is tested against a real, throwaway MySQL
database via `conftest.py`'s `mysql_config` fixture, seeded with real
`SpaceSector`/`save_sector` calls (same path `sectorGen.py`/`galaxyGen.py`
use) rather than hand-built rows, same convention as `test_api.py`.
"""

import math

import pytest

from queryDb import NavUnavailable, nav_between
from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.navGraph import build_knn_adjacency, shortest_path
from stellarObjects.navigation import course_between, warp_travel_times
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


# ---------------------------------------------------------------------
# stellarObjects.navigation -- pure math, no database.
# ---------------------------------------------------------------------

def test_course_between_planar_triangle():
    course = course_between((0, 0, 0), (3, 4, 0))
    assert course.distance_ly == pytest.approx(5.0)
    assert course.azimuth_deg == pytest.approx(math.degrees(math.atan2(4, 3)))
    assert course.altitude_deg == pytest.approx(0.0)


def test_course_between_azimuth_wraps_to_0_360():
    # Destination behind and to the left -> azimuth in the third quadrant,
    # not a negative angle.
    course = course_between((0, 0, 0), (-1, -1, 0))
    assert 180 < course.azimuth_deg < 270


def test_course_between_straight_up_and_down():
    up = course_between((0, 0, 0), (0, 0, 10))
    assert up.altitude_deg == pytest.approx(90.0)

    down = course_between((0, 0, 0), (0, 0, -10))
    assert down.altitude_deg == pytest.approx(-90.0)


def test_course_between_coincident_points_is_defined():
    course = course_between((5, 5, 5), (5, 5, 5))
    assert course == (0.0, 0.0, 0.0)


def test_warp_travel_times_default_factors_and_formula():
    legs = warp_travel_times(10.0)
    assert [leg.warp_factor for leg in legs] == [1, 3, 6, 9]

    for leg in legs:
        expected_velocity = leg.warp_factor ** (10 / 3)
        assert leg.velocity_multiple_of_c == pytest.approx(expected_velocity)
        assert leg.years == pytest.approx(10.0 / expected_velocity)
        assert leg.formatted  # non-empty human-readable string

    # Warp 1 is always exactly c -- 10 ly takes 10 years.
    assert legs[0].years == pytest.approx(10.0)
    # Higher warp factors must always be faster (less travel time).
    assert legs[0].years > legs[1].years > legs[2].years > legs[3].years


def test_warp_travel_times_custom_factors():
    legs = warp_travel_times(1.0, warp_factors=(2, 4))
    assert [leg.warp_factor for leg in legs] == [2, 4]


# ---------------------------------------------------------------------
# stellarObjects.navGraph -- pure math, no database.
# ---------------------------------------------------------------------

def test_build_knn_adjacency_is_symmetric():
    positions = {
        "A": (0, 0, 0),
        "B": (1, 0, 0),
        "C": (2, 0, 0),
        "D": (100, 100, 100),  # far outlier
    }
    graph = build_knn_adjacency(positions, k=1)

    for node, neighbors in graph.items():
        for neighbor, distance in neighbors.items():
            assert node in graph[neighbor], f"{node}->{neighbor} edge isn't symmetric"
            assert graph[neighbor][node] == pytest.approx(distance)


def test_build_knn_adjacency_includes_every_node_even_if_isolated():
    positions = {"only": (0, 0, 0)}
    graph = build_knn_adjacency(positions, k=5)
    assert graph == {"only": {}}


def test_shortest_path_prefers_the_true_shortest_route():
    # A straight line A-B-C-D; with k=1 each system only links to its
    # single nearest neighbor, so the path must hop through all of them.
    positions = {"A": (0, 0, 0), "B": (1, 0, 0), "C": (2, 0, 0), "D": (3, 0, 0)}
    graph = build_knn_adjacency(positions, k=1)

    path, distance = shortest_path(graph, "A", "D")
    assert path == ["A", "B", "C", "D"]
    assert distance == pytest.approx(3.0)


def test_shortest_path_same_start_and_end():
    graph = build_knn_adjacency({"A": (0, 0, 0), "B": (1, 0, 0)}, k=1)
    assert shortest_path(graph, "A", "A") == (["A"], 0.0)


def test_shortest_path_returns_none_when_unreachable():
    graph = {"A": {}, "B": {}}
    assert shortest_path(graph, "A", "B") is None


# ---------------------------------------------------------------------
# queryDb.nav_between -- real MySQL, seeded via SpaceSector/save_sector.
# ---------------------------------------------------------------------

def _make_system(star_type="G2V"):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.PLANETS = False
    return StarSystem(system_config=cfg), cfg


def _system_ids_by_name(conn, sector_id):
    rows = conn.execute(
        "SELECT id, name FROM star_systems WHERE sector_id = ? ORDER BY id", (sector_id,)
    ).fetchall()
    return [row["id"] for row in rows]


@pytest.fixture
def two_sector_galaxy(mysql_config):
    """
    Seeds two galaxy-placed sectors (20 ly apart on the galaxy-frame X
    axis) plus one standalone (non-galaxy-placed) sector, each with a
    small straight-line chain of systems -- enough to exercise same-
    sector routing, cross-sector routing, and every "unavailable" rule
    `nav_between` enforces.

    Returns:
        tuple: `(mysql_config, ids)` where `ids` is a dict with
              `sector_a`, `sector_b`, `sector_c` (the sectors' ids) and
              `a` (list of 3 system ids in sector A, in a straight line
              0/2/4 ly apart), `b` (2 system ids in sector B), `c` (1
              system id in the non-galaxy-placed sector), and
              `unplaced` (a system id with no sector at all).
    """
    sector_a = SpaceSector("Sector A", edge_ly=11.5)
    for i, position in enumerate([(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (4.0, 0.0, 0.0)]):
        system, cfg = _make_system()
        sector_a.add_system(system, position=position, system_config=cfg)

    sector_b = SpaceSector("Sector B", edge_ly=11.5)
    for position in [(-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)]:
        system, cfg = _make_system()
        sector_b.add_system(system, position=position, system_config=cfg)

    sector_c = SpaceSector("Sector C", edge_ly=11.5)
    system, cfg = _make_system()
    sector_c.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)

    empty_vertices = {"inner": [], "outer": []}
    sector_a_id = _db.save_sector(sector_a, config=mysql_config, galaxy_position={
        "center_x_pc": 0.0, "center_y_pc": 0.0, "center_z_pc": 0.0,
        "galactic_radius_pc": 0.0, "vertices_pc": empty_vertices,
    })
    sector_b_id = _db.save_sector(sector_b, config=mysql_config, galaxy_position={
        "center_x_pc": 6.1320300975969335, "center_y_pc": 0.0, "center_z_pc": 0.0,  # ~20 ly
        "galactic_radius_pc": 6.1320300975969335, "vertices_pc": empty_vertices,
    })
    sector_c_id = _db.save_sector(sector_c, config=mysql_config)  # no galaxy_position at all

    conn = _db.get_connection(mysql_config)
    try:
        ids = {
            "sector_a": sector_a_id,
            "sector_b": sector_b_id,
            "sector_c": sector_c_id,
            "a": _system_ids_by_name(conn, sector_a_id),
            "b": _system_ids_by_name(conn, sector_b_id),
            "c": _system_ids_by_name(conn, sector_c_id),
        }
        unplaced_cur = conn.execute("INSERT INTO system_configs (markdown) VALUES (0)")
        unplaced_config_id = unplaced_cur.lastrowid
        unplaced_cur = conn.execute(
            "INSERT INTO star_systems (sector_id, system_config_id, name, quadrant) VALUES (NULL, ?, ?, ?)",
            (unplaced_config_id, "Unplaced", "I"),
        )
        ids["unplaced"] = unplaced_cur.lastrowid
        conn.commit()
    finally:
        conn.close()

    return mysql_config, ids


def test_nav_between_same_sector(two_sector_galaxy):
    config, ids = two_sector_galaxy
    conn = _db.get_connection(config)
    try:
        result = nav_between(conn, ids["a"][0], ids["a"][2])
    finally:
        conn.close()

    assert result["scope"] == "sector"
    assert result["direct"].distance_ly == pytest.approx(4.0)
    assert result["route"]["path"][0] == ids["a"][0]
    assert result["route"]["path"][-1] == ids["a"][2]
    assert result["route"]["distance_ly"] == pytest.approx(4.0)
    assert [leg.warp_factor for leg in result["warp_times"]] == [1, 3, 6, 9]

    assert result["origin_position"] == pytest.approx((0.0, 0.0, 0.0))
    assert result["destination_position"] == pytest.approx((4.0, 0.0, 0.0))
    assert set(result["route"]["positions"]) == set(result["route"]["path"])
    assert result["route"]["positions"][ids["a"][0]] == pytest.approx(result["origin_position"])
    assert result["route"]["positions"][ids["a"][2]] == pytest.approx(result["destination_position"])


def test_nav_between_cross_sector_galaxy_scope(two_sector_galaxy):
    config, ids = two_sector_galaxy
    conn = _db.get_connection(config)
    try:
        result = nav_between(conn, ids["a"][0], ids["b"][0])
    finally:
        conn.close()

    assert result["scope"] == "galaxy"
    # Sector A's system 0 sits at galaxy x=0; sector B's system 0 sits at
    # galaxy x=~19 (20 ly center - 1 ly local offset).
    assert result["direct"].distance_ly == pytest.approx(19.0, abs=1e-6)
    assert result["route"] is not None
    assert result["route"]["path"][0] == ids["a"][0]
    assert result["route"]["path"][-1] == ids["b"][0]

    assert result["origin_position"] == pytest.approx((0.0, 0.0, 0.0))
    assert set(result["route"]["positions"]) == set(result["route"]["path"])


def test_nav_between_same_system_has_no_route(two_sector_galaxy):
    config, ids = two_sector_galaxy
    conn = _db.get_connection(config)
    try:
        result = nav_between(conn, ids["a"][0], ids["a"][0])
    finally:
        conn.close()

    assert result["direct"].distance_ly == pytest.approx(0.0)
    assert result["route"] is None


def test_nav_between_unavailable_when_either_system_unplaced(two_sector_galaxy):
    config, ids = two_sector_galaxy
    conn = _db.get_connection(config)
    try:
        with pytest.raises(NavUnavailable):
            nav_between(conn, ids["a"][0], ids["unplaced"])
    finally:
        conn.close()


def test_nav_between_unavailable_across_non_galaxy_sector(two_sector_galaxy):
    config, ids = two_sector_galaxy
    conn = _db.get_connection(config)
    try:
        with pytest.raises(NavUnavailable):
            nav_between(conn, ids["a"][0], ids["c"][0])
    finally:
        conn.close()


def test_nav_between_raises_value_error_for_missing_system(two_sector_galaxy):
    config, ids = two_sector_galaxy
    conn = _db.get_connection(config)
    try:
        with pytest.raises(ValueError):
            nav_between(conn, ids["a"][0], 999999999)
    finally:
        conn.close()
