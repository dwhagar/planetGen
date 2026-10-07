# tests/test_navigation.py

"""
NAV feature tests: the pure math modules (`planetgen.galaxy.navigation`,
`planetgen.galaxy.nav_graph`) need no database and run unconditionally;
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
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.galaxy.nav_graph import build_knn_adjacency, shortest_path
from planetgen.galaxy.navigation import (
    FRAME_GALACTIC, FRAME_SECTOR, FRAME_SYSTEM, course_between, fold_speed_c, format_course, fold_travel_times, warp_speed_c, warp_travel_times,
)
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.generation.system import StarSystem


# ---------------------------------------------------------------------
# planetgen.galaxy.navigation -- pure math, no database.
# ---------------------------------------------------------------------

def test_course_toward_the_center_is_0_mark_0():
    # Ship out on +X, target at the galactic core: dead ahead is North.
    course = course_between((1000, 0, 0), (0, 0, 0))
    assert course.distance_ly == pytest.approx(1000.0)
    assert course.bearing_deg == pytest.approx(0.0)
    assert course.mark_deg == pytest.approx(0.0)
    assert course.frame == FRAME_GALACTIC
    assert format_course(course.bearing_deg, course.mark_deg) == "000 mark 000"


def test_course_bearing_is_measured_from_the_center_not_plus_x():
    # The same displacement (+Y) reads differently depending on where the
    # ship sits relative to the center: North is always toward it.
    from_plus_x = course_between((10, 0, 0), (10, 5, 0))
    from_plus_y = course_between((0, 10, 0), (0, 15, 0))
    assert from_plus_y.bearing_deg == pytest.approx(180.0)  # straight away from center
    # East = North x Up: with North = -X and Up = +Z, East is +Y.
    assert from_plus_x.bearing_deg == pytest.approx(90.0)
    assert course_between((10, 0, 0), (10, -5, 0)).bearing_deg == pytest.approx(270.0)


def test_course_bearing_ignores_height_above_the_plane():
    # North is the ship-to-center vector flattened onto the plane, so a
    # ship far above the plane still has North pointing at the axis.
    course = course_between((10, 0, 50), (0, 0, 50))
    assert course.bearing_deg == pytest.approx(0.0)
    assert course.mark_deg == pytest.approx(0.0)


def test_course_mark_is_elevation_mod_360():
    up = course_between((10, 0, 0), (10, 0, 10))
    assert up.elevation_deg == pytest.approx(90.0)
    assert up.mark_deg == pytest.approx(90.0)

    down = course_between((10, 0, 0), (10, 0, -10))
    assert down.elevation_deg == pytest.approx(-90.0)
    assert down.mark_deg == pytest.approx(270.0)

    shallow = course_between((10, 0, 0), (0, 0, -10))
    assert shallow.elevation_deg == pytest.approx(-45.0)
    assert shallow.mark_deg == pytest.approx(315.0)
    assert format_course(shallow.bearing_deg, shallow.mark_deg) == "000 mark 315"


def test_course_frame_center_and_up_vector():
    # A sector frame centered elsewhere, and a system frame whose ecliptic
    # is tilted so "up" is +X: North still points at the center.
    sector = course_between((0, 0, 0), (5, 0, 0), frame=FRAME_SECTOR, center=(5, 0, 0))
    assert (sector.bearing_deg, sector.mark_deg, sector.frame) == (pytest.approx(0.0), pytest.approx(0.0), FRAME_SECTOR)

    tilted = course_between((0, 0, 10), (0, 0, 0), frame=FRAME_SYSTEM, center=(0, 0, 0), up=(2, 0, 0))
    assert tilted.bearing_deg == pytest.approx(0.0)
    assert tilted.mark_deg == pytest.approx(0.0)
    above = course_between((0, 0, 10), (5, 0, 10), frame=FRAME_SYSTEM, up=(1, 0, 0))
    assert above.elevation_deg == pytest.approx(90.0)


def test_course_over_the_pole_falls_back_to_plus_x():
    # Directly over the center, North is undefined: +X stands in.
    course = course_between((0, 0, 10), (5, 0, 10))
    assert course.bearing_deg == pytest.approx(0.0)
    # Up along +X: +Y stands in instead.
    course = course_between((10, 0, 0), (10, 5, 0), up=(1, 0, 0))
    assert course.bearing_deg == pytest.approx(0.0)


def test_course_rejects_a_zero_up_vector():
    with pytest.raises(ValueError):
        course_between((1, 0, 0), (2, 0, 0), up=(0, 0, 0))


def test_course_between_coincident_points_is_defined():
    course = course_between((5, 5, 5), (5, 5, 5))
    assert course == (0.0, 0.0, 0.0, 0.0, FRAME_GALACTIC)


@pytest.mark.parametrize("bearing,mark,text", [
    (0.0, 0.0, "000 mark 000"),
    (45.4, 330.0, "045 mark 330"),
    (359.6, 359.7, "000 mark 000"),
    (7.0, 90.0, "007 mark 090"),
])
def test_format_course(bearing, mark, text):
    assert format_course(bearing, mark) == text


# Boss's tables (docs/design/navigation-frames.md, "Travel speeds"): factor,
# speed in c (1 decimal), and days per kpc (1 kpc = 3,261.56 ly, 1 ly per
# 365.25 days at 1c; the table rounded with a slightly longer kpc, so the
# day counts match to 1e-5).
WARP_TABLE = [
    (1, 1.0, 1_191_286),
    (2, 10.1, 118_191),
    (4, 101.6, 11_726),
    (8, 1024.0, 1_163),
    (9, 1520.1, 784),
    (9.5, 1936.0, 615),
    (9.9, 2822.7, 422),
    (9.995, 12201.9, 98),
]
FOLD_TABLE = [
    (4, 256.0, 4_653),
    (5, 750.0, 1_588),
    (6, 1944.0, 613),
    (6.5, 3060.1, 389),
    (7, 4802.0, 248),
    (7.5, 7593.8, 157),
    (8, 12288.0, 97),
    (8.5, 20880.2, 57),
]
_LY_PER_KPC = 3261.56
_DAYS_PER_YEAR = 365.25


@pytest.mark.parametrize("warp_factor,speed_c,days_per_kpc", WARP_TABLE)
def test_warp_speed_matches_boss_table(warp_factor, speed_c, days_per_kpc):
    speed = warp_speed_c(warp_factor)
    assert round(speed, 1) == speed_c
    assert _LY_PER_KPC / speed * _DAYS_PER_YEAR == pytest.approx(days_per_kpc, rel=1e-5, abs=0.5)


@pytest.mark.parametrize("fold_factor,speed_c,days_per_kpc", FOLD_TABLE)
def test_fold_speed_matches_boss_table(fold_factor, speed_c, days_per_kpc):
    speed = fold_speed_c(fold_factor)
    assert round(speed, 1) == speed_c
    assert _LY_PER_KPC / speed * _DAYS_PER_YEAR == pytest.approx(days_per_kpc, rel=1e-5, abs=0.5)


@pytest.mark.parametrize("speed_c", [warp_speed_c, fold_speed_c])
@pytest.mark.parametrize("factor", [0, -1, 10, 10.5])
def test_speed_curves_reject_factors_outside_0_to_10(speed_c, factor):
    with pytest.raises(ValueError):
        speed_c(factor)


def test_warp_curve_always_increases():
    factors = [0.5 + i * 0.01 for i in range(950)]
    speeds = [warp_speed_c(w) for w in factors]
    assert all(b > a for a, b in zip(speeds, speeds[1:]))


def test_warp_travel_times_default_factors_and_formula():
    legs = warp_travel_times(10.0)
    assert [leg.warp_factor for leg in legs] == [row[0] for row in WARP_TABLE]

    for leg in legs:
        assert leg.velocity_multiple_of_c == pytest.approx(warp_speed_c(leg.warp_factor))
        assert leg.years == pytest.approx(10.0 / leg.velocity_multiple_of_c)
        assert leg.formatted  # non-empty human-readable string

    # Warp 1 is c -- 10 ly takes 10 years.
    assert legs[0].years == pytest.approx(10.0)
    # Higher warp factors must always be faster (less travel time).
    assert all(a.years > b.years for a, b in zip(legs, legs[1:]))


def test_warp_travel_times_custom_factors():
    legs = warp_travel_times(1.0, warp_factors=(2, 4))
    assert [leg.warp_factor for leg in legs] == [2, 4]


def test_fold_travel_times_default_factors():
    legs = fold_travel_times(3261.56)
    assert [leg.fold_factor for leg in legs] == [row[0] for row in FOLD_TABLE]
    for leg, (_factor, _speed, days_per_kpc) in zip(legs, FOLD_TABLE):
        assert leg.years * _DAYS_PER_YEAR == pytest.approx(days_per_kpc, rel=1e-5, abs=0.5)
        assert leg.formatted


# ---------------------------------------------------------------------
# planetgen.galaxy.nav_graph -- pure math, no database.
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
        "center_x_pc": 6.132027875711011, "center_y_pc": 0.0, "center_z_pc": 0.0,  # ~20 ly
        "galactic_radius_pc": 6.132027875711011, "vertices_pc": empty_vertices,
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
    assert result["direct"].frame == FRAME_SECTOR
    assert result["route"]["path"][0] == ids["a"][0]
    assert result["route"]["path"][-1] == ids["a"][2]
    assert result["route"]["distance_ly"] == pytest.approx(4.0)
    assert [leg.warp_factor for leg in result["warp_times"]] == [row[0] for row in WARP_TABLE]
    assert [leg.fold_factor for leg in result["fold_times"]] == [row[0] for row in FOLD_TABLE]

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
    assert result["direct"].frame == FRAME_GALACTIC
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


# ---------------------------------------------------------------------
# queryDb.nav_between -- phenomenon endpoints (nebula, standing in for
# every _PHENOMENON_TYPE_TO_TABLE type -- they all share the same
# galaxy-frame-placement/no-sector-local-position shape this exercises).
# ---------------------------------------------------------------------

def _insert_nebula(config, center_x_pc=None, sector_id=None):
    """Inserts one nebula row, galaxy-placed at `(center_x_pc, 0, 0)` pc
    when given, else left unplaced (never generated into the galaxy) --
    mirrors `two_sector_galaxy`'s own hand-picked-coordinates convention
    rather than running the real placement algorithm, since only the
    resulting position/placement state matters for these tests."""
    nebula = Nebula(SystemConfig())
    placement = None
    if center_x_pc is not None:
        placement = {
            "center_x_pc": center_x_pc, "center_y_pc": 0.0, "center_z_pc": 0.0,
            "galactic_radius_pc": abs(center_x_pc),
        }
    conn = _db.get_connection(config)
    try:
        nebula_id = _db.insert_nebula(conn, nebula, sector_id=sector_id, placement=placement)
        conn.commit()
    finally:
        conn.close()
    return nebula_id


def test_nav_between_system_to_placed_phenomenon_is_galaxy_scope(two_sector_galaxy):
    config, ids = two_sector_galaxy
    # 10 pc =~ 32.6 ly from the galaxy origin, same axis sector A's system
    # 0 already sits at (galaxy x=0).
    nebula_id = _insert_nebula(config, center_x_pc=10.0)

    conn = _db.get_connection(config)
    try:
        result = nav_between(conn, ids["a"][0], nebula_id, to_kind="phenomenon", to_type="nebula")
    finally:
        conn.close()

    assert result["scope"] == "galaxy"
    assert result["destination_position"][0] > 0
    assert result["direct"].distance_ly == pytest.approx(result["destination_position"][0], rel=1e-6)
    assert result["route"] is not None
    assert result["route"]["path"][0] == ids["a"][0]
    assert result["route"]["path"][-1] == "phenomenon:nebula:" + str(nebula_id)
    assert result["route"]["positions"]["phenomenon:nebula:" + str(nebula_id)] == pytest.approx(
        result["destination_position"]
    )


def test_nav_between_phenomenon_to_phenomenon_is_galaxy_scope(two_sector_galaxy):
    config, ids = two_sector_galaxy
    nebula_a = _insert_nebula(config, center_x_pc=5.0)
    nebula_b = _insert_nebula(config, center_x_pc=-5.0)

    conn = _db.get_connection(config)
    try:
        result = nav_between(
            conn, nebula_a, nebula_b,
            from_kind="phenomenon", from_type="nebula", to_kind="phenomenon", to_type="nebula",
        )
    finally:
        conn.close()

    assert result["scope"] == "galaxy"
    assert result["direct"].distance_ly > 0
    # Both endpoints are real, distinct nodes -- a route must exist (the
    # kNN graph over every system plus these two one-off phenomenon nodes
    # is always at least this reachable in a fixture this small).
    assert result["route"] is not None


def test_nav_between_unplaced_phenomenon_is_unavailable(two_sector_galaxy):
    config, ids = two_sector_galaxy
    nebula_id = _insert_nebula(config, center_x_pc=None)

    conn = _db.get_connection(config)
    try:
        with pytest.raises(NavUnavailable):
            nav_between(conn, ids["a"][0], nebula_id, to_kind="phenomenon", to_type="nebula")
    finally:
        conn.close()


def test_nav_between_phenomenon_never_qualifies_for_sector_scope(two_sector_galaxy):
    # Even when a phenomenon's own "nearest sector" convenience link
    # (sector_id) matches a system's real sector, that's not real
    # containment -- NAV between them must still resolve at galaxy scope,
    # never sector scope (see _load_nav_phenomenon_endpoint's docstring).
    config, ids = two_sector_galaxy
    nebula_id = _insert_nebula(config, center_x_pc=1.0, sector_id=ids["sector_a"])

    conn = _db.get_connection(config)
    try:
        result = nav_between(conn, ids["a"][0], nebula_id, to_kind="phenomenon", to_type="nebula")
    finally:
        conn.close()

    assert result["scope"] == "galaxy"


def test_nav_between_raises_value_error_for_unrecognized_phenomenon_type(two_sector_galaxy):
    config, ids = two_sector_galaxy
    conn = _db.get_connection(config)
    try:
        with pytest.raises(ValueError):
            nav_between(conn, ids["a"][0], 1, to_kind="phenomenon", to_type="not_a_real_type")
    finally:
        conn.close()


def test_nav_between_reaches_a_placed_supernova_remnant(two_sector_galaxy):
    # supernova_remnants gained galaxy-frame placement columns in v28 (see
    # schema.sql's "v28" header note), so a placed one is a NAV endpoint
    # like any other phenomenon; an unplaced one is unavailable, not a
    # raw SQL error.
    config, ids = two_sector_galaxy
    remnant = SupernovaRemnant(SystemConfig())
    placement = {"center_x_pc": 10.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 10.0}
    conn = _db.get_connection(config)
    try:
        placed_id = _db.insert_supernova_remnant(conn, remnant, placement=placement)
        unplaced_id = _db.insert_supernova_remnant(conn, SupernovaRemnant(SystemConfig()))
        conn.commit()
        result = nav_between(conn, ids["a"][0], placed_id, to_kind="phenomenon", to_type="supernova_remnant")
        with pytest.raises(NavUnavailable):
            nav_between(conn, ids["a"][0], unplaced_id, to_kind="phenomenon", to_type="supernova_remnant")
    finally:
        conn.close()

    assert result["scope"] == "galaxy"
    assert result["route"]["path"][-1] == f"phenomenon:supernova_remnant:{placed_id}"


def test_nav_between_raises_value_error_for_missing_phenomenon(two_sector_galaxy):
    config, ids = two_sector_galaxy
    conn = _db.get_connection(config)
    try:
        with pytest.raises(ValueError):
            nav_between(conn, ids["a"][0], 999999999, to_kind="phenomenon", to_type="nebula")
    finally:
        conn.close()
