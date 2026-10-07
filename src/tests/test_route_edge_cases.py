# tests/test_route_edge_cases.py

"""
Route edge cases (TEST.79), written before NAV.12 from the hop-length
study (2026-10-02, `nav-hop-length/report.md` in the project's shared
files; design in docs/design/course-routing.md).

The cases: an isolated system; empty and unfilled sectors between the
endpoints; the galaxy edge and the halo (a lone system 2 kpc above the
disk); both endpoints in one sector when the best route leaves it; and a
route graph in separate pieces (the study's 714 islands from 2,000
sectors).

The graph-level cases run on `navGraph.build_route_graph` with the same
`k` and island links `queryDb.nav_between` uses and need no database;
the `nav_between` cases need MySQL (`conftest.py`'s `mysql_config`) and
are skipped without it. The island cases failed before NAV.34 (the plain
6-nearest graph had no route) and pass now that the islands are joined.
The cases marked `xfail(strict=True)` are NAV.12's: it turns them green
and removes the marks (the result keys they read, `longest_hop_ly` and a
per-hop `unknown_space` flag in `route["hops"]`, are this file's reading
of NAV.12's "the longest hop is shown" and "`/api/nav` returns the flag
per hop"; NAV.12 may rename them here).
"""

import math
import random

import pytest

from planetgen.db.query import nav_between
from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.geometry import sector_position_pc
from planetgen.galaxy.nav_graph import (
    build_knn_adjacency, build_route_graph, connected_components, join_islands, shortest_path,
)
from planetgen.tuning import (
    DEFAULT_SECTOR_EDGE_LY, NAV_ADJACENCY_K, NAV_ISLAND_LINKS,
)
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem
from stellarObjects.utils import ly_to_pc, pc_to_ly

EDGE_LY = DEFAULT_SECTOR_EDGE_LY
SUN_LY = (pc_to_ly(7896.0), 0.0, 0.0)  # the study's home: the Sun's radius, in the plane


# ---------------------------------------------------------------------
# Graph-level cases -- pure math, no database.
# ---------------------------------------------------------------------

def _sector_points(rng, center, count):
    """`count` systems scattered through one sector-sized cube."""
    return [tuple(c + rng.uniform(-EDGE_LY / 2, EDGE_LY / 2) for c in center) for _ in range(count)]


def _neighbourhood(rng, center, sectors_a_side=3, per_sector=7):
    """A generated area: a cube of `sectors_a_side`^3 sectors, each with
    `per_sector` systems -- well over the 7 that make an island."""
    points = []
    half = (sectors_a_side - 1) / 2
    for i in range(sectors_a_side):
        for j in range(sectors_a_side):
            for k in range(sectors_a_side):
                sector_center = (
                    center[0] + (i - half) * EDGE_LY,
                    center[1] + (j - half) * EDGE_LY,
                    center[2] + (k - half) * EDGE_LY,
                )
                points += _sector_points(rng, sector_center, per_sector)
    return points


def _positions(*groups):
    """`{id: point}` over every group, plus each group's ids in order."""
    positions, ids = {}, []
    for group in groups:
        group_ids = []
        for point in group:
            positions[len(positions)] = point
            group_ids.append(len(positions) - 1)
        ids.append(group_ids)
    return positions, ids


def _route(positions, start_id, end_id):
    """The route `nav_between` would find over these positions."""
    graph = build_route_graph(positions, NAV_ADJACENCY_K, NAV_ISLAND_LINKS)
    return shortest_path(graph, start_id, end_id)


def _hops(positions, path):
    return [math.dist(positions[a], positions[b]) for a, b in zip(path, path[1:])]


def _closest_pair_distance(positions, ids_a, ids_b):
    return min(math.dist(positions[a], positions[b]) for a in ids_a for b in ids_b)


def test_separate_generated_areas_were_islands_and_now_have_a_route():
    rng = random.Random(79)
    home = _neighbourhood(rng, SUN_LY)
    far = _neighbourhood(rng, (SUN_LY[0] + pc_to_ly(1000.0), 0.0, 0.0))  # 1 kpc out
    positions, (home_ids, far_ids) = _positions(home, far)

    # The bug NAV.34 fixes: the plain 6-nearest graph is two islands.
    plain = build_knn_adjacency(positions, NAV_ADJACENCY_K)
    assert len(connected_components(plain)) == 2
    assert shortest_path(plain, home_ids[0], far_ids[0]) is None

    found = _route(positions, home_ids[0], far_ids[0])
    assert found is not None
    path, distance = found
    # The gap is crossed once, by the closest pair the two areas have.
    gap = _closest_pair_distance(positions, home_ids, far_ids)
    hops = _hops(positions, path)
    assert max(hops) == pytest.approx(gap)
    assert sum(1 for hop in hops if hop > 100 * EDGE_LY) == 1
    assert distance == pytest.approx(sum(hops))


def test_many_islands_over_the_disk_are_all_joined():
    # The study's 714 islands, scaled down: generated sectors scattered
    # over the whole disk, most with 7 or more systems.
    rng = random.Random(714)
    groups = []
    for _ in range(300):
        radius = pc_to_ly(15000.0) * math.sqrt(rng.random())
        angle = rng.uniform(0.0, 2 * math.pi)
        center = (radius * math.cos(angle), radius * math.sin(angle), pc_to_ly(rng.gauss(0.0, 150.0)))
        groups.append(_sector_points(rng, center, rng.randint(1, 14)))
    positions, ids = _positions(*groups)

    plain = build_knn_adjacency(positions, NAV_ADJACENCY_K)
    assert len(connected_components(plain)) > 100

    graph = build_route_graph(positions, NAV_ADJACENCY_K, NAV_ISLAND_LINKS)
    assert len(connected_components(graph)) == 1

    west = min(positions, key=lambda system_id: positions[system_id][0])
    east = max(positions, key=lambda system_id: positions[system_id][0])
    path, distance = shortest_path(graph, west, east)
    direct = math.dist(positions[west], positions[east])
    # The study's cross-galaxy route was 1.29 times direct.
    assert direct < distance < 1.6 * direct
    assert path[0] == west and path[-1] == east


def test_an_isolated_system_is_reached():
    rng = random.Random(1)
    home = _neighbourhood(rng, SUN_LY)
    lone = [(SUN_LY[0], pc_to_ly(300.0), 0.0)]
    positions, (home_ids, (lone_id,)) = _positions(home, lone)

    # A lone system is no island: its 6 nearest link it in (today's long
    # hidden hop), so the last hop is one of those, not a joining link.
    path, _distance = _route(positions, home_ids[0], lone_id)
    assert path[-1] == lone_id
    last_hop = _hops(positions, path)[-1]
    assert _closest_pair_distance(positions, home_ids, [lone_id]) <= last_hop < pc_to_ly(300.0) + 3 * EDGE_LY


def test_a_lone_system_in_the_halo_is_reached():
    # The study's S5: one system 2 kpc above the disk.
    rng = random.Random(5)
    home = _neighbourhood(rng, SUN_LY)
    halo = [(SUN_LY[0], 0.0, pc_to_ly(2000.0))]
    positions, (home_ids, (halo_id,)) = _positions(home, halo)

    path, _distance = _route(positions, home_ids[0], halo_id)
    assert path[-1] == halo_id
    assert _hops(positions, path)[-1] > pc_to_ly(1900.0)


def test_a_system_past_the_galaxy_edge_is_reached():
    rng = random.Random(15)
    edge_area = _neighbourhood(rng, (pc_to_ly(15000.0), 0.0, 0.0))
    beyond = _sector_points(rng, (pc_to_ly(20000.0), 0.0, 0.0), 8)  # an island of its own
    positions, (edge_ids, beyond_ids) = _positions(edge_area, beyond)

    found = _route(positions, edge_ids[0], beyond_ids[0])
    assert found is not None
    assert found[0][-1] == beyond_ids[0]


def test_empty_sectors_between_the_endpoints_are_crossed_in_one_hop():
    # Two generated areas with a row of empty sectors between them: the
    # route crosses the empty run in a single hop, never stopping in it.
    rng = random.Random(3)
    west = _neighbourhood(rng, SUN_LY)
    east = _neighbourhood(rng, (SUN_LY[0] + 10 * EDGE_LY, 0.0, 0.0))  # 7 empty sectors between
    positions, (west_ids, east_ids) = _positions(west, east)

    west_end = min(west_ids, key=lambda system_id: positions[system_id][0])
    east_end = max(east_ids, key=lambda system_id: positions[system_id][0])
    path, _distance = _route(positions, west_end, east_end)
    stops_x = [positions[system_id][0] for system_id in path]
    empty_low, empty_high = SUN_LY[0] + 1.5 * EDGE_LY, SUN_LY[0] + 8.5 * EDGE_LY
    assert not any(empty_low < x < empty_high for x in stops_x)
    assert max(_hops(positions, path)) >= empty_high - empty_low


def test_join_islands_leaves_a_connected_graph_alone():
    positions = {i: (float(i), 0.0, 0.0) for i in range(10)}
    graph = build_knn_adjacency(positions, NAV_ADJACENCY_K)
    before = {system_id: dict(edges) for system_id, edges in graph.items()}
    assert join_islands(graph, positions, NAV_ISLAND_LINKS) == before


def test_one_island_link_still_joins_everything():
    # Each round links each island to its nearest only; rounds repeat
    # until one island is left (the study's 1-link pass left 208).
    rng = random.Random(208)
    groups = [_sector_points(rng, (rng.uniform(0, 5000), rng.uniform(0, 5000), 0.0), 8) for _ in range(40)]
    positions, _ids = _positions(*groups)
    graph = join_islands(build_knn_adjacency(positions, NAV_ADJACENCY_K), positions, 1)
    assert len(connected_components(graph)) == 1


def test_phenomenon_endpoint_key_joins_like_a_system():
    # nav_between adds a phenomenon endpoint under a string key beside
    # int system ids; joining must not compare the two kinds of key.
    rng = random.Random(9)
    positions, (home_ids,) = _positions(_neighbourhood(rng, SUN_LY))
    positions["phenomenon:nebula:1"] = (SUN_LY[0] + 500.0, 0.0, 0.0)
    path, _distance = _route(positions, home_ids[0], "phenomenon:nebula:1")
    assert path[-1] == "phenomenon:nebula:1"


# ---------------------------------------------------------------------
# queryDb.nav_between cases -- real MySQL, seeded via save_sector.
# ---------------------------------------------------------------------

EDGE_PC = ly_to_pc(EDGE_LY)
RING = 1980  # about the Sun's radius at the default sector edge


def _make_system():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    return StarSystem(system_config=cfg), cfg


def _save_grid_sector(config, name, slot, local_positions, layer=0):
    """Saves one sector at its grid address with systems at the given
    sector-local light-year positions; returns its system ids in order."""
    sector = SpaceSector(name, edge_ly=EDGE_LY)
    for position in local_positions:
        system, cfg = _make_system()
        sector.add_system(system, position=position, system_config=cfg)
    center = sector_position_pc(RING, layer, slot, EDGE_PC)
    sector_id = store.save_sector(sector, config=config, galaxy_position={
        "center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
        "galactic_radius_pc": math.hypot(center[0], center[1]),
        "vertices_pc": {"inner": [], "outer": []},
        "ring_index": RING, "layer_index": layer, "ring_slot_index": slot,
    })
    conn = store.get_connection(config)
    try:
        rows = conn.execute(
            "SELECT id FROM star_systems WHERE sector_id = ? ORDER BY id", (sector_id,)
        ).fetchall()
    finally:
        conn.close()
    return [row["id"] for row in rows]


def _cluster(center, count=8, spread=1.0):
    """`count` sector-local positions in a tight ball -- 7 or more make an
    island in the 6-nearest graph."""
    rng = random.Random(count * 1000 + int(center[0] * 10))
    return [tuple(c + rng.uniform(-spread, spread) for c in center) for _ in range(count)]


def _nav(config, from_id, to_id):
    conn = store.get_connection(config)
    try:
        return nav_between(conn, from_id, to_id)
    finally:
        conn.close()


def test_nav_between_routes_between_two_separately_generated_areas(mysql_config):
    # Two generated sectors 20 slots apart, each an island of its own.
    west = _save_grid_sector(mysql_config, "West", 0, _cluster((0.0, 0.0, 0.0)))
    east = _save_grid_sector(mysql_config, "East", 20, _cluster((0.0, 0.0, 0.0)))

    result = _nav(mysql_config, west[0], east[0])
    assert result["scope"] == "galaxy"
    assert result["route"] is not None
    assert result["route"]["path"][0] == west[0]
    assert result["route"]["path"][-1] == east[0]


def test_nav_between_reaches_a_lone_system_in_the_halo(mysql_config):
    home = _save_grid_sector(mysql_config, "Home", 0, _cluster((0.0, 0.0, 0.0)))
    # Layer 500 is about 2 kpc above the disk at the default edge.
    halo = _save_grid_sector(mysql_config, "Halo", 0, [(0.0, 0.0, 0.0)], layer=500)

    result = _nav(mysql_config, home[0], halo[0])
    assert result["route"] is not None
    assert result["route"]["path"][-1] == halo[0]


@pytest.mark.xfail(strict=True, reason="NAV.12: the longest hop is shown with the route")
def test_nav_between_reports_the_longest_hop(mysql_config):
    home = _save_grid_sector(mysql_config, "Home", 0, _cluster((0.0, 0.0, 0.0)))
    halo = _save_grid_sector(mysql_config, "Halo", 0, [(0.0, 0.0, 0.0)], layer=500)

    route = _nav(mysql_config, home[0], halo[0])["route"]
    hops = [
        math.dist(route["positions"][a], route["positions"][b])
        for a, b in zip(route["path"], route["path"][1:])
    ]
    assert route["longest_hop_ly"] == pytest.approx(max(hops))


@pytest.mark.xfail(strict=True, reason="NAV.12: a hop through unfilled sectors is flagged as unknown space")
def test_nav_between_flags_a_hop_through_unfilled_sectors(mysql_config):
    # Slots 1-4 between the two generated sectors are never generated.
    west = _save_grid_sector(mysql_config, "West", 0, _cluster((0.0, 0.0, 0.0)))
    east = _save_grid_sector(mysql_config, "East", 5, _cluster((0.0, 0.0, 0.0)))

    route = _nav(mysql_config, west[0], east[0])["route"]
    flags = [hop["unknown_space"] for hop in route["hops"]]
    assert len(flags) == len(route["path"]) - 1
    assert flags.count(True) == 1


@pytest.mark.xfail(strict=True, reason="NAV.12: same-sector routes may leave the sector")
def test_nav_between_same_sector_route_uses_nearer_stars_next_door(mysql_config):
    # Both ends in one sector, 10 ly apart near its top face, with
    # nothing else in it; the sector above holds a line of stars just
    # over that face, every one nearer to the ends than they are to each
    # other. Routing by nearest stars, the course runs along that line.
    near_top, near_bottom = EDGE_LY / 2 - 0.3, -EDGE_LY / 2 + 0.3
    home = _save_grid_sector(mysql_config, "Home", 0, [(-5.0, 0.0, near_top), (5.0, 0.0, near_top)])
    above = _save_grid_sector(
        mysql_config, "Above", 0, [(float(x), 0.0, near_bottom) for x in range(-5, 6)], layer=1,
    )

    result = _nav(mysql_config, home[0], home[1])
    assert result["route"] is not None
    assert set(result["route"]["path"]) & set(above)
