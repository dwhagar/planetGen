"""NAV.10: routing through the systems near the direct line, A*, and the
line-segment object query."""

import math
import random

import pytest

from planetgen.db import corridor, query, store
from planetgen.db.query import nav_between
from planetgen.galaxy.nav_graph import build_knn_adjacency, shortest_path
from planetgen.physics.units import pc_to_ly
from tests.test_route_edge_cases import EDGE_PC, _cluster, _save_grid_sector


def test_segment_distance_is_to_the_nearest_point_of_the_segment():
    assert corridor.segment_distance((5, 3, 0), (0, 0, 0), (10, 0, 0)) == (3.0, 5.0)
    # Past an end, the end is the nearest point.
    assert corridor.segment_distance((14, 3, 0), (0, 0, 0), (10, 0, 0)) == (5.0, 10.0)
    assert corridor.segment_distance((3, 4, 0), (0, 0, 0), (0, 0, 0)) == (5.0, 0.0)


def test_a_star_finds_the_same_shortest_path_as_dijkstra():
    rng = random.Random(7)
    positions = {i: (rng.uniform(0, 100), rng.uniform(0, 100), rng.uniform(0, 100)) for i in range(150)}
    graph = build_knn_adjacency(positions, 4)
    for start, end in ((0, 149), (10, 77), (3, 90)):
        plain = shortest_path(graph, start, end)
        guided = shortest_path(graph, start, end, positions)
        assert (plain is None) == (guided is None)
        if plain is not None:
            assert guided[1] == pytest.approx(plain[1])


def test_a_phenomenon_key_can_be_the_goal():
    positions = {1: (0.0, 0.0, 0.0), 2: (1.1, 0.0, 0.0), "black-hole:7": (2.3, 0.2, 0.0)}
    graph = build_knn_adjacency(positions, 2)
    path, distance = shortest_path(graph, 1, "black-hole:7", positions)
    assert path == [1, 2, "black-hole:7"] or path == [1, "black-hole:7"]
    assert distance == pytest.approx(math.dist(positions[1], positions["black-hole:7"]), rel=0.2)


def _centers(config):
    conn = store.get_connection(config)
    try:
        rows = conn.execute("SELECT name, center_x_pc, center_y_pc, center_z_pc FROM sectors").fetchall()
    finally:
        conn.close()
    return {row["name"]: (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]) for row in rows}


def test_the_corridor_holds_the_systems_near_the_line_and_not_the_others(mysql_config):
    west = _save_grid_sector(mysql_config, "West", 0, _cluster((0.0, 0.0, 0.0)))
    middle = _save_grid_sector(mysql_config, "Middle", 3, _cluster((0.0, 0.0, 0.0)))
    east = _save_grid_sector(mysql_config, "East", 6, _cluster((0.0, 0.0, 0.0)))
    aside = _save_grid_sector(mysql_config, "Aside", 3, _cluster((0.0, 0.0, 0.0)), layer=40)
    centers = _centers(mysql_config)
    a, b = (tuple(pc_to_ly(v) for v in centers[name]) for name in ("West", "East"))
    conn = store.get_connection(mysql_config)
    try:
        loaded = corridor.positions_near_segment(conn, a, b, 40.0)
        elsewhere = corridor.positions_near_segment(conn, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1.0)
    finally:
        conn.close()
    assert set(west + middle + east) <= set(loaded)
    assert not set(aside) & set(loaded)
    assert elsewhere == {}


def test_objects_near_a_segment_come_in_order_along_it(mysql_config):
    west = _save_grid_sector(mysql_config, "West", 0, [(0.0, 0.0, 0.0)])
    east = _save_grid_sector(mysql_config, "East", 4, [(0.0, 0.0, 0.0)])
    centers = _centers(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        found = corridor.objects_near_segment(conn, centers["West"], centers["East"], 2.0, kinds=("system", "star"))
        none_far = corridor.objects_near_segment(conn, centers["West"], centers["East"], 2.0, kinds=("black_hole",))
    finally:
        conn.close()
    systems = [row for row in found if row["kind"] == "system"]
    assert [row["id"] for row in systems] == [west[0], east[0]]
    assert systems[0]["along_pc"] < systems[1]["along_pc"]
    assert any(row["kind"] == "star" and row["system_id"] == west[0] for row in found)
    assert none_far == []
    with pytest.raises(ValueError):
        corridor.objects_near_segment(None, centers["West"], centers["East"], -1.0)


def test_nav_between_reads_the_corridor_not_the_whole_galaxy(mysql_config, monkeypatch):
    west = _save_grid_sector(mysql_config, "West", 0, _cluster((0.0, 0.0, 0.0)))
    middle = _save_grid_sector(mysql_config, "Middle", 3, _cluster((0.0, 0.0, 0.0)))
    east = _save_grid_sector(mysql_config, "East", 6, _cluster((0.0, 0.0, 0.0)))
    aside = _save_grid_sector(mysql_config, "Aside", 3, _cluster((0.0, 0.0, 0.0)), layer=40)
    seen = []
    real = query.positions_near_segment
    monkeypatch.setattr(query, "positions_near_segment",
                        lambda conn, a, b, width: seen.append(real(conn, a, b, width)) or seen[-1])
    conn = store.get_connection(mysql_config)
    try:
        result = nav_between(conn, west[0], east[0])
    finally:
        conn.close()
    assert result["route"]["path"][0] == west[0] and result["route"]["path"][-1] == east[0]
    assert seen and not set(aside) & set(seen[0])
    assert set(result["route"]["path"]) <= set(west + middle + east)
