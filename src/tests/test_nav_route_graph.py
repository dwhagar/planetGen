# tests/test_nav_route_graph.py

"""
NAV.65: the graph round a route (`route.graph` of `db.query.nav_between` and `GET /api/nav?graph=1`).

The pure part (`nav_graph.relative_neighbours`, `nav_graph.route_graph`) is checked against a brute-force relative
neighbourhood graph and needs no database; the rest runs on a real MySQL database with systems placed in the galaxy
and compares the graph with a brute-force one over every system of the database.
"""

import itertools
import math
import random

import pytest

from planetgen.db import query, route_neighbours, store
from planetgen.galaxy.nav_graph import relative_neighbours, route_graph
from planetgen.physics.units import pc_to_ly
from planetgen.util import detmath
from tests.publicids import pid
from tests.test_api import client  # noqa: F401  (the fixture)
from planetgen.galaxy.sector import SpaceSector
from tests.test_route_edge_cases import EDGE_LY, _make_system, _save_grid_sector


def _cloud(rng, count, spread=10.0):
    return {index: tuple(rng.uniform(-spread, spread) for _ in range(3)) for index in range(count)}


def _brute_links(points):
    """Every relative-neighbourhood link of `points`, as a set of `frozenset({a, b})`, checked against every point."""
    links = set()
    for a, b in itertools.combinations(points, 2):
        between = detmath.dist(points[a], points[b])
        if not any(detmath.dist(points[a], points[c]) < between and detmath.dist(points[b], points[c]) < between
                   for c in points if c not in (a, b)):
            links.add(frozenset((a, b)))
    return links


# --- the rule --------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("seed", range(6))
def test_relative_neighbours_match_a_brute_force_graph(seed):
    points = _cloud(random.Random(seed), 40)
    expected = _brute_links(points)
    for origin, point in points.items():
        found = relative_neighbours(point, {key: other for key, other in points.items() if key != origin})
        assert {frozenset((origin, key)) for key, _distance in found} == {link for link in expected if origin in link}
        assert [distance for _key, distance in found] == sorted(distance for _key, distance in found)
        for key, distance in found:
            assert distance == pytest.approx(detmath.dist(point, points[key]))


def test_a_point_between_two_others_blocks_the_link_across_it():
    points = {"b": (1.0, 0.0, 0.0), "c": (2.0, 0.0, 0.0)}
    assert [key for key, _d in relative_neighbours((0.0, 0.0, 0.0), points)] == ["b"]
    # A point closer to one end only blocks nothing: both links stand.
    triangle = {"b": (1.0, 0.0, 0.0), "c": (0.4, 0.9, 0.0)}
    assert {key for key, _d in relative_neighbours((0.0, 0.0, 0.0), triangle)} == {"b", "c"}


def test_relative_neighbours_of_nothing_and_of_a_bad_point():
    assert relative_neighbours((0.0, 0.0, 0.0), {}) == []
    assert relative_neighbours((0.0, 0.0, 0.0), {1: (math.nan, 0.0, 0.0), 2: (3.0, 0.0, 0.0)}) == [(2, 3.0)]


# --- the graph round a route, no database ------------------------------------------------------------------------------

def _graph_of(points, path, cap=6, reach=math.inf):
    """`route_graph` over `points`, the systems round every stop being all the others."""
    def nearby(stop, point):
        return {key: other for key, other in points.items() if key != stop}, reach
    return route_graph(path, {stop: points[stop] for stop in path}, nearby, cap)


def _pairs(graph, **flags):
    return {frozenset((edge["from"], edge["to"])) for edge in graph["edges"]
            if all(edge[key] == value for key, value in flags.items())}


def test_the_hops_of_the_route_are_always_in_and_marked_on_route():
    points = {0: (0.0, 0.0, 0.0), 1: (1.0, 0.0, 0.0), 2: (2.0, 0.0, 0.0), 3: (3.0, 0.0, 0.0)}
    graph = _graph_of(points, [0, 2, 3])        # 0 to 2 is no relative-neighbourhood link (1 lies between)
    assert _pairs(graph, on_route=True) == {frozenset((0, 2)), frozenset((2, 3))}
    hop = next(edge for edge in graph["edges"] if edge["on_route"] and edge["from"] == 0)
    assert (hop["from"], hop["to"], hop["distance_ly"]) == (0, 2, 2.0)
    # System 1 is off the route and linked to the stops next to it; the hops come first.
    assert [node["id"] for node in graph["nodes"]] == [1]
    assert _pairs(graph, on_route=False) == {frozenset((0, 1)), frozenset((1, 2))}
    assert [edge["on_route"] for edge in graph["edges"]] == [True, True, False, False]


def test_a_link_between_two_stops_that_is_no_hop_is_listed_once_and_off_route():
    # A route 0 -> 2 -> 1 leaves out the short link 0-1.
    points = {0: (0.0, 0.0, 0.0), 1: (1.0, 0.0, 0.0), 2: (0.5, 40.0, 0.0)}
    graph = _graph_of(points, [0, 2, 1])
    links = [edge for edge in graph["edges"] if frozenset((edge["from"], edge["to"])) == frozenset((0, 1))]
    assert len(links) == 1 and links[0]["on_route"] is False
    assert (links[0]["from"], links[0]["to"]) == (0, 1)             # the earlier stop first
    assert graph["nodes"] == [] and graph["more"] == {}


def _icosahedron(radius, growth):
    """Twelve points round the origin, 63 degrees apart, at slightly growing radii: none shades another from the origin."""
    phi = (1 + math.sqrt(5)) / 2
    directions = [(0, s1, s2 * phi) for s1 in (1, -1) for s2 in (1, -1)]
    directions += [(s1, s2 * phi, 0) for s1 in (1, -1) for s2 in (1, -1)]
    directions += [(s1 * phi, 0, s2) for s1 in (1, -1) for s2 in (1, -1)]
    return {index: tuple(c / math.sqrt(1 + phi * phi) * (radius + growth * index) for c in direction)
            for index, direction in enumerate(directions)}


def test_a_stop_lists_at_most_the_cap_of_off_route_systems_and_counts_the_rest():
    points = {"stop": (0.0, 0.0, 0.0), "end": (500.0, 0.0, 0.0), **_icosahedron(5.0, 0.05)}
    whole = _graph_of(points, ["stop", "end"], cap=99)
    around = [edge for edge in whole["edges"] if not edge["on_route"] and edge["from"] == "stop"]
    assert len(around) == 12 and "stop" not in whole["more"]

    graph = _graph_of(points, ["stop", "end"], cap=4)
    off = [edge for edge in graph["edges"] if not edge["on_route"] and edge["from"] == "stop"]
    assert [edge["to"] for edge in off] == [edge["to"] for edge in around][:4]        # the nearest four
    assert graph["more"]["stop"] == 8
    assert [edge["distance_ly"] for edge in off] == sorted(edge["distance_ly"] for edge in off)
    assert not [edge for edge in _graph_of(points, ["stop", "end"], cap=0)["edges"] if edge["from"] == "stop" and not edge["on_route"]]
    assert _graph_of(points, ["stop", "end"], cap=0)["more"]["stop"] == 12


def test_a_system_next_to_two_stops_is_one_node_with_two_edges():
    points = {0: (0.0, 0.0, 0.0), 1: (4.0, 0.0, 0.0), 9: (2.0, 1.0, 0.0)}
    graph = _graph_of(points, [0, 1])
    assert [node["id"] for node in graph["nodes"]] == [9]
    assert _pairs(graph, on_route=False) == {frozenset((0, 9)), frozenset((1, 9))}


def test_stops_in_unknown_space_and_phenomenon_ends_take_part_like_any_stop():
    points = {"stone:black-hole:3": (10.0, 0.0, 0.0), "phenomenon:nebula:2": (20.0, 0.0, 0.0), 7: (0.0, 0.0, 0.0), 8: (10.0, 2.0, 0.0)}
    graph = _graph_of(points, [7, "stone:black-hole:3", "phenomenon:nebula:2"])
    assert [edge["on_route"] for edge in graph["edges"][:2]] == [True, True]
    assert frozenset((8, "stone:black-hole:3")) in _pairs(graph, on_route=False)
    assert [node["id"] for node in graph["nodes"]] == [8]


def test_the_same_route_gives_the_same_graph():
    points = _cloud(random.Random(3), 30)
    assert _graph_of(points, [0, 5, 9, 14]) == _graph_of(points, [0, 5, 9, 14])


# --- on a database -----------------------------------------------------------------------------------------------------

def _all_systems_ly(config):
    """`{id: (x, y, z)}` of every placed system, galaxy frame, light-years."""
    conn = store.get_connection(config)
    try:
        return {row["id"]: tuple(pc_to_ly(row[axis]) for axis in ("galaxy_x_pc", "galaxy_y_pc", "galaxy_z_pc"))
                for row in conn.execute("SELECT id, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc FROM star_systems").fetchall()}
    finally:
        conn.close()


def _nav(config, west, east, **kwargs):
    conn = store.get_connection(config)
    try:
        return query.nav_between(conn, west, east, **kwargs)
    finally:
        conn.close()


def _cloud_local(seed, count=16, spread=6.0):
    rng = random.Random(seed)
    return [tuple(rng.uniform(-spread, spread) for _ in range(3)) for _ in range(count)]


def _expected_links(config, route):
    """The links a stop of `route` has by the rule over every system of the database, and the hops."""
    everyone = _all_systems_ly(config)
    stops = [node for node in route["path"]]
    links = set()
    for stop in stops:
        for key, _distance in relative_neighbours(everyone[stop], {k: p for k, p in everyone.items() if k != stop}):
            links.add(frozenset((stop, key)))
    links |= {frozenset((a, b)) for a, b in zip(stops, stops[1:])}
    return links, everyone


@pytest.mark.parametrize("west_wedge, east_wedge", [(0, None), (0, 20), (2000, None), (2000, 2020)])
def test_the_graph_holds_the_hops_and_the_links_the_rule_gives(mysql_config, monkeypatch, west_wedge, east_wedge):
    """`None`: both ends in one sector (the route is turned back into the sector's frame, NAV.12); else two sectors apart.
    The wedge numbered 2000 lies at a longitude where the sector's axes are turned well away from the galaxy's (NAV.73)."""
    monkeypatch.setattr(query, "NAV_GRAPH_NEIGHBOURS", 99)
    west_ids = _save_grid_sector(mysql_config, "West", west_wedge, _cloud_local(1))
    if east_wedge is None:
        first, last = west_ids[0], west_ids[-1]
    else:
        east_ids = _save_grid_sector(mysql_config, "East", east_wedge, _cloud_local(2))
        first, last = west_ids[0], east_ids[-1]

    route = _nav(mysql_config, first, last, with_graph=True)["route"]
    graph = route["graph"]
    expected, everyone = _expected_links(mysql_config, route)

    assert {frozenset((edge["from"], edge["to"])) for edge in graph["edges"]} == expected
    hops = {frozenset(pair) for pair in zip(route["path"], route["path"][1:])}
    assert {frozenset((edge["from"], edge["to"])) for edge in graph["edges"] if edge["on_route"]} == hops
    assert len(graph["edges"]) == len({frozenset((edge["from"], edge["to"])) for edge in graph["edges"]})
    assert graph["more"] == {}
    hop_lengths = {frozenset((hop["from"], hop["to"])): hop["distance_ly"] for hop in route["hops"]}
    assert all(position == pytest.approx(everyone[node], abs=1e-6) for node, position in route["positions"].items()
               if east_wedge is not None)
    for edge in graph["edges"]:
        pair = frozenset((edge["from"], edge["to"]))
        if edge["on_route"]:
            assert edge["distance_ly"] == hop_lengths[pair]          # the route's own length, as `hops` lists it
            assert edge["distance_ly"] == pytest.approx(detmath.dist(everyone[edge["from"]], everyone[edge["to"]]), abs=1e-6)
        else:
            assert edge["distance_ly"] == pytest.approx(detmath.dist(everyone[edge["from"]], everyone[edge["to"]]), abs=1e-6)
            assert edge["from"] in route["path"] and edge["to"] not in route["path"] or (
                route["path"].index(edge["from"]) < route["path"].index(edge["to"]))
    # A node is no stop and sits where the index puts it, in the galaxy frame.
    for node in graph["nodes"]:
        assert node["id"] not in route["path"]
        assert node["position_ly"] == pytest.approx(everyone[node["id"]], abs=1e-6)
        assert any(edge["to"] == node["id"] for edge in graph["edges"])


def test_a_route_inside_a_sector_with_no_galaxy_placement_has_its_graph_from_the_sector(mysql_config, monkeypatch):
    monkeypatch.setattr(query, "NAV_GRAPH_NEIGHBOURS", 99)
    sector = SpaceSector("Loose", edge_ly=EDGE_LY)
    for position in _cloud_local(4):
        system, cfg = _make_system()
        sector.add_system(system, position=position, system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        positions = query._sector_local_positions(conn, sector_id)
        ids = sorted(positions)
        route = query.nav_between(conn, ids[0], ids[-1], with_graph=True)["route"]
    finally:
        conn.close()

    expected = {frozenset((a, b)) for a, b in zip(route["path"], route["path"][1:])}
    for stop in route["path"]:
        expected |= {frozenset((stop, key)) for key, _d in
                     relative_neighbours(positions[stop], {k: p for k, p in positions.items() if k != stop})}
    graph = route["graph"]
    assert {frozenset((edge["from"], edge["to"])) for edge in graph["edges"]} == expected
    for node in graph["nodes"]:
        assert node["position_ly"] == pytest.approx(positions[node["id"]], abs=1e-9)


def test_the_cap_cuts_the_systems_off_the_route_and_counts_what_it_cut(mysql_config, monkeypatch):
    ids = _save_grid_sector(mysql_config, "West", 0, _cloud_local(5, count=24, spread=8.0))
    monkeypatch.setattr(query, "NAV_GRAPH_NEIGHBOURS", 99)
    whole = _nav(mysql_config, ids[0], ids[-1], with_graph=True)["route"]["graph"]
    monkeypatch.setattr(query, "NAV_GRAPH_NEIGHBOURS", 1)
    cut = _nav(mysql_config, ids[0], ids[-1], with_graph=True)["route"]["graph"]

    off_by_stop = lambda graph: {stop: [e for e in graph["edges"] if not e["on_route"] and e["from"] == stop
                                        and e["to"] in {n["id"] for n in graph["nodes"]}]
                                 for stop in {e["from"] for e in graph["edges"]}}
    for stop, edges in off_by_stop(cut).items():
        assert len(edges) <= 1
        full = off_by_stop(whole).get(stop, [])
        assert len(edges) + cut["more"].get(stop, 0) == len(full)
        assert [e["to"] for e in edges] == [e["to"] for e in full][:len(edges)]       # the nearest are kept
    assert cut["more"] and all(count > 0 for count in cut["more"].values())
    assert {e["to"] for e in cut["edges"] if e["on_route"]} == {e["to"] for e in whole["edges"] if e["on_route"]}


def test_no_graph_is_read_unless_asked_for(mysql_config):
    ids = _save_grid_sector(mysql_config, "West", 0, _cloud_local(1))
    assert "graph" not in _nav(mysql_config, ids[0], ids[-1])["route"]


def test_the_nearest_systems_of_a_point_are_the_brute_force_nearest(mysql_config):
    ids = _save_grid_sector(mysql_config, "West", 0, _cloud_local(7, count=30, spread=9.0))
    everyone = _all_systems_ly(mysql_config)
    point = (everyone[ids[3]][0] + 0.7, everyone[ids[3]][1] - 0.4, everyone[ids[3]][2] + 0.2)
    conn = store.get_connection(mysql_config)
    try:
        found, reach = route_neighbours.nearby_systems(conn, point, 8, exclude_id=ids[0])
        sparse, sparse_reach = route_neighbours.nearby_systems(conn, point, 500)
    finally:
        conn.close()
    brute = sorted((detmath.dist(point, other), key) for key, other in everyone.items() if key != ids[0])[:8]
    assert list(found) == [key for _d, key in brute]
    assert reach == pytest.approx(brute[-1][0], rel=1e-6)
    assert found[ids[3]] == pytest.approx(everyone[ids[3]], abs=1e-6)
    assert set(sparse) == set(everyone) and sparse_reach >= max(detmath.dist(point, p) for p in everyone.values())


def test_the_api_carries_the_graph_only_when_asked(mysql_config, client):  # noqa: F811
    ids = _save_grid_sector(mysql_config, "West", 0, _cloud_local(1))
    base = f"/api/nav?from={pid('system', ids[0])}&to={pid('system', ids[-1])}"
    assert "graph" not in client.get(base).get_json()["route"]

    route = client.get(base + "&graph=1").get_json()["route"]
    graph = route["graph"]
    assert set(graph) == {"nodes", "edges", "more"}
    assert [e["from"] for e in graph["edges"] if e["on_route"]] == route["path"][:-1]
    assert [e["to"] for e in graph["edges"] if e["on_route"]] == route["path"][1:]
    assert all(set(e) == {"from", "to", "distance_ly", "on_route"} for e in graph["edges"])
    assert all(set(n) == {"id", "position_ly"} and len(n["position_ly"]) == 3 for n in graph["nodes"])
    assert all(isinstance(stop, str) for stop in graph["more"])
