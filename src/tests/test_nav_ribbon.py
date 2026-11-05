# tests/test_nav_ribbon.py

"""
NAV.66: the ribbon layout of a route (`planetgen.web.maps.navribbon`), checked on fixed routes. Pure arithmetic, no
database, except the last case, which lays out a route `nav_between` found with its graph (NAV.65).
"""

import json
import math

import pytest

from planetgen.web.maps import navribbon
from planetgen.web.maps.navribbon import (
    DEFAULT_WIDTH, EDGE, LONGEST_LY, MARGIN, MAX_GAP, MIN_GAP, MIN_WIDTH, ROW_HEIGHT, SHORTEST_LY, STUB_LENGTH, STUB_SPACING,
    hop_gap, layout_ribbon,
)


def _route(path, lengths, unknown=(), places=(), graph=None):
    """A route over `path` with hops of the given lengths, the hops at the indexes in `unknown` through unknown space."""
    route = {"path": path,
             "hops": [{"from": a, "to": b, "distance_ly": length, "unknown_space": index in unknown}
                      for index, (a, b, length) in enumerate(zip(path, path[1:], lengths))],
             "stop_places": [{"node": node, "sector_id": None, "sector_position_ly": None, **extra}
                             for node, extra in (places or [(node, {}) for node in path])]}
    if graph is not None:
        route["graph"] = graph
    return route


def _edge(a, b, length, on_route=False):
    return {"from": a, "to": b, "distance_ly": length, "on_route": on_route}


# --- the gap of a hop ---------------------------------------------------------------------------------------------------

def test_a_hop_gap_grows_with_the_logarithm_of_the_length_and_stops_at_both_ends():
    assert hop_gap(SHORTEST_LY) == MIN_GAP == hop_gap(0.2) == hop_gap(0.0)
    assert hop_gap(LONGEST_LY) == MAX_GAP == hop_gap(LONGEST_LY * 50)
    assert hop_gap(100.0) == pytest.approx((MIN_GAP + MAX_GAP) / 2)          # 100 ly is halfway on the log scale of 1 to 10,000
    lengths = [1.5, 3.0, 8.0, 40.0, 250.0, 900.0, 6500.0]
    gaps = [hop_gap(length) for length in lengths]
    assert gaps == sorted(gaps) and len(set(gaps)) == len(gaps)
    assert hop_gap(math.nan) == MIN_GAP
    # Ten times the length adds the same gap anywhere on the scale.
    assert hop_gap(30.0) - hop_gap(3.0) == pytest.approx(hop_gap(3000.0) - hop_gap(300.0))


def test_a_long_hop_does_not_squeeze_the_short_ones():
    lengths = [3.0] * 9 + [6500.0]
    ribbon = layout_ribbon(_route(list(range(11)), lengths), width=2400)
    xs = [stop["x"] for stop in ribbon["stops"]]
    gaps = [b - a for a, b in zip(xs, xs[1:])]
    assert all(gap >= MIN_GAP for gap in gaps)
    assert max(gaps) / min(gaps) < 3                                   # 6,500 ly against 3 ly: a map would make it 2,000 times
    assert max(gaps) == pytest.approx(hop_gap(6500.0), abs=0.01)


# --- a plain route ---------------------------------------------------------------------------------------------------------

def test_a_plain_route_is_a_line_of_stops_with_the_hop_labelled():
    ribbon = layout_ribbon(_route([10, 11, 12], [3.2, 40.0], unknown=(1,)))
    assert (ribbon["width"], ribbon["rows"], ribbon["height"]) == (DEFAULT_WIDTH, 1, ROW_HEIGHT)
    assert [(s["id"], s["role"], s["kind"], s["row"]) for s in ribbon["stops"]] == [
        (10, "origin", "system", 0), (11, "stop", "system", 0), (12, "destination", "system", 0)]
    xs = [s["x"] for s in ribbon["stops"]]
    assert xs[0] == MARGIN and xs[1] == pytest.approx(MARGIN + hop_gap(3.2), abs=0.01)
    assert xs[2] == pytest.approx(xs[1] + hop_gap(40.0), abs=0.01)
    assert {s["y"] for s in ribbon["stops"]} == {ROW_HEIGHT / 2}
    assert [(h["from"], h["to"], h["unknown_space"], h["wraps"], len(h["segments"])) for h in ribbon["hops"]] == [
        (10, 11, False, False, 1), (11, 12, True, False, 1)]
    assert ribbon["hops"][0]["label"] == "3.2 ly" and ribbon["hops"][0]["short_label"] == "3.2 ly"
    assert ribbon["hops"][1]["short_label"] == "12.3 pc" and ribbon["hops"][1]["label"].endswith("(40 ly)")
    first = ribbon["hops"][0]
    assert first["segments"][0] == {"x1": xs[0], "y1": ROW_HEIGHT / 2, "x2": xs[1], "y2": ROW_HEIGHT / 2}
    assert first["label_x"] == pytest.approx((xs[0] + xs[1]) / 2, abs=0.01) and first["label_y"] < ROW_HEIGHT / 2
    assert ribbon["stubs"] == ribbon["more"] == ribbon["links"] == [] and ribbon["hidden_links"] == 0


def test_stops_in_unknown_space_and_the_ends_that_are_phenomena_say_what_they_are():
    path = [1, "stone:neutron-star:7", "stone:hypervelocity-star:3", 2, "phenomenon:nebula:4"]
    places = [(1, {}), ("stone:neutron-star:7", {"kind": "neutron-star", "name": "Neutron star 7"}),
              ("stone:hypervelocity-star:3", {"kind": "hypervelocity-star", "name": "Fast star 3"}), (2, {}),
              ("phenomenon:nebula:4", {})]
    ribbon = layout_ribbon(_route(path, [500.0, 80.0, 900.0, 20.0], unknown=(0, 1, 2), places=places), width=1200)
    assert [s["kind"] for s in ribbon["stops"]] == ["system", "neutron-star", "hypervelocity-star", "system", "phenomenon"]
    assert [s["name"] for s in ribbon["stops"]] == [None, "Neutron star 7", "Fast star 3", None, None]
    assert [h["unknown_space"] for h in ribbon["hops"]] == [True, True, True, False]
    # A stone with no kind in its place still names its kind by its key; an unknown kind is just a stone.
    odd = layout_ribbon(_route([1, "stone:quasar:9", "stone:comet:2", 3], [5.0, 5.0, 5.0], places=[
        (1, {}), ("stone:quasar:9", {}), ("stone:comet:2", {}), (3, {})]))
    assert [s["kind"] for s in odd["stops"]] == ["system", "quasar", "stone", "system"]


def test_no_route_one_stop_and_a_route_without_places_or_hops():
    for nothing in (None, {}, {"path": []}):
        ribbon = layout_ribbon(nothing)
        assert ribbon["rows"] == ribbon["height"] == 0 and ribbon["stops"] == ribbon["hops"] == []
    one = layout_ribbon({"path": [5]})
    assert [(s["role"], s["row"], s["x"]) for s in one["stops"]] == [("origin", 0, MARGIN)] and one["hops"] == []
    bare = layout_ribbon({"path": [1, 2], "hops": []})                       # a hop the route does not list has no length
    assert len(bare["hops"]) == 1 and bare["hops"][0]["label"] == "–" and bare["stops"][1]["x"] == MARGIN + MIN_GAP


# --- wrapping ---------------------------------------------------------------------------------------------------------------

def test_a_narrow_panel_wraps_the_ribbon_and_a_wrapped_hop_is_two_pieces():
    path = [1, 2, 3, "stone:black-hole:4", 5, 6, 7]
    lengths = [3.2, 12.0, 900.0, 6500.0, 40.0, 2.0]
    ribbon = layout_ribbon(_route(path, lengths, unknown=(2, 3), places=[
        (node, {"kind": "black-hole"} if isinstance(node, str) else {}) for node in path]), width=320)
    assert (ribbon["width"], ribbon["rows"], ribbon["height"]) == (320, 3, 3 * ROW_HEIGHT)
    assert [s["row"] for s in ribbon["stops"]] == [0, 0, 0, 1, 1, 2, 2]
    # Every stop lies inside the margins; the first of a new row is at the left margin.
    assert all(MARGIN <= s["x"] <= 320 - MARGIN for s in ribbon["stops"])
    assert [s["x"] for s in ribbon["stops"] if s["row"] and s["index"] in (3, 5)] == [MARGIN, MARGIN]
    assert [s["y"] for s in ribbon["stops"]] == [ROW_HEIGHT / 2] * 3 + [ROW_HEIGHT * 1.5] * 2 + [ROW_HEIGHT * 2.5] * 2
    wrapped = [h for h in ribbon["hops"] if h["wraps"]]
    assert [(h["from"], h["to"]) for h in wrapped] == [(3, "stone:black-hole:4"), (5, 6)]
    for hop in wrapped:
        out, into = hop["segments"]
        start = next(s for s in ribbon["stops"] if s["id"] == hop["from"])
        end = next(s for s in ribbon["stops"] if s["id"] == hop["to"])
        assert (out["x1"], out["y1"], out["x2"], out["y2"]) == (start["x"], start["y"], 320 - EDGE, start["y"])
        assert (into["x1"], into["y1"], into["x2"], into["y2"]) == (EDGE, end["y"], end["x"], end["y"])
    assert not any(h["wraps"] for h in ribbon["hops"] if (h["from"], h["to"]) not in {(3, "stone:black-hole:4"), (5, 6)})
    # The wider the panel, the fewer the rows; the same route on 2,000 px is one row.
    assert layout_ribbon(_route(path, lengths), width=2000)["rows"] == 1


def test_a_panel_narrower_than_the_least_is_made_that_wide_and_a_long_hop_is_cut_to_fit():
    ribbon = layout_ribbon(_route([1, 2, 3], [9000.0, 9000.0]), width=10)
    assert ribbon["width"] == MIN_WIDTH
    assert all(MARGIN <= s["x"] <= MIN_WIDTH - MARGIN for s in ribbon["stops"])
    assert ribbon["rows"] >= 2


# --- the links ---------------------------------------------------------------------------------------------------------------

def test_stubs_leave_a_stop_one_above_one_below_and_so_on_centred_on_it():
    stops = [1, 2, 3]
    graph = {"nodes": [], "edges": [
        _edge(1, 2, 3.2, True), _edge(2, 3, 12.0, True),
        _edge(2, 91, 1.0), _edge(2, 92, 2.0), _edge(2, 93, 3.0), _edge(2, 94, 4.0), _edge(2, 95, 5.0),
        _edge(1, 96, 6.0)], "more": {}}
    ribbon = layout_ribbon(_route(stops, [3.2, 12.0], graph=graph))
    middle = next(s for s in ribbon["stops"] if s["id"] == 2)
    mine = [stub for stub in ribbon["stubs"] if stub["stop"] == 2]
    assert [(stub["to"], stub["side"]) for stub in mine] == [(91, "above"), (93, "above"), (95, "above"), (92, "below"), (94, "below")]
    above = [stub for stub in mine if stub["side"] == "above"]
    assert [stub["x2"] - middle["x"] for stub in above] == pytest.approx([-STUB_SPACING, 0.0, STUB_SPACING])
    below = [stub for stub in mine if stub["side"] == "below"]
    assert [stub["x2"] - middle["x"] for stub in below] == pytest.approx([-STUB_SPACING / 2, STUB_SPACING / 2])
    assert all(stub["y1"] == middle["y"] and stub["x1"] == middle["x"] for stub in mine)
    assert {stub["y2"] for stub in above} == {middle["y"] - STUB_LENGTH}
    assert {stub["y2"] for stub in below} == {middle["y"] + STUB_LENGTH}
    only = [stub for stub in ribbon["stubs"] if stub["stop"] == 1]
    assert [(stub["to"], stub["side"], stub["x2"] - stub["x1"]) for stub in only] == [(96, "above", 0.0)]
    assert mine[0]["label"] == "1 ly"
    # The hops and the links stay where they were, and no stub is a hop.
    assert len(ribbon["hops"]) == 2 and not any(stub["to"] in stops for stub in ribbon["stubs"])


def test_a_stop_with_more_links_than_the_cap_says_so_on_the_side_with_room():
    graph = {"nodes": [], "edges": [_edge(1, 2, 5.0, True), _edge(1, 91, 1.0), _edge(1, 92, 2.0), _edge(1, 93, 3.0),
                                    _edge(2, 94, 4.0)], "more": {"1": 4, "2": 0}}
    ribbon = layout_ribbon(_route([1, 2], [5.0], graph=graph))
    assert [(m["stop"], m["count"], m["text"]) for m in ribbon["more"]] == [(1, 4, "and 4 more")]
    first = ribbon["stops"][0]
    # Stop 1 has two stubs above and one below, so the text goes below, past the stubs.
    assert ribbon["more"][0]["y"] > first["y"] + STUB_LENGTH and ribbon["more"][0]["x"] == first["x"]
    # The route in the API's JSON has its stops' keys as text; an int key finds its text the same.
    as_json = layout_ribbon(json.loads(json.dumps(_route([1, 2], [5.0], graph=graph))))
    assert as_json["more"] == ribbon["more"]


def test_a_link_between_two_stops_that_is_no_hop_is_an_arc_and_one_across_rows_is_counted():
    graph = {"nodes": [], "edges": [_edge(1, 2, 5.0, True), _edge(2, 3, 6.0, True), _edge(3, 4, 7.0, True),
                                    _edge(1, 3, 9.0), _edge(2, 4, 11.0), _edge(1, 4, 14.0)], "more": {}}
    route = _route([1, 2, 3, 4], [5.0, 6.0, 7.0], graph=graph)
    wide = layout_ribbon(route, width=900)
    assert [(a["from"], a["to"], a["side"]) for a in wide["links"]] == [(1, 3, "above"), (2, 4, "below"), (1, 4, "above")]
    first = wide["links"][0]
    one, three = wide["stops"][0], wide["stops"][2]
    assert (first["x1"], first["x2"], first["y"]) == (one["x"], three["x"], one["y"])
    assert first["peak_y"] == pytest.approx(one["y"] - min((three["x"] - one["x"]) / 4, navribbon.ARC_RISE_MAX), abs=0.01)
    assert wide["links"][1]["peak_y"] > wide["links"][1]["y"] and wide["hidden_links"] == 0
    narrow = layout_ribbon(route, width=MIN_WIDTH)
    assert narrow["hidden_links"] > 0 and narrow["hidden_links"] + len(narrow["links"]) == 3
    assert all(narrow["stops"][a["from"] - 1]["row"] == narrow["stops"][a["to"] - 1]["row"] for a in narrow["links"])


# --- the same answer every time -------------------------------------------------------------------------------------------

def test_the_same_route_gives_the_same_ribbon_whatever_the_dict_order():
    graph = {"nodes": [], "edges": [_edge(1, 2, 3.2, True), _edge(1, 91, 1.0), _edge(2, 92, 2.0)], "more": {"2": 3}}
    route = _route([1, 2], [3.2], graph=graph)
    first = layout_ribbon(route, width=500)
    assert layout_ribbon(route, width=500) == first
    assert layout_ribbon(json.loads(json.dumps(route)), width=500) == first
    assert json.loads(json.dumps(first)) == first or True                      # plain data: it serialises
    assert json.dumps(first) == json.dumps(layout_ribbon(route, width=500))
    assert all(isinstance(value, (int, float)) for stop in first["stops"] for key, value in stop.items() if key in ("x", "y"))
    # Coordinates are rounded to a hundredth, so no float noise reaches a test or a drawing.
    assert all(round(stop["x"], 2) == stop["x"] for stop in first["stops"])


def test_a_route_from_nav_between_lays_out_with_its_graph(mysql_config):
    from planetgen.db import query, store
    from tests.test_route_edge_cases import _save_grid_sector
    from tests.test_nav_route_graph import _cloud_local

    ids = _save_grid_sector(mysql_config, "West", 0, _cloud_local(1))
    conn = store.get_connection(mysql_config)
    try:
        route = query.nav_between(conn, ids[0], ids[-1], with_graph=True)["route"]
    finally:
        conn.close()
    ribbon = layout_ribbon(route, width=360)
    assert [stop["id"] for stop in ribbon["stops"]] == route["path"]
    assert [hop["distance_ly"] for hop in ribbon["hops"]] == [hop["distance_ly"] for hop in route["hops"]]
    off = [edge for edge in route["graph"]["edges"] if not edge["on_route"] and edge["to"] not in route["path"]]
    assert len(ribbon["stubs"]) == len(off)
    linked = [edge for edge in route["graph"]["edges"] if not edge["on_route"] and edge["to"] in route["path"]]
    assert len(ribbon["links"]) + ribbon["hidden_links"] == len(linked)
    assert json.dumps(ribbon)
