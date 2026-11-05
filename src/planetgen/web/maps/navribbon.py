# planetgen/web/maps/navribbon.py

"""
The ribbon layout of a route (NAV.66)
=====================================

A schematic route map (layout L1 of reports/nav-schematic-route-map-2026-11-04): the route as one line of stops, left to
right, the links round it as short stubs above and below. Pure arithmetic over a route, no database and no I/O, so a
route always gives the same ribbon; `layout_ribbon` returns plain numbers and text that a drawing (NAV.67, an SVG with
the stop glyphs and the colours) places as they stand.

What the layout does and does not keep true:

- The order of the stops, the length of every hop (written on it), the flag for a hop through unknown space and the
  kind of each stop are kept. Where a stop lies in space is not: the ribbon is not a map.
- A hop's gap on the page follows the logarithm of its length (`hop_gap`), so one 6,500 ly hop does not squeeze nine
  hops of 3 ly to dots. The gap is the same for the same length on every route, so two ribbons compare.
- A ribbon wraps to a new row where the next stop would run past the panel's width (`width`, at least `MIN_WIDTH`); a hop
  that wraps is drawn as two pieces, one running off the right edge and one running in from the left edge, and says so
  (`wraps`).
- The links from `route.graph` (NAV.65) are the neighbour relation, not travel lanes: a link from a stop to a system off
  the route is a stub, alternating above and below the stop, at most as many as the graph lists (it is capped per stop),
  and the stops with more get an "and N more" text. A link between two stops of the route that is not a hop is an arc
  over or under the line when both stops are on one row, and is counted in `hidden_links` when they are not.

The input is `route` as `db.query.nav_between` returns it or `GET /api/nav` sends it (`path`, `hops`, `stop_places` and
optionally `graph`). A route without a `graph` has no stubs.
"""

import math

from planetgen.util.format import format_distance_ly

MIN_WIDTH = 240
"""int: The narrowest ribbon, pixels; a smaller `width` is raised to this (two stops with a hop and the margins)."""

DEFAULT_WIDTH = 720
"""int: The width of a ribbon when none is asked for, pixels."""

MARGIN = 28
"""int: Space at the left and right of the stops, pixels (the node and its label)."""

EDGE = 6
"""int: Where a hop that wraps leaves and enters the panel, pixels from its side."""

MIN_GAP = 68
MAX_GAP = 200
"""int: The least and the most pixels between two stops: a hop of `SHORTEST_LY` or less has the first, one of `LONGEST_LY`
or more the second."""

SHORTEST_LY = 1.0
LONGEST_LY = 10000.0
"""float: The hop lengths, light-years, the gap is scaled between (logarithmically)."""

STUB_LENGTH = 30
STUB_SPACING = 16
"""int: How far a stub reaches from the line and how far apart the stubs on one side of a stop lie, pixels."""

LABEL_ROOM = 22
"""int: Space above and below the stubs for their texts, pixels."""

ROW_HEIGHT = 2 * (STUB_LENGTH + LABEL_ROOM)
"""int: The height of one row of the ribbon, pixels; the line runs through its middle."""

ARC_RISE_MAX = STUB_LENGTH + 10
"""int: The most an arc over a link between two stops rises from the line, pixels."""

STONE_KINDS = ("bright-star", "black-hole", "neutron-star", "quasar", "hypervelocity-star")
"""tuple: The kinds of a stop in unknown space (`db.stones`); the kind of any other stop is `system` or `phenomenon`."""

_DIGITS = 2


def hop_gap(distance_ly):
    """
    The pixels between two stops joined by a hop of `distance_ly`: `MIN_GAP` for a hop of `SHORTEST_LY` or less, `MAX_GAP`
    for one of `LONGEST_LY` or more, between them in step with the logarithm of the length.
    """
    if not (distance_ly > SHORTEST_LY):          # also a missing or nan length
        return float(MIN_GAP)
    share = math.log(min(distance_ly, LONGEST_LY) / SHORTEST_LY) / math.log(LONGEST_LY / SHORTEST_LY)
    return MIN_GAP + (MAX_GAP - MIN_GAP) * share


def _round(value):
    return round(float(value), _DIGITS)


def _kind(node, place):
    """What a stop is: `system`, one of `STONE_KINDS`, or `phenomenon` (an end that is a nebula, a black hole and so on)."""
    if isinstance(node, str):
        if node.startswith("stone:"):
            kind = place.get("kind") or node.split(":")[1]
            return kind if kind in STONE_KINDS else "stone"
        if node.startswith("phenomenon:"):
            return "phenomenon"
    return "system"


def _labels(distance_ly):
    """`(label, short_label)` of a length: the site's text for it (`util.format.format_distance_ly`, the unit ladder, such as
    `12.3 pc (40 ly)`) and its first unit alone (`12.3 pc`), for a gap too narrow for the whole."""
    label = format_distance_ly(distance_ly)
    return label, label.split(" (")[0]


def _stub_slots(count):
    """The x offsets, pixels, of `count` stubs on one side of a stop, centered on it."""
    return [_round((slot - (count - 1) / 2) * STUB_SPACING) for slot in range(count)]


def layout_ribbon(route, width=DEFAULT_WIDTH):
    """
    Lays a route out as a ribbon.

    Args:
        route (dict or None): `path`, `hops` (`from`, `to`, `distance_ly`, `unknown_space`), `stop_places` (`node`, and
            `kind` and `name` on a stop in unknown space) and optionally `graph` (`nodes`, `edges`, `more`, as
            `nav_graph.route_graph` makes them); `None` or an empty path gives an empty ribbon.
        width (int): The panel's width, pixels; at least `MIN_WIDTH`.

    Returns:
        dict: `width` and `height` (pixels), `rows`, `stops`, `hops`, `stubs`, `more`, `links` and `hidden_links`:

        - `stops`: `{id, index, kind, role, name, row, x, y}` in path order; `role` is `origin`, `stop` or `destination`.
        - `hops`: `{from, to, distance_ly, label, short_label, unknown_space, wraps, segments, label_x, label_y}`; `segments` is one
          `{x1, y1, x2, y2}`, or two when the hop wraps (the first runs to the right edge, the second in from the left edge).
        - `stubs`: `{stop, to, distance_ly, label, short_label, side, x1, y1, x2, y2}`, `side` `above` or `below`.
        - `more`: `{stop, count, text, x, y}`, the stops whose links the graph's cap cut.
        - `links`: `{from, to, distance_ly, label, short_label, side, x1, x2, y, peak_y}` for each link between two stops of one row that
          is no hop; `hidden_links` counts those on two rows.
    """
    width = max(int(width), MIN_WIDTH)
    route = route or {}
    path = list(route.get("path") or [])
    places = {place["node"]: place for place in route.get("stop_places") or []}
    hop_by_start = {index: hop for index, hop in enumerate(route.get("hops") or [])}
    ribbon = {"width": width, "height": 0, "rows": 0, "stops": [], "hops": [], "stubs": [], "more": [], "links": [],
              "hidden_links": 0}
    if not path:
        return ribbon

    # The stops: left to right, on to the next row where the next one would run past the right margin.
    right = width - MARGIN
    row, x = 0, float(MARGIN)
    placed = []
    for index, node in enumerate(path):
        place = places.get(node, {})
        placed.append({"id": node, "index": index, "kind": _kind(node, place), "role": "stop", "name": place.get("name"),
                       "row": row, "x": x, "y": row * ROW_HEIGHT + ROW_HEIGHT / 2})
        if index + 1 < len(path):
            gap = min(hop_gap(hop_by_start.get(index, {}).get("distance_ly", math.nan)), right - MARGIN)
            if x + gap > right:
                row, x = row + 1, float(MARGIN)
            else:
                x += gap
    placed[0]["role"] = "origin"
    if len(placed) > 1:
        placed[-1]["role"] = "destination"
    where = {stop["id"]: stop for stop in placed}
    rows = placed[-1]["row"] + 1

    for stop in placed:
        ribbon["stops"].append({**stop, "x": _round(stop["x"]), "y": _round(stop["y"])})

    # The hops, a wrapped one as two pieces.
    for index, (a, b) in enumerate(zip(path, path[1:])):
        hop = hop_by_start.get(index, {"distance_ly": math.nan, "unknown_space": False})
        start, end = where[a], where[b]
        if start["row"] == end["row"]:
            segments = [(start["x"], start["y"], end["x"], end["y"])]
            label_x, label_y = (start["x"] + end["x"]) / 2, start["y"] - 8
        else:
            segments = [(start["x"], start["y"], width - EDGE, start["y"]), (EDGE, end["y"], end["x"], end["y"])]
            first_length, second_length = width - EDGE - start["x"], end["x"] - EDGE
            piece = segments[0] if first_length >= second_length else segments[1]
            label_x, label_y = (piece[0] + piece[2]) / 2, piece[1] - 8
        label, short_label = _labels(hop["distance_ly"])
        ribbon["hops"].append({
            "from": a, "to": b, "distance_ly": hop["distance_ly"], "label": label, "short_label": short_label,
            "unknown_space": bool(hop.get("unknown_space")), "wraps": start["row"] != end["row"],
            "segments": [{"x1": _round(x1), "y1": _round(y1), "x2": _round(x2), "y2": _round(y2)}
                         for x1, y1, x2, y2 in segments],
            "label_x": _round(label_x), "label_y": _round(label_y)})

    # The links round the stops.
    graph = route.get("graph") or {}
    on_route = set(path)
    stubs_of, between = {}, []
    for edge in graph.get("edges") or []:
        if edge.get("on_route"):
            continue
        if edge["to"] in on_route:
            between.append(edge)
        elif edge["from"] in where:
            stubs_of.setdefault(edge["from"], []).append(edge)
    for stop in placed:
        edges = stubs_of.get(stop["id"], [])
        sides = {"above": edges[0::2], "below": edges[1::2]}
        for side, chosen in sides.items():
            sign = -1 if side == "above" else 1
            for edge, dx in zip(chosen, _stub_slots(len(chosen))):
                label, short_label = _labels(edge["distance_ly"])
                ribbon["stubs"].append({
                    "stop": stop["id"], "to": edge["to"], "distance_ly": edge["distance_ly"],
                    "label": label, "short_label": short_label, "side": side,
                    "x1": _round(stop["x"]), "y1": _round(stop["y"]),
                    "x2": _round(stop["x"] + dx), "y2": _round(stop["y"] + sign * STUB_LENGTH)})
    ribbon["stubs"].sort(key=lambda stub: (where[stub["stop"]]["index"], stub["side"] != "above", stub["x2"]))

    left_over = {str(key): count for key, count in (graph.get("more") or {}).items()}
    for stop in placed:
        count = left_over.get(str(stop["id"]))
        if count:
            below = len(stubs_of.get(stop["id"], [])[1::2])
            above = len(stubs_of.get(stop["id"], [])[0::2])
            sign = 1 if below < above else -1          # the side with the fewer stubs has room
            ribbon["more"].append({"stop": stop["id"], "count": count, "text": f"and {count} more",
                                   "x": _round(stop["x"]), "y": _round(stop["y"] + sign * (STUB_LENGTH + LABEL_ROOM / 2))})

    arcs = 0
    for edge in between:
        a, b = where[edge["from"]], where[edge["to"]]
        if a["row"] != b["row"]:
            ribbon["hidden_links"] += 1
            continue
        side = "above" if arcs % 2 == 0 else "below"
        arcs += 1
        low, high = sorted((a, b), key=lambda stop: stop["x"])
        rise = min((high["x"] - low["x"]) / 4, ARC_RISE_MAX)
        label, short_label = _labels(edge["distance_ly"])
        ribbon["links"].append({
            "from": edge["from"], "to": edge["to"], "distance_ly": edge["distance_ly"],
            "label": label, "short_label": short_label, "side": side,
            "x1": _round(low["x"]), "x2": _round(high["x"]), "y": _round(a["y"]),
            "peak_y": _round(a["y"] + (-rise if side == "above" else rise))})

    ribbon["rows"] = rows
    ribbon["height"] = rows * ROW_HEIGHT
    return ribbon
