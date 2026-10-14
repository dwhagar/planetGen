# planetgen/web/nav_page.py

"""
The NAV page, `/nav` (was `nav.py`): course, distance and optimal route
between two endpoints, each a star system or a standalone phenomenon.
The course and route come from `GET /api/nav` (`queryDb.nav_between`);
this page only picks the endpoints and renders the result.

URL scheme (all GET, so every step is bookmarkable):

    /nav                                   choose a starting sector
    /nav?from_sector=5                     choose a starting system in it
    /nav?from=system:12                    choose a destination
    /nav?from=system:12&to_sector=9        ... a system in another sector
    /nav?from=system:12&to=system:40       the course and route
    /nav?from=nebula:3&to=system:40        a phenomenon endpoint

An endpoint is an object reference (`planetgen.galaxy.objectref`):
`system:<id>`, a body in a system (`star`, `planet`, `moon`, `belt`,
`comet`, NAV.16), or `<phenomenon type>:<id>` (`nebula`, `asteroid_field`,
`black_hole`, `neutron_star`, `supernova_remnant`, `rogue_planet`,
`interstellar_comet`, `quasar`). A bare number means a system. A course
to or from a body adds legs inside its system (see `queryDb.nav_course`).
`to` may be given without `from` ("navigate to here"): the origin picker
then carries it along, so choosing an origin lands on the course.

Build links with `nav_url(origin, destination)` / `endpoint(kind, id)`
below, or in a template `page_url('nav', **{'from': 'system:12'})`. The
older parameter style (`from_id`/`to_id`, or `from`/`to` with
`from_kind`/`to_kind="phenomenon"` and `from_type`/`to_type`, what
`nav.py` read and `page_url("nav", from_id=...)` produced) still works:
it redirects to the canonical URL.

A phenomenon or body is never offered in the pickers (there is no bounded,
dropdown-friendly list of every phenomenon); it becomes an endpoint only
through a link from its own page.
"""

import re
from urllib.parse import urlencode

from flask import redirect, request, url_for

from planetgen.web.lib import apiclient
from planetgen.web.lib.fmt import format_distance_ly
from planetgen.web.maps.navmap import render_nav_map_panel
from planetgen.galaxy import objectref
from planetgen.galaxy.navigation import format_course
from planetgen.physics.units import ly_to_pc

from . import bp
from .sector_page import PHENOMENON_TYPE_LABELS
from .helpers import crumb, current_admin, db_name, page_url, render_page, trusted_html

SECTOR_PICKER_LIMIT = 500
"""int: How many sectors the pickers offer: `GET /api/sectors`'s own
maximum page size (`docs/api.md`, "Pagination")."""

FRAME_LABELS = {
    "galactic": "Galactic Standard Frame: 000 mark 000 points at the galactic core",
    "sector": "Sector Local Frame: 000 mark 000 points at the sector's center",
    "system": "System Local Frame: 000 mark 000 points at the star",
}
"""dict: Each `navigation.Course.frame` value's line on the course panel."""

_LEGACY_PARAMS = ("from_id", "to_id", "from_kind", "to_kind", "from_type", "to_type")


def endpoint(kind, entity_id):
    """
    The `from`/`to` value for an endpoint: `endpoint("system", "FE81000A2B-0000005-000")` is
    `"system:FE81000A2B-0000005-000"`, `endpoint("planet", 3)` is `"planet:3"` (a phenomenon's
    kind is its type).
    """
    return objectref.format_public(kind, entity_id)


def nav_url(origin=None, destination=None):
    """`/nav` with `from=origin` and/or `to=destination` (both
    `endpoint(...)` strings, or `None`)."""
    return page_url("nav", **{"from": origin, "to": destination})


def parse_endpoint(raw):
    """
    Parses a `from`/`to` value.

    Returns:
        tuple: `(kind, id)` -- any `objectref` kind but a sector.

    Raises:
        apiclient.NotFoundError: For a value that is not an object
            reference (a bad id, the same 404 every page gives one).
    """
    try:
        kind, entity_id = objectref.parse_public(raw)
    except ValueError:
        kind = None
    if kind in (None, "sector"):
        raise apiclient.NotFoundError(f"No such system, body or phenomenon: {raw!r}")
    return kind, entity_id


def _legacy_redirect(args):
    """
    The canonical URL for a request using the old parameter style, or
    `None` when the request already uses the new one.
    """
    if not any(name in args for name in _LEGACY_PARAMS):
        return None
    params = []
    for prefix in ("from", "to"):
        raw = (args.get(f"{prefix}_id") or args.get(prefix) or "").strip()
        if not raw:
            continue
        if raw.isdecimal():  # an old row number names nothing since objects have IDs (API.23)
            continue
        params.append((prefix, raw))
    for name in ("from_sector", "to_sector"):
        if args.get(name):
            params.append((name, args[name]))
    return url_for("web.nav") + (f"?{urlencode(params, safe=':')}" if params else "")


def _sector_id(raw):
    """A sector's printed ID from a `from_sector`/`to_sector` value."""
    text = raw.strip().upper()
    if not re.fullmatch(r"[0-9A-F]{1,16}", text):
        raise apiclient.NotFoundError(f"No such sector: {raw!r}")
    return text


def _resolve(kind, entity_id):
    """
    One endpoint's display info: `ref`, `kind` (`"system"`, a body kind or
    `"phenomenon"`), `type` (phenomenon type or `None`), `id`, `name`,
    `url` (its page), `key` (the node id `GET /api/nav` uses for its
    system or phenomenon in `route.path`), `anchor_kind`/`anchor_id` (that
    system or phenomenon, which the map draws), `system_id`, `sector_id`
    (a system's or body's; `None` for a phenomenon) and `placed` (a
    phenomenon's galaxy position; `None` otherwise).
    """
    if kind == "system":
        system = apiclient.get_system(db_name(), entity_id)
        return {
            "ref": endpoint("system", system["id"]), "kind": "system", "type": None,
            "id": system["id"], "name": system["name"],
            "url": page_url("system", system_id=system["id"]), "key": system["id"],
            "anchor_kind": "system", "anchor_id": system["id"], "system_id": system["id"],
            "sector_id": system["sector_id"], "placed": None,
        }
    if kind in objectref.BODY_KINDS:
        body = apiclient.get_object(db_name(), endpoint(kind, entity_id))
        parents = {parent["kind"]: parent["ref"] for parent in body["parents"]}
        system_id = objectref.parse_public(parents["system"])[1]
        sector_ref = parents.get("sector")
        return {
            "ref": body["ref"], "kind": kind, "type": None, "id": entity_id, "name": body["name"],
            "url": page_url("system", system_id=system_id), "key": system_id,
            "anchor_kind": "system", "anchor_id": system_id, "system_id": system_id,
            "sector_id": objectref.parse_public(sector_ref)[1] if sector_ref else None, "placed": None,
        }
    detail = apiclient.get_phenomenon(db_name(), kind, entity_id)
    return {
        "ref": endpoint(kind, detail["id"]), "kind": "phenomenon", "type": kind, "id": detail["id"],
        "name": detail["name"],
        "url": page_url("phenomenon", phenomenon_type=kind, phenomenon_id=detail["id"]),
        # queryDb._phenomenon_nav_key's format.
        "key": f"phenomenon:{kind}:{detail['id']}",
        "anchor_kind": "phenomenon", "anchor_id": detail["id"], "system_id": None,
        "sector_id": None, "placed": detail.get("galactic_radius_pc") is not None,
    }


def _param_of(point):
    """The `from`/`to` value for a resolved endpoint."""
    return point["ref"]


def _sector_options():
    return [{"value": row["id"], "label": row["name"]}
            for row in apiclient.get_sectors(db_name(), limit=SECTOR_PICKER_LIMIT)["items"]]


def _system_options(systems, exclude=None):
    return [{"value": endpoint("system", row["id"]), "label": row["name"]}
            for row in systems if row["id"] != exclude]


def _picker(heading, field, label, options, hidden, button="Continue", start_over=None, empty=None):
    """One step of a picker: a GET form with one `<select>`."""
    return {
        "heading": heading, "field": field, "label": label, "options": options,
        "hidden": [(name, value) for name, value in hidden.items() if value],
        "button": button, "start_over": start_over, "empty": empty,
    }


def _origin_pickers(from_sector_raw, to_raw):
    """The sector-then-system origin picker, carrying `to` along."""
    hidden = {"to": to_raw}
    if not from_sector_raw:
        return [_picker("Choose a starting sector", "from_sector", "Sector", _sector_options(), hidden,
                        empty="No sectors have been generated yet.")]
    sector = apiclient.get_sector(db_name(), _sector_id(from_sector_raw))
    return [_picker(
        f"Choose a starting system in {sector['name']}", "from", "Starting system",
        _system_options(sector["systems"]), hidden, start_over=nav_url(destination=to_raw or None),
        empty="No systems are placed in this sector yet.",
    )]


def _destination_pickers(origin, to_sector_raw):
    """
    The destination pickers: for a system origin, every other system in
    its sector, plus (when that sector is galaxy-placed) the cross-sector
    sector-then-system picker; for a phenomenon origin, the cross-sector
    picker alone.
    """
    origin_param = _param_of(origin)
    hidden = {"from": origin_param}
    pickers = []
    cross_sector = True
    if origin["kind"] != "phenomenon":
        sector = apiclient.get_sector(db_name(), origin["sector_id"])
        cross_sector = sector["placed"]
        pickers.append(_picker(
            "Same-sector destination", "to", "Destination (same sector)",
            _system_options(sector["systems"], exclude=origin["system_id"]), hidden, button="Plot course",
            empty="No other systems are placed in this sector yet.",
        ))
    cross_heading = "Cross-sector destination" if origin["kind"] != "phenomenon" else "Destination"
    if not cross_sector:
        pickers.append(_picker(cross_heading, None, None, [], {}, empty=(
            "This sector has no galaxy placement, so NAV is only available to other systems in this "
            "same sector.")))
    elif not to_sector_raw:
        options = [row for row in _sector_options() if row["value"] != origin["sector_id"]]
        pickers.append(_picker(cross_heading + ": choose a sector", "to_sector", "Sector", options, hidden,
                               empty="No other galaxy-placed sectors exist yet."))
    else:
        sector = apiclient.get_sector(db_name(), _sector_id(to_sector_raw))
        pickers.append(_picker(
            f"{cross_heading}: choose a system in {sector['name']}", "to", "Destination system",
            _system_options(sector["systems"]), hidden, button="Plot course",
            start_over=nav_url(origin_param), empty="No systems are placed in this sector yet.",
        ))
    return pickers


def _map_picks(pick, other):
    """
    The "pick on a map" links beside the pickers (MAP.22, design doc
    section 9): choose the `pick` endpoint (`"from"`/`"to"`) on the
    Galaxy Map (`/galaxy?pick=...`, whose stage 8 hands a sector click to
    the Sector Map's pick mode) or, once the `other` endpoint is a known
    system in a sector, straight in that sector's Sector Map pick mode
    (`/sector/<id>?pick=...`). Both carry `other` along.

    Args:
        pick (str): `"from"` or `"to"`.
        other (dict or None): The other endpoint, `_resolve`d.
    """
    other_param = _param_of(other) if other else None
    keep = {("to" if pick == "from" else "from"): other_param}
    links = []
    if other is not None and other["sector_id"] is not None:
        links.append({"label": "Pick in this sector",
                      "url": page_url("sector", sector_id=other["sector_id"], pick=pick, **keep)})
    links.append({"label": "Pick on Galaxy Map", "url": page_url("galaxy", pick=pick, **keep)})
    return links


def _bookmark_pick(pick, other_param, label="Bookmarks"):
    """
    The NAV page's Bookmarks select (MAP.22, MAP.23, NAV.40; design doc
    section 9): `static/bookmarks.js` fills it with the browser's system
    and phenomenon bookmarks (and sectors, which open their system
    picker, and saved map views, which open the Galaxy Map there to pick
    on) and goes on with the `pick` endpoint set and `other_param`, the
    other endpoint's value, kept. `label` names the select.
    """
    keep = "to" if pick == "from" else "from"
    return {"pick": pick, "keep_name": keep, "keep_value": other_param or "", "db": db_name(),
            "nav_url": url_for("web.nav"), "label": label}


def _route_names(route, known):
    """`{node key: name}` for every stop on the route, looking up only
    intermediate systems (a phenomenon is only ever the first or last
    stop, already in `known`)."""
    names = dict(known)
    for node in (route or {}).get("path") or []:
        if node not in names and not (isinstance(node, str) and node.startswith("phenomenon:")):
            names[node] = apiclient.get_system(db_name(), node)["name"]
    return names


def _route_stops(route, names):
    """The route's stops as `{name, url}`, in order."""
    stops = []
    for node in route["path"]:
        if isinstance(node, str) and node.startswith("phenomenon:"):
            _, phenomenon_type, phenomenon_id = node.split(":", 2)
            url = page_url("phenomenon", phenomenon_type=phenomenon_type, phenomenon_id=int(phenomenon_id))
        else:
            url = page_url("system", system_id=node)
        stops.append({"name": names.get(node, str(node)), "url": url})
    for stop, hop in zip(stops, route.get("hops", [])):
        # UX.35: the hop to the next stop, shown between the two; NAV.11: and how long it takes.
        stop["hop_text"] = format_distance_ly(hop["distance_ly"])
        if "bearing_deg" in hop:  # NAV.42: the course to the next stop, "045 mark 012"
            stop["hop_course"] = format_course(hop["bearing_deg"], hop["mark_deg"])
        stop["hop_unknown"] = bool(hop.get("unknown_space"))
        stop["hop_times"] = [f"Warp {leg['warp_factor']:g}: {leg['formatted']}"
                             for leg in hop.get("warp_times", []) if leg["warp_factor"] in READOUT_WARP_FACTORS]
    return stops


MAX_STAY_MINUTES = 10_000_000.0
"""float: The longest stay per stop the NAV page asks the API for (`api.routes.NAV_MAX_STAY_MINUTES`)."""


def _stay_minutes(args):
    """
    NAV.11: the `stay` URL parameter, minutes spent at each stop of a
    route. Returns `(minutes, text, error)`: `text` is what the field
    shows, `error` a message when the value was unusable (the stay is then
    0, the default).
    """
    raw = (args.get("stay") or "").strip()
    if not raw:
        return 0.0, "", None
    try:
        minutes = float(raw)
    except ValueError:
        minutes = -1.0
    if not 0.0 <= minutes <= MAX_STAY_MINUTES:
        return 0.0, raw, f"The stay per stop must be a number of minutes from 0 to {MAX_STAY_MINUTES:,.0f}; 0 was used."
    return minutes, raw, None


ROUTE_COLLAPSE_AT = 9
"""int: A route with more stops than this shows a short list and keeps the whole one in a `<details>` (UX.35)."""


def _route_short_list(stops, route):
    """
    UX.35: the stops of a long route worth showing at once -- the first, the
    last three, both ends of the longest hop and of every hop through
    unknown space -- each with `gap_before` set when stops were left out
    ahead of it. A route of `ROUTE_COLLAPSE_AT` stops or fewer gives `None`
    (the whole list is shown).
    """
    if len(stops) <= ROUTE_COLLAPSE_AT:
        return None
    keep = {0, *range(len(stops) - 3, len(stops))}
    hops = route.get("hops", [])
    for index, hop in enumerate(hops):
        if hop.get("unknown_space"):
            keep.update((index, index + 1))
    if hops:
        longest = max(range(len(hops)), key=lambda index: hops[index]["distance_ly"])
        keep.update((longest, longest + 1))
    short, previous = [], -1
    for index in sorted(keep):
        short.append({**stops[index], "gap_before": index != previous + 1})
        previous = index
    return short


def _waypoints(origin, destination, result, names):
    """The origin, the route's intermediate hops, then the destination,
    in `navmap.render_nav_map_panel`'s input shape."""
    points = [{"id": origin["anchor_id"], "kind": origin["anchor_kind"], "type": origin["type"], "name": origin["name"],
               "position": result["origin_position"], "role": "origin"}]
    route = result["route"]
    if route is not None:
        hops = route.get("hops", [])

        def to_next(index):  # NAV.42: the course to the next stop, for the map's tooltip
            hop = hops[index] if index < len(hops) else None
            return format_course(hop["bearing_deg"], hop["mark_deg"]) if hop and "bearing_deg" in hop else None

        points[0]["course"] = to_next(0)
        for index, node in enumerate(route["path"][1:-1], start=1):
            points.append({"id": node, "kind": "system", "type": None, "name": names.get(node, str(node)),
                           "position": route["positions"][str(node)], "role": "hop", "course": to_next(index)})
    points.append({"id": destination["anchor_id"], "kind": destination["anchor_kind"], "type": destination["type"],
                   "name": destination["name"], "position": result["destination_position"],
                   "role": "destination"})
    return points


LEG_LABELS = {
    "within": "Within the system",
    "out": "Out of the system",
    "between": "Between systems",
    "into": "Into the system",
}
"""dict: Each `nav_course` leg kind's row label."""


def _leg_rows(legs, origin, destination, names):
    """The legs table's rows: label, the two ends' names, distance and course."""
    def name(ref):
        if ref == origin["ref"]:
            return origin["name"]
        if ref == destination["ref"]:
            return destination["name"]
        kind, entity_id = objectref.parse_public(ref)
        if kind == "system":
            return names.get(entity_id) or apiclient.get_system(db_name(), entity_id)["name"]
        return ref

    rows = []
    for leg in legs:
        direct = leg["direct"]
        rows.append({
            "label": LEG_LABELS[leg["kind"]], "from": name(leg["from"]), "to": name(leg["to"]),
            "distance": format_distance_ly(direct["distance_ly"]),
            "course": format_course(direct["bearing_deg"], direct["mark_deg"]),
            "frame": FRAME_LABELS.get(direct["frame"], direct["frame"]),
        })
    return rows


def galaxy_map_url(origin, destination):
    """`/galaxy` showing this course (`?course=<from>,<to>`, read by
    `web/galaxy_views.py`)."""
    return page_url("galaxy", course=f"{origin},{destination}")


def galaxy_course(from_raw, to_raw):
    """
    A course as the Galaxy Map draws it (the drill-down's section 9.4):
    its waypoints in galaxy-frame parsecs, so the map can open the
    smallest stage holding them and draw the line.

    Args:
        from_raw (str): The origin's `from` value (`endpoint(...)`).
        to_raw (str): The destination's `to` value.

    Returns:
        dict or None: `scope` (`"galaxy"` or `"sector"`), `points`
            (`[{name, url, role, x, y, z}]`, parsecs; only for
            `"galaxy"`), `sector` (`{id, name, ring, layer, slot}`, only
            for `"sector"`, where the whole course sits inside one
            sector) and `navUrl` (back to the course on the NAV page).
            `None` when either endpoint can't be navigated, or the two
            can't be navigated together.

    Raises:
        apiclient.NotFoundError: For an endpoint that doesn't exist.
    """
    origin = _resolve(*parse_endpoint(from_raw))
    destination = _resolve(*parse_endpoint(to_raw))
    if _unavailable_reason(origin) or _unavailable_reason(destination):
        return None
    try:
        result = apiclient.get_nav(db_name(), origin["ref"], destination["ref"])
    except apiclient.ApiError as exc:
        if isinstance(exc, apiclient.NotFoundError) or exc.status_code != 400:
            raise
        return None  # the two endpoints exist but can't be navigated together
    nav_url_here = nav_url(_param_of(origin), _param_of(destination))
    if result["scope"] == "system":
        return None  # inside one system: nothing to draw on the Galaxy Map
    if result["scope"] != "galaxy":
        # One sector holds the whole course: the map shows that sector.
        sector = apiclient.get_sector(db_name(), origin["sector_id"])
        if sector.get("ring_index") is None:
            return None
        return {
            "scope": "sector", "points": [], "navUrl": nav_url_here,
            "sector": {
                "id": sector["id"], "name": sector["name"], "ring": sector["ring_index"],
                "layer": sector["layer_index"], "slot": sector["ring_slot_index"],
            },
        }
    names = _route_names(result["route"], {origin["key"]: origin["name"], destination["key"]: destination["name"]})
    points = []
    for point in _waypoints(origin, destination, result, names):
        x, y, z = point["position"]
        points.append({
            "name": point["name"], "role": point["role"], "url": _point_url(point),
            "x": ly_to_pc(x), "y": ly_to_pc(y), "z": ly_to_pc(z),
        })
    direct = result["direct"]
    return {
        "scope": "galaxy", "points": points, "sector": None, "navUrl": nav_url_here,
        # NAV.20: the straight line apart from the route, and the readout beside the map.
        "direct": [points[0], points[-1]] if len(points) > 2 else None,
        "readout": _course_readout(result, direct, route=result["route"], stops=_route_stops(result["route"], names)
                                   if result["route"] and len(points) > 2 else []),
    }


READOUT_WARP_FACTORS = (1, 9)
"""tuple: The warp factors whose travel times the Galaxy Map's course
readout shows."""


CHART_ERROR_PREFIX = re.compile(r"^planetGen API error \(\d+\): ")


@bp.route("/nav/chart")
def nav_chart():
    """
    NAV.48 (admins): the uncharted sectors that block a course. Shows what charting would take (the cells
    of the route's hops through unknown space, whether a route through charted space exists, a one-cell
    border to tick), and a form that sends the job to the Generate page, which shows the size and time and
    asks to confirm before it starts. A visitor who is not logged in is sent to log in first.
    """
    admin = current_admin()
    here = request.full_path.rstrip("?")
    if admin is None:
        return redirect(page_url("login", next=here), code=302)
    if admin.get("must_change_credentials"):
        return redirect(page_url("account", next=here), code=302)

    from_raw = (request.args.get("from") or "").strip()
    to_raw = (request.args.get("to") or "").strip()
    if not from_raw or not to_raw:
        return redirect(page_url("nav"), code=302)
    origin = _resolve(*parse_endpoint(from_raw))
    destination = _resolve(*parse_endpoint(to_raw))
    border = request.args.get("border") == "1"
    title = f"Chart the course: {origin['name']} \u2192 {destination['name']}"
    crumbs = [crumb("Nav", "nav"), crumb(f"{origin['name']} \u2192 {destination['name']}", "nav",
                                          **{"from": _param_of(origin), "to": _param_of(destination)}),
              crumb("Chart")]
    course_url = nav_url(_param_of(origin), _param_of(destination))
    page = {"section": "nav", "description": "Generate the uncharted sectors that block a course."}
    try:
        plan = apiclient.get_nav_chart(db_name(), origin["ref"], destination["ref"], border=border)
    except apiclient.NotFoundError:
        raise
    except apiclient.ApiError as exc:
        if exc.status_code != 400:
            raise
        return render_page("nav_chart.html", title=title, breadcrumbs=crumbs, origin=origin, destination=destination,
                           error=CHART_ERROR_PREFIX.sub("", str(exc)), course_url=course_url, **page)
    bypass = plan.get("bypass")
    return render_page(
        "nav_chart.html", title=title, breadcrumbs=crumbs, origin=origin, destination=destination, plan=plan,
        border=border, course_url=course_url, needs_confirm=plan["count"] > plan["confirm_over"],
        bypass_text=_bypass_text(bypass, plan), generate_url=page_url("generate"),
        from_param=_param_of(origin), to_param=_param_of(destination), **page)


def _bypass_text(bypass, plan):
    """NAV.48: what the bypass test found, in a sentence (empty when no hop is unknown)."""
    if not bypass:
        return ""
    if bypass["found"]:
        return (f"A route through charted space exists, {format_distance_ly(bypass['distance_ly'])} against "
                f"{format_distance_ly(plan['route_distance_ly'])} for the shortest route, so charting is optional.")
    if not bypass["checked"]:
        return "Whether a route through charted space exists could not be worked out in time."
    return "No route through charted space exists, so these sectors have to be generated before the course can avoid unknown space."


def _course_readout(result, direct, route, stops):
    """
    The course readout beside the Galaxy Map (NAV.20): `distance`,
    `course` ("045 mark 000"), `frame`, `route_distance` (or `None`),
    `times` (`[{label, text}]`) and `stops` (`[{name, url}]`, the route's
    hops between the two ends).
    """
    return {
        "distance": format_distance_ly(direct["distance_ly"]),
        "course": format_course(direct["bearing_deg"], direct["mark_deg"]),
        "frame": FRAME_LABELS.get(direct["frame"], direct["frame"]),
        "route_distance": format_distance_ly(route["distance_ly"]) if route and stops else None,
        "times": [{"label": f"Warp {leg['warp_factor']:g}", "text": leg["formatted"]}
                  for leg in result["warp_times"] if leg["warp_factor"] in READOUT_WARP_FACTORS],
        "stops": stops[1:-1],
    }


def _point_url(point):
    """A waypoint's own page."""
    if point.get("kind") == "phenomenon":
        return page_url("phenomenon", phenomenon_type=point["type"], phenomenon_id=point["id"])
    return page_url("system", system_id=point["id"])


def _unavailable_reason(origin):
    if origin["kind"] == "phenomenon" and not origin["placed"]:
        return "NAV is not available: this phenomenon has not been placed in the galaxy."
    if origin["kind"] != "phenomenon" and origin["sector_id"] is None:
        return "NAV is not available: this system isn't assigned to a sector."
    return None


@bp.route("/nav")
def nav():
    """The NAV page; see the module docstring for its URL scheme."""
    canonical = _legacy_redirect(request.args)
    if canonical is not None:
        return redirect(canonical, code=301)

    from_raw = (request.args.get("from") or "").strip()
    to_raw = (request.args.get("to") or "").strip()
    page = {"section": "nav", "description": "Plot a course between two star systems or phenomena."}

    if not from_raw:
        destination = _resolve(*parse_endpoint(to_raw)) if to_raw else None  # a bad `to` is a 404 now
        return render_page(
            "nav.html", title="Nav", breadcrumbs=[crumb("Nav")],
            pickers=_origin_pickers(request.args.get("from_sector", ""), to_raw),
            map_picks=_map_picks("from", destination), map_picks_heading="Or pick a start on a map",
            landing=not request.args.get("from_sector"),
            bookmark_pick=_bookmark_pick("from", to_raw), **page,
        )

    origin = _resolve(*parse_endpoint(from_raw))
    crumbs = [crumb("Nav", "nav"), crumb(origin["name"])]
    title = f"Nav: {origin['name']}"
    reason = _unavailable_reason(origin)
    if reason:
        return render_page("nav.html", title=title, breadcrumbs=crumbs, origin=origin, error=reason, **page)

    if not to_raw:
        return render_page(
            "nav.html", title=title, breadcrumbs=crumbs, origin=origin,
            pickers=_destination_pickers(origin, request.args.get("to_sector", "")),
            map_picks=_map_picks("to", origin), map_picks_heading="Or pick a destination on a map",
            bookmark_pick=_bookmark_pick("to", _param_of(origin)), **page,
        )

    to_kind, to_id = parse_endpoint(to_raw)
    destination = _resolve(to_kind, to_id)
    stay, stay_text, stay_error = _stay_minutes(request.args)
    try:
        result = apiclient.get_nav(db_name(), origin["ref"], destination["ref"], stay=stay or None)
    except apiclient.NotFoundError:
        raise
    except apiclient.ApiError as exc:
        if exc.status_code != 400:
            raise
        # The two endpoints exist but do not support NAV together (e.g.
        # different sectors without galaxy placement): a message on the
        # page, not an error page.
        message = re.sub(r"^planetGen API error \(\d+\): ", "", str(exc))
        return render_page("nav.html", title=title, breadcrumbs=crumbs, origin=origin, error=message,
                           start_over=nav_url(_param_of(origin)), **page)

    # UX.70: the page is the whole course, not just where it starts.
    title = f"Course: {origin['name']} \u2192 {destination['name']}"
    crumbs = [crumb("Nav", "nav"), crumb(f"{origin['name']} \u2192 {destination['name']}")]
    route = result["route"]
    names = _route_names(route, {origin["key"]: origin["name"], destination["key"]: destination["name"]})
    direct = result["direct"]
    in_system = result["scope"] == "system"
    map_html = "" if in_system else render_nav_map_panel(
        page_url, _waypoints(origin, destination, result, names), has_route=route is not None)
    return render_page(
        "nav.html", title=title, breadcrumbs=crumbs, origin=origin, destination=destination,
        direct=direct, distance_text=format_distance_ly(direct["distance_ly"]),
        route_distance_text=format_distance_ly(route["distance_ly"]) if route else None,
        course=format_course(direct["bearing_deg"], direct["mark_deg"]),
        frame_label=FRAME_LABELS.get(direct["frame"], direct["frame"]), warp_times=result["warp_times"], fold_times=result["fold_times"],
        scope_label={"system": "Inside one system", "sector": "Same sector"}.get(
            result["scope"], "Cross-sector (galaxy)"),
        legs=_leg_rows(result["legs"], origin, destination, names), in_system=in_system, note=result["note"],
        route=route, stops=_route_stops(route, names) if route and route["path"] else [],
        stay_text=stay_text, stay_error=stay_error, stay_minutes=stay,
        from_param=_param_of(origin), to_param=_param_of(destination),
        short_stops=_route_short_list(_route_stops(route, names), route) if route and route["path"] else None,
        longest_hop_text=format_distance_ly(route["longest_hop_ly"]) if route and "longest_hop_ly" in route else None,
        unknown_jumps=sum(1 for hop in route.get("hops", []) if hop.get("unknown_space")) if route else 0,
        chart_url=page_url("nav_chart", **{"from": _param_of(origin), "to": _param_of(destination)})
        if current_admin() is not None and route and any(hop.get("unknown_space") for hop in route.get("hops", []))
        else None,
        map_html=trusted_html(map_html),
        reverse_url=nav_url(_param_of(destination), _param_of(origin)),
        bookmark_changes=[_bookmark_pick("from", _param_of(destination), "New start"),
                          _bookmark_pick("to", _param_of(origin), "New destination")],
        galaxy_map_url=galaxy_map_url(_param_of(origin), _param_of(destination)),
        start_over=nav_url(_param_of(origin)),
        **page,
    )
