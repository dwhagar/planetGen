# html/web/nav_page.py

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

An endpoint is `<kind>:<id>`: `system:<id>` for a star system, or
`<phenomenon type>:<id>` (`nebula`, `asteroid_field`, `black_hole`,
`neutron_star`, `supernova_remnant`, `rogue_planet`,
`interstellar_comet`, `quasar`) for a phenomenon. A bare number means a system.
`to` may be given without `from` ("navigate to here"): the origin picker
then carries it along, so choosing an origin lands on the course.

Build links with `nav_url(origin, destination)` / `endpoint(kind, id)`
below, or in a template `page_url('nav', **{'from': 'system:12'})`. The
older parameter style (`from_id`/`to_id`, or `from`/`to` with
`from_kind`/`to_kind="phenomenon"` and `from_type`/`to_type`, what
`nav.py` read and `page_url("nav", from_id=...)` produced) still works:
it redirects to the canonical URL.

A phenomenon is never offered in the pickers (there is no bounded,
dropdown-friendly list of every phenomenon); it becomes an endpoint only
through a link from its own page.
"""

import re
from urllib.parse import urlencode

from flask import redirect, request, url_for

import apiclient
from fmt import format_distance_ly
from navmap import render_nav_map_panel
from stellarObjects.navigation import format_course
from stellarObjects.utils import ly_to_pc

from . import bp
from .sector_page import PHENOMENON_TYPE_LABELS
from .helpers import crumb, db_name, page_url, render_page, trusted_html

SECTOR_PICKER_LIMIT = 500
"""int: How many sectors the pickers offer: `GET /api/sectors`'s own
maximum page size (`docs/api.md`, "Pagination")."""

FRAME_LABELS = {
    "galactic": "Galactic Standard Frame: 000 mark 000 points at the galactic core",
    "sector": "Sector Local Frame: 000 mark 000 points at the sector's center",
    "system": "System Local Frame: 000 mark 000 points at the star",
}
"""dict: Each `navigation.Course.frame` value's line on the course panel."""

_ENDPOINT_RE = re.compile(r"^(?:([a-z_]+):)?(\d+)$")

_LEGACY_PARAMS = ("from_id", "to_id", "from_kind", "to_kind", "from_type", "to_type")


def endpoint(kind, entity_id):
    """
    The `from`/`to` value for an endpoint: `endpoint("system", 12)` is
    `"system:12"`, `endpoint("nebula", 3)` is `"nebula:3"` (a phenomenon's
    kind is its type).
    """
    return f"{kind}:{int(entity_id)}"


def nav_url(origin=None, destination=None):
    """`/nav` with `from=origin` and/or `to=destination` (both
    `endpoint(...)` strings, or `None`)."""
    return page_url("nav", **{"from": origin, "to": destination})


def parse_endpoint(raw):
    """
    Parses a `from`/`to` value.

    Returns:
        tuple: `(kind, id)` -- kind `"system"` or a phenomenon type.

    Raises:
        apiclient.NotFoundError: For a value that is not `<kind>:<id>`
            (a bad id, the same 404 every page gives one).
    """
    match = _ENDPOINT_RE.match(raw.strip())
    if not match or match.group(1) not in (None, "system", *PHENOMENON_TYPE_LABELS):
        raise apiclient.NotFoundError(f"No such system or phenomenon: {raw!r}")
    return (match.group(1) or "system"), int(match.group(2))


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
        if raw.isdecimal():  # not isdigit(): "²" is a digit int() can't read
            kind = args.get(f"{prefix}_type") if args.get(f"{prefix}_kind") == "phenomenon" else None
            raw = endpoint(kind or "system", raw)
        params.append((prefix, raw))
    for name in ("from_sector", "to_sector"):
        if args.get(name):
            params.append((name, args[name]))
    return url_for("web.nav") + (f"?{urlencode(params, safe=':')}" if params else "")


def _sector_id(raw):
    try:
        return int(raw)
    except ValueError:
        raise apiclient.NotFoundError(f"No such sector: {raw!r}")


def _resolve(kind, entity_id):
    """
    One endpoint's display info: `kind`, `type` (phenomenon type or
    `None`), `id`, `name`, `url` (its page), `key` (the node id
    `GET /api/nav` uses for it in `route.path`), `sector_id` (a system's;
    `None` for a phenomenon) and `placed` (a phenomenon's galaxy
    position; `None` for a system).
    """
    if kind == "system":
        system = apiclient.get_system(db_name(), entity_id)
        return {
            "kind": "system", "type": None, "id": system["id"], "name": system["name"],
            "url": page_url("system", system_id=system["id"]), "key": system["id"],
            "sector_id": system["sector_id"], "placed": None,
        }
    detail = apiclient.get_phenomenon(db_name(), kind, entity_id)
    return {
        "kind": "phenomenon", "type": kind, "id": detail["id"], "name": detail["name"],
        "url": page_url("phenomenon", phenomenon_type=kind, phenomenon_id=detail["id"]),
        # queryDb._phenomenon_nav_key's format.
        "key": f"phenomenon:{kind}:{detail['id']}",
        "sector_id": None, "placed": detail.get("galactic_radius_pc") is not None,
    }


def _param_of(point):
    """The `from`/`to` value for a resolved endpoint."""
    return endpoint(point["type"] or "system", point["id"])


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
    if origin["kind"] == "system":
        sector = apiclient.get_sector(db_name(), origin["sector_id"])
        cross_sector = sector["placed"]
        pickers.append(_picker(
            "Same-sector destination", "to", "Destination (same sector)",
            _system_options(sector["systems"], exclude=origin["id"]), hidden, button="Plot course",
            empty="No other systems are placed in this sector yet.",
        ))
    cross_heading = "Cross-sector destination" if origin["kind"] == "system" else "Destination"
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


def _bookmark_pick(pick, other_param):
    """
    The NAV page's Bookmarks select (MAP.22, MAP.23; design doc section
    9): `static/bookmarks.js` fills it with the browser's system and
    phenomenon bookmarks (and sectors, which open their system picker)
    and goes to `/nav` with the `pick` endpoint set and `other_param`, the
    other endpoint's value, kept.
    """
    keep = "to" if pick == "from" else "from"
    return {"pick": pick, "keep_name": keep, "keep_value": other_param or "", "db": db_name(),
            "nav_url": url_for("web.nav")}


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
    return stops


def _waypoints(origin, destination, result, names):
    """The origin, the route's intermediate hops, then the destination,
    in `navmap.render_nav_map_panel`'s input shape."""
    points = [{"id": origin["id"], "kind": origin["kind"], "type": origin["type"], "name": origin["name"],
               "position": result["origin_position"], "role": "origin"}]
    route = result["route"]
    if route is not None:
        for node in route["path"][1:-1]:
            points.append({"id": node, "kind": "system", "type": None, "name": names.get(node, str(node)),
                           "position": route["positions"][str(node)], "role": "hop"})
    points.append({"id": destination["id"], "kind": destination["kind"], "type": destination["type"],
                   "name": destination["name"], "position": result["destination_position"],
                   "role": "destination"})
    return points


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
        result = apiclient.get_nav(
            db_name(), origin["id"], destination["id"],
            from_kind=origin["kind"], to_kind=destination["kind"],
            from_type=origin["type"], to_type=destination["type"],
        )
    except apiclient.ApiError as exc:
        if isinstance(exc, apiclient.NotFoundError) or exc.status_code != 400:
            raise
        return None  # the two endpoints exist but can't be navigated together
    nav_url_here = nav_url(_param_of(origin), _param_of(destination))
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
    return {"scope": "galaxy", "points": points, "sector": None, "navUrl": nav_url_here}


def _point_url(point):
    """A waypoint's own page."""
    if point.get("kind") == "phenomenon":
        return page_url("phenomenon", phenomenon_type=point["type"], phenomenon_id=point["id"])
    return page_url("system", system_id=point["id"])


def _unavailable_reason(origin):
    if origin["kind"] == "phenomenon" and not origin["placed"]:
        return "NAV is not available: this phenomenon has not been placed in the galaxy."
    if origin["kind"] == "system" and origin["sector_id"] is None:
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
    try:
        result = apiclient.get_nav(
            db_name(), origin["id"], to_id,
            from_kind=origin["kind"], to_kind="system" if to_kind == "system" else "phenomenon",
            from_type=origin["type"], to_type=None if to_kind == "system" else to_kind,
        )
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

    destination = _resolve(to_kind, to_id)
    route = result["route"]
    names = _route_names(route, {origin["key"]: origin["name"], destination["key"]: destination["name"]})
    direct = result["direct"]
    map_html = render_nav_map_panel(page_url, _waypoints(origin, destination, result, names),
                                    has_route=route is not None)
    return render_page(
        "nav.html", title=title, breadcrumbs=crumbs, origin=origin, destination=destination,
        direct=direct, distance_text=format_distance_ly(direct["distance_ly"]),
        route_distance_text=format_distance_ly(route["distance_ly"]) if route else None,
        course=format_course(direct["bearing_deg"], direct["mark_deg"]),
        frame_label=FRAME_LABELS.get(direct["frame"], direct["frame"]), warp_times=result["warp_times"], fold_times=result["fold_times"],
        scope_label="Same sector" if result["scope"] == "sector" else "Cross-sector (galaxy)",
        route=route, stops=_route_stops(route, names) if route and route["path"] else [],
        map_html=trusted_html(map_html),
        reverse_url=nav_url(_param_of(destination), _param_of(origin)),
        galaxy_map_url=galaxy_map_url(_param_of(origin), _param_of(destination)),
        start_over=nav_url(_param_of(origin)),
        **page,
    )
