# planetgen/web/galaxy_views.py

"""
The Galaxy Map page (`/galaxy`, was `galaxy.py`) and the JSON tile
endpoint its script fetches as the camera moves (`/galaxy/tiles`, was
`galaxy_tiles.py`).

The page renders the 3D map panel (`planetgen/web/maps/galaxymap3d.py`, drawn by
`static/galaxymap3d.js`) with the zoomed-all-the-way-out starting view's
tiles embedded, plus a plain-text table under it: a per-Quadrant summary,
or with `?quadrant=I|II|III|IV` that Quadrant's placed sectors nearest
the core first (paged with `?page=N`).

Tiles go through `planetgen/web/lib/tilecache.py`'s disk cache in this (the WSGI)
process: only tiles not already cached reach the API (in-process), and
whatever the API returns is written to the cache. The browser keeps its
own copy in `localStorage`, keyed as before, so a visitor's cached tiles
survive the move.
"""


from flask import jsonify, request, url_for

from planetgen.web.lib import apiclient
from planetgen.web.lib.fmt import format_distance_ly
from planetgen.web.maps.galaxymap import QUADRANT_LABELS, sector_quadrant, sector_zone, zone_bounds_ly
from planetgen.web.maps.galaxymap3d import initial_tile_request, render_galaxy_map3d_panel, view_radius_bounds
from planetgen.web.lib.pagination import page_slice, parse_page
from planetgen.util import log
from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY
from planetgen.physics.units import ly_to_pc, pc_to_ly
from planetgen.web.lib.tilecache import TileRequestError, fetch_stage, fetch_tiles

from planetgen.api.limiter import page_limit

from . import bp
from .helpers import crumb, current_admin, db_name, generate_target, page_url, pager, render_page, trusted_html

_ID_PLACEHOLDER = 987654321987
"""int: Stands in for a sector id while building the map's sector-link
template (`page_url` needs a real int for the `/sector/<int>` route), and
is then replaced by `{id}` for the script to fill in."""


def sector_url_template():
    """
    The URL of a sector page with `{id}` where the id goes, for
    `static/galaxymap3d.js`'s "View sector" link, built with `page_url`.
    """
    return page_url("sector", sector_id=_ID_PLACEHOLDER).replace(str(_ID_PLACEHOLDER), "{id}")


def system_url_template():
    """
    The URL of a star system page with `{id}` where the id goes, for
    `static/galaxymap3d.js`'s "View system" link on a filled bright star.
    """
    return page_url("system", system_id=_ID_PLACEHOLDER).replace(str(_ID_PLACEHOLDER), "{id}")


def phenomenon_url_template():
    """
    The URL of a phenomenon page with `{type}` and `{id}` where they go,
    for `static/galaxymap3d.js`'s "View phenomenon" link on a cloud.
    """
    url = page_url("phenomenon", phenomenon_type="nebula", phenomenon_id=_ID_PLACEHOLDER)
    return url.replace("/nebula/", "/{type}/").replace(str(_ID_PLACEHOLDER), "{id}")


def _parse_quadrant(raw):
    quadrant = (raw or "").strip().upper()
    return quadrant if quadrant in QUADRANT_LABELS else None


def _quadrant_summary_rows(sectors):
    """One row per Quadrant: placed sectors, total systems, Zone extent."""
    by_quadrant = {label: [] for label in QUADRANT_LABELS}
    for sector in sectors:
        by_quadrant[sector_quadrant(sector["x"], sector["y"])].append(sector)
    rows = []
    for label in QUADRANT_LABELS:
        members = by_quadrant[label]
        extent = None
        rings = [s["ring_index"] for s in members if s["ring_index"] is not None]
        if rings:
            _inner, outer_ly = zone_bounds_ly(sector_zone(max(rings)))
            extent = f"out to ~{format_distance_ly(outer_ly)}"
        rows.append({
            "label": label,
            "url": page_url("galaxy", quadrant=label, _anchor="galaxy-table"),
            "sector_count": len(members),
            "system_count": sum(s["system_count"] or 0 for s in members),
            "extent": extent,
        })
    return rows


def _quadrant_sector_rows(sectors, quadrant, page):
    """One page of the Quadrant's placed sectors, nearest the core first.
    Returns `(rows, page, total)`."""
    members = [s for s in sectors if sector_quadrant(s["x"], s["y"]) == quadrant]
    members.sort(key=lambda s: s["galactic_radius_pc"])
    page_members, page = page_slice(members, page)
    rows = [{
        "name": sector["name"],
        "url": page_url("sector", sector_id=sector["id"]),
        "zone": sector_zone(sector["ring_index"]) if sector["ring_index"] is not None else None,
        "distance_ly": pc_to_ly(sector["galactic_radius_pc"]),
        "distance": format_distance_ly(pc_to_ly(sector["galactic_radius_pc"])),
        "system_count": sector["system_count"] or 0,
    } for sector in page_members]
    return rows, page, len(members)


def _pick_from_args():
    """
    The NAV page's "Pick on Galaxy Map" (`?pick=from&to=...` or
    `?pick=to&from=...`, design doc section 9): the Sector Map's own pick
    mode (`sector_page._pick_mode`), whose `query` a sector link carries
    so the pick continues on that sector's page. `None` outside pick
    mode.
    """
    from .sector_page import _pick_mode  # a page module like this one; imported where used

    return _pick_mode(request.args)


def _has_territories(db):
    """Whether population has made any polities yet. Without one the map
    leaves out its Territories button, as the site does every population
    page; an API failure counts as none."""
    try:
        return apiclient.get_polities(db, limit=1)["total"] > 0
    except (apiclient.ApiError, apiclient.NotFoundError):
        return False


@bp.route("/galaxy")
@page_limit("galaxy")
def galaxy():
    """The 3D Galaxy Map plus the Quadrant table (`?quadrant=`, `?page=`)."""
    db = db_name()
    quadrant = _parse_quadrant(request.args.get("quadrant"))

    sectors = apiclient.get_galaxy_sectors(db)
    galaxy_shape = apiclient.get_galaxy_shape(db)
    edge_pc = galaxy_shape["edge_pc"] if galaxy_shape else ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    _min_radius, max_radius = view_radius_bounds(edge_pc, galaxy_shape)
    initial_view = fetch_tiles(db, initial_tile_request(max_radius))
    course = _course_from_args()
    pick = _pick_from_args()
    map_html = render_galaxy_map3d_panel(
        db, galaxy_shape, edge_pc, initial_view,
        fetch_path=url_for("web.galaxy_tiles"),
        stage_path=url_for("web.galaxy_stage"),
        locate_path=url_for("web.galaxy_locate"),
        territory_path=url_for("web.galaxy_territories") if _has_territories(db) else None,
        course=course,
        sector_url=sector_url_template() + (pick["query"] if pick else ""),
        generate=generate_target(current_admin()),
        pick=pick,
        phenomenon_url=phenomenon_url_template(),
        system_url=system_url_template(),
        nav_url=page_url("nav"),
    )

    context = {"quadrant": quadrant, "placed_count": len(sectors), "map_html": trusted_html(map_html)}
    if quadrant:
        rows, page, total = _quadrant_sector_rows(sectors, quadrant, parse_page(request.args.get("page")))
        context.update(
            sector_rows=rows,
            pager=pager("page", page, total, anchor="galaxy-table", label="Sector pages",
                        keep={"quadrant": quadrant}),
        )
        title = f"Galaxy Map: Quadrant {quadrant}"
        breadcrumbs = [crumb("Galaxy", "galaxy"), crumb(f"Quadrant {quadrant}")]
    else:
        context["summary_rows"] = _quadrant_summary_rows(sectors)
        title = "Galaxy Map"
        breadcrumbs = [crumb("Galaxy")]

    return render_page(
        "galaxy.html",
        title=title,
        section="galaxy",
        breadcrumbs=breadcrumbs,
        description="An interactive 3D map of every sector placed in this generated galaxy.",
        **context,
    )


def _course_from_args():
    """
    The NAV course the map should draw, from `?course=<from>,<to>` (the
    NAV result's "Show on Galaxy Map" link, the drill-down's section
    9.4), or `None`. A course that can't be drawn is simply not drawn;
    an endpoint that doesn't exist is the usual 404.
    """
    raw = (request.args.get("course") or "").strip()
    if not raw or "," not in raw:
        return None
    from_raw, to_raw = (part.strip() for part in raw.split(",", 1))
    if not from_raw or not to_raw:
        return None
    from .nav_page import galaxy_course  # imported here: both modules register routes on `bp`

    try:
        return galaxy_course(from_raw, to_raw)
    except apiclient.ApiError as exc:
        if isinstance(exc, apiclient.NotFoundError):
            raise
        log.exception(f"API error while plotting a course for the galaxy map: {exc}")
        return None


def _json_error(message, status):
    response = jsonify({"error": message})
    response.status_code = status
    return response


@bp.route("/galaxy/tiles")
@page_limit("galaxy_tiles")
def galaxy_tiles():
    """
    JSON for the map's script: `?tiles=<level/ix/iy/iz,...>` and optional
    `&stamp=<the browser cache's stamp>`. Returns
    `tilecache.fetch_tiles`' payload; a malformed request is a 400, an
    API failure a 502, both as `{"error": ...}` JSON.
    """
    tile_keys = [key for key in (request.args.get("tiles") or "").split(",") if key]
    known_stamp = request.args.get("stamp") or None
    try:
        payload = fetch_tiles(db_name(), tile_keys, known_stamp)
    except TileRequestError as exc:
        return _json_error(str(exc), 400)
    except apiclient.NotFoundError as exc:
        return _json_error(str(exc) or "Not found.", 404)
    except apiclient.ApiError as exc:
        log.exception(f"API error while fetching galaxy tiles: {exc}")
        return _json_error("The tiles could not be loaded. Please try again shortly.", 502)
    response = jsonify(payload)
    # Freshness is the tile cache's job (stamps); never let a shared
    # cache hand one visitor's stale answer to another.
    response.headers["Cache-Control"] = "no-store"
    return response


galaxy_tiles.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/galaxy/stage")
@page_limit("galaxy_tiles")
def galaxy_stage():
    """
    JSON for the map's drill-down: `?at=<m.ring.wedge.slab>` (or none, for
    the galaxy). Returns `tilecache.fetch_stage`'s payload; a malformed
    key is a 400, an API failure a 502, both as `{"error": ...}` JSON.
    """
    try:
        payload = fetch_stage(db_name(), request.args.get("at") or None)
    except TileRequestError as exc:
        return _json_error(str(exc), 400)
    except apiclient.NotFoundError as exc:
        return _json_error(str(exc) or "Not found.", 404)
    except apiclient.ApiError as exc:
        log.exception(f"API error while fetching a galaxy stage: {exc}")
        return _json_error("The map could not be loaded. Please try again shortly.", 502)
    response = jsonify(payload)
    response.headers["Cache-Control"] = "no-store"
    return response


galaxy_stage.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/galaxy/locate")
@page_limit("search")
def galaxy_locate():
    """
    JSON for the map's address bar: `?q=<part of a name>`. Returns
    `{"matches": [...]}` (`queryDb.galaxy_locate`); a 404 from the API is
    a 404 and any other API failure a 502, both as `{"error": ...}` JSON.
    """
    q = (request.args.get("q") or "").strip()[:200]
    if not q:
        return jsonify({"matches": []})
    try:
        matches = apiclient.get_galaxy_locate(db_name(), q)
    except apiclient.NotFoundError as exc:
        return _json_error(str(exc) or "Not found.", 404)
    except apiclient.ApiError as exc:
        log.exception(f"API error while looking up a name for the galaxy map: {exc}")
        return _json_error("The lookup failed. Please try again shortly.", 502)
    return jsonify({"matches": matches})


galaxy_locate.json_only = True  # not a page: tests/test_web_a11y.py skips it


TERRITORY_POLITY_LIMIT = 200
"""int: Most polities the territory overlay names. A galaxy has one
polity per spacefaring species, so this is far past any real count; the
rest are still drawn, just without a name in the legend."""


@bp.route("/galaxy/territories")
@page_limit("galaxy_tiles")
def galaxy_territories():
    """
    JSON for the map's territory overlay: `GET /api/territories`' owned
    systems and capitals, with each polity's name, color and system count
    folded in from `GET /api/polities` (the territory endpoint carries
    only ids). An API failure is a 502 as `{"error": ...}` JSON.
    """
    db = db_name()
    try:
        territories = apiclient.get_territories(db)
        polities = apiclient.get_polities(db, limit=TERRITORY_POLITY_LIMIT)["items"]
    except apiclient.ApiError as exc:
        log.exception(f"API error while fetching territories: {exc}")
        return _json_error("The territories could not be loaded. Please try again shortly.", 502)
    named = {polity["id"]: polity for polity in polities}
    merged = []
    for polity in territories["polities"]:
        detail = named.get(polity["id"], {})
        merged.append({
            "id": polity["id"], "capital_pc": polity["capital_pc"], "reach_ly": polity["reach_ly"],
            "name": detail.get("name"), "color": detail.get("color"),
            "government": detail.get("government"), "system_count": detail.get("system_count"),
        })
    response = jsonify({"points": territories["points"], "polities": merged})
    response.headers["Cache-Control"] = "no-store"
    return response


galaxy_territories.json_only = True  # not a page: tests/test_web_a11y.py skips it
