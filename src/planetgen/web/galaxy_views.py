# planetgen/web/galaxy_views.py

"""
The Galaxy Map page (`/galaxy`, was `galaxy.py`) and the endpoint its script fetches as the camera moves
(`/galaxy/nodes`: the sector-block nodes of the map's cut).

The page renders the 3D map panel (`planetgen/web/maps/galaxymap3d.py`, drawn by
`static/galaxymap3d.js`) with the zoomed-all-the-way-out starting view's
nodes embedded, plus a plain-text table under it: a per-Quadrant summary,
or with `?quadrant=I|II|III|IV` that Quadrant's placed sectors nearest
the core first (paged with `?page=N`).

Nodes go through `planetgen/web/lib/tilecache.py`'s disk cache in this (the WSGI)
process: only nodes not already cached reach the API (in-process), and
whatever the API returns is written to the cache. The browser keeps its
own copy in IndexedDB, so a visitor's cached nodes survive a reload.
"""


from flask import Response, current_app, jsonify, make_response, redirect, request, url_for

from planetgen.web.lib import apiclient, galaxystatus
from planetgen.web.lib.fmt import format_distance_ly
from planetgen.web.maps.galaxymap import QUADRANT_LABELS, sector_quadrant, sector_zone, zone_bounds_ly
from planetgen.web.maps.galaxymap3d import render_galaxy_map3d_panel
from planetgen.web.maps.starmap import uncharted_scene_data
from planetgen.web.warmup import opening_request
from planetgen.web.lib.datatable import Column, Facet, Table, in_memory, plain
from planetgen.util import log
from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY
from planetgen.physics.units import ly_to_pc, pc_to_ly
from planetgen.web.lib import tilepack
from planetgen.galaxy.window import DEFAULT_WINDOW
from planetgen.web.lib import tileprefetch
from planetgen.web.lib import tilecache
from planetgen.web.lib.tilecache import (
    TileRequestError, fetch_near, fetch_nodes, fetch_nodes_wire, fetch_stage, fetch_window)

from planetgen.api.limiter import page_limit

from planetgen.admin import maintenance

from . import bp, glance, tables
from .helpers import crumb, current_admin, db_name, generate_target, page_url, render_page, trusted_html

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


def uncharted_scene_url_template():
    """
    The URL of an uncharted sector's scene JSON (`galaxy_uncharted_scene`)
    with `{shell}`, `{band}` and `{wedge}` where the address goes, for the
    Galaxy Map to open a cell nothing was generated in (MAP.162).
    """
    return url_for("web.galaxy_uncharted_scene", shell=_ID_PLACEHOLDER, band=_ID_PLACEHOLDER + 1,
                   wedge=_ID_PLACEHOLDER + 2).replace(
        str(_ID_PLACEHOLDER + 1), "{band}").replace(str(_ID_PLACEHOLDER + 2), "{wedge}").replace(
        str(_ID_PLACEHOLDER), "{shell}")


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


def uncharted_phenomenon_url_template():
    """
    The URL of an uncharted phenomenon's page with `{id}` where the scatter row id goes, for
    `static/galaxymap3d.js`'s "View phenomenon" link on a scattered point (MAP.211).
    """
    return page_url("uncharted_phenomenon", scatter_id=_ID_PLACEHOLDER).replace(str(_ID_PLACEHOLDER), "{id}")


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
        shells = [s["shell_index"] for s in members if s["shell_index"] is not None]
        if shells:
            _inner, outer_ly = zone_bounds_ly(sector_zone(max(shells)))
            extent = f"out to ~{format_distance_ly(outer_ly)}"
        rows.append({
            "label": label,
            "url": page_url("galaxy", quadrant=label, _anchor="galaxy-table"),
            "sector_count": len(members),
            "system_count": sum(s["system_count"] or 0 for s in members),
            "extent": extent,
        })
    return rows


def _zone_label(sector):
    return f"Zone {sector_zone(sector['shell_index'])}" if sector["shell_index"] is not None else None


def _quadrant_cells(sector):
    radius_ly = pc_to_ly(sector["galactic_radius_pc"])
    return [
        {"text": sector["name"], "href": page_url("sector", sector_id=sector["id"])},
        {"text": _zone_label(sector) or "–"},
        {"text": plain(format_distance_ly(radius_ly))},
        {"text": str(sector["system_count"] or 0)},
    ]


def _quadrant_load(state, limit, offset, want_facets):
    """The placed sectors of the Quadrant in `?quadrant=` (the table route's `quadrant`)."""
    quadrant = _parse_quadrant(request.args.get("quadrant"))
    sectors = apiclient.get_galaxy_sectors(db_name()) if quadrant else []
    members = [s for s in sectors if sector_quadrant(s["x"], s["y"]) == quadrant]
    return in_memory(
        members, state, limit, offset, want_facets,
        {"name": lambda s, descending: s["name"].casefold(),
         "zone": lambda s, descending: s["shell_index"],
         "radius": lambda s, descending: s["galactic_radius_pc"],
         "systems": lambda s, descending: s["system_count"] or 0},
        {"zone": _zone_label}, _quadrant_cells)


QUADRANT_TABLE = tables.register(Table(
    "galaxy-quadrant", "Sectors",
    [Column("name", "Name"), Column("zone", "Zone"), Column("radius", "From core"), Column("systems", "Systems")],
    _quadrant_load, facets=[Facet("zone", "Zone")], noun=("sector", "sectors"), default_sort="radius",
))


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
    # UX.112: the map is being rebuilt; one page says so. The cache warm-up's own visits (`WARMUP`) are let through: they are what
    # builds the map's caches while the update holds the site closed.
    if maintenance.is_on() and not current_app.config.get("WARMUP"):
        return redirect(url_for("web.maintenance_page"), code=302)
    db = db_name()
    running = galaxystatus.status(db)
    if running["active"]:   # PERF.149: nothing is drawn from a half-built galaxy
        return _generating_page(running)
    quadrant = _parse_quadrant(request.args.get("quadrant"))

    sectors = apiclient.get_galaxy_sectors(db)
    galaxy_shape = apiclient.get_galaxy_shape(db)
    edge_pc = galaxy_shape["edge_pc"] if galaxy_shape else ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    initial_view = fetch_nodes(db, opening_request(galaxy_shape))
    course = _course_from_args()
    made = _made_from_args()
    pick = _pick_from_args()
    map_html = render_galaxy_map3d_panel(
        db, galaxy_shape, edge_pc, initial_view,
        fetch_path=url_for("web.galaxy_nodes"),
        regions_path=url_for("web.galaxy_regions"),
        stage_path=url_for("web.galaxy_stage"),
        window_path=url_for("web.galaxy_window"),
        near_path=url_for("web.galaxy_near"),
        locate_path=url_for("web.galaxy_locate"),
        territory_path=url_for("web.galaxy_territories") if _has_territories(db) else None,
        course=course,
        made=made,
        sector_url=sector_url_template(),
        generate=generate_target(current_admin()),
        pick=pick,
        phenomenon_url=phenomenon_url_template(),
        system_url=system_url_template(),
        nav_url=page_url("nav"),
        nearby_url=page_url("nearby"),
        uncharted_scene_url=uncharted_scene_url_template(),
        uncharted_phenomenon_url=uncharted_phenomenon_url_template(),
        superbubbles=apiclient.get_superbubbles(db),
    )

    context = {"quadrant": quadrant, "placed_count": len(sectors), "map_html": trusted_html(map_html)}
    if quadrant:
        context["sector_table"] = tables.render(QUADRANT_TABLE, request.path, anchor="galaxy-table",
                                                keep=("quadrant",), quadrant=quadrant)
        title = f"Galaxy Map: Quadrant {quadrant}"
        breadcrumbs = [crumb("Galaxy", "galaxy"), crumb(f"Quadrant {quadrant}")]
    else:
        context["summary_rows"] = _quadrant_summary_rows(sectors)
        context["glance"] = glance.galaxy_rows(sectors, apiclient.get_containers(db))
        title = "Galaxy Map"
        breadcrumbs = [crumb("Galaxy")]

    return render_page(
        "galaxy.html",
        title=title,
        section="galaxy",
        breadcrumbs=breadcrumbs,
        description="An interactive 3D map of every sector placed in this charted galaxy.",
        **context,
    )


GENERATING_REFRESH_SECONDS = 30
"""int: How often the page that stands in for the Galaxy Map while a run changes the galaxy asks again."""


def _generating_page(running):
    """
    The page the Galaxy Map answers with while a run is changing the galaxy (PERF.149): what the banner says, the stage the
    run is at, a link to the Queue page for an admin, and a `503` (the map is back when the run ends) that is not kept.
    """
    admin = current_admin()
    html, status = render_page(
        "galaxy_generating.html", title="Galaxy Map", section="galaxy", breadcrumbs=[crumb("Galaxy")],
        description="The Galaxy Map is closed while the galaxy is generated.", status=503, status_info=running,
        status_text=running["text"], admin=admin is not None and not admin.get("must_change_credentials"),
        refresh=GENERATING_REFRESH_SECONDS)
    response = make_response(html, status)
    response.headers["Cache-Control"] = "no-store"
    return response


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


def _made_from_args():
    """
    The sectors a generate job made, from `?made=<since>,<until>` (Unix
    seconds, the job page's "Show on Galaxy Map", ADM.31), or `None`. The
    map highlights them and opens on the smallest stage holding them.
    """
    raw = (request.args.get("made") or "").strip()
    if not raw:
        return None
    parts = raw.split(",")
    try:
        since = float(parts[0])
        until = float(parts[1]) if len(parts) > 1 and parts[1] else None
    except ValueError:
        return None
    try:
        made = apiclient.get_galaxy_made(db_name(), since, until)
    except apiclient.ApiError as exc:
        log.exception(f"API error while finding the sectors a job made: {exc}")
        return None
    if not made["items"]:
        return {"total": 0, "points": [], "sectors": []}
    return {
        "total": made["total"],
        "points": [{"x": i["x"], "y": i["y"], "z": i["z"]} for i in made["items"]],
        "sectors": [[i["shell_index"], i["band_index"], i["wedge_index"]] for i in made["items"]],
    }


def _json_error(message, status):
    response = jsonify({"error": message})
    response.status_code = status
    return response


@bp.route("/galaxy/nodes")
@page_limit("galaxy_tiles")
def galaxy_nodes():
    """
    The nodes for the map's script (MAP.260): `?keys=<m.shell.band.wedge,...>` and optional `&stamp=<the browser cache's
    stamp>`. Returns `tilecache.fetch_nodes_wire`'s payload packed (`tilepack`) and sent as the stored bytes; a malformed
    request is a 400, an API failure a 502, both as `{"error": ...}` JSON.
    """
    node_keys = [key for key in (request.args.get("keys") or "").split(",") if key]
    known_stamp = request.args.get("stamp") or None
    db_name_asked = db_name()
    try:
        body = fetch_nodes_wire(db_name_asked, node_keys, known_stamp)
    except TileRequestError as exc:
        return _json_error(str(exc), 400)
    except apiclient.NotFoundError as exc:
        return _json_error(str(exc) or "Not found.", 404)
    except apiclient.ApiError as exc:
        log.exception(f"API error while fetching galaxy nodes: {exc}")
        return _json_error("The nodes could not be loaded. Please try again shortly.", 502)
    response = Response(body, mimetype=tilepack.MIMETYPE)
    response.headers["Cache-Control"] = "no-store"
    app = current_app._get_current_object()
    response.call_on_close(lambda: _prefetch_around(app, db_name_asked, node_keys))
    return response


@bp.route("/galaxy/regions")
@page_limit("galaxy_tiles")
def galaxy_regions():
    """
    The charted regions for the map's script (MAP.316): `{"regions": [...], "edge_pc"}` from the API's `charted-regions`
    (every region with its node key, extent, share charted and name, by share charted). Read on every request, never
    cached here: it is the list the map's regions and the Menu's list are made from and changes with the sectors. An API
    failure is a 502 as `{"error": ...}` JSON.
    """
    try:
        return jsonify(apiclient.get_galaxy_charted_regions(db_name()))
    except apiclient.ApiError as exc:
        log.exception(f"API error while fetching the charted regions: {exc}")
        return _json_error("The charted regions could not be loaded. Please try again shortly.", 502)


galaxy_regions.json_only = True  # not a page: tests/test_web_a11y.py skips it


def _prefetch_around(app, db, node_keys):
    """Queues the nodes around a request that has been answered (the response has gone, so the request's context
    with it: the cache reads its stamp through the API in process, which needs one); never an error for the visitor."""
    try:
        with app.test_request_context("/"):
            generation, changed = tilecache.take_rebuilds(db)
            if changed:
                tileprefetch.rebuild(db, generation, changed)      # what a change touched, built again at once
            tileprefetch.schedule(db, node_keys)
    except Exception as exc:  # noqa: BLE001 -- a background nicety, never a failure of the request
        log.debug(f"Node prefetch: not scheduled ({exc}).")


galaxy_nodes.json_only = True  # not a page: tests/test_web_a11y.py skips it


def _with_nebula_cover(payload, db, at, g=None):
    """
    `payload` (a stage or a window) with the nebula cover of its cells (MAP.174): `nebula_cover`, a list of `[shell,
    band, wedge, cover]` for every cell with a cover above zero, filled or not, and `stats.nebula_cover` (0 when
    none) on each child that holds a generated sector. Cells of a sector have none, and their payload carries an empty
    list. Read from the API on every request, never from the stage cache, because the aggregates change without a
    sector changing. When the API cannot answer, `nebula_cover` is `null` and the map draws without it.
    """
    try:
        found = apiclient.get_galaxy_nebula_cover(db, at, g)
    except apiclient.ApiError as exc:
        log.error(f"Galaxy Map: the nebula cover could not be read ({exc}).")
        return {**payload, "nebula_cover": None}
    cells = found.get("cells") or []
    covers = {(cell[0], cell[1], cell[2]): cell[3] for cell in cells}
    children = payload.get("children")
    if children and found.get("child_m", 1) > 1:
        children = [{**child, "stats": {**child["stats"], "nebula_cover": covers.get(
            (child["shell"], child["band"], child["wedge"]), 0.0)}} for child in children]
    return {**payload, **({"children": children} if children is not None else {}), "nebula_cover": cells}


@bp.route("/galaxy/stage")
@page_limit("galaxy_tiles")
def galaxy_stage():
    """
    JSON for the map's drill-down: `?at=<m.shell.band.wedge>` (or none, for
    the galaxy). Returns `tilecache.fetch_stage`'s payload; a malformed
    key is a 400, an API failure a 502, both as `{"error": ...}` JSON.
    """
    try:
        payload = fetch_stage(db_name(), request.args.get("at") or None)
        payload = _with_nebula_cover(payload, db_name(), request.args.get("at") or None)
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


@bp.route("/galaxy/window")
@page_limit("galaxy_tiles")
def galaxy_window():
    """
    JSON for the map's data window (MAP.151): `?at=<m.shell.band.wedge>` the
    middle cell (any level, a sector as `1.shell.band.wedge`) and `?g=` the
    cells across (3 to 11, odd; 9 by default). Returns
    `tilecache.fetch_window`'s payload; a malformed key or width is a 400, an
    API failure a 502, both as `{"error": ...}` JSON.
    """
    try:
        payload = fetch_window(db_name(), request.args.get("at") or "", request.args.get("g") or DEFAULT_WINDOW)
        payload = _with_nebula_cover(payload, db_name(), request.args.get("at") or "",
                                     request.args.get("g") or DEFAULT_WINDOW)
    except TileRequestError as exc:
        return _json_error(str(exc), 400)
    except apiclient.NotFoundError as exc:
        return _json_error(str(exc) or "Not found.", 404)
    except apiclient.ApiError as exc:
        log.exception(f"API error while fetching a galaxy window: {exc}")
        return _json_error("The map could not be loaded. Please try again shortly.", 502)
    response = jsonify(payload)
    response.headers["Cache-Control"] = "no-store"
    return response


galaxy_window.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/galaxy/near")
@page_limit("galaxy_tiles")
def galaxy_near():
    """
    JSON for the map's camera coming near sectors (MAP.335): `?x=&y=&z=` a point in parsecs (the eye) and `?r=` how far round
    it to look, parsecs. Returns `tilecache.fetch_near`'s payload, the generated sectors with the ids a scene is fetched by; a
    malformed number is a 400, an API failure a 502, both as `{"error": ...}` JSON.
    """
    try:
        payload = fetch_near(db_name(), request.args.get("x"), request.args.get("y"), request.args.get("z"), request.args.get("r") or 0)
    except TileRequestError as exc:
        return _json_error(str(exc), 400)
    except apiclient.NotFoundError as exc:
        return _json_error(str(exc) or "Not found.", 404)
    except apiclient.ApiError as exc:
        log.exception(f"API error while fetching the sectors near a point: {exc}")
        return _json_error("The map could not be loaded. Please try again shortly.", 502)
    response = jsonify(payload)
    response.headers["Cache-Control"] = "no-store"
    return response


galaxy_near.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/galaxy/locate")
@page_limit("search")
def galaxy_locate():
    """
    JSON for the map's address bar: `?q=<part of a name>`. Returns
    `{"matches": [...]}` (`db.query.galaxy_locate`); a 404 from the API is
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


@bp.route("/galaxy/uncharted/<int:shell>/<int:band>/<int:wedge>/scene")
@page_limit("galaxy_tiles")
def galaxy_uncharted_scene(shell, band, wedge):
    """
    The scene JSON of a sector cell nothing was generated in
    (`starmap.uncharted_scene_data`, MAP.162): the waiting bright stars and
    unbuilt scattered phenomena the scatters left there, for the Galaxy Map
    to open in place with an "Uncharted" mark. A cell that has been
    generated redirects to its sector's scene; a cell outside the grid is a
    404 and an API failure a 502, both as `{"error": ...}` JSON.
    """
    try:
        cell = apiclient.get_uncharted_sector(db_name(), shell, band, wedge)
    except apiclient.NotFoundError as exc:
        return _json_error(str(exc) or "Not found.", 404)
    except apiclient.ApiError as exc:
        if exc.status_code == 400:  # an address the grid doesn't have
            return _json_error("There is no such sector cell.", 404)
        log.exception(f"API error while fetching an uncharted sector's scene: {exc}")
        return _json_error("The sector could not be loaded. Please try again shortly.", 502)
    if cell["sector_id"] is not None:
        return redirect(url_for("web.sector_scene", sector_id=cell["sector_id"]), code=302)
    scene = uncharted_scene_data(
        (shell, band, wedge), cell["designation"], cell["center_pc"], cell["edge_pc"],
        {"stars": cell["stars"], "scattered": cell["scattered"]}, generate=generate_target(current_admin()),
        link_url=page_url)
    response = jsonify(scene)
    response.headers["Cache-Control"] = "no-store"
    return response


galaxy_uncharted_scene.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/galaxy/nebula/<uid:nebula_id>/shape")
@page_limit("galaxy_tiles")
def galaxy_nebula_shape(nebula_id):
    """
    JSON for the maps' nebula meshes (MAP.103): `?lod=low|full`. Returns
    `GET /api/nebulae/<id>/shape`'s payload; an unknown nebula is a 404
    and an API failure a 502, both as `{"error": ...}` JSON.
    """
    lod = request.args.get("lod", "low")
    try:
        payload = apiclient.get_nebula_shape(db_name(), nebula_id, lod)
    except apiclient.NotFoundError as exc:
        return _json_error(str(exc) or "Not found.", 404)
    except apiclient.ApiError as exc:
        log.exception(f"API error while fetching a nebula's shape: {exc}")
        return _json_error("The nebula could not be loaded. Please try again shortly.", 502)
    response = jsonify(payload)
    # A nebula's shape is drawn from its own properties, so it doesn't change.
    response.headers["Cache-Control"] = "private, max-age=3600"
    return response


galaxy_nebula_shape.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/galaxy/system/<uid:system_id>/nearest")
@page_limit("galaxy_tiles")
def galaxy_system_nearest(system_id):
    """
    JSON for the Galaxy Map's star system card (MAP.296): the system's nearest neighbours,
    `{"nearest": [{"id", "name", "distance_ly", "url"}]}` nearest first, from the same stored neighbours the system
    page's "Nearest:" line names (`db.query.system_detail`). An unknown system is a 404 and an API failure a 502, both
    as `{"error": ...}` JSON.
    """
    try:
        detail = apiclient.get_system(db_name(), system_id)
    except apiclient.NotFoundError as exc:
        return _json_error(str(exc) or "Not found.", 404)
    except apiclient.ApiError as exc:
        log.exception(f"API error while fetching a system's nearest neighbours: {exc}")
        return _json_error("The nearest systems could not be loaded. Please try again shortly.", 502)
    nearest = [
        {"id": near["id"], "name": near["name"], "distance_ly": near["distance_ly"], "url": page_url("system", system_id=near["id"])}
        for near in detail.get("nearest_neighbors") or []
    ]
    response = jsonify({"nearest": nearest})
    response.headers["Cache-Control"] = "private, max-age=300"
    return response


galaxy_system_nearest.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/galaxy/nebula/<uid:nebula_id>/surroundings")
@page_limit("galaxy_tiles")
def galaxy_nebula_surroundings(nebula_id):
    """
    JSON for a nebula page's 3D view (MAP.105): the brightest stars round
    it. An unknown nebula is a 404 and an API failure a 502, both as
    `{"error": ...}` JSON.
    """
    try:
        payload = apiclient.get_nebula_surroundings(db_name(), nebula_id)
    except apiclient.NotFoundError as exc:
        return _json_error(str(exc) or "Not found.", 404)
    except apiclient.ApiError as exc:
        log.exception(f"API error while fetching a nebula's surroundings: {exc}")
        return _json_error("The nebula's surroundings could not be loaded. Please try again shortly.", 502)
    response = jsonify(payload)
    response.headers["Cache-Control"] = "no-store"
    return response


galaxy_nebula_surroundings.json_only = True  # not a page: tests/test_web_a11y.py skips it


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
