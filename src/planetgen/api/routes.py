# planetgen/api/routes.py

"""
JSON endpoints over the planetGen database.

The read endpoints each open a strictly-read-only connection
(`queryDb.open_readonly` -- same `file:...?mode=ro`-equivalent guarantee
the `planetgen.db.query` CLI already relies on, see that module's docstring for
how "read-only" is enforced by the configured account's grants). Listing
endpoints delegate straight to `planetgen.db.query`'s existing `list_sectors`/
`list_systems` functions rather than
re-implementing the same SQL a third time; the two detail endpoints go
through `planetgen.db.store.load_sector`/`load_star_system` for the full
nested object graph, serialized via each class's own `to_dict()` (Phase 1
serialization, `planetgen/util/serialization.py`).

Listing endpoints (`/sectors`, `/systems`) are paginated -- this project's
own roadmap (docs/TODO.md, Phase 4) plans galaxy-scale generation, so an
unbounded `SELECT *` here would eventually return an unbounded response.
`_paginate` centralizes parsing/validating `limit`/`offset` so both routes
apply the same defaults, cap, and error behavior.

The write endpoints (`POST`/`PATCH`/`DELETE` on `/sectors`/`/systems`,
near the bottom of this file) do real inserts/updates/deletes now, gated
behind admin authentication (`authz.require_admin`) and run against a
separate, write-capable database account -- see `docs/api.md`'s
"Write endpoints" section and `planetgen/admin/auth.py`.
"""

import math

import pymysql
from flask import Blueprint, current_app, g, jsonify, request, url_for

from planetgen.web.maps.systemscene import build_scene
from planetgen.generation import run_galaxy
from planetgen.queue import api_jobs
from planetgen.util import draw

from planetgen.db.query import (
    galaxy_layer_specs,
    facilities_for_system,
    facilities_in_sector,
    facility_detail,
    NO_SECTOR,
    bright_star_scatter_status,
    bright_stars_in_sector,
    uncharted_sector_contents,
    count_uncharted_systems,
    list_uncharted_systems,
    UNCHARTED_SYSTEM_SORTS,
    NavUnavailable,
    SEARCH_RESULT_LIMIT,
    SEARCH_RESULT_PANELS,
    SEARCH_TAG_FACETS,
    PHENOMENON_SORTS,
    SECTOR_SORTS,
    SYSTEM_SORTS,
    count_phenomena,
    count_sectors,
    is_capped,
    count_systems,
    MAX_TILES_PER_REQUEST,
    galaxy_changes,
    galaxy_density_shape,
    galaxy_placed_phenomena,
    galaxy_placed_sectors,
    sectors_made,
    galaxy_locate,
    galaxy_stage,
    galaxy_tiles,
    list_phenomena,
    phenomena_facets,
    sectors_facets,
    systems_facets,
    list_sectors,
    list_systems,
    nav_chart_plan,
    nav_course,
    open_readonly,
    phenomenon_detail as query_phenomenon_detail,
    resolve_object,
    nebula_shape as query_nebula_shape,
    nebula_surroundings as query_nebula_surroundings,
    search as run_search,
    sector_detail as query_sector_detail,
    system_detail as query_system_detail,
)
from planetgen.db import near as near_search
from planetgen.db import store
from planetgen.generation import bright_stars as brightStars, limits as generationLimits
from planetgen import tuning
from planetgen.galaxy import objectref as object_ref, version_check
from planetgen.db.store import (get_galaxy_bounds, get_galaxy_shape, get_sector_id_at, list_databases, resolve_database)
from planetgen.util.settings import get_settings
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.geometry import describe_sector_cell, sector_address_at
from planetgen.generation.system import StarSystem
from planetgen.db.render import FORMATS as SYSTEM_TEXT_FORMATS
from planetgen.db.render import render_system_sections, render_system_text
from planetgen.physics.units import ly_to_milliparsecs, ly_to_pc, pc_to_ly
from planetgen.util.format import format_distance_ly
from planetgen.wiki import WikiClient, WikiClientAuthError, WikiClientPageExistsError, WikiClientRequestError

from . import ids
from .authz import audit, require_admin
from .common import ApiError, request_json, require_json_body
from .version import API_VERSION
from .schemas import (
    BodyRename, FacilityCreate, NeighborhoodRequest, SectorCreate, SectorUpdate, SystemCreate, SystemPatch,
    SystemRename, WikiUpload, given, parse_body,
)
from .limiter import limiter, page_limit

bp = Blueprint("api", __name__, url_prefix="/api")

DEFAULT_PAGE_LIMIT = 100
MAX_PAGE_LIMIT = 500
MAX_PAGE_OFFSET = 2 ** 63 - 1
"""int: The largest `offset` MySQL's LIMIT/OFFSET takes as a plain
integer; past it the query itself fails, so it's a 400 instead."""
MAX_NEIGHBORHOOD_RADIUS_LY = generationLimits.MAX_GENERATE_RADIUS_LY
"""float: An upper bound on generate-neighborhood's `radius_ly` (about
652 ly, the same 200 pc cap the Generate page and `planetgen` use), so
no request can ask for an unbounded run."""

WRITE_RATE_LIMIT = "10 per minute"
"""str: Applied to every write route, on top of the app-wide default
limits (`config.Config.RATELIMIT_DEFAULT`) -- both apply together
(Flask-Limiter doesn't replace the default with a route's own `@limit`
unless told to), so a mutating endpoint is always throttled at least
this tightly regardless of how the global default is configured. Tighter
than any read endpoint's effective limit is today, on the general
principle that a request able to change data deserves more headroom
against abuse than one that only reads it."""


def _resolve_requested_db_config():
    """
    Resolves the connection config for this request: `current_app.config["MYSQL_CONFIG"]`
    (the process-wide default -- what every route used exclusively before
    this API learned about `?db=`) unless the request supplies a `db`
    query parameter, in which case it's validated against every schema on
    that same server whose name matches the configured prefix (see
    `planetgen.db.store.list_databases`/`resolve_database`) -- the same
    validation `html/lib/dbutil.resolve_db_name` has always applied for
    the CGI browser's own `?db=` picker, now shared by this API so a
    request can't select a schema this deployment never meant to expose.

    Returns:
        MySQLConfig: Ready to pass to `open_readonly`.

    Raises:
        ApiError: 404, if `db` is given but doesn't match a listed schema.
    """
    return _resolve_database_param(current_app.config["MYSQL_CONFIG"])


def _resolve_database_param(base_config):
    """`base_config`, or the `?db=` schema on its server (404 when it
    isn't one `list_databases` offers; 503, detail logged, when the
    server can't be asked)."""
    requested = request.args.get("db")
    if requested is None:
        return base_config
    try:
        return resolve_database(base_config, requested)
    except ValueError as exc:
        raise ApiError(str(exc), status_code=404)
    except pymysql.MySQLError as exc:
        _log_database_error(exc)
        raise ApiError(DATABASE_UNAVAILABLE, status_code=503)


def _statement_timeout_s():
    """`config.json`'s `mysql.statement_timeout_seconds` (PERF.17), or
    `None` when it's 0 or unset."""
    seconds = current_app.config.get("STATEMENT_TIMEOUT_S")
    if seconds is None:
        seconds = get_settings().mysql.statement_timeout_seconds
    return float(seconds) if seconds and float(seconds) > 0 else None


def get_db():
    """
    Returns the request-scoped read-only connection, opening one on first
    use (against `_resolve_requested_db_config`'s result -- the
    process-wide default database, or `?db=` when given). Reused for the
    lifetime of the request instead of one connection per query, then
    closed by `close_db` in the app's teardown handler.

    `queryDb.open_readonly` raises a bare `SystemExit` when the
    configured MySQL server is unreachable -- correct for its own CLI
    callers (a clean process exit with a message), but `SystemExit` is a
    `BaseException`, not an `Exception`: uncaught, it would propagate
    straight through Flask's request dispatch (and the WSGI worker
    handling it) instead of being turned into any HTTP response at all,
    from *every* route that calls this, not just `/health`'s own explicit
    check. Converted here into an `ApiError` (503 -- the same "database
    unreachable" status `/health` already wants to report) so it's just
    an ordinary exception from this point on: the app's registered
    `ApiError` handler turns it into the usual JSON error response for
    every other route, and `/health`'s own `except Exception` catches it
    directly.

    The response only says `DATABASE_UNAVAILABLE`: `open_readonly`'s
    message names the MySQL user, host and port, which are no business of
    a public client. The detail goes to the app's log.
    """
    if "db" not in g:
        try:
            g.db = open_readonly(_resolve_requested_db_config(), statement_timeout_s=_statement_timeout_s())
        except SystemExit as exc:
            _log_database_error(exc)
            raise ApiError(DATABASE_UNAVAILABLE, status_code=503)
    return g.db


DATABASE_UNAVAILABLE = "database unavailable"
"""str: The whole error message a client sees when the database can't be
opened or queried (a 503); the real reason is only logged."""


def _log_database_error(exc):
    """Logs why the database couldn't be used, for the operator only."""
    detail = str(exc) or type(exc).__name__
    # app.logger reaches Apache's error log and, when it's on, the debug
    # log (which listens on the root logger).
    current_app.logger.error("Database unavailable on %s %s: %s", request.method, request.path, detail)


def close_db(exception=None):
    db = g.pop("db", None)
    if db is not None:
        db.close()


def _parse_sort(query_args, allowed):
    """
    Parses the `sort`/`order` query parameters a sortable listing takes.

    Args:
        query_args (werkzeug.datastructures.MultiDict): `request.args`.
        allowed (Iterable[str]): The sort keys the listing offers; the first
            is the default.

    Returns:
        tuple[str, bool]: `(sort key, descending)`.

    Raises:
        ApiError: If `sort` is not one of `allowed` or `order` is not
            `asc`/`desc`.
    """
    allowed = list(allowed)
    sort = query_args.get("sort") or allowed[0]
    if sort not in allowed:
        raise ApiError(f"sort must be one of: {', '.join(allowed)}")
    order = (query_args.get("order") or "asc").lower()
    if order not in ("asc", "desc"):
        raise ApiError("order must be 'asc' or 'desc'")
    return sort, order == "desc"


def _parse_yes_no(raw, name):
    """`True` for "yes", `False` for "no", `None` when absent; anything else is an `ApiError`."""
    if raw is None or raw == "":
        return None
    if raw not in ("yes", "no"):
        raise ApiError(f"{name} must be 'yes' or 'no'")
    return raw == "yes"


def _paginate(query_args):
    """
    Parses and validates the `limit`/`offset` query parameters shared by
    every listing endpoint.

    `limit` defaults to `DEFAULT_PAGE_LIMIT` and is capped at
    `MAX_PAGE_LIMIT` (silently clamped, not rejected -- a client asking for
    "too much" isn't an error, just more than this API will hand back in
    one response). `offset` defaults to 0. Both must be non-negative
    integers when given at all.

    Args:
        query_args (werkzeug.datastructures.MultiDict): `request.args`.

    Returns:
        tuple[int, int]: `(limit, offset)`.

    Raises:
        ApiError: If `limit`/`offset` is present but not a non-negative
            integer.
    """
    raw_limit = query_args.get("limit")
    raw_offset = query_args.get("offset")

    if raw_limit is None:
        limit = DEFAULT_PAGE_LIMIT
    else:
        try:
            limit = int(raw_limit)
        except ValueError:
            raise ApiError(f"limit must be an integer, got {raw_limit!r}")
        if limit < 1:
            raise ApiError("limit must be at least 1")
        limit = min(limit, MAX_PAGE_LIMIT)

    if raw_offset is None:
        offset = 0
    else:
        try:
            offset = int(raw_offset)
        except ValueError:
            raise ApiError(f"offset must be an integer, got {raw_offset!r}")
        if offset < 0:
            raise ApiError("offset must be at least 0")
        if offset > MAX_PAGE_OFFSET:
            raise ApiError(f"offset must be at most {MAX_PAGE_OFFSET}")

    return limit, offset


def _parse_size_range(query_args, prefix):
    """
    Parses a `<prefix>_min_radius_km`/`<prefix>_max_radius_km` query
    parameter pair into a `(min_km, max_km)` tuple for `queryDb.search`'s
    `sizes` argument -- either bound may be omitted (`None`) for "no
    lower/upper limit"; giving neither returns `None` altogether (no size
    filter for this entity at all), the same "absent means inactive"
    convention every other search filter here already uses.

    Args:
        query_args (werkzeug.datastructures.MultiDict): `request.args`.
        prefix (str): `"star"`, `"planet"`, or `"moon"`.

    Returns:
        tuple[float or None, float or None] or None.

    Raises:
        ApiError: If either bound is present but not a valid non-negative
            number, or if both are given with min > max.
    """
    def _parse_bound(name):
        raw = query_args.get(name)
        if raw is None or raw.strip() == "":
            return None
        try:
            value = float(raw)
        except ValueError:
            raise ApiError(f"{name} must be a number, got {raw!r}")
        if not math.isfinite(value):
            raise ApiError(f"{name} must be a finite number, got {raw!r}")
        if value < 0:
            raise ApiError(f"{name} must be at least 0")
        return value

    min_km = _parse_bound(f"{prefix}_min_radius_km")
    max_km = _parse_bound(f"{prefix}_max_radius_km")
    if min_km is None and max_km is None:
        return None
    if min_km is not None and max_km is not None and min_km > max_km:
        raise ApiError(f"{prefix}_min_radius_km must not exceed {prefix}_max_radius_km")
    return (min_km, max_km)


@bp.route("/health")
@page_limit("health")
def health():
    """
    Liveness/readiness check for monitoring -- confirms the process is up
    and the configured database can actually be opened, not just that
    Flask is responding. Returns 503 (rather than letting the connection
    error propagate to a generic 500) so a monitor can tell "the API is
    running but its database is unreachable" apart from "the API itself is
    down".

    Also reports `schema_version`/`schema_current`: this API's read-only
    connection (`open_readonly`'s `ensure_schema=False`) never applies
    schema DDL itself, and even a read-write connection's `_ensure_schema`
    only runs `CREATE TABLE IF NOT EXISTS` -- a no-op against a table that
    already exists, so it never retroactively adds a migration's `ALTER
    TABLE` (e.g. schema v25's `idx_sectors_center`) to an existing
    database. Only `planetgen.cli.migrate` (run directly, or via `update.sh`/
    `install.sh`) actually advances an existing database's schema.
    Restarting this process alone -- a natural thing to try after pulling
    in a schema-fixing code change -- does *not* apply a pending
    migration, and previously wasn't surfaced anywhere: a stale schema
    silently kept e.g. the pre-v25 full-table-scan behind `GET
    /api/galaxy/view` that took the whole site down (`docs/apache-
    deployment.md`'s single `planetgen-api` process/thread pool serializes
    all API traffic, so one slow endpoint stalls every page). Surfaced
    here instead of only in `planetgen.cli.migrate`'s own output, so a live
    deployment that's fallen behind is visible without having to
    separately remember to go check.
    """
    try:
        get_db().execute("SELECT 1")
    except ApiError as exc:
        if exc.status_code != 503:
            raise  # e.g. an unknown `?db=` (404) -- the request's fault, not an outage
        # Already logged by get_db.
        return jsonify({"status": "error", "detail": DATABASE_UNAVAILABLE}), 503
    except Exception as exc:
        _log_database_error(exc)
        return jsonify({"status": "error", "detail": DATABASE_UNAVAILABLE}), 503

    # A separate try/except from the reachability check above: a database
    # that answers `SELECT 1` fine but has never had `schema.sql`/
    # `planetgen.cli.migrate` applied to it at all has no `schema_migrations` table
    # yet either -- one step further back than "some migrations pending"
    # (schema_row would simply come back empty for that), not an
    # unreachable database. Reported the same way a stale-but-present
    # schema_migrations row is below (200, schema_current False), not
    # folded into the 503 case above.
    try:
        schema_row = get_db().execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()
        schema_version = schema_row["version"] if schema_row else None
    except Exception:
        schema_version = None

    body = {
        "status": "ok",
        "api_version": API_VERSION,
        "schema_version": schema_version,
        "schema_current": schema_version == store.SCHEMA_VERSION,
    }
    if schema_version is not None and schema_version > store.SCHEMA_VERSION:
        body["detail"] = (
            f"Database schema is at v{schema_version}, newer than this code's v{store.SCHEMA_VERSION} -- "
            f"update planetGen (update.sh); don't run this older code's planetgen.cli.migrate against it."
        )
    elif schema_version != store.SCHEMA_VERSION:
        body["detail"] = (
            f"Database schema is at v{schema_version}, code expects v{store.SCHEMA_VERSION} -- "
            f"run planetgen.cli.migrate (or update.sh/install.sh) against this database."
            if schema_version is not None else
            f"Database schema has not been initialized yet (no schema_migrations table), code expects "
            f"v{store.SCHEMA_VERSION} -- run planetgen.cli.migrate (or update.sh/install.sh) against this database."
        )
    return jsonify(body)


@bp.route("/databases")
def databases():
    """
    Lists every MySQL schema on the configured server whose name matches
    this deployment's prefix (`planetgen.db.store.list_databases`), each
    with its size/last-modified stats plus a quick-glance sector/system
    count. Every other endpoint's own
    `?db=` selects among these same names (see `get_db`).
    """
    return jsonify({"items": list_databases(current_app.config["MYSQL_CONFIG"], with_counts=True)})


@bp.route("/sectors")
def sectors():
    """
    The sectors, nearest the galactic core first, paginated. The Sectors
    table (UX.41) also takes `sort` (`name`, `systems`, `density`,
    `position` or `distance`) with `order=asc|desc`, the repeatable filter
    `quadrant` (`I`-`IV`, `unplaced`) and `facets=1` for the Quadrant menu's
    option counts (`facets`). `total` counts the sectors that pass the filter (`total_capped`: it stopped at 10,000, PERF.77).
    """
    limit, offset = _paginate(request.args)
    sort = None
    descending = False
    if request.args.get("sort") or request.args.get("order"):
        sort, descending = _parse_sort(request.args, SECTOR_SORTS)
    quadrants = [v for v in request.args.getlist("quadrant") if v]
    db = get_db()
    total = count_sectors(db, quadrants=quadrants)
    body = {
        "items": list_sectors(db, limit=limit, offset=offset, sort=sort, descending=descending,
                              quadrants=quadrants),
        "total": total,
        "total_capped": is_capped(total),
        "limit": limit,
        "offset": offset,
    }
    if request.args.get("facets") == "1":
        body["facets"] = sectors_facets(db, quadrants=quadrants)
    return jsonify(body)


@bp.route("/sectors/<uid:sector_id>")
def sector_detail(sector_id):
    try:
        detail = query_sector_detail(get_db(), sector_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(detail)


def _parse_sector_id_filter(raw_sector_id):
    """
    Parses the `sector_id` query parameter `/api/systems` accepts: a
    sector's ID (one sector), the literal `"none"` (standalone systems --
    `queryDb.NO_SECTOR`, `html/browse.py`'s own table of systems generated
    with no sector), or absent (no filter).

    Raises:
        ApiError: If present and neither `"none"` nor a sector ID.
    """
    if raw_sector_id is None:
        return None
    if raw_sector_id.strip().lower() == "none":
        return NO_SECTOR
    try:
        ids.parse(ids.SECTOR, raw_sector_id)
    except ids.IdError:
        raise ApiError(f"sector_id must be a sector ID or 'none', got {raw_sector_id!r}")
    found = ids.row_id(get_db(), ids.SECTOR, raw_sector_id)
    return 0 if found is None else found  # a sector that is not there holds no systems


@bp.route("/systems")
def systems():
    """
    The systems, by name, paginated; `star_type` and `sector_id` filter them.
    The Systems tables (UX.41) also take `sort` (`name`, `sector`, `octant`
    or `binary`) with `order=asc|desc`, the filters `binary=yes|no`,
    `placement=sector|standalone` and `octant` (repeatable), and `facets=1`
    for the filter menus' option counts (`facets`: `placement`, `binary`,
    `octant`). `total` counts the systems that pass the filters. `after=<system ID>` starts the page after that system in the
    name sort instead of skipping `offset` rows (PERF.75: page 40,000 costs what page 1 does).
    """
    star_type = request.args.get("star_type")
    sector_id = _parse_sector_id_filter(request.args.get("sector_id"))
    sort, descending = _parse_sort(request.args, SYSTEM_SORTS)
    placement = request.args.get("placement")
    if placement not in (None, "", "sector", "standalone"):
        raise ApiError("placement must be 'sector' or 'standalone'")
    filters = {
        "binary": _parse_yes_no(request.args.get("binary"), "binary"),
        "octants": [v for v in request.args.getlist("octant") if v],
        "in_sector": None if not placement else placement == "sector",
    }

    limit, offset = _paginate(request.args)
    db = get_db()
    after = None
    if request.args.get("after"):
        try:
            ids.parse("system", request.args["after"])
        except ids.IdError:
            raise ApiError(f"after must be a system ID, got {request.args['after']!r}")
        after = ids.row_id(db, "system", request.args["after"])
    rows = list_systems(db, star_type_prefix=star_type, sector_id=sector_id, limit=limit, offset=offset,
                        sort=sort, descending=descending, after=after, **filters)
    total = count_systems(db, star_type_prefix=star_type, sector_id=sector_id, **filters)
    body = {
        "items": rows,
        "total": total,
        "total_capped": is_capped(total),
        "limit": limit,
        "offset": offset,
    }
    if request.args.get("facets") == "1":
        body["facets"] = systems_facets(db, sector_id=sector_id, **filters)
    return jsonify(body)


@bp.route("/systems/<uid:system_id>")
def system_detail(system_id):
    try:
        detail = query_system_detail(get_db(), system_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(detail)


@bp.route("/systems/<uid:system_id>/scene")
def system_scene(system_id):
    """
    `GET /api/systems/<id>/scene` -- the 3D system view's data (MAP.69):
    every star, planet, moon, belt and comet with its `ref`, radius, colour,
    orbit elements and position at `epoch` (see `web/maps/systemscene.py`).
    """
    try:
        scene = build_scene(get_db(), system_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(scene)


@bp.route("/systems/<uid:system_id>/text")
def system_text(system_id):
    """
    `GET /api/systems/<id>/text?format=wikitext|markdown` -- the system's
    full wiki page, rendered now from its database rows
    (`planetgen/db/render.py`; no page text is stored since
    schema v29). `format` defaults to `wikitext`.
    """
    fmt = request.args.get("format", "wikitext")
    if fmt not in SYSTEM_TEXT_FORMATS:
        raise ApiError(f"'format' must be one of: {', '.join(SYSTEM_TEXT_FORMATS)}")
    try:
        content = render_system_text(get_db(), system_id, fmt)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify({"id": system_id, "format": fmt, "content": content})


@bp.route("/systems/<uid:system_id>/sections")
def system_sections(system_id):
    """
    `GET /api/systems/<id>/sections` -- the same page as Markdown, split
    per star/planet/moon/belt/comet (keyed by row id) plus the system
    `overview`, for `html/system.py`'s expandable system list -- see
    `systemRender.render_system_sections`.
    """
    try:
        sections = render_system_sections(get_db(), system_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(sections)


@bp.route("/near")
def near():
    """
    Everything generated within a distance of a place (NAV.43). The place
    is `from` (an object reference such as `system:12`; a bare number is a
    system) or `point` (`x,y,z`, galaxy-frame parsecs); `distance` is in
    parsecs, up to 50. Optional: `kinds` (comma-separated), `limit`
    (default 50, at most 200) and `offset`.
    """
    args = request.args
    conn = get_db()
    if (args.get("from") is None) == (args.get("point") is None):
        raise ApiError("give exactly one of from (an object reference) or point (x,y,z in parsecs)")
    if args.get("distance") is None:
        raise ApiError("distance query parameter is required")
    try:
        if args.get("from") is not None:
            row_ref = _resolved_ref(args["from"], "from")
            try:
                place = near_search.place_from_reference(conn, row_ref)
            except near_search.NearError as exc:
                raise ApiError(str(exc).replace(row_ref, args["from"].strip()))
        else:
            place = near_search.place_from_point(args["point"].split(","))
        kinds = [kind for kind in args["kinds"].split(",") if kind] if args.get("kinds") else None
        result = near_search.objects_within(
            conn, place, args["distance"], kinds=kinds,
            limit=_int_arg(args, "limit", near_search.DEFAULT_LIMIT), offset=_int_arg(args, "offset", 0))
    except near_search.NearError as exc:
        raise ApiError(str(exc))
    return jsonify(result)


def _int_arg(query_args, name, default):
    raw = query_args.get(name)
    if raw is None:
        return default
    try:
        return int(raw)
    except ValueError:
        raise ApiError(f"{name} must be a whole number, got {raw!r}")


def _nav_ref_param(query_args, name):
    """
    A required `from`/`to` query parameter: an object reference
    (`planetgen.galaxy.objectref`; a bare number is a system).

    Raises:
        ApiError: If the parameter is missing or not a reference.
    """
    raw_value = query_args.get(name)
    if raw_value is None:
        raise ApiError(f"{name} query parameter is required")
    return _resolved_ref(raw_value, name)


def _resolved_ref(raw_value, name):
    """The row-id reference (`system:12`) of a public object reference (`system:FE81000A2B-0000005-000`).

    Raises:
        ApiError: 400 for text that is not a reference, 404 for one that names nothing.
    """
    try:
        resolved = ids.resolve_ref(get_db(), raw_value)
    except ids.IdError:
        raise ApiError(f"{name} must be an object reference such as system:FE81000A2B-0000005-000, got {raw_value!r}")
    if resolved is None:
        raise ApiError(f"{name} names nothing: {raw_value!r}", 404)
    return resolved


_HOST_KINDS = {"star": "star", "planet": "planet", "moon": "moon", "asteroid_belt": "belt",
               "asteroid_field": "asteroid_field", "space": "sector"}
"""dict: A facility's `host_type` -> the kind of object its `host_id` names."""


def _host_row(host_type, raw_value):
    """The row id of the host a facility request names by its printed ID.

    Raises:
        ApiError: 400 for a host type or ID that is not one, 404 for an ID that names nothing.
    """
    kind = _HOST_KINDS.get(host_type)
    if kind is None:
        raise ApiError(f"'host_type' must be one of: {', '.join(_HOST_KINDS)}")
    try:
        ids.parse(kind, raw_value)
        found = ids.row_id(get_db(), kind, raw_value)
    except ids.IdError:
        raise ApiError(f"host_id must be the ID of a {kind}, got {raw_value!r}")
    if found is None:
        raise ApiError(f"no {kind} with ID {raw_value}", 404)
    return found


NAV_MAX_STAY_MINUTES = 10_000_000.0
"""float: The longest stay per stop `?stay=` accepts (about 19 years)."""


def _nav_stay_param(query_args):
    """NAV.11: the optional `stay` query parameter, minutes spent at each
    stop of a route (default 0)."""
    raw_value = query_args.get("stay")
    if raw_value in (None, ""):
        return 0.0
    try:
        minutes = float(raw_value)
    except ValueError:
        raise ApiError(f"stay must be a number of minutes, got {raw_value!r}")
    if not 0.0 <= minutes <= NAV_MAX_STAY_MINUTES:
        raise ApiError(f"stay must be between 0 and {NAV_MAX_STAY_MINUTES:g} minutes, got {raw_value!r}")
    return minutes


@bp.route("/nav")
def nav():
    """
    Course, distance, and optimal route between two objects -- each any
    star system, body inside one, or standalone phenomenon, written as an
    object reference (`?from=planet:7&to=nebula:2`; a bare number is a
    system) -- see `queryDb.nav_course` for the legs and
    `queryDb.nav_between` for the availability rules, and docs/api.md for
    the response shape.
    """
    from_ref = _nav_ref_param(request.args, "from")
    to_ref = _nav_ref_param(request.args, "to")

    stay_minutes = _nav_stay_param(request.args)
    try:
        result = nav_course(get_db(), from_ref, to_ref, stay_minutes=stay_minutes)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    except NavUnavailable as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify({
        "scope": result["scope"],
        "origin": result["origin"],
        "destination": result["destination"],
        "direct": result["direct"]._asdict(),
        "warp_times": [leg._asdict() for leg in result["warp_times"]],
        "fold_times": [leg._asdict() for leg in result["fold_times"]],
        "origin_position": result["origin_position"],
        "destination_position": result["destination_position"],
        "route": _route_for_json(result["route"]),
        "legs": [{
            "kind": leg["kind"], "from": leg["from"], "to": leg["to"], "direct": leg["direct"]._asdict(),
            "warp_times": [t._asdict() for t in leg["warp_times"]],
            "fold_times": [t._asdict() for t in leg["fold_times"]],
        } for leg in result["legs"]],
        "note": result["note"],
    })


@bp.route("/nav/chart")
def nav_chart():
    """
    NAV.48: the uncharted sectors that block the course between two objects (`?from=...&to=...`, written
    as for `/api/nav`; `border=1` adds the uncharted cells next to them) -- see `queryDb.nav_chart_plan`.
    Answers `unknown_hops`, `cells` (`[ring, layer, slot]`, in the order the route enters them, inside the
    galaxy's outline), `count`, `outside_galaxy`, `route_distance_ly`, `bypass` (`{found, distance_ly,
    checked}`: whether a route through charted space exists, or `null` when no hop is unknown) and
    `confirm_over` (past this many cells the Generate page asks for a confirmation). Nothing is generated;
    the NAV chart page starts the job through the Generate page.
    """
    from_ref = _nav_ref_param(request.args, "from")
    to_ref = _nav_ref_param(request.args, "to")
    border = request.args.get("border", "") in ("1", "true", "yes", "on")
    try:
        plan = nav_chart_plan(get_db(), from_ref, to_ref, border=border)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    except NavUnavailable as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify({**plan, "cells": [list(cell) for cell in plan["cells"]], "count": len(plan["cells"]),
                    "confirm_over": tuning.NAV_CHART_CONFIRM_SECTORS})


def _route_for_json(route):
    """
    `route["positions"]` can have both plain-int keys (a system hop) and
    `queryDb._phenomenon_nav_key`-shaped string keys (a phenomenon
    origin/destination) since phenomenon endpoints were added -- Flask's
    `jsonify` sorts dict keys by default, and comparing an int key against
    a string key mid-sort raises `TypeError: '<' not supported between
    instances of 'int' and 'str'` (confirmed directly: a system-to-
    phenomenon `/api/nav` request 500'd here before this fix). JSON object
    keys are always strings regardless, so this stringifies every key at
    this JSON-serialization boundary specifically, rather than changing
    `nav_between`'s own return shape -- its direct Python callers (e.g.
    `test_navigation.py`) still index a system hop's own position by its
    real int id.
    """
    if route is None:
        return None
    return {
        "path": route["path"],
        "distance_ly": route["distance_ly"],
        "positions": {str(node_id): position for node_id, position in route["positions"].items()},
        "stop_places": route["stop_places"],
        "hops": route["hops"],
        "longest_hop_ly": route["longest_hop_ly"],
        "stops": route["stops"],
        "stay_minutes": route["stay_minutes"],
        "warp_times": [leg._asdict() for leg in route["warp_times"]],
        "fold_times": [leg._asdict() for leg in route["fold_times"]],
    }


@bp.route("/galaxy/sectors")
def galaxy_sectors():
    """
    Every galaxy-placed sector (`sectors.center_x/y/z_pc` not NULL), with
    its live system count -- the data the Galaxy Map (`/galaxy`) plots.
    Not paginated: bounded by how much of the galaxy has actually been
    generated so far (see `docs/TODO.md`'s Phase 4 lazy-generation design),
    not by the addressable galaxy's own astronomical scale.
    """
    return jsonify({"items": galaxy_placed_sectors(get_db())})


@bp.route("/galaxy/made")
def galaxy_made():
    """
    The galaxy-placed sectors created in a time window, for the Galaxy
    Map's "Show on Galaxy Map" after a generate job (ADM.31): `?since=`
    and optional `&until=` as Unix seconds. Returns `{"total", "items"}`
    (`queryDb.sectors_made`; at most `MADE_SECTOR_LIMIT` items).
    """
    try:
        since = float(request.args.get("since", ""))
        until = float(request.args["until"]) if request.args.get("until") else None
    except ValueError:
        return jsonify({"error": "since and until must be numbers (Unix seconds)."}), 400
    return jsonify(sectors_made(get_db(), since, until))


@bp.route("/galaxy/phenomena")
def galaxy_phenomena():
    """
    Every galaxy-placed standalone phenomenon (`center_x/y/z_pc` not
    NULL, in any `queryDb._PHENOMENON_TABLES` table -- see `schema.sql`'s
    "v18"/"v21"/"v29" header notes) -- the phenomenon counterpart to `/api/galaxy/sectors`, plotted
    as small dots on the same `/galaxy` Galaxy Map. Not paginated,
    for the same reason `/api/galaxy/sectors` isn't.
    """
    return jsonify({"items": galaxy_placed_phenomena(get_db())})


@bp.route("/galaxy/shape")
def galaxy_shape():
    """
    The galaxy's stored density-skeleton shape (`planetgen plan`'s own
    output, `queryDb.galaxy_density_shape`) -- the real spiral/disk/bulge
    model the Galaxy Map (`/galaxy`) shades its "expected density"
    cloud from, for whatever space hasn't actually been generated yet.
    `"shape"` is `null` when the skeleton has never been built (the map
    then falls back to its own generic illustrative gradient).
    `"bright_stars"` is `queryDb.bright_star_scatter_status`: whether the
    bright-star scatter has run, its threshold and seed, and the default
    threshold a plain plan uses. `"version_warning"` (DB.7) is a sentence
    when the galaxy holds sectors generated by another PlanetGen version
    than the running one, else `null`. `"layers"` (ADM.28) is
    `queryDb.galaxy_layer_specs`: the layer count, height, extent and how
    many sectors are charted, with a row per layer (a sample past 41).
    """
    conn = get_db()
    skeleton = get_galaxy_shape(conn)
    return jsonify({"shape": galaxy_density_shape(conn), "bright_stars": bright_star_scatter_status(conn),
                    "version_warning": version_check.galaxy_warning(conn),
                    "layers": galaxy_layer_specs(conn, skeleton.edge_pc) if skeleton else None})


@bp.route("/galaxy/cell")
def galaxy_cell():
    """
    One sector cell of the cylindrical grid, whether or not anything was
    ever generated there: `?ring=I&layer=J&slot=K` names it by address,
    `?x=&y=&z=` (parsecs, galaxy frame) by any point inside it. Returns
    `galaxyGeometry.describe_sector_cell` (center in Cartesian, cylindrical
    and spherical coordinates, bounds, 8 corners) plus `sector_id`, the
    generated sector there or `null`, and `in_galaxy`: whether the cell lies
    inside the planned galaxy's stored outline (`null` before any plan).
    Uses the stored skeleton's edge length, else the standard 4 pc.
    """
    def number(name, cast):
        value = request.args.get(name)
        if value is None:
            return None
        try:
            return cast(value)
        except ValueError:
            raise ApiError(f"{name} must be a number")

    conn = get_db()
    skeleton = get_galaxy_shape(conn)
    edge_pc = skeleton.edge_pc if skeleton else float(tuning.DEFAULT_SECTOR_EDGE_PC)
    ring, layer, slot = number("ring", int), number("layer", int), number("slot", int)
    x, y, z = number("x", float), number("y", float), number("z", float)
    # An address astronomically far out overflows the float geometry
    # (OverflowError) -- a bad request, not a server error.
    if None not in (ring, layer, slot):
        address = (ring, layer, slot)
    elif None not in (x, y, z) and all(math.isfinite(v) for v in (x, y, z)):
        try:
            address = sector_address_at((x, y, z), edge_pc)
        except OverflowError:
            raise ApiError("x, y and z are too far out")
    else:
        raise ApiError("give ring, layer and slot, or x, y and z")
    if address[0] < 0:
        raise ApiError("ring must be >= 0")
    try:
        cell = describe_sector_cell(*address, edge_pc)
    except OverflowError:
        raise ApiError("that cell is too far out")
    except ValueError as err:
        raise ApiError(str(err))
    cell["edge_pc"] = edge_pc
    cell["sector_id"] = get_sector_id_at(conn, *address)
    bounds = get_galaxy_bounds(conn)
    cell["in_galaxy"] = bounds.contains(address[0], address[1]) if bounds is not None else None
    return jsonify(cell)


@bp.route("/galaxy/tiles")
def galaxy_tiles_route():
    """
    The 3D Galaxy Map's cube tiles -- `tiles` is a comma-separated list of
    `level/ix/iy/iz` keys (at most `MAX_TILES_PER_REQUEST`). See
    `queryDb.galaxy_tiles` and `planetgen.galaxy.viewport`'s "Cube
    tiles" section. Each tile's work is bounded, so no request can scan an
    unbounded region (the removed `/galaxy/view` route could). Called by
    `/galaxy/tiles` (`planetgen/web/galaxy_views.py`), which caches every tile on disk and only forwards
    the ones it doesn't already have.
    """
    tile_keys = [key for key in (request.args.get("tiles") or "").split(",") if key]
    if len(tile_keys) > MAX_TILES_PER_REQUEST:
        raise ApiError(f"at most {MAX_TILES_PER_REQUEST} tiles per request")
    try:
        return jsonify(galaxy_tiles(get_db(), tile_keys))
    except ValueError as exc:
        raise ApiError(str(exc))


@bp.route("/galaxy/locate")
def galaxy_locate_route():
    """
    The Galaxy Map address bar's name lookup: `?q=<part of a name>`.
    Returns `{"matches": [...]}`, sectors and systems whose name contains
    it, each with its sector address -- see `queryDb.galaxy_locate`.
    Called by `/galaxy/locate` (`planetgen/web/galaxy_views.py`).
    """
    return jsonify({"matches": galaxy_locate(get_db(), request.args.get("q") or "")})


@bp.route("/galaxy/stage")
def galaxy_stage_route():
    """
    One Galaxy Map drill-down stage: `?at=m.ring.wedge.slab` (a block key,
    `planetgen.galaxy.drill`), or no `at` for the galaxy. Returns how
    many generated sectors each child block holds, and at a level-3 block
    the generated sectors themselves -- see `queryDb.galaxy_stage`.
    Called by `/galaxy/stage` (`planetgen/web/galaxy_views.py`), which caches
    each stage on disk.
    """
    try:
        return jsonify(galaxy_stage(get_db(), request.args.get("at") or None))
    except ValueError as exc:
        raise ApiError(str(exc))


@bp.route("/galaxy/stamp")
def galaxy_stamp_route():
    """
    `{"stamp": "<16 hex>", "state": "<token>"}` -- `stamp` changes whenever
    the galaxy's tile contents could change (see
    `queryDb.galaxy_content_stamp`), and `state` is what
    `/galaxy/changes?since=` takes to say which tiles did.
    """
    changes = galaxy_changes(get_db())
    return jsonify({"stamp": changes["stamp"], "state": changes["state"]})


@bp.route("/galaxy/changes")
def galaxy_changes_route():
    """
    `?since=<state>` -- which cube tiles changed since that state (see
    `queryDb.galaxy_changes`): `{"stamp", "state", "full", "tiles"}`.
    `planetgen/web/lib/tilecache.py` calls this about once a minute and deletes
    only the listed tiles, or all of them when `full`. A missing or
    unreadable `since` just answers `full`.
    """
    return jsonify(galaxy_changes(get_db(), request.args.get("since")))

@bp.route("/uncharted-systems")
def uncharted_systems():
    """
    `GET /api/uncharted-systems` -- the scattered stars no system was built
    around yet (UX.87), brightest first: `queryDb.list_uncharted_systems`'s
    rows (the star, its coordinates, its sector address, designation and
    position in the sector). Query parameters: `sort` (`luminosity` -- the
    default -- `temperature`, `type`, `sector` or `age`) with `order`
    (`asc` or `desc`; `desc` for luminosity when omitted, else `asc`), `limit`
    and `offset`. `total` counts every star waiting.
    """
    limit, offset = _paginate(request.args)
    sort, descending = _parse_sort(request.args, UNCHARTED_SYSTEM_SORTS)
    if "order" not in request.args:
        descending = sort == "luminosity"
    db = get_db()
    return jsonify({
        "items": list_uncharted_systems(db, limit=limit, offset=offset, sort=sort, descending=descending),
        "total": count_uncharted_systems(db), "limit": limit, "offset": offset,
    })


@bp.route("/uncharted-systems/<int:bright_star_id>/generate", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def generate_uncharted_system(bright_star_id):
    """
    `POST /api/uncharted-systems/<id>/generate` -- builds a system around
    one scattered star by itself (UX.87): the star's own data, a companion,
    planets and moons from its stored seed, saved as a standalone system and
    linked to the scattered star. Nothing else in its sector is generated
    (generate the whole sector instead to fill it). 404 for an unknown star;
    409 when a system already exists for it. Runs on the queue and answers
    `{"id": <system id>}` (201), or `202` with a job when it takes longer.
    """
    job_id, made = run_queued(generate_uncharted_system_job, _resolve_requested_write_db_config(), bright_star_id)
    if made is None:
        audit("system.generate_uncharted", target=f"job:{job_id}", detail=f"queued bright_star={bright_star_id}")
        return accepted(job_id)
    audit("system.generate_uncharted", target=f"system:{made['id']}", detail=f"bright_star={bright_star_id}")
    return jsonify(made), 201


def generate_uncharted_system_job(config, bright_star_id):
    """The queued body of `POST /api/uncharted-systems/<id>/generate`:
    builds the scattered star's system and links the two. Returns `{"id"}`."""
    conn = store.open_write(config)
    try:
        with conn:
            row = conn.execute("SELECT * FROM bright_stars WHERE id = ? FOR UPDATE", (bright_star_id,)).fetchone()
            if row is None:
                raise ApiError(f"no such uncharted star: {bright_star_id}", status_code=404)
            if row["star_system_id"] is not None:
                raise ApiError(f"star {bright_star_id} already has a system ({row['star_system_id']})",
                               status_code=409)
            system_config = SystemConfig()
            system_config.POPULATION = row["population"]
            x, y, z = (row[f"position_{axis}_mpc"] / brightStars.MPC_PER_PC for axis in "xyz")
            with draw.bound(row["seed"]):
                system = _generate_system(system_config, pc_to_ly(math.sqrt(x * x + y * y + z * z)),
                                          primary_star_params=brightStars.star_params(row))
            system_id = store.insert_star_system(conn, system, system_config)
            store.assign_uids(conn, system_ids=[system_id])
            store.mark_bright_star_filled(conn, bright_star_id, system_id)
    finally:
        conn.close()
    return {"id": system_id}


@bp.route("/phenomena")
def phenomena():
    """
    Every exotic phenomenon (nebula/asteroid field/black hole/neutron
    star/supernova remnant/rogue planet/interstellar comet), across every
    sector and regardless of galaxy placement -- the flat, paginated
    counterpart to `/api/galaxy/phenomena` (which only returns the
    galaxy-placed subset, for the Galaxy Map). `html/phenomena.py`'s own
    listing page.

    Query parameters: `sort` (`name` -- the default -- `type`,
    `descriptor`, `radius`, `sector` or `placed`) with `order` (`asc` or
    `desc`); the filters `type` and `descriptor` (each repeatable, any of),
    and `placed` (`yes` or `no`); and `facets=1` to add the option counts
    for the type and descriptor filter menus (`queryDb.phenomena_facets`).
    `total` counts the rows that pass the filters.
    """
    limit, offset = _paginate(request.args)
    sort, descending = _parse_sort(request.args, PHENOMENON_SORTS)
    filters = {
        "types": [v for v in request.args.getlist("type") if v],
        "descriptors": [v for v in request.args.getlist("descriptor") if v],
        "placed": _parse_yes_no(request.args.get("placed"), "placed"),
    }
    db = get_db()
    body = {
        "items": list_phenomena(db, limit=limit, offset=offset, sort=sort, descending=descending, **filters),
        "total": count_phenomena(db, **filters),
        "limit": limit,
        "offset": offset,
    }
    if request.args.get("facets") == "1":
        body["facets"] = phenomena_facets(db, **filters)
    return jsonify(body)


@bp.route("/phenomena/<phenomenon_type>/<uid:phenomenon_id>")
def phenomenon(phenomenon_type, phenomenon_id):
    """
    One phenomenon's full detail -- `html/phenomenon.py`'s info page.
    `phenomenon_type` is one of `queryDb._PHENOMENON_TYPE_TO_TABLE`'s keys
    (`nebula`/`asteroid_field`/`black_hole`/`neutron_star`/
    `supernova_remnant`/`rogue_planet`/`interstellar_comet`); anything
    else, or an id that doesn't exist under it, is a 404.
    """
    try:
        detail = query_phenomenon_detail(get_db(), phenomenon_type, phenomenon_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(detail)


@bp.route("/objects/<ref>")
def object_ref_route(ref):
    """
    Resolves one object reference (`<kind>:<id>`, NAV.7): name, parent
    chain, sibling references and positions in each frame -- see
    `queryDb.resolve_object` and docs/api.md.
    """
    try:
        resolved = ids.resolve_ref(get_db(), ref)
        if resolved is None:
            raise ValueError(f"nothing has the ID {ref}")
        kind, object_id = object_ref.parse(resolved)
        return jsonify(resolve_object(get_db(), kind, object_id))
    except ValueError as exc:  # includes ids.IdError
        return jsonify({"error": str(exc)}), 404


@bp.route("/nebulae/<uid:nebula_id>/shape")
def nebula_shape_route(nebula_id):
    """
    One nebula's shape as a triangle mesh (GEN.75): `?lod=low` (the default,
    a coarse mesh for the Galaxy Map) or `?lod=full`. Answers `{"id",
    "radius_ly", "center_pc" (or null for a nebula never placed), "lod",
    "vertices", "faces"}`: the vertices are `[x, y, z]` in units of the
    nebula's radius from its center (multiply by `radius_ly` for light-years),
    the faces triples of vertex indexes. A 404 for an unknown nebula or a
    `lod` other than `low` or `full`.
    """
    lod = request.args.get("lod", "low")
    if lod not in ("low", "full"):
        return jsonify({"error": f"lod must be 'low' or 'full', got {lod!r}"}), 404
    try:
        row, shape = query_nebula_shape(get_db(), nebula_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    vertices, faces = shape.mesh(lod)
    center = None if row["center_x_pc"] is None else [row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]]
    return jsonify({
        "id": nebula_id, "radius_ly": row["radius_ly"], "center_pc": center, "lod": lod,
        "vertices": [[round(float(v), 5) for v in vertex] for vertex in vertices],
        "faces": [[int(i) for i in face] for face in faces],
    })


@bp.route("/nebulae/<uid:nebula_id>/surroundings")
def nebula_surroundings_route(nebula_id):
    """
    The brightest stars round one nebula (MAP.105): `{"radius_pc",
    "half_width_pc", "stars": [{"x", "y", "z" (parsecs from the nebula's
    centre), "luminosity_sol", "temperature_k"}]}`, the most luminous first.
    A 404 for an unknown nebula.
    """
    try:
        return jsonify(query_nebula_surroundings(get_db(), nebula_id))
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404


@bp.route("/search")
def search():
    """
    Faceted search over the chosen database: click-to-filter attribute
    tags (object type; star spectral/luminosity class; planet/moon class,
    body type, supported life chemistry; asteroid belt density), a
    min/max size range per entity, plus a per-entity name search -- see
    `queryDb.search`'s own docstring for the full request/response shape.
    `html/search.py` is a thin renderer over this endpoint's JSON.

    Query parameters: `sector_q`/`system_q`/`star_q`/`planet_q`/`moon_q`
    (name search terms); every name in `queryDb.SEARCH_TAG_FACETS`
    (repeatable, e.g. `?class=M&class=K` for two active Class tags); and
    `star_min_radius_km`/`star_max_radius_km` (likewise `planet_`/
    `moon_`) for a size range, in km -- either bound may be omitted.

    Each result panel is paged independently: `limit` is the rows per
    panel (default `SEARCH_RESULT_LIMIT`, clamped to `MAX_PAGE_LIMIT`,
    validated like every other `limit`), and `sectors_offset`/
    `systems_offset`/`stars_offset`/`planets_offset`/`moons_offset`/
    `belts_offset` pick each panel's page. `panels` (comma-separated panel
    names) runs just those panels, the rest answering `null`.
    """
    args = request.args
    texts = {
        key: (args.get(key) or "").strip()
        for key in ("sector_q", "system_q", "star_q", "planet_q", "moon_q")
    }
    tags = {facet: {v.strip() for v in args.getlist(facet) if v.strip()} for facet in SEARCH_TAG_FACETS}
    tags["type"] &= {"star", "planet", "moon", "belt"}
    tags["body"] &= {"t", "g"}
    tags["moon_body"] &= {"t", "g"}
    sizes = {
        "star": _parse_size_range(args, "star"),
        "planet": _parse_size_range(args, "planet"),
        "moon": _parse_size_range(args, "moon"),
    }

    limit = SEARCH_RESULT_LIMIT
    if args.get("limit") is not None:
        limit, _offset = _paginate({"limit": args.get("limit")})
    offsets = {
        panel: _paginate({"offset": args.get(f"{panel}_offset")})[1]
        for panel in SEARCH_RESULT_PANELS
    }

    only = None
    if args.get("panels"):
        only = {panel for panel in args["panels"].split(",") if panel in SEARCH_RESULT_PANELS}
    return jsonify(run_search(get_db(), texts, tags, sizes=sizes, limit=limit, offsets=offsets, only=only))


@bp.route("/wiki-config")
def wiki_config():
    """`GET /api/wiki-config` -- `{"wikijs": bool, "mediawiki": bool}`,
    whether each backend has a usable `base_url` plus credentials
    configured (see `config.py`'s `_wiki_config`) and so is offered as an
    "Upload to Wiki" target at all. Read by the CGI browser
    (`html/system.py`/`html/sector.py`) to decide which backend choice(s)
    to show its upload form -- never leaks any of `WIKI_CONFIG`'s actual
    credential values, only these two booleans."""
    wiki = current_app.config["WIKI_CONFIG"]
    return jsonify({"wikijs": wiki["wikijs"]["configured"], "mediawiki": wiki["mediawiki"]["configured"]})


# ---------------------------------------------------------------------
# Write endpoints. Every one requires an authenticated admin
# (`authz.require_admin`) whose credentials aren't still the seeded
# default (`fresh=True` -- see `authz.py`/`adminAuth.py`), runs against
# `WRITE_MYSQL_CONFIG` (a separate, less-privileged-than-`CREATE`
# account -- see `config.py`) rather than the request-scoped read-only
# `get_db()` every read route above uses, and records one
# `admin_audit_log` row (`authz.audit`) after it actually succeeds, never
# speculatively. See docs/api.md's "Write endpoints" section for the
# request/response shapes.
# ---------------------------------------------------------------------

def _resolve_requested_write_db_config():
    """
    Same as `_resolve_requested_db_config` above, but resolved against
    `WRITE_MYSQL_CONFIG` (see `config.py`) instead of the read-only
    `MYSQL_CONFIG` every read route uses -- every write route's `?db=`
    selection (or lack of one, falling back to the write account's own
    configured default database) needs the write-capable account's
    credentials, not the `SELECT`-only one.
    """
    return _resolve_database_param(current_app.config["WRITE_MYSQL_CONFIG"])


def _write_conn():
    """Opens a write-capable connection against the content database
    `?db=` (or the configured default) selects -- `ensure_schema=False`,
    same reasoning as `queryDb.open_readonly`/`store.open_write`: the
    write-capable account has no `CREATE` grant (see `config.py`), so the
    schema must already exist."""
    return store.open_write(_resolve_requested_write_db_config())


@bp.route("/sectors", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def create_sector():
    """`POST /api/sectors` `{"name": str, "edge_ly": number > 0}` (both
    required) -- inserts a new, empty `sectors` row (no systems -- use
    `sectorGen.py`/`galaxyGen.py` to generate a populated one) and
    returns its id."""
    body = given(parse_body(SectorCreate, require_json_body()))

    conn = _write_conn()
    try:
        with conn:
            cur = conn.execute(
                "INSERT INTO sectors (name, edge_mpc) VALUES (?, ?)",
                (body["name"], ly_to_milliparsecs(body["edge_ly"])),
            )
            sector_id = cur.lastrowid
            store.set_sector_uid(conn, sector_id)
            printed = ids.printed(conn, ids.SECTOR, sector_id)
    finally:
        conn.close()

    audit("sector.create", target=f"sector:{printed}", detail=f"name={body['name']!r} edge_ly={body['edge_ly']}")
    return jsonify({"id": sector_id, "name": body["name"], "edge_ly": body["edge_ly"]}), 201


@bp.route("/sectors/<uid:sector_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def update_sector(sector_id):
    """`PATCH /api/sectors/<id>` -- any non-empty subset of `{"name": str,
    "edge_ly": number > 0, "wiki_url": str or null}`. `wiki_url` is the
    admin "manually set the wiki link" affordance (`html/admin.py`) --
    `null` clears it back to "no page yet" (see `schema.sql`'s "v22"
    header note); the same column is also written automatically by
    `POST /api/sectors/<id>/wiki` on a successful upload."""
    body = given(parse_body(SectorUpdate, require_json_body()))

    set_clauses, params = [], []
    if "name" in body:
        set_clauses.append("name = ?")
        params.append(body["name"])
    if "edge_ly" in body:
        set_clauses.append("edge_mpc = ?")
        params.append(ly_to_milliparsecs(body["edge_ly"]))
    if "wiki_url" in body:
        set_clauses.append("wiki_url = ?")
        params.append(body["wiki_url"])
    params.append(sector_id)

    conn = _write_conn()
    try:
        with conn:
            # Existence is checked explicitly rather than relying on the
            # UPDATE's own affected-row count: MySQL (and so pymysql)
            # reports 0 affected rows for a matched row whose new values
            # equal its current ones, which would otherwise misreport a
            # real, unchanged sector as "not found".
            if conn.execute("SELECT id FROM sectors WHERE id = ?", (sector_id,)).fetchone() is None:
                raise ApiError(f"no such sector: {sector_id}", status_code=404)
            conn.execute(f"UPDATE sectors SET {', '.join(set_clauses)} WHERE id = ?", tuple(params))
    finally:
        conn.close()

    audit("sector.update", target=f"sector:{sector_id}", detail=str(body))
    return jsonify({"status": "ok"})


@bp.route("/sectors/<uid:sector_id>", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def delete_sector(sector_id):
    """`DELETE /api/sectors/<id>` -- no request body. Its systems are
    detached, not deleted (`star_systems.sector_id` is `ON DELETE SET
    NULL`, per `schema.sql`), matching how a system can already exist
    with no sector at all."""
    conn = _write_conn()
    try:
        with conn:
            address = conn.execute("SELECT ring_index, layer_index, ring_slot_index FROM sectors WHERE id = ?",
                                   (sector_id,)).fetchone()
            deleted = conn.execute("DELETE FROM sectors WHERE id = ?", (sector_id,)).rowcount > 0
            if deleted:
                store.note_sector_deleted(conn)
            if address is not None:
                # GEN.44: the slot goes back to the bright-star level it had unfilled.
                store.forget_sector_fill(conn, (address["ring_index"], address["layer_index"],
                                              address["ring_slot_index"]))
    finally:
        conn.close()

    if not deleted:
        raise ApiError(f"no such sector: {sector_id}", status_code=404)
    audit("sector.delete", target=f"sector:{sector_id}")
    return jsonify({"status": "ok"})


@bp.route("/sectors/<uid:sector_id>/generate-neighborhood", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def generate_sector_neighborhood_route(sector_id):
    """`POST /api/sectors/<id>/generate-neighborhood` -- generates every
    not-yet-generated sector within `radius_ly` (optional JSON body
    field; defaults to `tuning.DEFAULT_GENERATE_RADIUS_PC`,
    12 pc, the sphere `planetgen galaxy`'s own random-start mode uses)
    of this already galaxy-placed sector -- see
    `run_galaxy.generate_sector_neighborhood`. Every new sector also gets
    the bright stars within 100 ly of it (GEN.23). **A large radius is
    genuinely large** (100 ly is ~2,000-3,000 candidate sector slots), so
    this can run for minutes to hours, not seconds. Runs synchronously like every
    other write route regardless (there's no background job queue in
    this project to hand it off to) -- `apiclient.py`'s own caller uses a
    much longer timeout than its other calls for exactly this reason, but
    a production deployment's own reverse-proxy/gateway timeout (Apache,
    etc.) may still need raising for this one route to ever complete over
    HTTP at all.

    PERF.3: `"estimate_only": true` returns the counts and the size and
    time `estimate` without writing anything; a run the database disk
    can't hold (over a quarter of it, or under 5 GB left) is refused with
    507 and nothing written."""
    body = request_json()
    if body is None:
        body = {}
    if not isinstance(body, dict):
        raise ApiError("request body must be a JSON object")
    request_body = parse_body(NeighborhoodRequest, body)
    radius_ly, estimate_only = request_body.radius_ly, request_body.estimate_only

    config = _resolve_requested_write_db_config()
    if not estimate_only:
        try:
            run_galaxy.require_math_check()
        except RuntimeError as exc:
            raise ApiError(str(exc), status_code=409)
    try:
        # A real run is estimated first, so an unknown sector (404), a
        # missing skeleton (409) and a database disk that can't hold it
        # (507) answer here, before anything is queued.
        result = run_galaxy.generate_sector_neighborhood(
            sector_id, radius_ly=radius_ly, config=config, estimate_only=True,
        )
    except ValueError as exc:
        raise ApiError(str(exc), status_code=404)
    except RuntimeError as exc:
        raise ApiError(str(exc), status_code=409)
    if estimate_only:
        return jsonify(result)
    if result["estimate"].get("refusal"):
        # PERF.3: the database disk can't hold it; nothing is queued.
        raise ApiError(result["estimate"]["refusal"], status_code=507)

    try:
        job_id = api_jobs.submit(api_jobs.generate_neighborhood, sector_id, radius_ly, config)
    except api_jobs.NoQueue as exc:
        raise ApiError(str(exc), status_code=503)
    audit(
        "sector.generate_neighborhood", target=f"sector:{sector_id}",
        detail=f"radius_ly={radius_ly!r} job={job_id}",
    )
    return accepted(job_id)


def accepted(job_id):
    """`202 Accepted` for work queued as job `job_id`: its state is at the
    `Location` header's `GET /api/jobs/<id>`."""
    url = f"/api/jobs/{job_id}"
    response = jsonify({"status": "accepted", "job_id": job_id, "status_url": url})
    response.status_code = 202
    response.headers["Location"] = url
    return response


def run_queued(function, *args):
    """
    Runs `function(*args)` on the queue (PERF.24) and waits a few seconds
    for it, so a quick edit still answers in the same response.

    Returns:
        tuple: `(job_id, result)`; `result` is `None` while the job is
            still running, and the route then answers `accepted(job_id)`.

    Raises:
        ApiError: The work refused (its own status: 404, 409 ...), failed
            (500), or no Redis server answers (503).
    """
    try:
        job_id = api_jobs.submit(function, *args)
        job = api_jobs.wait(job_id)
    except api_jobs.NoQueue as exc:
        raise ApiError(str(exc), status_code=503)
    if job is not None and job["state"] == "failed":
        raise ApiError(job["error"], status_code=job["error_status"] or 500)
    if job is not None and job["state"] == "succeeded":
        return job_id, job["result"]
    return job_id, None


@bp.route("/jobs/<job_id>", methods=["GET"])
@require_admin(scope="read")
def job_status_route(job_id):
    """`GET /api/jobs/<id>` -- where a queued API job stands: `state`
    (`queued`, `running`, `succeeded` or `failed`), its `result` once it
    succeeded, its `error` once it failed (PERF.24). 404 for an
    unknown or expired id; results are kept for a day. A finished job that
    generated sectors also gives `made_url`, the Galaxy Map fitted to them
    (ADM.31), else `null`."""
    try:
        job = api_jobs.status(job_id)
    except api_jobs.NoQueue as exc:
        raise ApiError(str(exc), status_code=503)
    if job is None:
        raise ApiError(f"no such job: {job_id}", status_code=404)
    made = job.pop("made")
    job["made_url"] = url_for("web.galaxy", made=f"{made['since']},{made['until']}") if made else None
    return jsonify(job)


def _sector_generation_context(conn, sector_id, system_config):
    """For a galaxy-placed sector: the system's distance from the galactic
    center in light-years, after giving `system_config` the sector's
    stellar population and (once a bright-star scatter ran) its dim-star
    cap (its own backfill level and mass, GEN.44, GEN.187), as a sector fill would
    (`brightStars.FillContext`). `None` for a
    sector outside the galaxy. 404 when the sector doesn't exist."""
    try:
        placement = store.get_sector_galaxy_position(conn, sector_id)
    except ValueError:
        raise ApiError(f"no such sector: {sector_id}", status_code=404)
    if placement is None:
        return None
    skeleton = get_galaxy_shape(conn)
    if skeleton is not None:
        address = (placement["ring_index"], placement["layer_index"], placement["ring_slot_index"])
        level = store.bright_star_fill_level(conn, *address)
        mass_limit = store.bright_star_fill_mass_limit(conn, *address)
        center = (placement["center_x_pc"], placement["center_y_pc"], placement["center_z_pc"])
        brightStars.FillContext(center, skeleton.shape, min_luminosity_sol=level,
                                star_mass_limit_sol=mass_limit).apply(system_config)
    return pc_to_ly(placement["galactic_radius_pc"])


def _generate_system(system_config, galactic_center_dist_ly=None, **kwargs):
    """Generates a `StarSystem` from a validated recipe (`kwargs`: the
    `StarSystem` arguments beyond it); a generation failure is the
    request's fault, so a 400, not a 500."""
    try:
        return StarSystem(system_config=system_config, galactic_center_dist_ly=galactic_center_dist_ly, **kwargs)
    except Exception as exc:
        # Generation can reject an internally-inconsistent recipe (e.g. an
        # impossible num_orbits/planet-class combination).
        raise ApiError(f"system generation failed: {exc}")


@bp.route("/systems", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def create_system():
    """
    `POST /api/systems` -- a generation "recipe" body (see
    `SYSTEM_CONFIG_ALLOWED_FIELDS`), the same shape `systemGen.py
    --system-file` already takes (`SystemConfig.from_dict`), plus optional
    `sector_id` and `position`. Without `sector_id` the new system is
    standalone. With it, the system joins that sector
    (`store.add_system_to_sector`): placed clear of every stored system's
    Hill sphere, or at `position` (`[x, y, z]` light-years from the
    sector's center, inside it), with the sector's stellar population
    when it is in the galaxy. Returns the new `star_systems.id` (and
    `position` when placed in a sector).
    """
    body = given(parse_body(SystemCreate, require_json_body()))
    sector_id, position = body.pop("sector_id", None), body.pop("position", None)
    if sector_id is not None:
        try:
            ids.parse(ids.SECTOR, sector_id)
        except ids.IdError:
            raise ApiError(f"sector_id must be a sector ID, got {sector_id!r}")
        found = ids.row_id(get_db(), ids.SECTOR, sector_id)
        if found is None:
            raise ApiError(f"no sector with ID {sector_id}", 404)
        sector_id = found
    if position is not None:
        position = tuple(position)

    job_id, made = run_queued(create_system_job, _resolve_requested_write_db_config(), body, sector_id, position)
    if made is None:
        audit("system.create", target=f"job:{job_id}", detail=f"queued star_type={body.get('star_type')!r}")
        return accepted(job_id)
    system_id = made["id"]
    audit("system.create", target=f"system:{system_id}",
          detail=f"star_type={body.get('star_type')!r} sector_id={sector_id!r}")
    return jsonify(made), 201


def create_system_job(config, body, sector_id, position):
    """The queued body of `POST /api/systems`: generates the system and
    stores it. Returns `{"id", ...}`, with `sector_id` and `position` when
    placed in a sector."""
    system_config = SystemConfig.from_dict(body)
    conn = store.open_write(config)
    try:
        with conn:
            if sector_id is None:
                system_id = store.insert_star_system(conn, _generate_system(system_config), system_config)
                store.assign_uids(conn, system_ids=[system_id])
                placed = None
            else:
                dist_ly = _sector_generation_context(conn, sector_id, system_config)
                star_system = _generate_system(system_config, dist_ly)
                try:
                    system_id, placed = store.add_system_to_sector(conn, sector_id, star_system, system_config,
                                                                 position=position)
                except ValueError as exc:
                    raise ApiError(str(exc), status_code=409 if position is None else 400)
    finally:
        conn.close()
    result = {"id": system_id}
    if placed is not None:
        result.update(sector_id=sector_id, position=list(placed))
    return result


_NAME_CLASH_LABELS = {
    "sectors": "a sector", "star_systems": "a star system", "stars": "a star",
}


def _require_unique_name(conn, name, exclude):
    """Raises a 409 when a sector, system or star other than the rows in
    `exclude` (`(table, id)` pairs) is already called `name`. Planet and
    moon names aren't checked (`store.name_in_use`)."""
    clash = store.name_in_use(conn, name, exclude=exclude)
    if clash is not None:
        raise ApiError(f"{_NAME_CLASH_LABELS[clash]} is already named {name!r}", status_code=409)


def _system_rename_exclusions(conn, system_id):
    """The system itself plus its single star, which shares its name --
    both are renamed together, so neither is a clash."""
    exclude = [("star_systems", system_id)]
    exclude.extend(
        ("stars", row["id"])
        for row in conn.execute(
            "SELECT id FROM stars WHERE star_system_id = ? AND role = 'single'", (system_id,)
        ).fetchall()
    )
    return exclude


@bp.route("/systems/<uid:system_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def update_system(system_id):
    """
    `PATCH /api/systems/<id>` -- `{"name": str}` renames a system, and
    every star, planet and moon still named after it
    (`store.rename_star_system`); 409 if a sector, system or star already has
    the name.

    `{"regenerate": recipe}` (the `POST /api/systems` recipe fields, `{}`
    for a fresh roll with defaults) replaces the system's stars, planets,
    moons, asteroid belts and comets with a newly generated set, keeping
    its id, name, sector, position and links
    (`store.replace_star_system_content`). A sector system keeps its
    sector's stellar population. 409 when the system was built around a
    pre-placed bright star, or hosts facilities unless `"drop_facilities":
    true` (they would be deleted with the bodies). Both may be sent
    together; the rename happens first.
    """
    parsed = parse_body(SystemPatch, require_json_body())
    body = given(parsed)
    name = parsed.name
    recipe = given(parsed.regenerate) if parsed.regenerate is not None else None
    drop_facilities = parsed.drop_facilities
    if recipe is not None:
        body["regenerate"] = recipe

    # A rename alone is one quick row write and stays in the request. A
    # regeneration generates a system, so it runs on the queue (PERF.24),
    # with any rename sent along, in the same transaction.
    config = _resolve_requested_write_db_config()
    if recipe is None:
        job_id, final_name = None, update_system_job(config, system_id, name, None, drop_facilities)
    else:
        job_id, final_name = run_queued(update_system_job, config, system_id, name, recipe, drop_facilities)
    detail = {key: body[key] for key in ("name", "drop_facilities") if key in body}
    if recipe is not None:
        detail["regenerate"] = recipe
    if final_name is None:
        audit("system.update", target=f"system:{system_id}", detail=f"queued job={job_id}: {detail}")
        return accepted(job_id)
    audit("system.update", target=f"system:{system_id}", detail=str(detail))
    return jsonify({"status": "ok", "id": system_id, "name": final_name, "regenerated": recipe is not None})


def update_system_job(config, system_id, name, recipe, drop_facilities):
    """The queued body of `PATCH /api/systems/<id>`: the rename, then the
    regeneration, in one transaction. Returns the system's final name."""
    conn = store.open_write(config)
    try:
        with conn:
            row = conn.execute("SELECT id, sector_id FROM star_systems WHERE id = ?", (system_id,)).fetchone()
            if row is None:
                raise ApiError(f"no such system: {system_id}", status_code=404)
            if name is not None:
                _require_unique_name(conn, name, _system_rename_exclusions(conn, system_id))
                store.rename_star_system(conn, system_id, name)
            if recipe is not None:
                blockers = store.system_content_blockers(conn, system_id)
                if blockers["bright_star"]:
                    raise ApiError("this system is built around a pre-placed bright star, so its content "
                                   "can't be regenerated", status_code=409)
                if blockers["facilities"] and not drop_facilities:
                    raise ApiError(f"this system hosts {blockers['facilities']} facilities, which regenerating "
                                   "would delete; send \"drop_facilities\": true to go ahead", status_code=409)
                system_config = SystemConfig.from_dict(recipe)
                dist_ly = (_sector_generation_context(conn, row["sector_id"], system_config)
                           if row["sector_id"] is not None else None)
                store.replace_star_system_content(conn, system_id, _generate_system(system_config, dist_ly),
                                                system_config)
            return conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"]
    finally:
        conn.close()


@bp.route("/stars/<uid:star_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def rename_star(star_id):
    """`PATCH /api/stars/<id>` `{"name": str}` -- renames a star. A single
    star shares its system's name, so this renames the system too; a
    binary's star is renamed on its own, with the planets and moons named
    after it (`store.rename_star`). 409 if a sector, system or star already has the
    name."""
    name = parse_body(SystemRename, require_json_body()).name

    conn = _write_conn()
    try:
        with conn:
            row = conn.execute("SELECT star_system_id, role FROM stars WHERE id = ?", (star_id,)).fetchone()
            if row is None:
                raise ApiError(f"no such star: {star_id}", status_code=404)
            if row["role"] == "single":
                exclude = _system_rename_exclusions(conn, row["star_system_id"])
            else:
                exclude = [("stars", star_id)]
            _require_unique_name(conn, name, exclude)
            store.rename_star(conn, star_id, name)
    finally:
        conn.close()

    audit("star.rename", target=f"star:{star_id}", detail=str({"name": name}))
    return jsonify({"status": "ok", "id": star_id, "star_system_id": row["star_system_id"], "name": name})


def _rename_planet_or_moon(table, kind, body_id):
    """Shared body of `rename_planet`/`rename_moon`."""
    name = parse_body(BodyRename, require_json_body()).name

    conn = _write_conn()
    try:
        with conn:
            row = conn.execute(f"SELECT star_system_id FROM {table} WHERE id = ?", (body_id,)).fetchone()
            if row is None:
                raise ApiError(f"no such {kind}: {body_id}", status_code=404)
            _require_unique_name(conn, name, [(table, body_id)])
            store.rename_body(conn, table, body_id, name)
    finally:
        conn.close()

    audit(f"{kind}.rename", target=f"{kind}:{body_id}", detail=str({"name": name}))
    return jsonify({"status": "ok", "id": body_id, "star_system_id": row["star_system_id"], "name": name})


@bp.route("/planets/<uid:planet_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def rename_planet(planet_id):
    """`PATCH /api/planets/<id>` `{"name": str}` -- renames one planet. Its
    moons keep their names. 409 if a sector, system or star already has the name."""
    return _rename_planet_or_moon("planets", "planet", planet_id)


@bp.route("/moons/<uid:moon_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def rename_moon(moon_id):
    """`PATCH /api/moons/<id>` `{"name": str}` -- renames one moon. 409 if
    a sector, system or star already has the name."""
    return _rename_planet_or_moon("moons", "moon", moon_id)


@bp.route("/systems/<uid:system_id>", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def delete_system(system_id):
    """`DELETE /api/systems/<id>` -- no request body. Its stars/planets/
    moons/belts are deleted with it (all `ON DELETE CASCADE`, per
    `schema.sql`). Bumps the system's sector's `modified_at`, since a
    deleted row leaves no timestamp of its own behind (`schema.sql`'s
    "v27" header note)."""
    conn = _write_conn()
    try:
        with conn:
            row = conn.execute("SELECT sector_id FROM star_systems WHERE id = ?", (system_id,)).fetchone()
            deleted = conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,)).rowcount > 0
            if deleted:
                store.touch_sector(conn, row["sector_id"])
    finally:
        conn.close()

    if not deleted:
        raise ApiError(f"no such system: {system_id}", status_code=404)
    audit("system.delete", target=f"system:{system_id}")
    return jsonify({"status": "ok"})


# ---------------------------------------------------------------------
# Facilities (schema v42): starbases, colonies and
# outposts. See planetgen/population/facilities.py for the placement rules.
# ---------------------------------------------------------------------

@bp.route("/facilities/<uid:facility_id>")
def facility(facility_id):
    """`GET /api/facilities/<id>` -- one facility."""
    found = facility_detail(get_db(), facility_id)
    if found is None:
        raise ApiError(f"no such facility: {facility_id}", status_code=404)
    return jsonify(found)


@bp.route("/systems/<uid:system_id>/facilities")
def system_facilities(system_id):
    """`GET /api/systems/<id>/facilities` -- every facility in a system.
    404 for an unknown system."""
    db = get_db()
    if db.execute("SELECT 1 FROM star_systems WHERE id = ?", (system_id,)).fetchone() is None:
        raise ApiError(f"no such system: {system_id}", status_code=404)
    return jsonify({"items": facilities_for_system(db, system_id)})


@bp.route("/sectors/<uid:sector_id>/facilities")
def sector_facilities(sector_id):
    """`GET /api/sectors/<id>/facilities` -- stand-alone facilities parked
    in a sector and those on its asteroid fields. 404 for an unknown
    sector."""
    db = get_db()
    if db.execute("SELECT 1 FROM sectors WHERE id = ?", (sector_id,)).fetchone() is None:
        raise ApiError(f"no such sector: {sector_id}", status_code=404)
    return jsonify({"items": facilities_in_sector(db, sector_id)})


@bp.route("/galaxy/bright-stars")
def galaxy_bright_stars_in_cell():
    """`GET /api/galaxy/bright-stars?ring=&layer=&slot=[&all=1]` -- the
    pre-placed bright stars in one sector cell, brightest first
    (`queryDb.bright_stars_in_sector`): only those not yet built into a
    system unless `all=1`. `{"items": []}` when no scatter has run (or
    the database predates them)."""
    try:
        ring, layer, slot = (int(request.args[key]) for key in ("ring", "layer", "slot"))
    except (KeyError, ValueError):
        raise ApiError("'ring', 'layer' and 'slot' must be integers")
    unfilled_only = request.args.get("all") not in ("1", "true")
    try:
        items = bright_stars_in_sector(get_db(), ring, layer, slot, unfilled_only=unfilled_only)
    except pymysql.err.ProgrammingError:  # no bright_stars table yet (before v43)
        items = []
    return jsonify({"items": items})


@bp.route("/galaxy/uncharted")
def galaxy_uncharted_sector():
    """`GET /api/galaxy/uncharted?ring=&layer=&slot=` -- one sector cell's
    place and what the scatters left in it, for opening a cell nothing was
    generated in (MAP.162): `designation`, `center_pc`, `edge_pc`,
    `sector_id` (the generated sector there, or `null`), `stars` (its
    waiting bright stars, `queryDb.bright_stars_in_sector`) and `scattered`
    (its unbuilt scattered phenomena, `queryDb.uncharted_sector_contents`).
    A cell outside the grid is a 400."""
    try:
        ring, layer, slot = (int(request.args[key]) for key in ("ring", "layer", "slot"))
    except (KeyError, ValueError):
        raise ApiError("'ring', 'layer' and 'slot' must be integers")
    if ring < 0:
        raise ApiError("ring must be >= 0")
    conn = get_db()
    skeleton = get_galaxy_shape(conn)
    edge_pc = skeleton.edge_pc if skeleton else float(tuning.DEFAULT_SECTOR_EDGE_PC)
    try:
        cell = describe_sector_cell(ring, layer, slot, edge_pc)
    except OverflowError:
        raise ApiError("that cell is too far out")
    except ValueError as err:
        raise ApiError(str(err))
    return jsonify({
        "ring_index": ring, "layer_index": layer, "ring_slot_index": slot, "designation": cell["designation"],
        "center_pc": cell["cartesian_pc"], "edge_pc": edge_pc, "sector_id": get_sector_id_at(conn, ring, layer, slot),
        **uncharted_sector_contents(conn, ring, layer, slot),
    })


@bp.route("/facilities/orbit")
def facility_orbit():
    """`GET /api/facilities/orbit?host_type=star|planet|moon&host_id=N[&distance_km=X]`
    -- the orbit (distance, period, speed) an orbital facility would get,
    without saving anything, so a form can show it first, plus
    `min_distance_km`/`max_distance_km`, the orbits the host allows (just
    above its surface to the edge of its sphere of influence)."""
    host_type = request.args.get("host_type", "")
    host_id = _host_row(host_type, request.args.get("host_id", ""))
    try:
        raw_distance = request.args.get("distance_km")
        distance_km = None if raw_distance in (None, "") else float(raw_distance)
    except ValueError:
        raise ApiError("distance_km must be a number")
    try:
        return jsonify(store.facility_orbit(get_db(), host_type, host_id, distance_km))
    except store.FacilityError as exc:
        raise ApiError(str(exc), status_code=404 if exc.not_found else 400)


@bp.route("/facilities", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def create_facility():
    """`POST /api/facilities` `{"name", "kind", "placement", "host_type",
    "host_id"}` (required) plus optional `distance_km`/`phase_deg`
    (orbital), `offset_ly` (stand-alone, `[x, y, z]` from the sector's
    center) and `description` -- adds a facility (`store.add_facility`).
    400 when the placement rules refuse it, 404 when the host is missing."""
    body = given(parse_body(FacilityCreate, require_json_body()))
    body["host_id"] = _host_row(body["host_type"], body["host_id"])
    conn = _write_conn()
    try:
        with conn:
            facility_id = store.add_facility(
                conn, body["name"].strip(), body["kind"], body["placement"], body["host_type"], body["host_id"],
                distance_km=body.get("distance_km"), phase_deg=body.get("phase_deg"),
                offset_ly=body.get("offset_ly"), description=body.get("description"),
            )
    except store.FacilityError as exc:
        raise ApiError(str(exc), status_code=404 if exc.not_found else 400)
    finally:
        conn.close()

    audit("facility.create", target=f"facility:{facility_id}",
          detail=str({key: body[key] for key in ("name", "kind", "placement", "host_type", "host_id")}))
    return jsonify({"id": facility_id}), 201


@bp.route("/facilities/<uid:facility_id>", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def delete_facility(facility_id):
    """`DELETE /api/facilities/<id>` -- no request body."""
    conn = _write_conn()
    try:
        with conn:
            deleted = store.delete_facility(conn, facility_id)
    finally:
        conn.close()
    if not deleted:
        raise ApiError(f"no such facility: {facility_id}", status_code=404)
    audit("facility.delete", target=f"facility:{facility_id}")
    return jsonify({"status": "ok"})


# ---------------------------------------------------------------------
# Wiki publishing (schema.sql's "v22" header note) -- both routes below
# share the same request shape and backend-resolution/error-mapping
# helpers, differing only in where their page content comes from
# (`systemRender.render_system_text` for a system, `_sector_wiki_content`
# for a sector's summary -- both built on the fly) and which column the
# resulting URL is written back to.
# ---------------------------------------------------------------------

def _require_wiki_configured(backend):
    """Raises a 501 when `backend` has no `base_url`/credentials
    configured at all -- rather than queuing an upload that can only
    fail on an empty URL."""
    if not current_app.config["WIKI_CONFIG"][backend]["configured"]:
        raise ApiError(f"wiki publishing is not configured for {backend!r}", status_code=501)


def wiki_settings(backend):
    """`backend`'s wiki settings as the worker resolves them itself
    (`config.py`'s `_wiki_config`: environment, then `config.json`), so
    the wiki's token or password never travels through Redis."""
    from planetgen.api.config import _wiki_config
    settings = _wiki_config(get_settings().wiki)[backend]
    if not settings["configured"]:
        raise ApiError(f"wiki publishing is not configured for {backend!r}", status_code=501)
    return settings


def _wiki_client(backend):
    """
    Builds a `planetgen.wiki.WikiClient` for `backend` from `wiki_settings`.

    Args:
        backend (str): `"wikijs"` or `"mediawiki"`.

    Returns:
        planetgen.wiki.WikiClient

    Raises:
        ApiError: 501 if `backend` isn't configured.
    """
    settings = wiki_settings(backend)
    if backend == "wikijs":
        return WikiClient(backend="wikijs", base_url=settings["base_url"], api_token=settings["api_token"])
    return WikiClient(
        backend="mediawiki", base_url=settings["base_url"],
        username=settings["username"], password=settings["password"],
    )


def _open_read_for_job(config):
    """A read-only connection for a queued job, as `get_db` opens one."""
    try:
        return open_readonly(config, statement_timeout_s=None)
    except SystemExit as exc:
        raise ApiError(f"{DATABASE_UNAVAILABLE} ({exc})", status_code=503)


def _wiki_upload_request(body):
    """A `POST .../wiki` body checked as `schemas.WikiUpload`: `{"backend":
    "wikijs" | "mediawiki", "path": str}`. `path` is the target page's
    path/slug for `"wikijs"` (which addresses a page separately from its
    title -- see `planetgen/wiki/wikijs.py`) and required for it;
    `"mediawiki"` has no such separate concept (its title *is* its address,
    see `planetgen/wiki/mediawiki.py`), so `path` is accepted but ignored.

    Returns:
        tuple[str, str or None]: `(backend, path)` -- `path` stripped, or
            `None` if not given/blank.
    """
    upload = parse_body(WikiUpload, body)
    return upload.backend, upload.path


def _create_wiki_page(client, path, title, content):
    """Wraps `client.create_page`, mapping `planetgen.wiki`'s own exception
    hierarchy onto this API's status codes -- shared by both routes below
    so neither duplicates the mapping.

    Raises:
        ApiError: 409 (a page already exists there) or 502 (auth
            rejected, or the wiki couldn't be reached/errored).
    """
    try:
        return client.create_page(path=path, title=title, content=content)
    except WikiClientPageExistsError as exc:
        raise ApiError(str(exc), status_code=409)
    except (WikiClientAuthError, WikiClientRequestError) as exc:
        raise ApiError(str(exc), status_code=502)


@bp.route("/systems/<uid:system_id>/wiki", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def upload_system_wiki(system_id):
    """
    `POST /api/systems/<id>/wiki` `{"backend": "wikijs"|"mediawiki",
    "path": str}` -- publishes this system's page -- rendered now
    from its database rows, Markdown for `wikijs` and wikitext for
    `mediawiki` (`planetgen/db/render.py`) -- to the chosen wiki,
    then records the
    new page's URL on `star_systems.wikijs_url`/`mediawiki_url` (see
    `schema.sql`'s "v22" header note) so `html/system.py` can swap its
    description section for a link to it afterward.

    A real database read (the system's own name/content, via `get_db()`'s
    request-scoped read-only connection) but not a database *write* until
    after the wiki call succeeds -- unlike every other write route above,
    which mutate the local database directly, this one's real effect is
    external, so `_write_conn()` is only opened once there is a URL worth
    persisting.
    """
    body = require_json_body()
    backend, path = _wiki_upload_request(body)
    _require_wiki_configured(backend)
    job_id, page = run_queued(upload_system_wiki_job, _resolve_requested_db_config(),
                              _resolve_requested_write_db_config(), system_id, backend, path)
    if page is None:
        audit("system.wiki_upload", target=f"system:{system_id}", detail=f"backend={backend!r} queued job={job_id}")
        return accepted(job_id)
    audit("system.wiki_upload", target=f"system:{system_id}", detail=f"backend={backend!r} path={page['path']!r}")
    return jsonify(page), 201


def upload_system_wiki_job(read_config, write_config, system_id, backend, path):
    """The queued body of `POST /api/systems/<id>/wiki`: renders the
    page, creates it on the wiki, records its URL. Returns the page."""
    client = _wiki_client(backend)
    db = _open_read_for_job(read_config)
    try:
        try:
            system = query_system_detail(db, system_id)
            content = render_system_text(db, system_id, "markdown" if backend == "wikijs" else "wikitext")
        except ValueError:
            raise ApiError(f"no such system: {system_id}", status_code=404)
    finally:
        db.close()

    page_path = path if backend == "wikijs" else system["name"]
    page = _create_wiki_page(client, page_path, system["name"], content)

    url_column = "wikijs_url" if backend == "wikijs" else "mediawiki_url"
    conn = store.open_write(write_config)
    try:
        with conn:
            conn.execute(f"UPDATE star_systems SET {url_column} = ? WHERE id = ?", (page.url, system_id))
    finally:
        conn.close()
    return {"id": page.id, "path": page.path, "title": page.title, "url": page.url}


def _sector_wiki_content(sector):
    """
    Builds a sector-summary page's content, fresh, in both Markdown (for
    `wikijs`) and wikitext (for `mediawiki`) -- like a system's own page
    (`systemRender.render_system_text`), nothing is stored
    (see `schema.sql`'s "v22" header note), so this is rendered on the
    fly from `sector`'s already-queried detail (`queryDb.sector_detail`'s
    shape) at upload time: the sector's name/size, and one table row per
    system placed in it, with the same star-type fallback
    `html/sector.py`'s own systems table already uses (a `'close'` pair's
    merged `binary_type` if it has one, else each of its own stars'
    `star_type` joined by `" / "`).

    Args:
        sector (dict): `queryDb.sector_detail`'s return shape.

    Returns:
        tuple[str, str]: `(markdown_content, wikitext_content)`.
    """
    def star_type_for(system):
        if system["is_binary"] and system.get("binary_type"):
            return system["binary_type"]
        return " / ".join(star["star_type"] for star in system["stars"]) if system["stars"] else ""

    systems = sector["systems"]

    md_rows = "\n".join(
        f'| {s["name"]} | {s["quadrant"] or ""} | {"Yes" if s["is_binary"] else "No"} | '
        f'{star_type_for(s)} | {s["location"] or ""} |'
        for s in systems
    ) or "| *(none)* | | | | |"
    markdown_content = (
        f'# {sector["name"]}\n\n'
        f'**Edge:** {format_distance_ly(sector["edge_ly"])}  \n'
        f'**Systems:** {sector["system_count"]}\n\n'
        "## Systems\n\n"
        "| Name | Octant | Binary | Star type | Location |\n"
        "|---|---|---|---|---|\n"
        f"{md_rows}\n"
    )

    wiki_rows = "\n|-\n".join(
        f'| {s["name"]} || {s["quadrant"] or ""} || {"Yes" if s["is_binary"] else "No"} || '
        f'{star_type_for(s)} || {s["location"] or ""}'
        for s in systems
    ) or "| ''(none)'' ||  ||  ||  || "
    wikitext_content = (
        f'= {sector["name"]} =\n\n'
        f"'''Edge:''' {format_distance_ly(sector['edge_ly'])}\n\n"
        f"'''Systems:''' {sector['system_count']}\n\n"
        "== Systems ==\n\n"
        '{| class="wikitable"\n'
        "! Name !! Octant !! Binary !! Star type !! Location\n"
        "|-\n"
        f"{wiki_rows}\n"
        "|}\n"
    )
    return markdown_content, wikitext_content


@bp.route("/sectors/<uid:sector_id>/wiki", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def upload_sector_wiki(sector_id):
    """
    `POST /api/sectors/<id>/wiki` `{"backend": "wikijs"|"mediawiki",
    "path": str}` -- publishes a freshly generated sector-summary page
    (see `_sector_wiki_content`) to the chosen wiki, then records the new
    page's URL on `sectors.wiki_url` (a single column, not one per
    backend the way `star_systems` has -- see `schema.sql`'s "v22" header
    note) so `html/sector.py` can link to it afterward. The same column
    can also be set directly, without uploading anything, via
    `PATCH /api/sectors/<id>`'s own `wiki_url` field (`html/admin.py`'s
    manual-link admin section).
    """
    body = require_json_body()
    backend, path = _wiki_upload_request(body)
    _require_wiki_configured(backend)
    job_id, page = run_queued(upload_sector_wiki_job, _resolve_requested_db_config(),
                              _resolve_requested_write_db_config(), sector_id, backend, path)
    if page is None:
        audit("sector.wiki_upload", target=f"sector:{sector_id}", detail=f"backend={backend!r} queued job={job_id}")
        return accepted(job_id)
    audit("sector.wiki_upload", target=f"sector:{sector_id}", detail=f"backend={backend!r} path={page['path']!r}")
    return jsonify(page), 201


def upload_sector_wiki_job(read_config, write_config, sector_id, backend, path):
    """The queued body of `POST /api/sectors/<id>/wiki`. Returns the page."""
    client = _wiki_client(backend)
    db = _open_read_for_job(read_config)
    try:
        try:
            sector = query_sector_detail(db, sector_id)
        except ValueError:
            raise ApiError(f"no such sector: {sector_id}", status_code=404)
    finally:
        db.close()

    markdown_content, wikitext_content = _sector_wiki_content(sector)
    content = markdown_content if backend == "wikijs" else wikitext_content
    page_path = path if backend == "wikijs" else sector["name"]
    page = _create_wiki_page(client, page_path, sector["name"], content)

    conn = store.open_write(write_config)
    try:
        with conn:
            conn.execute("UPDATE sectors SET wiki_url = ? WHERE id = ?", (page.url, sector_id))
    finally:
        conn.close()
    return {"id": page.id, "path": page.path, "title": page.title, "url": page.url}
