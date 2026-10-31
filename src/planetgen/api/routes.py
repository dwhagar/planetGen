# planetgen/api/routes.py

"""
JSON endpoints over the planetGen database.

The read endpoints each open a strictly-read-only connection
(`db.query.open_readonly` -- same `file:...?mode=ro`-equivalent guarantee
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
separate, write-capable database account -- see `docs/reference/api.md`'s
"Write endpoints" section and `planetgen/admin/auth.py`.
"""

import math

import pymysql
from flask import Blueprint, current_app, g, jsonify, request, url_for

from planetgen.db import search_values
from planetgen.db import children as children_db
from planetgen.describe.sheet import sheet as object_sheet
from planetgen.web.maps.systemscene import build_scene
from planetgen.generation import run_galaxy
from planetgen.queue import api_jobs, pagejobs
from planetgen.util import draw, generator_epoch

from planetgen.db.shell_specs import galaxy_shell_specs
from planetgen.db.query import (
    facilities_for_system,
    facilities_in_sector,
    facility_detail,
    FACILITY_SORTS,
    list_facilities,
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
    galaxy_changes,
    galaxy_density_shape,
    galaxy_placed_phenomena,
    galaxy_placed_sectors,
    sectors_made,
    galaxy_locate,
    galaxy_stage,
    galaxy_edge_pc,
    galaxy_nebula_cover,
    galaxy_window,
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
    scattered_phenomenon_detail as query_scattered_phenomenon_detail,
    uncharted_system_detail as query_uncharted_system_detail,
    resolve_object,
    nebula_shape as query_nebula_shape,
    nebula_surroundings as query_nebula_surroundings,
    search as run_search,
    parse_star_count,
    sector_detail as query_sector_detail,
    system_detail as query_system_detail,
)
from planetgen.db import near as near_search
from planetgen.db import region_stats
from planetgen.names import naming_key
from planetgen.db import cloudcontents, containers as containersDb, nodes, objectlist
from planetgen.db import stats as dbStats
from planetgen.db import scatter_systems, store
from planetgen.generation import bright_stars as brightStars, limits as generationLimits, stored_trees
from planetgen import tuning
from planetgen.admin import maintenance
from planetgen.population import facility_ratings, facility_suggest, facility_types
from planetgen.galaxy import objectref as object_ref, scatter_limits, superbubble_field, version_check
from planetgen.galaxy.profile import sector_profile
from planetgen.galaxy.profile_report import profile_report
from planetgen.galaxy.resolution import cell_resolution
from planetgen.db.store import (get_galaxy_bounds, get_galaxy_shape, get_sector_id_at, list_databases, resolve_database)
from planetgen.util.settings import get_settings
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.profile import profile_at
from planetgen.galaxy import drill, spherical
from planetgen.galaxy.spherical import describe_sector_cell, sector_address_at
from planetgen.galaxy.window import DEFAULT_WINDOW
from planetgen.generation.system import StarSystem
from planetgen.db.render import FORMATS as SYSTEM_TEXT_FORMATS
from planetgen.db.render import describe_system, render_system_text
from planetgen.physics.units import ly_to_milliparsecs, pc_to_ly
from planetgen.util import detmath

from . import ids
from .authz import audit, require_admin
from .common import ApiError, request_json, require_json_body
from .version import API_VERSION
from .schemas import (
    BlockGenerateRequest, BodyRename, FacilityCreate, FacilityUpdate, NeighborhoodRequest, SectorCreate, SectorUpdate, SystemCreate, SystemPatch,
    SystemRename, given, parse_body,
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

    `db.query.open_readonly` raises a bare `SystemExit` when the
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


def _parse_stars(query_args):
    """The `stars` query parameter (UX.239) as `(least, most)`, or `None` when absent; anything `parse_star_count` refuses is an `ApiError`."""
    try:
        return parse_star_count(query_args.get("stars"))
    except ValueError as exc:
        raise ApiError(str(exc))


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
    parameter pair into a `(min_km, max_km)` tuple for `db.query.search`'s
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
        "generator_version": generator_epoch.generator_epoch(),
    }
    state = maintenance.read()
    if state["enabled"]:
        # Open to anyone while the site is closed (`maintenance_gate.OPEN_PATHS`): the update's own reload check
        # asks while its maintenance mode is still on (OPS.75).
        body["maintenance"] = {"enabled": True, "cause": state["cause"], "reason": state["reason"]}
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
    `after=<sector ID>` starts the page after that sector in the default order and the `name` and `distance` sorts
    (PERF.137: any page costs the same, and `offset` is ignored); the other sorts ignore it.
    """
    limit, offset = _paginate(request.args)
    sort = None
    descending = False
    if request.args.get("sort") or request.args.get("order"):
        sort, descending = _parse_sort(request.args, SECTOR_SORTS)
    quadrants = [v for v in request.args.getlist("quadrant") if v]
    db = get_db()
    after = None
    if request.args.get("after"):
        after = ids.row_id(db, "sector", request.args["after"])
        if after is None:
            raise ApiError(f"after must be a sector ID, got {request.args['after']!r}")
    total = count_sectors(db, quadrants=quadrants)
    body = {
        "items": list_sectors(db, limit=limit, offset=offset, sort=sort, descending=descending,
                              quadrants=quadrants, after=after),
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
    conn = get_db()
    try:
        detail = query_sector_detail(conn, sector_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    detail["superbubble"] = superbubble_of(conn, detail)
    detail["region"] = region_of(conn, (detail["center_x_pc"], detail["center_y_pc"], detail["center_z_pc"])
                                 if detail.get("placed") else None)
    detail["unresolved"] = unresolved_of(conn, detail)
    detail["scatter_tier"] = scatter_tier_of(conn, detail)
    return jsonify(detail)


def _parse_sector_id_filter(raw_sector_id):
    """
    Parses the `sector_id` query parameter `/api/systems` accepts: a
    sector's ID (one sector), the literal `"none"` (standalone systems --
    `db.query.NO_SECTOR`, `html/browse.py`'s own table of systems generated
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


MAX_START_LENGTH = 80
"""int: Longest `from=` text kept (the jump-by-value box)."""


@bp.route("/systems")
def systems():
    """
    The systems, by name, paginated; `star_type` and `sector_id` filter them, and so does `stars` (UX.239): a count (`stars=3`, exactly
    three), a count and a plus (`stars=3+`, three or more) or a name (`single`, `binary`, `trinary`, `quaternary`: exactly that many;
    `quintenary`: five or more), anything else a 400.
    The Systems tables (UX.41) also take `sort` (`name`, `sector`, `octant` or `stars`, the number of stars) with `order=asc|desc`, the filters `binary=yes|no`,
    `placement=sector|standalone` and `octant` (repeatable), and `facets=1`
    for the filter menus' option counts (`facets`: `placement`, `binary`,
    `octant`). `total` counts the systems that pass the filters. `after=<system ID>` starts the page after that system in the
    sort instead of skipping `offset` rows (PERF.75, PERF.206: page 40,000 costs what page 1 does, in the name, sector
    and octant sorts; a list filtered by star type or sector reads by `offset`). `from=<text>` is the
    jump-by-value box (PERF.75): with the `name` sort the list begins at the first name from that text on (a letter, or
    the start of a name), with the `sector` sort at the first sector from it on, in the sort's direction; `total`
    counts the systems from there and `offset` and `after` count from there too. The other sorts ignore it.
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
        "stars": _parse_stars(request.args),
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
    start = (request.args.get("from") or "").strip()[:MAX_START_LENGTH] or None
    rows = list_systems(db, star_type_prefix=star_type, sector_id=sector_id, limit=limit, offset=offset,
                        sort=sort, descending=descending, after=after, start=start, **filters)
    total = count_systems(db, star_type_prefix=star_type, sector_id=sector_id, start=start, sort=sort,
                          descending=descending, **filters)
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
    conn = get_db()
    try:
        detail = query_system_detail(conn, system_id)
        detail["region"] = region_of(conn, resolve_object(conn, "system", system_id)["positions"]["galaxy_pc"])
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


@bp.route("/systems/<uid:system_id>/description")
def system_description(system_id):
    """
    `GET /api/systems/<id>/description` -- the system's Document (API.25,
    `planetgen/describe/`): sections, fact and rating tables, compositions
    and prose, every number a typed quantity, built now from its database
    rows. The web page renders it, splitting it per body (`describe.model.split_by_body`); `/text` writes it as
    text on request.
    """
    try:
        return jsonify(describe_system(get_db(), system_id))
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404


@bp.route("/systems/<uid:system_id>/text")
def system_text(system_id):
    """
    `GET /api/systems/<id>/text?format=markdown|text` -- the
    system's description written as Markdown or plain text,
    on request only (`planetgen/db/render.py`; no page text is stored since
    schema v29). `format` defaults to `markdown`.
    """
    fmt = request.args.get("format", "markdown")
    if fmt not in SYSTEM_TEXT_FORMATS:
        raise ApiError(f"'format' must be one of: {', '.join(SYSTEM_TEXT_FORMATS)}")
    try:
        content = render_system_text(get_db(), system_id, fmt)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify({"id": system_id, "format": fmt, "content": content})


def _life_air(table, body_id):
    """Shared body of `planet_life_air` and `moon_life_air`."""
    from planetgen.generation import evolution
    from planetgen.physics import life_air, life_classes

    conn = get_db()
    row = conn.execute(f"SELECT star_system_id, life_class, life_class_variant, life_stage FROM {table} WHERE id = ?", (body_id,)).fetchone()
    if row is None:
        raise ApiError(f"no such {table[:-1]}: {body_id}", status_code=404)
    class_id = life_classes.living_class_id(row["life_class"], row["life_class_variant"], row["life_stage"])
    history = None
    if class_id in life_classes.BY_ID:
        system = store.load_star_system(conn, row["star_system_id"])
        bodies = [planet for planet in system.planets if hasattr(planet, "moons")]
        bodies += [moon for planet in bodies for moon in planet.moons]
        body = next((b for b in bodies if getattr(b, "db_id", None) == body_id and (b.is_moon == (table == "moons"))), None)
        if body is not None:
            ages = dict(zip(evolution.MILESTONE_KEYS, (getattr(body, field, None) for field in evolution.SCHEDULE_FIELDS)))
            history = life_air.air_history(body, [life_air.Living(life_classes.BY_ID[class_id], 1.0, ages.get("photosynthesis"))], ages)
    return jsonify({"id": ids.printed(conn, table[:-1], body_id), "life_class_id": class_id, "life_air": history})


@bp.route("/planets/<uid:planet_id>/life-air")
def planet_life_air(planet_id):
    """
    `GET /api/planets/<id>/life-air` -- the Life and air chart's data (UX.237, `physics.life_air.air_history`): the planet's air through the
    age of its star as far as the life class that lives on it changed it, worked out from the stored air on request and stored nowhere. `{"id",
    "life_class_id", "life_air"}`; `life_air` is `null` for a world with no living class, a giant, or an air the life did not change.
    Otherwise it has `age_gyr` (the star's age, Gyr), `producers_gyr` and `rose_gyr` (the ages the first producers began and the oxygen rose,
    Gyr, the latter `null` when no class makes any), `classes` (`{"life_class", "variant", "share", "producers_gyr"}` for each class whose life
    the air sums, with its share of the food production), `gases` (`{"formula", "name", "steps": [[age_gyr, partial_pressure_kpa], ...]}`, each
    value holding until the next age), `stages` (`{"key", "age_gyr"}` for each milestone the schedule reaches) and `speculative` (the
    ledger's `{"tier", "tweak"}` for a class whose numbers are tuned, else `null`).
    """
    return _life_air("planets", planet_id)


@bp.route("/moons/<uid:moon_id>/life-air")
def moon_life_air(moon_id):
    """`GET /api/moons/<id>/life-air` -- the same for a moon (see `planet_life_air`)."""
    return _life_air("moons", moon_id)


def _life_readings(table, body_id):
    """Shared body of `planet_life_readings` and `moon_life_readings`."""
    from planetgen.physics import constants, life_classes, life_readings

    conn = get_db()
    rows = conn.execute(f"SELECT * FROM {table} WHERE id = ?", (body_id,)).fetchall()
    if not rows:
        raise ApiError(f"no such {table[:-1]}: {body_id}", status_code=404)
    values = store.with_life_inputs(conn, rows, table == "moons")[0]
    values["is_moon"] = table == "moons"
    star = None
    star_row = None if values.get("star_id") is None else conn.execute(
        "SELECT mass_kg, age_gy, metallicity_fe_h FROM stars WHERE id = ?", (values["star_id"],)).fetchone()
    if star_row is not None:
        star = life_readings.star_of({"mass_sol": star_row["mass_kg"] / constants.SOLAR_MASS_TO_KG, "age_gyr": star_row["age_gy"],
                                      "fe_h": star_row["metallicity_fe_h"] or 0.0})
    found = [{"kind": reading.kind, "label": reading.label, "score": reading.score, "abundance": reading.abundance,
              "trend": reading.trend, "security": reading.security, "years_left_gyr": reading.years_left_gyr,
              "text": life_readings.sentence(reading)} for reading in life_readings.readings(values, star)]
    class_id = life_classes.living_class_id(values.get("life_class"), values.get("life_class_variant"), values.get("life_stage"))
    return jsonify({"id": ids.printed(conn, table[:-1], body_id), "life_class_id": class_id, "readings": found})


@bp.route("/planets/<uid:planet_id>/life-readings")
def planet_life_readings(planet_id):
    """
    `GET /api/planets/<id>/life-readings` -- the three words for each kind of life on the planet (UX.251, `physics.life_readings`),
    worked out on request from the stored body and its star and stored nowhere. `{"id", "life_class_id", "readings"}`; `readings` is
    empty for a gas giant and for a body whose kind scores are not stored, otherwise one entry for each of microbes, flora, fauna and
    sentient life in that order: `kind`, `label`, `score` (the stored 0 to 9), `abundance` (Absent, Rare, Occasional, Frequent, Abundant or
    Dominant), `trend` (Growing, Stable or Declining), `security` (Secure, Marginal, Precarious or Dying out), `years_left_gyr` (Gyr until the
    living class stops scoring on the world as its star brightens, `null` when it lasts or is not known) and `text` ("Microbes: abundant,
    stable, secure."). The trend and security of a kind that is Absent are `null`.
    """
    return _life_readings("planets", planet_id)


@bp.route("/moons/<uid:moon_id>/life-readings")
def moon_life_readings(moon_id):
    """`GET /api/moons/<id>/life-readings` -- the same for a moon (see `planet_life_readings`)."""
    return _life_readings("moons", moon_id)


@bp.route("/near")
def near():
    """
    Everything generated within a distance of a place (NAV.43). The place
    is `from` (an object reference such as `system:12`; a bare number is a
    system) or `point` (`x,y,z`, galaxy-frame parsecs); `distance` is in
    parsecs, up to 50. Optional: `kinds` (comma-separated), `limit`
    (default 50, at most 200) and `offset`. `stars` (UX.239; a count such as `3`, `3+` for three or more, or `single`, `binary`,
    `trinary`, `quaternary`, `quintenary`) keeps the systems with that many stars and what lives in them; phenomena and facilities
    outside a system are left out, and each `system` row carries its `star_count`.

    `nearest=K` (NAV.61) asks for the nearest K rows instead: the search tries 2, 4, 8, 16 and 32 pc and stops at the
    first sphere that holds K, so it is quick at the dense core as well as at the Sun; `distance` is then optional and
    caps how far it looks. A sphere that holds more than 100,000 systems is refused with a 400 that says to shrink the
    distance. The answer carries `systems_in_range` (the systems in the sphere, counted no further than 10,000) and
    `systems_in_range_capped` (true when there are more).
    """
    args = request.args
    conn = get_db()
    stars = _parse_stars(args)
    if (args.get("from") is None) == (args.get("point") is None):
        raise ApiError("give exactly one of from (an object reference) or point (x,y,z in parsecs)")
    if args.get("distance") is None and args.get("nearest") is None:
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
        if args.get("nearest") is not None:
            result = near_search.nearest(
                conn, place, _int_arg(args, "nearest", 1), kinds=kinds,
                max_distance_pc=args.get("distance", near_search.MAX_DISTANCE_PC), stars=stars)
        else:
            result = near_search.objects_within(
                conn, place, args["distance"], kinds=kinds,
                limit=_int_arg(args, "limit", near_search.DEFAULT_LIMIT), offset=_int_arg(args, "offset", 0), stars=stars)
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


def _float_arg(query_args, name):
    """A required finite number from the query string; anything else is an `ApiError`."""
    raw = query_args.get(name)
    try:
        value = float(raw)
    except (TypeError, ValueError):
        raise ApiError(f"{name} must be a number, got {raw!r}")
    if value != value or value in (float("inf"), float("-inf")):
        raise ApiError(f"{name} must be a finite number, got {raw!r}")
    return value


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
    system) -- see `db.query.nav_course` for the legs and
    `db.query.nav_between` for the availability rules, and docs/reference/api.md for
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
    as for `/api/nav`; `border=1` adds the uncharted cells next to them) -- see `db.query.nav_chart_plan`.
    Answers `unknown_hops`, `cells` (`[shell, band, wedge]`, in the order the route enters them, inside the
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
    `db.query._phenomenon_nav_key`-shaped string keys (a phenomenon
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
    (`db.query.sectors_made`; at most `MADE_SECTOR_LIMIT` items).
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
    NULL, in any `db.query._PHENOMENON_TABLES` table -- see `schema.sql`'s
    "v18"/"v21"/"v29" header notes) -- the phenomenon counterpart to `/api/galaxy/sectors`, plotted
    as small dots on the same `/galaxy` Galaxy Map. Not paginated,
    for the same reason `/api/galaxy/sectors` isn't.
    """
    return jsonify({"items": galaxy_placed_phenomena(get_db())})


@bp.route("/galaxy/shape")
def galaxy_shape():
    """
    The galaxy's stored density-skeleton shape (`planetgen plan`'s own
    output, `db.query.galaxy_density_shape`) -- the real spiral/disk/bulge
    model the Galaxy Map (`/galaxy`) shades its "expected density"
    cloud from, for whatever space hasn't actually been generated yet.
    `"shape"` is `null` when the skeleton has never been built (the map
    then falls back to its own generic illustrative gradient).
    `"bright_stars"` is `db.query.bright_star_scatter_status`: whether the
    bright-star scatter has run, its threshold and seed, and the default
    threshold a plain plan uses. `"version_warning"` (DB.7) is a sentence
    when the galaxy holds sectors generated by another PlanetGen version
    than the running one, else `null`. `"shells"` (ADM.28) is
    `shell_specs.galaxy_shell_specs`: the shell count, the outermost shell's reach
    and how many sectors are charted, with a row per shell (a sample past 41).
    """
    conn = get_db()
    skeleton = get_galaxy_shape(conn)
    return jsonify({"shape": galaxy_density_shape(conn), "bright_stars": bright_star_scatter_status(conn),
                    "version_warning": version_check.galaxy_warning(conn),
                    "shells": galaxy_shell_specs(conn, skeleton.edge_pc) if skeleton else None})


@bp.route("/galaxy/superbubbles")
def galaxy_superbubbles():
    """
    The galaxy's superbubbles (GEN.254, `galaxy.superbubble_field`): the cavities the Galaxy Map draws as faint
    shells, `{"superbubbles": [{"x", "y", "z", "radius_pc"}, ...]}` in galaxy-frame parsecs, the biggest first;
    empty before a plan. Nothing is stored: they are drawn from the galaxy's seed.
    """
    skeleton = get_galaxy_shape(get_db())
    if skeleton is None or skeleton.galaxy_seed is None:
        return jsonify({"superbubbles": []})
    reach = (skeleton.outer_shell_index + 1) * skeleton.edge_pc
    bubbles = sorted(superbubble_field.all_bubbles(skeleton.galaxy_seed, skeleton.shape, reach),
                     key=lambda bubble: (-bubble.radius_pc, bubble.index, bubble.n))
    return jsonify({"superbubbles": [
        {"x": round(b.center_pc[0], 1), "y": round(b.center_pc[1], 1), "z": round(b.center_pc[2], 1),
         "radius_pc": round(b.radius_pc, 1)} for b in bubbles]})


def region_of(conn, position_pc):
    """
    The named region holding a point (UX.231, `region_stats.region_of`): `{"name": "Vang Zarh Region", "parent": "Vang Reach"}`, or `None`
    for a point with no position, one in a cell with no generated sector (no region exists there), or a galaxy with no naming key (the
    name cannot be spelled). A sector, system or object answer carries it as `region`.
    """
    if position_pc is None:
        return None
    edge_pc = galaxy_edge_pc(conn)
    found = region_stats.region_of(conn, position_pc, edge_pc)
    return None if found is None else region_stats.labels(found["block"], naming_key.active_key(), edge_pc)


def superbubble_of(conn, detail):
    """
    The superbubble a placed sector's centre lies in (`sector_detail`'s `superbubble`), or `None`:
    `{"radius_pc", "ambient_density_cm3", "temperature_k", "centre_distance_pc"}`.
    """
    skeleton = get_galaxy_shape(conn)
    if skeleton is None or skeleton.galaxy_seed is None or not detail.get("placed"):
        return None
    center = (detail["center_x_pc"], detail["center_y_pc"], detail["center_z_pc"])
    bubble = superbubble_field.containing(skeleton.galaxy_seed, skeleton.shape, center)
    if bubble is None:
        return None
    return {"radius_pc": round(bubble.radius_pc, 1), "ambient_density_cm3": float("%.3g" % bubble.ambient_density_cm3),
            "temperature_k": bubble.temperature_k,
            "centre_distance_pc": round(detmath.dist(bubble.center_pc, center), 1)}


def scatter_tier_of(conn, detail):
    """
    The tier of the sector's band and the limits the galaxy scatter drew its stars at there (GEN.333): `{"tier",
    "min_luminosity_sol", "mass_limit_sol"}` (`mass_limit_sol` `None` when the scatter had no mass pass), or `None`
    for a sector outside the galaxy or before any scatter ran.
    """
    address = (detail.get("shell_index"), detail.get("band_index"), detail.get("wedge_index"))
    if not detail.get("placed") or None in address:
        return None
    limits = store.scatter_limits_of(conn)
    found = limits.limit_at(address[0], address[1])
    if found is None:
        return None
    return {"tier": scatter_limits.TIER_NAMES[limits.tier_at(address[0], address[1])],
            "min_luminosity_sol": float("%.6g" % found.min_luminosity_sol),
            "mass_limit_sol": None if found.mass_limit_sol is None else float("%.6g" % found.mass_limit_sol)}


def unresolved_of(conn, detail):
    """
    What a placed sector leaves unresolved (`sector_detail`'s `unresolved`, GEN.281), or `None` when the cell expected
    no more systems than the cell star ceiling and was filled whole: `{"ceiling_systems", "luminosity_cap_sol",
    "expected_systems", "unresolved_systems", "unresolved_mass_solar"}`, worked out from the address.
    """
    skeleton = get_galaxy_shape(conn)
    address = (detail.get("shell_index"), detail.get("band_index"), detail.get("wedge_index"))
    if skeleton is None or not detail.get("placed") or None in address:
        return None
    found = cell_resolution(sector_profile(address, skeleton), skeleton.cell_star_ceiling,
                            skeleton.fill_luminosity_cap_sol,
                            store.scatter_limits_of(conn, skeleton).mass_limit_at(address[0], address[1]))
    if not found.capped:
        return None
    return {"ceiling_systems": found.ceiling_systems, "luminosity_cap_sol": found.luminosity_cap_sol,
            "expected_systems": float("%.4g" % found.expected_systems),
            "unresolved_systems": float("%.4g" % found.unresolved_systems),
            "unresolved_mass_solar": float("%.4g" % found.unresolved_mass_solar)}


def _cell_query_address(skeleton):
    """
    The `(shell, band, wedge)` and the sector edge a `/galaxy/cell`-style query names: `?shell=&band=&wedge=` by
    address, `?x=&y=&z=` (parsecs, galaxy frame) by any point inside the cell. Uses the stored skeleton's edge length,
    else the standard 4 pc.

    Raises:
        ApiError: 400 for a missing, non-numeric or non-finite part, a negative shell, or a point too far out.
    """
    def number(name, cast):
        value = request.args.get(name)
        if value is None:
            return None
        try:
            return cast(value)
        except ValueError:
            raise ApiError(f"{name} must be a number")

    edge_pc = skeleton.edge_pc if skeleton else float(tuning.DEFAULT_SECTOR_EDGE_PC)
    shell, band, wedge = number("shell", int), number("band", int), number("wedge", int)
    x, y, z = number("x", float), number("y", float), number("z", float)
    # An address astronomically far out overflows the float geometry
    # (OverflowError) -- a bad request, not a server error.
    if None not in (shell, band, wedge):
        address = (shell, band, wedge)
    elif None not in (x, y, z) and all(math.isfinite(v) for v in (x, y, z)):
        try:
            address = sector_address_at((x, y, z), edge_pc)
        except OverflowError:
            raise ApiError("x, y and z are too far out")
    else:
        raise ApiError("give shell, band and wedge, or x, y and z")
    if address[0] < 0:
        raise ApiError("shell must be >= 0")
    return address, edge_pc


@bp.route("/galaxy/cell")
def galaxy_cell():
    """
    One sector cell of the spherical grid, whether or not anything was
    ever generated there: `?shell=I&band=J&wedge=K` names it by address,
    `?x=&y=&z=` (parsecs, galaxy frame) by any point inside it. Returns
    `spherical.describe_sector_cell` (center in Cartesian, cylindrical
    and spherical coordinates, bounds, 8 corners) plus `sector_id`, the
    generated sector there or `null`, and `in_galaxy`: whether the cell lies
    inside the planned galaxy's stored outline (`null` before any plan).
    Uses the stored skeleton's edge length, else the standard 4 pc.
    """
    conn = get_db()
    skeleton = get_galaxy_shape(conn)
    address, edge_pc = _cell_query_address(skeleton)
    try:
        cell = describe_sector_cell(*address, edge_pc)
    except OverflowError:
        raise ApiError("that cell is too far out")
    except ValueError as err:
        raise ApiError(str(err))
    cell["edge_pc"] = edge_pc
    cell["sector_id"] = get_sector_id_at(conn, *address)
    bounds = get_galaxy_bounds(conn)
    cell["in_galaxy"] = bounds.qualifies(address[0], address[1]) if bounds is not None else None
    return jsonify(cell)


@bp.route("/galaxy/profile")
def galaxy_profile():
    """
    The Sector Profile of any cell of the grid (API.27), charted or not, so an uncharted sector can be previewed:
    `?shell=I&band=J&wedge=K` by address, `?x=&y=&z=` by any point inside the cell. Read only: nothing is generated,
    stored or drawn from the run's random numbers; `galaxy.profile.sector_profile` is closed form in the stored
    skeleton and the address, so two calls give the same answer. See `profile_report.profile_report` for the body.
    `404` before a galaxy is planned (there is no shape to read), `400` for a bad query or a cell past the grid's
    reach (`spherical.MAX_ADDRESSABLE_SHELL`).
    """
    conn = get_db()
    skeleton = get_galaxy_shape(conn)
    if skeleton is None:
        raise ApiError("no galaxy has been planned, so there is no profile to read", 404)
    address, _edge_pc = _cell_query_address(skeleton)
    if address[0] > spherical.MAX_ADDRESSABLE_SHELL:
        raise ApiError(f"shell must be <= {spherical.MAX_ADDRESSABLE_SHELL}, the outermost shell of the grid's reach")
    try:
        profile = sector_profile(address, skeleton)
    except OverflowError:
        raise ApiError("that cell is too far out")
    except ValueError as err:
        raise ApiError(str(err))
    body = profile_report(profile)
    body["sector_id"] = get_sector_id_at(conn, *address)
    bounds = get_galaxy_bounds(conn)
    body["in_galaxy"] = bounds.qualifies(address[0], address[1]) if bounds is not None else None
    body["unresolved"] = unresolved_of(conn, {"placed": True, "shell_index": address[0], "band_index": address[1],
                                              "wedge_index": address[2]})
    resolution = cell_resolution(profile, skeleton.cell_star_ceiling, skeleton.fill_luminosity_cap_sol,
                                 store.bright_star_mass_limit(conn))
    body["estimate"] = run_galaxy.estimate_cell(profile.density, resolution.resolved_systems,
                                                config=_resolve_requested_db_config())
    return jsonify(body)


@bp.route("/galaxy/nodes")
def galaxy_nodes_route():
    """
    The 3D Galaxy Map's sector-block nodes (MAP.260) -- `keys` is a comma-separated list of `m.shell.band.wedge`
    keys (a sector's is `1.shell.band.wedge`; at most `nodes.MAX_NODES_PER_REQUEST`). See `db.nodes.galaxy_nodes`: each
    node is the brightest objects inside a block and what is stored about the block as a whole. Called by
    `/galaxy/nodes` (`planetgen/web/galaxy_views.py`), which caches every node on disk and only forwards the ones it
    does not already have.
    """
    node_keys = [key for key in (request.args.get("keys") or "").split(",") if key]
    if len(node_keys) > nodes.MAX_NODES_PER_REQUEST:
        raise ApiError(f"at most {nodes.MAX_NODES_PER_REQUEST} nodes per request")
    try:
        return jsonify(nodes.galaxy_nodes(get_db(), node_keys))
    except ValueError as exc:
        raise ApiError(str(exc))


@bp.route("/galaxy/locate")
def galaxy_locate_route():
    """
    The Galaxy Map address bar's name lookup: `?q=<part of a name>`.
    Returns `{"matches": [...]}`, sectors and systems whose name contains
    it, each with its sector address -- see `db.query.galaxy_locate`.
    Called by `/galaxy/locate` (`planetgen/web/galaxy_views.py`).
    """
    return jsonify({"matches": galaxy_locate(get_db(), request.args.get("q") or "")})


@bp.route("/galaxy/stage")
def galaxy_stage_route():
    """
    One Galaxy Map drill-down stage: `?at=m.shell.band.wedge` (a block key,
    `planetgen.galaxy.drill`), or no `at` for the galaxy. Returns how
    many generated sectors each child block holds, and at a level-3 block
    the generated sectors themselves -- see `db.query.galaxy_stage`.
    Called by `/galaxy/stage` (`planetgen/web/galaxy_views.py`), which caches
    each stage on disk.
    """
    try:
        return jsonify(galaxy_stage(get_db(), request.args.get("at") or None))
    except ValueError as exc:
        raise ApiError(str(exc))


@bp.route("/galaxy/window")
def galaxy_window_route():
    """
    A Galaxy Map data window of cells bigger than a sector: `?at=m.shell.band.wedge`
    (the middle cell) and `?g=` (cells across), read from the stored block
    totals in one query -- see `db.query.galaxy_window`. Called by
    `/galaxy/window` (`planetgen/web/galaxy_views.py`).
    """
    try:
        conn = get_db()
        return jsonify(galaxy_window(conn, request.args.get("at") or "", request.args.get("g") or DEFAULT_WINDOW,
                                     galaxy_edge_pc(conn)))
    except ValueError as exc:
        raise ApiError(str(exc))


@bp.route("/galaxy/nebula-cover")
def galaxy_nebula_cover_route():
    """
    The share of each cell inside a nebula for a Galaxy Map stage (`?at=m.shell.band.wedge`, or none for the galaxy)
    or window (`?at=` the middle cell and `?g=` the cells across), filled cells and unfilled alike -- see
    `db.query.galaxy_nebula_cover`. Called by `/galaxy/stage` and `/galaxy/window` (`planetgen/web/galaxy_views.py`),
    which add it to their payloads; it is not cached with them because it changes without a sector changing.
    """
    try:
        conn = get_db()
        return jsonify(galaxy_nebula_cover(conn, request.args.get("at") or None, request.args.get("g") or None,
                                           galaxy_edge_pc(conn)))
    except ValueError as exc:
        raise ApiError(str(exc))


@bp.route("/galaxy/charted-regions")
def galaxy_charted_regions_route():
    """
    `{"regions": [...], "edge_pc"}` -- every charted region (a 9-block holding a generated sector, `region_stats`), by share
    charted, each with its node key, the centre and radius of its extent in parsecs, its sector count and capacity, its
    `share` (0 to 1) and, with a naming key in force, its `name` and `parent`, and its `population` once head counts are stored
    (absent before, MAP.318) -- see `db.region_stats.charted_regions`.
    Read by the Galaxy Map (MAP.316) for the regions it draws at galactic scale and its charted-regions list.
    """
    conn = get_db()
    edge_pc = galaxy_edge_pc(conn)
    return jsonify({"regions": region_stats.charted_regions(conn, naming_key.active_key(), edge_pc), "edge_pc": edge_pc})


@bp.route("/galaxy/stamp")
def galaxy_stamp_route():
    """
    `{"stamp": "<16 hex>", "state": "<token>"}` -- `stamp` changes whenever
    the galaxy's node contents could change (see
    `db.query.galaxy_content_stamp`), and `state` is what
    `/galaxy/changes?since=` takes to say which nodes did.
    """
    changes = galaxy_changes(get_db())
    return jsonify({"stamp": changes["stamp"], "state": changes["state"]})


@bp.route("/galaxy/changes")
def galaxy_changes_route():
    """
    `?since=<state>` -- which nodes and stages changed since that state (see
    `db.query.galaxy_changes`): `{"stamp", "tile_stamp", "state", "full", "tile_full", "stages", "nodes"}`, and, unless `full`, `sectors` (the printed IDs of the sectors changed) and
    `grew` (a sector or a system was added), which the page cache (PERF.104) reads.
    `planetgen/web/lib/tilecache.py` calls this about once a minute and deletes
    only the listed nodes, or all of them when `tile_full`. A missing or
    unreadable `since` just answers `full`.
    """
    conn = get_db()
    changes = galaxy_changes(conn, request.args.get("since"))
    if "sector_ids" in changes:
        # The page cache tags a sector by its printed ID (API.23), as the pages and the API speak it.
        printed = ids.printed_many(conn, ids.SECTOR, changes.pop("sector_ids"))
        changes["sectors"] = sorted(printed.values())
    return jsonify(changes)

@bp.route("/uncharted-systems")
def uncharted_systems():
    """
    `GET /api/uncharted-systems` -- the scattered stars no system was built
    around yet (UX.87), brightest first: `db.query.list_uncharted_systems`'s
    rows (the star, its coordinates, its sector address, designation and
    position in the sector). Query parameters: `sort` (`luminosity` -- the
    default -- `temperature`, `type`, `sector` or `age`) with `order`
    (`asc` or `desc`; `desc` for luminosity when omitted, else `asc`), `limit`
    and `offset`. `after=<star id>` starts the page after that star in the sort instead of skipping `offset`
    rows (PERF.75). `from=<value>` (PERF.135, the jump-by-value box) starts the list at the first star at or after
    that value in the sort and its direction: a luminosity in solar luminosities, a temperature in kelvin, an age in
    Gyr, a shell number, or the start of a star type; a value that is not one is a 400. `total` counts every star
    waiting (`total_capped`: it stopped at 10,000 before the first stored count was made, PERF.112), or from the
    start when `from` is given. A value past the last star answers the last page, `total` its length and `past_end`
    true.
    """
    limit, offset = _paginate(request.args)
    sort, descending = _parse_sort(request.args, UNCHARTED_SYSTEM_SORTS)
    if "order" not in request.args:
        descending = sort == "luminosity"
    after = None
    if request.args.get("after"):
        try:
            after = int(request.args["after"])
        except ValueError:
            raise ApiError(f"after must be a star id, got {request.args['after']!r}")
    start = (request.args.get("from") or "").strip()[:MAX_START_LENGTH] or None
    db = get_db()
    try:
        total = count_uncharted_systems(db, start=start, sort=sort, descending=descending)
        items = list_uncharted_systems(db, limit=limit, offset=offset, sort=sort, descending=descending, after=after,
                                       start=start)
        past_end = False
        if start and not items and not offset and after is None and limit is not None:
            past_end = True
            items = list_uncharted_systems(db, limit=limit, sort=sort, descending=descending, start=start, last=True)
            total = len(items)
    except ValueError as exc:
        raise ApiError(str(exc))
    return jsonify({
        "items": items, "total": total, "total_capped": is_capped(total), "limit": limit, "offset": offset,
        "past_end": past_end,
    })


@bp.route("/object-lists/<kind>")
def object_list(kind):
    """
    `GET /api/object-lists/<kind>` -- one keyset page of the objects of a kind from the object index (UX.141): `kind` is
    `system`, `star`, `planet`, `moon`, `belt`, `comet`, `facility` or a phenomenon type. Query parameters: `sort`
    (`name` -- the default --, `class`, `mass` or `radius`; the last three list the objects that have one) with
    `order` (`asc` or `desc`), `limit` (50 by default, 200 at most), `after=<object ID>` (the `next` of the previous
    page: no `offset`), `name` (names starting with it), `class` (exactly this class), `sector_id` and `from` (begin the list at the first
    name from this one on, with `sort=name` only; `past_end` is true when none lies at or past it and the last page
    is answered). Answers `items`
    (`kind`, `id`, `source_id`, `name`, `sector_id`, `system_id`, `parent_id`, `object_class`, `mass_kg`, `radius_km`,
    `placed`), `next` (null at the end), `total` and `total_capped`. 400 for an unknown kind, sort or filter, an `after`
    that is not an object ID, a `from` with another sort, or a `sector_id` that is not a number.
    """
    sort, descending = _parse_sort(request.args, objectlist.SORTS)
    try:
        limit = int(request.args.get("limit") or objectlist.DEFAULT_LIMIT)
        sector_id = int(request.args["sector_id"]) if request.args.get("sector_id") else None
    except ValueError:
        raise ApiError("limit and sector_id must be whole numbers")
    try:
        page = objectlist.list_objects(
            get_db(), kind, limit=limit, sort=sort, descending=descending, after=request.args.get("after") or None,
            name=(request.args.get("name") or "").strip()[:MAX_START_LENGTH] or None,
            object_class=request.args.get("class") or None, sector_id=sector_id,
            start=(request.args.get("from") or "").strip()[:MAX_START_LENGTH] or None)
    except ValueError as exc:
        raise ApiError(str(exc))
    return jsonify({**page, "total_capped": is_capped(page["total"]), "limit": min(limit, objectlist.MAX_LIMIT)})


@bp.route("/uncharted-systems/<int:bright_star_id>")
def uncharted_system(bright_star_id):
    """
    `GET /api/uncharted-systems/<id>` -- one scattered star (UX.114): the same fields as a `GET /api/uncharted-systems`
    item, `system_id` (the system built around it, `null` while it is uncharted) and, for a star with company
    that is not built yet, its `members` and orbit `nodes` (`db.query.uncharted_system_detail`). 404 for an id that
    is no row.
    """
    try:
        return jsonify(query_uncharted_system_detail(get_db(), bright_star_id))
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404


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
            # GEN.315: the system is the tree the row stored, and draws no companion of its own
            plan = stored_trees.plan(scatter_systems.read_stored_trees(conn, "bright_stars", [bright_star_id])
                                     .get(bright_star_id), store.get_orbit_epoch_unix(conn))
            with draw.bound(row["seed"]):
                system = _generate_system(system_config, pc_to_ly(math.sqrt(x * x + y * y + z * z)),
                                          primary_star_params=brightStars.star_params(row), member_tree=plan)
            system_id = store.insert_star_system(conn, system, system_config)
            store.assign_uids(conn, system_ids=[system_id])
            store.mark_bright_star_filled(conn, bright_star_id, system_id)
            if system.tree_built:
                scatter_systems.delete_member_trees(conn, "bright_stars", [bright_star_id])
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
    for the type and descriptor filter menus (`db.query.phenomena_facets`).
    `total` counts the rows that pass the filters (`total_capped`: it stopped at 10,000, PERF.69).
    `after=<scatter id>` (the `scatter_id` of a scattered row) starts the page after that row in the tail of
    unbuilt scattered phenomena instead of skipping `offset` rows (PERF.75). `after=<type>:<ID>` (the `type` and `id`
    of a built row, such as `nebula:0008000000-3000000-000`) starts it after that row in the order of `sort`, then
    the tail (PERF.205, PERF.211: any page costs the same, `offset` is ignored with it; the first page of any sort needs
    no `after` and costs the same too, while a page by `offset` past the first sorts the whole set of built rows). `from=<scatter id>` (PERF.135, the
    jump-by-value box) lists that tail from the row with that id on (the built phenomena come before it and are not
    listed); a value that is not a whole number from 1 up is a 400, and one past the last scattered row answers the last
    page of the tail, `total` its length and `past_end` true.
    """
    limit, offset = _paginate(request.args)
    sort, descending = _parse_sort(request.args, PHENOMENON_SORTS)
    filters = {
        "types": [v for v in request.args.getlist("type") if v],
        "descriptors": [v for v in request.args.getlist("descriptor") if v],
        "placed": _parse_yes_no(request.args.get("placed"), "placed"),
    }
    db = get_db()
    after = None
    if request.args.get("after"):
        type_label, colon, printed = request.args["after"].partition(":")
        if colon:
            row = ids.row_id(db, type_label, printed) if type_label in ids.PHENOMENON_KINDS else None
            if row is None:
                raise ApiError(f"after must be a scatter id or <type>:<ID> of a phenomenon, got {request.args['after']!r}")
            after = (type_label, row)
        else:
            try:
                after = int(request.args["after"])
            except ValueError:
                raise ApiError(f"after must be a scatter id or <type>:<ID> of a phenomenon, got {request.args['after']!r}")
    start = (request.args.get("from") or "").strip()[:MAX_START_LENGTH] or None
    try:
        total = count_phenomena(db, start=start, **filters)
        items = list_phenomena(db, limit=limit, offset=offset, sort=sort, descending=descending, after=after,
                               start=start, **filters)
        past_end = False
        if start and not items and not offset and after is None and limit is not None:
            past_end = True
            items = list_phenomena(db, limit=limit, sort=sort, descending=descending, start=start, last=True,
                                   **filters)
            total = len(items)
    except ValueError as exc:
        raise ApiError(str(exc))
    body = {
        "items": items,
        "total": total,
        "total_capped": is_capped(total),
        "limit": limit,
        "offset": offset,
        "past_end": past_end,
    }
    if request.args.get("facets") == "1":
        body["facets"] = phenomena_facets(db, **filters)
    return jsonify(body)


@bp.route("/phenomena/<phenomenon_type>/<uid:phenomenon_id>")
def phenomenon(phenomenon_type, phenomenon_id):
    """
    One phenomenon's full detail -- `html/phenomenon.py`'s info page.
    `phenomenon_type` is one of `db.query._PHENOMENON_TYPE_TO_TABLE`'s keys
    (`nebula`/`asteroid_field`/`black_hole`/`neutron_star`/
    `supernova_remnant`/`rogue_planet`/`interstellar_comet`); anything
    else, or an id that doesn't exist under it, is a 404.
    """
    try:
        detail = query_phenomenon_detail(get_db(), phenomenon_type, phenomenon_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(detail)


@bp.route("/phenomena/<phenomenon_type>/<uid:phenomenon_id>/contents")
def phenomenon_contents(phenomenon_type, phenomenon_id):
    """
    `GET /api/phenomena/<type>/<id>/contents` -- what sits inside a nebula or supernova remnant (UX.144): `counts` by
    kind, `total`, and up to `limit` (default 100, at most 500) `items` (`kind`, `id`, `name`) in kind, name order. 404
    for another type or an id that is no row.
    """
    conn = get_db()
    try:
        query_phenomenon_detail(conn, phenomenon_type, phenomenon_id)
        limit = min(max(_int_arg(request.args, "limit", cloudcontents.DEFAULT_LIMIT), 1), 500)
        return jsonify(cloudcontents.cloud_contents(conn, phenomenon_type, phenomenon_id, limit))
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404


@bp.route("/containers")
def containers_list():
    """
    `GET /api/containers[?kind=globular_cluster|superbubble]` -- the galaxy's stored containers (DB.85), by kind then
    id: `{"containers": [{"id", "uid", "kind", "name", "center_pc", "extent_pc", "region", "bubble"?}, ...]}`. 400 for
    another kind.
    """
    try:
        return jsonify({"containers": containersDb.listing(get_db(), request.args.get("kind"))})
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


@bp.route("/containers/at")
def containers_at():
    """
    `GET /api/containers/at?x=&y=&z=[&kind=]` -- the containers holding a galaxy-frame point (parsecs), deepest first,
    each with `depth` (distance from its centre over its extent) and `distance_pc`. 400 for a missing coordinate or
    another kind.
    """
    try:
        point = [_float_arg(request.args, axis) for axis in ("x", "y", "z")]
        kinds = [request.args["kind"]] if request.args.get("kind") else None
        return jsonify({"containers": containersDb.holding(get_db(), point, kinds)})
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


@bp.route("/containers/<printed_id>")
def container_detail(printed_id):
    """
    `GET /api/containers/<object ID>` -- one container by its printed object ID. 404 for an ID that is no container.
    """
    try:
        found = containersDb.get(get_db(), uid=printed_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    if found is None:
        return jsonify({"error": f"no container {printed_id}"}), 404
    return jsonify(found)


@bp.route("/containers/<printed_id>/contents")
def container_contents(printed_id):
    """
    `GET /api/containers/<object ID>/contents[?limit=]` -- what sits inside a container: `counts` by kind, `total`, and
    up to `limit` (default 100, at most 500) `items` (`kind`, `id`, `name`) in kind, name order. 404 for an ID that is
    no container.
    """
    conn = get_db()
    try:
        found = containersDb.get(conn, uid=printed_id)
        if found is None:
            raise ValueError(f"no container {printed_id}")
        limit = min(max(_int_arg(request.args, "limit", containersDb.DEFAULT_LIMIT), 1), 500)
        return jsonify(containersDb.contents(conn, found["id"], limit))
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404


@bp.route("/uncharted-phenomena/<int:scatter_id>")
def uncharted_phenomenon(scatter_id):
    """
    `GET /api/uncharted-phenomena/<id>` -- one scattered phenomenon by its `phenomenon_scatter` row (MAP.211):
    what is known of a black hole, neutron star, remnant, planetary nebula, quasar or hypervelocity star
    before its sector is generated, and once it is, the object it became (`db.query.scattered_phenomenon_detail`).
    404 for an id that is no row.
    """
    try:
        return jsonify(query_scattered_phenomenon_detail(get_db(), scatter_id))
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404


@bp.route("/objects/<ref>")
def object_ref_route(ref):
    """
    Resolves one object reference (`<kind>:<id>`, NAV.7): name, parent
    chain, sibling references and positions in each frame -- see
    `db.query.resolve_object` and docs/reference/api.md.
    """
    try:
        resolved = ids.resolve_ref(get_db(), ref)
        if resolved is None:
            raise ValueError(f"nothing has the ID {ref}")
        kind, object_id = object_ref.parse(resolved)
        found = resolve_object(get_db(), kind, object_id)
        found["region"] = region_of(get_db(), found["positions"]["galaxy_pc"])
        return jsonify(found)
    except ValueError as exc:  # includes ids.IdError
        return jsonify({"error": str(exc)}), 404


def _body_of(ref):
    """`(kind, row id)` of the star, planet, moon, belt or comet a public reference names (the sheet and the children
    of a page's object); `ValueError` for anything else, with the message the 404 carries."""
    resolved = ids.resolve_ref(get_db(), ref)
    if resolved is None:
        raise ValueError(f"nothing has the ID {ref}")
    return object_ref.parse(resolved)


@bp.route("/objects/<ref>/sheet")
def object_sheet_route(ref):
    """
    `GET /api/objects/<ref>/sheet` (UX.168, UX.175) -- every stored property of a star, planet, moon, belt, comet or
    phenomenon, grouped
    for its page (`describe/sheet.py`): `kind` and `groups`, each `{key, label, rows}` with a row
    `{key, label, value (the Document's inline list: quantities with their units), text, stored}`. A property the row
    holds no value for is `stored: false` and says "not stored". 404 for a reference that names nothing or a kind with
    no sheet.
    """
    try:
        kind, object_id = _body_of(ref)
        return jsonify(object_sheet(get_db(), kind, object_id))
    except ValueError as exc:  # includes ids.IdError
        return jsonify({"error": str(exc)}), 404


@bp.route("/objects/<ref>/children")
def object_children_route(ref):
    """
    `GET /api/objects/<ref>/children[?limit=&after=]` (UX.169) -- what belongs to an object, one keyset page at a time
    (`db/children.py`): a sector's systems then its phenomena, a system's stars then its circumbinary bodies, a star's
    planets and belts then its comets, a planet's moons. `ref`, `kind`, `total`, `items` (each `ref`, `kind`, `id`,
    `name`, `summary`, `counts` of everything below it by kind, `child_count`) and `next`, the `after` of the page
    after (`null` at the end). `limit` is 1 to 200 (default 50). A moon, belt, comet or phenomenon holds nothing.
    404 for a reference that names nothing, 400 for a bad `limit` or `after`.
    """
    try:
        resolved = ids.resolve_ref(get_db(), ref)
        if resolved is None:
            raise ValueError(f"nothing has the ID {ref}")
        kind, object_id = object_ref.parse(resolved)
        limit = _int_arg(request.args, "limit", children_db.DEFAULT_LIMIT)
        return jsonify(children_db.children(get_db(), kind, object_id, limit, request.args.get("after") or None))
    except children_db.PageError as exc:
        return jsonify({"error": str(exc)}), 400
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
    `db.query.search`'s own docstring for the full request/response shape.
    `html/search.py` is a thin renderer over this endpoint's JSON.

    Query parameters: `sector_q`/`system_q`/`star_q`/`planet_q`/`moon_q`
    (name search terms); every name in `db.query.SEARCH_TAG_FACETS`
    (repeatable, e.g. `?class=M&class=K` for two active Class tags); and
    `star_min_radius_km`/`star_max_radius_km` (likewise `planet_`/
    `moon_`) for a size range, in km -- either bound may be omitted.

    UX.214: the score filters, each repeatable and named for a planet (`habitability`, `domain`, `life_level`, `living`)
    or with `moon_` in front for a moon -- `habitability=blue|green|yellow|red` (the PHI-6 colour), `domain=<domain>:<colour>`
    (`suitability:green`, `water:red`; the domains are pressure, temperature, air, water, radiation and suitability),
    `life_level=<microbial|flora|fauna|sentient>:<1|5>` (that life score 1 or more, or 5 or more, of 9) and
    `living=alive|lifeless|food|civilization`. `planet_sort`/`moon_sort` order that panel by `phi4`, `suitability`, a domain
    name, `microbial`, `flora`, `fauna` or `sentient`: worst first, or `-` in front for best first (an unscored body
    comes last either way).

    Each result panel is paged independently: `limit` is the rows per
    panel (default `SEARCH_RESULT_LIMIT`, clamped to `MAX_PAGE_LIMIT`,
    validated like every other `limit`), and `sectors_offset`/
    `systems_offset`/`stars_offset`/`planets_offset`/`moons_offset`/
    `belts_offset` pick each panel's page. `panels` (comma-separated panel
    names; `none` for no panel at all, which still answers the facets and name suggestions) runs just those
    panels, the rest answering `null`; the Search page fetches each panel on its own so one slow panel
    times out alone (UX.97).

    UX.226: numeric range filters and sorts on the stored physical values (`db.search_values.KEYS`, listed with their units in the
    response's `values`). `<entity>_min_<key>` and `<entity>_max_<key>` (`planet_min_gravity=0.8&planet_max_gravity=1.2`,
    `planet_max_surface_dose=2`, `star_min_luminosity=0.5`) keep the bodies whose stored value is in the range, bounds in the key's
    default unit or in the unit `<entity>_<key>_unit` names (`planet_surface_pressure_unit=kPa`; it must measure the same kind of
    thing); a body with no value for the key is out of the range. The keys: gravity (g0), surface_temperature (K), surface_pressure
    (bar), orbit_distance (au for a planet, km for a moon), day_length and rotation_period (h), year_length (yr), mass (Earth masses
    for a planet or a moon, solar masses for a star), surface_dose (mSv/yr), magnetic_moment (A*m2), ocean_depth (km), moon_count
    (a planet), gas_<formula> (the partial pressure of a gas, bar; no row is 0), and for a star luminosity (Lsun), age (Gyr) and
    temperature (K). `star_sort`, `planet_sort` and `moon_sort` also take a key: smallest first, or `-` in front for largest first, a
    body with no value last. Rows then carry `values` (`{key: {value, unit, symbol}}`) for the keys in use, and `active_ranges`
    lists the filters for chips. A bound that is not a number, a minimum over a maximum or a unit of another kind is a 400.

    UX.239: `stars` keeps the systems with that many stars: a count (`stars=3`, exactly three), a count and a plus (`stars=3+`, three or more)
    or a name (`single`, `binary`, `trinary`, `quaternary`: exactly that many; `quintenary`: five or more). The Systems panel runs for it with
    or without a name; each row carries `star_count`. Anything else is a 400.
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
        only = {panel for panel in args["panels"].split(",") if panel in SEARCH_RESULT_PANELS}  # "none": no panel
    sorts = {entity: args.get(f"{entity}_sort") for entity in ("star", "planet", "moon") if args.get(f"{entity}_sort")}
    stars = _parse_stars(args)
    try:
        ranges = search_values.ranges_from_args(args)
    except ValueError as exc:
        raise ApiError(str(exc))
    return jsonify(run_search(get_db(), texts, tags, sizes=sizes, limit=limit, offsets=offsets, only=only, sorts=sorts, stars=stars,
                              ranges=ranges))


# ---------------------------------------------------------------------
# Write endpoints. Every one requires an authenticated admin
# (`authz.require_admin`) whose credentials aren't still the seeded
# default (`fresh=True` -- see `authz.py`/`admin/auth.py`), runs against
# `WRITE_MYSQL_CONFIG` (a separate, less-privileged-than-`CREATE`
# account -- see `config.py`) rather than the request-scoped read-only
# `get_db()` every read route above uses, and records one
# `admin_audit_log` row (`authz.audit`) after it actually succeeds, never
# speculatively. See docs/reference/api.md's "Write endpoints" section for the
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
    same reasoning as `db.query.open_readonly`/`store.open_write`: the
    write-capable account has no `CREATE` grant (see `config.py`), so the
    schema must already exist."""
    return store.open_write(_resolve_requested_write_db_config())


@bp.route("/sectors", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def create_sector():
    """`POST /api/sectors` `{"name": str, "edge_ly": number > 0}` (both
    required) -- inserts a new, empty `sectors` row (no systems -- use
    `generation/run_sector.py`/`generation/run_galaxy.py` to generate a populated one) and
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
    "edge_ly": number > 0}`."""
    body = given(parse_body(SectorUpdate, require_json_body()))

    set_clauses, params = [], []
    if "name" in body:
        set_clauses.append("name = ?")
        params.append(body["name"])
    if "edge_ly" in body:
        set_clauses.append("edge_mpc = ?")
        params.append(ly_to_milliparsecs(body["edge_ly"]))
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
            conn.execute(f"UPDATE sectors SET {', '.join(set_clauses)}, parity_dirty = 1 WHERE id = ?", tuple(params))
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
            address = conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors WHERE id = ?",
                                   (sector_id,)).fetchone()
            deleted = conn.execute("DELETE FROM sectors WHERE id = ?", (sector_id,)).rowcount > 0
            if deleted:
                store.note_sector_deleted(conn)
            if address is not None:
                # GEN.44: the cell goes back to the bright-star level it had unfilled.
                store.forget_sector_fill(conn, (address["shell_index"], address["band_index"],
                                              address["wedge_index"]))
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
    genuinely large** (100 ly is ~2,000-3,000 candidate sector cells), so
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
    507 and nothing written; a radius holding more than 2 million sectors
    (`run_common.WEB_ENUMERATION_LIMIT`, PERF.36) is refused with 413."""
    body = request_json()
    if body is None:
        body = {}
    if not isinstance(body, dict):
        raise ApiError("request body must be a JSON object")
    request_body = parse_body(NeighborhoodRequest, body)
    radius_ly, estimate_only = request_body.radius_ly, request_body.estimate_only
    disregard = request_body.disregard_missing_data

    config = _resolve_requested_write_db_config()
    if not estimate_only:
        try:
            run_galaxy.require_datasets(disregard)
        except RuntimeError as exc:
            raise ApiError(str(exc), status_code=409)
    try:
        # A real run is estimated first, so an unknown sector (404), a
        # missing skeleton (409) and a database disk that can't hold it
        # (507) answer here, before anything is queued.
        result = run_galaxy.generate_sector_neighborhood(
            sector_id, radius_ly=radius_ly, config=config, estimate_only=True,
        )
    except run_galaxy.EnumerationTooLarge as exc:
        raise ApiError(str(exc), status_code=413)
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
        job_id = api_jobs.submit(api_jobs.generate_neighborhood, sector_id, radius_ly, config, disregard,
                                 expected_seconds=result["estimate"].get("seconds"))
    except api_jobs.NoQueue as exc:
        raise ApiError(str(exc), status_code=503)
    audit(
        "sector.generate_neighborhood", target=f"sector:{sector_id}",
        detail=f"radius_ly={radius_ly!r} job={job_id}" + (" disregard_missing_data" if disregard else ""),
    )
    return accepted(job_id)


@bp.route("/galaxy/blocks/<key>/generate", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def generate_block_route(key):
    """`POST /api/galaxy/blocks/<key>/generate` (GEN.310) -- generates every not-yet-generated sector the stored outline
    holds inside the block of the Galaxy Map's ladder `key` names (`m.shell.band.wedge` of a stage level: 243, 81,
    27, 9 or 3), shell by shell outward, with the bright-star backfill and the settle after them -- see
    `run_galaxy.generate_block`. A queued job like generate-neighborhood: `202` with the job's `status_url`; the job
    reports its progress there and in the queue. `"estimate_only": true` returns the counts and the size and time
    `estimate` without writing anything. A key that names no block is refused with 400, a block holding more than
    2 million sectors (`run_common.WEB_ENUMERATION_LIMIT`) with 413, a galaxy never planned with 409, and a run the
    database disk can't hold with 507."""
    body = request_json()
    if body is None:
        body = {}
    if not isinstance(body, dict):
        raise ApiError("request body must be a JSON object")
    request_body = parse_body(BlockGenerateRequest, body)
    estimate_only, disregard = request_body.estimate_only, request_body.disregard_missing_data
    try:
        block = drill.parse_stage_key(key)
    except ValueError as exc:
        raise ApiError(str(exc), status_code=400)
    key = drill.format_drill_key(block)

    config = _resolve_requested_write_db_config()
    if not estimate_only:
        try:
            run_galaxy.require_datasets(disregard)
        except RuntimeError as exc:
            raise ApiError(str(exc), status_code=409)
    try:
        result = run_galaxy.generate_block(key, config=config, estimate_only=True)
    except run_galaxy.EnumerationTooLarge as exc:
        raise ApiError(str(exc), status_code=413)
    except ValueError as exc:
        raise ApiError(str(exc), status_code=400)
    except RuntimeError as exc:
        raise ApiError(str(exc), status_code=409)
    if estimate_only:
        return jsonify(result)
    if result["estimate"].get("refusal"):
        raise ApiError(result["estimate"]["refusal"], status_code=507)

    try:
        job_id = api_jobs.submit(api_jobs.generate_block, key, config, disregard,
                                 expected_seconds=result["estimate"].get("seconds"))
    except api_jobs.NoQueue as exc:
        raise ApiError(str(exc), status_code=503)
    audit("galaxy.generate_block", target=f"block:{key}",
          detail=f"job={job_id}" + (" disregard_missing_data" if disregard else ""))
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


@bp.route("/table-estimates")
def table_estimates_route():
    """`GET /api/table-estimates` -- `{"tables": {name: rows}}`, the database's own estimate (from
    `information_schema`, instant and within some tens of percent) of every table's row count. The web pages read it
    once a minute to decide which requests are slow (PERF.105, `web/lib/slowpages.py`)."""
    conn = get_db()
    return jsonify({"tables": {row["name"]: row["approx_rows"] for row in dbStats.table_stats(conn)}})


@bp.route("/page-jobs/<job_id>")
def page_job_status_route(job_id):
    """`GET /api/page-jobs/<id>` -- where the job building a slow page's answer stands (PERF.105): `state`
    (`queued`, `running`, `succeeded`, `failed` or `cancelled`), `done` of `total` rows (an estimate: `estimated`
    is true), `unit`, `eta_seconds`, `seconds`, `error`, `cancelling`, `cache_key` and `target`. Public, like the
    pages it serves; it names only page jobs, never the admin's other jobs. 404 for an unknown or expired id
    (a job is kept an hour); 503 without Redis."""
    try:
        job = pagejobs.status(job_id)
    except pagejobs.NoQueue as exc:
        raise ApiError(str(exc), status_code=503)
    if job is None:
        raise ApiError(f"no such page job: {job_id}", status_code=404)
    return jsonify(job)


@bp.route("/page-jobs/<job_id>/cancel", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def page_job_cancel_route(job_id):
    """`POST /api/page-jobs/<id>/cancel` -- admin only: stops the job building a slow page's answer if it has not
    started its query (a query already running finishes, and its answer is kept in the cache); answers `GET`'s body.
    Every visitor waiting on the same request shares the job, so a visitor's own "Cancel" only stops that page
    waiting (it stops polling) and the job goes on for the others. 404 and 503 as `GET`."""
    try:
        job = pagejobs.cancel(job_id)
    except pagejobs.NoQueue as exc:
        raise ApiError(str(exc), status_code=503)
    if job is None:
        raise ApiError(f"no such page job: {job_id}", status_code=404)
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
        address = (placement["shell_index"], placement["band_index"], placement["wedge_index"])
        level = store.bright_star_fill_level(conn, *address)
        mass_limit = store.bright_star_fill_mass_limit(conn, *address)
        center = (placement["center_x_pc"], placement["center_y_pc"], placement["center_z_pc"])
        profile = profile_at(center, skeleton.shape, skeleton.edge_pc, skeleton.expected_system_count_at_density_1,
                             address=address)
        brightStars.FillContext(profile, min_luminosity_sol=level, star_mass_limit_sol=mass_limit).apply(system_config)
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
    `SYSTEM_CONFIG_ALLOWED_FIELDS`), the same shape `generation/run_system.py
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
                store.refresh_edited_sectors(conn, [row["sector_id"]])    # DB.83: the sector's figures and the blocks above
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
    """`GET /api/galaxy/bright-stars?shell=&band=&wedge=[&all=1]` -- the
    pre-placed bright stars in one sector cell, brightest first
    (`db.query.bright_stars_in_sector`): only those not yet built into a
    system unless `all=1`. `{"items": []}` when no scatter has run (or
    the database predates them)."""
    try:
        shell, band, wedge = (int(request.args[key]) for key in ("shell", "band", "wedge"))
    except (KeyError, ValueError):
        raise ApiError("'shell', 'band' and 'wedge' must be integers")
    unfilled_only = request.args.get("all") not in ("1", "true")
    try:
        items = bright_stars_in_sector(get_db(), shell, band, wedge, unfilled_only=unfilled_only)
    except pymysql.err.ProgrammingError:  # no bright_stars table yet (before v43)
        items = []
    return jsonify({"items": items})


@bp.route("/galaxy/uncharted")
def galaxy_uncharted_sector():
    """`GET /api/galaxy/uncharted?shell=&band=&wedge=` -- one sector cell's
    place and what the scatters left in it, for opening a cell nothing was
    generated in (MAP.162): `designation`, `center_pc`, `edge_pc`,
    `sector_id` (the generated sector there, or `null`), `stars` (its
    waiting bright stars, `db.query.bright_stars_in_sector`) and `scattered`
    (its unbuilt scattered phenomena, `db.query.uncharted_sector_contents`).
    A cell outside the grid is a 400."""
    try:
        shell, band, wedge = (int(request.args[key]) for key in ("shell", "band", "wedge"))
    except (KeyError, ValueError):
        raise ApiError("'shell', 'band' and 'wedge' must be integers")
    if shell < 0:
        raise ApiError("shell must be >= 0")
    conn = get_db()
    skeleton = get_galaxy_shape(conn)
    edge_pc = skeleton.edge_pc if skeleton else float(tuning.DEFAULT_SECTOR_EDGE_PC)
    try:
        cell = describe_sector_cell(shell, band, wedge, edge_pc)
    except OverflowError:
        raise ApiError("that cell is too far out")
    except ValueError as err:
        raise ApiError(str(err))
    return jsonify({
        "shell_index": shell, "band_index": band, "wedge_index": wedge, "designation": cell["designation"],
        "center_pc": cell["cartesian_pc"], "edge_pc": edge_pc, "sector_id": get_sector_id_at(conn, shell, band, wedge),
        **uncharted_sector_contents(conn, shell, band, wedge),
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


def _facility_attrs(body):
    """The keys of a create or change body that `facility_attrs.resolve` reads."""
    attrs = {key: body[key] for key in ("type_id", "affiliation_polity_id", "affiliation_name", "ratings_note",
                                        "fields") if key in body}
    if body.get("ratings") is not None:
        attrs["ratings"] = body["ratings"]
    return attrs


@bp.route("/facilities")
def facility_list():
    """`GET /api/facilities` -- a page of facilities. Filters: `kind`, `type_id`, `affiliation_polity_id`,
    `affiliation_name`, `rating_crime`, `rating_housing`, `rating_resources`, `rating_maintenance`,
    `rating_health` (each `green`, `yellow` or `red`) and `overall` (the worst of the five). `sort` is `name`
    (the default), `kind`, `type` or `worst` (worst overall rating first, unassessed last; `order=desc` reverses
    the ratings); `limit`/`offset` page. `total` counts what passes the filters."""
    args = request.args
    kind = args.get("kind") or None
    if kind is not None and kind not in tuning.FACILITY_KINDS:
        raise ApiError(f"kind must be one of: {', '.join(sorted(tuning.FACILITY_KINDS))}")

    def whole(name):
        raw = args.get(name)
        if raw in (None, ""):
            return None
        try:
            return int(raw)
        except ValueError:
            raise ApiError(f"{name} must be an integer, got {raw!r}") from None

    def color(name):
        raw = args.get(name)
        if raw in (None, ""):
            return None
        if raw not in facility_ratings.COLORS:
            raise ApiError(f"{name} must be one of: {', '.join(facility_ratings.COLORS)}")
        return raw

    ratings = {dimension: color(f"rating_{dimension}") for dimension in facility_ratings.RATINGS}
    sort, descending = _parse_sort(args, FACILITY_SORTS)
    limit, offset = _paginate(args)
    items, total = list_facilities(
        get_db(), kind=kind, type_id=whole("type_id"), polity_id=whole("affiliation_polity_id"),
        affiliation_name=args.get("affiliation_name") or None,
        ratings={name: value for name, value in ratings.items() if value is not None}, overall=color("overall"),
        sort=sort, descending=descending, limit=limit, offset=offset)
    return jsonify({"items": items, "total": total, "limit": limit, "offset": offset})


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
                attrs=_facility_attrs(body),
            )
    except store.FacilityError as exc:
        raise ApiError(str(exc), status_code=404 if exc.not_found else 400)
    except facility_types.FacilityTypeError as exc:
        raise ApiError(str(exc), status_code=exc.status)
    finally:
        conn.close()

    audit("facility.create", target=f"facility:{facility_id}",
          detail=str({key: body[key] for key in ("name", "kind", "placement", "host_type", "host_id")}))
    return jsonify({"id": facility_id}), 201


@bp.route("/facilities/<uid:facility_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def update_facility(facility_id):
    """`PATCH /api/facilities/<id>` -- any of `name`, `description`, `type_id`, `affiliation_polity_id` or
    `affiliation_name`, `ratings`, `ratings_note` and `fields` (see `POST /api/facilities`); a key sent as `null`
    clears it. The place never changes. Returns the facility."""
    body = given(parse_body(FacilityUpdate, require_json_body()))
    conn = _write_conn()
    try:
        with conn:
            found = store.update_facility(conn, facility_id, name=body.get("name", "").strip() or None,
                                          description=body.get("description"), attrs=_facility_attrs(body))
            detail = facility_detail(conn, facility_id) if found else None   # not `get_db()`: its snapshot predates the write
    except facility_types.FacilityTypeError as exc:
        raise ApiError(str(exc), status_code=exc.status)
    finally:
        conn.close()
    if not found:
        raise ApiError(f"no such facility: {facility_id}", status_code=404)
    audit("facility.update", target=f"facility:{facility_id}", detail=str(sorted(body)))
    return jsonify(detail)


@bp.route("/facilities/<uid:facility_id>/suggest-ratings", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def suggest_facility_ratings(facility_id):
    """`POST /api/facilities/<id>/suggest-ratings` (no body) -- fills the ratings nobody set from the stored
    data (`facility_suggest`), marks them `generated` (unless some were set by hand) and never changes a rating
    that is set. Returns `{"filled": [dimensions], "facility": {...}}`."""
    conn = _write_conn()
    try:
        with conn:
            filled = facility_suggest.suggest_for_facility(conn, facility_id)
            detail = facility_detail(conn, facility_id) if filled is not None else None
    finally:
        conn.close()
    if filled is None:
        raise ApiError(f"no such facility: {facility_id}", status_code=404)
    audit("facility.suggest_ratings", target=f"facility:{facility_id}", detail=", ".join(filled) or "none")
    return jsonify({"filled": filled, "facility": detail})


@bp.route("/facilities/suggest-ratings", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def suggest_all_facility_ratings():
    """`POST /api/facilities/suggest-ratings` (no body) -- the same for every facility with an unset rating (at
    most `facility_suggest.BATCH_LIMIT`, lowest ids first). Returns `{"facilities", "ratings"}`, how many facilities
    changed and how many ratings were filled."""
    conn = _write_conn()
    try:
        with conn:
            facilities, filled = facility_suggest.suggest_for_unassessed(conn)
    finally:
        conn.close()
    audit("facility.suggest_ratings", detail=f"{facilities} facilities, {filled} ratings")
    return jsonify({"facilities": facilities, "ratings": filled})


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
