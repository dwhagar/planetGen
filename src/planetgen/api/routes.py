# planetgen/api/routes.py

"""
JSON endpoints over the planetGen database.

The read endpoints each open a strictly-read-only connection
(`queryDb.open_readonly` -- same `file:...?mode=ro`-equivalent guarantee
the `planetgen.db.query` CLI already relies on, see that module's docstring for
how "read-only" is enforced by the configured account's grants). Listing
endpoints delegate straight to `planetgen.db.query`'s existing `list_sectors`/
`list_systems`/`systems_within_radius` functions rather than
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
import os
import sys

import pymysql
from flask import Blueprint, current_app, g, jsonify, request

# generate.py lives at the repo root, three levels above src/planetgen/api/ (this
# file) -- src/ itself is already on sys.path (see html/wsgi.py's own
# docstring), but the repo root isn't, so it's added here specifically for
# this import. Only `generate_sector_neighborhood_route` below needs it.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import generate  # noqa: E402

from planetgen.db.query import (
    facilities_for_system,
    facilities_in_sector,
    facility_detail,
    NO_SECTOR,
    bright_star_scatter_status,
    bright_stars_in_sector,
    NavUnavailable,
    SEARCH_RESULT_LIMIT,
    SEARCH_RESULT_PANELS,
    SEARCH_TAG_FACETS,
    count_phenomena,
    count_sectors,
    count_systems,
    MAX_TILES_PER_REQUEST,
    galaxy_changes,
    galaxy_density_shape,
    galaxy_placed_phenomena,
    galaxy_placed_sectors,
    galaxy_locate,
    galaxy_stage,
    galaxy_tiles,
    list_phenomena,
    list_sectors,
    list_systems,
    nav_between,
    open_readonly,
    phenomenon_detail as query_phenomenon_detail,
    search as run_search,
    sector_detail as query_sector_detail,
    system_detail as query_system_detail,
    systems_within_radius,
)
from planetgen.db import store
from planetgen.generation import bright_stars as brightStars, limits as generationLimits
from planetgen import tuning
from planetgen.population import facilities as facility_rules
from planetgen.db.store import MySQLConfig, get_galaxy_bounds, get_galaxy_shape, get_sector_id_at, list_databases, resolve_database
from planetgen.util.appconfig import load_config
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.geometry import describe_sector_cell, sector_address_at
from planetgen.generation.system import StarSystem
from planetgen.db.render import FORMATS as SYSTEM_TEXT_FORMATS
from planetgen.db.render import render_system_sections, render_system_text
from stellarObjects.utils import format_distance_ly, ly_to_milliparsecs, ly_to_pc, pc_to_ly
from wikiClient import WikiClient, WikiClientAuthError, WikiClientPageExistsError, WikiClientRequestError

from .authz import audit, require_admin
from .common import ApiError, is_http_url, require_json_body
from .limiter import limiter, page_limit

bp = Blueprint("api", __name__, url_prefix="/api")

DEFAULT_PAGE_LIMIT = 100
MAX_PAGE_LIMIT = 500
MAX_PAGE_OFFSET = 2 ** 63 - 1
"""int: The largest `offset` MySQL's LIMIT/OFFSET takes as a plain
integer; past it the query itself fails, so it's a 400 instead."""
MAX_NAME_LENGTH = 255
"""int: `sectors.name`/`star_systems.name` are VARCHAR(255)."""
MAX_WIKI_URL_LENGTH = 2048
"""int: `sectors.wiki_url` is VARCHAR(2048)."""
MAX_SECTOR_EDGE_LY = 1e9
"""float: An upper bound on a sector's `edge_ly` -- far beyond any real
sector, well short of overflowing the unit conversion."""
MAX_NEIGHBORHOOD_RADIUS_LY = generationLimits.MAX_GENERATE_RADIUS_LY
"""float: An upper bound on generate-neighborhood's `radius_ly` (about
652 ly, the same 200 pc cap the Generate page and `generate.py` use), so
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
        seconds = load_config()["mysql"].get("statement_timeout_seconds") or 0
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
    base_config = current_app.config["MYSQL_CONFIG"]
    entries = list_databases(base_config)
    items = []
    for entry in entries:
        conn = None
        try:
            # Opened inside the try: a listed schema can still fail to open
            # (dropped since the listing, no grant), and `open_readonly`
            # signals that with SystemExit, which would otherwise escape
            # Flask altogether.
            conn = open_readonly(MySQLConfig(
                host=base_config.host, port=base_config.port,
                user=base_config.user, password=base_config.password, database=entry["name"],
            ))
            sector_count = conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"]
            system_count = conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]
        except (Exception, SystemExit):
            # A schema matching the configured prefix but missing this
            # project's own tables (e.g. mid-migration, or a stray
            # unrelated database sharing the prefix) shouldn't take down
            # the whole listing -- report it with unknown counts instead.
            sector_count = system_count = None
        finally:
            if conn is not None:
                conn.close()
        items.append({**entry, "sector_count": sector_count, "system_count": system_count})
    return jsonify({"items": items})


@bp.route("/sectors")
def sectors():
    limit, offset = _paginate(request.args)
    db = get_db()
    return jsonify({
        "items": list_sectors(db, limit=limit, offset=offset),
        "total": count_sectors(db),
        "limit": limit,
        "offset": offset,
    })


@bp.route("/sectors/<int:sector_id>")
def sector_detail(sector_id):
    try:
        detail = query_sector_detail(get_db(), sector_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(detail)


def _parse_sector_id_filter(raw_sector_id):
    """
    Parses the `sector_id` query parameter `/api/systems` accepts: an
    integer (one sector), the literal `"none"` (standalone systems --
    `queryDb.NO_SECTOR`, `html/browse.py`'s own table of systems generated
    with no sector), or absent (no filter).

    Raises:
        ApiError: If present and neither `"none"` nor a valid integer.
    """
    if raw_sector_id is None:
        return None
    if raw_sector_id.strip().lower() == "none":
        return NO_SECTOR
    try:
        return int(raw_sector_id)
    except ValueError:
        raise ApiError(f"sector_id must be an integer or 'none', got {raw_sector_id!r}")


@bp.route("/systems")
def systems():
    star_type = request.args.get("star_type")
    sector_id = _parse_sector_id_filter(request.args.get("sector_id"))

    limit, offset = _paginate(request.args)
    db = get_db()
    rows = list_systems(db, star_type_prefix=star_type, sector_id=sector_id, limit=limit, offset=offset)
    return jsonify({
        "items": rows,
        "total": count_systems(db, star_type_prefix=star_type, sector_id=sector_id),
        "limit": limit,
        "offset": offset,
    })


@bp.route("/systems/<int:system_id>")
def system_detail(system_id):
    try:
        detail = query_system_detail(get_db(), system_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(detail)


@bp.route("/systems/<int:system_id>/text")
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


@bp.route("/systems/<int:system_id>/sections")
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


@bp.route("/systems/<int:system_id>/near")
def systems_near(system_id):
    # This endpoint itself still returns bare JSON, just a list of ids and
    # distances -- but "no way to actually see the two points" is covered
    # now by /api/nav's web page (`html/nav.py`), which renders the
    # `navmap.py` top-down plot on top of the same course/route data.
    raw_radius = request.args.get("radius")
    if raw_radius is None:
        raise ApiError("radius query parameter is required")
    try:
        radius = float(raw_radius)
    except ValueError:
        raise ApiError(f"radius must be a number, got {raw_radius!r}")
    if not math.isfinite(radius) or radius <= 0:
        raise ApiError("radius must be a finite number greater than 0")

    try:
        matches = systems_within_radius(get_db(), system_id, radius)
    except SystemExit as exc:
        # systems_within_radius is shared with the planetgen.db.query CLI and raises
        # SystemExit (its CLI-appropriate error signal) for a missing/
        # unplaced system id -- caught here rather than changing its shared
        # behavior just for this one caller.
        return jsonify({"error": str(exc)}), 404
    return jsonify(matches)


def _parse_system_id_param(query_args, name):
    """
    Parses a required `<name>` query parameter as an id -- a
    `star_systems.id` when its matching `<name>_kind` is `"system"`
    (the default), or a phenomenon table's own row id when it's
    `"phenomenon"` (see `_parse_nav_endpoint_params`).

    Args:
        query_args (werkzeug.datastructures.MultiDict): `request.args`.
        name (str): The query parameter's name (`"from"` or `"to"`).

    Returns:
        int: The parsed id.

    Raises:
        ApiError: If the parameter is missing or not an integer.
    """
    raw_value = query_args.get(name)
    if raw_value is None:
        raise ApiError(f"{name} query parameter is required")
    try:
        return int(raw_value)
    except ValueError:
        raise ApiError(f"{name} must be an integer, got {raw_value!r}")


def _parse_nav_endpoint_params(query_args, name):
    """
    Parses one `/api/nav` endpoint's full reference: its id
    (`_parse_system_id_param`) plus `<name>_kind` (`"system"`, the
    default, or `"phenomenon"`) and, when it's a phenomenon,
    `<name>_type` (one of `queryDb._PHENOMENON_TYPE_TO_TABLE`'s keys --
    validated by `nav_between`/`_load_nav_phenomenon_endpoint` itself,
    via the same `ValueError` -> 404 handling `phenomenon()` above
    already uses for the same set of types, not re-validated here).

    Args:
        query_args (werkzeug.datastructures.MultiDict): `request.args`.
        name (str): `"from"` or `"to"`.

    Returns:
        tuple: `(id, kind, phenomenon_type)` -- `phenomenon_type` is
              `None` when `kind == "system"`.

    Raises:
        ApiError: If the id is missing/not an integer, `<name>_kind` is
            neither `"system"` nor `"phenomenon"`, or `kind ==
            "phenomenon"` with no `<name>_type` given.
    """
    endpoint_id = _parse_system_id_param(query_args, name)
    kind = query_args.get(f"{name}_kind", "system")
    if kind not in ("system", "phenomenon"):
        raise ApiError(f"{name}_kind must be 'system' or 'phenomenon', got {kind!r}")
    phenomenon_type = None
    if kind == "phenomenon":
        phenomenon_type = query_args.get(f"{name}_type")
        if not phenomenon_type:
            raise ApiError(f"{name}_type query parameter is required when {name}_kind is 'phenomenon'")
    return endpoint_id, kind, phenomenon_type


@bp.route("/nav")
def nav():
    """
    Course, distance, and optimal route between two endpoints -- each
    either a star system (the default) or a standalone phenomenon
    (`?from_kind=phenomenon&from_type=nebula&from=<id>`, and likewise for
    `to`) -- see `queryDb.nav_between` for the full availability rules
    and docs/api.md for the response shape.
    """
    from_id, from_kind, from_type = _parse_nav_endpoint_params(request.args, "from")
    to_id, to_kind, to_type = _parse_nav_endpoint_params(request.args, "to")

    try:
        result = nav_between(
            get_db(), from_id, to_id,
            from_kind=from_kind, to_kind=to_kind, from_type=from_type, to_type=to_type,
        )
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    except NavUnavailable as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify({
        "scope": result["scope"],
        "direct": result["direct"]._asdict(),
        "warp_times": [leg._asdict() for leg in result["warp_times"]],
        "fold_times": [leg._asdict() for leg in result["fold_times"]],
        "origin_position": result["origin_position"],
        "destination_position": result["destination_position"],
        "route": _route_for_json(result["route"]),
    })


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
    The galaxy's stored density-skeleton shape (`generate.py plan`'s own
    output, `queryDb.galaxy_density_shape`) -- the real spiral/disk/bulge
    model the Galaxy Map (`/galaxy`) shades its "expected density"
    cloud from, for whatever space hasn't actually been generated yet.
    `"shape"` is `null` when the skeleton has never been built (the map
    then falls back to its own generic illustrative gradient).
    `"bright_stars"` is `queryDb.bright_star_scatter_status`: whether the
    bright-star scatter has run, its threshold and seed, and the default
    threshold a plain plan uses.
    """
    conn = get_db()
    return jsonify({"shape": galaxy_density_shape(conn), "bright_stars": bright_star_scatter_status(conn)})


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

@bp.route("/phenomena")
def phenomena():
    """
    Every exotic phenomenon (nebula/asteroid field/black hole/neutron
    star/supernova remnant/rogue planet/interstellar comet), across every
    sector and regardless of galaxy placement -- the flat, paginated
    counterpart to `/api/galaxy/phenomena` (which only returns the
    galaxy-placed subset, for the Galaxy Map). `html/phenomena.py`'s own
    listing page.
    """
    limit, offset = _paginate(request.args)
    db = get_db()
    return jsonify({
        "items": list_phenomena(db, limit=limit, offset=offset),
        "total": count_phenomena(db),
        "limit": limit,
        "offset": offset,
    })


@bp.route("/phenomena/<phenomenon_type>/<int:phenomenon_id>")
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
    `belts_offset` pick each panel's page.
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

    return jsonify(run_search(get_db(), texts, tags, sizes=sizes, limit=limit, offsets=offsets))


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

SECTOR_FIELDS = {
    # field name -> (expected Python type(s), validator)
    "name": (str, lambda v: bool(v.strip()) and len(v) <= MAX_NAME_LENGTH),
    # `0 < v <= MAX`, not `v > 0`: NaN, infinity and absurd sizes (which
    # overflow the unit conversion) are all refused.
    "edge_ly": ((int, float), lambda v: 0 < v <= MAX_SECTOR_EDGE_LY),
}
"""dict: The `sectors` JSON object shape both `create_sector` and
`update_sector` validate against -- see docs/api.md's "Sectors" write
schema. Shared so the two can never validate the same field two
different ways."""

SECTOR_UPDATE_FIELDS = {
    **SECTOR_FIELDS,
    # `None` clears a manually-set/uploaded link back to "no page yet" --
    # see schema.sql's "v22" header note. `create_sector` deliberately
    # doesn't accept this (a brand-new, just-generated sector has never
    # been uploaded anywhere), so it's added only to `update_sector`'s own
    # allowed-fields set, not to SECTOR_FIELDS itself.
    # Only an absolute http/https URL with a host: the sector page links
    # to it, so `javascript:`/`data:` URLs are refused (`is_http_url`).
    "wiki_url": ((str, type(None)), lambda v: v is None or (len(v) <= MAX_WIKI_URL_LENGTH and is_http_url(v))),
}
"""dict: `SECTOR_FIELDS` plus `update_sector`-only fields -- see
`_validate_sector_fields`'s `allowed` parameter."""


def _validate_sector_fields(body, required, allowed=SECTOR_FIELDS):
    """
    Checks `body` against `allowed`: every key in `required` must be
    present, and every key actually present (required or not) must match
    its expected type and pass its validator.

    Args:
        body (dict): The parsed request body.
        required (set[str]): Field names that must be present.
        allowed (dict): The field shape to validate against -- `SECTOR_FIELDS`
            for `create_sector`, `SECTOR_UPDATE_FIELDS` for `update_sector`
            (see that dict's own docstring for why they differ).

    Raises:
        ApiError: On a missing required field, an unrecognized field, a
            wrong-typed field, or one that fails its validator.
    """
    unknown = set(body) - set(allowed)
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")

    missing = required - set(body)
    if missing:
        raise ApiError(f"missing required field(s): {', '.join(sorted(missing))}")

    for field, (expected_type, is_valid) in allowed.items():
        if field not in body:
            continue
        value = body[field]
        if not isinstance(value, expected_type) or isinstance(value, bool) or not is_valid(value):
            raise ApiError(f"'{field}' is invalid: {value!r}")


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
    body = require_json_body()
    _validate_sector_fields(body, required={"name", "edge_ly"})

    conn = _write_conn()
    try:
        with conn:
            cur = conn.execute(
                "INSERT INTO sectors (name, edge_mpc) VALUES (?, ?)",
                (body["name"], ly_to_milliparsecs(body["edge_ly"])),
            )
            sector_id = cur.lastrowid
    finally:
        conn.close()

    audit("sector.create", target=f"sector:{sector_id}", detail=f"name={body['name']!r} edge_ly={body['edge_ly']}")
    return jsonify({"id": sector_id, "name": body["name"], "edge_ly": body["edge_ly"]}), 201


@bp.route("/sectors/<int:sector_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def update_sector(sector_id):
    """`PATCH /api/sectors/<id>` -- any non-empty subset of `{"name": str,
    "edge_ly": number > 0, "wiki_url": str or null}`. `wiki_url` is the
    admin "manually set the wiki link" affordance (`html/admin.py`) --
    `null` clears it back to "no page yet" (see `schema.sql`'s "v22"
    header note); the same column is also written automatically by
    `POST /api/sectors/<id>/wiki` on a successful upload."""
    body = require_json_body()
    if not body:
        raise ApiError("body must include at least one field to update")
    _validate_sector_fields(body, required=set(), allowed=SECTOR_UPDATE_FIELDS)

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


@bp.route("/sectors/<int:sector_id>", methods=["DELETE"])
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


@bp.route("/sectors/<int:sector_id>/generate-neighborhood", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def generate_sector_neighborhood_route(sector_id):
    """`POST /api/sectors/<id>/generate-neighborhood` -- generates every
    not-yet-generated sector within `radius_ly` (optional JSON body
    field; defaults to `tuning.DEFAULT_GENERATE_RADIUS_PC`,
    12 pc, the sphere `generate.py galaxy`'s own random-start mode uses)
    of this already galaxy-placed sector -- see
    `generate.generate_sector_neighborhood`. Every new sector also gets
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
    body = request.get_json(silent=True)
    if body is None:
        body = {}
    if not isinstance(body, dict):
        raise ApiError("request body must be a JSON object")
    radius_ly = body.get("radius_ly")
    if radius_ly is not None and (
        not isinstance(radius_ly, (int, float)) or isinstance(radius_ly, bool)
        or not 0 < radius_ly <= MAX_NEIGHBORHOOD_RADIUS_LY
    ):
        raise ApiError(f"'radius_ly' is invalid: {radius_ly!r}")
    estimate_only = body.get("estimate_only", False)
    if not isinstance(estimate_only, bool):
        raise ApiError(f"'estimate_only' is invalid: {estimate_only!r}")

    try:
        result = generate.generate_sector_neighborhood(
            sector_id, radius_ly=radius_ly, config=_resolve_requested_write_db_config(),
            estimate_only=estimate_only,
        )
    except ValueError as exc:
        raise ApiError(str(exc), status_code=404)
    except generate.GenerationRefused as exc:
        # PERF.3: the database disk can't hold it; nothing was written.
        raise ApiError(str(exc), status_code=507)
    except RuntimeError as exc:
        # The galaxy's density skeleton (`generate.py plan`) has never
        # been built -- generate_sector_neighborhood needs it to gate each
        # candidate slot's own generation on local stellar density.
        raise ApiError(str(exc), status_code=409)

    if estimate_only:
        return jsonify(result)
    audit(
        "sector.generate_neighborhood", target=f"sector:{sector_id}",
        detail=f"radius_ly={radius_ly!r} generated={result['generated']}",
    )
    return jsonify(result)


SYSTEM_CONFIG_TRISTATE_FIELDS = {
    "habitable_world", "asteroid_belt", "large_star", "moons",
    "max_planets", "planets", "intelligent_life", "binary_system",
    "wide_binary",
}
"""set[str]: `SystemConfig` fields that are `bool` or `None` (`None` = let
the generator decide) -- see `config.py`'s own field docstrings."""

SYSTEM_CONFIG_ALLOWED_FIELDS = SYSTEM_CONFIG_TRISTATE_FIELDS | {"markdown", "star_type", "name", "age", "num_orbits"}
"""set[str]: Every recipe field `POST /api/systems` currently accepts --
`SystemConfig.SERIALIZABLE_FIELDS` minus `SLOTS` (a nested per-orbit
structure not validated in this pass; a request naming it is rejected as
an unrecognized field rather than silently accepted-but-ignored)."""


def _validate_system_config_body(body):
    """
    Validates a `POST /api/systems` recipe body against
    `SYSTEM_CONFIG_ALLOWED_FIELDS` -- deliberately stricter than
    `SystemConfig.from_dict`/`fields_from_dict` itself (which silently
    ignores an unrecognized key rather than rejecting it), since a typo
    in a request body should be a `400`, not a silently-ignored no-op.

    Raises:
        ApiError: On an unrecognized field or one with an invalid type/
            value.
    """
    unknown = set(body) - SYSTEM_CONFIG_ALLOWED_FIELDS
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")

    if "markdown" in body and not isinstance(body["markdown"], bool):
        raise ApiError("'markdown' must be a boolean")
    for field in SYSTEM_CONFIG_TRISTATE_FIELDS:
        if field in body and body[field] is not None and not isinstance(body[field], bool):
            raise ApiError(f"'{field}' must be a boolean or null")
    for field in ("star_type", "name"):
        if field in body and body[field] is not None and not isinstance(body[field], str):
            raise ApiError(f"'{field}' must be a string or null")
    if isinstance(body.get("name"), str) and len(body["name"]) > store.SYSTEM_NAME_MAX_LENGTH:
        raise ApiError(f"'name' must be at most {store.SYSTEM_NAME_MAX_LENGTH} characters")
    if "age" in body and body["age"] not in (None, "young", "old"):
        raise ApiError("'age' must be 'young', 'old', or null")
    if "num_orbits" in body:
        value = body["num_orbits"]
        if value is not None and (
            not isinstance(value, int) or isinstance(value, bool)
            or not 0 < value <= generationLimits.MAX_NUM_ORBITS
        ):
            raise ApiError(
                f"'num_orbits' must be a positive integer up to {generationLimits.MAX_NUM_ORBITS}, or null"
            )


def _placement_fields(body):
    """Pops and checks `POST /api/systems`' placement fields: `sector_id`
    (a positive integer) and `position` (three finite light-year numbers,
    sector-local, only with `sector_id`). Returns `(sector_id, position)`."""
    sector_id = body.pop("sector_id", None)
    position = body.pop("position", None)
    if sector_id is not None and (not isinstance(sector_id, int) or isinstance(sector_id, bool) or sector_id < 1):
        raise ApiError("'sector_id' must be a positive integer or null")
    if position is not None:
        if sector_id is None:
            raise ApiError("'position' needs 'sector_id'")
        if (not isinstance(position, list) or len(position) != 3
                or not all(isinstance(c, (int, float)) and not isinstance(c, bool) and math.isfinite(c)
                           for c in position)):
            raise ApiError("'position' must be three finite numbers (light-years from the sector's center)")
        position = tuple(float(c) for c in position)
    return sector_id, position


def _sector_generation_context(conn, sector_id, system_config):
    """For a galaxy-placed sector: the system's distance from the galactic
    center in light-years, after giving `system_config` the sector's
    stellar population and (once a bright-star scatter ran) its dim-star
    cap (its own backfill level, GEN.44), as a sector fill would
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
        level = store.bright_star_fill_level(conn, placement["ring_index"], placement["layer_index"],
                                           placement["ring_slot_index"])
        center = (placement["center_x_pc"], placement["center_y_pc"], placement["center_z_pc"])
        brightStars.FillContext(center, skeleton.shape, min_luminosity_sol=level).apply(system_config)
    return pc_to_ly(placement["galactic_radius_pc"])


def _generate_system(system_config, galactic_center_dist_ly=None):
    """Generates a `StarSystem` from a validated recipe; a generation
    failure is the request's fault, so a 400, not a 500."""
    try:
        return StarSystem(system_config=system_config, galactic_center_dist_ly=galactic_center_dist_ly)
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
    body = require_json_body()
    sector_id, position = _placement_fields(body)
    _validate_system_config_body(body)
    system_config = SystemConfig.from_dict(body)

    conn = _write_conn()
    try:
        with conn:
            if sector_id is None:
                system_id = store.insert_star_system(conn, _generate_system(system_config), system_config)
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

    audit("system.create", target=f"system:{system_id}",
          detail=f"star_type={body.get('star_type')!r} sector_id={sector_id!r}")
    result = {"id": system_id}
    if placed is not None:
        result.update(sector_id=sector_id, position=list(placed))
    return jsonify(result), 201


NAME_MAX_LENGTH = MAX_NAME_LENGTH
"""int: Every `name` column is `VARCHAR(255)` (see `schema.sql`)."""

_NAME_CLASH_LABELS = {
    "sectors": "a sector", "star_systems": "a star system", "stars": "a star",
}


def _rename_body(max_length=NAME_MAX_LENGTH):
    """
    Validates a rename request's body -- exactly `{"name": str}`, not
    blank once trimmed, at most `max_length` characters -- and
    returns the trimmed name.

    Raises:
        ApiError: 400 on any other shape.
    """
    body = require_json_body()
    unknown = set(body) - {"name"}
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")
    return _check_name(body.get("name"), max_length)


def _check_name(name, max_length=NAME_MAX_LENGTH):
    """A new name, trimmed and checked as `_rename_body` describes. A
    system or star passes `store.SYSTEM_NAME_MAX_LENGTH`: its planets and
    moons are named after it, so it needs room to grow."""
    if not isinstance(name, str) or not name.strip():
        raise ApiError("'name' must be a non-empty string")
    name = " ".join(name.split())
    if len(name) > max_length:
        raise ApiError(f"'name' must be at most {max_length} characters")
    return name


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


_SYSTEM_PATCH_FIELDS = {"name", "regenerate", "drop_facilities"}


@bp.route("/systems/<int:system_id>", methods=["PATCH"])
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
    body = require_json_body()
    unknown = set(body) - _SYSTEM_PATCH_FIELDS
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")
    if "name" not in body and "regenerate" not in body:
        raise ApiError("send 'name', 'regenerate', or both")
    name = _check_name(body["name"], store.SYSTEM_NAME_MAX_LENGTH) if "name" in body else None
    recipe = body.get("regenerate")
    drop_facilities = body.get("drop_facilities", False)
    if not isinstance(drop_facilities, bool):
        raise ApiError("'drop_facilities' must be a boolean")
    if "regenerate" in body:
        if not isinstance(recipe, dict):
            raise ApiError("'regenerate' must be an object (a recipe, or {} for defaults)")
        _validate_system_config_body(recipe)
        if "name" in recipe:
            raise ApiError("'regenerate' keeps the system's name; rename with 'name' instead")

    conn = _write_conn()
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
            final_name = conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"]
    finally:
        conn.close()

    detail = {key: body[key] for key in ("name", "drop_facilities") if key in body}
    if recipe is not None:
        detail["regenerate"] = recipe
    audit("system.update", target=f"system:{system_id}", detail=str(detail))
    return jsonify({"status": "ok", "id": system_id, "name": final_name, "regenerated": recipe is not None})


@bp.route("/stars/<int:star_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def rename_star(star_id):
    """`PATCH /api/stars/<id>` `{"name": str}` -- renames a star. A single
    star shares its system's name, so this renames the system too; a
    binary's star is renamed on its own, with the planets and moons named
    after it (`store.rename_star`). 409 if a sector, system or star already has the
    name."""
    name = _rename_body(store.SYSTEM_NAME_MAX_LENGTH)

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
    name = _rename_body()

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


@bp.route("/planets/<int:planet_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def rename_planet(planet_id):
    """`PATCH /api/planets/<id>` `{"name": str}` -- renames one planet. Its
    moons keep their names. 409 if a sector, system or star already has the name."""
    return _rename_planet_or_moon("planets", "planet", planet_id)


@bp.route("/moons/<int:moon_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def rename_moon(moon_id):
    """`PATCH /api/moons/<id>` `{"name": str}` -- renames one moon. 409 if
    a sector, system or star already has the name."""
    return _rename_planet_or_moon("moons", "moon", moon_id)


@bp.route("/systems/<int:system_id>", methods=["DELETE"])
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

FACILITY_FIELDS = {
    "name": (str, lambda v: bool(v.strip()) and len(v) <= MAX_NAME_LENGTH),
    "kind": (str, lambda v: v in tuning.FACILITY_KINDS),
    "placement": (str, lambda v: v in facility_rules.PLACEMENTS),
    "host_type": (str, lambda v: v in facility_rules.HOST_TYPES),
    "host_id": (int, lambda v: v > 0),
    "distance_km": ((int, float), lambda v: math.isfinite(v) and v > 0),
    "phase_deg": ((int, float), math.isfinite),
    "offset_ly": (list, lambda v: len(v) == 3 and all(
        isinstance(c, (int, float)) and not isinstance(c, bool) and math.isfinite(c) for c in v)),
    "description": (str, lambda v: len(v) <= 4000),
}
"""dict: The `POST /api/facilities` body shape."""


@bp.route("/facilities/<int:facility_id>")
def facility(facility_id):
    """`GET /api/facilities/<id>` -- one facility."""
    found = facility_detail(get_db(), facility_id)
    if found is None:
        raise ApiError(f"no such facility: {facility_id}", status_code=404)
    return jsonify(found)


@bp.route("/systems/<int:system_id>/facilities")
def system_facilities(system_id):
    """`GET /api/systems/<id>/facilities` -- every facility in a system.
    404 for an unknown system."""
    db = get_db()
    if db.execute("SELECT 1 FROM star_systems WHERE id = ?", (system_id,)).fetchone() is None:
        raise ApiError(f"no such system: {system_id}", status_code=404)
    return jsonify({"items": facilities_for_system(db, system_id)})


@bp.route("/sectors/<int:sector_id>/facilities")
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


@bp.route("/facilities/orbit")
def facility_orbit():
    """`GET /api/facilities/orbit?host_type=star|planet|moon&host_id=N[&distance_km=X]`
    -- the orbit (distance, period, speed) an orbital facility would get,
    without saving anything, so a form can show it first, plus
    `min_distance_km`/`max_distance_km`, the orbits the host allows (just
    above its surface to the edge of its sphere of influence)."""
    host_type = request.args.get("host_type", "")
    try:
        host_id = int(request.args.get("host_id", ""))
        raw_distance = request.args.get("distance_km")
        distance_km = None if raw_distance in (None, "") else float(raw_distance)
    except ValueError:
        raise ApiError("host_id must be a whole number and distance_km a number")
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
    body = require_json_body()
    _validate_sector_fields(body, required={"name", "kind", "placement", "host_type", "host_id"},
                            allowed=FACILITY_FIELDS)
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


@bp.route("/facilities/<int:facility_id>", methods=["DELETE"])
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

def _wiki_client_for(backend):
    """
    Builds a `wikiClient.WikiClient` for `backend`, from
    `current_app.config["WIKI_CONFIG"]` (see `config.py`'s `_wiki_config`).

    Args:
        backend (str): `"wikijs"` or `"mediawiki"`.

    Returns:
        wikiClient.WikiClient

    Raises:
        ApiError: 501 if `backend` has no `base_url`/credentials
            configured at all -- rather than trying, and failing, to
            reach an empty URL.
    """
    settings = current_app.config["WIKI_CONFIG"][backend]
    if not settings["configured"]:
        raise ApiError(f"wiki publishing is not configured for {backend!r}", status_code=501)

    if backend == "wikijs":
        return WikiClient(backend="wikijs", base_url=settings["base_url"], api_token=settings["api_token"])
    return WikiClient(
        backend="mediawiki", base_url=settings["base_url"],
        username=settings["username"], password=settings["password"],
    )


def _wiki_upload_request(body):
    """
    Validates a `POST .../wiki` request body.

    Args:
        body (dict): The parsed request body -- `{"backend": "wikijs" |
            "mediawiki", "path": str}`. `path` is the target page's
            path/slug for `"wikijs"` (which addresses a page separately
            from its title -- see `wikiClient/wikijs.py`) and required
            for it; `"mediawiki"` has no such separate concept (its
            title *is* its address, see `wikiClient/mediawiki.py`), so
            `path` is accepted but ignored for it.

    Returns:
        tuple[str, str or None]: `(backend, path)` -- `path` stripped, or
            `None` if not given/blank.

    Raises:
        ApiError: On an unrecognized field, a missing/invalid `backend`,
            or a missing/blank `path` when `backend == "wikijs"`.
    """
    unknown = set(body) - {"backend", "path"}
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")

    backend = body.get("backend")
    if backend not in ("wikijs", "mediawiki"):
        raise ApiError("'backend' must be 'wikijs' or 'mediawiki'")

    path = body.get("path")
    if path is not None and not isinstance(path, str):
        raise ApiError("'path' must be a string")
    path = path.strip() if path else None

    if backend == "wikijs" and not path:
        raise ApiError("'path' is required for the 'wikijs' backend")

    return backend, path


def _create_wiki_page(client, path, title, content):
    """Wraps `client.create_page`, mapping `wikiClient`'s own exception
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


@bp.route("/systems/<int:system_id>/wiki", methods=["POST"])
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
    client = _wiki_client_for(backend)

    try:
        system = query_system_detail(get_db(), system_id)
        content = render_system_text(get_db(), system_id, "markdown" if backend == "wikijs" else "wikitext")
    except ValueError:
        raise ApiError(f"no such system: {system_id}", status_code=404)

    page_path = path if backend == "wikijs" else system["name"]
    page = _create_wiki_page(client, page_path, system["name"], content)

    url_column = "wikijs_url" if backend == "wikijs" else "mediawiki_url"
    conn = _write_conn()
    try:
        with conn:
            conn.execute(f"UPDATE star_systems SET {url_column} = ? WHERE id = ?", (page.url, system_id))
    finally:
        conn.close()

    audit("system.wiki_upload", target=f"system:{system_id}", detail=f"backend={backend!r} path={page.path!r}")
    return jsonify({"id": page.id, "path": page.path, "title": page.title, "url": page.url}), 201


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
        f'**Cube edge:** {format_distance_ly(sector["edge_ly"])}  \n'
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
        f"'''Cube edge:''' {format_distance_ly(sector['edge_ly'])}\n\n"
        f"'''Systems:''' {sector['system_count']}\n\n"
        "== Systems ==\n\n"
        '{| class="wikitable"\n'
        "! Name !! Octant !! Binary !! Star type !! Location\n"
        "|-\n"
        f"{wiki_rows}\n"
        "|}\n"
    )
    return markdown_content, wikitext_content


@bp.route("/sectors/<int:sector_id>/wiki", methods=["POST"])
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
    client = _wiki_client_for(backend)

    try:
        sector = query_sector_detail(get_db(), sector_id)
    except ValueError:
        raise ApiError(f"no such sector: {sector_id}", status_code=404)

    markdown_content, wikitext_content = _sector_wiki_content(sector)
    content = markdown_content if backend == "wikijs" else wikitext_content
    page_path = path if backend == "wikijs" else sector["name"]
    page = _create_wiki_page(client, page_path, sector["name"], content)

    conn = _write_conn()
    try:
        with conn:
            conn.execute("UPDATE sectors SET wiki_url = ? WHERE id = ?", (page.url, sector_id))
    finally:
        conn.close()

    audit("sector.wiki_upload", target=f"sector:{sector_id}", detail=f"backend={backend!r} path={page.path!r}")
    return jsonify({"id": page.id, "path": page.path, "title": page.title, "url": page.url}), 201
