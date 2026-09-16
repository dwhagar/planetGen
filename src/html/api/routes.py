# html/api/routes.py

"""
JSON endpoints over the planetGen database.

The read endpoints each open a strictly-read-only connection
(`queryDb.open_readonly` -- same `file:...?mode=ro`-equivalent guarantee
the `queryDb.py` CLI already relies on, see that module's docstring for
how "read-only" is enforced by the configured account's grants). Listing
endpoints delegate straight to `queryDb.py`'s existing `list_sectors`/
`list_systems`/`systems_within_radius` functions rather than
re-implementing the same SQL a third time; the two detail endpoints go
through `stellarObjects._db.load_sector`/`load_star_system` for the full
nested object graph, serialized via each class's own `to_dict()` (Phase 1
serialization, `stellarObjects/serialization.py`).

Listing endpoints (`/sectors`, `/systems`) are paginated -- this project's
own roadmap (docs/TODO.md, Phase 4) plans galaxy-scale generation, so an
unbounded `SELECT *` here would eventually return an unbounded response.
`_paginate` centralizes parsing/validating `limit`/`offset` so both routes
apply the same defaults, cap, and error behavior.

The write endpoints (`POST`/`PATCH`/`DELETE` on `/sectors`/`/systems`,
near the bottom of this file) do real inserts/updates/deletes now, gated
behind admin authentication (`authz.require_admin`) and run against a
separate, write-capable database account -- see `docs/api.md`'s
"Write endpoints" section and `stellarObjects/adminAuth.py`.
"""

from flask import Blueprint, current_app, g, jsonify, request

from queryDb import (
    NO_SECTOR,
    NavUnavailable,
    SEARCH_TAG_FACETS,
    count_sectors,
    count_systems,
    galaxy_placed_phenomena,
    galaxy_placed_sectors,
    list_sectors,
    list_systems,
    nav_between,
    open_readonly,
    search as run_search,
    sector_detail as query_sector_detail,
    system_detail as query_system_detail,
    systems_within_radius,
)
from stellarObjects import _db
from stellarObjects._db import MySQLConfig, list_databases, resolve_database
from stellarObjects.config import SystemConfig
from stellarObjects.systemData import StarSystem
from stellarObjects.utils import ly_to_milliparsecs

from .authz import audit, require_admin
from .common import ApiError, require_json_body
from .limiter import limiter

bp = Blueprint("api", __name__, url_prefix="/api")

DEFAULT_PAGE_LIMIT = 100
MAX_PAGE_LIMIT = 500

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
    `stellarObjects._db.list_databases`/`resolve_database`) -- the same
    validation `html/lib/dbutil.resolve_db_name` has always applied for
    the CGI browser's own `?db=` picker, now shared by this API so a
    request can't select a schema this deployment never meant to expose.

    Returns:
        MySQLConfig: Ready to pass to `open_readonly`.

    Raises:
        ApiError: 404, if `db` is given but doesn't match a listed schema.
    """
    base_config = current_app.config["MYSQL_CONFIG"]
    requested = request.args.get("db")
    if requested is None:
        return base_config
    try:
        return resolve_database(base_config, requested)
    except ValueError as exc:
        raise ApiError(str(exc), status_code=404)


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
    """
    if "db" not in g:
        try:
            g.db = open_readonly(_resolve_requested_db_config())
        except SystemExit as exc:
            raise ApiError(str(exc), status_code=503)
    return g.db


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
@limiter.exempt
def health():
    """
    Liveness/readiness check for monitoring -- confirms the process is up
    and the configured database can actually be opened, not just that
    Flask is responding. Returns 503 (rather than letting the connection
    error propagate to a generic 500) so a monitor can tell "the API is
    running but its database is unreachable" apart from "the API itself is
    down".
    """
    try:
        get_db().execute("SELECT 1")
    except Exception as exc:
        return jsonify({"status": "error", "detail": str(exc)}), 503
    return jsonify({"status": "ok"})


@bp.route("/databases")
def databases():
    """
    Lists every MySQL schema on the configured server whose name matches
    this deployment's prefix (`stellarObjects._db.list_databases`), each
    with its size/last-modified stats plus a quick-glance sector/system
    count -- the data `html/index.py`'s database picker needs, and (via
    the sidenav's "Databases" link, shown on every page) whether that
    picker has anything to offer at all. Every other endpoint's own
    `?db=` selects among these same names (see `get_db`).
    """
    base_config = current_app.config["MYSQL_CONFIG"]
    entries = list_databases(base_config)
    items = []
    for entry in entries:
        conn = open_readonly(MySQLConfig(
            host=base_config.host, port=base_config.port,
            user=base_config.user, password=base_config.password, database=entry["name"],
        ))
        try:
            sector_count = conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"]
            system_count = conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]
        except Exception:
            # A schema matching the configured prefix but missing this
            # project's own tables (e.g. mid-migration, or a stray
            # unrelated database sharing the prefix) shouldn't take down
            # the whole listing -- report it with unknown counts instead.
            sector_count = system_count = None
        finally:
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
    if radius <= 0:
        raise ApiError("radius must be greater than 0")

    try:
        matches = systems_within_radius(get_db(), system_id, radius)
    except SystemExit as exc:
        # systems_within_radius is shared with the queryDb.py CLI and raises
        # SystemExit (its CLI-appropriate error signal) for a missing/
        # unplaced system id -- caught here rather than changing its shared
        # behavior just for this one caller.
        return jsonify({"error": str(exc)}), 404
    return jsonify(matches)


def _parse_system_id_param(query_args, name):
    """
    Parses a required `<name>` query parameter as a `star_systems.id`.

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


@bp.route("/nav")
def nav():
    """
    Course, distance, and optimal route between two systems -- see
    `queryDb.nav_between` for the full availability rules (a system not
    assigned to any sector, or two systems in different sectors where
    either sector lacks a galaxy placement, both mean NAV isn't
    available for that pair) and docs/api.md for the response shape.
    """
    from_id = _parse_system_id_param(request.args, "from")
    to_id = _parse_system_id_param(request.args, "to")

    try:
        result = nav_between(get_db(), from_id, to_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    except NavUnavailable as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify({
        "scope": result["scope"],
        "direct": result["direct"]._asdict(),
        "warp_times": [leg._asdict() for leg in result["warp_times"]],
        "origin_position": result["origin_position"],
        "destination_position": result["destination_position"],
        "route": result["route"],
    })


@bp.route("/galaxy/sectors")
def galaxy_sectors():
    """
    Every galaxy-placed sector (`sectors.center_x/y/z_pc` not NULL), with
    its live system count -- the data `html/galaxy.py`'s Galaxy Map plots.
    Not paginated: bounded by how much of the galaxy has actually been
    generated so far (see `docs/TODO.md`'s Phase 4 lazy-generation design),
    not by the addressable galaxy's own astronomical scale.
    """
    return jsonify({"items": galaxy_placed_sectors(get_db())})


@bp.route("/galaxy/phenomena")
def galaxy_phenomena():
    """
    Every galaxy-placed nebula/asteroid field (`nebulae`/`asteroid_fields`.
    `center_x/y/z_pc` not NULL -- see `schema.sql`'s "v18" header note) --
    the phenomenon counterpart to `/api/galaxy/sectors`, plotted as small
    dots on the same `html/galaxy.py` Galaxy Map. Not paginated, for the
    same reason `/api/galaxy/sectors` isn't.
    """
    return jsonify({"items": galaxy_placed_phenomena(get_db())})


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

    return jsonify(run_search(get_db(), texts, tags, sizes=sizes))


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
    "name": (str, lambda v: bool(v.strip())),
    "edge_ly": ((int, float), lambda v: v > 0),
}
"""dict: The `sectors` JSON object shape both `create_sector` and
`update_sector` validate against -- see docs/api.md's "Sectors" write
schema. Shared so the two can never validate the same field two
different ways."""


def _validate_sector_fields(body, required):
    """
    Checks `body` against `SECTOR_FIELDS`: every key in `required` must be
    present, and every key actually present (required or not) must match
    its expected type and pass its validator.

    Args:
        body (dict): The parsed request body.
        required (set[str]): Field names that must be present.

    Raises:
        ApiError: On a missing required field, an unrecognized field, a
            wrong-typed field, or one that fails its validator.
    """
    unknown = set(body) - set(SECTOR_FIELDS)
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")

    missing = required - set(body)
    if missing:
        raise ApiError(f"missing required field(s): {', '.join(sorted(missing))}")

    for field, (expected_type, is_valid) in SECTOR_FIELDS.items():
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
    base_config = current_app.config["WRITE_MYSQL_CONFIG"]
    requested = request.args.get("db")
    if requested is None:
        return base_config
    try:
        return resolve_database(base_config, requested)
    except ValueError as exc:
        raise ApiError(str(exc), status_code=404)


def _write_conn():
    """Opens a write-capable connection against the content database
    `?db=` (or the configured default) selects -- `ensure_schema=False`,
    same reasoning as `queryDb.open_readonly`/`_db.open_write`: the
    write-capable account has no `CREATE` grant (see `config.py`), so the
    schema must already exist."""
    return _db.open_write(_resolve_requested_write_db_config())


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
    "edge_ly": number > 0}`."""
    body = require_json_body()
    if not body:
        raise ApiError("body must include at least one field to update")
    _validate_sector_fields(body, required=set())

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
            deleted = conn.execute("DELETE FROM sectors WHERE id = ?", (sector_id,)).rowcount > 0
    finally:
        conn.close()

    if not deleted:
        raise ApiError(f"no such sector: {sector_id}", status_code=404)
    audit("sector.delete", target=f"sector:{sector_id}")
    return jsonify({"status": "ok"})


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
    if "age" in body and body["age"] not in (None, "young", "old"):
        raise ApiError("'age' must be 'young', 'old', or null")
    if "num_orbits" in body:
        value = body["num_orbits"]
        if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value <= 0):
            raise ApiError("'num_orbits' must be a positive integer or null")


@bp.route("/systems", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def create_system():
    """
    `POST /api/systems` -- a generation "recipe" body (see
    `SYSTEM_CONFIG_ALLOWED_FIELDS`), the same shape `systemGen.py
    --system-file` already takes (`SystemConfig.from_dict`). Generates a
    new **standalone** system (no `sector_id` -- attaching a newly
    generated system to an existing sector needs that sector's own
    placement/Hill-sphere separation logic, tracked as follow-up work in
    docs/TODO.md, not implemented here) and persists it via the same
    `insert_star_system` every other entry point uses. Returns the new
    `star_systems.id`.
    """
    body = require_json_body()
    _validate_system_config_body(body)
    system_config = SystemConfig.from_dict(body)

    try:
        star_system = StarSystem(system_config=system_config)
    except Exception as exc:
        # Generation can reject an internally-inconsistent recipe (e.g. an
        # impossible num_orbits/planet-class combination) -- surfaced as a
        # 400 (the request body was the problem), not an unhandled 500.
        raise ApiError(f"system generation failed: {exc}")

    conn = _write_conn()
    try:
        with conn:
            system_id = _db.insert_star_system(conn, star_system, system_config)
    finally:
        conn.close()

    audit("system.create", target=f"system:{system_id}", detail=f"star_type={body.get('star_type')!r}")
    return jsonify({"id": system_id}), 201


SYSTEM_UPDATE_FIELDS = {
    "name": (str, lambda v: bool(v.strip())),
}
"""dict: `PATCH /api/systems/<id>`'s allowed fields -- metadata only
(a rename) in this pass; editing a system's generated content (stars/
planets/moons/belts) is out of scope here, same reasoning as `POST
/api/systems` accepting only a generation recipe rather than a
hand-edited object graph (see docs/api.md)."""


@bp.route("/systems/<int:system_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def update_system(system_id):
    """`PATCH /api/systems/<id>` `{"name": str}` -- renames a system."""
    body = require_json_body()
    if not body:
        raise ApiError("body must include at least one field to update")
    unknown = set(body) - set(SYSTEM_UPDATE_FIELDS)
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")
    for field, (expected_type, is_valid) in SYSTEM_UPDATE_FIELDS.items():
        if field not in body:
            continue
        value = body[field]
        if not isinstance(value, expected_type) or not is_valid(value):
            raise ApiError(f"'{field}' is invalid: {value!r}")

    conn = _write_conn()
    try:
        with conn:
            if conn.execute("SELECT id FROM star_systems WHERE id = ?", (system_id,)).fetchone() is None:
                raise ApiError(f"no such system: {system_id}", status_code=404)
            conn.execute("UPDATE star_systems SET name = ? WHERE id = ?", (body["name"], system_id))
    finally:
        conn.close()

    audit("system.update", target=f"system:{system_id}", detail=str(body))
    return jsonify({"status": "ok"})


@bp.route("/systems/<int:system_id>", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def delete_system(system_id):
    """`DELETE /api/systems/<id>` -- no request body. Its stars/planets/
    moons/belts are deleted with it (all `ON DELETE CASCADE`, per
    `schema.sql`)."""
    conn = _write_conn()
    try:
        with conn:
            deleted = conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,)).rowcount > 0
    finally:
        conn.close()

    if not deleted:
        raise ApiError(f"no such system: {system_id}", status_code=404)
    audit("system.delete", target=f"system:{system_id}")
    return jsonify({"status": "ok"})
