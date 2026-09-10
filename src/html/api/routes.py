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
near the bottom of this file) are deliberately stubs: they validate their
request body against the JSON shape documented in `docs/api.md` and
return `501 Not Implemented`, but never touch the database. Wiring them
up to real inserts/updates/deletes is future work -- see each stub's own
docstring and `config.py`'s module docstring for what that will also
require (a write-capable database account, and some form of
authentication/authorization, neither of which exists yet).
"""

from flask import Blueprint, current_app, g, jsonify, request

from queryDb import (
    NO_SECTOR,
    NavUnavailable,
    SEARCH_TAG_FACETS,
    count_sectors,
    count_systems,
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
from stellarObjects._db import MySQLConfig, list_databases, resolve_database

from .limiter import limiter

bp = Blueprint("api", __name__, url_prefix="/api")

DEFAULT_PAGE_LIMIT = 100
MAX_PAGE_LIMIT = 500

WRITE_RATE_LIMIT = "10 per minute"
"""str: Applied to every write stub, on top of the app-wide default
limits (`config.Config.RATELIMIT_DEFAULT`) -- both apply together
(Flask-Limiter doesn't replace the default with a route's own `@limit`
unless told to), so a mutating endpoint is always throttled at least
this tightly regardless of how the global default is configured. Tighter
than any read endpoint's effective limit is today, on the general
principle that a request able to change data deserves more headroom
against abuse than one that only reads it -- revisit once these are real
writes rather than stubs, since the right number depends on actual usage
patterns this project doesn't have yet."""


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
    """
    if "db" not in g:
        g.db = open_readonly(_resolve_requested_db_config())
    return g.db


def close_db(exception=None):
    db = g.pop("db", None)
    if db is not None:
        db.close()


class ApiError(Exception):
    """
    Raised by a route to end the request with a JSON `{"error": ...}` body
    and a specific status code, handled by `app.py`'s error handler. Beats
    each route hand-rolling its own `return jsonify(...), status` for
    validation failures, so every 400 in this API is worded and shaped the
    same way.
    """

    def __init__(self, message, status_code=400):
        super().__init__(message)
        self.message = message
        self.status_code = status_code


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


def _require_json_body():
    """
    Parses the request body as JSON, for every write endpoint below.

    `request.get_json(silent=True)` returns `None` for a missing/empty
    body, a body that isn't valid JSON, *or* a `Content-Type` other than
    `application/json` -- this API doesn't need to tell those apart for a
    caller, they're all just "you didn't send a JSON object".

    Returns:
        dict: The parsed body.

    Raises:
        ApiError: If the body is missing, isn't valid JSON, or isn't a
            JSON *object* (a bare list/string/number is valid JSON but
            not a usable request body here).
    """
    body = request.get_json(silent=True)
    if not isinstance(body, dict):
        raise ApiError("request body must be a JSON object")
    return body


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


@bp.route("/search")
def search():
    """
    Faceted search over the chosen database: click-to-filter attribute
    tags (object type; star spectral/luminosity class; planet/moon class,
    body type, supported life chemistry; asteroid belt density) plus a
    per-entity name search -- see `queryDb.search`'s own docstring for the
    full request/response shape. `html/search.py` is a thin renderer over
    this endpoint's JSON.

    Query parameters: `sector_q`/`system_q`/`star_q`/`planet_q`/`moon_q`
    (name search terms) plus every name in `queryDb.SEARCH_TAG_FACETS`
    (repeatable, e.g. `?class=M&class=K` for two active Class tags).
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

    return jsonify(run_search(get_db(), texts, tags))


# ---------------------------------------------------------------------
# Write endpoints -- stubs. See the module docstring: these validate
# their request body against the shape documented in docs/api.md and
# always respond 501, without touching the database. `_require_json_body`
# doesn't get to know that, though, so it still runs against every
# request, malformed or not, exactly as it will once these are real.
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


@bp.route("/sectors", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
def create_sector():
    """
    Stub. Request body: `{"name": str, "edge_ly": number > 0}` (both
    required) -- see docs/api.md. Validates the body, then always
    responds 501; no `sectors` row is ever inserted.
    """
    body = _require_json_body()
    _validate_sector_fields(body, required={"name", "edge_ly"})
    raise ApiError("sector creation is not implemented yet", status_code=501)


@bp.route("/sectors/<int:sector_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
def update_sector(sector_id):
    """
    Stub. Request body: any non-empty subset of `{"name": str,
    "edge_ly": number > 0}` -- see docs/api.md. Validates the body, then
    always responds 501; no `sectors` row is ever modified. Doesn't
    check `sector_id` actually exists first (there's nothing to apply
    the update to yet either way), so a nonexistent id also gets 501,
    not 404 -- revisit once this does something real.
    """
    body = _require_json_body()
    if not body:
        raise ApiError("body must include at least one field to update")
    _validate_sector_fields(body, required=set())
    raise ApiError(f"sector {sector_id} modification is not implemented yet", status_code=501)


@bp.route("/sectors/<int:sector_id>", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
def delete_sector(sector_id):
    """Stub. No request body. Always responds 501; no `sectors` row is
    ever deleted."""
    raise ApiError(f"sector {sector_id} deletion is not implemented yet", status_code=501)


@bp.route("/systems", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
def create_system():
    """
    Stub. Request body must be a JSON object -- unlike sectors, this
    project hasn't yet decided a system's create schema (a generation
    "recipe" shaped like `SystemConfig`, e.g. `star_type`/`planets`/
    `moons`, closer to what `systemGen.py --star-type ...` takes; vs. a
    fully-specified object graph shaped like `StarSystem.to_dict()`, with
    every star/planet/moon/belt spelled out; or supporting both) -- see
    docs/api.md. Always responds 501; no `star_systems` row is ever
    inserted.
    """
    _require_json_body()
    raise ApiError("system creation is not implemented yet", status_code=501)


@bp.route("/systems/<int:system_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
def update_system(system_id):
    """Stub. Request body must be a JSON object -- field-level schema
    still TBD, same open question as `create_system`. Always responds
    501; no `star_systems` row is ever modified."""
    _require_json_body()
    raise ApiError(f"system {system_id} modification is not implemented yet", status_code=501)


@bp.route("/systems/<int:system_id>", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
def delete_system(system_id):
    """Stub. No request body. Always responds 501; no `star_systems` row
    (or its stars/planets/moons/belts) is ever deleted."""
    raise ApiError(f"system {system_id} deletion is not implemented yet", status_code=501)
