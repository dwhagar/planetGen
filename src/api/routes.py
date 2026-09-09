# api/routes.py

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
    count_sectors,
    count_systems,
    list_sectors,
    list_systems,
    open_readonly,
    systems_within_radius,
)
from stellarObjects._db import load_sector, load_star_system

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


def get_db():
    """
    Returns the request-scoped read-only connection, opening one on first
    use. Reused for the lifetime of the request instead of one connection
    per query, then closed by `close_db` in the app's teardown handler.
    """
    if "db" not in g:
        g.db = open_readonly(current_app.config["MYSQL_CONFIG"])
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
        sector = load_sector(get_db(), sector_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(sector.to_dict())


@bp.route("/systems")
def systems():
    star_type = request.args.get("star_type")

    raw_sector_id = request.args.get("sector_id")
    sector_id = None
    if raw_sector_id is not None:
        try:
            sector_id = int(raw_sector_id)
        except ValueError:
            raise ApiError(f"sector_id must be an integer, got {raw_sector_id!r}")

    limit, offset = _paginate(request.args)
    db = get_db()
    rows = list_systems(db, star_type_prefix=star_type, sector_id=sector_id, limit=limit, offset=offset)
    return jsonify({
        "items": [dict(row) for row in rows],
        "total": count_systems(db, star_type_prefix=star_type, sector_id=sector_id),
        "limit": limit,
        "offset": offset,
    })


@bp.route("/systems/<int:system_id>")
def system_detail(system_id):
    try:
        system = load_star_system(get_db(), system_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(system.to_dict())


@bp.route("/systems/<int:system_id>/near")
def systems_near(system_id):
    # TODO: this endpoint already computes the distance between two placed
    # systems (via `systems_within_radius`) but returns bare JSON -- there's
    # no way to actually see the two points, just a number. A rendered
    # image (even a simple 2D projection) or a small web-page
    # visualization of "here's system A, here's system B, here's the line
    # between them in galactic space" would build on this route's existing
    # data rather than needing new queries. See docs/TODO.md, "Investigate
    # Further".
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
