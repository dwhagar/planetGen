# api/routes.py

"""
Read-only JSON endpoints over the planetGen database.

Every route opens a strictly-read-only connection (`queryDb.open_readonly`
-- same `file:...?mode=ro` guarantee the `queryDb.py` CLI already relies on,
so this API can never accidentally write to a database another process is
using). Listing endpoints delegate straight to `queryDb.py`'s existing
`list_sectors`/`list_systems`/`systems_within_radius` functions rather than
re-implementing the same SQL a third time; the two detail endpoints go
through `stellarObjects._db.load_sector`/`load_star_system` for the full
nested object graph, serialized via each class's own `to_dict()` (Phase 1
serialization, `stellarObjects/serialization.py`).

Listing endpoints (`/sectors`, `/systems`) are paginated -- this project's
own roadmap (docs/TODO.md, Phase 4) plans galaxy-scale generation, so an
unbounded `SELECT *` here would eventually return an unbounded response.
`_paginate` centralizes parsing/validating `limit`/`offset` so both routes
apply the same defaults, cap, and error behavior.
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

bp = Blueprint("api", __name__, url_prefix="/api")

DEFAULT_PAGE_LIMIT = 100
MAX_PAGE_LIMIT = 500


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


@bp.route("/health")
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
