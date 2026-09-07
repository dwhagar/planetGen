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
"""

from flask import Blueprint, current_app, g, jsonify, request

from queryDb import list_sectors, list_systems, open_readonly, systems_within_radius
from stellarObjects._db import load_sector, load_star_system

bp = Blueprint("api", __name__, url_prefix="/api")


def get_db():
    """
    Returns the request-scoped read-only connection, opening one on first
    use. Reused for the lifetime of the request instead of one connection
    per query, then closed by `close_db` in the app's teardown handler.
    """
    if "db" not in g:
        g.db = open_readonly(current_app.config["DB_PATH"])
    return g.db


def close_db(exception=None):
    db = g.pop("db", None)
    if db is not None:
        db.close()


@bp.route("/sectors")
def sectors():
    return jsonify(list_sectors(get_db()))


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
    sector_id = request.args.get("sector_id", type=int)
    rows = list_systems(get_db(), star_type_prefix=star_type, sector_id=sector_id)
    return jsonify([dict(row) for row in rows])


@bp.route("/systems/<int:system_id>")
def system_detail(system_id):
    try:
        system = load_star_system(get_db(), system_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(system.to_dict())


@bp.route("/systems/<int:system_id>/near")
def systems_near(system_id):
    radius = request.args.get("radius", type=float)
    if radius is None:
        return jsonify({"error": "radius query parameter is required"}), 400
    try:
        matches = systems_within_radius(get_db(), system_id, radius)
    except SystemExit as exc:
        # systems_within_radius is shared with the queryDb.py CLI and raises
        # SystemExit (its CLI-appropriate error signal) for a missing/
        # unplaced system id -- caught here rather than changing its shared
        # behavior just for this one caller.
        return jsonify({"error": str(exc)}), 404
    return jsonify(matches)
