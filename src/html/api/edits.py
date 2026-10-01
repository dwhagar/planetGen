# html/api/edits.py

"""
Admin editing endpoints (TODO ADM.1): delete and regenerate one planet,
moon, asteroid belt, phenomenon or sector (ADM.8). Every route needs an
admin whose credentials are current, writes an audit-log row, and answers
with what the edit did.

A body edit loads its system (`_db.load_star_system`), changes it with
`stellarObjects.adminEdits`, re-validates it from the moons outward
(`stellarObjects.validation`) and writes it back in place
(`stellarObjects.editStore`), so the other bodies keep their rows. The
answer lists the bodies the re-validation moved and any warnings left.

Kept in its own blueprint, like `population.py`, so this work doesn't
touch `routes.py`.
"""

import random

from flask import Blueprint, jsonify, request

import generate
from stellarObjects import _db, adminEdits, editStore
from stellarObjects.config import SystemConfig

from .authz import audit, require_admin
from .common import ApiError
from .limiter import limiter
from .routes import WRITE_RATE_LIMIT, _resolve_requested_write_db_config, _write_conn

bp = Blueprint("edits", __name__, url_prefix="/api")

_BODY_TABLES = {"planet": "planets", "moon": "moons", "belt": "asteroid_belts"}
_FACILITY_COLUMNS = {"planet": "planet_id", "moon": "moon_id", "belt": "asteroid_belt_id"}


def _options():
    """The optional JSON body: `{"drop_facilities": bool}`."""
    body = request.get_json(silent=True)
    if body is None:
        body = {}
    if not isinstance(body, dict):
        raise ApiError("request body must be a JSON object")
    unknown = set(body) - {"drop_facilities"}
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")
    drop = body.get("drop_facilities", False)
    if not isinstance(drop, bool):
        raise ApiError("'drop_facilities' must be a boolean")
    return drop


def _facilities_lost(conn, kind, body, regenerate):
    """How many facilities the edit would delete: those on the body and
    its moons for a delete; those on a regenerated planet's moons (its
    moons are replaced; the planet's own row stays)."""
    ids = []
    if kind == "planet":
        ids = [("moon_id", m.db_id) for m in body.moons if getattr(m, "db_id", None) is not None]
        if not regenerate:
            ids.append(("planet_id", body.db_id))
    elif not regenerate:
        ids = [(_FACILITY_COLUMNS[kind], body.db_id)]
    return sum(conn.execute(f"SELECT COUNT(*) AS n FROM facilities WHERE {column} = ?", (value,)).fetchone()["n"]
               for column, value in ids)


def _result_json(result, **extra):
    return jsonify({"status": "ok", "summary": result.summary, "moved": result.moved,
                    "reclassified": result.reclassified, "removed": result.removed,
                    "warnings": result.warnings, **extra})


def _edit_body(kind, body_id, regenerate):
    """Shared body of the planet/moon/belt delete and regenerate routes."""
    drop_facilities = _options()
    conn = _write_conn()
    try:
        with conn:
            row = conn.execute(f"SELECT star_system_id FROM {_BODY_TABLES[kind]} WHERE id = ?",
                               (body_id,)).fetchone()
            if row is None:
                raise ApiError(f"no such {kind}: {body_id}", status_code=404)
            system_id = row["star_system_id"]
            system = _db.load_star_system(conn, system_id)
            body, owner = adminEdits.find_body(system, kind, body_id)
            lost = _facilities_lost(conn, kind, body, regenerate)
            if lost and not drop_facilities:
                raise ApiError(f"this would delete {lost} facilit{'y' if lost == 1 else 'ies'}; send "
                               "\"drop_facilities\": true to go ahead", status_code=409)
            try:
                if regenerate:
                    result = adminEdits.regenerate_body(system, kind, body, owner)
                else:
                    result = adminEdits.remove_body(system, body, owner)
            except ValueError as exc:
                raise ApiError(str(exc), status_code=409)
            editStore.save_system_edits(conn, system_id, system)
    finally:
        conn.close()
    audit(f"{kind}.{'regenerate' if regenerate else 'delete'}", target=f"{kind}:{body_id}",
          detail=result.summary)
    return _result_json(result, star_system_id=system_id)


def _body_routes(kind, plural):
    """Registers `DELETE /api/<plural>/<id>` and `POST
    /api/<plural>/<id>/regenerate` for one body kind."""
    def delete(body_id):
        return _edit_body(kind, body_id, regenerate=False)

    def regenerate(body_id):
        return _edit_body(kind, body_id, regenerate=True)

    delete.__name__ = f"delete_{kind}"
    regenerate.__name__ = f"regenerate_{kind}"
    delete.__doc__ = (f"`DELETE /api/{plural}/<id>` -- removes one {kind} from its system "
                      "(optional body `{\"drop_facilities\": true}`; 409 without it when facilities "
                      "would go too).")
    regenerate.__doc__ = (f"`POST /api/{plural}/<id>/regenerate` -- rolls the {kind} again at the same "
                          "orbit, keeping its name and row, then re-validates the system.")
    bp.route(f"/{plural}/<int:body_id>", methods=["DELETE"])(
        limiter.limit(WRITE_RATE_LIMIT)(require_admin(fresh=True)(delete)))
    bp.route(f"/{plural}/<int:body_id>/regenerate", methods=["POST"])(
        limiter.limit(WRITE_RATE_LIMIT)(require_admin(fresh=True)(regenerate)))


_body_routes("planet", "planets")
_body_routes("moon", "moons")
_body_routes("belt", "belts")


# ---------------------------------------------------------------------
# Phenomena
# ---------------------------------------------------------------------

@bp.route("/phenomena/<phenomenon_type>/<int:phenomenon_id>", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def delete_phenomenon(phenomenon_type, phenomenon_id):
    """`DELETE /api/phenomena/<type>/<id>` -- deletes one standalone
    phenomenon (a system's own black hole or neutron star goes with its
    system instead: 409)."""
    _options()
    conn = _write_conn()
    try:
        with conn:
            try:
                deleted = editStore.delete_phenomenon(conn, phenomenon_type, phenomenon_id)
            except editStore.EditError as exc:
                raise ApiError(str(exc), status_code=409 if phenomenon_type in editStore.PHENOMENON_TABLES else 404)
    finally:
        conn.close()
    if not deleted:
        raise ApiError(f"no such {phenomenon_type}: {phenomenon_id}", status_code=404)
    audit("phenomenon.delete", target=f"{phenomenon_type}:{phenomenon_id}")
    return jsonify({"status": "ok", "summary": "Deleted."})


@bp.route("/phenomena/<phenomenon_type>/<int:phenomenon_id>/regenerate", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def regenerate_phenomenon(phenomenon_type, phenomenon_id):
    """`POST /api/phenomena/<type>/<id>/regenerate` -- rolls a standalone
    phenomenon again as the same type, keeping its id, name, sector,
    galaxy position and what it sits inside."""
    _options()
    if phenomenon_type not in editStore.PHENOMENON_TABLES:
        raise ApiError(f"unknown phenomenon type: {phenomenon_type}", status_code=404)
    conn = _write_conn()
    try:
        with conn:
            row = editStore.phenomenon_row(conn, phenomenon_type, phenomenon_id)
            if row is None:
                raise ApiError(f"no such {phenomenon_type}: {phenomenon_id}", status_code=404)
            fresh = generate.generate_phenomenon(editStore.GENERATOR_TYPES[phenomenon_type], SystemConfig(),
                                                 anchor_system=False, name=row["name"])
            try:
                editStore.replace_phenomenon_content(conn, phenomenon_type, phenomenon_id, fresh)
            except editStore.EditError as exc:
                raise ApiError(str(exc), status_code=409)
    finally:
        conn.close()
    audit("phenomenon.regenerate", target=f"{phenomenon_type}:{phenomenon_id}")
    return jsonify({"status": "ok", "summary": f"Regenerated {row['name']}."})


# ---------------------------------------------------------------------
# Sectors
# ---------------------------------------------------------------------

@bp.route("/sectors/<int:sector_id>/contents", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def delete_sector_with_contents(sector_id):
    """`DELETE /api/sectors/<id>/contents` -- deletes the sector together
    with its star systems, the phenomena filed under it and facilities
    parked in it, leaving its slot in the galaxy unfilled so it can be
    generated again. (`DELETE /api/sectors/<id>` deletes only the sector
    row and keeps its systems as standalone ones.)"""
    _options()
    conn = _write_conn()
    try:
        with conn:
            counts = editStore.delete_sector_with_contents(conn, sector_id)
    finally:
        conn.close()
    if counts is None:
        raise ApiError(f"no such sector: {sector_id}", status_code=404)
    audit("sector.delete", target=f"sector:{sector_id}", detail=f"with contents {counts}")
    return jsonify({"status": "ok", **counts})


@bp.route("/sectors/<int:sector_id>/regenerate", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def regenerate_sector(sector_id):
    """`POST /api/sectors/<id>/regenerate` -- deletes a galaxy-placed
    sector with everything in it (as `DELETE .../contents`) and generates
    its slot again from the galaxy's density plan. The new sector gets a
    new id and name (`sector_id` in the answer; `null` if the slot no
    longer qualifies for a sector). 409 for a sector off the galaxy grid
    or before `generate.py plan` has run."""
    _options()
    config = _resolve_requested_write_db_config()
    conn = _write_conn()
    try:
        with conn:
            try:
                address = editStore.sector_address(conn, sector_id)
            except editStore.EditError:
                raise ApiError(f"no such sector: {sector_id}", status_code=404)
            if address is None:
                raise ApiError("this sector isn't placed in the galaxy, so it can't be generated again",
                               status_code=409)
            if _db.get_galaxy_shape(conn) is None:
                raise ApiError("the galaxy has no density plan yet; run 'generate.py plan' first", status_code=409)
            counts = editStore.delete_sector_with_contents(conn, sector_id)
    finally:
        conn.close()
    random.seed()
    result = generate.ensure_sector_generated(*address, config=config)
    audit("sector.regenerate", target=f"sector:{sector_id}",
          detail=f"deleted {counts}; new sector {result['sector_id']}")
    return jsonify({"status": "ok", "deleted": counts, "sector_id": result["sector_id"],
                    "sector_name": result["sector_name"]})
