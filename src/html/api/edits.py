# html/api/edits.py

"""
Admin editing endpoints (TODO ADM.1): delete and regenerate one planet,
moon, asteroid belt, phenomenon or sector (ADM.8), change a planet's or
moon's class (ADM.6) and change a single-star system's star (ADM.7).
Every write needs an admin whose credentials are current, writes an
audit-log row, and answers with what the edit did.

A body edit loads its system (`store.load_star_system`), changes it with
`planetgen.admin.edits`, re-validates it from the moons outward
(`planetgen.generation.validation`) and writes it back in place
(`planetgen.db.edits`), so the other bodies keep their rows. The
answer lists the bodies the re-validation moved and any warnings left.

Kept in its own blueprint, like `population.py`, so this work doesn't
touch `routes.py`.
"""

import random

from flask import Blueprint, jsonify, request

import generate
from planetgen.db import edits as editStore, store
from planetgen.admin import edits as adminEdits
from planetgen.generation import validation
from planetgen import tuning
from planetgen.generation.config import SystemConfig

from .authz import audit, require_admin
from .common import ApiError
from .limiter import limiter
from .routes import WRITE_RATE_LIMIT, _resolve_requested_write_db_config, _write_conn

bp = Blueprint("edits", __name__, url_prefix="/api")

_BODY_TABLES = {"planet": "planets", "moon": "moons", "belt": "asteroid_belts"}
_FACILITY_COLUMNS = {"planet": "planet_id", "moon": "moon_id", "belt": "asteroid_belt_id"}


def _json_object(allowed):
    """The optional JSON body, checked for unknown fields."""
    body = request.get_json(silent=True)
    if body is None:
        body = {}
    if not isinstance(body, dict):
        raise ApiError("request body must be a JSON object")
    unknown = set(body) - set(allowed)
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")
    return body


def _options():
    """The optional JSON body: `{"drop_facilities": bool}`."""
    body = _json_object({"drop_facilities"})
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
            system = store.load_star_system(conn, system_id)
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
# ADM.6 and ADM.7: a body's class, a system's star
# ---------------------------------------------------------------------

def _load_system_of(conn, kind, body_id):
    row = conn.execute(f"SELECT star_system_id FROM {_BODY_TABLES[kind]} WHERE id = ?", (body_id,)).fetchone()
    if row is None:
        raise ApiError(f"no such {kind}: {body_id}", status_code=404)
    return row["star_system_id"], store.load_star_system(conn, row["star_system_id"])


@bp.route("/systems/<int:system_id>/class-options")
@require_admin()
def system_class_options(system_id):
    """`GET /api/systems/<id>/class-options` -- for each planet and moon
    (`"planet:<id>"`, `"moon:<id>"`), the classes it could take without
    moving anything (`recommended`, most common first), plus every class
    code (`all`) for a forced change."""
    conn = _write_conn()
    try:
        try:
            system = store.load_star_system(conn, system_id)
        except ValueError:
            raise ApiError(f"no such system: {system_id}", status_code=404)
        options = adminEdits.class_options(system)
    finally:
        conn.close()
    return jsonify({"recommended": options, "all": sorted(tuning.PLANET_CLASSES)})


def _change_class(kind, body_id):
    body = _json_object({"class", "force"})
    planet_class = body.get("class")
    force = body.get("force", False)
    if not isinstance(planet_class, str) or planet_class not in tuning.PLANET_CLASSES:
        raise ApiError(f"'class' must be one of {', '.join(sorted(tuning.PLANET_CLASSES))}")
    if not isinstance(force, bool):
        raise ApiError("'force' must be a boolean")
    conn = _write_conn()
    try:
        with conn:
            system_id, system = _load_system_of(conn, kind, body_id)
            target, owner = adminEdits.find_body(system, kind, body_id)
            try:
                result = adminEdits.change_class(system, target, owner, planet_class, force=force)
            except ValueError as exc:
                raise ApiError(str(exc), status_code=409)
            editStore.save_system_edits(conn, system_id, system)
    finally:
        conn.close()
    audit(f"{kind}.class", target=f"{kind}:{body_id}", detail=f"{planet_class} force={force}: {result.summary}")
    return _result_json(result, star_system_id=system_id)


@bp.route("/planets/<int:body_id>/class", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def change_planet_class(body_id):
    """`POST /api/planets/<id>/class` `{"class": "M", "force": false}` --
    changes a planet's class (ADM.6), re-generating its surface
    conditions as that class and keeping its orbit, mass and name
    (ADM.27). A class its mass doesn't fit is refused (409) and nothing
    changes. Without `force` only a recommended class is accepted (409
    otherwise); with it any class its mass fits. The rest of the
    system is re-spaced to fit, never trimmed; what still doesn't validate
    comes back in `warnings`, and the change is saved either way."""
    return _change_class("planet", body_id)


@bp.route("/moons/<int:body_id>/class", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def change_moon_class(body_id):
    """`POST /api/moons/<id>/class` -- as for a planet, for a moon."""
    return _change_class("moon", body_id)


def _removed_rows(system_before_ids, system):
    """The planet, moon and belt rows an edit dropped:
    `[(facility column, row id)]`."""
    kept = set()
    for _star, planets in validation.star_lists(system):
        for body in planets:
            kind = "belt" if body.body_type == 'a' else "planet"
            kept.add((kind, getattr(body, "db_id", None)))
            for moon in getattr(body, "moons", ()):
                kept.add(("moon", getattr(moon, "db_id", None)))
    return [(_FACILITY_COLUMNS[kind], row_id) for kind, row_id in system_before_ids if (kind, row_id) not in kept]


def _body_ids(system):
    ids = []
    for _star, planets in validation.star_lists(system):
        for body in planets:
            ids.append(("belt" if body.body_type == 'a' else "planet", body.db_id))
            for moon in getattr(body, "moons", ()):
                ids.append(("moon", moon.db_id))
    return ids


@bp.route("/systems/<int:system_id>/star", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def change_system_star(system_id):
    """`POST /api/systems/<id>/star` `{"star_type": "K2V"}` -- replaces a
    single-star system's star with a new one of that type (ADM.7). Every
    planet, moon and belt keeps its class; orbits scale with the new
    star's light and are re-spaced, and bodies past its farthest stable
    orbit (or moons a planet moved inward can no longer hold) are removed
    and listed in `removed`. 409 for a binary, a black hole or neutron
    star, a system built around a pre-placed bright star, or when removed
    bodies host facilities and `"drop_facilities": true` isn't sent."""
    body = _json_object({"star_type", "drop_facilities"})
    star_type = body.get("star_type")
    drop_facilities = body.get("drop_facilities", False)
    if not isinstance(star_type, str):
        raise ApiError("'star_type' must be a spectral type such as \"K2V\"")
    if not isinstance(drop_facilities, bool):
        raise ApiError("'drop_facilities' must be a boolean")
    conn = _write_conn()
    try:
        with conn:
            try:
                system = store.load_star_system(conn, system_id)
            except ValueError:
                raise ApiError(f"no such system: {system_id}", status_code=404)
            if store.system_content_blockers(conn, system_id)["bright_star"]:
                raise ApiError("this system is built around a pre-placed bright star, so its star can't be changed",
                               status_code=409)
            before = _body_ids(system)
            try:
                result, new_star = adminEdits.change_star(system, star_type)
            except ValueError as exc:
                raise ApiError(str(exc), status_code=409 if "single star" in str(exc) or "black hole" in str(exc)
                               else 400)
            lost = sum(conn.execute(f"SELECT COUNT(*) AS n FROM facilities WHERE {column} = ?",
                                    (row_id,)).fetchone()["n"]
                       for column, row_id in _removed_rows(before, system))
            if lost and not drop_facilities:
                raise ApiError(f"this would delete {lost} facilit{'y' if lost == 1 else 'ies'} on the bodies "
                               "with no room left; send \"drop_facilities\": true to go ahead", status_code=409)
            editStore.save_system_edits(conn, system_id, system, stars=[new_star])
    finally:
        conn.close()
    audit("system.star", target=f"system:{system_id}", detail=f"{star_type}: {result.summary}")
    return _result_json(result, star_system_id=system_id, star_type=new_star.type.split()[0])


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
    new id and name (`sector_id` in the answer; `null` if the slot is
    outside the galaxy's outline). 409 for a sector off the galaxy grid
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
            if store.get_galaxy_shape(conn) is None:
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
