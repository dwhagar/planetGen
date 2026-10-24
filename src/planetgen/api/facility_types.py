"""
Facility types and their typed fields (POP.10, schema v89): anyone may read
them, an admin defines them. See `planetgen.population.facility_types` for
the rules and `docs/design/habitability/population-and-politics.md`.

Kept in its own blueprint, like `population.py`, so this work doesn't touch
`routes.py`.
"""

from flask import Blueprint, jsonify, request

from planetgen import tuning
from planetgen.population import facility_types as types

from .authz import audit, require_admin
from .common import ApiError, require_json_body
from .limiter import limiter
from .routes import WRITE_RATE_LIMIT, _write_conn, get_db
from .schemas import (
    FacilityFieldCreate, FacilityFieldUpdate, FacilityTypeCreate, FacilityTypeUpdate, given, parse_body,
)

bp = Blueprint("facility_types", __name__, url_prefix="/api")


def _truth(raw, name):
    if raw is None:
        return False
    if raw.strip().lower() in ("1", "true", "yes"):
        return True
    if raw.strip().lower() in ("0", "false", "no"):
        return False
    raise ApiError(f"{name} must be true or false, got {raw!r}")


def _write(action):
    """Runs `action(conn)` in a transaction on a write connection; a broken rule is the answer's status."""
    conn = _write_conn()
    try:
        with conn:
            return action(conn)
    except types.FacilityTypeError as exc:
        raise ApiError(str(exc), status_code=exc.status)
    finally:
        conn.close()


@bp.route("/facility-types")
def list_facility_types():
    """`GET /api/facility-types[?kind=][&active=1]` -- every type in sort order, each with its `fields` and the
    number of `facilities` that use it. `active=1` leaves out the types an admin retired."""
    kind = request.args.get("kind") or None
    if kind is not None and kind not in tuning.FACILITY_KINDS:
        raise ApiError(f"kind must be one of: {', '.join(sorted(tuning.FACILITY_KINDS))}")
    items = types.list_types(get_db(), kind=kind, active_only=_truth(request.args.get("active"), "active"))
    return jsonify({"items": items, "icons": list(types.ICON_NAMES)})


@bp.route("/facility-types/<int:type_id>")
def facility_type(type_id):
    """`GET /api/facility-types/<id>` -- one type with its fields."""
    found = types.get_type(get_db(), type_id)
    if found is None:
        raise ApiError(f"no such facility type: {type_id}", status_code=404)
    return jsonify(found)


@bp.route("/facility-types", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def create_facility_type():
    """`POST /api/facility-types` `{"name", "kind"}` plus `icon`, `description`, `sort_order`, `active` and
    `fields` (a list of field definitions, as `POST /api/facility-types/<id>/fields` takes) -- returns `{"id"}`."""
    body = given(parse_body(FacilityTypeCreate, require_json_body()))
    fields = body.pop("fields", [])
    type_id = _write(lambda conn: types.create_type(conn, body, fields))
    audit("facility_type.create", target=f"facility_type:{type_id}", detail=str({"name": body["name"], "kind": body["kind"]}))
    return jsonify({"id": type_id}), 201


@bp.route("/facility-types/<int:type_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def update_facility_type(type_id):
    """`PATCH /api/facility-types/<id>` -- any of `name`, `kind`, `icon`, `description`, `sort_order`, `active`.
    409 for a name already taken, and for a new `kind` while facilities use the type."""
    body = given(parse_body(FacilityTypeUpdate, require_json_body()))
    def change(conn):
        types.update_type(conn, type_id, body)
        return types.get_type(conn, type_id)

    found = _write(change)
    audit("facility_type.update", target=f"facility_type:{type_id}", detail=str(sorted(body)))
    return jsonify(found)


@bp.route("/facility-types/<int:type_id>", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def delete_facility_type(type_id):
    """`DELETE /api/facility-types/<id>` -- 409 while facilities use the type (set it inactive instead)."""
    _write(lambda conn: types.delete_type(conn, type_id))
    audit("facility_type.delete", target=f"facility_type:{type_id}")
    return jsonify({"status": "ok"})


@bp.route("/facility-types/examples", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def install_example_types():
    """`POST /api/facility-types/examples` (no body) -- adds the example types that are not there yet; returns
    `{"added": [names]}`."""
    added = _write(types.install_examples)
    audit("facility_type.examples", detail=", ".join(added) or "none")
    return jsonify({"added": added}), 201 if added else 200


@bp.route("/facility-types/<int:type_id>/fields", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def add_facility_type_field(type_id):
    """`POST /api/facility-types/<id>/fields` `{"key", "label", "datatype"}` (`int`, `decimal`, `text`, `bool`
    or `choice`) plus `unit`, `required`, `min_value`, `max_value`, `choices`, `default_text`, `sort_order` --
    returns `{"id"}`. 409 for a key the type has or a 25th field."""
    body = given(parse_body(FacilityFieldCreate, require_json_body()))
    field_id = _write(lambda conn: types.add_field(conn, type_id, body))
    audit("facility_type.field_add", target=f"facility_type:{type_id}", detail=body["key"])
    return jsonify({"id": field_id}), 201


@bp.route("/facility-types/<int:type_id>/fields/<int:field_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def update_facility_type_field(type_id, field_id):
    """`PATCH /api/facility-types/<id>/fields/<field_id>` -- any field definition column but `key`. 409 when a
    change would orphan a value a facility already holds."""
    body = given(parse_body(FacilityFieldUpdate, require_json_body()))
    def change(conn):
        types.update_field(conn, type_id, field_id, body)
        return types.get_type(conn, type_id)

    found = _write(change)
    audit("facility_type.field_update", target=f"facility_type:{type_id}", detail=f"{field_id}: {sorted(body)}")
    return jsonify(found)


@bp.route("/facility-types/<int:type_id>/fields/<int:field_id>", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def delete_facility_type_field(type_id, field_id):
    """`DELETE /api/facility-types/<id>/fields/<field_id>` -- 409 while any facility holds a value for it."""
    _write(lambda conn: types.delete_field(conn, type_id, field_id))
    audit("facility_type.field_delete", target=f"facility_type:{type_id}", detail=str(field_id))
    return jsonify({"status": "ok"})
