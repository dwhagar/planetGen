# planetgen/web/facility_type_pages.py

"""
The admin Facility types page (POP.10): `/admin/facility-types` lists every
facility type with its typed fields and lets an admin add, change, retire and
delete types and their fields. Each form is a plain POST carrying
`csrf_field()` and a `type_action`; the page calls the API (`api/facility_types.py`,
which holds the rules), then redirects (303) back to the GET page, which shows
what happened once through the admin flash cookie. A refusal shows the API's
reason in the same way.
"""

from flask import abort, request

from planetgen.population.facility_icons import ICON_NAMES
from planetgen.tuning import FACILITY_KINDS
from planetgen.population import facility_types
from planetgen.web.lib import apiclient

from . import bp
from .admin_pages import _api_message, _cookie_header, _flash, _render, _require_admin, _see_other, _take_flash
from .helpers import crumb, db_name, page_url

DATATYPE_LABELS = (("int", "Whole number"), ("decimal", "Decimal number"), ("text", "Text"), ("bool", "Yes or no"),
                   ("choice", "One of a list"))
"""tuple: `(datatype, label)` for a field's datatype list."""

ACTIONS = ("create_type", "update_type", "delete_type", "add_field", "update_field", "delete_field", "examples")
"""tuple: The `type_action` values the page's POST accepts."""


def _text(name):
    return (request.form.get(name) or "").strip()


def _whole(name, default=0):
    raw = _text(name)
    if not raw:
        return default
    try:
        return int(raw)
    except ValueError:
        raise ValueError(f"{name.replace('_', ' ').capitalize()} must be a whole number.") from None


def _number(name):
    raw = _text(name)
    if not raw:
        return None
    try:
        return float(raw)
    except ValueError:
        raise ValueError(f"{name.replace('_', ' ').capitalize()} must be a number.") from None


def _type_body():
    """The `POST`/`PATCH /api/facility-types` body of the type form."""
    return {"name": _text("name"), "kind": _text("kind"), "icon": _text("icon") or None,
            "description": _text("description") or None, "sort_order": _whole("sort_order"),
            "active": request.form.get("active") == "1"}


def _field_body(updating):
    """The body of the field form: a blank optional value is sent as `null` when changing a field (so it clears)
    and left out when adding one."""
    body = {"label": _text("label"), "datatype": _text("datatype"), "unit": _text("unit") or None,
            "required": request.form.get("required") == "1", "min_value": _number("min_value"),
            "max_value": _number("max_value"), "default_text": _text("default_text") or None,
            "sort_order": _whole("sort_order")}
    options = [line.strip() for line in request.form.get("choices", "").splitlines() if line.strip()]
    body["choices"] = options if body["datatype"] == "choice" else None
    if not updating:
        body["key"] = _text("key")
        body = {key: value for key, value in body.items() if value is not None}
    return body


def _do(action, cookie_header, db):
    """Runs one action; returns the message to flash. Raises `ValueError` (a form problem) or `ApiError`."""
    write = lambda method, path, body=None: apiclient.facility_type_write(method, cookie_header, db, path, body)  # noqa: E731
    type_id = _whole("type_id", default=None) if action not in ("create_type", "examples") else None
    if action == "create_type":
        created = write("POST", "", _type_body())
        return f"Added the type {_text('name')}.", f"type-{created['id']}"
    if action == "examples":
        added = write("POST", "/examples")["added"]
        return (f"Added {', '.join(added)}." if added else "The example type is there already."), None
    if type_id is None:
        raise ValueError("Choose a facility type.")
    if action == "update_type":
        write("PATCH", f"/{type_id}", _type_body())
        return f"Saved {_text('name')}.", f"type-{type_id}"
    if action == "delete_type":
        write("DELETE", f"/{type_id}")
        return "Deleted the type.", None
    if action == "add_field":
        write("POST", f"/{type_id}/fields", _field_body(updating=False))
        return f"Added the field {_text('key')}.", f"type-{type_id}"
    field_id = _whole("field_id", default=None)
    if field_id is None:
        raise ValueError("Choose a field.")
    if action == "update_field":
        write("PATCH", f"/{type_id}/fields/{field_id}", _field_body(updating=True))
        return "Saved the field.", f"type-{type_id}"
    write("DELETE", f"/{type_id}/fields/{field_id}")
    return "Deleted the field.", f"type-{type_id}"


@bp.route("/admin/facility-types", methods=["GET", "POST"])
def admin_facility_types():
    """The Facility types page: every type with its fields, and the forms that change them."""
    identity, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    db = db_name()

    if request.method == "POST":
        action = request.form.get("type_action")
        if action not in ACTIONS:
            abort(400)
        try:
            message, anchor = _do(action, cookie_header, db)
        except ValueError as exc:
            return _flash(_see_other(page_url("admin_facility_types")), error=str(exc))
        except (apiclient.ApiError, apiclient.NotFoundError) as exc:
            return _flash(_see_other(page_url("admin_facility_types")), error=_api_message(exc))
        return _flash(_see_other(page_url("admin_facility_types", _anchor=anchor)), message=message)

    flashed = _take_flash()
    return _render(
        "admin_facility_types.html", flashed=True, title="Facility types", section="admin_facility_types",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Facility types")], admin=identity,
        types=apiclient.get_facility_types(db, active_only=False), kinds=list(FACILITY_KINDS),
        icons=ICON_NAMES, datatypes=DATATYPE_LABELS, max_fields=facility_types.MAX_FIELDS,
        message=flashed.get("message"), error=flashed.get("error"),
    )
