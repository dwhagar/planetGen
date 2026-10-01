# html/web/system_facilities.py

"""
The system page's Facilities panel (`/system/<id>`, `web/system_pages.py`,
schema v42): the system's starbases, colonies and outposts, and for an
admin a form to place a new one and a Remove button on each.

The form is a plain POST to the system page carrying `csrf_field()` and a
`facility_action`:

- `preview`: checks the choice against the placement rules
  (`stellarObjects.facilities.check_facility`, the same check
  `POST /api/facilities` runs) and, for an orbital facility, shows the
  orbit `GET /api/facilities/orbit` works out from the host's mass, then
  re-renders the page with the form still filled in. Nothing is saved.
- `save`: `POST /api/facilities`. Success redirects (303) back to the GET
  page with `?facility=added`; a refusal re-renders the page with the
  API's reason next to the form.
- `remove`: `DELETE /api/facilities/<id>` for one of this system's own
  facilities, then the same redirect with `?facility=removed`. The button
  sits behind a `<details>` "Remove" disclosure, the confirm step, so it
  needs no script.
"""

import math
import re

from flask import abort, redirect, request

import apiclient
from fmt import format_distance_km
from systempage import facility_kind_label, facility_row
from tabledisplay import format_period

from stellarObjects import facilities as facility_rules
from stellarObjects.physical_constants import AU_TO_KM
from stellarObjects.program_constants import FACILITY_KINDS

from .helpers import current_admin, db_name, page_url, trusted_html

ACTIONS = ("preview", "save", "remove")
"""tuple: The `facility_action` values the system page's POST accepts."""

PLACEMENT_OPTIONS = (
    ("orbital", "In orbit around it"),
    ("terrestrial", "On its surface (a planet or moon)"),
    ("asteroid", "In the asteroid belt"),
)
"""tuple: `(placement, label)` for the form. Stand-alone facilities park
in a sector, not a system, so they aren't offered here."""

DISTANCE_UNITS = {"km": 1.0, "au": AU_TO_KM}
"""dict: The orbital distance field's units, in km."""

MESSAGES = {
    # ?facility=<code> after a save or remove -> (is_error, message)
    "added": (False, "Facility added."),
    "removed": (False, "Facility removed."),
    "gone": (True, "That facility was already removed."),
}
"""dict: Fixed messages for the redirect after a save or remove, so
nothing from the query string reaches the page as text."""

_BODY_TYPES = {"t": "terrestrial", "g": "gas giant"}


def _api_message(exc):
    """An `ApiError`'s message without `apiclient`'s "planetGen API error
    (NNN): " prefix."""
    message = re.sub(r"^planetGen API error \(\d+\): ", "", str(exc))
    if getattr(exc, "status_code", None) in (401, 403):
        return "Not allowed. Change the admin username and password the installer set first."
    return message


def host_options(system):
    """
    Every body in the system a facility can go on or around, as the host
    `<select>`'s options: `{"value": "<host_type>:<id>", "label", "type",
    "id", "body_type"}` -- its star(s), each planet with its moons right
    after it, then each asteroid belt.
    """
    options = [{"value": f"star:{star['id']}", "label": f"Star: {star['name']}", "type": "star",
                "id": star["id"], "body_type": None} for star in system["stars"]]
    for planet in system["planets"]:
        options.append({
            "value": f"planet:{planet['id']}",
            "label": f"Planet: {planet['name']} ({_BODY_TYPES.get(planet['body_type'], planet['body_type'])})",
            "type": "planet", "id": planet["id"], "body_type": planet["body_type"],
        })
        for moon in planet.get("moons") or []:
            options.append({
                "value": f"moon:{moon['id']}",
                "label": f"Moon: {moon['name']} of {planet['name']} "
                         f"({_BODY_TYPES.get(moon['body_type'], moon['body_type'])})",
                "type": "moon", "id": moon["id"], "body_type": moon["body_type"],
            })
    for belt in system["belts"]:
        options.append({
            "value": f"asteroid_belt:{belt['id']}",
            "label": f"Asteroid belt at {format_distance_km(belt['distance_km'])}",
            "type": "asteroid_belt", "id": belt["id"], "body_type": None,
        })
    return options


def kind_options():
    """`(kind, label)` for every `program_constants.FACILITY_KINDS` key,
    with its one-line meaning."""
    return [(kind, f"{facility_kind_label(kind)}: {meaning}") for kind, meaning in FACILITY_KINDS.items()]


def _host_label(facility):
    """Where a facility is, by its host's name: "Sol (star)"."""
    kind = facility["host_type"].replace("_", " ")
    if facility.get("host_name"):
        return f"{facility['host_name']} ({kind})"
    return kind.capitalize()


def facility_rows(facilities):
    """The Facilities table's rows: `id`, `name`, `kind`, `host`,
    `placement`, `distance`, `period`, `speed` (the HTML ones already
    escaped by `systempage.facility_row`, marked trusted here)."""
    rows = []
    for facility in facilities:
        row = {key: (trusted_html(value) if value is not None else None)
               for key, value in facility_row(facility).items()}
        row["id"] = facility["id"]
        row["host"] = _host_label(facility)
        rows.append(row)
    return rows


def _read_form(system):
    """
    The submitted form's values, the host it names, the distance in km and
    every problem found before asking the API anything.

    Returns:
        tuple: `(values, host, distance_km, errors)`.
    """
    values = {key: (request.form.get(key) or "").strip()
              for key in ("name", "host", "placement", "kind", "distance", "distance_unit", "description")}
    if values["distance_unit"] not in DISTANCE_UNITS:
        values["distance_unit"] = "km"
    errors = []
    host = next((option for option in host_options(system) if option["value"] == values["host"]), None)
    if host is None:
        errors.append("Choose where the facility goes.")
    if values["placement"] not in dict(PLACEMENT_OPTIONS):
        errors.append("Choose how it is placed: in orbit, on the surface or in the belt.")
    if values["kind"] not in FACILITY_KINDS:
        errors.append("Choose the kind of facility.")

    distance_km = None
    if values["placement"] == "orbital" and values["distance"]:
        try:
            distance_km = float(values["distance"].replace(",", "")) * DISTANCE_UNITS[values["distance_unit"]]
        except ValueError:
            distance_km = None
        if distance_km is None or not math.isfinite(distance_km) or distance_km <= 0:
            errors.append("The orbital distance must be a number above zero.")
            distance_km = None

    if not errors:
        problem = facility_rules.check_facility(values["kind"], values["placement"], host["type"], host["body_type"])
        if problem:
            errors.append(f"Not allowed: {problem}.")
    return values, host, distance_km, errors


def _preview_orbit(host, distance_km):
    """`GET /api/facilities/orbit` for the form's host and distance, as
    display values, or `(None, error)`."""
    try:
        orbit = apiclient.get_facility_orbit(db_name(), host["type"], host["id"], distance_km)
    except apiclient.NotFoundError as exc:
        return None, str(exc)
    except apiclient.ApiError as exc:
        return None, f"Not allowed: {_api_message(exc)}."
    return {
        "distance": trusted_html(format_distance_km(orbit["distance_km"])),
        "period": format_period(orbit["period_years"]),
        "speed": f'{orbit["orbital_speed_kms"]:,.2f} km/s',
    }, None


def handle_post(system_id, system, facilities):
    """
    Runs one `facility_action` POST (the CSRF token was already checked
    app-wide). Returns a redirect response after a save or remove, else
    the form state for the page to show: `{"fields", "errors",
    "preview"}`.
    """
    if current_admin() is None:
        abort(403)
    action = request.form.get("facility_action")
    cookie_header = request.headers.get("Cookie")
    done = lambda code: redirect(  # noqa: E731
        page_url("system", system_id=system_id, facility=code, _anchor="facilities"), code=303)

    if action == "remove":
        try:
            facility_id = int(request.form.get("facility_id", ""))
        except ValueError:
            facility_id = None
        if facility_id not in {facility["id"] for facility in facilities}:
            return {"fields": {}, "errors": ["That facility is not in this system."], "preview": None}
        try:
            apiclient.delete_facility(cookie_header, db_name(), facility_id)
        except apiclient.NotFoundError:
            return done("gone")
        except apiclient.ApiError as exc:
            return {"fields": {}, "errors": [_api_message(exc)], "preview": None}
        return done("removed")

    values, host, distance_km, errors = _read_form(system)
    if action == "save" and not values["name"]:
        errors.insert(0, "Give the facility a name.")
    state = {"fields": values, "errors": errors, "preview": None}
    if errors:
        return state

    if action == "preview":
        if values["placement"] == "orbital":
            state["preview"], error = _preview_orbit(host, distance_km)
            if error:
                errors.append(error)
        else:
            state["preview"] = {"distance": None}
        return state

    body = {"name": values["name"], "kind": values["kind"], "placement": values["placement"],
            "host_type": host["type"], "host_id": host["id"]}
    if distance_km is not None:
        body["distance_km"] = distance_km
    if values["description"]:
        body["description"] = values["description"]
    try:
        apiclient.create_facility(cookie_header, db_name(), body)
    except (apiclient.NotFoundError, apiclient.ApiError) as exc:
        errors.append(_api_message(exc))
        return state
    return done("added")
