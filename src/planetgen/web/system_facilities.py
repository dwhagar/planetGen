# planetgen/web/system_facilities.py

"""
The system page's Facilities panel (`/system/<id>`, `web/system_pages.py`,
schema v42): the system's starbases, colonies and outposts, and for an
admin a form to place a new one and a Remove button on each.

The form is a plain POST to the system page carrying `csrf_field()` and a
`facility_action`:

- `preview`: checks the choice against the placement rules
  (`planetgen.population.facilities.check_facility`, the same check
  `POST /api/facilities` runs) and, for an orbital facility, shows the
  orbit `GET /api/facilities/orbit` works out from the host's mass, then
  re-renders the page with the form still filled in. Nothing is saved.

The form asks for the name, then the placement, then the host; the host
list holds only the hosts that take that placement
(`facilities.host_placements`; `static/facilityform.js` hides the rest,
and the server refuses a mismatch). An orbital facility's distance is a
logarithmic slider (`orbit_step`, `facilities.distance_from_step`) from
just above the host's surface to the edge of its sphere of influence
(`facilities.orbit_limits`), shown only for "in orbit". A facility on a
surface takes no distance, and one in a belt gets a random spot in it
when saved (`_db.add_facility`).
- `save`: `POST /api/facilities`. Success redirects (303) back to the GET
  page with `?facility=added`; a refusal re-renders the page with the
  API's reason next to the form.
- `remove`: `DELETE /api/facilities/<id>` for one of this system's own
  facilities, then the same redirect with `?facility=removed`. The button
  sits behind a `<details>` "Remove" disclosure, the confirm step, so it
  needs no script.
"""

import re

from flask import abort, redirect, request

from planetgen.web.lib import apiclient, facility_display
from planetgen.web.lib.fmt import format_distance_km, format_speed_kms
from planetgen.web.lib.systempage import facility_kind_label, facility_row
from planetgen.web.lib.tabledisplay import format_period

from planetgen.admin import activity_log
from planetgen.population import facilities as facility_rules, facility_ratings
from planetgen.tuning import FACILITY_KINDS, FACILITY_ORBIT_STEPS

from .helpers import current_admin, db_name, page_url, trusted_html

ACTIONS = ("preview", "save", "remove", "edit", "suggest")
"""tuple: The `facility_action` values the system page's POST accepts."""

PLACEMENT_OPTIONS = (
    ("orbital", "In orbit around it"),
    ("terrestrial", "On its surface (a planet or moon)"),
    ("asteroid", "In the asteroid belt"),
)
"""tuple: `(placement, label)` for the form. Stand-alone facilities park
in a sector, not a system, so they aren't offered here."""

ORBIT_STEP_DEFAULT = FACILITY_ORBIT_STEPS // 2
"""int: Where the orbit slider starts: halfway along its log scale, the
geometric mean of the host's surface and the edge of its sphere of
influence."""

RATING_CHOICES = (("", "Not assessed"), ("green", "Green"), ("yellow", "Yellow"), ("red", "Red"))
"""tuple: `(value, label)` for a rating's `<select>`."""

MESSAGES = {
    # ?facility=<code> after a save or remove -> (is_error, message)
    "added": (False, "Facility added."),
    "removed": (False, "Facility removed."),
    "gone": (True, "That facility was already removed."),
    "edited": (False, "Facility changed."),
    "suggested": (False, "Ratings suggested for the ones that were not set."),
    "nothing": (False, "Every rating was already set, so nothing was suggested."),
}
"""dict: Fixed messages for the redirect after a save or remove, so
nothing from the query string reaches the page as text."""

_BODY_TYPES = {"t": "terrestrial", "g": "gas giant"}

_WORST = ("green", "yellow", "red")


def _api_message(exc):
    """An `ApiError`'s message without `apiclient`'s "planetGen API error
    (NNN): " prefix."""
    message = re.sub(r"^planetGen API error \(\d+\): ", "", str(exc))
    if getattr(exc, "status_code", None) in (401, 403):
        return "Not allowed. Change the admin username and password the installer set first."
    return message


def _star_host(system, star):
    """`facilities.star_host` for one of `system`'s stars."""
    return facility_rules.star_host(
        system["stars"], star["id"], system.get("binary_configuration"),
        system.get("binary_separation_km"), system.get("binary_heliosphere_radius_km"))


def _orbit_fields(mass_kg, radius_km, sphere_km):
    """An orbitable host's `mass_kg`, `min_km` and `max_km` (the orbit
    slider's ends, `facilities.orbit_limits`), or `None`s when its size
    or mass is unknown."""
    if not radius_km or not mass_kg:
        return {"mass_kg": None, "min_km": None, "max_km": None}
    lowest, highest = facility_rules.orbit_limits(radius_km, sphere_km)
    return {"mass_kg": mass_kg, "min_km": lowest, "max_km": highest}


def host_options(system):
    """
    Every body in the system a facility can go on or around, as the host
    `<select>`'s options: `{"value": "<host_type>:<id>", "label", "type",
    "id", "body_type", "placements", "mass_kg", "min_km", "max_km"}` --
    its star(s), each planet with its moons right after it, then each
    asteroid belt. `placements` are the placements the host takes
    (`facilities.host_placements`); `mass_kg`, `min_km` and `max_km` give
    an orbitable host's slider range (`None` for a belt).
    """
    def option(host_type, host_id, label, body_type, orbit):
        return {"value": f"{host_type}:{host_id}", "label": label, "type": host_type, "id": host_id,
                "body_type": body_type, "placements": facility_rules.host_placements(host_type, body_type),
                **orbit}

    options = [option("star", star["id"], f"Star: {star['name']}", None, _orbit_fields(*_star_host(system, star)))
               for star in system["stars"]]
    for planet in system["planets"]:
        options.append(option(
            "planet", planet["id"],
            f"Planet: {planet['name']} ({_BODY_TYPES.get(planet['body_type'], planet['body_type'])})",
            planet["body_type"],
            _orbit_fields(planet.get("mass_kg"), planet.get("radius_km"), planet.get("hill_radius_km"))))
        for moon in planet.get("moons") or []:
            options.append(option(
                "moon", moon["id"],
                f"Moon: {moon['name']} of {planet['name']} "
                f"({_BODY_TYPES.get(moon['body_type'], moon['body_type'])})",
                moon["body_type"],
                _orbit_fields(moon.get("mass_kg"), moon.get("radius_km"), moon.get("hill_radius_km"))))
    for belt in system["belts"]:
        options.append(option(
            "asteroid_belt", belt["id"],
            f"Asteroid belt from {format_distance_km(belt['lower_limit_km'])}"
            f" to {format_distance_km(belt['upper_limit_km'])}",
            None, {"mass_kg": None, "min_km": None, "max_km": None}))
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
    """The Facilities table's rows, the worst-rated first (POP.10): `id`, `name`, `kind`, `host`,
    `placement`, `distance`, `period`, `speed` (the HTML ones already
    escaped by `systempage.facility_row`, marked trusted here), and the type
    (`type`, with its icon), `affiliation`, the rating `chips` and the field
    values as `fields`."""
    rows = []
    for facility in sorted(facilities, key=lambda f: -_WORST.index(f["overall"]) if f.get("overall") else 1):
        row = {key: (trusted_html(value) if value is not None else None)
               for key, value in facility_row(facility).items()}
        row["id"] = facility["id"]
        row["host"] = _host_label(facility)
        row["type"] = trusted_html(facility_display.type_html(facility))
        row["affiliation"] = facility_display.affiliation_text(facility)
        row["chips"] = trusted_html(facility_display.rating_chips(facility))
        row["fields"] = facility_display.fields_text(facility)
        # What the edit form starts from.
        row["edit"] = {
            "affiliation": f"polity:{facility['affiliation_polity_id']}" if facility.get("affiliation_polity_id") else "",
            "affiliation_name": None if facility.get("affiliation_polity_id") else facility.get("affiliation_name"),
            "ratings_note": facility.get("ratings_note"),
            **{f"rating_{name}": facility["ratings"][name] or "" for name in facility_ratings.RATINGS},
        }
        row["polity"] = ({"id": facility["affiliation_polity_id"], "name": facility["affiliation_name"]}
                         if facility.get("affiliation_polity_id") else None)
        rows.append(row)
    return rows


def form_context(system_id):
    """
    What the placement and edit forms need beyond the system itself: `types` (the active facility types with their
    fields, `GET /api/facility-types`) and `owner` (the polity that owns the system, whose name a new facility is
    affiliated with by default, or `None`).
    """
    try:
        owner = apiclient.get_system_owner(db_name(), system_id)
    except apiclient.ApiError:
        owner = None
    return {"types": apiclient.get_facility_types(db_name()), "owner": owner}


def _convert_field(field, raw):
    """One field's text from the form as the API takes it, or raises `ValueError` with the problem."""
    if field["datatype"] == "int":
        try:
            return int(raw)
        except ValueError:
            raise ValueError(f"{field['label']} must be a whole number.") from None
    if field["datatype"] == "decimal":
        try:
            return float(raw)
        except ValueError:
            raise ValueError(f"{field['label']} must be a number.") from None
    if field["datatype"] == "bool":
        if raw not in ("yes", "no"):
            raise ValueError(f"{field['label']} must be yes or no.")
        return raw == "yes"
    return raw


def read_extras(types):
    """
    The type, field values, affiliation and ratings of a submitted facility form.

    Args:
        types (list): The active types, as `form_context` returns them.

    Returns:
        tuple: `(values, body, errors)` -- the raw values to show the form again, the keys to add to the API body,
            and the problems found before asking the API anything.
    """
    form = request.form
    values = {key: (form.get(key) or "").strip()
              for key in ("type_id", "affiliation", "affiliation_name", "ratings_note")}
    values.update({f"rating_{name}": (form.get(f"rating_{name}") or "").strip() for name in facility_ratings.RATINGS})
    values.update({key: value.strip() for key, value in form.items() if key.startswith("field_")})
    body, errors = {}, []
    chosen = None
    if values["type_id"]:
        chosen = next((t for t in types if str(t["id"]) == values["type_id"]), None)
        if chosen is None:
            errors.append("Choose one of the listed facility types.")
        else:
            body["type_id"] = chosen["id"]
            fields = {}
            for field in chosen["fields"]:
                raw = values.get(f"field_{field['key']}", "")
                if raw == "":
                    continue
                try:
                    fields[field["key"]] = _convert_field(field, raw)
                except ValueError as exc:
                    errors.append(str(exc))
            if fields:
                body["fields"] = fields
    if values["affiliation_name"]:
        body["affiliation_name"] = values["affiliation_name"]
    elif values["affiliation"].startswith("polity:") and values["affiliation"][7:].isdigit():
        body["affiliation_polity_id"] = int(values["affiliation"][7:])
    ratings = {name: values[f"rating_{name}"] for name in facility_ratings.RATINGS if values[f"rating_{name}"]}
    if any(color not in facility_ratings.COLORS for color in ratings.values()):
        errors.append("A rating is Green, Yellow or Red.")
    elif ratings:
        body["ratings"] = ratings
    if values["ratings_note"]:
        body["ratings_note"] = values["ratings_note"]
    return values, body, errors


def _read_form(system):
    """
    The submitted form's values, the host it names, the orbit distance in
    km (from the slider, for an orbital facility) and every problem found
    before asking the API anything.

    Returns:
        tuple: `(values, host, distance_km, errors)`.
    """
    values = {key: (request.form.get(key) or "").strip()
              for key in ("name", "host", "placement", "kind", "orbit_step", "description")}
    try:
        step = int(values["orbit_step"])
    except ValueError:
        step = ORBIT_STEP_DEFAULT
    values["orbit_step"] = min(max(step, 0), FACILITY_ORBIT_STEPS)
    errors = []
    if values["placement"] not in dict(PLACEMENT_OPTIONS):
        errors.append("Choose how it is placed: in orbit, on the surface or in the belt.")
    host = next((option for option in host_options(system) if option["value"] == values["host"]), None)
    if host is None:
        errors.append("Choose where the facility goes.")
    elif not errors and values["placement"] not in host["placements"]:
        errors.append(f"{host['label']} can't take a facility placed "
                      f"{dict(PLACEMENT_OPTIONS)[values['placement']].lower()}; choose another host.")
    if values["kind"] not in FACILITY_KINDS:
        errors.append("Choose the kind of facility.")

    distance_km = None
    if not errors and values["placement"] == "orbital" and host["min_km"] is not None:
        distance_km = facility_rules.distance_from_step(values["orbit_step"], host["min_km"], host["max_km"])

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
        "speed": format_speed_kms(orbit["orbital_speed_kms"]),
    }, None


def _edit_body():
    """The `PATCH /api/facilities/<id>` body of the edit form: the affiliation, the five ratings (a blank one
    clears it) and the note, replaced as a set."""
    form = request.form
    body = {"affiliation_polity_id": None, "affiliation_name": None,
            "ratings": {name: (form.get(f"rating_{name}") or "").strip() or None
                        for name in facility_ratings.RATINGS},
            "ratings_note": (form.get("ratings_note") or "").strip() or None}
    name, chosen = (form.get("affiliation_name") or "").strip(), (form.get("affiliation") or "").strip()
    if name:
        body["affiliation_name"] = name
    elif chosen.startswith("polity:") and chosen[7:].isdigit():
        body["affiliation_polity_id"] = int(chosen[7:])
    return body


def handle_post(system_id, system, facilities, context):
    """
    Runs one `facility_action` POST (the CSRF token was already checked
    app-wide). Returns a redirect response after a save, remove, edit or
    suggestion, else the form state for the page to show: `{"fields",
    "errors", "preview"}`. `context` is `form_context`'s answer.
    """
    if current_admin() is None:
        activity_log.event("AUTHZ", "admin.required", path=request.path)
        abort(403)
    action = request.form.get("facility_action")
    cookie_header = request.headers.get("Cookie")
    done = lambda code: redirect(  # noqa: E731
        page_url("system", system_id=system_id, facility=code, _anchor="facilities"), code=303)

    if action == "remove":
        facility_id = request.form.get("facility_id", "").strip().upper()
        if facility_id not in {str(facility["id"]) for facility in facilities}:
            return {"fields": {}, "errors": ["That facility is not in this system."], "preview": None}
        try:
            apiclient.delete_facility(cookie_header, db_name(), facility_id)
        except apiclient.NotFoundError:
            return done("gone")
        except apiclient.ApiError as exc:
            return {"fields": {}, "errors": [_api_message(exc)], "preview": None}
        return done("removed")

    if action in ("edit", "suggest"):
        facility_id = request.form.get("facility_id", "").strip().upper()
        if facility_id not in {str(facility["id"]) for facility in facilities}:
            return {"fields": {}, "errors": ["That facility is not in this system."], "preview": None}
        try:
            if action == "suggest":
                return done("suggested" if apiclient.suggest_facility_ratings(cookie_header, db_name(), facility_id)
                            else "nothing")
            apiclient.update_facility(cookie_header, db_name(), facility_id, _edit_body())
        except apiclient.NotFoundError:
            return done("gone")
        except apiclient.ApiError as exc:
            return {"fields": {}, "errors": [_api_message(exc)], "preview": None}
        return done("edited")

    values, host, distance_km, errors = _read_form(system)
    extras, extra_body, extra_errors = read_extras(context["types"])
    values.update(extras)
    errors.extend(extra_errors)
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
        elif host["type"] == "asteroid_belt":
            belt = next(belt for belt in system["belts"] if belt["id"] == host["id"])
            state["preview"] = {"distance": None, "belt": trusted_html(
                f"{format_distance_km(belt['lower_limit_km'])} and {format_distance_km(belt['upper_limit_km'])}")}
        else:
            state["preview"] = {"distance": None}
        return state

    body = {"name": values["name"], "kind": values["kind"], "placement": values["placement"],
            "host_type": host["type"], "host_id": host["id"]}
    if distance_km is not None:
        body["distance_km"] = distance_km
    if values["description"]:
        body["description"] = values["description"]
    body.update(extra_body)
    try:
        apiclient.create_facility(cookie_header, db_name(), body)
    except (apiclient.NotFoundError, apiclient.ApiError) as exc:
        errors.append(_api_message(exc))
        return state
    return done("added")
