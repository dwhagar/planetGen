# html/web/system_pages.py

"""
The system and phenomenon pages (moved from the CGI `system.py`,
`phenomena.py` and `phenomenon.py`):

- `/system/<id>`: one star system -- the system map, the expandable body
  list with each body's generated description, the Wikitext/Markdown
  code views (`?code=wikitext|markdown`), the Stars/Planets/Belts/Comets
  tables, the system's facilities, and for an admin an "Upload to Wiki"
  form (POST, CSRF-checked, then redirected back to the GET page) and the
  facility form and Remove buttons (`web/system_facilities.py`).
- `/phenomena`: every exotic phenomenon, paged with `?page=N`.
- `/phenomenon/<type>/<id>`: one phenomenon's data and its AU-scale
  diagram.

Nothing here reads a database name from the request (`db_name()`), and
every link is a plain GET link (`page_url`).
"""

from flask import abort, redirect, request
import apiclient
from fmt import (
    format_number, runaway_text,
    format_distance_km, format_distance_ly, format_distance_pc, linkify_location, nearest_neighbors_location,
    nearest_systems_html,
)
from pagination import fetch_page, parse_page
from phenomenonmap import render_phenomenon_map_panel
from phenomenonrender import render_phenomenon_view_panel, view_kind
from systemmap import render_system_map_panel
from systempage import stars_html, system_list_html
from tabledisplay import format_star_radius, to_plain_text

from classref import ROGUE_MASS_CLASS_NAMES
from stellarObjects import activitylog
from stellarObjects.program_constants import NEBULA_CLASSES

from . import bp
from . import system_facilities
from .class_pages import class_url
from .helpers import (
    crumb, current_admin, db_name, page_url, pager, population_status, render_page, trusted_html,
)
from .nav_page import endpoint, nav_url

# ---------------------------------------------------------------------
# /system/<id>
# ---------------------------------------------------------------------

CODE_FORMATS = {"wikitext": "Wikitext", "markdown": "Markdown"}
"""dict: `?code=` values the system page accepts, with their labels."""

WIKI_BACKENDS = {"wikijs": "Wiki.js", "mediawiki": "MediaWiki"}

WIKI_MESSAGES = {
    # ?wiki=<code> after an upload POST -> (is_error, message)
    "uploaded": (False, "Uploaded to the wiki."),
    "exists": (True, "A page already exists at that location."),
    "invalid": (True, "The upload was rejected. Wiki.js needs a path; check the backend and path."),
    "unconfigured": (True, "That wiki is not configured for this site."),
    "forbidden": (True, "Not allowed. Change the admin username and password the installer set first."),
    "failed": (True, "The upload failed. The wiki could not be reached or returned an error."),
}
"""dict: The upload outcome shown after the POST-redirect-GET. Only these
fixed messages are shown, so nothing from the query string reaches the
page as text."""

_LOCATION_MARKER = " -- nearest: "


def _system_url(system_id):
    """`fmt`'s location `system_url` hook: a plain GET link."""
    return page_url("system", system_id=system_id)


def nav_links(kind, entity_id):
    """
    "Navigate from here" / "Navigate to here" URLs for a system
    (`kind="system"`) or a phenomenon (`kind` = its type):
    `/nav?from=<kind>:<id>` and `/nav?to=<kind>:<id>` (`web/nav_page.py`).
    """
    point = endpoint(kind, entity_id)
    return {"from": nav_url(origin=point), "to": nav_url(destination=point)}


def _location_html(system):
    """The "Location:" line: the sector name and nearest neighbours, each
    linked, or `None` for a system with no stored location."""
    if not system["location"]:
        return None
    neighbors = system.get("nearest_neighbors")
    if neighbors:
        return trusted_html(nearest_neighbors_location(system["location"], neighbors, _system_url))
    name_to_id = {row["name"]: row["id"] for row in system.get("sector_siblings") or []}
    return trusted_html(linkify_location(system["location"], name_to_id, _system_url))


def _system_crumbs(system):
    """Home > Sectors > <sector> > <system>, or Home > Systems > <system>
    for a standalone system. The sector's name is the stored location's
    prefix (the system detail carries no separate sector name)."""
    if system["sector_id"] is None:
        return "systems", [crumb("Systems", "systems"), crumb(system["name"])]
    sector_name = (system["location"] or "").split(_LOCATION_MARKER, 1)[0].strip() or "Sector"
    return "sectors", [
        crumb("Sectors", "sectors"),
        crumb(sector_name, "sector", sector_id=system["sector_id"]),
        crumb(system["name"]),
    ]


def _inside_link(phenomenon_type, phenomenon_id, name):
    """The "Inside: <cloud>" badge's `{name, url}`, linking the nebula or
    supernova remnant something sits in (schema v39)."""
    return {"name": name, "url": page_url("phenomenon", phenomenon_type=phenomenon_type,
                                          phenomenon_id=phenomenon_id)}


def _phenomenon_inside(detail):
    """`_inside_link` for a phenomenon's own `inside_nebula_id`/
    `inside_remnant_id`, or `None` in open space. The cloud's name costs
    one more lookup."""
    for column, kind in (("inside_nebula_id", "nebula"), ("inside_remnant_id", "supernova_remnant")):
        cloud_id = detail.get(column)
        if cloud_id is not None:
            cloud = apiclient.get_phenomenon(db_name(), kind, cloud_id)
            return _inside_link(kind, cloud_id, cloud["name"])
    return None


def _badges(system):
    badges = []
    if system["quadrant"]:
        badges.append(f"Octant {system['quadrant']}")
    if system["is_binary"]:
        suffix = {"close": " (close)", "wide": " (wide)"}.get(system.get("binary_configuration"), "")
        badges.append(f"Binary system{suffix}")
    else:
        badges.append("Single star")
    runaway = runaway_text(system)
    if runaway:
        badges.append(runaway)
    return badges


def _code_buttons(system_id, code_fmt):
    """The Wikitext/Markdown toggle links: `?code=<fmt>` shows that code
    box, and the active one becomes "Hide <label>" (back to no code)."""
    buttons = []
    for fmt, label in CODE_FORMATS.items():
        active = fmt == code_fmt
        buttons.append({
            "label": f"Hide {label}" if active else label,
            "url": page_url("system", system_id=system_id, code=None if active else fmt,
                            _anchor="system-panel"),
            "active": active,
        })
    return buttons


def _wiki_upload_options(system, wiki_config):
    """The backends an admin may still upload to: configured site-wide
    and not yet uploaded for this system."""
    uploaded = {"wikijs": system["wikijs_url"], "mediawiki": system["mediawiki_url"]}
    return [{"value": value, "label": label} for value, label in WIKI_BACKENDS.items()
            if wiki_config.get(value) and not uploaded[value]]


def _species_by_planet(db, system):
    """`{planet_id: {"name", "url"}}`: each life world's dominant species
    (only planets with life can have one), or `{}` before population has
    made any species."""
    if not population_status()["species"]:
        return {}
    found = {}
    for planet in system["planets"]:
        if not planet.get("life_chemical"):
            continue
        try:
            species = apiclient.get_planet_species(db, planet["id"])
        except (apiclient.ApiError, apiclient.NotFoundError):
            species = None
        if species:
            found[planet["id"]] = {"name": species["name"],
                                   "url": page_url("species_page", species_id=species["id"])}
    return found


def _territory(db, system_id):
    """`{"name", "url", "color"}` of the polity whose territory this
    system is in, or `None` (also before population has made any)."""
    if not population_status()["territories"]:
        return None
    try:
        owner = apiclient.get_system_owner(db, system_id)
    except (apiclient.ApiError, apiclient.NotFoundError):
        return None
    if not owner:
        return None
    return {"name": owner["polity_name"], "color": owner.get("color"),
            "url": page_url("polity_page", polity_id=owner["polity_id"])}


@bp.route("/system/<int:system_id>", methods=["GET", "POST"])
def system(system_id):
    """One star system. A POST with a `facility_action` is the admin
    facility form (`web/system_facilities.py`); any other POST is the
    admin "Upload to Wiki" form."""
    facility_action = request.form.get("facility_action") if request.method == "POST" else None
    if request.method == "POST" and facility_action not in system_facilities.ACTIONS:
        return _upload_to_wiki(system_id)

    db = db_name()
    code_fmt = request.args.get("code")
    if code_fmt not in CODE_FORMATS:
        code_fmt = None

    detail = apiclient.get_system(db, system_id)
    facilities = apiclient.get_system_facilities(db, system_id)
    facility_form = None
    if facility_action:
        facility_form = system_facilities.handle_post(system_id, detail, facilities)
        if not isinstance(facility_form, dict):
            return facility_form
    sections = apiclient.get_system_sections(db, system_id)
    code_content = apiclient.get_system_text(db, system_id, code_fmt)["content"] if code_fmt else None

    admin = current_admin()
    wiki_options = _wiki_upload_options(detail, apiclient.get_wiki_config()) if admin else []
    wiki_status = WIKI_MESSAGES.get(request.args.get("wiki")) if admin else None

    section, breadcrumbs = _system_crumbs(detail)
    inside = detail.get("inside")
    map_html = ""
    if detail["stars"]:
        map_html = render_system_map_panel(detail, detail["stars"], detail["planets"], detail["belts"],
                                           facilities=facilities)
    # NAV measures from a sector position, so a standalone system gets no
    # links; the NAV page itself works out whether cross-sector NAV applies.
    links = nav_links("system", system_id) if detail["sector_id"] is not None else None

    return render_page(
        "system.html",
        title=detail["name"],
        section=section,
        breadcrumbs=breadcrumbs,
        description=f"The {detail['name']} star system: its stars, planets, moons, belts and comets.",
        system=detail,
        badges=_badges(detail),
        inside=_inside_link(inside["type"], inside["id"], inside["name"]) if inside else None,
        nav_links=links,
        location_html=_location_html(detail),
        map_html=trusted_html(map_html),
        code_fmt=code_fmt,
        code_label=CODE_FORMATS.get(code_fmt),
        code_content=code_content,
        code_rows=min(code_content.count("\n") + 3, 30) if code_content else 0,
        code_buttons=_code_buttons(system_id, code_fmt),
        system_list_html=trusted_html(system_list_html(detail, sections, class_url, facilities,
                                                       species=_species_by_planet(db, detail))),
        territory=_territory(db, system_id),
        stars_html=trusted_html(stars_html(detail["stars"], class_url)),
        wiki_options=wiki_options,
        wiki_status=wiki_status,
        admin=admin,
        facilities=system_facilities.facility_rows(facilities),
        facility_status=system_facilities.MESSAGES.get(request.args.get("facility")) if admin else None,
        facility_form=facility_form or {"fields": {}, "errors": [], "preview": None},
        facility_hosts=system_facilities.host_options(detail) if admin else [],
        facility_placements=system_facilities.PLACEMENT_OPTIONS,
        facility_kinds=system_facilities.kind_options(),
        facility_orbit_steps=system_facilities.FACILITY_ORBIT_STEPS,
        facility_orbit_default=system_facilities.ORBIT_STEP_DEFAULT,
    )


def _upload_to_wiki(system_id):
    """
    The admin "Upload to Wiki" POST (the CSRF token was already checked
    app-wide). Always answers with a 303 back to the GET page carrying a
    fixed `?wiki=<outcome>` code, so a reload never re-posts.
    """
    if current_admin() is None:
        activitylog.event("AUTHZ", "admin.required", path=request.path)
        abort(403)
    backend = request.form.get("backend", "")
    path = (request.form.get("path") or "").strip() or None
    outcome = "uploaded"
    if backend not in WIKI_BACKENDS:
        outcome = "invalid"
    else:
        try:
            apiclient.upload_system_to_wiki(request.headers.get("Cookie"), db_name(), system_id, backend, path)
        except apiclient.NotFoundError:
            raise
        except apiclient.ApiError as exc:
            outcome = {
                409: "exists", 400: "invalid", 401: "forbidden", 403: "forbidden", 501: "unconfigured",
            }.get(exc.status_code, "failed")
    return redirect(page_url("system", system_id=system_id, wiki=outcome, _anchor="wiki-upload"), code=303)


# ---------------------------------------------------------------------
# /phenomena and /phenomenon/<type>/<id>
# ---------------------------------------------------------------------

TYPE_LABELS = {
    "nebula": "Nebula", "asteroid_field": "Asteroid Field",
    "black_hole": "Black Hole", "neutron_star": "Neutron Star",
    "supernova_remnant": "Supernova Remnant",
    "rogue_planet": "Rogue Planet", "interstellar_comet": "Interstellar Comet",
    "quasar": "Quasar",
}
"""dict: Every phenomenon type (`queryDb._PHENOMENON_TABLES`) and its label."""


def _title_case(value):
    return (value or "").replace("_", " ").replace("-", " ").capitalize()


@bp.route("/phenomena")
def phenomena():
    """Every exotic phenomenon, placed on the galaxy map or not, paged
    with `?page=N`."""
    envelope, page = fetch_page(
        lambda limit, offset: apiclient.get_phenomena(db_name(), limit=limit, offset=offset),
        parse_page(request.args.get("page")),
    )
    rows = [{
        "name": row["name"],
        "url": page_url("phenomenon", phenomenon_type=row["type"], phenomenon_id=row["id"]),
        "type": TYPE_LABELS.get(row["type"], row["type"]),
        "descriptor": _title_case(row["descriptor"]),
        "radius": format_distance_ly(row["radius_ly"]) if row["radius_ly"] else None,
        "sector_name": row["sector_name"],
        "sector_url": page_url("sector", sector_id=row["sector_id"]) if row["sector_id"] is not None else None,
        "placed": bool(row["placed"]),
    } for row in envelope["items"]]
    return render_page(
        "phenomena.html",
        title="Phenomena",
        section="phenomena",
        breadcrumbs=[crumb("Phenomena")],
        description="Every nebula, black hole, neutron star and other exotic phenomenon in this generated galaxy.",
        rows=rows,
        total=envelope["total"],
        pager=pager("page", page, envelope["total"], anchor="phenomena-list", label="Phenomena pages"),
    )


# Every distance and non-body radius goes through the distance ladder
# (fmt.format_distance_*); a rogue planet's radius is a body radius, so it
# is km in scientific notation like every planet's.
_ly = format_distance_ly


def _km(value):
    return format_distance_km(value)


def _body_radius(value):
    return to_plain_text(format_star_radius(value))


_DAYS_PER_YEAR = 365.25


def _light_days(value):
    return format_distance_ly(value / _DAYS_PER_YEAR)


def _bool_text(value):
    return "Yes" if value else "No"


def _progenitor_text(value):
    # "Type Ia" is stored correctly cased; _title_case would lowercase "Ia".
    return value if value == "Type Ia" else _title_case(value)


def _rogue_planet_type_text(value):
    return {"t": "Terrestrial", "g": "Gas Giant"}.get(value, value)


def _cloud_class_text(value):
    """A nebula or remnant class letter with its name, e.g. "D: Classical
    H II region" (`program_constants.NEBULA_CLASSES`, schema v38)."""
    entry = NEBULA_CLASSES.get(value)
    return f"{value}: {entry['name']}" if entry else value


# The contents a nebula or supernova remnant carries (schema v38). An
# unset value is stored as 0 or ''; those rows are left out.
_CLOUD_CONTENTS = [
    ("dominant_species", "Dominant Species", lambda v: v or None),
    ("density_cm3", "Density", lambda v: f"{format_number(v, ',.3g')} particles/cm\u00b3" if v else None),
    ("temperature_k", "Gas Temperature", lambda v: f"{format_number(v, ',.0f')} K" if v else None),
    ("extinction_av", "Extinction", lambda v: f"{format_number(v, ',.2f')} magnitudes (visual)" if v else None),
]


def _rogue_mass_bin_text(value):
    return ROGUE_MASS_CLASS_NAMES.get(value, value)


def _optional(fmt):
    """A formatter that leaves a `None` (not applicable) row out."""
    return lambda v: fmt(v) if v is not None else None


_SPEED = ("galactic_orbital_speed_kms", "Galactic Orbital Speed", lambda v: f"{format_number(v, ',.1f')} km/s")
_PERIOD = ("galactic_orbital_period_gy", "Galactic Orbital Period", lambda v: f"{format_number(v, ',.2f')} Gy")

FIELD_SPECS = {
    "nebula": [
        ("nebula_class", "Class", _cloud_class_text),
        ("nebula_type", "Nebula Type", _title_case),
        ("radius_ly", "Radius", _ly),
        ("composition", "Composition", str),
        ("formation_cause", "Formation", str),
        *_CLOUD_CONTENTS,
        _SPEED, _PERIOD,
    ],
    "asteroid_field": [
        ("field_class", "Class", str),
        ("composition_family", "Composition Family", _title_case),
        ("density", "Density", _title_case),
        ("radius_ly", "Radius", _ly),
        ("composition_summary", "Composition", str),
        _SPEED, _PERIOD,
    ],
    "black_hole": [
        ("mass_class", "Class", _title_case),
        ("mass_solar", "Mass", lambda v: f"{format_number(v, ',.2f')} solar masses"),
        ("event_horizon_radius_km", "Event Horizon Radius", _km),
        ("spin", "Spin (dimensionless)", lambda v: f"{v:.3f}"),
        ("has_accretion_disk", "Accretion Disk", _bool_text),
        ("temperature_k", "Hawking Temperature", lambda v: f"{v:.2e} K"),
        ("luminosity_w", "Luminosity", lambda v: f"{format_number(v, '.2e')} W"),
        ("age_gy", "Age", lambda v: f"{format_number(v, ',.2f')} Gy"),
        (_SPEED[0], _SPEED[1], _optional(_SPEED[2])),
        (_PERIOD[0], _PERIOD[1], _optional(_PERIOD[2])),
    ],
    "neutron_star": [
        ("mass_solar", "Mass", lambda v: f"{format_number(v, ',.2f')} solar masses"),
        ("radius_km", "Radius", _km),
        ("spin_period_ms", "Spin Period", lambda v: f"{format_number(v, ',.2f')} ms"),
        ("magnetic_field_gauss", "Magnetic Field", lambda v: f"{v:.2e} G"),
        ("pulsar_type", "Pulsar Type", _title_case),
        ("surface_temperature_k", "Surface Temperature", lambda v: f"{format_number(v, ',.0f')} K"),
        ("luminosity_w", "Luminosity", lambda v: f"{format_number(v, '.2e')} W"),
        ("age_gy", "Age", lambda v: f"{format_number(v, ',.2f')} Gy"),
        (_SPEED[0], _SPEED[1], _optional(_SPEED[2])),
        (_PERIOD[0], _PERIOD[1], _optional(_PERIOD[2])),
    ],
    "supernova_remnant": [
        ("remnant_class", "Class", _cloud_class_text),
        ("morphology", "Morphology", _title_case),
        ("progenitor_type", "Progenitor Type", _progenitor_text),
        ("age_years", "Age", lambda v: f"{format_number(v, ',.0f')} years"),
        ("radius_ly", "Radius", _ly),
        ("compact_remnant_kind", "Compact Remnant Left Behind", _title_case),
        *_CLOUD_CONTENTS,
        _SPEED, _PERIOD,
    ],
    "rogue_planet": [
        ("planet_type", "Type", _rogue_planet_type_text),
        ("mass_bin", "Mass Class", _rogue_mass_bin_text),
        ("mass_kg", "Mass", lambda v: f"{v:.2e} kg"),
        ("radius_km", "Radius", _body_radius),
        ("composition", "Composition", str),
        ("has_internal_heat", "Internal Heat", _bool_text),
        ("has_moons", "Has Moons", _bool_text),
        _SPEED, _PERIOD,
    ],
    "interstellar_comet": [
        ("nucleus_diameter_km", "Nucleus Diameter", _km),
        ("velocity_kms", "Velocity", lambda v: f"{format_number(v, ',.1f')} km/s"),
        ("is_active", "Active", _bool_text),
        ("composition_summary", "Composition", str),
        _SPEED, _PERIOD,
    ],
    # No galactic-orbit rows: a quasar sits at the galactic center.
    "quasar": [
        ("black_hole_mass_solar", "Black Hole Mass", lambda v: f"{format_number(v, '.2e')} solar masses"),
        ("event_horizon_radius_km", "Event Horizon Radius", _km),
        ("luminosity_w", "Luminosity", lambda v: f"{format_number(v, '.2e')} W"),
        ("eddington_ratio", "Eddington Ratio", lambda v: f"{v:.0%}"),
        ("accretion_rate_solar_per_year", "Accretion Rate", lambda v: f"{format_number(v, ',.2f')} solar masses/year"),
        ("broad_line_region_light_days", "Broad-Line Region Radius", _light_days),
        ("is_radio_loud", "Radio-Loud (Jets)", _bool_text),
        ("jet_length_ly", "Jet Length", _ly),
        ("active_age_years", "Active For", lambda v: f"{format_number(v, ',.0f')} years"),
    ],
}
"""dict: (column, label, formatter) per type. A formatter returns plain
text (the template escapes it) or `None` to leave the row out."""

CLASS_COLUMNS = {
    "nebula_class": "nebula",
    "remnant_class": "supernova-remnant",
    "field_class": "asteroid-field",
    "mass_class": "black-hole",
    "mass_bin": "rogue-planet",
}
"""dict: The columns whose value is a class, and the class type
(`lib/classref.py`) whose page it links to. An asteroid field's "C3"
links by its letter."""


def phenomenon_fields(phenomenon_type, detail):
    """`[(label, text, url)]` for the data table; `url` is the value's
    class page, or `None`. Needs a request context."""
    fields = []
    for column, label, formatter in FIELD_SPECS.get(phenomenon_type, []):
        raw = detail.get(column)
        if raw is None:
            continue
        text = formatter(raw)
        if text is not None:
            url = class_url(CLASS_COLUMNS[column], raw) if column in CLASS_COLUMNS else None
            fields.append((label, text, url))
    return fields


@bp.route("/phenomenon/<phenomenon_type>/<int:phenomenon_id>")
def phenomenon(phenomenon_type, phenomenon_id):
    """One phenomenon: its view (a render, the AU-scale diagram, or none
    for an asteroid field; see `phenomenonrender.view_kind`) and its data
    table."""
    if phenomenon_type not in TYPE_LABELS:
        abort(404)
    detail = apiclient.get_phenomenon(db_name(), phenomenon_type, phenomenon_id)
    type_label = TYPE_LABELS[phenomenon_type]

    distance = None
    radius_pc = detail.get("galactic_radius_pc")
    if radius_pc is not None:
        distance = format_distance_pc(radius_pc)
    sector = None
    if detail.get("sector_id") is not None:
        sector = {"name": detail.get("sector_name") or "Sector",
                  "url": page_url("sector", sector_id=detail["sector_id"])}

    kind = view_kind(phenomenon_type)
    if kind == "render":
        map_html = render_phenomenon_view_panel(phenomenon_type, detail)
    elif kind == "map":
        map_html = render_phenomenon_map_panel(
            phenomenon_type, detail["name"], detail.get("radius_ly") or 0,
        )
    else:
        map_html = ""
    return render_page(
        "phenomenon.html",
        title=detail["name"],
        section="phenomena",
        breadcrumbs=[crumb("Phenomena", "phenomena"), crumb(detail["name"])],
        description=f"{detail['name']}, a {type_label.lower()} in this generated galaxy.",
        type_label=type_label,
        distance=distance,
        sector=sector,
        octant=detail.get("quadrant"),
        nearest_html=trusted_html(nearest_systems_html(detail["nearest"], _system_url))
        if detail.get("nearest") else None,
        inside=_phenomenon_inside(detail),
        nav_links=nav_links(phenomenon_type, detail["id"]),
        view_kind=kind,
        map_html=trusted_html(map_html),
        fields=phenomenon_fields(phenomenon_type, detail),
    )
