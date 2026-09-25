# html/web/system_pages.py

"""
The system and phenomenon pages (moved from the CGI `system.py`,
`phenomena.py` and `phenomenon.py`):

- `/system/<id>`: one star system -- the system map, the expandable body
  list with each body's generated description, the Wikitext/Markdown
  code views (`?code=wikitext|markdown`), the Stars/Planets/Belts/Comets
  tables, and for an admin an "Upload to Wiki" form (POST, CSRF-checked,
  then redirected back to the GET page).
- `/phenomena`: every exotic phenomenon, paged with `?page=N`.
- `/phenomenon/<type>/<id>`: one phenomenon's data and its AU-scale
  diagram.

Nothing here reads a database name from the request (`db_name()`), and
every link is a plain GET link (`page_url`).
"""

from flask import abort, redirect, request
from markupsafe import escape

import apiclient
from fmt import nearest_neighbors_location, linkify_location
from pagination import fetch_page, parse_page
from phenomenonmap import render_phenomenon_map_panel
from systemmap import render_system_map_panel
from systempage import bodies_html, stars_html, system_list_html

from . import bp
from .helpers import crumb, current_admin, db_name, page_url, pager, render_page, trusted_html

try:
    from stellarObjects.utils import pc_to_ly
except ImportError:  # pragma: no cover -- the package is always importable next to the API
    pc_to_ly = None

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
    "forbidden": (True, "Not allowed. Change the default admin username and password first."),
    "failed": (True, "The upload failed. The wiki could not be reached or returned an error."),
}
"""dict: The upload outcome shown after the POST-redirect-GET. Only these
fixed messages are shown, so nothing from the query string reaches the
page as text."""

_LOCATION_MARKER = " -- nearest: "


def _system_link(system_id, label_html):
    """`fmt`'s location `link` hook: a plain GET link to a system page.
    `label_html` is already escaped by `fmt`."""
    return f'<a href="{escape(page_url("system", system_id=system_id))}">{label_html}</a>'


def _location_html(system):
    """The "Location:" line: the sector name and nearest neighbours, each
    linked, or `None` for a system with no stored location."""
    if not system["location"]:
        return None
    neighbors = system.get("nearest_neighbors")
    if neighbors:
        return trusted_html(nearest_neighbors_location(None, system["location"], neighbors, link=_system_link))
    name_to_id = {row["name"]: row["id"] for row in system.get("sector_siblings") or []}
    return trusted_html(linkify_location(None, system["location"], name_to_id, link=_system_link))


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


def _badges(system):
    badges = []
    if system["quadrant"]:
        badges.append(f"Octant {system['quadrant']}")
    if system["is_binary"]:
        suffix = {"close": " (close)", "wide": " (wide)"}.get(system.get("binary_configuration"), "")
        badges.append(f"Binary system{suffix}")
    else:
        badges.append("Single star")
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


@bp.route("/system/<int:system_id>", methods=["GET", "POST"])
def system(system_id):
    """One star system. POST is the admin "Upload to Wiki" form."""
    if request.method == "POST":
        return _upload_to_wiki(system_id)

    db = db_name()
    code_fmt = request.args.get("code")
    if code_fmt not in CODE_FORMATS:
        code_fmt = None

    detail = apiclient.get_system(db, system_id)
    sections = apiclient.get_system_sections(db, system_id)
    code_content = apiclient.get_system_text(db, system_id, code_fmt)["content"] if code_fmt else None

    admin = current_admin()
    wiki_options = _wiki_upload_options(detail, apiclient.get_wiki_config()) if admin else []
    wiki_status = WIKI_MESSAGES.get(request.args.get("wiki")) if admin else None

    section, breadcrumbs = _system_crumbs(detail)
    map_html = ""
    if detail["stars"]:
        map_html = render_system_map_panel(detail, detail["stars"], detail["planets"], detail["belts"])
    nav_links = None
    if detail["sector_id"] is not None:
        # NAV measures from a sector position, so a standalone system gets
        # no links; nav.py itself works out whether cross-sector NAV applies.
        nav_links = {
            "from": page_url("nav", from_id=system_id),
            "to": page_url("nav", to_id=system_id),
        }

    return render_page(
        "system.html",
        title=detail["name"],
        section=section,
        breadcrumbs=breadcrumbs,
        description=f"The {detail['name']} star system: its stars, planets, moons, belts and comets.",
        system=detail,
        badges=_badges(detail),
        nav_links=nav_links,
        location_html=_location_html(detail),
        map_html=trusted_html(map_html),
        code_fmt=code_fmt,
        code_label=CODE_FORMATS.get(code_fmt),
        code_content=code_content,
        code_rows=min(code_content.count("\n") + 3, 30) if code_content else 0,
        code_buttons=_code_buttons(system_id, code_fmt),
        system_list_html=trusted_html(system_list_html(detail, sections)),
        stars_html=trusted_html(stars_html(detail["stars"])),
        bodies_html=trusted_html(bodies_html(
            detail["planets"], detail["belts"], detail["comets"], detail["stars"],
            detail.get("binary_configuration"),
        )),
        wiki_options=wiki_options,
        wiki_status=wiki_status,
    )


def _upload_to_wiki(system_id):
    """
    The admin "Upload to Wiki" POST (the CSRF token was already checked
    app-wide). Always answers with a 303 back to the GET page carrying a
    fixed `?wiki=<outcome>` code, so a reload never re-posts.
    """
    if current_admin() is None:
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
        "radius": f"{row['radius_ly']:,.2f} ly" if row["radius_ly"] else None,
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


def _ly(value):
    return f"{value:,.2f} ly"


def _bool_text(value):
    return "Yes" if value else "No"


def _progenitor_text(value):
    # "Type Ia" is stored correctly cased; _title_case would lowercase "Ia".
    return value if value == "Type Ia" else _title_case(value)


def _rogue_planet_type_text(value):
    return {"t": "Terrestrial", "g": "Gas Giant"}.get(value, value)


def _optional(fmt):
    """A formatter that leaves a `None` (not applicable) row out."""
    return lambda v: fmt(v) if v is not None else None


_SPEED = ("galactic_orbital_speed_kms", "Galactic Orbital Speed", lambda v: f"{v:,.1f} km/s")
_PERIOD = ("galactic_orbital_period_gy", "Galactic Orbital Period", lambda v: f"{v:,.2f} Gy")

FIELD_SPECS = {
    "nebula": [
        ("nebula_type", "Nebula Type", _title_case),
        ("radius_ly", "Radius", _ly),
        ("composition", "Composition", str),
        ("formation_cause", "Formation", str),
        _SPEED, _PERIOD,
    ],
    "asteroid_field": [
        ("density", "Density", _title_case),
        ("radius_ly", "Radius", _ly),
        ("composition_summary", "Composition", str),
        _SPEED, _PERIOD,
    ],
    "black_hole": [
        ("mass_solar", "Mass", lambda v: f"{v:,.2f} solar masses"),
        ("event_horizon_radius_km", "Event Horizon Radius", lambda v: f"{v:,.1f} km"),
        ("spin", "Spin (dimensionless)", lambda v: f"{v:.3f}"),
        ("has_accretion_disk", "Accretion Disk", _bool_text),
        ("temperature_k", "Hawking Temperature", lambda v: f"{v:.2e} K"),
        ("luminosity_w", "Luminosity", lambda v: f"{v:.2e} W"),
        ("age_gy", "Age", lambda v: f"{v:,.2f} Gy"),
        (_SPEED[0], _SPEED[1], _optional(_SPEED[2])),
        (_PERIOD[0], _PERIOD[1], _optional(_PERIOD[2])),
    ],
    "neutron_star": [
        ("mass_solar", "Mass", lambda v: f"{v:,.2f} solar masses"),
        ("radius_km", "Radius", lambda v: f"{v:,.2f} km"),
        ("spin_period_ms", "Spin Period", lambda v: f"{v:,.2f} ms"),
        ("magnetic_field_gauss", "Magnetic Field", lambda v: f"{v:.2e} G"),
        ("pulsar_type", "Pulsar Type", _title_case),
        ("surface_temperature_k", "Surface Temperature", lambda v: f"{v:,.0f} K"),
        ("luminosity_w", "Luminosity", lambda v: f"{v:.2e} W"),
        ("age_gy", "Age", lambda v: f"{v:,.2f} Gy"),
        (_SPEED[0], _SPEED[1], _optional(_SPEED[2])),
        (_PERIOD[0], _PERIOD[1], _optional(_PERIOD[2])),
    ],
    "supernova_remnant": [
        ("morphology", "Morphology", _title_case),
        ("progenitor_type", "Progenitor Type", _progenitor_text),
        ("age_years", "Age", lambda v: f"{v:,.0f} years"),
        ("radius_ly", "Radius", _ly),
        ("compact_remnant_kind", "Compact Remnant Left Behind", _title_case),
        _SPEED, _PERIOD,
    ],
    "rogue_planet": [
        ("planet_type", "Type", _rogue_planet_type_text),
        ("mass_kg", "Mass", lambda v: f"{v:.2e} kg"),
        ("radius_km", "Radius", lambda v: f"{v:,.0f} km"),
        ("composition", "Composition", str),
        ("has_internal_heat", "Internal Heat", _bool_text),
        ("has_moons", "Has Moons", _bool_text),
        _SPEED, _PERIOD,
    ],
    "interstellar_comet": [
        ("nucleus_diameter_km", "Nucleus Diameter", lambda v: f"{v:,.2f} km"),
        ("velocity_kms", "Velocity", lambda v: f"{v:,.1f} km/s"),
        ("is_active", "Active", _bool_text),
        ("composition_summary", "Composition", str),
        _SPEED, _PERIOD,
    ],
}
"""dict: (column, label, formatter) per type. A formatter returns plain
text (the template escapes it) or `None` to leave the row out."""


def phenomenon_fields(phenomenon_type, detail):
    """`[(label, text)]` for the data table."""
    fields = []
    for column, label, formatter in FIELD_SPECS.get(phenomenon_type, []):
        raw = detail.get(column)
        if raw is None:
            continue
        text = formatter(raw)
        if text is not None:
            fields.append((label, text))
    return fields


@bp.route("/phenomenon/<phenomenon_type>/<int:phenomenon_id>")
def phenomenon(phenomenon_type, phenomenon_id):
    """One phenomenon: its data table and AU-scale diagram."""
    if phenomenon_type not in TYPE_LABELS:
        abort(404)
    detail = apiclient.get_phenomenon(db_name(), phenomenon_type, phenomenon_id)
    type_label = TYPE_LABELS[phenomenon_type]

    distance = None
    radius_pc = detail.get("galactic_radius_pc")
    if radius_pc is not None:
        distance = _ly(pc_to_ly(radius_pc) if pc_to_ly else radius_pc * 3.2616)
    sector = None
    if detail.get("sector_id") is not None:
        sector = {"name": detail.get("sector_name") or "Sector",
                  "url": page_url("sector", sector_id=detail["sector_id"])}

    # Offered for every type: nav.py itself says when a phenomenon was
    # never placed in the galaxy.
    nav_links = {
        "from": page_url("nav", from_id=detail["id"], from_kind="phenomenon", from_type=phenomenon_type),
        "to": page_url("nav", to_id=detail["id"], to_kind="phenomenon", to_type=phenomenon_type),
    }
    map_html = render_phenomenon_map_panel(
        phenomenon_type, detail["name"], detail.get("radius_ly") or 0, include_scripts=False,
    )
    return render_page(
        "phenomenon.html",
        title=detail["name"],
        section="phenomena",
        breadcrumbs=[crumb("Phenomena", "phenomena"), crumb(detail["name"])],
        description=f"{detail['name']}, a {type_label.lower()} in this generated galaxy.",
        type_label=type_label,
        distance=distance,
        sector=sector,
        nav_links=nav_links,
        map_html=trusted_html(map_html),
        fields=phenomenon_fields(phenomenon_type, detail),
    )
