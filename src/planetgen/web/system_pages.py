# planetgen/web/system_pages.py

"""
The system and phenomenon pages (moved from the CGI `system.py`,
`phenomena.py` and `phenomenon.py`):

- `/system/<id>`: one star system -- the system map, the expandable body
  list with each body's generated description, the Markdown code view
  (`?code=markdown`), the Stars/Planets/Belts/Comets tables, the
  system's facilities, and for an admin the facility form and Remove
  buttons (`web/system_facilities.py`).
- `/phenomena`: every exotic phenomenon, paged with `?page=N`.
- `/phenomenon/<type>/<id>`: one phenomenon's data and its AU-scale
  diagram.

Nothing here reads a database name from the request (`db_name()`), and
every link is a plain GET link (`page_url`).
"""

import math

from flask import abort, jsonify, redirect, request
from markupsafe import Markup
from planetgen.describe.html import blocks_html
from planetgen.describe.model import sections, split_by_body
from planetgen.describe.stored_phenomena import stored_phenomenon_document
from planetgen.physics.quantity import Quantity
from planetgen.physics.lensing import einstein_radius_against_galactic_centre_m
from planetgen.web.lib import apiclient, facility_display
from planetgen.web.lib.slowroute import slow_route
from planetgen.web.lib.fmt import (
    distance_km_html, distance_ly_html, format_members, format_number, linkify_nearest,
    measured_html, nearest_systems_html, period_years_html, runaway_text, speed_kms_html,
)
from planetgen.web.lib.datatable import Column, Facet, Result, Table, read_from_start
from planetgen.web.maps.starmap import HYPERVELOCITY_EXPECTED_STAR
from planetgen.web import system_orbits_view
from planetgen.web.scatter_members_view import members_view
from planetgen.web.maps.phenomenonmap import render_phenomenon_map_panel
from planetgen.web.maps.phenomenonrender import render_nebula_view_panel, render_phenomenon_view_panel, view_kind
from planetgen.web.maps.systemmap import render_system_map_panel
from planetgen.web.lib.tabledisplay import format_star_radius, to_plain_text

from planetgen.web.lib.classref import ROGUE_MASS_CLASS_NAMES
from planetgen.generation.belt import format_composition_summary
from planetgen.generation.phenomena.compact_remnant import hawking_luminosity_w, hawking_temperature_k
from planetgen.generation.phenomena.rogue import format_comet_composition_summary
from planetgen.tuning import NEBULA_CLASSES
from planetgen.physics.rogue_surface import SURFACE_REGIME_LABELS
from planetgen.util import chem

from . import bp, habitability_symbol_row, tables
from . import edit_actions, glance, item_pages, system_facilities, system_species
from .class_pages import class_url
from .helpers import (
    bookmark, crumb, current_admin, db_name, page_url, population_status, render_page, trusted_html,
)
from .nav_page import endpoint, nav_url
from .nearby_page import DEFAULT_DISTANCE_PC, nearby_url

# ---------------------------------------------------------------------
# /system/<id>
# ---------------------------------------------------------------------

CODE_FORMATS = {"markdown": "Markdown"}
"""dict: `?code=` values the system page accepts, with their labels."""

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


def pick_nav(kind, entity_id, args):
    """
    `(nav_links, pick)` for a system or phenomenon page (NAV.15): the
    Navigate links, or while a NAV start or destination is picked
    (`?pick=from|to`, `sector_page._pick_mode`) `{"pick", "pickLabel"}`,
    the one button that returns to the NAV page with this object, and the
    pick mode for the banner and the bookmarks that keep it.
    """
    from .sector_page import _pick_mode  # a page module like this one; imported where used

    links = nav_links(kind, entity_id)
    pick = _pick_mode(args)
    if pick is None:
        return links, None
    here = endpoint(kind, entity_id)
    back = (nav_url(origin=here, destination=pick["other"]) if pick["pick"] == "from"
            else nav_url(origin=pick["other"], destination=here))
    return {"pick": back, "pickLabel": pick["label"]}, pick


def body_nav(args, course=True):
    """
    The links a body's info panel on the System Map gets, as URL templates
    with `{ref}` where the body's object reference goes (`static/systemnav.js`):
    the NAV buttons (NAV.50), `take` and `label` while a NAV start or
    destination is being picked (`?pick=from|to`), else `start` and `end`,
    the "What's within N pc" link (`nearby`, `nearby-label`, NAV.45) and the
    View button to the body's own page (`view`, UX.143, MAP.254), which a course
    being picked leaves out. `course=False` (a system with no
    sector, which NAV cannot measure from) gives the nearby link alone.
    """
    from .sector_page import _pick_mode

    token = "__REF__"
    pick = _pick_mode(args) if course else None
    nearby = {"nearby": nearby_url(token, DEFAULT_DISTANCE_PC).replace(token, "{ref}"),
              "nearby-label": f"What’s within {DEFAULT_DISTANCE_PC} pc",
              "view": page_url("object_page", ref=token).replace(token, "{ref}")}
    if pick is None:
        nav = ({"start": nav_url(origin=token).replace(token, "{ref}"),
                "end": nav_url(destination=token).replace(token, "{ref}")} if course else {})
        return {**nav, **nearby}
    url = (nav_url(origin=token, destination=pick["other"]) if pick["pick"] == "from"
           else nav_url(origin=pick["other"], destination=token))
    return {"take": url.replace(token, "{ref}"), "label": pick["label"]}


def _nearest_html(system):
    """The "Nearest:" line: the three nearest neighbours, each linked with its
    distance (the breadcrumb carries the sector, UX.52), or `None` for a
    system with none."""
    neighbors = system.get("nearest_neighbors")
    if neighbors:
        return trusted_html(nearest_systems_html(neighbors, _system_url))
    name_to_id = {row["name"]: row["id"] for row in system.get("sector_siblings") or []}
    text = linkify_nearest(system["location"], name_to_id, _system_url)
    return trusted_html(text) if text else None


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


def _phenomenon_inside_cloud(detail):
    """A phenomenon's own `inside_nebula_id`/`inside_remnant_id` as the Location block's `inside` (`type`, `id`, `name`),
    or `None` in open space."""
    for column, kind in (("inside_nebula_id", "nebula"), ("inside_remnant_id", "supernova_remnant")):
        cloud_id = detail.get(column)
        if cloud_id is not None:
            return {"type": kind, "id": cloud_id, "name": apiclient.get_phenomenon(db_name(), kind, cloud_id)["name"]}
    return None


def _phenomenon_contents(phenomenon_type, phenomenon_id):
    """The "What is inside it" block of a nebula or supernova remnant (UX.144), or `None` for any other type."""
    if phenomenon_type not in ("nebula", "supernova_remnant"):
        return None
    try:
        return item_pages.contents(apiclient.get_phenomenon_contents(db_name(), phenomenon_type, phenomenon_id))
    except (apiclient.ApiError, apiclient.NotFoundError):
        return None   # the list is extra: the page opens without it


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
    """The Markdown toggle link: `?code=<fmt>` shows that code
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


def _with_keep_out(context, ref):
    """A page's sheet context with the object's keep-out row first in its At a glance panel (UX.192). The page opens
    without the row when the object cannot be read, as it does without its contents list."""
    try:
        found = apiclient.get_object(db_name(), ref)
    except (apiclient.ApiError, apiclient.NotFoundError):
        return context
    row, region = glance.keep_out_row(found), item_pages.region_row(found)
    region = region and {"label": region["label"], "html": trusted_html(region["html"])}
    return {**context, "glance": ([row] if row else []) + context["glance"] + ([region] if region else [])}


@bp.route("/system/<uid:system_id>/scene")
def system_scene(system_id):
    """The 3D system view's scene JSON (`GET /api/systems/<id>/scene`,
    MAP.69), for `static/systempage3d.js`."""
    response = jsonify(apiclient.get_system_scene(db_name(), system_id))
    response.headers["Cache-Control"] = "no-store"
    return response


system_scene.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/system/<uid:system_id>", methods=["GET", "POST"])
def system(system_id):
    """One star system. A POST with a `facility_action` is the admin
    facility form (`web/system_facilities.py`); any other POST is refused."""
    facility_action = request.form.get("facility_action") if request.method == "POST" else None
    if request.method == "POST" and "edit_action" in request.form:
        return _edit_post(system_id)
    if request.method == "POST" and facility_action not in system_facilities.ACTIONS:
        abort(400)

    db = db_name()
    code_fmt = request.args.get("code")
    if code_fmt not in CODE_FORMATS:
        code_fmt = None

    detail = apiclient.get_system(db, system_id)
    facilities = apiclient.get_system_facilities(db, system_id)
    facility_form = None
    facility_context = {"types": [], "owner": None}
    if facility_action or current_admin() is not None:
        facility_context = system_facilities.form_context(system_id)
    if facility_action:
        facility_form = system_facilities.handle_post(system_id, detail, facilities, facility_context)
        if not isinstance(facility_form, dict):
            return facility_form
    sections = split_by_body(apiclient.get_system_description(db, system_id))
    code_content = apiclient.get_system_text(db, system_id, code_fmt)["content"] if code_fmt else None

    admin = current_admin()

    section, breadcrumbs = _system_crumbs(detail)
    inside = detail.get("inside")
    map_html = ""
    if detail["stars"]:
        map_html = render_system_map_panel(detail, detail["stars"], detail["planets"], detail["belts"],
                                           facilities=facilities,
                                           scene_url=page_url("system_scene", system_id=system_id),
                                           nav=body_nav(request.args, course=detail["sector_id"] is not None))
    # NAV measures from a sector position, so a standalone system gets no
    # links; the NAV page itself works out whether cross-sector NAV applies.
    links, pick = pick_nav("system", system_id, request.args) if detail["sector_id"] is not None else (None, None)

    tree, tree_marks = item_pages.system_tree(db, system_id, detail["name"])
    return render_page(
        "system.html",
        title=detail["name"],
        section=section,
        breadcrumbs=breadcrumbs,
        description=f"The {detail['name']} star system: its stars, planets, moons, belts and comets.",
        system=detail,
        badges=_badges(detail),
        inside=_inside_link(inside["type"], inside["id"], inside["name"]) if inside else None,
        location=item_pages.sector_location(detail.get("sector_id"), detail.get("sector_name"), inside),
        nav_links=links,
        nearby_url=nearby_url(endpoint("system", system_id)),
        pick=pick,
        bookmark=bookmark("system", endpoint("system", system_id), detail["name"],
                          page_url("system", system_id=system_id)),
        nearest_html=_nearest_html(detail),
        map_html=trusted_html(map_html),
        code_fmt=code_fmt,
        code_label=CODE_FORMATS.get(code_fmt),
        code_content=code_content,
        code_rows=min(code_content.count("\n") + 3, 30) if code_content else 0,
        code_buttons=_code_buttons(system_id, code_fmt),
        overview=trusted_html(blocks_html(sections["overview"], item_pages._ref_link) if sections["overview"] else ""),
        glance=glance.system_rows(detail, facilities, habitability_symbol_row),
        tree=tree,
        tree_marks=tree_marks,
        orbits=system_orbits_view.orbit_rows(detail),
        territory=_territory(db, system_id),
        admin=admin,
        facilities=system_facilities.facility_rows(facilities),
        habitable_species=system_species.panel_rows(db, system_id),
        facility_status=system_facilities.MESSAGES.get(request.args.get("facility")) if admin else None,
        facility_form=facility_form or {"fields": {}, "errors": [], "preview": None},
        facility_hosts=system_facilities.host_options(detail) if admin else [],
        facility_context=facility_context if admin else {"types": [], "owner": None},
        facility_ratings=system_facilities.RATING_CHOICES,
        facility_dimensions=facility_display.DIMENSION_LABELS,
        facility_placements=system_facilities.PLACEMENT_OPTIONS,
        facility_kinds=system_facilities.kind_options(),
        facility_orbit_steps=system_facilities.FACILITY_ORBIT_STEPS,
        facility_orbit_default=system_facilities.ORBIT_STEP_DEFAULT,
        can_edit=edit_actions.can_edit(admin),
        star_editable=not detail.get("is_binary") and len(detail["stars"]) == 1,
        star_type=detail["stars"][0]["star_type"] if detail["stars"] else "",
    )


def _edit_post(system_id):
    """An admin Delete or Regenerate POST (`web/edit_actions.py`) for the system itself; a body's own are posted to its
    page (`web/item_pages.py`)."""
    detail = apiclient.get_system(db_name(), system_id)
    allowed = {("system", system_id)}
    gone = (page_url("sector", sector_id=detail["sector_id"]) if detail["sector_id"] is not None
            else page_url("systems"))
    return edit_actions.handle_post(allowed, page_url("system", system_id=system_id, _anchor="admin-menu"),
                                    after_delete={"system": gone})


# ---------------------------------------------------------------------
# /phenomena and /phenomenon/<type>/<id>
# ---------------------------------------------------------------------

TYPE_LABELS = {
    "nebula": "Nebula", "asteroid_field": "Asteroid Field",
    "black_hole": "Black Hole", "neutron_star": "Neutron Star",
    "supernova_remnant": "Supernova Remnant",
    "rogue_planet": "Rogue Planet", "interstellar_comet": "Interstellar Comet",
    "quasar": "Quasar", "hypervelocity_star": "Hypervelocity Star",
}
"""dict: Every phenomenon type (`planetgen.db.query._PHENOMENON_TABLES`) and its label."""


def _title_case(value):
    return (value or "").replace("_", " ").replace("-", " ").capitalize()


def _phenomena_load(state, limit, offset, want_facets, after=None):
    """The Phenomena table's page of rows (`lib/datatable.Result`); `next` is the last scattered row's id when the
    page ends inside the scattered tail, which every page after it is read from by key (PERF.75), and otherwise the
    `<type>:<id>` of the last built row, whatever the sort (PERF.205, PERF.211)."""
    filters = {"types": state.filters["type"], "descriptors": state.filters["descriptor"]}
    envelope, note = read_from_start(
        lambda start: apiclient.get_phenomena(
            db_name(), limit=limit, offset=offset, sort=state.sort, descending=state.descending,
            facets=want_facets, after=after, start=start, **filters),
        state.start)
    rows = [[
        {"text": row["name"], "href": page_url("uncharted_phenomenon", scatter_id=row["scatter_id"])}
        if row["scattered"] else
        {"text": row["name"], "href": page_url("phenomenon", phenomenon_type=row["type"], phenomenon_id=row["id"])},
        {"text": TYPE_LABELS.get(row["type"], row["type"])},
        {"text": _title_case(row["descriptor"])},
        {"text": distance_ly_html(row["radius_ly"])} if row["radius_ly"] else {"text": "\u2013"},
        {"text": row["sector_name"], "href": page_url("sector", sector_id=row["sector_id"])}
        if row["sector_id"] is not None else {"text": "None", "muted": True},
        {"text": "Yes" if row["placed"] else "No"},
        {"text": format_members(row.get("member_count"), row.get("system_mass_solar"))},
    ] for row in envelope["items"]]
    facets = None
    if want_facets:
        facets = {
            "type": [{"value": option["value"], "label": TYPE_LABELS.get(option["value"], option["value"]),
                      "count": option["count"]} for option in envelope["facets"]["type"]],
            "descriptor": [{"value": option["value"], "label": _title_case(option["value"]),
                            "count": option["count"]} for option in envelope["facets"]["descriptor"]],
        }
    last = envelope["items"][-1] if envelope["items"] else None
    following = None
    if last and last["scattered"]:
        following = last["scatter_id"]
    elif last:      # the built rows page by key too, in the order of any sort (PERF.205, PERF.211)
        following = f"{last['type']}:{last['id']}"
    return Result(rows, envelope["total"], facets, envelope.get("total_capped", False), following, note)


PHENOMENA_TABLE = tables.register(Table(
    "phenomena", "Phenomena",
    [Column("name", "Name"), Column("type", "Type"), Column("descriptor", "Descriptor"),
     Column("radius", "Radius"), Column("sector", "Sector"), Column("placed", "On Galaxy Map"),
     Column("members", "Members", sortable=False)],
    _phenomena_load,
    facets=[Facet("type", "Type"), Facet("descriptor", "Descriptor")],
    noun=("phenomenon", "phenomena"), keyset=True,
    jump={"name": "A scattered phenomenon's number, such as 5000"},
))


@bp.route("/phenomena")
@slow_route("Phenomena", section="phenomena", breadcrumbs=lambda: [crumb("Phenomena")],
            description="Every nebula, black hole, neutron star and other exotic phenomenon in this charted galaxy.")
def phenomena():
    """Every exotic phenomenon, placed on the galaxy map or not, as a data
    table (UX.41): sort by any column, filter by type and descriptor, and
    scroll through the pages."""
    table = tables.render(PHENOMENA_TABLE, request.path, anchor="phenomena-list")
    return render_page(
        "phenomena.html",
        title="Phenomena",
        section="phenomena",
        breadcrumbs=[crumb("Phenomena")],
        description="Every nebula, black hole, neutron star and other exotic phenomenon in this charted galaxy.",
        table=table,
        total=table["total"],
    )


_SCATTER_SUBTYPE_NOTES = {
    "supermassive": "a supermassive black hole",
    "intermediate": "an intermediate-mass black hole",
    "stellar": "a stellar-mass black hole",
}

_SCATTER_ORIGINS = {
    "hypervelocity_star": "Ejected from the central black hole, flying outward",
}


def _uncharted_fields(detail):
    """`[(label, text, url)]` for an uncharted phenomenon's data table (MAP.211): what the scatter drew and
    nothing it has not. A star's own spectral class, mass and luminosity are drawn when its sector is made."""
    kind = detail["type"]
    subtype = detail.get("subtype")
    label = TYPE_LABELS.get(kind, kind)
    fields = [("Type", f"{label} ({_title_case(subtype)})" if subtype else label, None),
              ("Status", "Uncharted: its sector has not been generated, so it has no record of its own yet", None)]
    if kind == "hypervelocity_star":
        star = HYPERVELOCITY_EXPECTED_STAR
        fields.append((
            "Expected star",
            Markup(f"{star['star_type']}-type, a young main-sequence star of about "
                   f"{_kelvin(star['temperature_k'])} and {format_number(star['luminosity_sol'], ',.0f')} "
                   "solar luminosities; the star itself is drawn when its sector is generated"), None))
    if detail.get("mass_solar") is not None:
        fields.append(("Mass", _solar_mass(detail["mass_solar"]), None))
    if detail.get("member_count"):
        fields.append(("Members", f"{detail['member_count']} (see below)", None))
        fields.append(("System mass", _solar_mass(detail["system_mass_solar"]), None))
    if detail.get("speed_kms") is not None:
        fields.append(("Speed", speed_kms_html(detail["speed_kms"]), None))
        if detail.get("flight_years"):
            fields.append(("Time in flight", period_years_html(detail["flight_years"]), None))
    if kind in _SCATTER_ORIGINS:
        fields.append(("Origin", _SCATTER_ORIGINS[kind], None))
    shell, band, wedge = detail["address"]
    fields.append(("Sector cell", f"Shell {shell}, band {band}, wedge {wedge}", None))
    fields.append(("Position x, y, z", ", ".join(format_number(detail[axis], ",.1f") for axis in "xyz") + " pc", None))
    return fields


@bp.route("/uncharted-phenomenon/<int:scatter_id>")
def uncharted_phenomenon(scatter_id):
    """
    A phenomenon the scatter placed (MAP.211). While its sector is not generated this is its page: what is known
    of it, with the kind of star a hypervelocity star will be. Once the sector is generated it is the object
    itself, so the page sends the visitor to that object's page: the star system for a hypervelocity star, the
    phenomenon page for anything else.
    """
    try:
        detail = apiclient.get_uncharted_phenomenon(db_name(), scatter_id)
    except apiclient.NotFoundError:
        abort(404)
    target = detail.get("target")
    if detail["built"] and target:
        if target["kind"] == "system":
            return redirect(page_url("system", system_id=target["id"]), code=302)
        return redirect(page_url("phenomenon", phenomenon_type=target["kind"], phenomenon_id=target["id"]), code=302)
    label = TYPE_LABELS.get(detail["type"], detail["type"])
    name = f"{label} (uncharted)" if not detail["built"] else f"{label} (generated)"
    sector = None
    if detail.get("sector_id") is not None:
        sector = {"name": detail.get("sector_name") or "Sector", "url": page_url("sector", sector_id=detail["sector_id"])}
    return render_page(
        "uncharted_phenomenon.html",
        title=name,
        section="phenomena",
        breadcrumbs=[crumb("Phenomena", "phenomena"), crumb(name)],
        description=f"{name}, scattered across this galaxy and not yet generated.",
        type_label=label,
        distance=measured_html("pc", detail["distance_from_core_pc"]),
        sector=sector,
        built_note=None if not detail["built"] else
        "Its sector has been generated, but this record could not be matched to what was made from it.",
        fields=_uncharted_fields(detail),
        members=members_view(detail),
    )


# Every distance and non-body radius goes through the distance ladder
# (fmt.format_distance_*); a rogue planet's radius is a body radius, so it
# is km in scientific notation like every planet's.
_ly = distance_ly_html
_km = distance_km_html


def _quantity(unit, companions=True):
    """A field formatter for a measurement (UX.78): the value in the registry unit `unit`, in a span the reader's unit
    system re-renders (`fmt.measured_html`)."""
    return lambda value: measured_html(unit, value, companions)


_solar_mass = _quantity("Msun")
_power = _quantity("W")
_kelvin = _quantity("K", companions=False)


def _body_radius(value):
    return to_plain_text(format_star_radius(value))


_light_days = _quantity("ld")


def _bool_text(value):
    return "Yes" if value else "No"


def _progenitor_text(value):
    # "Type Ia" is stored correctly cased; _title_case would lowercase "Ia".
    return value if value == "Type Ia" else _title_case(value)


def _pressure_text(value):
    """The pressure on the ladder, "None" for no air at all, a trace in its own unit alone."""
    if not value:
        return "None"
    if value < 100:
        return Markup(f"{measured_html('Pa', value, companions=False)} (a trace)")
    return measured_html("Pa", value)


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
    ("dominant_species", "Dominant Species", lambda v: Markup(chem.mark_html(v)) if v else None),
    ("density_cm3", "Density", lambda v: measured_html("cm-3", v) if v else None),
    ("temperature_k", "Gas Temperature", lambda v: f"{format_number(v, ',.0f')} K" if v else None),
    ("extinction_av", "Extinction", lambda v: f"{format_number(v, ',.2f')} magnitudes (visual)" if v else None),
]


def _rogue_mass_bin_text(value):
    return ROGUE_MASS_CLASS_NAMES.get(value, value)


def _optional(fmt):
    """A formatter that leaves a `None` (not applicable) row out."""
    return lambda v: fmt(v) if v is not None else None


def _field_composition_text(rows):
    """An asteroid field's composition rows (DB.2) as a phrase: "high
    concentrations of iron and trace amounts of platinum"."""
    return format_composition_summary([(row["component"], row["concentration"]) for row in rows]) if rows else None


def _comet_composition_text(components):
    """An interstellar comet's composition rows (DB.2) as a phrase."""
    return format_comet_composition_summary(components) if components else None


_SPEED = ("galactic_orbital_speed_kms", "Galactic Orbital Speed", speed_kms_html)
_PERIOD = ("galactic_orbital_period_gy", "Galactic Orbital Period", _quantity("Gyr"))

_EINSTEIN_RADIUS = ("einstein_radius_m", "Einstein Radius (seen from the reference radius, against the galactic centre)",
                    _quantity("m"))
"""The derived row of a compact object's page (GEN.154): see `_derived_fields`."""

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
        ("composition", "Composition", _field_composition_text),
        _SPEED, _PERIOD,
    ],
    "black_hole": [
        ("mass_class", "Class", _title_case),
        ("mass_solar", "Mass", _solar_mass),
        ("event_horizon_radius_km", "Event Horizon Radius", _km),
        ("spin", "Spin (dimensionless)", lambda v: f"{v:.3f}"),
        ("has_accretion_disk", "Accretion Disk", _bool_text),
        ("temperature_k", "Hawking Temperature", _kelvin),
        ("luminosity_w", "Luminosity", _power),
        ("age_gy", "Age", lambda v: f"{format_number(v, ',.2f')} Gyr"),
        _EINSTEIN_RADIUS,
        (_SPEED[0], _SPEED[1], _optional(_SPEED[2])),
        (_PERIOD[0], _PERIOD[1], _optional(_PERIOD[2])),
    ],
    "neutron_star": [
        ("mass_solar", "Mass", _solar_mass),
        ("radius_km", "Radius", _km),
        ("spin_period_ms", "Spin Period", _quantity("ms")),
        ("magnetic_field_gauss", "Magnetic Field", _quantity("G")),
        ("pulsar_type", "Pulsar Type", _title_case),
        ("surface_temperature_k", "Surface Temperature", _kelvin),
        ("luminosity_w", "Luminosity", _power),
        ("age_gy", "Age", lambda v: f"{format_number(v, ',.2f')} Gyr"),
        _EINSTEIN_RADIUS,
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
        ("planet_class", "Planet Class", str),
        ("mass_bin", "Mass Class", _rogue_mass_bin_text),
        ("mass_kg", "Mass", _quantity("kg")),
        ("radius_km", "Radius", _body_radius),
        ("gravity_g", "Surface Gravity", _quantity("g0")),
        ("composition", "Composition", str),
        ("has_internal_heat", "Geologically Active", _bool_text),
        ("has_moons", "Has Moons", _bool_text),
        # Surface conditions (schema v48, rogueSurface): no star, so only
        # its own heat. A giant's temperature and pressure are at 1 bar.
        ("age_gy", "Age", lambda v: f"{format_number(v, ',.2f')} Gyr"),
        ("internal_heat_flux_w_m2", "Internal Heat Flow", _quantity("W_m2")),
        ("effective_temperature_k", "Effective Temperature", _kelvin),
        ("surface_regime", "Surface", lambda v: SURFACE_REGIME_LABELS.get(v, v)),
        ("surface_temperature_k", "Surface Temperature", _quantity("K")),
        ("surface_pressure_pa", "Surface Pressure", _pressure_text),
        ("ice_shell_thickness_km", "Ice Thickness", lambda v: f"{format_number(v, ',.1f')} km" if v else None),
        ("ocean_depth_km", "Liquid Ocean Depth", lambda v: f"{format_number(v, ',.0f')} km" if v else None),
        ("hp_ice_km", "High-Pressure Ice Below", lambda v: f"{format_number(v, ',.0f')} km" if v else None),
        # Habitability (schema v153): the PHI-6 score of the surface a visitor would find
        ("phi4", "PHI-6", lambda v: f"{v:.2f}"),
        _EINSTEIN_RADIUS,
        _SPEED, _PERIOD,
    ],
    "interstellar_comet": [
        ("nucleus_diameter_km", "Nucleus Diameter", _km),
        ("velocity_kms", "Velocity", speed_kms_html),
        ("is_active", "Active", _bool_text),
        ("composition", "Composition", _comet_composition_text),
        _SPEED, _PERIOD,
    ],
    # No galactic-orbit rows: a quasar sits at the galactic center.
    "quasar": [
        ("black_hole_mass_solar", "Black Hole Mass", _solar_mass),
        ("event_horizon_radius_km", "Event Horizon Radius", _km),
        ("luminosity_w", "Luminosity", _power),
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
    "planet_class": "planet",
}
"""dict: The columns whose value is a class, and the class type
(`planetgen/web/lib/classref.py`) whose page it links to. An asteroid field's "C3"
links by its letter."""


def _with_galactic_speed(detail):
    """`detail` plus `galactic_orbital_speed_kms`, the length of the movement vector it stores (GEN.336): the page shows
    a speed, the table keeps the vector."""
    parts = [detail.get(f"velocity_{axis}_kms") for axis in "xyz"]
    if any(part is None for part in parts):
        return detail
    return {**detail, "galactic_orbital_speed_kms": math.sqrt(sum(part * part for part in parts))}


def _with_derived_fields(phenomenon_type, detail):
    """`detail` plus the values worked out from it rather than stored: the
    Einstein radius (GEN.154) of a black hole, neutron star or rogue planet
    that has a galactic position, left out for a lens beyond the galaxy's reference radius from the centre
    (`physics/lensing.py`, GEN.404)."""
    if phenomenon_type in ("black_hole", "neutron_star"):
        mass_kg = Quantity("mass", detail.get("mass_solar") or 0, "Msun").kg
    elif phenomenon_type == "rogue_planet":
        mass_kg = detail.get("mass_kg") or 0
    else:
        return detail
    radius_pc = detail.get("galactic_radius_pc")
    observer_pc = detail.get("reference_radius_pc")
    if radius_pc is None or observer_pc is None:
        return detail
    einstein_m = einstein_radius_against_galactic_centre_m(mass_kg, radius_pc, observer_pc)
    if einstein_m is None:
        return detail
    return {**detail, "einstein_radius_m": einstein_m}


def phenomenon_fields(phenomenon_type, detail):
    """`[(label, text, url)]` for the data table; `url` is the value's
    class page, or `None`. Needs a request context."""
    fields = []
    if phenomenon_type == "black_hole" and not detail.get("has_accretion_disk") and detail.get("mass_solar"):
        # Stored before GEN.82 as a flat 0: show the Hawking values.
        detail = dict(detail)
        if not detail.get("temperature_k"):
            detail["temperature_k"] = hawking_temperature_k(detail["mass_solar"])
        if not detail.get("luminosity_w"):
            detail["luminosity_w"] = hawking_luminosity_w(detail["mass_solar"])
    detail = _with_derived_fields(phenomenon_type, _with_galactic_speed(detail))
    for column, label, formatter in FIELD_SPECS.get(phenomenon_type, []):
        raw = detail.get(column)
        if raw is None:
            continue
        text = formatter(raw)
        if text is not None:
            url = class_url(CLASS_COLUMNS[column], raw) if column in CLASS_COLUMNS else None
            fields.append((label, text, url))
    return fields


BODY_REF_KINDS = ("star", "planet", "moon", "belt", "comet")
"""tuple: The reference kinds `/object/<ref>` sends to a row on the system page (NAV.8)."""


def object_url(kind, object_id):
    """
    Where an object reference opens (NAV.8): a sector, system or phenomenon on its own page, and a star, planet, moon,
    belt or comet on the page of its own (`/object/planet:<ID>`, `web/item_pages.py`).
    """
    if kind == "sector":
        return page_url("sector", sector_id=object_id)
    if kind == "system":
        return page_url("system", system_id=object_id)
    if kind in BODY_REF_KINDS:
        return item_pages.item_url(f"{kind}:{object_id}")
    return page_url("phenomenon", phenomenon_type=kind, phenomenon_id=object_id)


@bp.route("/object/<ref>", methods=["GET", "POST"])
def object_page(ref):
    """Any object reference (`planet:12`, `nebula:3`, `system:` and a printed ID): the page the site has for the kind, or,
    for a star, planet, moon, belt or comet, its own page (UX.143, `web/item_pages.py`). A POST is an admin Regenerate or
    Delete of that body."""
    return item_pages.object_page_response(ref)


@bp.route("/object/<ref>/tree")
def object_tree(ref):
    """One page of an object's children as the rows of a summary tree level (UX.172, `static/objecttree.js`)."""
    return item_pages.tree_response(ref)


object_tree.json_only = True  # not a page: tests/test_web_a11y.py skips it




def _phenomenon_facts(phenomenon_type, detail):
    """The paragraph under "About" on a phenomenon page (UX.185): the sentences the generator's description says about the
    kind, built from the stored row (`describe/stored_phenomena.py`); empty for a kind with none."""
    document = stored_phenomenon_document(phenomenon_type, detail)
    if document is None:
        return ""
    return "".join(blocks_html([block for block in piece["blocks"] if block["type"] == "prose"])
                   for piece in sections(document))


def _phenomenon_peer(phenomenon_type, peer):
    """`{"name", "url"}` of the previous or next phenomenon of the type (`phenomenon_detail`'s `previous`/`next`), or `None`."""
    if not peer:
        return None
    return {"name": peer["name"], "url": page_url("phenomenon", phenomenon_type=phenomenon_type, phenomenon_id=peer["id"])}


@bp.route("/phenomenon/<phenomenon_type>/<uid:phenomenon_id>", methods=["GET", "POST"])
def phenomenon(phenomenon_type, phenomenon_id):
    """One phenomenon: its view (a render, the AU-scale diagram, or none
    for an asteroid field; see `phenomenonrender.view_kind`) and its data
    table. A POST is an admin Delete or Regenerate (`web/edit_actions.py`)."""
    if phenomenon_type not in TYPE_LABELS:
        abort(404)
    if request.method == "POST":
        return edit_actions.handle_post(
            {(phenomenon_type, phenomenon_id)},
            page_url("phenomenon", phenomenon_type=phenomenon_type, phenomenon_id=phenomenon_id),
            after_delete={phenomenon_type: page_url("phenomena")})
    detail = apiclient.get_phenomenon(db_name(), phenomenon_type, phenomenon_id)
    type_label = TYPE_LABELS[phenomenon_type]

    distance = None
    radius_pc = detail.get("galactic_radius_pc")
    if radius_pc is not None:
        distance = measured_html("pc", radius_pc)
    sector = None
    if detail.get("sector_id") is not None:
        sector = {"name": detail.get("sector_name") or "Sector",
                  "url": page_url("sector", sector_id=detail["sector_id"])}

    kind = view_kind(phenomenon_type)
    if kind == "render":
        map_html = render_phenomenon_view_panel(phenomenon_type, detail)
    elif kind == "nebula":
        map_html = render_nebula_view_panel(detail)
    elif kind == "map":
        map_html = render_phenomenon_map_panel(
            phenomenon_type, detail["name"], detail.get("radius_ly") or 0,
        )
    else:
        map_html = ""
    phenomenon_nav = pick_nav(phenomenon_type, detail["id"], request.args)
    return render_page(
        "phenomenon.html",
        title=detail["name"],
        section="phenomena",
        breadcrumbs=[crumb("Phenomena", "phenomena"), crumb(detail["name"])],
        description=f"{detail['name']}, a {type_label.lower()} in this charted galaxy.",
        type_label=type_label,
        distance=distance,
        sector=sector,
        octant=detail.get("quadrant"),
        nearest_html=trusted_html(nearest_systems_html(detail["nearest"], _system_url))
        if detail.get("nearest") else None,
        inside=_phenomenon_inside(detail),
        location=item_pages.sector_location(detail.get("sector_id"), detail.get("sector_name"),
                                            _phenomenon_inside_cloud(detail)),
        contents=_phenomenon_contents(phenomenon_type, detail["id"]),
        nav_links=phenomenon_nav[0],
        nearby_url=nearby_url(endpoint(phenomenon_type, phenomenon_id)),
        pick=phenomenon_nav[1],
        galaxy_url=page_url("sector_on_galaxy_map", sector_id=detail["sector_id"])
        if detail.get("sector_id") is not None else None,
        bookmark=bookmark(phenomenon_type, endpoint(phenomenon_type, detail["id"]), detail["name"],
                          page_url("phenomenon", phenomenon_type=phenomenon_type, phenomenon_id=detail["id"])),
        view_kind=kind,
        map_html=trusted_html(map_html),
        fields=phenomenon_fields(phenomenon_type, detail),
        facts=trusted_html(_phenomenon_facts(phenomenon_type, detail)),
        previous=_phenomenon_peer(phenomenon_type, detail.get("previous")),
        following=_phenomenon_peer(phenomenon_type, detail.get("next")),
        **_with_keep_out(item_pages.sheet_context(apiclient.get_object_sheet(db_name(), f"{phenomenon_type}:{phenomenon_id}")),
                         f"{phenomenon_type}:{phenomenon_id}"),
        # A system's own black hole or neutron star is edited with its system.
        can_edit=edit_actions.can_edit(current_admin()) and detail.get("star_id") is None,
        edit_target=f"{phenomenon_type}:{detail['id']}",
    )
