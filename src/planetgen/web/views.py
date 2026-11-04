# planetgen/web/views.py

"""
The routes of the Flask-served pages. First pages moved from CGI: the
home page (was `index.py`/`browse.py`), `/sectors` and `/systems` (the
two halves of `browse.py`), and `/search` (was `search.py`).
"""

import re

from flask import abort, flash, get_flashed_messages, redirect, request, url_for

from planetgen.web.lib import apiclient
from planetgen.web.lib.slowroute import slow_route
from planetgen.web.lib.datatable import Column, Facet, Result, Table, plain, read_from_start
from planetgen.web.lib.fmt import (
    format_density, format_distance_ly, format_luminosity_sol, format_members, format_number, format_mass_solar,
    measured_html)
from planetgen.web.scatter_members_view import members_view
from planetgen.web.maps.galaxymap import sector_quadrant

from planetgen.api.limiter import page_limit

from . import bp, csrf, tables
from . import searchpage
from .helpers import crumb, current_admin, db_name, page_url, render_page


def _binary_filter(values):
    """The `binary` argument for `yes`/`no` menu choices: one choice filters, none or both don't."""
    chosen = {value for value in values if value in ("yes", "no")}
    return None if len(chosen) != 1 else chosen == {"yes"}


def _sector_row(sector):
    """One Sectors table row: a list of cells (`lib/datatable.py`)."""
    quadrant = sector_quadrant(sector["center_x_pc"], sector["center_y_pc"]) if sector["placed"] else None
    return [
        {"text": sector["name"], "href": page_url("sector", sector_id=sector["id"])},
        {"text": str(sector["system_count"])},
        {"text": plain(format_density(sector["edge_ly"], sector["system_count"]))},
        {"text": f"Quadrant {quadrant}", "href": page_url("galaxy", quadrant=quadrant)} if quadrant
        else {"text": "Unplaced", "muted": True},
        {"text": plain(format_distance_ly(sector.get("galactic_radius_ly")))},
    ]


def _sectors_load(state, limit, offset, want_facets, after=None):
    envelope = apiclient.get_sectors(
        db_name(), limit=limit, offset=offset, after=after, sort=state.sort, descending=state.descending,
        quadrants=state.filters["quadrant"], facets=want_facets)
    facets = None
    if want_facets:
        facets = {"quadrant": [
            {"value": option["value"],
             "label": "Unplaced" if option["value"] == "unplaced" else f"Quadrant {option['value']}",
             "count": option["count"]} for option in envelope["facets"]["quadrant"]]}
    items = envelope["items"]
    return Result([_sector_row(sector) for sector in items], envelope["total"], facets,
                  envelope.get("total_capped", False), items[-1]["id"] if items else None)


SECTORS_TABLE = tables.register(Table(
    "sectors", "Sectors",
    [Column("name", "Name"), Column("systems", "Systems"), Column("density", "Density"),
     Column("position", "Galaxy Position"), Column("distance", "Distance from core")],
    _sectors_load, facets=[Facet("quadrant", "Quadrant")], prefix="sectors_", default_sort="distance",
    noun=("sector", "sectors"), keyset=True,
))

_BINARY_LABELS = {"yes": "Binary", "no": "Single star"}
_PLACEMENT_LABELS = {"sector": "In a sector", "standalone": "Standalone"}


def _system_row(system, with_sector):
    """One Systems table row (`with_sector`: the All Systems columns, else the Standalone ones)."""
    row = [{"text": system["name"], "href": page_url("system", system_id=system["id"])}]
    if with_sector:
        row.append({"text": system["sector_name"], "href": page_url("sector", sector_id=system["sector_id"])}
                   if system.get("sector_id") else {"text": "Standalone"})
        row.append({"text": system.get("quadrant") or "\u2013"})
    row.append({"text": str(system["star_count"])})
    row.append({"text": str(system["object_count"])})
    row.append({"text": "Yes" if system["has_habitable_world"] else "No"})
    row.append({"text": system["star_summary"] or "\u2013"})
    return row


def _facet_options(options, labels=None):
    return [{"value": option["value"], "label": (labels or {}).get(option["value"], option["value"]),
             "count": option["count"]} for option in options]


def _all_systems_load(state, limit, offset, want_facets, after=None):
    placement = state.filters["placement"]
    envelope = apiclient.get_systems(
        db_name(), limit=limit, offset=offset, after=after, sort=state.sort, descending=state.descending,
        binary=_binary_filter(state.filters["binary"]),
        placement=placement[0] if len(placement) == 1 and placement[0] in _PLACEMENT_LABELS else None,
        octants=state.filters["octant"], facets=want_facets, start=state.start)
    facets = None
    if want_facets:
        facets = {
            "placement": _facet_options(envelope["facets"]["placement"], _PLACEMENT_LABELS),
            "binary": _facet_options(envelope["facets"]["binary"], _BINARY_LABELS),
            "octant": _facet_options(envelope["facets"]["octant"]),
        }
    items = envelope["items"]
    return Result([_system_row(system, True) for system in items], envelope["total"], facets,
                  envelope.get("total_capped", False), items[-1]["id"] if items else None)


ALL_SYSTEMS_TABLE = tables.register(Table(
    "systems", "All systems",
    [Column("name", "Name"), Column("sector", "Sector"), Column("octant", "Octant"),
     Column("stars", "Stars"), Column("objects", "Objects", sortable=False),
     Column("habitable", "Habitable world", sortable=False), Column("star_type", "Star type", sortable=False)],
    _all_systems_load,
    facets=[Facet("placement", "Where"), Facet("binary", "Stars"), Facet("octant", "Octant")],
    prefix="systems_", noun=("system", "systems"), keyset=True, jump={"name": "A letter or the start of a name", "sector": "A sector name or its first letters"},
))


_FLASH_UNCHARTED = "uncharted"

_POPULATION_LABELS = {"young": "Young disk", "intermediate": "Disk", "old": "Old disk", "bulge": "Bulge"}


def _uncharted_row(star, can_generate):
    """One Uncharted stars table row (UX.87): the scattered star's own data, where it is and a
    Generate button for an admin."""
    label = f"Uncharted {star['star_type']} star {star['id']}"
    coords = ", ".join(f"{star[axis]:.2f}" for axis in ("x", "y", "z"))
    local = ", ".join(f"{star[f'local_{axis}']:.2f}" for axis in ("x", "y", "z"))
    generate = {"text": ""}
    if can_generate:
        generate = {"text": "", "form": {
            "action": url_for("web.generate_uncharted_system", bright_star_id=star["id"]), "button": "Generate",
            "label": f"Generate {label} by itself",
            "fields": [[csrf.FIELD_NAME, csrf.csrf_token()]]}}
    return [
        {"text": label, "href": page_url("uncharted_star", bright_star_id=star["id"])},
        {"text": f"{format_luminosity_sol(star['luminosity_sol'])}"},
        {"text": f"{star['star_type']} ({_POPULATION_LABELS.get(star['population'], star['population'])})"},
        {"text": f"{format_number(round(star['temperature_k']))} K"},
        {"text": f"{star['age_gy']:.2f} Gyr"},
        {"text": f"{star['designation']} (shell {star['shell_index']}, band {star['band_index']}, "
                 f"wedge {star['wedge_index']})",
         "href": page_url("galaxy", sector=star["designation"], open="1")},
        {"text": f"{coords} pc"},
        {"text": f"{local} pc"},
        {"text": format_members(star.get("member_count"), star.get("system_mass_solar"))},
        generate,
    ]


def _uncharted_load(state, limit, offset, want_facets, after=None):
    # Brightest first is the table's default, so its ascending order is the API's descending one.
    descending = state.descending if state.sort != "luminosity" else not state.descending
    envelope, note = read_from_start(
        lambda start: apiclient.get_uncharted_systems(
            db_name(), limit=limit, offset=offset, sort=state.sort, descending=descending, after=after, start=start),
        state.start)
    can_generate = current_admin() is not None
    items = envelope["items"]
    return Result([_uncharted_row(star, can_generate) for star in items], envelope["total"], None,
                  envelope.get("total_capped", False), items[-1]["id"] if items else None, note)


UNCHARTED_SYSTEMS_TABLE = tables.register(Table(
    "uncharted_systems", "Uncharted stars",
    [Column("name", "Name", sortable=False), Column("luminosity", "Luminosity"), Column("type", "Star"),
     Column("temperature", "Temperature"), Column("age", "Age"), Column("sector", "Sector and cell"),
     Column("position", "Galaxy x, y, z", sortable=False), Column("local", "In the sector x, y, z", sortable=False),
     Column("members", "Members", sortable=False), Column("generate", "Generate", sortable=False)],
    _uncharted_load, prefix="uncharted_", noun=("uncharted star", "uncharted stars"), keyset=True,
    jump={"luminosity": "Solar luminosities, such as 25", "temperature": "Kelvin, such as 6000",
          "type": "A star type, such as G2", "sector": "A shell number, such as 12", "age": "Gyr, such as 3.5"},
))


def _standalone_url():
    """The Systems list filtered to the standalone systems (UX.55)."""
    placement = ALL_SYSTEMS_TABLE.facets[0].param
    return page_url("systems", _anchor="all-systems", **{placement: "standalone"})


def _systems_total(**filters):
    return apiclient.get_systems(db_name(), limit=1, offset=0, **filters)["total"]


@bp.route("/")
def index():
    """Home, the front door (UX.55): the counts of what is here, each a link
    to its list, and the Galaxy Map, NAV and Search."""
    db = db_name()
    sectors = apiclient.get_sectors(db, limit=1, offset=0)
    systems = apiclient.get_systems(db, limit=1, offset=0)
    phenomena = apiclient.get_phenomena(db, limit=1, offset=0)
    standalone = apiclient.get_systems(db, limit=1, offset=0, sector_id="none")
    # A count the API stopped at its cap (`total_capped`) is shown "9,999+", as the lists show it (UX.110).
    counts = [
        (sectors, ("sector", "sectors"), page_url("sectors")),
        (systems, ("system", "systems"), page_url("systems")),
        (phenomena, ("phenomenon", "phenomena"), page_url("phenomena")),
        (standalone, ("standalone system", "standalone systems"), _standalone_url()),
    ]
    return render_page(
        "index.html",
        title="Home",
        description="Every sector, star system and phenomenon in this charted galaxy.",
        counts=[{"total": envelope["total"], "capped": bool(envelope.get("total_capped")), "url": url,
                 "noun": nouns[0] if envelope["total"] == 1 and not envelope.get("total_capped") else nouns[1]}
                for envelope, nouns, url in counts],
        tools=[{"label": "Galaxy Map", "hint": "Fly through the whole galaxy.", "url": page_url("galaxy")},
               {"label": "Navigate", "hint": "Plan a route between two places.", "url": page_url("nav")},
               {"label": "Search", "hint": "Find a sector, system, star, planet or moon by name.",
                "url": page_url("search")}],
    )


@bp.route("/sectors")
def sectors():
    """Every sector, nearest the galactic core first."""
    return render_page(
        "sectors.html",
        title="Sectors",
        section="sectors",
        breadcrumbs=[crumb("Sectors")],
        description="Every sector in this charted galaxy.",
        sectors=tables.render(SECTORS_TABLE, request.path, anchor="sectors"),
    )


@bp.route("/systems")
def systems():
    """Every system with its sector and octant; the standalone ones (generated
    outside any sector) are a filter on the list (UX.55). `?uncharted=1`
    lists the scattered stars no system was built around instead (UX.87)."""
    show_uncharted = request.args.get("uncharted") == "1"
    if show_uncharted:
        table = tables.render(UNCHARTED_SYSTEMS_TABLE, request.path, anchor="all-systems", keep=("uncharted",))
    else:
        table = tables.render(ALL_SYSTEMS_TABLE, request.path, anchor="all-systems")
    return render_page(
        "systems.html",
        title="Uncharted stars" if show_uncharted else "Systems",
        section="systems",
        breadcrumbs=[crumb("Systems", "systems"), crumb("Uncharted stars")] if show_uncharted else [crumb("Systems")],
        description="Every star system in this charted galaxy.",
        all_systems=table,
        show_uncharted=show_uncharted,
        uncharted_url=page_url("systems", uncharted="1", _anchor="all-systems"),
        systems_url=page_url("systems", _anchor="all-systems"),
        can_generate=current_admin() is not None,
        messages=get_flashed_messages(category_filter=[_FLASH_UNCHARTED]),
        standalone_url=_standalone_url(),
        standalone_total=_systems_total(sector_id="none"),
        uncharted_total=apiclient.get_uncharted_systems(db_name(), limit=1, offset=0)["total"],
    )


@bp.route("/systems/uncharted/<int:bright_star_id>")
def uncharted_star(bright_star_id):
    """
    One scattered star no system was built around yet (UX.114): its own data, where it is and, for a star with company,
    its members and total mass. Once a system is built around it the page sends the visitor to that system.
    """
    try:
        detail = apiclient.get_uncharted_system(db_name(), bright_star_id)
    except apiclient.NotFoundError:
        abort(404)
    if detail.get("system_id") is not None:
        return redirect(page_url("system", system_id=detail["system_id"]), code=302)
    name = f"Uncharted {detail['star_type']} star {detail['id']}"
    population = _POPULATION_LABELS.get(detail["population"], detail["population"])
    fields = [
        ("Star", f"{detail['star_type']} ({population})"),
        ("Luminosity", format_luminosity_sol(detail["luminosity_sol"])),
        ("Temperature", f"{format_number(round(detail['temperature_k']))} K"),
        ("Mass", format_mass_solar(detail["mass_solar"])),
        ("Age", f"{detail['age_gy']:.2f} Gyr"),
        *([("Rotation period", measured_html("h", detail["rotation_period_hours"]))]
          if detail.get("rotation_period_hours") is not None else []),
        ("Sector cell", f"{detail['designation']} (shell {detail['shell_index']}, band {detail['band_index']}, "
                        f"wedge {detail['wedge_index']})"),
        ("Position x, y, z", ", ".join(f"{detail[axis]:.2f}" for axis in ("x", "y", "z")) + " pc"),
    ]
    nebula = detail.get("hosted_nebula")
    if nebula:
        fields.append(("Lights a nebula", f"{nebula['name']}, a {nebula['nebula_type']} nebula of class "
                                          f"{nebula['nebula_class']} ({nebula['class_name']}), "
                                          f"radius {format_distance_ly(nebula['radius_ly'])}"))
    if detail.get("member_count"):
        fields.insert(4, ("Members", f"{detail['member_count']} (see below)"))
        fields.insert(5, ("System mass", format_mass_solar(detail["system_mass_solar"])))
    return render_page(
        "uncharted_star.html",
        title=name,
        section="systems",
        breadcrumbs=[crumb("Systems", "systems"), crumb("Uncharted stars", "systems", uncharted="1"), crumb(name)],
        description=f"{name}, scattered across this galaxy and not yet generated.",
        fields=fields,
        members=members_view(detail),
        sector_url=page_url("galaxy", sector=detail["designation"], open="1"),
        designation=detail["designation"],
    )


@bp.route("/systems/uncharted/<int:bright_star_id>/generate", methods=["POST"])
def generate_uncharted_system(bright_star_id):
    """An admin's Generate button on an uncharted star (UX.87): builds that one system and opens it;
    the rest of its sector stays ungenerated."""
    back = redirect(page_url("systems", uncharted="1", _anchor="all-systems"), code=303)
    admin = current_admin()
    if admin is None or admin["must_change_credentials"]:
        return back
    try:
        result = apiclient.admin_edit(request.headers.get("Cookie"), db_name(), "POST",
                                      f"/uncharted-systems/{bright_star_id}/generate")
    except apiclient.NotFoundError:
        flash("That star is no longer waiting: it may have been generated already.", _FLASH_UNCHARTED)
        return back
    except apiclient.ApiError as exc:
        if exc.status_code is None or exc.status_code >= 500:
            raise
        flash(re.sub(r"^planetGen API error \(\d+\): ", "", str(exc)), _FLASH_UNCHARTED)
        return back
    if result.get("id"):
        return redirect(page_url("system", system_id=result["id"]), code=303)
    flash("The system is being generated; reload in a moment to see it.", _FLASH_UNCHARTED)
    return back


@bp.route("/search")
@page_limit("search")
@slow_route("Search", breadcrumbs=lambda: [crumb("Search")],
            description="Search this charted galaxy by name, size and type.")
def search():
    """
    Faceted search (was `search.py`): `?q=` searches every kind of name
    at once (the header search box), plus per-object name fields, size
    ranges and click-to-filter tags, each a GET parameter, so any search
    is a bookmarkable URL. See `web/searchpage.py`.
    """
    if searchpage.needs_canonical_redirect(request.args):
        state = searchpage.SearchState.from_args(request.args)
        return redirect(state.url(with_pages=True), code=302)
    state = searchpage.SearchState.from_args(request.args)
    data = searchpage.fetch(state)
    return render_page(
        "search.html",
        title="Search",
        breadcrumbs=[crumb("Search")],
        description="Search this charted galaxy by name, size and type.",
        state=state,
        lookup=searchpage.lookup_notice(state.q),
        name_fields=searchpage.NAME_FIELDS,
        size_entities=searchpage.SIZE_ENTITIES,
        value_groups=searchpage.value_form(state),
        star_choices=searchpage.STAR_COUNT_CHOICES,
        autocomplete=data["autocomplete"],
        tag_groups=searchpage.tag_groups(state, data["facets"]),
        chips=searchpage.active_filters(state, data["facet_labels"]),
        panels=searchpage.result_panels(state, data["results"]) if state.active() else None,
    )
