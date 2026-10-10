# planetgen/web/views.py

"""
The routes of the Flask-served pages. First pages moved from CGI: the
home page (was `index.py`/`browse.py`), `/sectors` and `/systems` (the
two halves of `browse.py`), and `/search` (was `search.py`).
"""

import re

from flask import flash, get_flashed_messages, redirect, request, url_for

from planetgen.web.lib import apiclient
from planetgen.web.lib.datatable import Column, Facet, Result, Table, plain
from planetgen.web.lib.fmt import format_density, format_distance_ly, format_number
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


def _sectors_load(state, limit, offset, want_facets):
    envelope = apiclient.get_sectors(
        db_name(), limit=limit, offset=offset, sort=state.sort, descending=state.descending,
        quadrants=state.filters["quadrant"], facets=want_facets)
    facets = None
    if want_facets:
        facets = {"quadrant": [
            {"value": option["value"],
             "label": "Unplaced" if option["value"] == "unplaced" else f"Quadrant {option['value']}",
             "count": option["count"]} for option in envelope["facets"]["quadrant"]]}
    return Result([_sector_row(sector) for sector in envelope["items"]], envelope["total"], facets)


SECTORS_TABLE = tables.register(Table(
    "sectors", "Sectors",
    [Column("name", "Name"), Column("systems", "Systems"), Column("density", "Density"),
     Column("position", "Galaxy Position"), Column("distance", "Distance from core")],
    _sectors_load, facets=[Facet("quadrant", "Quadrant")], prefix="sectors_", default_sort="distance",
    noun=("sector", "sectors"),
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
    row.append({"text": "Yes" if system["is_binary"] else "No"})
    row.append({"text": system["star_summary"] or "\u2013"})
    return row


def _facet_options(options, labels=None):
    return [{"value": option["value"], "label": (labels or {}).get(option["value"], option["value"]),
             "count": option["count"]} for option in options]


def _all_systems_load(state, limit, offset, want_facets):
    placement = state.filters["placement"]
    envelope = apiclient.get_systems(
        db_name(), limit=limit, offset=offset, sort=state.sort, descending=state.descending,
        binary=_binary_filter(state.filters["binary"]),
        placement=placement[0] if len(placement) == 1 and placement[0] in _PLACEMENT_LABELS else None,
        octants=state.filters["octant"], facets=want_facets)
    facets = None
    if want_facets:
        facets = {
            "placement": _facet_options(envelope["facets"]["placement"], _PLACEMENT_LABELS),
            "binary": _facet_options(envelope["facets"]["binary"], _BINARY_LABELS),
            "octant": _facet_options(envelope["facets"]["octant"]),
        }
    return Result([_system_row(system, True) for system in envelope["items"]], envelope["total"], facets)


ALL_SYSTEMS_TABLE = tables.register(Table(
    "systems", "All systems",
    [Column("name", "Name"), Column("sector", "Sector"), Column("octant", "Octant"), Column("binary", "Binary"),
     Column("star_type", "Star type", sortable=False)],
    _all_systems_load,
    facets=[Facet("placement", "Where"), Facet("binary", "Stars"), Facet("octant", "Octant")],
    prefix="systems_", noun=("system", "systems"),
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
        {"text": label},
        {"text": f"{format_number(star['luminosity_sol'])} L\u2609"},
        {"text": f"{star['star_type']} ({_POPULATION_LABELS.get(star['population'], star['population'])})"},
        {"text": f"{format_number(round(star['temperature_k']))} K"},
        {"text": f"{star['age_gy']:.2f} Gy"},
        {"text": f"{star['designation']} (ring {star['ring_index']}, layer {star['layer_index']}, "
                 f"slot {star['ring_slot_index']})",
         "href": page_url("galaxy", sector=star["designation"], open="1")},
        {"text": f"{coords} pc"},
        {"text": f"{local} pc"},
        generate,
    ]


def _uncharted_load(state, limit, offset, want_facets):
    # Brightest first is the table's default, so its ascending order is the API's descending one.
    descending = state.descending if state.sort != "luminosity" else not state.descending
    envelope = apiclient.get_uncharted_systems(
        db_name(), limit=limit, offset=offset, sort=state.sort, descending=descending)
    can_generate = current_admin() is not None
    return Result([_uncharted_row(star, can_generate) for star in envelope["items"]], envelope["total"], None)


UNCHARTED_SYSTEMS_TABLE = tables.register(Table(
    "uncharted_systems", "Uncharted stars",
    [Column("name", "Name", sortable=False), Column("luminosity", "Luminosity"), Column("type", "Star"),
     Column("temperature", "Temperature"), Column("age", "Age"), Column("sector", "Sector and cell"),
     Column("position", "Galaxy x, y, z", sortable=False), Column("local", "In the sector x, y, z", sortable=False),
     Column("generate", "Generate", sortable=False)],
    _uncharted_load, prefix="uncharted_", noun=("uncharted star", "uncharted stars"),
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
    counts = [
        (apiclient.get_sectors(db, limit=1, offset=0)["total"], ("sector", "sectors"), page_url("sectors")),
        (_systems_total(), ("system", "systems"), page_url("systems")),
        (apiclient.get_phenomena(db, limit=1, offset=0)["total"], ("phenomenon", "phenomena"),
         page_url("phenomena")),
        (_systems_total(sector_id="none"), ("standalone system", "standalone systems"), _standalone_url()),
    ]
    return render_page(
        "index.html",
        title="Home",
        description="Every sector, star system and phenomenon in this generated galaxy.",
        counts=[{"total": total, "noun": nouns[0] if total == 1 else nouns[1], "url": url}
                for total, nouns, url in counts],
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
        description="Every sector in this generated galaxy.",
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
        description="Every star system in this generated galaxy.",
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
    data = apiclient.get_search(
        db_name(), state.name_terms(), state.tags, sizes=state.sizes(), limit=searchpage.PAGE_SIZE,
        offsets=searchpage.api_offsets(state),
    )
    return render_page(
        "search.html",
        title="Search",
        breadcrumbs=[crumb("Search")],
        description="Search this generated galaxy by name, size and type.",
        state=state,
        name_fields=searchpage.NAME_FIELDS,
        size_entities=searchpage.SIZE_ENTITIES,
        autocomplete=data["autocomplete"],
        tag_groups=searchpage.tag_groups(state, data["facets"]),
        chips=searchpage.active_filters(state, data["facet_labels"]),
        panels=searchpage.result_panels(state, data["results"]) if state.active() else None,
    )
