# planetgen/web/views.py

"""
The routes of the Flask-served pages. First pages moved from CGI: the
home page (was `index.py`/`browse.py`), `/sectors` and `/systems` (the
two halves of `browse.py`), and `/search` (was `search.py`).
"""

from flask import redirect, request

from planetgen.web.lib import apiclient
from planetgen.web.lib.datatable import Column, Facet, Result, Table, plain
from planetgen.web.lib.fmt import format_density, format_distance_ly
from planetgen.web.maps.galaxymap import sector_quadrant

from planetgen.api.limiter import page_limit

from . import bp, tables
from . import searchpage
from .helpers import crumb, db_name, page_url, render_page


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


def _standalone_load(state, limit, offset, want_facets):
    envelope = apiclient.get_systems(
        db_name(), sector_id="none", limit=limit, offset=offset, sort=state.sort, descending=state.descending,
        binary=_binary_filter(state.filters["standalone_binary"]), facets=want_facets)
    facets = None
    if want_facets:
        facets = {"standalone_binary": _facet_options(envelope["facets"]["binary"], _BINARY_LABELS)}
    return Result([_system_row(system, False) for system in envelope["items"]], envelope["total"], facets)


STANDALONE_TABLE = tables.register(Table(
    "standalone-systems", "Standalone systems",
    [Column("name", "Name"), Column("binary", "Binary"), Column("star_type", "Star type", sortable=False)],
    _standalone_load, facets=[Facet("standalone_binary", "Stars")], prefix="standalone_",
    noun=("standalone system", "standalone systems"),
))


@bp.route("/")
def index():
    """Home: every sector and every standalone system, each a table of its own
    (`sectors_*` and `standalone_*` query parameters)."""
    sectors = tables.render(SECTORS_TABLE, request.path, anchor="sectors")
    systems = tables.render(STANDALONE_TABLE, request.path, anchor="standalone-systems")
    return render_page(
        "index.html",
        title="Home",
        description="Every sector and standalone star system in this generated galaxy.",
        sectors=sectors,
        systems=systems,
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
    """Every system with its sector and octant, and the standalone ones
    (generated outside any sector) in their own table."""
    return render_page(
        "systems.html",
        title="Systems",
        section="systems",
        breadcrumbs=[crumb("Systems")],
        description="Every star system in this generated galaxy.",
        all_systems=tables.render(ALL_SYSTEMS_TABLE, request.path, anchor="all-systems"),
        systems=tables.render(STANDALONE_TABLE, request.path, anchor="standalone-systems"),
    )


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
