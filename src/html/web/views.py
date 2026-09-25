# html/web/views.py

"""
The routes of the Flask-served pages. First pages moved from CGI: the
home page (was `index.py`/`browse.py`), `/sectors` and `/systems` (the
two halves of `browse.py`), and `/search` (was `search.py`).
"""

from flask import redirect, request

import apiclient
from fmt import format_density, format_distance_ly
from galaxymap import sector_quadrant
from pagination import fetch_page, parse_page

from . import bp
from . import searchpage
from .helpers import crumb, db_name, page_url, pager, render_page, trusted_html


def _sector_rows(sectors):
    """Template-ready rows for the Sectors table."""
    rows = []
    for sector in sectors:
        quadrant = None
        if sector["placed"]:
            quadrant = sector_quadrant(sector["center_x_pc"], sector["center_y_pc"])
        rows.append({
            "name": sector["name"],
            "url": page_url("sector", sector_id=sector["id"]),
            "system_count": sector["system_count"],
            "density": trusted_html(format_density(sector["edge_ly"], sector["system_count"])),
            "quadrant": quadrant,
            "quadrant_url": page_url("galaxy", quadrant=quadrant) if quadrant else None,
            "distance": trusted_html(format_distance_ly(sector.get("galactic_radius_ly"))),
        })
    return rows


def _system_rows(systems):
    """Template-ready rows for the Standalone Systems table."""
    return [{
        "name": system["name"],
        "url": page_url("system", system_id=system["id"]),
        "is_binary": bool(system["is_binary"]),
        "star_summary": system["star_summary"],
    } for system in systems]


def _sectors_panel(keep=None):
    """One page of sectors plus its pager (`?sectors_page=N`)."""
    db = db_name()
    envelope, page = fetch_page(
        lambda limit, offset: apiclient.get_sectors(db, limit=limit, offset=offset),
        parse_page(request.args.get("sectors_page")),
    )
    return {
        "rows": _sector_rows(envelope["items"]),
        "total": envelope["total"],
        "page": page,
        "pager": pager("sectors_page", page, envelope["total"], anchor="sectors",
                       label="Sector pages", keep=keep),
    }


def _systems_panel(keep=None):
    """One page of standalone systems plus its pager
    (`?standalone_page=N`)."""
    db = db_name()
    envelope, page = fetch_page(
        lambda limit, offset: apiclient.get_systems(db, sector_id="none", limit=limit, offset=offset),
        parse_page(request.args.get("standalone_page")),
    )
    return {
        "rows": _system_rows(envelope["items"]),
        "total": envelope["total"],
        "page": page,
        "pager": pager("standalone_page", page, envelope["total"], anchor="standalone-systems",
                       label="Standalone system pages", keep=keep),
    }


@bp.route("/")
def index():
    """Home: every sector and every standalone system, each paged on its
    own (`?sectors_page=N&standalone_page=M`)."""
    sectors_page = parse_page(request.args.get("sectors_page"))
    standalone_page = parse_page(request.args.get("standalone_page"))
    sectors = _sectors_panel(keep={"standalone_page": standalone_page})
    systems = _systems_panel(keep={"sectors_page": sectors_page})
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
        sectors=_sectors_panel(),
    )


@bp.route("/systems")
def systems():
    """Every standalone system (one generated outside any sector)."""
    return render_page(
        "systems.html",
        title="Systems",
        section="systems",
        breadcrumbs=[crumb("Systems")],
        description="Every standalone star system in this generated galaxy.",
        systems=_systems_panel(),
    )


@bp.route("/search")
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
