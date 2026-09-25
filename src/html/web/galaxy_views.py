# html/web/galaxy_views.py

"""
The Galaxy Map page (`/galaxy`, was `galaxy.py`) and the JSON tile
endpoint its script fetches as the camera moves (`/galaxy/tiles`, was
`galaxy_tiles.py`).

The page renders the 3D map panel (`lib/galaxymap3d.py`, drawn by
`static/galaxymap3d.js`) with the zoomed-all-the-way-out starting view's
tiles embedded, plus a plain-text table under it: a per-Quadrant summary,
or with `?quadrant=I|II|III|IV` that Quadrant's placed sectors nearest
the core first (paged with `?page=N`).

Tiles go through `lib/tilecache.py`'s disk cache in this (the WSGI)
process: only tiles not already cached reach the API (in-process), and
whatever the API returns is written to the cache. The browser keeps its
own copy in `localStorage`, keyed as before, so a visitor's cached tiles
survive the move.
"""

from flask import jsonify, request, url_for

import apiclient
from galaxymap import QUADRANT_LABELS, ring_bounds_ly, sector_quadrant, sector_ring
from galaxymap3d import initial_tile_request, render_galaxy_map3d_panel, view_radius_bounds
from pagination import page_slice, parse_page
from stellarObjects import log
from stellarObjects.program_constants import DEFAULT_SECTOR_EDGE_LY
from stellarObjects.utils import ly_to_pc, pc_to_ly
from tilecache import TileRequestError, fetch_tiles

from . import bp
from .helpers import crumb, db_name, page_url, pager, render_page, trusted_html

_ID_PLACEHOLDER = 987654321987
"""int: Stands in for a sector id while building the map's sector-link
template (`page_url` needs a real int for the `/sector/<int>` route), and
is then replaced by `{id}` for the script to fill in."""


def sector_url_template():
    """
    The URL of a sector page with `{id}` where the id goes, for
    `static/galaxymap3d.js`'s "View sector" link. Goes through `page_url`,
    so it follows the sector page wherever it lives (CGI or Flask).
    """
    return page_url("sector", sector_id=_ID_PLACEHOLDER).replace(str(_ID_PLACEHOLDER), "{id}")


def _parse_quadrant(raw):
    quadrant = (raw or "").strip().upper()
    return quadrant if quadrant in QUADRANT_LABELS else None


def _quadrant_summary_rows(sectors):
    """One row per Quadrant: placed sectors, total systems, Ring extent."""
    by_quadrant = {label: [] for label in QUADRANT_LABELS}
    for sector in sectors:
        by_quadrant[sector_quadrant(sector["x"], sector["y"])].append(sector)
    rows = []
    for label in QUADRANT_LABELS:
        members = by_quadrant[label]
        extent = None
        if members:
            _inner, outer_ly = ring_bounds_ly(max(sector_ring(s["shell_index"]) for s in members))
            extent = f"out to ~{outer_ly:,.0f} ly"
        rows.append({
            "label": label,
            "url": page_url("galaxy", quadrant=label, _anchor="galaxy-table"),
            "sector_count": len(members),
            "system_count": sum(s["system_count"] or 0 for s in members),
            "extent": extent,
        })
    return rows


def _quadrant_sector_rows(sectors, quadrant, page):
    """One page of the Quadrant's placed sectors, nearest the core first.
    Returns `(rows, page, total)`."""
    members = [s for s in sectors if sector_quadrant(s["x"], s["y"]) == quadrant]
    members.sort(key=lambda s: s["galactic_radius_pc"])
    page_members, page = page_slice(members, page)
    rows = [{
        "name": sector["name"],
        "url": page_url("sector", sector_id=sector["id"]),
        "ring": sector_ring(sector["shell_index"]),
        "distance_ly": pc_to_ly(sector["galactic_radius_pc"]),
        "system_count": sector["system_count"] or 0,
    } for sector in page_members]
    return rows, page, len(members)


@bp.route("/galaxy")
def galaxy():
    """The 3D Galaxy Map plus the Quadrant table (`?quadrant=`, `?page=`)."""
    db = db_name()
    quadrant = _parse_quadrant(request.args.get("quadrant"))

    sectors = apiclient.get_galaxy_sectors(db)
    galaxy_shape = apiclient.get_galaxy_shape(db)
    edge_pc = galaxy_shape["edge_pc"] if galaxy_shape else ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    _min_radius, max_radius = view_radius_bounds(edge_pc, galaxy_shape)
    tile_keys, density_key = initial_tile_request(max_radius, galaxy_shape is not None)
    initial_view = fetch_tiles(db, tile_keys, density_key)
    map_html = render_galaxy_map3d_panel(
        db, galaxy_shape, edge_pc, initial_view,
        fetch_path=url_for("web.galaxy_tiles"),
        sector_url=sector_url_template(),
    )

    context = {"quadrant": quadrant, "placed_count": len(sectors), "map_html": trusted_html(map_html)}
    if quadrant:
        rows, page, total = _quadrant_sector_rows(sectors, quadrant, parse_page(request.args.get("page")))
        context.update(
            sector_rows=rows,
            pager=pager("page", page, total, anchor="galaxy-table", label="Sector pages",
                        keep={"quadrant": quadrant}),
        )
        title = f"Galaxy Map: Quadrant {quadrant}"
        breadcrumbs = [crumb("Galaxy", "galaxy"), crumb(f"Quadrant {quadrant}")]
    else:
        context["summary_rows"] = _quadrant_summary_rows(sectors)
        title = "Galaxy Map"
        breadcrumbs = [crumb("Galaxy")]

    return render_page(
        "galaxy.html",
        title=title,
        section="galaxy",
        breadcrumbs=breadcrumbs,
        description="An interactive 3D map of every sector placed in this generated galaxy.",
        **context,
    )


def _json_error(message, status):
    response = jsonify({"error": message})
    response.status_code = status
    return response


@bp.route("/galaxy/tiles")
def galaxy_tiles():
    """
    JSON for the map's script: `?tiles=<level/ix/iy/iz,...>`, optional
    `&density=<key>` and `&stamp=<the browser cache's stamp>`. Returns
    `tilecache.fetch_tiles`' payload; a malformed request is a 400, an
    API failure a 502, both as `{"error": ...}` JSON.
    """
    tile_keys = [key for key in (request.args.get("tiles") or "").split(",") if key]
    density_key = request.args.get("density") or None
    known_stamp = request.args.get("stamp") or None
    try:
        payload = fetch_tiles(db_name(), tile_keys, density_key, known_stamp)
    except TileRequestError as exc:
        return _json_error(str(exc), 400)
    except apiclient.NotFoundError as exc:
        return _json_error(str(exc) or "Not found.", 404)
    except apiclient.ApiError as exc:
        log.exception(f"API error while fetching galaxy tiles: {exc}")
        return _json_error("The tiles could not be loaded. Please try again shortly.", 502)
    response = jsonify(payload)
    # Freshness is the tile cache's job (stamps); never let a shared
    # cache hand one visitor's stale answer to another.
    response.headers["Cache-Control"] = "no-store"
    return response
