#!/usr/bin/env python3
# html/galaxy.py

"""
Galaxy Map page: a real perspective-camera WebGL scene
(`lib/galaxymap3d.py` + `static/galaxymap3d.js`, three.js) a visitor can
rotate, dolly, and click through -- replaced the previous flat, face-on
SVG projection (`lib/galaxymap.py`'s own former `render_galaxy_map_panel`,
now removed) entirely, rather than living alongside it as a separate
page: a flat SVG's fixed-radius markers necessarily grow relative to the
view as you zoom in (there's no camera to shrink them the way a real
perspective projection does for free), which read as "the star icon gets
bigger and bigger, hiding sectors" and made scroll/+/-/click zoom feel
broken well before it actually was.

This page itself only fetches what the map's *first* paint needs
(`get_galaxy_shape`, then the zoomed-all-the-way-out starting view's
tiles through `lib/tilecache.py`'s disk cache) -- every later view, as
the visitor's camera moves, is fetched directly by the page's own
client-side JS from `galaxy_tiles.py`, never through this handler again.

Below the map, this page still keeps its own two data tables --
independent of how the map is drawn, and still useful as a plain-text
overview: every sector actually placed in the galaxy (`GET
/api/galaxy/sectors`), grouped into four azimuthal Quadrants (I-IV) and
concentric Rings (fixed-width `shell_index` bands, `lib/galaxymap.py`).
`?quadrant=I|II|III|IV` swaps the page's own table from a per-Quadrant
summary (the full-galaxy default, since a flat list of every placed
sector at once doesn't scale) to a full sector list for that one
Quadrant, sorted by distance from the core -- it no longer crops the map
itself (the 3D map has its own free-flying camera, not a Quadrant-cropped
viewBox).

Sectors never placed in the galaxy (made via `sectorGen.py`'s own
standalone CLI -- most of what exists in a database today) have no
position to plot here at all; they stay in `browse.py`'s own flat sector
table unchanged, which now also links each *placed* sector's row into
this page (see `browse.py`).
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import get_galaxy_sectors, get_galaxy_shape
from fmt import esc, post_link, static_url
from galaxymap import QUADRANT_LABELS, ring_bounds_ly, sector_quadrant, sector_ring
from galaxymap3d import initial_tile_request, render_galaxy_map3d_panel, view_radius_bounds
from page import nav_params, run
from pagination import page_slice, parse_page, render_pagination
from tilecache import fetch_tiles

try:
    from stellarObjects.program_constants import DEFAULT_SECTOR_EDGE_LY
    from stellarObjects.utils import ly_to_pc, pc_to_ly
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # duplicated fallback, matching every other page's identical pattern.
    DEFAULT_SECTOR_EDGE_LY = 11.5
    pc_to_ly = None
    def ly_to_pc(ly):
        return ly / 3.2616


def _display_ly(galactic_radius_pc):
    if pc_to_ly:
        return pc_to_ly(galactic_radius_pc)
    return galactic_radius_pc * 3.2616


def _quadrant_summary_table(db_name, sectors):
    """Full-galaxy default view: one row per Quadrant (sector/system
    counts, Ring span) rather than every placed sector at once -- the
    per-Quadrant drill-down (`?quadrant=`) is where an actual sector list
    shows up (`_quadrant_sector_table`)."""
    by_quadrant = {label: [] for label in QUADRANT_LABELS}
    for sector in sectors:
        by_quadrant[sector_quadrant(sector["x"], sector["y"])].append(sector)

    rows = []
    for label in QUADRANT_LABELS:
        members = by_quadrant[label]
        system_total = sum(s["system_count"] or 0 for s in members)
        if members:
            max_ring = max(sector_ring(s["shell_index"]) for s in members)
            _inner, outer_ly = ring_bounds_ly(max_ring)
            span = f"out to ~{outer_ly:,.0f} ly"
        else:
            span = "&ndash;"
        rows.append(
            "<tr>"
            f'<td>{post_link("galaxy.py", {"db": db_name, "quadrant": label}, f"Quadrant {label}")}</td>'
            f"<td>{len(members)}</td>"
            f"<td>{system_total}</td>"
            f"<td>{span}</td>"
            "</tr>"
        )
    return "".join(rows)


def _quadrant_sector_table(db_name, sectors, quadrant, page):
    """One page (`lib/pagination.py`) of the Quadrant's placed sectors,
    nearest the core first. Returns `(rows_html, pager_html)`."""
    members = [s for s in sectors if sector_quadrant(s["x"], s["y"]) == quadrant]
    members.sort(key=lambda s: s["galactic_radius_pc"])
    page_members, page = page_slice(members, page)
    pager_html = render_pagination(
        "galaxy.py", {"db": db_name, "quadrant": quadrant}, "page", page, len(members),
        anchor="galaxy-table", label="Sector pages",
    )

    rows = []
    for sector in page_members:
        ring = sector_ring(sector["shell_index"])
        distance_ly = _display_ly(sector["galactic_radius_pc"])
        rows.append(
            "<tr>"
            f'<td>{post_link("sector.py", {"db": db_name, "id": sector["id"]}, esc(sector["name"]))}</td>'
            f"<td>Ring {ring}</td>"
            f"<td>{distance_ly:,.1f} ly</td>"
            f'<td>{sector["system_count"] or 0}</td>'
            "</tr>"
        )
    rows_html = "".join(rows) or '<tr><td colspan="4"><em>No sectors placed in this Quadrant yet.</em></td></tr>'
    return rows_html, pager_html


def handler():
    params = nav_params()
    db_name = params.get("db", "")
    quadrant = (params.get("quadrant") or "").upper() or None
    if quadrant not in QUADRANT_LABELS:
        quadrant = None

    sectors = get_galaxy_sectors(db_name)
    galaxy_shape = get_galaxy_shape(db_name)
    edge_pc = galaxy_shape["edge_pc"] if galaxy_shape else ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    _min_radius, max_radius = view_radius_bounds(edge_pc, galaxy_shape)
    tile_keys, density_key = initial_tile_request(max_radius, galaxy_shape is not None)
    initial_view = fetch_tiles(db_name, tile_keys, density_key)

    map_html = render_galaxy_map3d_panel(db_name, galaxy_shape, edge_pc, initial_view)

    if quadrant:
        table_title = f"Sectors in Quadrant {quadrant}"
        table_head = "<tr><th>Name</th><th>Ring</th><th>Distance from core</th><th>Systems</th></tr>"
        table_rows, pager_html = _quadrant_sector_table(db_name, sectors, quadrant, parse_page(params.get("page")))
    else:
        table_title = "Quadrants"
        table_head = "<tr><th>Quadrant</th><th>Placed sectors</th><th>Total systems</th><th>Extent</th></tr>"
        table_rows = _quadrant_summary_table(db_name, sectors)
        pager_html = ""

    placed_count = len(sectors)
    badges_html = (
        f'<p class="badges"><span class="badge">{placed_count} placed sector'
        f'{"s" if placed_count != 1 else ""}</span></p>'
    )

    title = f"Galaxy Map: Quadrant {quadrant}" if quadrant else "Galaxy Map"
    breadcrumb_html = f'<p class="breadcrumb">{post_link("browse.py", {"db": db_name}, esc(db_name))} &rarr; {esc(title)}</p>'
    body = f"""
<div class="page-subhead">{breadcrumb_html}{badges_html}</div>
{map_html}
<section class="panel" id="galaxy-table">
<h2>{esc(table_title)}</h2>
<div class="table-scroll"><table>
  <thead>{table_head}</thead>
  <tbody>{table_rows}</tbody>
</table></div>
{pager_html}
</section>
<script type="module" src="{static_url("galaxymap3d.js")}"></script>
"""
    return title, body


run(handler)
