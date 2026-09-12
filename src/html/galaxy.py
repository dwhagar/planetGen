#!/usr/bin/env python3
# html/galaxy.py

"""
Galaxy Map page: every sector actually placed in the galaxy (a non-NULL
`sectors.center_x/y/z_pc` -- see `docs/design/galaxy-coordinate-system.md`
and `galaxyGen.py`) plotted by its real position, grouped into four
azimuthal Quadrants (I-IV) and concentric Rings (fixed-width `shell_index`
bands) -- see `lib/galaxymap.py` for the map geometry/rendering this page
just supplies data and drill-down tables around.

`?quadrant=I|II|III|IV` zooms the map into that Quadrant and swaps the
page's own table from a per-Quadrant summary (the full-galaxy default,
since a flat list of every placed sector at once doesn't scale) to a full
sector list for that one Quadrant, sorted by distance from the core.

Sectors never placed in the galaxy (made via `sectorGen.py`'s own
standalone CLI -- most of what exists in a database today) have no
position to plot here at all; they stay in `browse.py`'s own flat sector
table unchanged, which now also links each *placed* sector's row into this
page (see `browse.py`).

Every galaxy-placed nebula/asteroid field (`phenomenonGen.py --sector-id`,
schema.sql's "v18" header note) is also plotted here, as a small fixed-size
dot -- unlike a sector, a phenomenon has no real "how much is here"
quantity to size a dot by, and at this scale its own physical extent
(which can itself span several sectors) would be a misleading dot size;
that real extent is instead depicted where it belongs, as a translucent
cloud on the Sector Map (`sector.py`) of any sector it reaches into.
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import get_galaxy_phenomena, get_galaxy_sectors
from fmt import esc
from galaxymap import QUADRANT_LABELS, render_galaxy_map_panel, ring_bounds_ly, sector_quadrant, sector_ring
from page import query_params, run

try:
    from stellarObjects.utils import pc_to_ly
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # fall back to showing the raw stored parsec unit rather than failing.
    pc_to_ly = None


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
            f'<td><a href="galaxy.py?db={esc(db_name)}&quadrant={label}">Quadrant {label}</a></td>'
            f"<td>{len(members)}</td>"
            f"<td>{system_total}</td>"
            f"<td>{span}</td>"
            "</tr>"
        )
    return "".join(rows)


def _quadrant_sector_table(db_name, sectors, quadrant):
    members = [s for s in sectors if sector_quadrant(s["x"], s["y"]) == quadrant]
    members.sort(key=lambda s: s["galactic_radius_pc"])

    rows = []
    for sector in members:
        ring = sector_ring(sector["shell_index"])
        distance_ly = _display_ly(sector["galactic_radius_pc"])
        rows.append(
            "<tr>"
            f'<td><a href="sector.py?db={esc(db_name)}&id={sector["id"]}">{esc(sector["name"])}</a></td>'
            f"<td>Ring {ring}</td>"
            f"<td>{distance_ly:,.1f} ly</td>"
            f'<td>{sector["system_count"] or 0}</td>'
            "</tr>"
        )
    return "".join(rows) or '<tr><td colspan="4"><em>No sectors placed in this Quadrant yet.</em></td></tr>'


def handler():
    params = query_params()
    db_name = params.get("db", "")
    quadrant = (params.get("quadrant") or "").upper() or None
    if quadrant not in QUADRANT_LABELS:
        quadrant = None

    sectors = get_galaxy_sectors(db_name)
    phenomena = get_galaxy_phenomena(db_name)

    map_html = render_galaxy_map_panel(db_name, sectors, quadrant=quadrant, phenomena=phenomena)

    if quadrant:
        table_title = f"Sectors in Quadrant {quadrant}"
        table_head = "<tr><th>Name</th><th>Ring</th><th>Distance from core</th><th>Systems</th></tr>"
        table_rows = _quadrant_sector_table(db_name, sectors, quadrant)
    else:
        table_title = "Quadrants"
        table_head = "<tr><th>Quadrant</th><th>Placed sectors</th><th>Total systems</th><th>Extent</th></tr>"
        table_rows = _quadrant_summary_table(db_name, sectors)

    placed_count = len(sectors)
    phenomenon_count = len(phenomena)
    badge_bits = [f"{placed_count} placed sector{'s' if placed_count != 1 else ''}"]
    if phenomenon_count:
        badge_bits.append(f"{phenomenon_count} placed nebula/asteroid field{'s' if phenomenon_count != 1 else ''}")
    badges_html = "<p class=\"badges\">" + "".join(
        f'<span class="badge">{bit}</span>' for bit in badge_bits
    ) + "</p>"

    title = f"Galaxy Map: Quadrant {quadrant}" if quadrant else "Galaxy Map"
    body = f"""
<p class="breadcrumb"><a href="index.py">Databases</a> &rarr; <a href="browse.py?db={esc(db_name)}">{esc(db_name)}</a> &rarr; {esc(title)}</p>
{badges_html}
{map_html}
<section class="panel">
<h2>{esc(table_title)}</h2>
<div class="table-scroll"><table>
  <thead>{table_head}</thead>
  <tbody>{table_rows}</tbody>
</table></div>
</section>
"""
    return title, body


run(handler)
