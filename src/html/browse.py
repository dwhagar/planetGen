#!/usr/bin/env python3
# html/browse.py

"""
Per-database overview: every sector (with its system count) and every
standalone system (one generated with no sector, `sector_id IS NULL`).

Both tables are capped at the API's own max page size (`GET /api/sectors`/
`GET /api/systems`, 500 rows -- see docs/api.md's "Pagination") rather
than paginated here: a flat HTML table was always meant for a browsable,
human-scale database, and this project's own galaxy-scale roadmap
(docs/TODO.md Phase 4) is exactly why that cap exists on the API side in
the first place -- `search.py`/`galaxy.py` are the tools for a database
too big for this page's own tables to show in full.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from apiclient import get_sectors, get_systems
from fmt import esc, format_density
from galaxymap import sector_quadrant
from page import query_params, run

_MAX_ROWS = 500


def _galaxy_position_cell(db_name, sector):
    """Links a placed sector's row straight into its Quadrant on
    `galaxy.py`; an unplaced sector (no galaxy position -- everything
    made via `sectorGen.py`'s own standalone CLI, likely most rows in an
    existing database) has no galaxy position to link to at all."""
    if not sector["placed"]:
        return '<span class="hint">Unplaced</span>'
    quadrant = sector_quadrant(sector["center_x_pc"], sector["center_y_pc"])
    return f'<a href="galaxy.py?db={esc(db_name)}&quadrant={quadrant}">Quadrant {quadrant}</a>'


def _truncated_note(total, shown):
    if total <= shown:
        return ""
    return f'<p class="hint">Showing the first {shown} of {total} -- try Search for a more targeted view.</p>'


def handler():
    params = query_params()
    db_name = params.get("db", "")

    sectors_page = get_sectors(db_name, limit=_MAX_ROWS)
    sectors, sectors_total = sectors_page["items"], sectors_page["total"]

    standalone_page = get_systems(db_name, sector_id="none", limit=_MAX_ROWS)
    standalone, standalone_total = standalone_page["items"], standalone_page["total"]

    sector_rows = "".join(
        "<tr>"
        f'<td><a href="sector.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
        f'<td>{row["system_count"]}</td>'
        f'<td>{format_density(row["edge_ly"], row["system_count"])}</td>'
        f'<td>{_galaxy_position_cell(db_name, row)}</td>'
        "</tr>"
        for row in sectors
    ) or '<tr><td colspan="4"><em>None</em></td></tr>'

    standalone_rows = "".join(
        "<tr>"
        f'<td><a href="system.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
        f'<td>{_esc_bool(row["is_binary"])}</td>'
        f'<td>{esc(row["star_summary"])}</td>'
        "</tr>"
        for row in standalone
    ) or '<tr><td colspan="3"><em>None</em></td></tr>'

    badges_html = "<p class=\"badges\">" + "".join(
        f'<span class="badge">{bit}</span>' for bit in (
            f"{sectors_total} sector{'s' if sectors_total != 1 else ''}",
            f"{standalone_total} standalone system{'s' if standalone_total != 1 else ''}",
        )
    ) + "</p>"

    body = f"""
<p class="breadcrumb"><a href="index.py">Databases</a> &rarr; {esc(db_name)} &middot; <a href="search.py?db={esc(db_name)}">Search</a></p>
{badges_html}
<section class="panel" id="sectors">
<div class="panel-header">
  <h2>Sectors</h2>
  <a href="galaxy.py?db={esc(db_name)}">Galaxy Map &rarr;</a>
</div>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Systems</th><th>Density</th><th>Galaxy Position</th></tr></thead>
  <tbody>{sector_rows}</tbody>
</table></div>
{_truncated_note(sectors_total, len(sectors))}
</section>

<section class="panel" id="standalone-systems">
<h2>Standalone Systems</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Binary</th><th>Star type</th></tr></thead>
  <tbody>{standalone_rows}</tbody>
</table></div>
{_truncated_note(standalone_total, len(standalone))}
</section>
"""
    return f"Browse: {db_name}", body


def _esc_bool(value):
    return "Yes" if value else "No"


run(handler)
