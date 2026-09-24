#!/usr/bin/env python3
# html/browse.py

"""
Per-database overview: every sector (with its system count) and every
standalone system (one generated with no sector, `sector_id IS NULL`).

Both tables are paginated through the API (`GET /api/sectors`/`GET
/api/systems`'s own `limit`/`offset` -- see docs/api.md's "Pagination"),
`pagination.PAGE_SIZE` rows at a time with the site's shared pager
(`lib/pagination.py`). Each table's own `sectors_page`/`standalone_page`
parameter is independent, so paging one table never resets the other's
current page.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from apiclient import get_sectors, get_systems
from fmt import esc, format_density, post_link
from galaxymap import sector_quadrant
from page import nav_params, run
from pagination import fetch_page, parse_page, render_pagination


def _galaxy_position_cell(db_name, sector):
    """Links a placed sector's row straight into its Quadrant on
    `galaxy.py`; an unplaced sector (no galaxy position -- everything
    made via `sectorGen.py`'s own standalone CLI, likely most rows in an
    existing database) has no galaxy position to link to at all."""
    if not sector["placed"]:
        return '<span class="hint">Unplaced</span>'
    quadrant = sector_quadrant(sector["center_x_pc"], sector["center_y_pc"])
    return post_link("galaxy.py", {"db": db_name, "quadrant": quadrant}, f"Quadrant {quadrant}")


def handler(db_name=None):
    params = nav_params()
    if db_name is None:
        db_name = params.get("db", "")
    sectors_page, sectors_page_no = fetch_page(
        lambda limit, offset: get_sectors(db_name, limit=limit, offset=offset),
        parse_page(params.get("sectors_page")),
    )
    sectors, sectors_total = sectors_page["items"], sectors_page["total"]

    standalone_page, standalone_page_no = fetch_page(
        lambda limit, offset: get_systems(db_name, sector_id="none", limit=limit, offset=offset),
        parse_page(params.get("standalone_page")),
    )
    standalone, standalone_total = standalone_page["items"], standalone_page["total"]

    page_state = {"db": db_name, "sectors_page": sectors_page_no, "standalone_page": standalone_page_no}
    sectors_pager = render_pagination(
        "browse.py", page_state, "sectors_page", sectors_page_no, sectors_total,
        anchor="sectors", label="Sector pages",
    )
    standalone_pager = render_pagination(
        "browse.py", page_state, "standalone_page", standalone_page_no, standalone_total,
        anchor="standalone-systems", label="Standalone system pages",
    )

    sector_rows = "".join(
        "<tr>"
        f'<td>{post_link("sector.py", {"db": db_name, "id": row["id"]}, esc(row["name"]))}</td>'
        f'<td>{row["system_count"]}</td>'
        f'<td>{format_density(row["edge_ly"], row["system_count"])}</td>'
        f'<td>{_galaxy_position_cell(db_name, row)}</td>'
        "</tr>"
        for row in sectors
    ) or '<tr><td colspan="4"><em>None</em></td></tr>'

    standalone_rows = "".join(
        "<tr>"
        f'<td>{post_link("system.py", {"db": db_name, "id": row["id"]}, esc(row["name"]))}</td>'
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
<p class="breadcrumb">{esc(db_name)} &middot; {post_link("search.py", {"db": db_name}, "Search")}</p>
{badges_html}
<section class="panel" id="sectors">
<div class="panel-header">
  <h2>Sectors</h2>
  {post_link("galaxy.py", {"db": db_name}, "Galaxy Map &rarr;")}
</div>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Systems</th><th>Density</th><th>Galaxy Position</th></tr></thead>
  <tbody>{sector_rows}</tbody>
</table></div>
{sectors_pager}
</section>

<section class="panel" id="standalone-systems">
<h2>Standalone Systems</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Binary</th><th>Star type</th></tr></thead>
  <tbody>{standalone_rows}</tbody>
</table></div>
{standalone_pager}
</section>
"""
    return f"Browse: {db_name}", body


def _esc_bool(value):
    return "Yes" if value else "No"


if __name__ == "__main__":
    run(handler)
