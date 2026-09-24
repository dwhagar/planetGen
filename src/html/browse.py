#!/usr/bin/env python3
# html/browse.py

"""
Per-database overview: every sector (with its system count) and every
standalone system (one generated with no sector, `sector_id IS NULL`).

Both tables are genuinely paginated (`GET /api/sectors`/`GET /api/systems`'s
own `limit`/`offset` -- see docs/api.md's "Pagination"), `_PAGE_SIZE` rows
at a time with Prev/Next controls, rather than the single-page "first 500,
then go use Search instead" cap this page used to have -- that cap meant a
database with more sectors/standalone systems than fit on one page had no
way to ever reach the rest of them from this listing at all, only via
Search's own targeted filters. Each table's own `sector_offset`/
`standalone_offset` query/POST parameter is independent, so paginating one
table never resets the other's current page.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from apiclient import get_sectors, get_systems
from fmt import esc, format_density, post_link
from galaxymap import sector_quadrant
from page import nav_params, run

_PAGE_SIZE = 100
"""int: Rows per page for both tables -- matches the API's own default
`limit` (see docs/api.md's "Pagination"), well under its 500-row max."""


def _galaxy_position_cell(db_name, sector):
    """Links a placed sector's row straight into its Quadrant on
    `galaxy.py`; an unplaced sector (no galaxy position -- everything
    made via `sectorGen.py`'s own standalone CLI, likely most rows in an
    existing database) has no galaxy position to link to at all."""
    if not sector["placed"]:
        return '<span class="hint">Unplaced</span>'
    quadrant = sector_quadrant(sector["center_x_pc"], sector["center_y_pc"])
    return post_link("galaxy.py", {"db": db_name, "quadrant": quadrant}, f"Quadrant {quadrant}")


def _parse_offset(raw):
    """Parses a `sector_offset`/`standalone_offset` nav param back into a
    non-negative int -- defensively: a hand-edited/stale URL carrying a
    negative or non-numeric value falls back to 0 (page 1) rather than
    raising, since this only ever controls which page is shown, never
    anything security-sensitive."""
    try:
        value = int(raw)
    except (TypeError, ValueError):
        return 0
    return max(0, value)


def _pagination_controls(db_name, offset_param, total, offset, other_offset_param, other_offset):
    """
    Builds a "showing X-Y of Z" line plus Prev/Next controls for one
    table -- `post_link`s carrying the OTHER table's own current offset
    unchanged (`other_offset_param`/`other_offset`), so paginating this
    table never resets the other one back to page 1.

    Args:
        db_name (str): The current `?db=` value.
        offset_param (str): This table's own offset param name
            (`"sector_offset"` or `"standalone_offset"`).
        total (int): Total row count (`get_sectors`/`get_systems`'s own
            `total`).
        offset (int): This table's own current offset.
        other_offset_param (str): The other table's offset param name.
        other_offset (int): The other table's own current offset --
            carried through unchanged on every link this builds.

    Returns:
        str: A `<div class="pagination">` block, or `""` if everything
            fits on one page already (nothing to paginate).
    """
    if total <= _PAGE_SIZE:
        return ""

    shown_from = offset + 1
    shown_to = min(offset + _PAGE_SIZE, total)
    summary = f"{shown_from:,}&ndash;{shown_to:,} of {total:,}"

    prev_offset = max(0, offset - _PAGE_SIZE)
    next_offset = offset + _PAGE_SIZE
    prev_params = {"db": db_name, offset_param: prev_offset, other_offset_param: other_offset}
    next_params = {"db": db_name, offset_param: next_offset, other_offset_param: other_offset}
    prev_html = (
        post_link("browse.py", prev_params, "&larr; Prev") if offset > 0 else '<span class="hint">&larr; Prev</span>'
    )
    next_html = (
        post_link("browse.py", next_params, "Next &rarr;")
        if next_offset < total
        else '<span class="hint">Next &rarr;</span>'
    )
    return f'<div class="pagination">{prev_html}<span class="hint">{summary}</span>{next_html}</div>'


def handler(db_name=None):
    params = nav_params()
    if db_name is None:
        db_name = params.get("db", "")
    sector_offset = _parse_offset(params.get("sector_offset"))
    standalone_offset = _parse_offset(params.get("standalone_offset"))

    sectors_page = get_sectors(db_name, limit=_PAGE_SIZE, offset=sector_offset)
    sectors, sectors_total = sectors_page["items"], sectors_page["total"]

    standalone_page = get_systems(db_name, sector_id="none", limit=_PAGE_SIZE, offset=standalone_offset)
    standalone, standalone_total = standalone_page["items"], standalone_page["total"]

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
{_pagination_controls(db_name, "sector_offset", sectors_total, sector_offset, "standalone_offset", standalone_offset)}
</section>

<section class="panel" id="standalone-systems">
<h2>Standalone Systems</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Binary</th><th>Star type</th></tr></thead>
  <tbody>{standalone_rows}</tbody>
</table></div>
{_pagination_controls(db_name, "standalone_offset", standalone_total, standalone_offset, "sector_offset", sector_offset)}
</section>
"""
    return f"Browse: {db_name}", body


def _esc_bool(value):
    return "Yes" if value else "No"


if __name__ == "__main__":
    run(handler)
