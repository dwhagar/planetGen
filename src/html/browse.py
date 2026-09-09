#!/usr/bin/env python3
# html/browse.py

"""
Per-database overview: every sector (with its system count) and every
standalone system (one generated with no sector, `sector_id IS NULL`).
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from dbutil import esc, fetch_all, format_density, open_readonly, resolve_db_path
from galaxymap import sector_quadrant
from page import query_params, run


def _star_summary(conn, star_system_id, is_binary):
    """One-line star description for a system row: the single/primary
    star's descriptive type, or a binary's combined type."""
    if is_binary:
        row = conn.execute(
            "SELECT binary_type FROM star_systems WHERE id = ?", (star_system_id,)
        ).fetchone()
        return row["binary_type"] if row else ""
    row = conn.execute(
        "SELECT star_type FROM stars WHERE star_system_id = ? AND role = 'single'",
        (star_system_id,),
    ).fetchone()
    return row["star_type"] if row else ""


def _galaxy_position_cell(db_name, row):
    """Links a placed sector's row straight into its Quadrant on
    `galaxy.py`; an unplaced sector (`center_x_pc IS NULL` -- everything
    made via `sectorGen.py`'s own standalone CLI, likely most rows in an
    existing database) has no galaxy position to link to at all."""
    if row["center_x_pc"] is None:
        return '<span class="hint">Unplaced</span>'
    quadrant = sector_quadrant(row["center_x_pc"], row["center_y_pc"])
    return f'<a href="galaxy.py?db={esc(db_name)}&quadrant={quadrant}">Quadrant {quadrant}</a>'


def handler():
    params = query_params()
    db_name = params.get("db", "")
    path = resolve_db_path(db_name)
    conn = open_readonly(path)
    try:
        sectors = fetch_all(
            conn,
            """
            SELECT s.id, s.name, s.edge_mpc, s.center_x_pc, s.center_y_pc,
                   COUNT(ss.id) AS system_count
            FROM sectors s
            LEFT JOIN star_systems ss ON ss.sector_id = s.id
            GROUP BY s.id
            ORDER BY s.name
            """,
        )
        standalone = fetch_all(
            conn,
            """
            SELECT id, name, is_binary
            FROM star_systems
            WHERE sector_id IS NULL
            ORDER BY name
            """,
        )

        sector_rows = "".join(
            "<tr>"
            f'<td><a href="sector.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
            f'<td>{row["system_count"]}</td>'
            f'<td>{format_density(row["edge_mpc"], row["system_count"])}</td>'
            f'<td>{_galaxy_position_cell(db_name, row)}</td>'
            "</tr>"
            for row in sectors
        ) or '<tr><td colspan="4"><em>None</em></td></tr>'

        standalone_rows = "".join(
            "<tr>"
            f'<td><a href="system.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
            f'<td>{_esc_bool(row["is_binary"])}</td>'
            f'<td>{esc(_star_summary(conn, row["id"], row["is_binary"]))}</td>'
            "</tr>"
            for row in standalone
        ) or '<tr><td colspan="3"><em>None</em></td></tr>'
    finally:
        conn.close()

    badges_html = "<p class=\"badges\">" + "".join(
        f'<span class="badge">{bit}</span>' for bit in (
            f"{len(sectors)} sector{'s' if len(sectors) != 1 else ''}",
            f"{len(standalone)} standalone system{'s' if len(standalone) != 1 else ''}",
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
</section>

<section class="panel" id="standalone-systems">
<h2>Standalone Systems</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Binary</th><th>Star type</th></tr></thead>
  <tbody>{standalone_rows}</tbody>
</table></div>
</section>
"""
    return f"Browse: {db_name}", body


def _esc_bool(value):
    return "Yes" if value else "No"


run(handler)
