#!/usr/bin/env python3
# html/sector.py

"""
Sector detail page: the sector's name/size and every system placed in it,
with quadrant and star-type info, linking to `system.py` for each -- plus
an isometric "Sector Map" (see `lib/starmap.py`) of the same systems
plotted by position within the sector's cube.
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))
# Falls back to src/ (stellarObjects now lives at src/stellarObjects/, src
# layout) so `stellarObjects` is importable even when it hasn't been
# `pip install`-ed system-wide -- true for the default deployment layout
# (`html/` and `src/` as siblings under /var/lib/planetGen).
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(_HTML_DIR)), "src"))

from dbutil import NotFoundError, esc, fetch_all, fetch_one, linkify_location, open_readonly, resolve_db_path
from page import query_params, run
from starmap import render_map_panel

try:
    from stellarObjects.utils import milliparsecs_to_ly
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # fall back to showing the raw stored unit rather than failing outright.
    milliparsecs_to_ly = None


def handler():
    params = query_params()
    db_name = params.get("db", "")
    sector_id = params.get("id", "")
    path = resolve_db_path(db_name)
    conn = open_readonly(path)
    try:
        sector = fetch_one(conn, "SELECT * FROM sectors WHERE id = ?", (sector_id,))
        if sector is None:
            raise NotFoundError(f"No such sector: {sector_id!r}")

        systems = fetch_all(
            conn,
            """
            SELECT id, name, is_binary, quadrant, location,
                   position_x_mpc, position_y_mpc, position_z_mpc, binary_type
            FROM star_systems
            WHERE sector_id = ?
            ORDER BY name
            """,
            (sector_id,),
        )

        name_to_id = {row["name"]: row["id"] for row in systems}

        rows = []
        map_systems = []
        for row in systems:
            if row["is_binary"]:
                star_type = row["binary_type"] or ""
            else:
                star_row = conn.execute(
                    "SELECT star_type FROM stars WHERE star_system_id = ? AND role = 'single'",
                    (row["id"],),
                ).fetchone()
                star_type = star_row["star_type"] if star_row else ""
            rows.append(
                "<tr>"
                f'<td><a href="system.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
                f'<td>{esc(row["quadrant"])}</td>'
                f'<td>{"Yes" if row["is_binary"] else "No"}</td>'
                f'<td>{esc(star_type)}</td>'
                f'<td>{linkify_location(db_name, row["location"], name_to_id)}</td>'
                "</tr>"
            )
            if row["position_x_mpc"] is not None:
                # One row for a single star, two (primary then secondary) for
                # a binary -- each with its own star_type/temperature/
                # radius/luminosity, unlike `binary_type`/`binary_temperature_k`
                # on `star_systems` (which have no per-component color info,
                # only a combined summary), so the map can draw and color
                # each component of a binary independently.
                star_rows = fetch_all(
                    conn,
                    "SELECT role, star_type, temperature_k, radius_km, luminosity_w, table_temp"
                    " FROM stars WHERE star_system_id = ?"
                    " ORDER BY CASE role WHEN 'secondary' THEN 1 ELSE 0 END",
                    (row["id"],),
                )
                if not star_rows:
                    continue
                map_systems.append({
                    "id": row["id"],
                    "name": row["name"],
                    "quadrant": row["quadrant"],
                    "location": row["location"],
                    "x": row["position_x_mpc"],
                    "y": row["position_y_mpc"],
                    "z": row["position_z_mpc"],
                    "stars": [
                        {
                            "star_type": star_row["star_type"],
                            "temperature_k": star_row["temperature_k"],
                            "radius_km": star_row["radius_km"],
                            "luminosity_w": star_row["luminosity_w"],
                            "temp_display": star_row["table_temp"],
                        }
                        for star_row in star_rows
                    ],
                })
        rows_html = "".join(rows) or '<tr><td colspan="5"><em>None</em></td></tr>'
        map_html = render_map_panel(db_name, sector["edge_mpc"], map_systems)

        edge_ly = milliparsecs_to_ly(sector["edge_mpc"]) if milliparsecs_to_ly else None
        edge_text = f"{edge_ly:,.2f} ly" if edge_ly is not None else f"{sector['edge_mpc']:,.2f} mpc"
    finally:
        conn.close()

    system_count = len(systems)
    badges_html = "<p class=\"badges\">" + "".join(
        f'<span class="badge">{bit}</span>' for bit in (
            f"Cube edge {esc(edge_text)}",
            f"{system_count} system{'s' if system_count != 1 else ''}",
        )
    ) + "</p>"

    body = f"""
<p class="breadcrumb"><a href="index.py">Databases</a> &rarr; <a href="browse.py?db={esc(db_name)}">{esc(db_name)}</a> &rarr; {esc(sector['name'])}</p>
{badges_html}
{map_html}
<section class="panel">
<h2>Systems</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Quadrant</th><th>Binary</th><th>Star type</th><th>Location</th></tr></thead>
  <tbody>{rows_html}</tbody>
</table></div>
</section>
<script src="static/sectormap.js" defer></script>
"""
    return f"Sector: {sector['name']}", body


run(handler)
