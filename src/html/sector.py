#!/usr/bin/env python3
# html/sector.py

"""
Sector detail page: the sector's name/size and every system placed in it,
with quadrant and star-type info, linking to `system.py` for each -- plus
an interactive 3D "Sector Map" (see `lib/starmap.py`) of the same systems
plotted by position within the sector, outlined by the sector's real
on-shell wedge shape when it has a galaxy placement (else a plain cube).
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import get_sector
from fmt import esc, linkify_location
from galaxymap import sector_quadrant
from page import query_params, run
from starmap import render_map_panel


def handler():
    params = query_params()
    db_name = params.get("db", "")
    sector_id = params.get("id", "")

    sector = get_sector(db_name, sector_id)
    systems = sector["systems"]
    name_to_id = {row["name"]: row["id"] for row in systems}

    rows = []
    map_systems = []
    for row in systems:
        star_type = row["binary_type"] if row["is_binary"] else (row["stars"][0]["star_type"] if row["stars"] else "")
        rows.append(
            "<tr>"
            f'<td><a href="system.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
            f'<td>{esc(row["quadrant"])}</td>'
            f'<td>{"Yes" if row["is_binary"] else "No"}</td>'
            f'<td>{esc(star_type or "")}</td>'
            f'<td>{linkify_location(db_name, row["location"], name_to_id)}</td>'
            "</tr>"
        )
        if row["position_x_mpc"] is not None and row["stars"]:
            # One entry for a single star, two (primary, then secondary)
            # for a binary -- each with its own star_type/temperature/
            # radius/luminosity, so the map can draw and color each
            # component of a binary independently.
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
                        "star_type": star["star_type"],
                        "temperature_k": star["temperature_k"],
                        "radius_km": star["radius_km"],
                        "luminosity_w": star["luminosity_w"],
                        "temp_display": f'{int(star["temperature_k"])} K',
                    }
                    for star in row["stars"]
                ],
            })
    rows_html = "".join(rows) or '<tr><td colspan="5"><em>None</em></td></tr>'

    center_pc = (
        (sector["center_x_pc"], sector["center_y_pc"], sector["center_z_pc"])
        if sector["placed"]
        else None
    )
    map_html = render_map_panel(
        db_name, sector["edge_mpc"], sector["shell_index"], sector["shell_slot_index"], center_pc, map_systems
    )

    edge_text = f"{sector['edge_ly']:,.2f} ly"

    system_count = len(systems)
    badge_bits = [
        f"Cube edge {esc(edge_text)}",
        f"{system_count} system{'s' if system_count != 1 else ''}",
    ]
    if sector["placed"]:
        quadrant = sector_quadrant(sector["center_x_pc"], sector["center_y_pc"])
        badge_bits.append(
            f'<a href="galaxy.py?db={esc(db_name)}&amp;quadrant={quadrant}">View on Galaxy Map (Quadrant {quadrant})</a>'
        )
    badges_html = "<p class=\"badges\">" + "".join(
        f'<span class="badge">{bit}</span>' for bit in badge_bits
    ) + "</p>"

    body = f"""
<p class="breadcrumb"><a href="index.py">Databases</a> &rarr; <a href="browse.py?db={esc(db_name)}">{esc(db_name)}</a> &rarr; {esc(sector['name'])}</p>
{badges_html}
{map_html}
<section class="panel">
<h2>Systems</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Octant</th><th>Binary</th><th>Star type</th><th>Location</th></tr></thead>
  <tbody>{rows_html}</tbody>
</table></div>
</section>
<script src="static/sectormap.js" defer></script>
"""
    return f"Sector: {sector['name']}", body


run(handler)
