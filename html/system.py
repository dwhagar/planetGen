#!/usr/bin/env python3
# html/system.py

"""
System detail page: stars, planets/moons, asteroid belts, and the full
rendered wiki page (wikitext or Markdown, toggled via `?format=`) stored
for this system.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from dbutil import NotFoundError, esc, fetch_all, fetch_one, open_readonly, resolve_db_path
from page import query_params, run


def _stars_html(conn, system_id):
    stars = fetch_all(
        conn,
        "SELECT * FROM stars WHERE star_system_id = ? ORDER BY CASE role WHEN 'primary' THEN 0 WHEN 'single' THEN 0 ELSE 1 END",
        (system_id,),
    )
    rows = "".join(
        "<tr>"
        f'<td>{esc(row["role"])}</td>'
        f'<td>{esc(row["name"])}</td>'
        f'<td>{esc(row["star_type"])}</td>'
        f'<td>{esc(row["table_mass"])}</td>'
        f'<td>{esc(row["table_radius"])}</td>'
        f'<td>{esc(row["table_temp"])}</td>'
        f'<td>{esc(row["table_lum"])}</td>'
        "</tr>"
        for row in stars
    )
    if not rows:
        return ""
    return f"""
<h2>Stars</h2>
<table>
  <thead><tr><th>Role</th><th>Name</th><th>Type</th><th>Mass</th><th>Radius</th><th>Temp</th><th>Luminosity</th></tr></thead>
  <tbody>{rows}</tbody>
</table>
"""


def _planet_row(conn, planet, indent=""):
    rows = [
        "<tr>"
        f'<td>{indent}{esc(planet["name"])}</td>'
        f'<td>{esc(planet["planet_class"])}</td>'
        f'<td>{"Gas Giant" if planet["body_type"] == "g" else "Terrestrial"}</td>'
        f'<td>{esc(planet["zone"])}</td>'
        f'<td>{esc(planet["table_distance"])}</td>'
        f'<td>{esc(planet["table_period"])}</td>'
        f'<td>{esc(planet["table_gravity"])}</td>'
        "</tr>"
    ]
    moons = fetch_all(
        conn,
        "SELECT * FROM planets WHERE parent_planet_id = ? ORDER BY orbital_index",
        (planet["id"],),
    )
    for moon in moons:
        rows.extend(_planet_row(conn, moon, indent="&nbsp;&nbsp;&nbsp;&nbsp;└ "))
    return rows


def _bodies_html(conn, system_id):
    planets = fetch_all(
        conn,
        "SELECT * FROM planets WHERE star_system_id = ? AND parent_planet_id IS NULL ORDER BY orbital_index",
        (system_id,),
    )
    belts = fetch_all(
        conn,
        "SELECT * FROM asteroid_belts WHERE star_system_id = ? ORDER BY orbital_index",
        (system_id,),
    )

    planet_rows = []
    for planet in planets:
        planet_rows.extend(_planet_row(conn, planet))
    planet_html = ""
    if planet_rows:
        planet_html = f"""
<h2>Planets &amp; Moons</h2>
<table>
  <thead><tr><th>Name</th><th>Class</th><th>Type</th><th>Zone</th><th>Distance</th><th>Period</th><th>Gravity</th></tr></thead>
  <tbody>{''.join(planet_rows)}</tbody>
</table>
"""

    belt_html = ""
    if belts:
        belt_rows = "".join(
            "<tr>"
            f'<td>{esc(row["density"])}</td>'
            f'<td>{row["distance_km"]:,.0f} km</td>'
            f'<td>{esc(row["composition_summary"])}</td>'
            "</tr>"
            for row in belts
        )
        belt_html = f"""
<h2>Asteroid Belts</h2>
<table>
  <thead><tr><th>Density</th><th>Distance</th><th>Composition</th></tr></thead>
  <tbody>{belt_rows}</tbody>
</table>
"""

    return planet_html + belt_html


def handler():
    params = query_params()
    db_name = params.get("db", "")
    system_id = params.get("id", "")
    fmt = params.get("format", "wikitext")
    if fmt not in ("wikitext", "markdown"):
        fmt = "wikitext"

    path = resolve_db_path(db_name)
    conn = open_readonly(path)
    try:
        system = fetch_one(conn, "SELECT * FROM star_systems WHERE id = ?", (system_id,))
        if system is None:
            raise NotFoundError(f"No such system: {system_id!r}")

        back_html = f'<p><a href="index.py">Databases</a>'
        if system["sector_id"] is not None:
            back_html += f' &rarr; <a href="sector.py?db={esc(db_name)}&id={system["sector_id"]}">Sector</a>'
        back_html += f" &rarr; {esc(system['name'])}</p>"

        summary_bits = []
        if system["quadrant"]:
            summary_bits.append(f"Quadrant {esc(system['quadrant'])}")
        summary_bits.append("Binary system" if system["is_binary"] else "Single star")
        summary_html = f"<p>{' &middot; '.join(summary_bits)}</p>"

        content = system["wikitext_content"] if fmt == "wikitext" else system["markdown_content"]
        other_fmt = "markdown" if fmt == "wikitext" else "wikitext"
        toggle_html = (
            f'<p><a href="system.py?db={esc(db_name)}&id={system_id}&format={other_fmt}">'
            f'Switch to {other_fmt}</a></p>'
        )

        stars_html = _stars_html(conn, system_id)
        bodies_html = _bodies_html(conn, system_id)
    finally:
        conn.close()

    rows = (content or "").count("\n") + 3
    body = f"""
{back_html}
{summary_html}
{stars_html}
{bodies_html}
<h2>Rendered Page ({esc(fmt)})</h2>
{toggle_html}
<textarea readonly rows="{min(rows, 40)}" class="content-box">{esc(content)}</textarea>
"""
    return f"System: {system['name']}", body


run(handler)
