#!/usr/bin/env python3
# html/system.py

"""
System detail page: stars, planets/moons, asteroid belts, and the
system's full description. Defaults to rendering `markdown_content` as
actual HTML (via `mdconvert.markdown_to_html`) so the description reads
like a normal page instead of a wall of raw Markdown; `?view=source`
switches to the original raw-text view (wikitext or Markdown, toggled via
`&format=`), which is what you want when copy-pasting into a wiki.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from dbutil import NotFoundError, esc, fetch_all, fetch_one, open_readonly, resolve_db_path
from mdconvert import markdown_to_html
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
<section class="panel">
<h2>Stars</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Role</th><th>Name</th><th>Type</th><th>Mass</th><th>Radius</th><th>Temp</th><th>Luminosity</th></tr></thead>
  <tbody>{rows}</tbody>
</table></div>
</section>
"""


def _body_row(row, indent=""):
    return (
        "<tr>"
        f'<td>{indent}{esc(row["name"])}</td>'
        f'<td>{esc(row["planet_class"])}</td>'
        f'<td>{"Gas Giant" if row["body_type"] == "g" else "Terrestrial"}</td>'
        f'<td>{esc(row["zone"])}</td>'
        f'<td>{esc(row["table_distance"])}</td>'
        f'<td>{esc(row["table_period"])}</td>'
        f'<td>{esc(row["table_gravity"])}</td>'
        "</tr>"
    )


def _planet_rows(conn, planet):
    """One row for `planet` plus one indented row per moon -- moons live
    in their own table as of schema v2 and never have moons of their own
    (`Planet.__init__` only calls `generate_moons` `if not self.is_moon`),
    so this never needs to recurse further."""
    rows = [_body_row(planet)]
    moons = fetch_all(conn, "SELECT * FROM moons WHERE planet_id = ? ORDER BY orbital_index", (planet["id"],))
    rows.extend(_body_row(moon, indent="&nbsp;&nbsp;&nbsp;&nbsp;└ ") for moon in moons)
    return rows


def _bodies_html(conn, system_id):
    planets = fetch_all(
        conn,
        "SELECT * FROM planets WHERE star_system_id = ? ORDER BY orbital_index",
        (system_id,),
    )
    belts = fetch_all(
        conn,
        "SELECT * FROM asteroid_belts WHERE star_system_id = ? ORDER BY orbital_index",
        (system_id,),
    )

    planet_rows = []
    for planet in planets:
        planet_rows.extend(_planet_rows(conn, planet))
    planet_html = ""
    if planet_rows:
        planet_html = f"""
<section class="panel">
<h2>Planets &amp; Moons</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Class</th><th>Type</th><th>Zone</th><th>Distance</th><th>Period</th><th>Gravity</th></tr></thead>
  <tbody>{''.join(planet_rows)}</tbody>
</table></div>
</section>
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
<section class="panel">
<h2>Asteroid Belts</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Density</th><th>Distance</th><th>Composition</th></tr></thead>
  <tbody>{belt_rows}</tbody>
</table></div>
</section>
"""

    return planet_html + belt_html


def _description_html(db_name, system_id, view, fmt, markdown_content, wikitext_content):
    """Builds the description section: rendered HTML by default, or the
    raw wikitext/Markdown source (for copy-pasting into a wiki) when
    `view=source`."""
    base_url = f"system.py?db={esc(db_name)}&id={system_id}"

    if view == "source":
        content = wikitext_content if fmt == "wikitext" else markdown_content
        other_fmt = "markdown" if fmt == "wikitext" else "wikitext"
        rows = min((content or "").count("\n") + 3, 40)
        return f"""
<section class="panel">
<div class="panel-header">
  <h2>Description</h2>
  <div class="view-toggle">
    <a href="{base_url}">Rendered</a>
    <span class="view-toggle-current">{esc(fmt.capitalize())} source</span>
    <a href="{base_url}&view=source&format={other_fmt}">{esc(other_fmt.capitalize())} source</a>
  </div>
</div>
<textarea readonly rows="{rows}" class="content-box" aria-label="{esc(fmt)} source">{esc(content)}</textarea>
</section>
"""

    rendered = markdown_to_html(markdown_content)
    return f"""
<section class="panel">
<div class="panel-header">
  <h2>Description</h2>
  <div class="view-toggle">
    <span class="view-toggle-current">Rendered</span>
    <a href="{base_url}&view=source">Wikitext source</a>
    <a href="{base_url}&view=source&format=markdown">Markdown source</a>
  </div>
</div>
<article class="prose">
{rendered}
</article>
</section>
"""


def handler():
    params = query_params()
    db_name = params.get("db", "")
    system_id = params.get("id", "")
    view = params.get("view", "rendered")
    if view not in ("rendered", "source"):
        view = "rendered"
    fmt = params.get("format", "wikitext")
    if fmt not in ("wikitext", "markdown"):
        fmt = "wikitext"

    path = resolve_db_path(db_name)
    conn = open_readonly(path)
    try:
        system = fetch_one(conn, "SELECT * FROM star_systems WHERE id = ?", (system_id,))
        if system is None:
            raise NotFoundError(f"No such system: {system_id!r}")

        back_html = '<p class="breadcrumb"><a href="index.py">Databases</a>'
        if system["sector_id"] is not None:
            back_html += f' &rarr; <a href="sector.py?db={esc(db_name)}&id={system["sector_id"]}">Sector</a>'
        back_html += f" &rarr; {esc(system['name'])}</p>"

        summary_bits = []
        if system["quadrant"]:
            summary_bits.append(f"Quadrant {esc(system['quadrant'])}")
        summary_bits.append("Binary system" if system["is_binary"] else "Single star")
        summary_html = "<p class=\"badges\">" + "".join(
            f'<span class="badge">{bit}</span>' for bit in summary_bits
        ) + "</p>"

        stars_html = _stars_html(conn, system_id)
        bodies_html = _bodies_html(conn, system_id)
        description_html = _description_html(
            db_name, system_id, view, fmt, system["markdown_content"], system["wikitext_content"]
        )
    finally:
        conn.close()

    body = f"""
{back_html}
{summary_html}
{description_html}
{stars_html}
{bodies_html}
"""
    return f"System: {system['name']}", body


run(handler)
