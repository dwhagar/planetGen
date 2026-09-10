#!/usr/bin/env python3
# html/system.py

"""
System detail page: stars, planets/moons, asteroid belts, and the
system's full description. Defaults to rendering `markdown_content` as
actual HTML (via `mdconvert.markdown_to_html_with_headings`) so the
description reads like a normal page instead of a wall of raw Markdown,
with a table-of-contents linking to each heading (a collapsed pulldown on
narrow windows, a fixed sidebar on wide ones); `?view=source`
switches to the original raw-text view (wikitext or Markdown, toggled via
`&format=`), which is what you want when copy-pasting into a wiki.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from apiclient import get_system
from fmt import esc, linkify_location
from mdconvert import markdown_to_html_with_headings
from page import query_params, run
from systemmap import render_system_map_panel
from tabledisplay import (
    format_body_distance, format_period, format_star_luminosity, format_star_mass, format_star_radius,
)


def _stars_html(stars):
    rows = "".join(
        "<tr>"
        f'<td>{esc(star["role"])}</td>'
        f'<td>{esc(star["name"])}</td>'
        f'<td>{esc(star["star_type"])}</td>'
        # Not esc()'d: these are built entirely from floats and fixed unit
        # literals (never database TEXT), and legitimately contain a raw
        # `<sup>exponent</sup>` -- see `tabledisplay.py`.
        f'<td>{format_star_mass(star["mass_kg"])}</td>'
        f'<td>{format_star_radius(star["radius_km"])}</td>'
        f'<td>{int(star["temperature_k"])} K</td>'
        f'<td>{format_star_luminosity(star["luminosity_w"])}</td>'
        "</tr>"
        for star in stars
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


def _body_row(row, is_moon, indent=""):
    return (
        "<tr>"
        f'<td>{indent}{esc(row["name"])}</td>'
        f'<td>{esc(row["planet_class"])}</td>'
        f'<td>{"Gas Giant" if row["body_type"] == "g" else "Terrestrial"}</td>'
        f'<td>{esc(row["zone"])}</td>'
        # Not esc()'d -- see the same note in `_stars_html`.
        f'<td>{format_body_distance(row["distance_km"], is_moon)}</td>'
        f'<td>{format_period(row["period_years"])}</td>'
        f'<td>{round(row["gravity_g"], 3) if row["gravity_g"] is not None else ""} g</td>'
        "</tr>"
    )


def _planet_rows(planet):
    """One row for `planet` plus one indented row per moon -- moons never
    have moons of their own, so this never needs to recurse further."""
    rows = [_body_row(planet, is_moon=False)]
    rows.extend(
        _body_row(moon, is_moon=True, indent="&nbsp;&nbsp;&nbsp;&nbsp;└ ")
        for moon in planet["moons"]
    )
    return rows


def _bodies_html(planets, belts):
    planet_rows = []
    for planet in planets:
        planet_rows.extend(_planet_rows(planet))
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


def _toc_html(headings):
    """
    Builds the table-of-contents linking to each heading
    `markdown_to_html_with_headings` found in the rendered description --
    skipped entirely when there's nothing worth a contents list for (just
    the system's own top-level heading, or no description at all).

    A checkbox-driven disclosure: collapsed into a pulldown above the
    prose by default (so it doesn't eat vertical space next to the
    reading column), or -- once the window is wide enough to have real
    margin space beyond the centered content column -- fixed in the
    right-hand margin and always shown open (see `.toc` in style.css,
    which also explains why this isn't a native <details>/<summary>).
    Only one of these ever exists on a page, so the fixed `id` below
    never collides.
    """
    if len(headings) <= 1:
        return ""
    items = "".join(
        f'<li class="toc-level-{heading["level"]}"><a href="#{heading["id"]}">{esc(heading["text"])}</a></li>'
        for heading in headings
    )
    return f"""
<nav class="toc" aria-label="Table of contents">
<input type="checkbox" id="toc-toggle" class="toc-toggle">
<label for="toc-toggle" class="toc-title">Contents</label>
<ul>{items}</ul>
</nav>
"""


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

    rendered, headings = markdown_to_html_with_headings(markdown_content)
    toc_html = _toc_html(headings)
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
{toc_html}
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

    system = get_system(db_name, system_id)

    back_html = '<p class="breadcrumb"><a href="index.py">Databases</a>'
    if system["sector_id"] is not None:
        back_html += f' &rarr; <a href="sector.py?db={esc(db_name)}&id={system["sector_id"]}">Sector</a>'
    back_html += f" &rarr; {esc(system['name'])}</p>"

    summary_bits = []
    if system["quadrant"]:
        summary_bits.append(f"Octant {esc(system['quadrant'])}")
    summary_bits.append("Binary system" if system["is_binary"] else "Single star")
    summary_html = "<p class=\"badges\">" + "".join(
        f'<span class="badge">{bit}</span>' for bit in summary_bits
    ) + "</p>"

    nav_html = ""
    if system["sector_id"] is not None:
        # NAV needs a sector to measure a position from at all -- see
        # queryDb.nav_between's own availability rules, which this only
        # pre-checks the first (cheapest) condition of. A system in a
        # non-galaxy-placed sector still gets the link: same-sector NAV
        # is always available once that much is true, nav.py itself
        # works out whether cross-sector NAV also applies.
        nav_html = (
            f'<p><a class="btn" href="nav.py?db={esc(db_name)}&from={system_id}">Navigate from here</a></p>'
        )

    location_html = ""
    if system["location"]:
        name_to_id = {row["name"]: row["id"] for row in system["sector_siblings"]}
        location_html = (
            f'<p class="location">Location: '
            f'{linkify_location(db_name, system["location"], name_to_id)}</p>'
        )

    map_html = ""
    if system["stars"]:
        map_html = render_system_map_panel(system, system["stars"], system["planets"], system["belts"])
    stars_html = _stars_html(system["stars"])
    bodies_html = _bodies_html(system["planets"], system["belts"])
    description_html = _description_html(
        db_name, system_id, view, fmt, system["markdown_content"], system["wikitext_content"]
    )

    body = f"""
{back_html}
{summary_html}
{nav_html}
{location_html}
{map_html}
{description_html}
{stars_html}
{bodies_html}
<script src="static/systemmap.js" defer></script>
"""
    return f"System: {system['name']}", body


run(handler)
