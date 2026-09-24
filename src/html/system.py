#!/usr/bin/env python3
# html/system.py

"""
System detail page: the system rendered natively as an expandable list --
its stars, planets (moons nested under each), asteroid belts and comets,
each row showing compact stats (class, type, habitable, inhabited) and
opening onto that body's own generated description -- plus the stars/
bodies tables and the system map.

Nothing about the page text is stored (schema v28): the list's
descriptions come from `GET /api/systems/<id>/sections`, and the
Wikitext/Markdown buttons (`code=wikitext|markdown`) show the full wiki
page from `GET /api/systems/<id>/text` in a code box with a Copy button
(`static/copycode.js`), both rendered from the database rows on demand.

Once this system has been uploaded to a wiki (`star_systems.wikijs_url`/
`mediawiki_url` -- see `schema.sql`'s "v22" header note), the System
panel links to that page (opening in a new tab). An admin session
(`auth_me`) additionally gets an "Upload to Wiki" form offering whichever
backend(s) are both configured deployment-wide (`GET /api/wiki-config`)
and not yet uploaded to.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from apiclient import (
    ApiError, auth_me, get_system, get_system_sections, get_system_text, get_wiki_config, upload_system_to_wiki,
)
from fmt import esc, linkify_location, nearest_neighbors_location, post_link
from mdconvert import markdown_to_html
from page import form_params, incoming_cookie_header, nav_params, run
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


def _planets_table_html(planets, heading="Planets &amp; Moons"):
    planet_rows = []
    for planet in planets:
        planet_rows.extend(_planet_rows(planet))
    if not planet_rows:
        return ""
    return f"""
<section class="panel">
<h2>{heading}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Class</th><th>Type</th><th>Zone</th><th>Distance</th><th>Period</th><th>Gravity</th></tr></thead>
  <tbody>{''.join(planet_rows)}</tbody>
</table></div>
</section>
"""


def _belts_table_html(belts, heading="Asteroid Belts"):
    if not belts:
        return ""
    belt_rows = "".join(
        "<tr>"
        f'<td>{esc(row["density"])}</td>'
        f'<td>{row["distance_km"]:,.0f} km</td>'
        f'<td>{esc(row["composition_summary"])}</td>'
        "</tr>"
        for row in belts
    )
    return f"""
<section class="panel">
<h2>{heading}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Density</th><th>Distance</th><th>Composition</th></tr></thead>
  <tbody>{belt_rows}</tbody>
</table></div>
</section>
"""


def _comets_table_html(comets, heading="Comets"):
    if not comets:
        return ""
    comet_rows = "".join(
        "<tr>"
        f'<td>{esc(row["name"])}</td>'
        f'<td>{"Elliptical" if row["orbit_type"] == "elliptical" else "Parabolic"}</td>'
        # Not esc()'d -- see the same note in `_stars_html`.
        f'<td>{format_body_distance(row["perihelion_distance_km"], is_moon=False)}</td>'
        f'<td>{row["eccentricity"]:.3f}</td>'
        f'<td>{format_period(row["orbital_period_years"]) if row["orbital_period_years"] is not None else "single apparition"}</td>'
        f'<td>{"Active" if row["is_active"] else "Dormant"}</td>'
        f'<td>{esc(row["composition_summary"])}</td>'
        "</tr>"
        for row in comets
    )
    return f"""
<section class="panel">
<h2>{heading}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Orbit</th><th>Perihelion</th><th>Eccentricity</th><th>Period</th><th>Activity</th><th>Composition</th></tr></thead>
  <tbody>{comet_rows}</tbody>
</table></div>
</section>
"""


def _bodies_html(planets, belts, comets, stars, binary_configuration):
    """
    Builds the "Planets & Moons"/"Asteroid Belts"/"Comets" section(s).

    For a `'wide'` (S-type) binary, `planets`/`belts`/`comets` belong to
    two different, independent stars (disambiguated by each row's own
    `star_id`, matched against `stars`' own `id` -- see
    `queryDb.system_detail`'s docstring) -- rendering them in one flat
    table the way a single star's or a `'close'` (P-type) pair's bodies
    already are would misrepresent which star each one actually orbits
    (a `'close'` pair's planets genuinely have no single owning star --
    they orbit the merged pair together -- so that case is unaffected).
    Grouped into one labeled section per star instead, in the same order
    `stars` already comes in (primary first).
    """
    if binary_configuration != "wide":
        return _planets_table_html(planets) + _belts_table_html(belts) + _comets_table_html(comets)

    sections = []
    for star in stars:
        star_planets = [p for p in planets if p["star_id"] == star["id"]]
        star_belts = [b for b in belts if b["star_id"] == star["id"]]
        star_comets = [c for c in comets if c["star_id"] == star["id"]]
        label = f'{esc(star["name"])} ({esc(star["role"])})'
        sections.append(_planets_table_html(star_planets, heading=f"Planets &amp; Moons — {label}"))
        sections.append(_belts_table_html(star_belts, heading=f"Asteroid Belts — {label}"))
        sections.append(_comets_table_html(star_comets, heading=f"Comets — {label}"))
    return "".join(sections)


_BODY_TYPE_LABELS = {"t": "Terrestrial", "g": "Gas giant"}


def _flag_html(label, value):
    """One yes/no chip in a system-list row, e.g. "Habitable: Yes"."""
    state = "yes" if value else "no"
    return f'<span class="flag flag-{state}">{label}: {"Yes" if value else "No"}</span>'


def _row_html(title, stats, markdown, children_html="", children_visible=False):
    """
    One clickable row of the system list: a native `<details>` whose
    summary line is the body's name plus its compact stats, opening onto
    that body's own generated description. Needs no script.

    `children_html` (nested rows) sits inside the `<details>` -- shown only
    once the row is opened, as a planet's moons are -- or, with
    `children_visible`, right after it, always shown, as a wide pair's
    star's own bodies are.
    """
    stats_html = "".join(stats)
    description = markdown_to_html(markdown) if markdown else ""
    inside, after = ("", children_html) if children_visible else (children_html, "")
    return f"""
<li><details class="body-row">
<summary><span class="body-name">{title}</span><span class="body-stats">{stats_html}</span></summary>
<div class="body-detail prose">{description}</div>
{inside}
</details>{after}</li>"""


def _star_row_html(star, sections, children_html=""):
    stats = [f'<span class="stat">{esc(star["star_type"])}</span>']
    if star["role"] != "single":
        stats.append(f'<span class="stat">{esc(star["role"].capitalize())}</span>')
    return _row_html(
        esc(star["name"]), stats, sections["stars"].get(str(star["id"])), children_html, children_visible=True,
    )


def _planet_row_html(body, sections, is_moon=False):
    stats = [
        f'<span class="stat">Class {esc(body["planet_class"])}</span>' if body["planet_class"] else "",
        f'<span class="stat">{_BODY_TYPE_LABELS.get(body["body_type"], "")}</span>',
        _flag_html("Habitable", body["habitable"]),
        _flag_html("Inhabited", body["inhabited"]),
    ]
    children_html = ""
    if is_moon:
        markdown = sections["moons"].get(str(body["id"]))
    else:
        markdown = sections["planets"].get(str(body["id"]))
        moons = body.get("moons") or []
        if moons:
            stats.append(f'<span class="stat">{len(moons)} moon{"s" if len(moons) != 1 else ""}</span>')
            children_html = '<ul class="system-list">' + "".join(
                _planet_row_html(moon, sections, is_moon=True) for moon in moons
            ) + "</ul>"
    return _row_html(esc(body["name"]), stats, markdown, children_html)


def _belt_row_html(belt, sections):
    stats = [f'<span class="stat">{esc(belt["density"]).capitalize()}</span>'] if belt.get("density") else []
    return _row_html("Asteroid Belt", stats, sections["belts"].get(str(belt["id"])))


def _comet_row_html(comet, sections):
    kind = "Elliptical" if comet["orbit_type"] == "elliptical" else "Parabolic"
    stats = [
        f'<span class="stat">{kind} comet</span>',
        f'<span class="stat">{"Active" if comet["is_active"] else "Dormant"}</span>',
    ]
    return _row_html(esc(comet["name"]), stats, sections["comets"].get(str(comet["id"])))


def _orbiting_rows_html(planets, belts, comets, sections):
    """A star's (or a close pair's) own bodies, in orbital order -- planets
    and belts share one `orbital_index` space per star -- then comets."""
    ordered = sorted(
        [("planet", p) for p in planets] + [("belt", b) for b in belts],
        key=lambda item: item[1]["orbital_index"],
    )
    rows = [
        _planet_row_html(body, sections) if kind == "planet" else _belt_row_html(body, sections)
        for kind, body in ordered
    ]
    rows.extend(_comet_row_html(comet, sections) for comet in comets)
    return "".join(rows)


def _system_list_html(system, sections):
    """
    The system rendered natively: the page's overview (a binary pair's
    own data, the system summary, any flavor text) above an expandable
    list of its stars, planets (with their moons nested under them),
    asteroid belts and comets. Each row shows compact stats and opens onto
    that body's own generated description -- see
    `stellarObjects.systemRender.render_system_sections`.

    A `'wide'` (S-type) pair's bodies each orbit one of its two stars, so
    they're nested under that star's own row; a single star's or a
    `'close'` (P-type) pair's bodies orbit the whole system and follow
    the star rows at the top level.
    """
    stars, planets, belts, comets = system["stars"], system["planets"], system["belts"], system["comets"]
    if system.get("binary_configuration") == "wide":
        rows = []
        for star in stars:
            children = _orbiting_rows_html(
                [p for p in planets if p["star_id"] == star["id"]],
                [b for b in belts if b["star_id"] == star["id"]],
                [c for c in comets if c["star_id"] == star["id"]],
                sections,
            )
            children_html = f'<ul class="system-list">{children}</ul>' if children else ""
            rows.append(_star_row_html(star, sections, children_html))
        rows_html = "".join(rows)
    else:
        rows_html = "".join(_star_row_html(star, sections) for star in stars)
        rows_html += _orbiting_rows_html(planets, belts, comets, sections)

    overview_html = markdown_to_html(sections["overview"]) if sections["overview"] else ""
    return f"""
<div class="prose system-overview">{overview_html}</div>
<ul class="system-list system-list-root">{rows_html}</ul>
"""


def _code_html(code_fmt, content):
    """The generated wiki page in a read-only code box, with a Copy button
    (`static/copycode.js` -- the page's Content-Security-Policy allows no
    inline script)."""
    label = "Wikitext" if code_fmt == "wikitext" else "Markdown"
    rows = min(content.count("\n") + 3, 30)
    return f"""
<div class="code-view">
  <div class="code-view-header">
    <span class="view-toggle-current">{label}</span>
    <button type="button" class="btn copy-btn" data-copy-target="system-code">Copy</button>
  </div>
  <textarea readonly rows="{rows}" id="system-code" class="content-box" aria-label="{label} source">{esc(content)}</textarea>
</div>
"""


def _system_section_html(db_name, system_id, system, sections, code_fmt, code_content):
    """
    The "System" panel: the natively rendered system list, with Wikitext/
    Markdown buttons that show the page's generated code (rendered on
    demand -- nothing is stored, see `schema.sql`'s "v28" header note) in
    a code box above it, plus links to the wiki copies once uploaded.
    """
    base_params = {"db": db_name, "id": system_id}
    buttons = []
    for fmt, label in (("wikitext", "Wikitext"), ("markdown", "Markdown")):
        if fmt == code_fmt:
            buttons.append(post_link("system.py", base_params, f"Hide {label}", css_class="btn btn-active"))
        else:
            buttons.append(post_link("system.py", {**base_params, "code": fmt}, label, css_class="btn"))

    wiki_links = []
    if system["wikijs_url"]:
        wiki_links.append(
            f'<a href="{esc(system["wikijs_url"])}" target="_blank" rel="noopener noreferrer">View on Wiki.js</a>'
        )
    if system["mediawiki_url"]:
        wiki_links.append(
            f'<a href="{esc(system["mediawiki_url"])}" target="_blank" rel="noopener noreferrer">View on MediaWiki</a>'
        )
    wiki_html = f'<p class="hint">Published on the wiki: {" &middot; ".join(wiki_links)}</p>' if wiki_links else ""
    code_html = _code_html(code_fmt, code_content) if code_fmt else ""

    return f"""
<section class="panel" id="system-panel">
<div class="panel-header">
  <h2>System</h2>
  <div class="view-toggle">{"".join(buttons)}</div>
</div>
{wiki_html}
{code_html}
{_system_list_html(system, sections)}
</section>
"""


def _wiki_upload_section_html(db_name, system_id, system, wiki_config, wiki_message, wiki_error):
    """
    Builds the "Upload to Wiki" form -- one radio option per backend that
    is both configured deployment-wide (`wiki_config`, `GET
    /api/wiki-config`) and not yet uploaded to for this system (this
    system's own `wikijs_url`/`mediawiki_url`). Returns just the
    message/error (no form at all) once every configured backend has
    already been uploaded to, or none are configured -- the caller only
    reaches this for an authenticated admin session in the first place
    (see the module-level POST handling below `handler`).
    """
    message_html = f'<p class="hint">{esc(wiki_message)}</p>' if wiki_message else ""
    error_html = f'<p class="error">{esc(wiki_error)}</p>' if wiki_error else ""

    options = []
    if wiki_config.get("wikijs") and not system["wikijs_url"]:
        options.append(("wikijs", "Wiki.js"))
    if wiki_config.get("mediawiki") and not system["mediawiki_url"]:
        options.append(("mediawiki", "MediaWiki"))
    if not options:
        return f"{message_html}{error_html}"

    radios = " ".join(
        f'<label><input type="radio" name="backend" value="{value}"{" checked" if i == 0 else ""}> {label}</label>'
        for i, (value, label) in enumerate(options)
    )
    return f"""
{message_html}{error_html}
<section class="panel">
<h2>Upload to Wiki</h2>
<form method="post" action="system.py" class="search-form">
  <input type="hidden" name="action" value="upload_wiki">
  <input type="hidden" name="db" value="{esc(db_name)}">
  <input type="hidden" name="id" value="{esc(system_id)}">
  <div class="search-fields">
    <div class="search-field">{radios}</div>
    <label class="search-field">Path (Wiki.js only -- MediaWiki uses this system's name)
      <input type="text" name="path" placeholder="e.g. systems/{esc(system['name'])}">
    </label>
  </div>
  <div class="search-actions">
    <button type="submit" class="btn">Upload</button>
  </div>
</form>
</section>
"""


def handler():
    params = nav_params()
    db_name = params.get("db", "")
    system_id = params.get("id", "")
    code_fmt = params.get("code")
    if code_fmt not in ("wikitext", "markdown"):
        code_fmt = None

    system = get_system(db_name, system_id)
    sections = get_system_sections(db_name, system_id)
    code_content = get_system_text(db_name, system_id, code_fmt)["content"] if code_fmt else ""

    back_html = f'<p class="breadcrumb">{post_link("browse.py", {"db": db_name}, esc(db_name))}'
    if system["sector_id"] is not None:
        sector_link = post_link("sector.py", {"db": db_name, "id": system["sector_id"]}, "Sector")
        back_html += f" &rarr; {sector_link}"
    back_html += f" &rarr; {esc(system['name'])}</p>"

    summary_bits = []
    if system["quadrant"]:
        summary_bits.append(f"Octant {esc(system['quadrant'])}")
    if system["is_binary"]:
        config_label = {"close": " (close)", "wide": " (wide)"}.get(system.get("binary_configuration"), "")
        summary_bits.append(f"Binary system{config_label}")
    else:
        summary_bits.append("Single star")
    summary_html = "<p class=\"badges\">" + "".join(
        f'<span class="badge">{bit}</span>' for bit in summary_bits
    ) + "</p>"

    wiki_upload_html = ""
    if identity is not None:
        wiki_config = get_wiki_config()
        wiki_upload_html = _wiki_upload_section_html(db_name, system_id, system, wiki_config, wiki_message, wiki_error)

    nav_html = ""
    if system["sector_id"] is not None:
        # NAV needs a sector to measure a position from at all -- see
        # queryDb.nav_between's own availability rules, which this only
        # pre-checks the first (cheapest) condition of. A system in a
        # non-galaxy-placed sector still gets the links: same-sector NAV
        # is always available once that much is true, nav.py itself
        # works out whether cross-sector NAV also applies. "From" sets
        # nav.py's origin directly; "To" sets only its destination and
        # lets nav.py prompt for an origin -- the symmetric entry point
        # nav.py's own origin picker now supports.
        nav_html = (
            post_link("nav.py", {"db": db_name, "from": system_id}, "Navigate from here", css_class="btn")
            + post_link("nav.py", {"db": db_name, "to": system_id}, "Navigate to here", css_class="btn")
        )

    location_html = ""
    neighbors = system.get("nearest_neighbors")
    if system["location"] and neighbors:
        location_html = (
            f'<span class="location">Location: '
            f'{nearest_neighbors_location(db_name, system["location"], neighbors)}</span>'
        )
    elif system["location"]:
        # An API without `nearest_neighbors`, or a system with no stored
        # position: fall back to linking whichever names in the stored
        # string still match a system in the same sector.
        name_to_id = {row["name"]: row["id"] for row in system["sector_siblings"]}
        location_html = (
            f'<span class="location">Location: '
            f'{linkify_location(db_name, system["location"], name_to_id)}</span>'
        )

    # One compact flex row (breadcrumb + Octant/binary badges + the Navigate
    # button + nearest-neighbor location) instead of four separately
    # stacked, vertically spread-out blocks -- see static/style.css's
    # `.page-subhead` rule.
    subhead_html = f'<div class="page-subhead">{back_html}{summary_html}{nav_html}{location_html}</div>'

    map_html = ""
    if system["stars"]:
        map_html = render_system_map_panel(system, system["stars"], system["planets"], system["belts"])
    stars_html = _stars_html(system["stars"])
    bodies_html = _bodies_html(
        system["planets"], system["belts"], system["comets"], system["stars"], system.get("binary_configuration")
    )
    system_html = _system_section_html(db_name, system_id, system, sections, code_fmt, code_content)

    body = f"""
{subhead_html}
{map_html}
{wiki_upload_html}
{system_html}
{stars_html}
{bodies_html}
<script type="module" src="static/systemmap.js"></script>
<script type="module" src="static/copycode.js"></script>
"""
    return f"System: {system['name']}", body


cookie_header = incoming_cookie_header()
try:
    identity = auth_me(cookie_header)
except ApiError:
    # Fails quiet, same as _sidenav_html's own admin-session check --
    # this only decides whether the "Upload to Wiki" section renders at
    # all; handler()'s own read of the system itself will hit (and report)
    # the same API-unreachable failure a moment later.
    identity = None

wiki_message = None
wiki_error = None
if identity is not None and os.environ.get("REQUEST_METHOD", "GET").upper() == "POST":
    fields = form_params()
    if fields.get("action") == "upload_wiki":
        db_name = fields.get("db", "")
        system_id = fields.get("id", "")
        backend = fields.get("backend", "")
        path = fields.get("path", "").strip() or None
        try:
            page = upload_system_to_wiki(cookie_header, db_name, system_id, backend, path)
            wiki_message = f"Uploaded to the wiki: {page['url']}"
        except ApiError as exc:
            wiki_error = "A page already exists at that location." if exc.status_code == 409 else str(exc)

run(handler)
