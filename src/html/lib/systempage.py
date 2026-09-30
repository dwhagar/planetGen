# html/lib/systempage.py

"""
The HTML pieces of the system page (`/system/<id>`, `web/system_pages.py`)
built in Python rather than in its Jinja template: the natively rendered
system list (stars, planets with their moons nested under them, asteroid
belts and comets, each a native `<details>` row opening onto that body's
own generated description) and the Stars / Planets & Moons / Asteroid
Belts / Comets tables.

Moved here unchanged from the old CGI `system.py`. Every value taken from
the database is escaped with `fmt.esc`; the only raw markup comes from
`tabledisplay.py`'s unit formatters (built from floats and fixed unit
literals, e.g. `<sup>`) and from `mdconvert.markdown_to_html`, which
escapes its own input. So the page wraps what these return in
`trusted_html`.
"""

from fmt import esc, format_distance_km
from mdconvert import markdown_to_html
from tabledisplay import (
    format_body_distance, format_period, format_star_luminosity, format_star_mass, format_star_radius,
)


def stars_html(stars):
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
<div class="table-scroll" tabindex="0"><table>
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
        # Not esc()'d -- see the same note in `stars_html`.
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
<div class="table-scroll" tabindex="0"><table>
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
        f'<td>{format_distance_km(row["distance_km"])}</td>'
        f'<td>{esc(row["composition_summary"])}</td>'
        "</tr>"
        for row in belts
    )
    return f"""
<section class="panel">
<h2>{heading}</h2>
<div class="table-scroll" tabindex="0"><table>
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
        # Not esc()'d -- see the same note in `stars_html`.
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
<div class="table-scroll" tabindex="0"><table>
  <thead><tr><th>Name</th><th>Orbit</th><th>Perihelion</th><th>Eccentricity</th><th>Period</th><th>Activity</th><th>Composition</th></tr></thead>
  <tbody>{comet_rows}</tbody>
</table></div>
</section>
"""


# TODO(web-pages #48): the Planets & Moons, Asteroid Belts and Comets
# tables repeat the system list; remove them once the list rows carry their
# columns (class, type, zone, distance, period, gravity) as stats.
def bodies_html(planets, belts, comets, stars, binary_configuration):
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


# TODO(system-list #2): a non-terrestrial planet is never habitable, and a
# habitable one is always terrestrial, so show one type chip: "Gas Giant",
# "Terrestrial" or "Habitable" (never Terrestrial and Habitable together).
# Add a new chip, shown only when one of the planet's moons is habitable
# (e.g. "Habitable moon"); queryDb._with_life_fields already marks each
# moon. _body_row above (the table view) follows the same rule.
# TODO(web-pages #48): moons are buried after the planet's description
# inside its <details>; give them their own expandable group right under
# the planet's row (a nested "N moons" <details>), no script needed.
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


# TODO(system-list #2): the list shows only "Asteroid Belt" and density;
# add its distance from the star in the most meaningful unit (distances
# #1).
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


# TODO(web-pages #48): comets go in their relative order from the star too,
# not appended last: sort them in with planets and belts by distance (a
# comet by its semi-major axis, perihelion_distance_km / (1 -
# eccentricity); parabolic ones last, by perihelion).
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


def system_list_html(system, sections):
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
