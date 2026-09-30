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


_BODY_TYPE_LABELS = {"t": "Terrestrial", "g": "Gas Giant"}


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


def _type_chip(body):
    """One chip for what kind of world this is: "Habitable" (always
    terrestrial), else "Terrestrial" or "Gas Giant", never two of them."""
    if body["habitable"]:
        return '<span class="flag flag-yes">Habitable</span>'
    return f'<span class="stat">{_BODY_TYPE_LABELS.get(body["body_type"], "")}</span>'


def _stat(text):
    return f'<span class="stat">{text}</span>' if text else ""


def _gravity_text(gravity_g):
    return f"{round(gravity_g, 3)} g" if gravity_g is not None else ""


def _planet_row_html(body, sections, is_moon=False):
    """
    A planet's (or moon's) row: its class, one type chip, a "Habitable
    moon" chip when one of its moons is habitable, "Inhabited" when it is,
    and the zone, distance, period and gravity the old Planets & Moons
    table showed. A planet's moons follow in their own "N moons" group
    right under its row, collapsed, needing no script.
    """
    moons = [] if is_moon else (body.get("moons") or [])
    stats = [
        _stat(f'Class {esc(body["planet_class"])}' if body["planet_class"] else ""),
        _type_chip(body),
    ]
    if any(moon["habitable"] for moon in moons):
        stats.append('<span class="flag flag-yes">Habitable moon</span>')
    if body["inhabited"]:
        stats.append('<span class="flag flag-yes">Inhabited</span>')
    stats += [
        _stat(esc(body.get("zone") or "")),
        # Not esc()'d: built from floats and fixed unit literals.
        _stat(format_body_distance(body["distance_km"], is_moon)),
        _stat(format_period(body["period_years"]) if body.get("period_years") is not None else ""),
        _stat(_gravity_text(body.get("gravity_g"))),
    ]
    section = sections["moons" if is_moon else "planets"]
    after_html = ""
    if moons:
        label = f'{len(moons)} moon{"s" if len(moons) != 1 else ""}'
        moon_rows = "".join(_planet_row_html(moon, sections, is_moon=True) for moon in moons)
        after_html = (
            f'<details class="moon-group"><summary>{label} of {esc(body["name"])}</summary>'
            f'<ul class="system-list">{moon_rows}</ul></details>'
        )
    return _row_html(esc(body["name"]), stats, section.get(str(body["id"])), after_html, children_visible=True)


def _belt_row_html(belt, sections):
    stats = [
        _stat(esc(belt["density"]).capitalize() if belt.get("density") else ""),
        _stat(format_distance_km(belt["distance_km"])),
    ]
    if belt.get("lower_limit_km") is not None and belt.get("upper_limit_km") is not None:
        stats.append(_stat(
            f'{format_distance_km(belt["lower_limit_km"])} to {format_distance_km(belt["upper_limit_km"])}'
        ))
    return _row_html("Asteroid Belt", stats, sections["belts"].get(str(belt["id"])))


def comet_orbit_key_km(comet):
    """
    Where a comet sorts among a star's planets and belts: its semi-major
    axis, `perihelion_distance_km / (1 - eccentricity)`. A parabolic (or
    hyperbolic) comet has no finite axis, so it sorts after every bound
    body, by perihelion.

    Returns:
        tuple: `(0, axis_km)` for a bound orbit, `(1, perihelion_km)`
               otherwise.
    """
    eccentricity = comet.get("eccentricity") or 0.0
    perihelion_km = comet["perihelion_distance_km"] or 0.0
    if comet.get("orbit_type") == "elliptical" and eccentricity < 1.0:
        return (0, perihelion_km / (1.0 - eccentricity))
    return (1, perihelion_km)


def _comet_row_html(comet, sections):
    kind = "Elliptical" if comet["orbit_type"] == "elliptical" else "Parabolic"
    stats = [
        _stat(f"{kind} comet"),
        _stat("Active" if comet["is_active"] else "Dormant"),
        _stat(f'Perihelion {format_distance_km(comet["perihelion_distance_km"])}'),
        _stat(format_period(comet["orbital_period_years"]) if comet.get("orbital_period_years") is not None
              else "Single apparition"),
    ]
    return _row_html(esc(comet["name"]), stats, sections["comets"].get(str(comet["id"])))


def _orbiting_rows_html(planets, belts, comets, sections):
    """
    A star's (or a close pair's) own bodies in their order out from the
    star: planets and belts by their shared `orbital_index`, and each comet
    slotted in before the first planet or belt farther out than its
    semi-major axis (`comet_orbit_key_km`); unbound comets last.
    """
    ordered = sorted(
        [("planet", p) for p in planets] + [("belt", b) for b in belts],
        key=lambda item: item[1]["orbital_index"],
    )
    pending = sorted(comets, key=comet_orbit_key_km)
    rows = []
    for kind, body in ordered:
        while pending and comet_orbit_key_km(pending[0]) < (0, body["distance_km"] or 0.0):
            rows.append(_comet_row_html(pending.pop(0), sections))
        rows.append(_planet_row_html(body, sections) if kind == "planet" else _belt_row_html(body, sections))
    rows.extend(_comet_row_html(comet, sections) for comet in pending)
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
