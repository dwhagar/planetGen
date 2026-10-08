# planetgen/web/lib/systempage.py

"""
The HTML pieces of the system page (`/system/<id>`, `web/system_pages.py`)
built in Python rather than in its Jinja template: the natively rendered
system list (stars, planets with their moons nested under them, asteroid
belts and comets, each a native `<details>` row opening onto that body's
own generated description) and the Stars / Planets & Moons / Asteroid
Belts / Comets tables. A body's facilities (schema v42: starbases,
colonies, outposts) are listed in its own row's detail, and `facility_row`
formats one for the page's Facilities panel.

Moved here unchanged from the old CGI `system.py`. Every value taken from
the database is escaped with `fmt.esc`; the only raw markup comes from
`tabledisplay.py`'s unit formatters (built from floats and fixed unit
literals, e.g. `<sup>`) and from `mdconvert.markdown_to_html`, which
escapes its own input. So the page wraps what these return in
`trusted_html`.

Class labels (a star's type, a planet's class) link to their class
reference pages (`/classes/...`) through an optional `class_url(type_slug,
code)` hook, the way `fmt`'s location helpers take a `system_url` hook:
it returns the page's URL or `None`, and without it the labels stay
plain text. A link never goes in a row's `<summary>` (a link inside the
disclosure button would be a control nested in a control); it opens the
row's detail instead.
"""

from planetgen.web.lib.classref import PARABOLIC_COMET_CLASS, class_entry, star_type_classes
from planetgen.web.lib.fmt import esc, format_distance_km, format_number, format_speed_kms
from planetgen.web.lib.mdconvert import markdown_to_html
from planetgen.web.lib.tabledisplay import (
    format_body_distance, format_period, format_star_luminosity, format_star_mass, format_star_radius,
)


def _link(url, text_html):
    """`text_html` (already escaped) as a link to `url`, or unchanged
    when `url` is `None`."""
    return f'<a href="{esc(url)}">{text_html}</a>' if url else text_html


def _star_class_urls(star_type, class_url):
    """`(spectral class URL, luminosity class URL)` for a star type, each
    `None` when it can't be linked."""
    classes = star_type_classes(star_type) if class_url else None
    if classes is None:
        return None, None
    return class_url("star-spectral", classes[0]), class_url("star-luminosity", classes[1])


def stars_html(stars, class_url=None, notes_html=""):
    """The Stars table. With `class_url` (see the module docstring), a
    star's type links to its spectral class page. `notes_html`
    (`star_notes_html`) goes under the table."""
    # UX.53: "single" says nothing; the Role column is for binaries and multiples.
    show_role = any(star["role"] != "single" for star in stars)
    rows = "".join(
        "<tr>"
        + (f'<td data-label="Role">{esc(star["role"])}</td>' if show_role else "")
        + f'<td class="stack-title">{esc(star["name"])}</td>'
        f'<td data-label="Type">{_link(_star_class_urls(star["star_type"], class_url)[0], esc(star["star_type"]))}</td>'
        # Not esc()'d: these are built entirely from floats and fixed unit
        # literals (never database TEXT), and legitimately contain a raw
        # `<sup>exponent</sup>` -- see `tabledisplay.py`.
        f'<td data-label="Mass">{format_star_mass(star["mass_kg"])}</td>'
        f'<td data-label="Radius">{format_star_radius(star["radius_km"])}</td>'
        f'<td data-label="Temp">{format_number(star["temperature_k"])} K</td>'
        f'<td data-label="Luminosity">{format_star_luminosity(star["luminosity_w"])}</td>'
        "</tr>"
        for star in stars
    )
    if not rows:
        return ""
    role_head = "<th>Role</th>" if show_role else ""
    return f"""
<section class="panel">
<h2>Stars</h2>
<div class="table-scroll" tabindex="0"><table class="table-stack">
  <thead><tr>{role_head}<th>Name</th><th>Type</th><th>Mass</th><th>Radius</th><th>Temp</th><th>Luminosity</th></tr></thead>
  <tbody>{rows}</tbody>
</table></div>
{notes_html}</section>
"""


def _star_has_list_row(system, star, by_host):
    """Whether `star` still has its own row in the System list: a wide pair's
    stars do (their bodies nest under them), and so does a star that hosts a
    facility. Any other star is only in the Stars table (UX.53)."""
    return system.get("binary_configuration") == "wide" or bool((by_host or {}).get(("star", star["id"])))


def star_notes_html(system, sections, class_url=None, facilities=()):
    """What the System list's star rows used to say, for the stars that have
    no row any more (UX.53): each one's description and its spectral and
    luminosity class links, under the Stars table."""
    by_host = _facilities_by_host(facilities)
    notes = []
    for star in system["stars"]:
        if _star_has_list_row(system, star, by_host):
            continue
        text = sections["stars"].get(str(star["id"]))
        links = []
        classes = star_type_classes(star["star_type"]) if class_url else None
        if classes:
            spectral_url, luminosity_url = _star_class_urls(star["star_type"], class_url)
            links = [_link(spectral_url, f"Spectral class {esc(classes[0])}") if spectral_url else "",
                     _link(luminosity_url, f"Luminosity class {esc(classes[1])}") if luminosity_url else ""]
        links = [link for link in links if link]
        if not text and not links:
            continue
        heading = f"<h3>{esc(star['name'])}</h3>" if len(system["stars"]) > 1 else ""
        link_line = f'<p class="class-links">{" &middot; ".join(links)}</p>' if links else ""
        body = markdown_to_html(text) if text else ""
        notes.append(f'<div class="star-note prose">{heading}{body}{link_line}</div>')
    return "".join(notes)


_BODY_TYPE_LABELS = {"t": "Terrestrial", "g": "Gas Giant"}

_FACILITY_PLACEMENT_LABELS = {
    "terrestrial": "On the surface",
    "orbital": "In orbit",
    "asteroid": "Among the asteroids",
    "standalone": "Parked in space",
}
"""dict: How each `facilities.placement` reads on the page."""


def facility_kind_label(kind):
    """A `facilities.kind` as a label: "mining-colony" -> "Mining colony"."""
    return (kind or "").replace("-", " ").capitalize()


def facility_row(facility):
    """
    One facility's display values, every one already escaped or built from
    numbers (the page wraps them in `trusted_html`): `name`, `kind`,
    `placement` and, for an orbital facility, `distance`, `period` and
    `speed` (else `None`). The orbit is the stored one the API returns,
    worked out from the host's mass when it was saved.
    """
    orbital = facility.get("orbit_distance_km") is not None
    return {
        "name": esc(facility["name"]),
        "kind": esc(facility_kind_label(facility["kind"])),
        "placement": _FACILITY_PLACEMENT_LABELS.get(facility["placement"], esc(facility["placement"])),
        "distance": format_distance_km(facility["orbit_distance_km"]) if orbital else None,
        "period": esc(format_period(facility["orbit_period_years"]))
        if orbital and facility.get("orbit_period_years") is not None else None,
        "speed": format_speed_kms(facility["orbital_speed_kms"])
        if orbital and facility.get("orbital_speed_kms") is not None else None,
    }


def _facilities_html(facilities):
    """A body's own facilities as a list at the top of its row's detail,
    or `""` for none."""
    if not facilities:
        return ""
    items = []
    for facility in facilities:
        row = facility_row(facility)
        where = row["placement"]
        if row["distance"]:
            where += f' at {row["distance"]}'
            where += "".join(f", {part}" for part in (row["period"], row["speed"]) if part)
        items.append(f'<li><span class="facility-name">{row["name"]}</span> '
                     f'<span class="stat">{row["kind"]}</span> {where}</li>')
    return f'<ul class="facility-list" aria-label="Facilities">{"".join(items)}</ul>'


def _facility_stat(facilities):
    """The summary chip counting a body's facilities, or `""`."""
    if not facilities:
        return ""
    count = len(facilities)
    return _stat(f'{count} facilit{"y" if count == 1 else "ies"}')


def _facilities_by_host(facilities):
    """`{(host_type, host_id): [facility, ...]}`."""
    by_host = {}
    for facility in facilities or ():
        by_host.setdefault((facility["host_type"], facility["host_id"]), []).append(facility)
    return by_host


def _row_admin_html(admin_rows, kind, body_id):
    """The Admin menu at the end of a planet's, moon's or belt's row (UX.68):
    one item per action of that body (`admin_rows`, `{"planet:5": {"id":
    "body-2", "items": [(key, label), ...]}}`), each opening the dialog the
    system page's template (`partials/admin_menu.html`) draws for it. `""`
    where the body has none (not an admin, or a comet)."""
    row = (admin_rows or {}).get(f"{kind}:{body_id}")
    if not row:
        return ""
    items = "".join(
        f'<sl-menu-item data-dialog="admin-{esc(row["id"])}-{esc(key)}">{esc(label)}</sl-menu-item>'
        for key, label in row["items"]
    )
    return (f'<sl-dropdown class="row-admin" placement="bottom-end" distance="4">'
            f'<sl-button slot="trigger" size="small" caret>Admin</sl-button>'
            f'<sl-menu aria-label="Admin actions for {esc(row["label"])}">{items}</sl-menu></sl-dropdown>')


def _row_html(title, stats, markdown, children_html="", children_visible=False, links=(), facilities=(), admin=""):
    """
    One clickable row of the system list: a native `<details>` whose
    summary line is the body's name plus its compact stats, opening onto
    that body's own generated description. Needs no script.

    `children_html` (nested rows) sits inside the `<details>` -- shown only
    once the row is opened, as a planet's moons are -- or, with
    `children_visible`, right after it, always shown, as a wide pair's
    star's own bodies are.

    `links` (the body's class links, HTML) open the detail, above the
    description; the body's `facilities` are counted in the summary and
    listed under the links; `admin` (an Admin menu, `_row_admin_html`) sits
    at the row's right end.
    """
    stats_html = "".join(stats) + _facility_stat(facilities)
    description = _facilities_html(facilities) + (markdown_to_html(markdown) if markdown else "")
    links = [link for link in links if link]
    if links:
        description = f'<p class="class-links">{" &middot; ".join(links)}</p>{description}'
    inside, after = ("", children_html) if children_visible else (children_html, "")
    return f"""
<li{' class="has-admin"' if admin else ""}><details class="body-row">
<summary><span class="body-name">{title}</span><span class="body-stats">{stats_html}</span></summary>
<div class="body-detail prose">{description}</div>
{inside}
</details>{admin}{after}</li>"""


def _star_row_html(star, sections, children_html="", class_url=None, by_host=None):
    stats = [f'<span class="stat">{esc(star["star_type"])}</span>']
    if star["role"] != "single":
        stats.append(f'<span class="stat">{esc(star["role"].capitalize())}</span>')
    links = []
    classes = star_type_classes(star["star_type"]) if class_url else None
    if classes:
        spectral_url, luminosity_url = _star_class_urls(star["star_type"], class_url)
        links = [_link(spectral_url, f"Spectral class {esc(classes[0])}") if spectral_url else "",
                 _link(luminosity_url, f"Luminosity class {esc(classes[1])}") if luminosity_url else ""]
    return _row_html(
        esc(star["name"]), stats, sections["stars"].get(str(star["id"])), children_html, children_visible=True,
        links=links, facilities=(by_host or {}).get(("star", star["id"])),
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


def _planet_row_html(body, sections, is_moon=False, class_url=None, by_host=None, species=None, admin_rows=None):
    """
    A planet's (or moon's) row: its class, one type chip, a "Habitable
    moon" chip when one of its moons is habitable, "Inhabited" when it is,
    and its distance, period and gravity. The zone is left to the body's
    own page. A planet's moons follow in their own "N moons" group
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
        # Not esc()'d: built from floats and fixed unit literals.
        _stat(format_body_distance(body["distance_km"], is_moon)),
        _stat(format_period(body["period_years"]) if body.get("period_years") is not None else ""),
        _stat(_gravity_text(body.get("gravity_g"))),
    ]
    section = sections["moons" if is_moon else "planets"]
    after_html = ""
    if moons:
        label = f'{len(moons)} moon{"s" if len(moons) != 1 else ""}'
        moon_rows = "".join(_planet_row_html(moon, sections, is_moon=True, class_url=class_url, by_host=by_host,
                                              admin_rows=admin_rows)
                            for moon in moons)
        after_html = (
            f'<details class="moon-group"><summary>{label} of {esc(body["name"])}</summary>'
            f'<ul class="system-list">{moon_rows}</ul></details>'
        )
    class_link = ""
    planet_url = class_url("planet", body["planet_class"]) if class_url and body["planet_class"] else None
    if planet_url:
        class_link = _link(planet_url, f'Planet class {esc(body["planet_class"])}')
    dominant = None if is_moon else (species or {}).get(body["id"])
    species_link = ""
    if dominant:
        stats.insert(2, f'<span class="stat">Species: {esc(dominant["name"])}</span>')
        species_link = f'Dominant species: {_link(dominant["url"], esc(dominant["name"]))}'
    return _row_html(esc(body["name"]), stats, section.get(str(body["id"])), after_html, children_visible=True,
                     links=[class_link, species_link], facilities=(by_host or {}).get(("moon" if is_moon else "planet", body["id"])),
                     admin=_row_admin_html(admin_rows, "moon" if is_moon else "planet", body["id"]))


BELT_TOP_MINERALS = 3
"""int: How many of a belt's components its row names (largest share first)."""


def _belt_minerals_text(composition):
    """A belt's top components as "Iron, nickel, olivine", or `""`."""
    names = [part["component"] for part in (composition or [])[:BELT_TOP_MINERALS]]
    return esc(", ".join(names).capitalize()) if names else ""


def _belt_row_html(belt, sections, by_host=None, admin_rows=None):
    """
    A belt's row: its density, its range (the nominal distance only when
    the range is missing) and its top minerals.
    """
    if belt.get("lower_limit_km") is not None and belt.get("upper_limit_km") is not None:
        where = f'{format_distance_km(belt["lower_limit_km"])} to {format_distance_km(belt["upper_limit_km"])}'
    else:
        where = format_distance_km(belt["distance_km"])
    stats = [
        _stat(esc(belt["density"]).capitalize() if belt.get("density") else ""),
        _stat(where),
        _stat(_belt_minerals_text(belt.get("composition"))),
    ]
    return _row_html("Asteroid Belt", stats, sections["belts"].get(str(belt["id"])),
                     facilities=(by_host or {}).get(("asteroid_belt", belt["id"])),
                     admin=_row_admin_html(admin_rows, "belt", belt["id"]))


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


def _comet_row_html(comet, sections, class_url=None):
    kind = "Elliptical" if comet["orbit_type"] == "elliptical" else "Parabolic"
    stats = [
        _stat(f"{kind} comet"),
        _stat("Active" if comet["is_active"] else "Dormant"),
        _stat(f'Perihelion {format_distance_km(comet["perihelion_distance_km"])}'),
        _stat(format_period(comet["orbital_period_years"]) if comet.get("orbital_period_years") is not None
              else "Single apparition"),
    ]
    # UX.29: every comet links its class, a parabolic one too.
    code = comet.get("period_class") if comet["orbit_type"] == "elliptical" else PARABOLIC_COMET_CLASS
    period_url = class_url("comet", code) if class_url and code else None
    links = [_link(period_url, esc(class_entry("comet", code)["name"]))] if period_url else []
    return _row_html(esc(comet["name"]), stats, sections["comets"].get(str(comet["id"])), links=links)


def _orbiting_rows_html(planets, belts, comets, sections, class_url=None, by_host=None, species=None, admin_rows=None):
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
            rows.append(_comet_row_html(pending.pop(0), sections, class_url))
        if kind == "planet":
            rows.append(_planet_row_html(body, sections, class_url=class_url, by_host=by_host, species=species,
                                         admin_rows=admin_rows))
        else:
            rows.append(_belt_row_html(body, sections, by_host, admin_rows))
    rows.extend(_comet_row_html(comet, sections, class_url) for comet in pending)
    return "".join(rows)


def system_list_html(system, sections, class_url=None, facilities=(), species=None, admin_rows=None):
    """
    The system rendered natively: the page's overview (a binary pair's
    own data, the system summary, any flavor text) above an expandable
    list of its stars, planets (with their moons nested under them),
    asteroid belts and comets. Each row shows compact stats and opens onto
    that body's own generated description -- see
    `planetgen.db.render.render_system_sections`.

    A `'wide'` (S-type) pair's bodies each orbit one of its two stars, so
    they're nested under that star's own row; a single star's or a
    `'close'` (P-type) pair's bodies orbit the whole system and follow
    the star rows at the top level.

    `facilities` (`GET /api/systems/<id>/facilities`'s items) are listed
    in their host's own row, and `species` (`{planet_id: {"name", "url"}}`)
    names each life world's dominant species in its row. `admin_rows`
    (`_row_admin_html`) gives each planet, moon and belt row its Admin menu.
    """
    stars, planets, belts, comets = system["stars"], system["planets"], system["belts"], system["comets"]
    by_host = _facilities_by_host(facilities)
    if system.get("binary_configuration") == "wide":
        rows = []
        for star in stars:
            children = _orbiting_rows_html(
                [p for p in planets if p["star_id"] == star["id"]],
                [b for b in belts if b["star_id"] == star["id"]],
                [c for c in comets if c["star_id"] == star["id"]],
                sections,
                class_url,
                by_host,
                species,
                admin_rows,
            )
            children_html = f'<ul class="system-list">{children}</ul>' if children else ""
            rows.append(_star_row_html(star, sections, children_html, class_url, by_host))
        rows_html = "".join(rows)
    else:
        # UX.53: the Stars table already shows these stars, so the list
        # starts at the first planet or belt (a star that hosts a facility
        # keeps its row, which is where the facility is listed).
        rows_html = "".join(_star_row_html(star, sections, class_url=class_url, by_host=by_host)
                            for star in stars if _star_has_list_row(system, star, by_host))
        rows_html += _orbiting_rows_html(planets, belts, comets, sections, class_url, by_host, species, admin_rows)

    overview_html = markdown_to_html(sections["overview"]) if sections["overview"] else ""
    return f"""
<div class="prose system-overview">{overview_html}</div>
<ul class="system-list system-list-root">{rows_html}</ul>
"""
