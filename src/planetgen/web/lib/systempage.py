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
literals, e.g. `<sup>`) and from `describe.html`, which escapes all of the
Document's text. So the page wraps what these return in `trusted_html`.

Class labels (a star's type, a planet's class) link to their class
reference pages (`/classes/...`) through an optional `class_url(type_slug,
code)` hook, the way `fmt`'s location helpers take a `system_url` hook:
it returns the page's URL or `None`, and without it the labels stay
plain text. A class's own letter or name in a row's `<summary>` is the link
(UX.182): there is no separate "Planet class M" line. The symbols and chips
beside it link into the habitability page the same way (UX.181).
"""

from planetgen.describe.html import blocks_html
from planetgen.physics.habitability_world import DOMAINS, TIER_NAMES
from planetgen.web.lib.classref import PARABOLIC_COMET_CLASS, class_entry, star_type_classes, star_type_split
from planetgen.web.lib.fmt import (
    distance_km_html, esc, format_gravity_g_short, format_number, measured_html, period_years_html, speed_kms_html,
)
from planetgen.web.lib.habitability_symbols import CHIP_ANCHORS, PHI6_CHIP_ANCHOR, equipment_anchor, help_link, symbol_row
from planetgen.web.lib.tabledisplay import format_star_radius


def _ref_url(ref):
    """The page fragment of the row a Document reference (`"planet:12"`) names."""
    kind, _, db_id = ref.partition(":")
    return f"#{body_anchor(kind, db_id)}" if kind in ("star", "planet", "moon", "belt", "comet") else None


def _document_html(blocks):
    """A row's Document blocks (`planetgen.describe.model.split_by_body`) as HTML, or `""` for none."""
    return blocks_html(blocks, _ref_url) if blocks else ""


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


def _star_type_html(star_type, class_url):
    """A star type as HTML (UX.182): with `class_url`, its spectral part ("G2") links to its spectral class page and its
    luminosity part ("V") to its luminosity class page, the words after them stay plain; unlinked where either page
    is missing or the type does not parse."""
    parts = star_type_split(star_type) if class_url else None
    if parts is None:
        return esc(star_type)
    spectral_url, luminosity_url = _star_class_urls(star_type, class_url)
    spectral, luminosity, rest = parts
    return f"{_link(spectral_url, esc(spectral))}{_link(luminosity_url, esc(luminosity))}{esc(rest)}"


def stars_html(stars, class_url=None, notes_html=""):
    """The Stars table. With `class_url` (see the module docstring), a
    star's type links to its spectral and luminosity class pages. `notes_html`
    (`star_notes_html`) goes under the table."""
    # UX.53: "single" says nothing; the Role column is for binaries and multiples.
    show_role = any(star["role"] != "single" for star in stars)
    rows = "".join(
        "<tr>"
        + (f'<td data-label="Role">{esc(star["role"])}</td>' if show_role else "")
        + f'<td class="stack-title">{esc(star["name"])}</td>'
        f'<td data-label="Type">{_star_type_html(star["star_type"], class_url)}</td>'
        # Not esc()'d: these are built entirely from floats and fixed unit
        # literals (never database TEXT), and legitimately contain a raw
        # `<sup>exponent</sup>` -- see `tabledisplay.py`.
        f'<td data-label="Mass">{measured_html("kg", star["mass_kg"])}</td>'
        f'<td data-label="Radius">{format_star_radius(star["radius_km"])}</td>'
        f'<td data-label="Temp">{format_number(star["temperature_k"])} K</td>'
        f'<td data-label="Luminosity">{measured_html("W", star["luminosity_w"])}</td>'
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
    no row any more (UX.53): each one's description, under the Stars table."""
    by_host = _facilities_by_host(facilities)
    notes = []
    for star in system["stars"]:
        if _star_has_list_row(system, star, by_host):
            continue
        blocks = sections["stars"].get(str(star["id"]))
        if not blocks:
            continue
        heading = f"<h3>{esc(star['name'])}</h3>" if len(system["stars"]) > 1 else ""
        notes.append(f'<div class="star-note prose">{heading}{_document_html(blocks)}</div>')
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
        "distance": distance_km_html(facility["orbit_distance_km"]) if orbital else None,
        "period": period_years_html(facility["orbit_period_years"])
        if orbital and facility.get("orbit_period_years") is not None else None,
        "speed": speed_kms_html(facility["orbital_speed_kms"])
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


def body_anchor(kind, body_id):
    """The fragment (without `#`) of a body's row on the system page (NAV.8): `planet-12`, `moon-7`,
    `star-3`, `belt-2`, `comet-9`. `static/bodyanchor.js` opens and highlights the row it names."""
    return f"{kind}-{body_id}"


def _row_html(title, stats, blocks, children_html="", children_visible=False, links=(), facilities=(), admin="",
              anchor=None):
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
    at the row's right end. `anchor` (`body_anchor`) is the row's `id`.
    """
    stats_html = "".join(stats) + _facility_stat(facilities)
    description = _facilities_html(facilities) + _document_html(blocks)
    links = [link for link in links if link]
    if links:
        description = f'<p class="class-links">{" &middot; ".join(links)}</p>{description}'
    inside, after = ("", children_html) if children_visible else (children_html, "")
    return f"""
<li{f' id="{anchor}"' if anchor else ""}{' class="has-admin"' if admin else ""}><details class="body-row">
<summary><span class="body-name">{title}</span><span class="body-stats">{stats_html}</span></summary>
<div class="body-detail prose">{description}</div>
{inside}
</details>{admin}{after}</li>"""


def _star_row_html(star, sections, children_html="", class_url=None, by_host=None):
    stats = [f'<span class="stat">{_star_type_html(star["star_type"], class_url)}</span>']
    if star["role"] != "single":
        stats.append(f'<span class="stat">{esc(star["role"].capitalize())}</span>')
    return _row_html(
        esc(star["name"]), stats, sections["stars"].get(str(star["id"])), children_html, children_visible=True,
        facilities=(by_host or {}).get(("star", star["id"])), anchor=body_anchor("star", star["id"]),
    )


def _flag_chip(text, help_url):
    """A "flag-yes" chip (UX.181), linked to its entry on the habitability page when `help_url` is given."""
    return help_link(f'<span class="flag flag-yes">{text}</span>', help_url, CHIP_ANCHORS[text.lower()])


def _habitable_chip(body, help_url=None):
    """The "Habitable class" chip (the class can carry life, which says nothing of its PHI-6 rating), or `""`. The
    row names the class and no Terrestrial or Gas Giant label (UX.32)."""
    return _flag_chip("Habitable class", help_url) if body["habitable"] else ""


def _stat(text):
    return f'<span class="stat">{text}</span>' if text else ""


def _habitability_chip(body, help_url=None):
    """GEN.89, GEN.655, UX.154: the equipment a human needs on this body as one chip per item, solid for a required item
    and outlined for a recommended one, each in the colour of the PHI-6 domain it comes from; a body that needs nothing
    shows one chip of its worst PHI-6 domain's colour; empty for a gas giant or a body not scored. The ratings
    themselves are the symbols beside it (UX.100). With `help_url` each chip links to the habitability page's entry for
    its item (UX.181)."""
    levels = body.get("equipment")
    if levels is None:
        return ""
    chips = []
    for item in levels["items"]:
        tier = max(body.get(f"tier_{domain.lower()}") or 0 for domain in item.get("domains") or [item["domain"]])
        text = item["item"][:1].upper() + item["item"][1:]
        hidden = ' <span class="sr-only">(recommended)</span>' if item["severity"] == "recommended" else ""
        chip = f'<span class="flag flag-eq flag-eq-{item["severity"]} flag-tier flag-tier-{tier}">{esc(text)}{hidden}</span>'
        chips.append(help_link(chip, help_url, equipment_anchor(item["item"])))
    if not chips:
        worst = max(body.get(f"tier_{domain}") or 0 for domain in DOMAINS)
        return help_link(f'<span class="flag flag-tier flag-tier-{worst}">PHI-6 {TIER_NAMES[worst]}</span>', help_url,
                         PHI6_CHIP_ANCHOR)
    return f'<span class="eq-chips">{"".join(chips)}</span>'


def _gravity_text(gravity_g):
    return format_gravity_g_short(gravity_g) if gravity_g is not None else ""


def _planet_row_html(body, sections, is_moon=False, class_url=None, by_host=None, species=None, admin_rows=None,
                     sprite_url=None, help_url=None):
    """
    A planet's (or moon's) row: its class, a "Habitable class" chip for a
    habitable class, a "Habitable-class moon" chip when one of its moons is
    habitable, "Inhabited" when it is, its moon count, its distance, period
    and gravity, and the equipment chip and the PHI-6 symbols (UX.100; none
    where `sprite_url` is not given or the body is unscored). The class letter
    links to its class page (UX.182); the chips and symbols link into the
    habitability page at `help_url` (UX.181). The zone is
    left to the body's own page. A planet's moons follow in their own "N
    moons" group right under its row, collapsed, needing no script.
    """
    moons = [] if is_moon else (body.get("moons") or [])
    class_code = body["planet_class"]
    class_text = _link(class_url("planet", class_code) if class_url and class_code else None, esc(class_code))
    stats = [
        _stat(f"Class {class_text}" if class_code else ""),
        _habitable_chip(body, help_url),
    ]
    if any(moon["habitable"] for moon in moons):
        stats.append(_flag_chip("Habitable-class moon", help_url))
    if body["inhabited"]:
        stats.append(_flag_chip("Inhabited", help_url))
    if moons:
        stats.append(_stat(f'{len(moons)} moon{"s" if len(moons) != 1 else ""}'))
    stats += [
        # Not esc()'d: built from floats and fixed unit literals.
        _stat(distance_km_html(body["distance_km"])),
        _stat(period_years_html(body["period_years"]) if body.get("period_years") is not None else ""),
        _stat(_gravity_text(body.get("gravity_g"))),
        _habitability_chip(body, help_url),
        symbol_row(body, sprite_url, help_url=help_url) if sprite_url else "",
    ]
    section = sections["moons" if is_moon else "planets"]
    after_html = ""
    if moons:
        label = f'{len(moons)} moon{"s" if len(moons) != 1 else ""}'
        moon_rows = "".join(_planet_row_html(moon, sections, is_moon=True, class_url=class_url, by_host=by_host,
                                              admin_rows=admin_rows, sprite_url=sprite_url, help_url=help_url)
                            for moon in moons)
        after_html += (
            f'<details class="moon-group"><summary>{label} of {esc(body["name"])}</summary>'
            f'<ul class="system-list">{moon_rows}</ul></details>'
        )
    dominant = None if is_moon else (species or {}).get(body["id"])
    species_link = ""
    if dominant:
        stats.insert(2, f'<span class="stat">Species: {esc(dominant["name"])}</span>')
        species_link = f'Dominant species: {_link(dominant["url"], esc(dominant["name"]))}'
    return _row_html(esc(body["name"]), stats, section.get(str(body["id"])), after_html, children_visible=True,
                     links=[species_link], facilities=(by_host or {}).get(("moon" if is_moon else "planet", body["id"])),
                     admin=_row_admin_html(admin_rows, "moon" if is_moon else "planet", body["id"]),
                     anchor=body_anchor("moon" if is_moon else "planet", body["id"]))


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
        where = f'{distance_km_html(belt["lower_limit_km"])} to {distance_km_html(belt["upper_limit_km"])}'
    else:
        where = distance_km_html(belt["distance_km"])
    stats = [
        _stat(esc(belt["density"]).capitalize() if belt.get("density") else ""),
        _stat(where),
        _stat(_belt_minerals_text(belt.get("composition"))),
    ]
    return _row_html("Asteroid Belt", stats, sections["belts"].get(str(belt["id"])),
                     facilities=(by_host or {}).get(("asteroid_belt", belt["id"])),
                     admin=_row_admin_html(admin_rows, "belt", belt["id"]), anchor=body_anchor("belt", belt["id"]))


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
    # UX.29: every comet has a class, a parabolic one too; UX.182: its name is the link.
    code = comet.get("period_class") if comet["orbit_type"] == "elliptical" else PARABOLIC_COMET_CLASS
    period_url = class_url("comet", code) if class_url and code else None
    stats = [
        _stat(_link(period_url, esc(class_entry("comet", code)["name"])) if period_url else f"{kind} comet"),
        _stat("Active" if comet["is_active"] else "Dormant"),
        _stat(f'Perihelion {distance_km_html(comet["perihelion_distance_km"])}'),
        _stat(period_years_html(comet["orbital_period_years"]) if comet.get("orbital_period_years") is not None
              else "Single apparition"),
    ]
    return _row_html(esc(comet["name"]), stats, sections["comets"].get(str(comet["id"])),
                     anchor=body_anchor("comet", comet["id"]))


def _orbiting_rows_html(planets, belts, comets, sections, class_url=None, by_host=None, species=None, admin_rows=None,
                        sprite_url=None, help_url=None):
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
                                         admin_rows=admin_rows, sprite_url=sprite_url, help_url=help_url))
        else:
            rows.append(_belt_row_html(body, sections, by_host, admin_rows))
    rows.extend(_comet_row_html(comet, sections, class_url) for comet in pending)
    return "".join(rows)


def system_list_html(system, sections, class_url=None, facilities=(), species=None, admin_rows=None,
                     habitability_url=None, sprite_url=None):
    """
    The system rendered natively: the page's overview (a binary pair's
    own data, the system summary, any flavor text) above an expandable
    list of its stars, planets (with their moons nested under them),
    asteroid belts and comets. Each row shows compact stats and opens onto
    that body's own generated description -- see
    `planetgen.describe.model.split_by_body`.

    A `'wide'` (S-type) pair's bodies each orbit one of its two stars, so
    they're nested under that star's own row; a single star's or a
    `'close'` (P-type) pair's bodies orbit the whole system and follow
    the star rows at the top level.

    `facilities` (`GET /api/systems/<id>/facilities`'s items) are listed
    in their host's own row, and `species` (`{planet_id: {"name", "url"}}`)
    names each life world's dominant species in its row. `admin_rows`
    (`_row_admin_html`) gives each planet, moon and belt row its Admin menu.
    `habitability_url` (UX.181) is the habitability page: each chip and PHI-6
    symbol in a row links to its section there. `sprite_url` (the icon
    sprite, UX.100) draws each scored planet's and moon's PHI-6 symbols in its row.
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
                sprite_url,
                habitability_url,
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
        rows_html += _orbiting_rows_html(planets, belts, comets, sections, class_url, by_host, species, admin_rows,
                                         sprite_url, habitability_url)

    overview_html = _document_html(sections["overview"])
    return f"""
<div class="prose system-overview">{overview_html}</div>
<ul class="system-list system-list-root">{rows_html}</ul>
"""
