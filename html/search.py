#!/usr/bin/env python3
# html/search.py

"""
Cross-database search.

Two complementary ways in, per one chosen database (`?db=`, same
convention as every other page here):

  1. Click-to-filter attribute tags -- object type, star spectral/
     luminosity class, planet class/body type, and supported life
     chemistry -- built from exactly the *distinct values actually
     present* in the chosen database (an empty facet, e.g. no binary
     stars generated, simply renders no button for it). Clicking a tag
     toggles it on/off via a plain link that rewrites the query string,
     so this works with JavaScript disabled, same as every other page in
     `html/`.
  2. A name search, one field per nameable entity (sector, star system,
     star, planet/moon), each with its own HTML5 `<datalist>` for
     autocomplete -- no JavaScript, just the browser's native
     suggestion UI, populated from that entity's own distinct names.
     Asteroid belts have no name of their own (see `schema.sql`'s own
     comment on this), so the only way to reach them here is the
     "Asteroid Belt" object-type tag.

Tag facets and name fields combine as independent AND'd filters *within*
the object type they apply to (e.g. a star spectral tag plus a star name
term both narrow the Stars panel); a result panel for a given object type
is only shown when there's a specific reason to query it -- one of its
own tags/name field is active, or its type is explicitly selected via the
"Object Type" tag group. Picking an explicit Object Type tag acts as a
master filter: e.g. selecting only "Stars" hides the Planets panel even
if a planet-class tag happens to also be selected.
"""

import os
import sys
from urllib.parse import parse_qs, urlencode

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))
# Falls back to the project root (html/'s parent) so `stellarObjects` is
# importable even when it hasn't been `pip install`-ed system-wide -- true
# for the default deployment layout (`html/` and `stellarObjects/` as
# siblings under /var/lib/planetGen).
sys.path.append(os.path.dirname(_HTML_DIR))

from dbutil import esc, fetch_all, fetch_one, open_readonly, resolve_db_path
from page import run

try:
    from stellarObjects.physical_constants import SPECTRAL_CLASS_COLORS
    from stellarObjects.program_constants import PLANET_CLASSES
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # tag labels fall back to their raw codes rather than failing outright.
    SPECTRAL_CLASS_COLORS = {}
    PLANET_CLASSES = {}

RESULT_LIMIT = 300
AUTOCOMPLETE_LIMIT = 500

TAG_FACETS = ("type", "spectral", "luminosity", "class", "body", "life", "density")

# starData.py's own Yerkes-class-to-descriptive-label mapping (Star.__init__),
# duplicated here (not imported) since it's a plain literal there, not an
# importable constant. "D" is an alternate white-dwarf code seen elsewhere
# in starData.py alongside "VII".
YERKES_LABELS = {
    "0": "Hypergiant",
    "IA": "Supergiant",
    "IAB": "Intermediate-size Luminous Supergiant",
    "IB": "Less Luminous Supergiant",
    "II": "Bright Giant",
    "III": "Giant",
    "IV": "Subgiant",
    "V": "Main Sequence",
    "VII": "White Dwarf",
    "D": "White Dwarf",
}
YERKES_ORDER = ["0", "IA", "IAB", "IB", "II", "III", "IV", "V", "VII", "D"]

BODY_LABELS = {"t": "Terrestrial", "g": "Gas Giant"}


def _like_pattern(term):
    """Escapes `%`/`_`/`\\` in a user-supplied substring so it's safe to
    use as a SQL LIKE pattern (paired with `ESCAPE '\\'` in the query)."""
    escaped = term.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")
    return f"%{escaped}%"


def _count_suffix(n, truncated):
    return f" ({n}{'+' if truncated else ''})"


def _truncated_note(truncated):
    if not truncated:
        return ""
    return f'<p class="hint">Showing the first {RESULT_LIMIT} matches -- refine your search for more precise results.</p>'


def _sector_link(db_name, sector_id):
    if sector_id is None:
        return "Standalone"
    return f'<a href="sector.py?db={esc(db_name)}&id={sector_id}">View</a>'


def _system_star_summary(conn, star_system_id, is_binary):
    """One-line star description for a system row, same logic as
    browse.py's own helper of the same purpose (small enough, and
    specific enough to this page's query shape, not to share a module
    over)."""
    if is_binary:
        row = fetch_one(conn, "SELECT binary_type FROM star_systems WHERE id = ?", (star_system_id,))
        return row["binary_type"] if row else ""
    row = fetch_one(
        conn,
        "SELECT star_type FROM stars WHERE star_system_id = ? AND role = 'single'",
        (star_system_id,),
    )
    return row["star_type"] if row else ""


# ---------------------------------------------------------------------
# Facet option discovery -- each returns a list of
# (value, label, count, tooltip_or_None) tuples, one per distinct value
# actually present in the database (never a fixed/static enumeration).
# ---------------------------------------------------------------------

def _facet_options_type(conn):
    star_c = fetch_one(conn, "SELECT COUNT(*) AS c FROM stars")["c"]
    planet_c = fetch_one(conn, "SELECT COUNT(*) AS c FROM planets WHERE is_moon = 0")["c"]
    moon_c = fetch_one(conn, "SELECT COUNT(*) AS c FROM planets WHERE is_moon = 1")["c"]
    belt_c = fetch_one(conn, "SELECT COUNT(*) AS c FROM asteroid_belts")["c"]
    opts = []
    for value, label, count in (
        ("star", "Stars", star_c),
        ("planet", "Planets", planet_c),
        ("moon", "Moons", moon_c),
        ("belt", "Asteroid Belts", belt_c),
    ):
        if count:
            opts.append((value, label, count, None))
    return opts


def _facet_options_spectral(conn):
    rows = fetch_all(conn, "SELECT SUBSTR(star_type, 1, 1) AS v, COUNT(*) AS c FROM stars GROUP BY v ORDER BY v")
    opts = []
    for row in rows:
        letter = row["v"]
        color = SPECTRAL_CLASS_COLORS.get(letter)
        tip = f"{color} star" if color else None
        opts.append((letter, f"{letter}-Type Star", row["c"], tip))
    return opts


def _facet_options_luminosity(conn):
    rows = fetch_all(
        conn,
        "SELECT yerkes_class AS v, COUNT(*) AS c FROM stars WHERE yerkes_class IS NOT NULL GROUP BY v",
    )
    rows = sorted(rows, key=lambda r: YERKES_ORDER.index(r["v"]) if r["v"] in YERKES_ORDER else len(YERKES_ORDER))
    return [
        (row["v"], YERKES_LABELS.get(row["v"], row["v"]), row["c"], f"Yerkes class {row['v']}")
        for row in rows
    ]


def _facet_options_class(conn):
    rows = fetch_all(
        conn,
        "SELECT planet_class AS v, COUNT(*) AS c FROM planets WHERE planet_class IS NOT NULL GROUP BY v ORDER BY v",
    )
    opts = []
    for row in rows:
        info = PLANET_CLASSES.get(row["v"], {})
        opts.append((row["v"], f"Class {row['v']} Planet", row["c"], info.get("description")))
    return opts


def _facet_options_body(conn):
    rows = fetch_all(conn, "SELECT body_type AS v, COUNT(*) AS c FROM planets GROUP BY v ORDER BY v")
    return [(row["v"], BODY_LABELS.get(row["v"], row["v"]), row["c"], None) for row in rows]


def _facet_options_life(conn):
    rows = fetch_all(
        conn,
        "SELECT life_chemical AS v, COUNT(*) AS c FROM planets WHERE life_chemical IS NOT NULL GROUP BY v ORDER BY v",
    )
    return [(row["v"], row["v"], row["c"], None) for row in rows]


def _facet_options_density(conn):
    rows = fetch_all(conn, "SELECT density AS v, COUNT(*) AS c FROM asteroid_belts GROUP BY v ORDER BY v")
    return [(row["v"], row["v"].capitalize(), row["c"], None) for row in rows]


def _name_list(conn, table, limit=AUTOCOMPLETE_LIMIT):
    rows = fetch_all(conn, f"SELECT DISTINCT name FROM {table} ORDER BY name LIMIT ?", (limit,))
    return [row["name"] for row in rows]


def _datalist_html(list_id, names):
    options = "".join(f'<option value="{esc(n)}"></option>' for n in names)
    return f'<datalist id="{list_id}">{options}</datalist>'


# ---------------------------------------------------------------------
# Query-string state and link building -- every tag button and "remove
# filter" chip is a plain <a> that re-renders this same page with one
# value toggled, preserving every other currently active filter.
# ---------------------------------------------------------------------

def _build_url(db_name, texts, tags):
    params = [("db", db_name)]
    for key in ("sector_q", "system_q", "star_q", "planet_q"):
        value = texts.get(key)
        if value:
            params.append((key, value))
    for facet in TAG_FACETS:
        for value in sorted(tags.get(facet, ())):
            params.append((facet, value))
    # Query-string "&" separators are ours (not user data), but still
    # need entity-escaping to be strictly valid inside an href attribute.
    return "search.py?" + urlencode(params).replace("&", "&amp;")


def _toggle_url(state, facet, value):
    new_tags = {key: set(values) for key, values in state["tags"].items()}
    values = new_tags.setdefault(facet, set())
    if value in values:
        values.discard(value)
    else:
        values.add(value)
    return _build_url(state["db"], state["texts"], new_tags)


def _remove_text_url(state, key):
    new_texts = dict(state["texts"])
    new_texts[key] = ""
    return _build_url(state["db"], new_texts, state["tags"])


def _tag_group_html(title, facet, options, selected, state):
    if not options:
        return ""
    buttons = []
    for value, label, count, tip in options:
        css_class = "tag active" if value in selected else "tag"
        title_attr = f' title="{esc(tip)}"' if tip else ""
        buttons.append(
            f'<a class="{css_class}" href="{_toggle_url(state, facet, value)}"{title_attr}>'
            f'{esc(label)} <span class="tag-count">{count}</span></a>'
        )
    return f"""
<div class="tag-group">
<h3>{esc(title)}</h3>
<div class="tag-list">{''.join(buttons)}</div>
</div>
"""


def _active_filters_html(state, facet_labels):
    chips = []
    text_labels = {"sector_q": "Sector", "system_q": "System", "star_q": "Star", "planet_q": "Planet"}
    for key, label in text_labels.items():
        value = state["texts"].get(key)
        if value:
            chips.append((f'{esc(label)}: &ldquo;{esc(value)}&rdquo;', _remove_text_url(state, key)))
    for facet in TAG_FACETS:
        for value in sorted(state["tags"].get(facet, ())):
            label = facet_labels.get((facet, value), value)
            chips.append((esc(label), _toggle_url(state, facet, value)))
    if not chips:
        return ""
    items = "".join(
        f'<span class="filter-chip">{text}<a href="{url}" aria-label="Remove filter">&times;</a></span>'
        for text, url in chips
    )
    return f'<div class="active-filters">{items}</div>'


# ---------------------------------------------------------------------
# Result panels
# ---------------------------------------------------------------------

def _sectors_panel(conn, db_name, term):
    rows = fetch_all(
        conn,
        "SELECT id, name, edge_mpc FROM sectors WHERE name LIKE ? ESCAPE '\\' ORDER BY name LIMIT ?",
        (_like_pattern(term), RESULT_LIMIT + 1),
    )
    truncated = len(rows) > RESULT_LIMIT
    rows = rows[:RESULT_LIMIT]
    body_rows = "".join(
        "<tr>"
        f'<td><a href="sector.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
        f'<td>{row["edge_mpc"]:,.2f} mpc</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="2"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Sectors{_count_suffix(len(rows), truncated)}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Cube Edge</th></tr></thead>
  <tbody>{body_rows}</tbody>
</table></div>
{_truncated_note(truncated)}
</section>
"""


def _systems_panel(conn, db_name, term):
    rows = fetch_all(
        conn,
        """
        SELECT id, name, sector_id, is_binary
        FROM star_systems
        WHERE name LIKE ? ESCAPE '\\'
        ORDER BY name
        LIMIT ?
        """,
        (_like_pattern(term), RESULT_LIMIT + 1),
    )
    truncated = len(rows) > RESULT_LIMIT
    rows = rows[:RESULT_LIMIT]
    body_rows = "".join(
        "<tr>"
        f'<td><a href="system.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
        f'<td>{_sector_link(db_name, row["sector_id"])}</td>'
        f'<td>{"Yes" if row["is_binary"] else "No"}</td>'
        f'<td>{esc(_system_star_summary(conn, row["id"], row["is_binary"]))}</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="4"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Systems{_count_suffix(len(rows), truncated)}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Sector</th><th>Binary</th><th>Star type</th></tr></thead>
  <tbody>{body_rows}</tbody>
</table></div>
{_truncated_note(truncated)}
</section>
"""


def _stars_panel(conn, db_name, spectral_tags, luminosity_tags, term):
    clauses, params = [], []
    if spectral_tags:
        clauses.append(f"SUBSTR(s.star_type, 1, 1) IN ({','.join('?' * len(spectral_tags))})")
        params.extend(sorted(spectral_tags))
    if luminosity_tags:
        clauses.append(f"s.yerkes_class IN ({','.join('?' * len(luminosity_tags))})")
        params.extend(sorted(luminosity_tags))
    if term:
        clauses.append("s.name LIKE ? ESCAPE '\\'")
        params.append(_like_pattern(term))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    params.append(RESULT_LIMIT + 1)
    rows = fetch_all(
        conn,
        f"""
        SELECT s.name, s.role, s.star_type, s.star_system_id, ss.name AS system_name, ss.sector_id
        FROM stars s
        JOIN star_systems ss ON ss.id = s.star_system_id
        WHERE 1=1{where}
        ORDER BY ss.name, s.name
        LIMIT ?
        """,
        params,
    )
    truncated = len(rows) > RESULT_LIMIT
    rows = rows[:RESULT_LIMIT]
    body_rows = "".join(
        "<tr>"
        f'<td>{esc(row["name"])}</td>'
        f'<td>{esc(row["role"])}</td>'
        f'<td>{esc(row["star_type"])}</td>'
        f'<td><a href="system.py?db={esc(db_name)}&id={row["star_system_id"]}">{esc(row["system_name"])}</a></td>'
        f'<td>{_sector_link(db_name, row["sector_id"])}</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="5"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Stars{_count_suffix(len(rows), truncated)}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Role</th><th>Type</th><th>System</th><th>Sector</th></tr></thead>
  <tbody>{body_rows}</tbody>
</table></div>
{_truncated_note(truncated)}
</section>
"""


def _planets_panel(conn, db_name, class_tags, body_tags, life_tags, term, is_moon_filter):
    clauses, params = [], []
    if is_moon_filter is not None:
        clauses.append("p.is_moon = ?")
        params.append(is_moon_filter)
    if class_tags:
        clauses.append(f"p.planet_class IN ({','.join('?' * len(class_tags))})")
        params.extend(sorted(class_tags))
    if body_tags:
        clauses.append(f"p.body_type IN ({','.join('?' * len(body_tags))})")
        params.extend(sorted(body_tags))
    if life_tags:
        clauses.append(f"p.life_chemical IN ({','.join('?' * len(life_tags))})")
        params.extend(sorted(life_tags))
    if term:
        clauses.append("p.name LIKE ? ESCAPE '\\'")
        params.append(_like_pattern(term))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    params.append(RESULT_LIMIT + 1)
    rows = fetch_all(
        conn,
        f"""
        SELECT p.name, p.planet_class, p.body_type, p.life_chemical, p.is_moon,
               p.star_system_id, ss.name AS system_name, ss.sector_id
        FROM planets p
        JOIN star_systems ss ON ss.id = p.star_system_id
        WHERE 1=1{where}
        ORDER BY ss.name, p.orbital_index
        LIMIT ?
        """,
        params,
    )
    truncated = len(rows) > RESULT_LIMIT
    rows = rows[:RESULT_LIMIT]
    body_rows = "".join(
        "<tr>"
        f'<td>{esc(row["name"])}</td>'
        f'<td>{"Moon" if row["is_moon"] else "Planet"}</td>'
        f'<td>{esc(row["planet_class"]) or "&mdash;"}</td>'
        f'<td>{BODY_LABELS.get(row["body_type"], esc(row["body_type"]))}</td>'
        f'<td>{esc(row["life_chemical"]) or "&mdash;"}</td>'
        f'<td><a href="system.py?db={esc(db_name)}&id={row["star_system_id"]}">{esc(row["system_name"])}</a></td>'
        f'<td>{_sector_link(db_name, row["sector_id"])}</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="7"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Planets &amp; Moons{_count_suffix(len(rows), truncated)}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Type</th><th>Class</th><th>Body</th><th>Life Chemistry</th><th>System</th><th>Sector</th></tr></thead>
  <tbody>{body_rows}</tbody>
</table></div>
{_truncated_note(truncated)}
</section>
"""


def _belts_panel(conn, db_name, density_tags):
    clauses, params = [], []
    if density_tags:
        clauses.append(f"ab.density IN ({','.join('?' * len(density_tags))})")
        params.extend(sorted(density_tags))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    params.append(RESULT_LIMIT + 1)
    rows = fetch_all(
        conn,
        f"""
        SELECT ab.density, ab.composition_summary, ab.star_system_id, ss.name AS system_name, ss.sector_id
        FROM asteroid_belts ab
        JOIN star_systems ss ON ss.id = ab.star_system_id
        WHERE 1=1{where}
        ORDER BY ss.name, ab.orbital_index
        LIMIT ?
        """,
        params,
    )
    truncated = len(rows) > RESULT_LIMIT
    rows = rows[:RESULT_LIMIT]
    body_rows = "".join(
        "<tr>"
        f'<td>{esc(row["density"]).capitalize()}</td>'
        f'<td>{esc(row["composition_summary"])}</td>'
        f'<td><a href="system.py?db={esc(db_name)}&id={row["star_system_id"]}">{esc(row["system_name"])}</a></td>'
        f'<td>{_sector_link(db_name, row["sector_id"])}</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="4"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Asteroid Belts{_count_suffix(len(rows), truncated)}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Density</th><th>Composition</th><th>System</th><th>Sector</th></tr></thead>
  <tbody>{body_rows}</tbody>
</table></div>
{_truncated_note(truncated)}
</section>
"""


def _render_results(conn, db_name, texts, tags):
    type_tags = tags["type"]
    spectral_tags, luminosity_tags = tags["spectral"], tags["luminosity"]
    class_tags, body_tags, life_tags = tags["class"], tags["body"], tags["life"]
    density_tags = tags["density"]

    star_has_reason = bool(spectral_tags or luminosity_tags or texts["star_q"])
    planet_has_reason = bool(class_tags or body_tags or life_tags or texts["planet_q"])
    belt_has_reason = bool(density_tags)

    if type_tags:
        # An explicit Object Type selection is a master filter: it wins
        # over an attribute tag from another object type that also
        # happens to be selected (see the module docstring).
        stars_included = "star" in type_tags
        planets_included = ("planet" in type_tags) or ("moon" in type_tags)
        belts_included = "belt" in type_tags
    else:
        stars_included = star_has_reason
        planets_included = planet_has_reason
        belts_included = belt_has_reason

    sections = []
    if texts["sector_q"]:
        sections.append(_sectors_panel(conn, db_name, texts["sector_q"]))
    if texts["system_q"]:
        sections.append(_systems_panel(conn, db_name, texts["system_q"]))
    if stars_included:
        sections.append(_stars_panel(conn, db_name, spectral_tags, luminosity_tags, texts["star_q"]))
    if planets_included:
        is_moon_filter = None
        if type_tags:
            has_planet, has_moon = "planet" in type_tags, "moon" in type_tags
            if has_planet and not has_moon:
                is_moon_filter = 0
            elif has_moon and not has_planet:
                is_moon_filter = 1
        sections.append(
            _planets_panel(conn, db_name, class_tags, body_tags, life_tags, texts["planet_q"], is_moon_filter)
        )
    if belts_included:
        sections.append(_belts_panel(conn, db_name, density_tags))

    return "".join(sections) if sections else '<p class="hint">No matching objects.</p>'


def handler():
    raw = parse_qs(os.environ.get("QUERY_STRING", ""))

    def _first(key):
        values = raw.get(key)
        return values[0].strip() if values else ""

    db_name = _first("db")
    path = resolve_db_path(db_name)
    conn = open_readonly(path)
    try:
        texts = {
            "sector_q": _first("sector_q"),
            "system_q": _first("system_q"),
            "star_q": _first("star_q"),
            "planet_q": _first("planet_q"),
        }
        tags = {facet: {v.strip() for v in raw.get(facet, []) if v.strip()} for facet in TAG_FACETS}
        tags["type"] &= {"star", "planet", "moon", "belt"}
        tags["body"] &= {"t", "g"}

        state = {"db": db_name, "texts": texts, "tags": tags}

        facet_defs = (
            ("Object Type", "type", _facet_options_type(conn)),
            ("Star: Spectral Class", "spectral", _facet_options_spectral(conn)),
            ("Star: Luminosity Class", "luminosity", _facet_options_luminosity(conn)),
            ("Planet Class", "class", _facet_options_class(conn)),
            ("Planet Body Type", "body", _facet_options_body(conn)),
            ("Supported Life Chemistry", "life", _facet_options_life(conn)),
            ("Asteroid Belt Density", "density", _facet_options_density(conn)),
        )

        facet_labels = {
            (facet, value): label
            for _, facet, options in facet_defs
            for value, label, _count, _tip in options
        }

        tag_browser_html = "".join(
            _tag_group_html(title, facet, options, tags[facet], state) for title, facet, options in facet_defs
        )

        hidden_tag_inputs = "".join(
            f'<input type="hidden" name="{facet}" value="{esc(v)}">'
            for facet in TAG_FACETS
            for v in sorted(tags[facet])
        )
        datalists = "".join(
            [
                _datalist_html("dl-sector", _name_list(conn, "sectors")),
                _datalist_html("dl-system", _name_list(conn, "star_systems")),
                _datalist_html("dl-star", _name_list(conn, "stars")),
                _datalist_html("dl-planet", _name_list(conn, "planets")),
            ]
        )

        form_html = f"""
<form method="get" action="search.py" class="search-form">
  <input type="hidden" name="db" value="{esc(db_name)}">
  {hidden_tag_inputs}
  <div class="search-fields">
    <label class="search-field">Sector name
      <input type="text" name="sector_q" value="{esc(texts['sector_q'])}" list="dl-sector" autocomplete="off" placeholder="e.g. Voranthis Kelmoor">
    </label>
    <label class="search-field">System name
      <input type="text" name="system_q" value="{esc(texts['system_q'])}" list="dl-system" autocomplete="off" placeholder="e.g. Kepler-42">
    </label>
    <label class="search-field">Star name
      <input type="text" name="star_q" value="{esc(texts['star_q'])}" list="dl-star" autocomplete="off" placeholder="e.g. Kepler-42 A">
    </label>
    <label class="search-field">Planet / moon name
      <input type="text" name="planet_q" value="{esc(texts['planet_q'])}" list="dl-planet" autocomplete="off" placeholder="e.g. Kepler-42 b">
    </label>
  </div>
  {datalists}
  <div class="search-actions">
    <button type="submit" class="btn">Search</button>
    <a href="search.py?db={esc(db_name)}">Clear all filters</a>
  </div>
</form>
"""

        active_filters_html = _active_filters_html(state, facet_labels)

        any_active = bool(
            any(tags[facet] for facet in TAG_FACETS) or any(texts.values())
        )
        if any_active:
            results_html = _render_results(conn, db_name, texts, tags)
        else:
            results_html = '<p class="hint">Select a tag below, or enter a name above and press Search, to see matching results.</p>'
    finally:
        conn.close()

    breadcrumb = (
        '<p class="breadcrumb"><a href="index.py">Databases</a> &rarr; '
        f'<a href="browse.py?db={esc(db_name)}">{esc(db_name)}</a> &rarr; Search</p>'
    )

    body = f"""
{breadcrumb}
<section class="panel">
<h2>Search</h2>
{form_html}
{active_filters_html}
</section>
<section class="panel">
<h2>Browse by Tag</h2>
{tag_browser_html}
</section>
{results_html}
"""
    return f"Search: {db_name}", body


run(handler)
