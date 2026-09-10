#!/usr/bin/env python3
# html/search.py

"""
Cross-database search.

Two complementary ways in, per one chosen database (`?db=`, same
convention as every other page here):

  1. Click-to-filter attribute tags -- object type, star spectral/
     luminosity class, planet class/body type/life chemistry, and the
     same three for moons -- built from exactly the *distinct values
     actually present* in the chosen database (an empty facet, e.g. no
     binary stars generated, simply renders no button for it). Planets
     and moons live in separate tables (`planets`/`moons`, schema v2 --
     see `schema.sql`'s "v2" header note) with their own separate facets,
     so a "Class D" tag only ever means one or the other, never both at
     once the way a shared `is_moon`-flagged table would. Clicking a tag
     toggles it on/off via a plain link that rewrites the query string,
     so this works with JavaScript disabled, same as every other page in
     `html/`.
  2. A name search, one field per nameable entity (sector, star system,
     star, planet, moon), each with its own HTML5 `<datalist>` for
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

This page is a thin renderer over `GET /api/search` -- every query
(facet-option discovery, name search, autocomplete) lives once, in
`queryDb.search` (shared with the API itself), not duplicated here.
"""

import os
import sys
from urllib.parse import parse_qs, urlencode

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import get_search
from fmt import esc

from page import run

# Mirrors queryDb.SEARCH_TAG_FACETS -- duplicated, not imported, since
# this page only ever talks to the database through the API (see the
# module docstring); it needs to know the set of valid facet names to
# parse the incoming query string before it has a response to derive
# them from.
TAG_FACETS = (
    "type", "spectral", "luminosity",
    "class", "body", "life",
    "moon_class", "moon_body", "moon_life",
    "density",
)


def _count_suffix(n, truncated):
    return f" ({n}{'+' if truncated else ''})"


def _truncated_note(result_limit, truncated):
    if not truncated:
        return ""
    return f'<p class="hint">Showing the first {result_limit} matches -- refine your search for more precise results.</p>'


def _sector_link(db_name, sector_id):
    if sector_id is None:
        return "Standalone"
    return f'<a href="sector.py?db={esc(db_name)}&id={sector_id}">View</a>'


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
    for key in ("sector_q", "system_q", "star_q", "planet_q", "moon_q"):
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
    for option in options:
        value, label, count, tip = option["value"], option["label"], option["count"], option["tooltip"]
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
    text_labels = {
        "sector_q": "Sector", "system_q": "System", "star_q": "Star",
        "planet_q": "Planet", "moon_q": "Moon",
    }
    for key, label in text_labels.items():
        value = state["texts"].get(key)
        if value:
            chips.append((f'{esc(label)}: &ldquo;{esc(value)}&rdquo;', _remove_text_url(state, key)))
    for facet in TAG_FACETS:
        for value in sorted(state["tags"].get(facet, ())):
            label = facet_labels.get(f"{facet}:{value}", value)
            chips.append((esc(label), _toggle_url(state, facet, value)))
    if not chips:
        return ""
    items = "".join(
        f'<span class="filter-chip">{text}<a href="{url}" aria-label="Remove filter">&times;</a></span>'
        for text, url in chips
    )
    return f'<div class="active-filters">{items}</div>'


# ---------------------------------------------------------------------
# Result panels -- each renders one `results[panel]` entry
# (`{"rows": [...], "truncated": bool}`) from GET /api/search.
# ---------------------------------------------------------------------

def _sectors_panel(db_name, result):
    rows = result["rows"]
    body_rows = "".join(
        "<tr>"
        f'<td><a href="sector.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
        f'<td>{row["edge_mpc"]:,.2f} mpc</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="2"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Sectors{_count_suffix(len(rows), result["truncated"])}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Cube Edge</th></tr></thead>
  <tbody>{body_rows}</tbody>
</table></div>
{_truncated_note(len(rows), result["truncated"])}
</section>
"""


def _systems_panel(db_name, result):
    rows = result["rows"]
    body_rows = "".join(
        "<tr>"
        f'<td><a href="system.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
        f'<td>{_sector_link(db_name, row["sector_id"])}</td>'
        f'<td>{"Yes" if row["is_binary"] else "No"}</td>'
        f'<td>{esc(row["star_summary"])}</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="4"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Systems{_count_suffix(len(rows), result["truncated"])}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Sector</th><th>Binary</th><th>Star type</th></tr></thead>
  <tbody>{body_rows}</tbody>
</table></div>
{_truncated_note(len(rows), result["truncated"])}
</section>
"""


def _stars_panel(db_name, result):
    rows = result["rows"]
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
<h2>Stars{_count_suffix(len(rows), result["truncated"])}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Role</th><th>Type</th><th>System</th><th>Sector</th></tr></thead>
  <tbody>{body_rows}</tbody>
</table></div>
{_truncated_note(len(rows), result["truncated"])}
</section>
"""


def _planets_panel(db_name, result):
    rows = result["rows"]
    body_rows = "".join(
        "<tr>"
        f'<td>{esc(row["name"])}</td>'
        f'<td>{esc(row["planet_class"]) or "&mdash;"}</td>'
        f'<td>{"Gas Giant" if row["body_type"] == "g" else "Terrestrial"}</td>'
        f'<td>{esc(row["life_chemical"]) or "&mdash;"}</td>'
        f'<td><a href="system.py?db={esc(db_name)}&id={row["star_system_id"]}">{esc(row["system_name"])}</a></td>'
        f'<td>{_sector_link(db_name, row["sector_id"])}</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="6"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Planets{_count_suffix(len(rows), result["truncated"])}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Class</th><th>Body</th><th>Life Chemistry</th><th>System</th><th>Sector</th></tr></thead>
  <tbody>{body_rows}</tbody>
</table></div>
{_truncated_note(len(rows), result["truncated"])}
</section>
"""


def _moons_panel(db_name, result):
    rows = result["rows"]
    body_rows = "".join(
        "<tr>"
        f'<td>{esc(row["name"])}</td>'
        f'<td>{esc(row["planet_class"]) or "&mdash;"}</td>'
        f'<td>{"Gas Giant" if row["body_type"] == "g" else "Terrestrial"}</td>'
        f'<td>{esc(row["life_chemical"]) or "&mdash;"}</td>'
        f'<td>{esc(row["planet_name"])}</td>'
        f'<td><a href="system.py?db={esc(db_name)}&id={row["star_system_id"]}">{esc(row["system_name"])}</a></td>'
        f'<td>{_sector_link(db_name, row["sector_id"])}</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="7"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Moons{_count_suffix(len(rows), result["truncated"])}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Class</th><th>Body</th><th>Life Chemistry</th><th>Orbits</th><th>System</th><th>Sector</th></tr></thead>
  <tbody>{body_rows}</tbody>
</table></div>
{_truncated_note(len(rows), result["truncated"])}
</section>
"""


def _belts_panel(db_name, result):
    rows = result["rows"]
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
<h2>Asteroid Belts{_count_suffix(len(rows), result["truncated"])}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Density</th><th>Composition</th><th>System</th><th>Sector</th></tr></thead>
  <tbody>{body_rows}</tbody>
</table></div>
{_truncated_note(len(rows), result["truncated"])}
</section>
"""


_PANEL_RENDERERS = {
    "sectors": _sectors_panel,
    "systems": _systems_panel,
    "stars": _stars_panel,
    "planets": _planets_panel,
    "moons": _moons_panel,
    "belts": _belts_panel,
}


def _render_results(db_name, results):
    sections = []
    for panel in ("sectors", "systems", "stars", "planets", "moons", "belts"):
        result = results[panel]
        if result is not None:
            sections.append(_PANEL_RENDERERS[panel](db_name, result))
    return "".join(sections) if sections else '<p class="hint">No matching objects.</p>'


def handler():
    raw = parse_qs(os.environ.get("QUERY_STRING", ""))

    def _first(key):
        values = raw.get(key)
        return values[0].strip() if values else ""

    db_name = _first("db")
    texts = {
        "sector_q": _first("sector_q"),
        "system_q": _first("system_q"),
        "star_q": _first("star_q"),
        "planet_q": _first("planet_q"),
        "moon_q": _first("moon_q"),
    }
    tags = {facet: {v.strip() for v in raw.get(facet, []) if v.strip()} for facet in TAG_FACETS}
    tags["type"] &= {"star", "planet", "moon", "belt"}
    tags["body"] &= {"t", "g"}
    tags["moon_body"] &= {"t", "g"}

    state = {"db": db_name, "texts": texts, "tags": tags}

    data = get_search(db_name, texts, tags)

    facet_titles = {
        "type": "Object Type",
        "spectral": "Star: Spectral Class",
        "luminosity": "Star: Luminosity Class",
        "class": "Planet Class",
        "body": "Planet Body Type",
        "life": "Planet Supported Life Chemistry",
        "moon_class": "Moon Class",
        "moon_body": "Moon Body Type",
        "moon_life": "Moon Supported Life Chemistry",
        "density": "Asteroid Belt Density",
    }
    tag_browser_html = "".join(
        _tag_group_html(facet_titles[facet], facet, data["facets"][facet], tags[facet], state)
        for facet in TAG_FACETS
    )

    hidden_tag_inputs = "".join(
        f'<input type="hidden" name="{facet}" value="{esc(v)}">'
        for facet in TAG_FACETS
        for v in sorted(tags[facet])
    )
    datalists = "".join([
        _datalist_html("dl-sector", data["autocomplete"]["sectors"]),
        _datalist_html("dl-system", data["autocomplete"]["systems"]),
        _datalist_html("dl-star", data["autocomplete"]["stars"]),
        _datalist_html("dl-planet", data["autocomplete"]["planets"]),
        _datalist_html("dl-moon", data["autocomplete"]["moons"]),
    ])

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
    <label class="search-field">Planet name
      <input type="text" name="planet_q" value="{esc(texts['planet_q'])}" list="dl-planet" autocomplete="off" placeholder="e.g. Kepler-42 b">
    </label>
    <label class="search-field">Moon name
      <input type="text" name="moon_q" value="{esc(texts['moon_q'])}" list="dl-moon" autocomplete="off" placeholder="e.g. Kepler-42 b I">
    </label>
  </div>
  {datalists}
  <div class="search-actions">
    <button type="submit" class="btn">Search</button>
    <a href="search.py?db={esc(db_name)}">Clear all filters</a>
  </div>
</form>
"""

    active_filters_html = _active_filters_html(state, data["facet_labels"])

    any_active = bool(any(tags[facet] for facet in TAG_FACETS) or any(texts.values()))
    if any_active:
        results_html = _render_results(db_name, data["results"])
    else:
        results_html = '<p class="hint">Select a tag below, or enter a name above and press Search, to see matching results.</p>'

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
