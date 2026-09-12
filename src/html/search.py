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

# Min/max size (radius_km) query-string field names, one pair per entity
# with a size of its own -- threaded through `_build_url`/`handler()`
# exactly like the plain text fields (`sector_q` etc.) below: opaque
# strings preserved verbatim across a tag toggle, parsed to a float only
# at the `get_search()` call site (see `_size_range_from_texts`).
SIZE_ENTITIES = ("star", "planet", "moon")
SIZE_FIELDS = tuple(f"{entity}_{bound}_radius_km" for entity in SIZE_ENTITIES for bound in ("min", "max"))


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
    for key in ("sector_q", "system_q", "star_q", "planet_q", "moon_q") + SIZE_FIELDS:
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


def _remove_size_url(state, entity):
    new_texts = dict(state["texts"])
    new_texts[f"{entity}_min_radius_km"] = ""
    new_texts[f"{entity}_max_radius_km"] = ""
    return _build_url(state["db"], new_texts, state["tags"])


def _size_range_from_texts(texts, entity):
    """
    Parses `texts["<entity>_min_radius_km"]`/`"_max_radius_km"` (plain
    strings from the query string, same as every other search field) into
    a `(min_km, max_km)` tuple for `apiclient.get_search`'s `sizes`
    argument -- `None` for a bound that's empty or, from a hand-edited
    URL, not actually a valid number (silently ignored rather than
    erroring this read-only browse page over it; `GET /api/search`
    itself still validates strictly). Returns `None` altogether when
    neither bound is set, matching `get_search`'s own "absent means no
    filter" convention.
    """
    def _parse(key):
        raw = texts.get(key, "")
        if not raw:
            return None
        try:
            return float(raw)
        except ValueError:
            return None

    min_km = _parse(f"{entity}_min_radius_km")
    max_km = _parse(f"{entity}_max_radius_km")
    if min_km is None and max_km is None:
        return None
    return (min_km, max_km)


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
    size_labels = {"star": "Star size", "planet": "Planet size", "moon": "Moon size"}
    for entity, label in size_labels.items():
        size_range = _size_range_from_texts(state["texts"], entity)
        if size_range is not None:
            min_km, max_km = size_range
            if min_km is not None and max_km is not None:
                range_text = f"{min_km:,.0f}&ndash;{max_km:,.0f} km"
            elif min_km is not None:
                range_text = f"&ge; {min_km:,.0f} km"
            else:
                range_text = f"&le; {max_km:,.0f} km"
            chips.append((f'{esc(label)}: {range_text}', _remove_size_url(state, entity)))
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


def _radius_km(row):
    return f'{row["radius_km"]:,.0f} km' if row.get("radius_km") is not None else "&mdash;"


def _stars_panel(db_name, result):
    rows = result["rows"]
    body_rows = "".join(
        "<tr>"
        f'<td>{esc(row["name"])}</td>'
        f'<td>{esc(row["role"])}</td>'
        f'<td>{esc(row["star_type"])}</td>'
        f'<td>{_radius_km(row)}</td>'
        f'<td><a href="system.py?db={esc(db_name)}&id={row["star_system_id"]}">{esc(row["system_name"])}</a></td>'
        f'<td>{_sector_link(db_name, row["sector_id"])}</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="6"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Stars{_count_suffix(len(rows), result["truncated"])}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Role</th><th>Type</th><th>Radius</th><th>System</th><th>Sector</th></tr></thead>
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
        f'<td>{_radius_km(row)}</td>'
        f'<td>{esc(row["life_chemical"]) or "&mdash;"}</td>'
        f'<td><a href="system.py?db={esc(db_name)}&id={row["star_system_id"]}">{esc(row["system_name"])}</a></td>'
        f'<td>{_sector_link(db_name, row["sector_id"])}</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="7"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Planets{_count_suffix(len(rows), result["truncated"])}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Class</th><th>Body</th><th>Radius</th><th>Life Chemistry</th><th>System</th><th>Sector</th></tr></thead>
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
        f'<td>{_radius_km(row)}</td>'
        f'<td>{esc(row["life_chemical"]) or "&mdash;"}</td>'
        f'<td>{esc(row["planet_name"])}</td>'
        f'<td><a href="system.py?db={esc(db_name)}&id={row["star_system_id"]}">{esc(row["system_name"])}</a></td>'
        f'<td>{_sector_link(db_name, row["sector_id"])}</td>'
        "</tr>"
        for row in rows
    ) or '<tr><td colspan="8"><em>None</em></td></tr>'
    return f"""
<section class="panel">
<h2>Moons{_count_suffix(len(rows), result["truncated"])}</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Class</th><th>Body</th><th>Radius</th><th>Life Chemistry</th><th>Orbits</th><th>System</th><th>Sector</th></tr></thead>
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
    for field in SIZE_FIELDS:
        texts[field] = _first(field)
    tags = {facet: {v.strip() for v in raw.get(facet, []) if v.strip()} for facet in TAG_FACETS}
    tags["type"] &= {"star", "planet", "moon", "belt"}
    tags["body"] &= {"t", "g"}
    tags["moon_body"] &= {"t", "g"}

    state = {"db": db_name, "texts": texts, "tags": tags}

    sizes = {entity: _size_range_from_texts(texts, entity) for entity in SIZE_ENTITIES}
    data = get_search(db_name, texts, tags, sizes=sizes)

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

    def _size_field_html(entity, label, placeholder_min, placeholder_max):
        min_key, max_key = f"{entity}_min_radius_km", f"{entity}_max_radius_km"
        return f"""
    <div class="search-field search-field-size">
      <span>{esc(label)} radius (km)</span>
      <input type="number" name="{min_key}" value="{esc(texts[min_key])}" min="0" step="any" placeholder="{placeholder_min}" aria-label="Minimum {esc(label)} radius, km">
      <span>&ndash;</span>
      <input type="number" name="{max_key}" value="{esc(texts[max_key])}" min="0" step="any" placeholder="{placeholder_max}" aria-label="Maximum {esc(label)} radius, km">
    </div>"""

    size_fields_html = "".join([
        _size_field_html("star", "Star", "e.g. 0", "e.g. 696000 (Sun)"),
        _size_field_html("planet", "Planet", "e.g. 0", "e.g. 6371 (Earth)"),
        _size_field_html("moon", "Moon", "e.g. 0", "e.g. 1737 (Moon)"),
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
  <div class="search-fields search-fields-size">
    {size_fields_html}
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
