# planetgen/web/searchpage.py

"""
The search page (`/search`, was `search.py`): request parsing, URL
building and template data. The route itself is `views.search`.

Every filter is a GET query parameter, so a search is an ordinary,
bookmarkable URL:

- `q`: the main name search (the header search box sends it). It
  searches the names of sectors, star systems, stars, planets and moons
  at once. A per-object field below overrides it for that object; when an
  Object Type tag is active, `q` searches only the selected object types
  (not sectors and systems).
- `sector_q`, `system_q`, `star_q`, `planet_q`, `moon_q`: a name search
  for one kind of object.
- `<star|planet|moon>_<min|max>_radius_km`: a size range.
- One repeated parameter per tag facet (`type=star&spectral=G&spectral=K`),
  see `TAG_FACETS`.
- `<panel>_page`: each result panel's page (`sectors_page`, ...,
  `belts_page`).

The query logic itself lives in `queryDb.search`, behind `GET /api/search`
(`apiclient.get_search`, called in-process).
"""

import math
from urllib.parse import urlencode

from flask import g, request, url_for

from planetgen.physics.habitability_world import EQUIPMENT_LABELS
from planetgen.web.lib import apiclient
from planetgen.web.lib.datatable import Column, Result, Table
from planetgen.web.lib.fmt import format_number
from planetgen.web.lib.pagination import PAGE_SIZE, page_offset, parse_page

from . import tables
from .helpers import db_name, page_url
from .sector_page import PHENOMENON_TYPE_LABELS

# Mirrors queryDb.SEARCH_TAG_FACETS (this layer talks to the database only
# through the API, so it doesn't import planetgen.db.query).
TAG_FACETS = (
    "type", "spectral", "luminosity",
    "class", "body", "life", "equipment",
    "moon_class", "moon_body", "moon_life", "moon_equipment",
    "density",
    "phenomenon", "phenomenon_class",
)

FACET_TITLES = {
    "type": "Object Type",
    "spectral": "Star: Spectral Class",
    "luminosity": "Star: Luminosity Class",
    "class": "Planet Class",
    "body": "Planet Body Type",
    "life": "Planet Supported Life Chemistry",
    "equipment": "Planet: Equipment a Human Needs",
    "moon_class": "Moon Class",
    "moon_body": "Moon Body Type",
    "moon_life": "Moon Supported Life Chemistry",
    "moon_equipment": "Moon: Equipment a Human Needs",
    "density": "Asteroid Belt Density",
    "phenomenon": "Phenomenon",
    "phenomenon_class": "Phenomenon Class",
}

# Facets whose values are a fixed set; anything else from a hand-edited
# URL is dropped.
_FACET_ALLOWED = {
    "type": {"star", "planet", "moon", "belt"},
    "body": {"t", "g"},
    "moon_body": {"t", "g"},
    "equipment": {"0", "1", "2", "3", "4"},
    "moon_equipment": {"0", "1", "2", "3", "4"},
    "phenomenon": {"nebula", "asteroid_field", "black_hole", "neutron_star", "supernova_remnant",
                   "rogue_planet", "interstellar_comet", "quasar"},
}

NAME_FIELDS = (
    # (parameter, label, autocomplete key, placeholder)
    ("sector_q", "Sector", "sectors", "e.g. Voranthis Kelmoor"),
    ("system_q", "System", "systems", "e.g. Kepler-42"),
    ("star_q", "Star", "stars", "e.g. Kepler-42 A"),
    ("planet_q", "Planet", "planets", "e.g. Kepler-42 b"),
    ("moon_q", "Moon", "moons", "e.g. Kepler-42 b I"),
)

# The name fields `q` fills in when no Object Type tag is active; with one,
# only the object fields (the type tags pick which of those panels show).
_Q_ALL = ("sector_q", "system_q", "star_q", "planet_q", "moon_q")
_Q_TYPED = ("star_q", "planet_q", "moon_q")

SIZE_ENTITIES = (
    # (entity, label, max placeholder)
    ("star", "Star", "e.g. 696000 (Sun)"),
    ("planet", "Planet", "e.g. 6371 (Earth)"),
    ("moon", "Moon", "e.g. 1737 (Moon)"),
)
SIZE_FIELDS = tuple(f"{entity}_{bound}_radius_km" for entity, _l, _p in SIZE_ENTITIES for bound in ("min", "max"))

RESULT_PANELS = (
    # (panel, heading) -- mirrors queryDb.SEARCH_RESULT_PANELS
    ("sectors", "Sectors"),
    ("systems", "Systems"),
    ("stars", "Stars"),
    ("planets", "Planets"),
    ("moons", "Moons"),
    ("belts", "Asteroid Belts"),
    ("phenomena", "Phenomena"),
)
PANEL_NAMES = tuple(panel for panel, _heading in RESULT_PANELS)

_ALL_PARAMS = ("q",) + tuple(f[0] for f in NAME_FIELDS) + SIZE_FIELDS + TAG_FACETS + tuple(
    f"{panel}_page" for panel in PANEL_NAMES)


class SearchState:
    """One search: the request's filters, parsed and cleaned."""

    def __init__(self, q="", texts=None, tags=None, pages=None):
        self.q = q
        self.texts = {key: "" for key in tuple(f[0] for f in NAME_FIELDS) + SIZE_FIELDS}
        self.texts.update(texts or {})
        self.tags = {facet: set() for facet in TAG_FACETS}
        for facet, values in (tags or {}).items():
            self.tags[facet] = set(values)
        self.pages = {panel: 1 for panel in PANEL_NAMES}
        self.pages.update(pages or {})

    @classmethod
    def from_args(cls, args):
        """Parses a request's query (`request.args`)."""
        def first(key):
            return (args.get(key) or "").strip()

        texts = {key: first(key) for key in tuple(f[0] for f in NAME_FIELDS) + SIZE_FIELDS}
        tags = {}
        for facet in TAG_FACETS:
            values = {value.strip() for value in args.getlist(facet) if value.strip()}
            if facet in _FACET_ALLOWED:
                values &= _FACET_ALLOWED[facet]
            tags[facet] = values
        pages = {panel: parse_page(args.get(f"{panel}_page")) for panel in PANEL_NAMES}
        return cls(first("q"), texts, tags, pages)

    def copy(self, **changes):
        state = SearchState(self.q, dict(self.texts), {k: set(v) for k, v in self.tags.items()},
                            dict(self.pages))
        for key, value in changes.items():
            setattr(state, key, value)
        return state

    # -- what gets searched ------------------------------------------------

    def name_terms(self):
        """The five name terms sent to the API: each per-object field, or
        `q` where that field is empty (see the module docstring)."""
        fills = _Q_TYPED if self.tags["type"] else _Q_ALL
        terms = {}
        for key, _label, _ac, _ph in NAME_FIELDS:
            terms[key] = self.texts[key] or (self.q if key in fills else "")
        return terms

    def size_range(self, entity):
        """`(min_km, max_km)` for `apiclient.get_search`'s `sizes`, or
        `None` when neither bound is a number (a non-number from a
        hand-edited URL is ignored, not an error)."""
        def parse(key):
            raw = self.texts.get(key, "")
            try:
                value = float(raw) if raw else None
            except ValueError:
                return None
            # Only a finite, non-negative size is a size: NaN, infinity
            # and negatives from a hand-edited URL are ignored like any
            # other non-number (the API would refuse them).
            if value is None or not math.isfinite(value) or value < 0:
                return None
            return value

        low, high = parse(f"{entity}_min_radius_km"), parse(f"{entity}_max_radius_km")
        if low is None and high is None:
            return None
        if low is not None and high is not None and low > high:
            return None  # an inverted range, likewise ignored
        return (low, high)

    def sizes(self):
        return {entity: self.size_range(entity) for entity, _l, _p in SIZE_ENTITIES}

    def active(self):
        """True when anything at all is being searched for."""
        return bool(self.q or any(self.texts.values()) or any(self.tags.values()))

    def advanced_open(self):
        """Open the per-object fields when one of them is in use."""
        return any(self.texts.values())

    # -- URLs ----------------------------------------------------------------

    def params(self, with_pages=False):
        """The query parameters for this search, in a stable order, empty
        values left out."""
        pairs = []
        if self.q:
            pairs.append(("q", self.q))
        for key in tuple(f[0] for f in NAME_FIELDS) + SIZE_FIELDS:
            if self.texts.get(key):
                pairs.append((key, self.texts[key]))
        for facet in TAG_FACETS:
            for value in sorted(self.tags[facet]):
                pairs.append((facet, value))
        if with_pages:
            for panel in PANEL_NAMES:
                if self.pages.get(panel, 1) > 1:
                    pairs.append((f"{panel}_page", self.pages[panel]))
        return pairs

    def url(self, with_pages=False):
        pairs = self.params(with_pages)
        base = url_for("web.search")
        return f"{base}?{urlencode(pairs)}" if pairs else base


def needs_canonical_redirect(args):
    """True when the query carries empty or unknown parameters (a
    submitted form sends every field, filled or not), so the view can
    redirect to the short form of the same search."""
    for key, value in args.items(multi=True):
        if key not in _ALL_PARAMS or not value.strip():
            return True
    return False


# ---------------------------------------------------------------------
# Template data
# ---------------------------------------------------------------------

def _toggle(state, facet, value):
    tags = {k: set(v) for k, v in state.tags.items()}
    tags[facet] ^= {value}
    return state.copy(tags=tags).url()


def tag_groups(state, facets):
    """The "Browse by Tag" groups: one per facet with any values."""
    groups = []
    for facet in TAG_FACETS:
        options = facets.get(facet) or []
        if not options:
            continue
        groups.append({
            "title": FACET_TITLES[facet],
            "options": [{
                "label": option["label"],
                "count": option["count"],
                "tooltip": option.get("tooltip"),
                "active": option["value"] in state.tags[facet],
                "url": _toggle(state, facet, option["value"]),
            } for option in options],
        })
    return groups


# Body radii are the distance ladder's exception (always km), and these
# chips echo the km the visitor typed, so they stay plain km.
def _km(value):
    return f"{value:,.0f} km"


def active_filters(state, facet_labels):
    """The removable filter chips: `[{"label", "url"}]`, each URL being
    this search without that one filter."""
    chips = []  # the main name is not echoed: the input above says it (UX.54)
    for key, label, _ac, _ph in NAME_FIELDS:
        if state.texts[key]:
            texts = dict(state.texts, **{key: ""})
            chips.append({"label": f"{label}: “{state.texts[key]}”", "url": state.copy(texts=texts).url()})
    for entity, label, _ph in SIZE_ENTITIES:
        size_range = state.size_range(entity)
        if size_range is None:
            continue
        low, high = size_range
        if low is not None and high is not None:
            text = f"{low:,.0f}–{_km(high)}"
        elif low is not None:
            text = f"≥ {_km(low)}"
        else:
            text = f"≤ {_km(high)}"
        texts = dict(state.texts, **{f"{entity}_min_radius_km": "", f"{entity}_max_radius_km": ""})
        chips.append({"label": f"{label} size: {text}", "url": state.copy(texts=texts).url()})
    for facet in TAG_FACETS:
        for value in sorted(state.tags[facet]):
            chips.append({"label": facet_labels.get(f"{facet}:{value}", value), "url": _toggle(state, facet, value)})
    return chips


def _sector_cell(sector_id):
    return None if sector_id is None else page_url("sector", sector_id=sector_id)


def _rows(panel, rows):
    """Template-ready rows for one result panel."""
    out = []
    for row in rows:
        item = dict(row)
        if panel == "sectors":
            item["url"] = page_url("sector", sector_id=row["id"])
            item["galaxy_url"] = page_url("sector_on_galaxy_map", sector_id=row["id"])
        elif panel == "phenomena":
            item["url"] = page_url("phenomenon", phenomenon_type=row["type"], phenomenon_id=row["id"])
            item["type_label"] = PHENOMENON_TYPE_LABELS.get(row["type"], row["type"])
            item["sector_url"] = _sector_cell(row["sector_id"])
        elif panel == "systems":
            item["url"] = page_url("system", system_id=row["id"])
            item["sector_url"] = _sector_cell(row["sector_id"])
            if row["sector_id"] is not None:
                item["galaxy_url"] = page_url("sector_on_galaxy_map", sector_id=row["sector_id"])
        else:
            item["system_url"] = page_url("system", system_id=row["star_system_id"])
            item["sector_url"] = _sector_cell(row["sector_id"])
        if "radius_km" in row:
            item["radius"] = _km(row["radius_km"]) if row["radius_km"] is not None else None
        if "body_type" in row:
            item["body"] = "Gas Giant" if row["body_type"] == "g" else "Terrestrial"
        out.append(item)
    return out


def _link(text, href):
    return {"text": text, "href": href}


def _sector_link(row):
    return _link("View", row["sector_url"]) if row["sector_url"] else {"text": "Standalone"}


def _radius(row):
    return {"text": row["radius"]} if row.get("radius") else {"text": "—"}


def _show_on_map(row):
    return (_link("Show", row["galaxy_url"]) | {"label": f"Show {row['name']} on the Galaxy Map"}
            if row.get("galaxy_url") else {"text": "—"})


_PANEL_COLUMNS = {
    "sectors": ["Name", "Edge", "Galaxy Map"],
    "systems": ["Name", "Sector", "Binary", "Star type", "Galaxy Map"],
    "stars": ["Name", "Role", "Type", "Radius", "System", "Sector"],
    "planets": ["Name", "Class", "Body", "Radius", "Life Chemistry", "Equipment", "System", "Sector"],
    "moons": ["Name", "Class", "Body", "Radius", "Life Chemistry", "Equipment", "Orbits", "System", "Sector"],
    "belts": ["Density", "Composition", "System", "Sector"],
    "phenomena": ["Name", "Type", "Class", "Sector"],
}
_PANEL_NOUNS = {
    "sectors": ("sector", "sectors"), "systems": ("system", "systems"), "stars": ("star", "stars"),
    "planets": ("planet", "planets"), "moons": ("moon", "moons"), "belts": ("asteroid belt", "asteroid belts"),
    "phenomena": ("phenomenon", "phenomena"),
}


def _cells(panel, row):
    """One search result row as the data table's cells (the columns of `_PANEL_COLUMNS`)."""
    if panel == "sectors":
        return [_link(row["name"], row["url"]), {"text": f"{format_number(row['edge_mpc'], ',.2f')} mpc"},
                _show_on_map(row)]
    if panel == "systems":
        return [_link(row["name"], row["url"]), _sector_link(row), {"text": "Yes" if row["is_binary"] else "No"},
                {"text": row["star_summary"]}, _show_on_map(row)]
    if panel == "stars":
        return [{"text": row["name"]}, {"text": row["role"]}, {"text": row["star_type"]}, _radius(row),
                _link(row["system_name"], row["system_url"]), _sector_link(row)]
    if panel in ("planets", "moons"):
        cells = [{"text": row["name"]}, {"text": row["planet_class"] or "—"}, {"text": row["body"]}, _radius(row),
                 {"text": row["life_chemical"] or "—"},
                 {"text": EQUIPMENT_LABELS[row["equipment_tier"]] if row.get("equipment_tier") is not None else "—"}]
        if panel == "moons":
            cells.append({"text": row["planet_name"]})
        return cells + [_link(row["system_name"], row["system_url"]), _sector_link(row)]
    if panel == "phenomena":
        return [_link(row["name"], row["url"]), {"text": row["type_label"]},
                {"text": row["phenomenon_class"] or "—"}, _sector_link(row)]
    return [{"text": row["density"].capitalize()}, {"text": row["composition_summary"]},
            _link(row["system_name"], row["system_url"]), _sector_link(row)]


def _panel_loader(panel):
    def load(_table_state, limit, offset, _want_facets):
        cached = (getattr(g, "search_results", None) or {}).get(panel)
        if cached is None or cached["offset"] != offset or cached["limit"] != limit:
            state = SearchState.from_args(request.args)
            data = apiclient.get_search(
                db_name(), state.name_terms(), state.tags, sizes=state.sizes(), limit=limit,
                offsets={panel: offset}, panels=[panel])
            cached = data["results"][panel]
        if cached is None:
            return Result([], 0)
        rows = _rows(panel, cached["rows"])
        return Result([_cells(panel, row) for row in rows], cached["total"], None,
                      cached.get("total_capped", False))
    return load


PANEL_TABLES = {
    panel: tables.register(Table(
        f"search-{panel}", heading,
        [Column(label.lower().replace(" ", "_"), label, sortable=False) for label in _PANEL_COLUMNS[panel]],
        _panel_loader(panel), prefix=f"{panel}_", noun=_PANEL_NOUNS[panel]))
    for panel, heading in RESULT_PANELS
}


def result_panels(state, results):
    """The result panels that ran, each `{"panel", "heading", "total",
    "total_capped", "table"}` (`total_capped`: more matches than `total`,
    shown as "300+"), each a data table (UX.41) that keeps the search and
    every other panel's page. The API's results are handed to the tables
    so the page asks it once."""
    g.search_results = results
    path = url_for("web.search")
    source = {}
    for key, value in state.params():
        source.setdefault(key, []).append(value)
    panels = []
    for panel, heading in RESULT_PANELS:
        result = results.get(panel)
        if result is None or not result["total"]:
            continue  # a group with no matches is left out (UX.54)
        keep = tuple(key for key in _ALL_PARAMS if key != f"{panel}_page")
        panels.append({
            "panel": panel,
            "heading": heading,
            "total": result["total"],
            "total_capped": result.get("total_capped", False),
            "table": tables.render(PANEL_TABLES[panel], path, anchor=f"search-{panel}", keep=keep, **source),
        })
    return panels


def api_offsets(state):
    return {panel: page_offset(page) for panel, page in state.pages.items()}


__all__ = [
    "PAGE_SIZE", "NAME_FIELDS", "SIZE_ENTITIES", "SearchState", "needs_canonical_redirect",
    "tag_groups", "active_filters", "result_panels", "api_offsets",
]
