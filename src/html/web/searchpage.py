# html/web/searchpage.py

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

from urllib.parse import urlencode

from flask import request, url_for

from pagination import PAGE_SIZE, page_offset, parse_page, render_pagination

from .helpers import page_url, trusted_html

# Mirrors queryDb.SEARCH_TAG_FACETS (this layer talks to the database only
# through the API, so it doesn't import queryDb).
TAG_FACETS = (
    "type", "spectral", "luminosity",
    "class", "body", "life",
    "moon_class", "moon_body", "moon_life",
    "density",
)

FACET_TITLES = {
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

# Facets whose values are a fixed set; anything else from a hand-edited
# URL is dropped.
_FACET_ALLOWED = {
    "type": {"star", "planet", "moon", "belt"},
    "body": {"t", "g"},
    "moon_body": {"t", "g"},
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
            return value if value is None or value == value else None  # drop NaN

        low, high = parse(f"{entity}_min_radius_km"), parse(f"{entity}_max_radius_km")
        if low is None and high is None:
            return None
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


def _km(value):
    return f"{value:,.0f} km"


def active_filters(state, facet_labels):
    """The removable filter chips: `[{"label", "url"}]`, each URL being
    this search without that one filter."""
    chips = []
    if state.q:
        chips.append({"label": f"Name: “{state.q}”", "url": state.copy(q="").url()})
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
        elif panel == "systems":
            item["url"] = page_url("system", system_id=row["id"])
            item["sector_url"] = _sector_cell(row["sector_id"])
        else:
            item["system_url"] = page_url("system", system_id=row["star_system_id"])
            item["sector_url"] = _sector_cell(row["sector_id"])
        if "radius_km" in row:
            item["radius"] = _km(row["radius_km"]) if row["radius_km"] is not None else None
        if "body_type" in row:
            item["body"] = "Gas Giant" if row["body_type"] == "g" else "Terrestrial"
        out.append(item)
    return out


def result_panels(state, results):
    """The result panels that ran, each `{"panel", "heading", "total",
    "rows", "pager"}`. Every pager keeps the search and every other
    panel's page as the API actually returned it (a past-the-end page
    comes back as the last one)."""
    shown = [(panel, heading, results.get(panel)) for panel, heading in RESULT_PANELS
             if results.get(panel) is not None]
    pages = {panel: result["offset"] // result["limit"] + 1 for panel, _h, result in shown}
    paged = state.copy(pages=dict(state.pages, **pages))
    base = paged.params(with_pages=True)
    action = url_for("web.search")
    panels = []
    for panel, heading, result in shown:
        panels.append({
            "panel": panel,
            "heading": heading,
            "total": result["total"],
            "rows": _rows(panel, result["rows"]),
            "pager": trusted_html(render_pagination(
                action, base, f"{panel}_page", pages[panel], result["total"], page_size=result["limit"],
                anchor=f"search-{panel}", label=f"{heading} result pages", method="get",
            )),
        })
    return panels


def api_offsets(state):
    return {panel: page_offset(page) for panel, page in state.pages.items()}


__all__ = [
    "PAGE_SIZE", "NAME_FIELDS", "SIZE_ENTITIES", "SearchState", "needs_canonical_redirect",
    "tag_groups", "active_filters", "result_panels", "api_offsets",
]
