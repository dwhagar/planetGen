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
- `<star|planet|moon>_<min|max>_<key>` and `<entity>_<key>_unit`: a range on a physical value in the unit named (UX.227, the keys of
  `planetgen.db.search_values`: gravity, surface temperature, pressure, dose, orbit distance, mass, the partial pressure of a gas, ...), and
  `<entity>_value_sort` (`gravity` or `-gravity`) to order that kind of object by one. `planet_other_gas` with `planet_other_gas_min`,
  `_max` and `_unit` (likewise `moon_`) is the form's way to range any gas outside its short list: it turns into the key of that gas.
- `stars`: how many stars a system holds (UX.239): `3` exactly three, `3+` three or more, or `single`, `binary`, `trinary`, `quaternary`,
  `quintenary` (five or more). A `stars:<value>` word typed into the main search `q` (the header search box) is the same filter.
- One repeated parameter per tag facet (`type=star&spectral=G&spectral=K`),
  see `TAG_FACETS`.
- `<panel>_page`: each result panel's page (`sectors_page`, ...,
  `belts_page`).

The query logic itself lives in `planetgen.db.query.search`, behind `GET /api/search`
(`apiclient.get_search`, called in-process).
"""

import functools
import math
import re
from urllib.parse import urlencode

from flask import g, request, url_for

from planetgen.db import search_values
from planetgen.physics import equipment, score_filters
from planetgen.physics.habitability_world import TIER_NAMES
from planetgen.util import star_count
from planetgen.web.lib import apiclient, datatable, slowpages
from planetgen.web.lib.datatable import Column, Result, Table
from planetgen.web.lib.fmt import format_number
from planetgen.web.lib.habitability_cells import equipment_cell, habitability_cell
from planetgen.web.lib.habitability_symbols import LIFE_NAMES
from planetgen.web.lib.pagination import PAGE_SIZE, page_offset, parse_page

from . import tables
from .helpers import db_name, page_url
from .item_pages import item_url
from .sector_page import PHENOMENON_TYPE_LABELS

# Mirrors planetgen.db.query.SEARCH_TAG_FACETS (this layer talks to the database only
# through the API, so it doesn't import planetgen.db.query).
TAG_FACETS = (
    "type", "spectral", "luminosity",
    "class", "body", "life", "equipment",
    "moon_class", "moon_body", "moon_life", "moon_equipment",
    *score_filters.FACETS, *(f"moon_{name}" for name in score_filters.FACETS),
    "density",
    "phenomenon", "phenomenon_class",
)

FACET_TITLES = {
    "type": "Object Type",
    "spectral": "Star: Spectral Class",
    "luminosity": "Star: Luminosity Class",
    "class": "Planet Class",
    "body": "Planet Body Type",
    "life": "Planet Life Class",
    "equipment": "Planet: Equipment a Human Needs",
    "moon_class": "Moon Class",
    "moon_body": "Moon Body Type",
    "moon_life": "Moon Life Class",
    "moon_equipment": "Moon: Equipment a Human Needs",
    **{name: f"Planet: {title}" for name, (title, _ticks) in score_filters.FACETS.items()},
    **{f"moon_{name}": f"Moon: {title}" for name, (title, _ticks) in score_filters.FACETS.items()},
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
    "equipment": equipment.TICK_VALUES,
    "moon_equipment": equipment.TICK_VALUES,
    **{name: score_filters.values(name) for name in score_filters.FACETS},
    **{f"moon_{name}": score_filters.values(name) for name in score_filters.FACETS},
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

VALUE_ENTITIES = (("planet", "Planet"), ("moon", "Moon"), ("star", "Star"))
"""tuple: The kinds of object that take ranges on physical values (UX.227), `(entity, label)`, in the order the form shows them."""
VALUE_KEYS = {entity: tuple(key for key, value in search_values.KEYS.items() if entity in value.entities) for entity, _label in VALUE_ENTITIES}
"""dict: Entity -> the keys of `planetgen.db.search_values` it has."""
DEFAULT_GASES = ("N2", "O2", "CO2", "H2O", "CH4", "Ar", "H2", "He", "NH3", "SO2")
"""tuple: The gases the form shows a row for; any other gas is ranged through the "another gas" row (or the URL)."""
OTHER_GAS = "other_gas"
"""str: The form's "another gas" row is `<entity>_other_gas` (the formula) with `_min`, `_max` and `_unit`."""
GAS_ENTITIES = ("planet", "moon")


def value_field_names(entity, key):
    """`(minimum, maximum, unit)` parameter names of `key` for `entity`."""
    return (f"{entity}_min_{key}", f"{entity}_max_{key}", f"{entity}_{key}_unit")


VALUE_FIELDS = tuple(name for entity, _label in VALUE_ENTITIES for key in VALUE_KEYS[entity] for name in value_field_names(entity, key))
"""tuple: Every range parameter of a physical value (UX.227)."""
VALUE_SORT_FIELDS = tuple(f"{entity}_value_sort" for entity, _label in VALUE_ENTITIES)
"""tuple: The parameters that order a kind of object by a physical value."""
_OTHER_GAS_FIELDS = tuple(f"{entity}_{OTHER_GAS}{part}" for entity in GAS_ENTITIES for part in ("", "_min", "_max", "_unit"))
STAR_COUNT_CHOICES = star_count.STAR_COUNT_CHOICES
"""tuple: The stars filter's choices, `(value, label)`, for the form's menu."""
STARS_FIELD = "stars"
"""str: The stars filter's parameter (UX.239); it is kept with the other text fields of a state."""
TEXT_FIELDS = tuple(f[0] for f in NAME_FIELDS) + SIZE_FIELDS + VALUE_FIELDS + VALUE_SORT_FIELDS + (STARS_FIELD,)
"""tuple: Every parameter a state holds as text."""
_Q_STARS = re.compile(r"(?:^|\s)stars:(\S+)", re.IGNORECASE)

RESULT_PANELS = (
    # (panel, heading) -- mirrors planetgen.db.query.SEARCH_RESULT_PANELS
    ("sectors", "Sectors"),
    ("systems", "Systems"),
    ("stars", "Stars"),
    ("planets", "Planets"),
    ("moons", "Moons"),
    ("belts", "Asteroid Belts"),
    ("phenomena", "Phenomena"),
)
PANEL_NAMES = tuple(panel for panel, _heading in RESULT_PANELS)

_ALL_PARAMS = ("q",) + TEXT_FIELDS + TAG_FACETS + tuple(
    f"{panel}_{part}" for panel in PANEL_NAMES for part in ("page", "sort", "order"))


def _valid_stars(text):
    """The stars filter as typed when it is one (a count, a count and a plus, or a name), else empty: a hand-edited URL is ignored, not an error."""
    try:
        star_count.parse_star_count(text)
    except ValueError:
        return ""
    return text.strip().lower()


def _finite(text):
    """The number a range box holds, or `None` for an empty box or anything that is not a finite number (a hand-edited URL is ignored, not an error)."""
    try:
        value = float(text) if text else None
    except ValueError:
        return None
    return value if value is not None and math.isfinite(value) else None


def _plain(number):
    """A number as the range boxes and chips write it: the shortest form that reads back the same to six figures."""
    return f"{number:.6g}"


def _key_label(entity, key):
    """What a chip or a sort menu calls `key`: the registry's label with its first letter lowered ("surface gravity")."""
    label = search_values.KEYS[key].label.split(" (")[0]
    return label[:1].lower() + label[1:]


@functools.lru_cache(maxsize=None)
def _unit_symbol(unit):
    return search_values.symbol(unit)


def _tidy_units(texts):
    """Drops a unit that says nothing: the form sends every row's unit, so one with no bound beside it or the key's own unit is left out of the URL."""
    for entity, _label in VALUE_ENTITIES:
        for key in VALUE_KEYS[entity]:
            least_name, most_name, unit_name = value_field_names(entity, key)
            if texts.get(unit_name) and (not (texts.get(least_name) or texts.get(most_name))
                                         or texts[unit_name] == search_values.KEYS[key].default_unit(entity)):
                texts[unit_name] = ""


def _unit_ok(entity, key, unit):
    """True when `unit` measures what `key` measures for `entity` (any unit of its kind in the registry, as the API takes)."""
    try:
        search_values.convert(1.0, unit, search_values.KEYS[key].default_unit(entity))
    except (ValueError, KeyError, TypeError):
        return False
    return True


class SearchState:
    """One search: the request's filters, parsed and cleaned."""

    def __init__(self, q="", texts=None, tags=None, pages=None):
        self.q = q
        self.texts = {key: "" for key in TEXT_FIELDS}
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

        texts = {key: first(key) for key in TEXT_FIELDS}
        q = first("q")
        found = _Q_STARS.search(q)
        if found:       # `stars:3` typed in the search box is the stars filter, and the rest of the box stays a name
            q = _Q_STARS.sub("", q, count=1).strip()
            texts[STARS_FIELD] = texts[STARS_FIELD] or found.group(1)
        texts[STARS_FIELD] = _valid_stars(texts[STARS_FIELD])
        for entity in GAS_ENTITIES:      # the form's "another gas" row is the range of that gas's own key
            key = search_values.GAS_PREFIX + first(f"{entity}_{OTHER_GAS}")
            if key in search_values.KEYS and entity in search_values.KEYS[key].entities:
                for part, name in zip(("min", "max", "unit"), value_field_names(entity, key)):
                    texts[name] = texts[name] or first(f"{entity}_{OTHER_GAS}_{part}")
        _tidy_units(texts)
        tags = {}
        for facet in TAG_FACETS:
            values = {value.strip() for value in args.getlist(facet) if value.strip()}
            if facet in _FACET_ALLOWED:
                values &= _FACET_ALLOWED[facet]
            tags[facet] = values
        pages = {panel: parse_page(args.get(f"{panel}_page")) for panel in PANEL_NAMES}
        return cls(q, texts, tags, pages)

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

    def stars(self):
        """The stars filter as `apiclient.get_search` sends it (`stars=`), or `None` for no filter."""
        return self.texts.get(STARS_FIELD) or None

    def value_ranges(self):
        """
        The ranges on physical values as `apiclient.get_search` sends them (UX.227): `{entity: {key: (least, most, unit)}}`, each bound
        a number or `None` and `unit` the unit key or `None` for the key's own. A box that is not a finite number, an inverted range and
        a unit that does not measure the key are ignored, like a size range from a hand-edited URL.
        """
        ranges = {}
        for entity, _label in VALUE_ENTITIES:
            for key in VALUE_KEYS[entity]:
                least_name, most_name, unit_name = value_field_names(entity, key)
                least, most = _finite(self.texts.get(least_name)), _finite(self.texts.get(most_name))
                if least is None and most is None or (least is not None and most is not None and least > most):
                    continue
                unit = self.texts.get(unit_name) or None
                if unit is not None and (unit == search_values.KEYS[key].default_unit(entity) or not _unit_ok(entity, key, unit)):
                    unit = None
                ranges.setdefault(entity, {})[key] = (least, most, unit)
        return ranges

    def value_sorts(self):
        """The order a kind of object is asked for by a physical value, `{entity: "gravity" or "-gravity"}`; a key it does not have is ignored."""
        sorts = {}
        for entity, _label in VALUE_ENTITIES:
            sort = self.texts.get(f"{entity}_value_sort", "")
            if sort and search_values.sorts_in(sort, entity) is not None:
                sorts[entity] = sort
        return sorts

    def active(self):
        """True when anything at all is being searched for (an order by a value alone is not)."""
        return bool(self.q or any(value for key, value in self.texts.items() if key not in VALUE_SORT_FIELDS) or any(self.tags.values()))

    def advanced_open(self):
        """Open the per-object fields when one of them is in use (the stars filter sits in the main form)."""
        return any(value for key, value in self.texts.items() if key != STARS_FIELD and key not in VALUE_FIELDS and key not in VALUE_SORT_FIELDS)

    def values_open(self):
        """Open the physical values when one range or order is in use."""
        return any(self.texts.get(key) for key in VALUE_FIELDS + VALUE_SORT_FIELDS)

    # -- URLs ----------------------------------------------------------------

    def params(self, with_pages=False):
        """The query parameters for this search, in a stable order, empty
        values left out."""
        pairs = []
        if self.q:
            pairs.append(("q", self.q))
        for key in TEXT_FIELDS:
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
    sent = {key: (args.get(key) or "").strip() for key in VALUE_FIELDS}
    tidy = dict(sent)
    _tidy_units(tidy)
    return tidy != sent


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
    if state.stars():
        texts = dict(state.texts, **{STARS_FIELD: ""})
        chips.append({"label": f"Stars: {star_count.star_count_label(star_count.parse_star_count(state.stars()))}",
                      "url": state.copy(texts=texts).url()})
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
    entity_labels = dict(VALUE_ENTITIES)
    for entity, by_key in state.value_ranges().items():
        for key, (least, most, unit) in by_key.items():
            symbol = _unit_symbol(unit or search_values.KEYS[key].default_unit(entity))
            tail = f" {symbol}" if symbol else ""
            if least is not None and most is not None:
                text = f"{_plain(least)}–{_plain(most)}{tail}"
            elif least is not None:
                text = f"≥ {_plain(least)}{tail}"
            else:
                text = f"≤ {_plain(most)}{tail}"
            texts = dict(state.texts, **{name: "" for name in value_field_names(entity, key)})
            chips.append({"label": f"{entity_labels[entity]} {_key_label(entity, key)}: {text}", "url": state.copy(texts=texts).url()})
    for entity, sort in state.value_sorts().items():
        texts = dict(state.texts, **{f"{entity}_value_sort": ""})
        chips.append({"label": f"{entity_labels[entity]}s ordered by {_key_label(entity, sort.lstrip('-'))}, "
                               f"{'high to low' if sort.startswith('-') else 'low to high'}",
                      "url": state.copy(texts=texts).url()})
    for facet in TAG_FACETS:
        for value in sorted(state.tags[facet]):
            chips.append({"label": facet_labels.get(f"{facet}:{value}", value), "url": _toggle(state, facet, value)})
    return chips


def value_form(state):
    """
    The physical-value rows of the form (UX.227), one group per kind of object: `[{"entity", "label", "rows", "other_gas", "sort"}]`. A row is a key with its
    two boxes, its unit menu and what the state holds; the gas rows are the short list (`DEFAULT_GASES`) and any gas in use, and
    `other_gas` is the row for the rest. `sort` is the menu that orders that kind of object by a value.
    """
    ranges, sorts = state.value_ranges(), state.value_sorts()
    groups = []
    for entity, label in VALUE_ENTITIES:
        in_use = set(ranges.get(entity, {})) | ({sorts[entity].lstrip("-")} if entity in sorts else set())
        keys = [key for key in VALUE_KEYS[entity]
                if not key.startswith(search_values.GAS_PREFIX) or key[len(search_values.GAS_PREFIX):] in DEFAULT_GASES or key in in_use]
        rows = []
        for key in keys:
            least_name, most_name, unit_name = value_field_names(entity, key)
            default = search_values.KEYS[key].default_unit(entity)
            chosen = (ranges.get(entity, {}).get(key) or (None, None, None))[2] or default
            rows.append({
                "key": key, "label": search_values.KEYS[key].label, "min_name": least_name, "max_name": most_name, "unit_name": unit_name,
                "min": state.texts.get(least_name, ""), "max": state.texts.get(most_name, ""), "symbol": _unit_symbol(default),
                "units": [{"value": unit, "label": _unit_symbol(unit), "selected": unit == chosen} for unit in search_values.unit_choices(entity, key)],
            })
        other = None
        if entity in GAS_ENTITIES:
            shown = set(keys)
            other = {
                "name": f"{entity}_{OTHER_GAS}", "min_name": f"{entity}_{OTHER_GAS}_min", "max_name": f"{entity}_{OTHER_GAS}_max",
                "unit_name": f"{entity}_{OTHER_GAS}_unit",
                "gases": [key[len(search_values.GAS_PREFIX):] for key in VALUE_KEYS[entity]
                          if key.startswith(search_values.GAS_PREFIX) and key not in shown],
                "units": [{"value": unit, "label": _unit_symbol(unit), "selected": unit == "bar"}
                          for unit in search_values.unit_choices(entity, search_values.GAS_PREFIX + DEFAULT_GASES[0])],
            }
        options = [{"value": "", "label": "The page's own order", "selected": entity not in sorts}]
        for key in keys:
            name = search_values.KEYS[key].label.split(" (")[0]
            options.append({"value": key, "label": f"{name}, low to high", "selected": sorts.get(entity) == key})
            options.append({"value": f"-{key}", "label": f"{name}, high to low", "selected": sorts.get(entity) == f"-{key}"})
        groups.append({"entity": entity, "label": label, "rows": rows, "other_gas": other,
                       "sort": {"name": f"{entity}_value_sort", "options": options}})
    return groups


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
        elif row.get("star_system_id") is None:    # DB.94: a moon of a free planet has no system page
            item["system_url"] = None
            item["url"] = item_url(row["ref"])
            item["sector_url"] = _sector_cell(row["sector_id"])
        else:
            item["system_url"] = page_url("system", system_id=row["star_system_id"])
            # UX.232: the body's own page by the reference the API gave it (its printed ID), never by the row's number.
            item["url"] = item_url(row["ref"])
            item["sector_url"] = _sector_cell(row["sector_id"])
        if "radius_km" in row:
            item["radius"] = _km(row["radius_km"]) if row["radius_km"] is not None else None
        if "body_type" in row:
            item["body"] = "Gas Giant" if row["body_type"] == "g" else "Terrestrial"
        out.append(item)
    return out


def _link(text, href):
    return {"text": text, "href": href} if href else {"text": text or "—"}


def _name_link(row):
    """The Name cell of a result row (GEN.1140): the link, and a note "formerly X" when the word searched is a name the object used to have."""
    cell = _link(row["name"], row["url"])
    if row.get("formerly"):
        cell["note"] = "(formerly " + ", ".join(row["formerly"]) + ")"
    return cell


def _sector_link(row):
    return _link("View", row["sector_url"]) if row["sector_url"] else {"text": "Standalone"}


def _radius(row):
    return {"text": row["radius"]} if row.get("radius") else {"text": "—"}


def _show_on_map(row):
    return (_link("Show", row["galaxy_url"]) | {"label": f"Show {row['name']} on the Galaxy Map"}
            if row.get("galaxy_url") else {"text": "—"})


SCORE_COLUMNS = ("Suitability", "Microbes", "Flora", "Fauna", "Sentient")
"""tuple: The planet and moon columns that show one score each (UX.214)."""

SORT_KEYS = {"Habitability": "phi4", "Suitability": "suitability", "Microbes": "microbial", "Flora": "flora",
             "Fauna": "fauna", "Sentient": "sentient", "System": "system"}
"""dict: The planet and moon columns that sort, to the API's `planet_sort`/`moon_sort` key (`system` is the page's own order)."""

_PANEL_COLUMNS = {
    "sectors": ["Name", "Edge", "Galaxy Map"],
    "systems": ["Name", "Sector", "Stars", "Objects", "Habitable world", "Star type", "Galaxy Map"],
    "stars": ["Name", "Role", "Type", "Radius", "Values", "System", "Sector"],
    "planets": ["Name", "Class", "Body", "Radius", "Life Class", "Equipment", "Habitability", *SCORE_COLUMNS, "Values",
                "System", "Sector"],
    "moons": ["Name", "Class", "Body", "Radius", "Life Class", "Equipment", "Habitability", *SCORE_COLUMNS, "Values", "Orbits",
              "System", "Sector"],
    "belts": ["Density", "Composition", "System", "Sector"],
    "phenomena": ["Name", "Type", "Class", "Sector"],
}
_PANEL_NOUNS = {
    "sectors": ("sector", "sectors"), "systems": ("system", "systems"), "stars": ("star", "stars"),
    "planets": ("planet", "planets"), "moons": ("moon", "moons"), "belts": ("asteroid belt", "asteroid belts"),
    "phenomena": ("phenomenon", "phenomena"),
}


def _score_cells(row):
    """The Suitability colour and the four life scores of a planet or moon row, one cell each (`SCORE_COLUMNS`); a muted
    dash for a body that was not scored."""
    if row.get("phi4") is None:
        return [{"text": "—", "muted": True}] * len(SCORE_COLUMNS)
    cells = [{"text": TIER_NAMES[row["tier_suitability"]]}]
    for name in LIFE_NAMES:
        score = row.get(f"life_{name}")
        cells.append({"text": "—", "muted": True} if score is None else {"text": f"{score} of 9"})
    return cells


def _values_cell(row):
    """The physical values the search asked for (UX.227), each with its unit, in one cell; a dash when it asked for none."""
    shown = [f"{search_values.KEYS[key].label.split(' (')[0]} {item['value']:.4g}" + (f" {item['symbol']}" if item["symbol"] else "")
             for key, item in (row.get("values") or {}).items() if item.get("value") is not None]
    return {"text": "; ".join(shown)} if shown else {"text": "—", "muted": True}


def _cells(panel, row):
    """One search result row as the data table's cells (the columns of `_PANEL_COLUMNS`)."""
    if panel == "sectors":
        return [_name_link(row), {"text": f"{format_number(row['edge_mpc'], ',.2f')} mpc"},
                _show_on_map(row)]
    if panel == "systems":
        return [_name_link(row), _sector_link(row), {"text": str(row["star_count"])},
                {"text": str(row["object_count"])}, {"text": "Yes" if row["has_habitable_world"] else "No"},
                {"text": row["star_summary"]}, _show_on_map(row)]
    if panel == "stars":
        return [_name_link(row), {"text": row["role"]}, {"text": row["star_type"]}, _radius(row), _values_cell(row),
                _link(row["system_name"], row["system_url"]), _sector_link(row)]
    if panel in ("planets", "moons"):
        cells = [_name_link(row), {"text": row["planet_class"] or "—"}, {"text": row["body"]}, _radius(row),
                 {"text": row["life_class_id"] or "—"},
                 equipment_cell(row),
                 habitability_cell(row, row["url"]), *_score_cells(row), _values_cell(row)]
        if panel == "moons":
            cells.append({"text": row["planet_name"]})
        return cells + [_link(row["system_name"], row["system_url"]), _sector_link(row)]
    if panel == "phenomena":
        return [_link(row["name"], row["url"]), {"text": row["type_label"]},
                {"text": row["phenomenon_class"] or "—"}, _sector_link(row)]
    return [_link(row["density"].capitalize(), row["url"]), {"text": row["composition_summary"]},
            _link(row["system_name"], row["system_url"]), _sector_link(row)]


TIMED_OUT = {"timed_out": True}
"""dict: Stands in for a panel's result when its query hit the statement time limit (UX.97)."""

TIMED_OUT_MESSAGE = "This took too long; try again or narrow the search."
"""str: What a panel that timed out says."""


def api_sorts(args, state=None):
    """The `planet_sort`/`moon_sort` the request's table sorts ask the API for (UX.214): the column's key, `-` in front
    when the table runs best first (descending); the page's own order sends none. Where a table is not sorted, the order the form
    asked for by a physical value (`SearchState.value_sorts`, UX.227) stands, for the stars too."""
    sorts = {}
    for panel, entity in (("planets", "planet"), ("moons", "moon")):
        table_state = datatable.parse_state(PANEL_TABLES[panel], args)
        key = table_state.sort
        if key and key != "system":
            sorts[entity] = ("-" if table_state.descending else "") + key
    for entity, sort in (state or SearchState.from_args(args)).value_sorts().items():
        sorts.setdefault(entity, sort)
    return sorts


def lookup_notice(text):
    """
    The note above the results when the search box holds an Object ID, or a name an object used to have (GEN.1140): where it is now, or that its
    object was deleted. `None` for anything else, a name an object has now (the results below show it) and a lookup that fails (the search
    must not depend on it).

    Returns:
        dict or None: `{"text", "url", "link"}`; `url` is `None` for a deleted object.
    """
    from planetgen.web import item_pages
    text = (text or "").strip()
    if not text:
        return None
    try:
        found = apiclient.get_lookup(db_name(), text)
    except apiclient.ApiError:
        return None
    if not found:
        return None
    if found["status"] == "gone":
        return {"text": f"{text} belonged to an object that has been deleted.", "url": None, "link": None}
    is_id = found.get("id") is not None and found["id"].upper() == text.upper()
    if found["status"] == "current" and not is_id:
        return None
    url = item_pages.item_url(found["ref"])
    if found["status"] == "redirect":
        return {"text": f"{text} is now", "url": url, "link": found["name"] if not is_id else f"{found['name']} ({found['id']})"}
    return {"text": f"{text} is", "url": url, "link": found["name"]}


def fetch(state):
    """
    Everything the Search page needs from the API, one request per result panel (UX.97): the facets and name
    suggestions come from a request that runs no panel, then each panel is asked for on its own, so each runs
    under its own statement time limit and a panel that times out (a 504) is `TIMED_OUT` while the others show.
    Any other API failure raises as before.

    A request that is slow is built by a page job (PERF.105) and raises `slowpages.Pending`. Every request is still
    asked for before the first of them is raised, so all their jobs run side by side and the visitor waits for the
    slowest, not for each in turn (UX.234: answers waited for one at a time went stale before the last was built).

    Returns:
        dict: `get_search`'s answer, its `results` holding every panel.

    Raises:
        slowpages.Pending: At least one request is being built; the first one.
    """
    waiting = []
    try:
        data = apiclient.get_search(db_name(), state.name_terms(), state.tags, sizes=state.sizes(), limit=PAGE_SIZE,
                                    offsets=api_offsets(state), panels=[], sorts=api_sorts(request.args, state), stars=state.stars(),
                                    ranges=state.value_ranges())
    except slowpages.Pending as pending:
        data = None
        waiting.append(pending)
    results = {}
    for panel in PANEL_NAMES:
        try:
            answer = apiclient.get_search(
                db_name(), state.name_terms(), state.tags, sizes=state.sizes(), limit=PAGE_SIZE,
                offsets={panel: api_offsets(state).get(panel, 0)}, panels=[panel], sorts=api_sorts(request.args, state), stars=state.stars(),
                ranges=state.value_ranges())
            results[panel] = answer["results"].get(panel)
        except slowpages.Pending as pending:
            waiting.append(pending)
        except apiclient.ApiError as exc:
            if exc.status_code != 504:
                raise
            results[panel] = TIMED_OUT
    if waiting:
        raise waiting[0]
    data["results"] = results
    return data


def _panel_loader(panel):
    def load(_table_state, limit, offset, _want_facets):
        cached = (getattr(g, "search_results", None) or {}).get(panel)
        if cached is None or cached["offset"] != offset or cached["limit"] != limit:
            state = SearchState.from_args(request.args)
            data = apiclient.get_search(
                db_name(), state.name_terms(), state.tags, sizes=state.sizes(), limit=limit,
                offsets={panel: offset}, panels=[panel], sorts=api_sorts(request.args, state), stars=state.stars(),
                ranges=state.value_ranges())
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
        [Column(SORT_KEYS[label] if panel in ("planets", "moons") and label in SORT_KEYS
                else label.lower().replace(" ", "_"), label,
                sortable=panel in ("planets", "moons") and label in SORT_KEYS)
         for label in _PANEL_COLUMNS[panel]],
        _panel_loader(panel), prefix=f"{panel}_", noun=_PANEL_NOUNS[panel],
        default_sort="system" if panel in ("planets", "moons") else None))
    for panel, heading in RESULT_PANELS
}


def result_panels(state, results):
    """The result panels that ran, each `{"panel", "heading", "total",
    "total_capped", "table"}` (`total_capped`: more matches than `total`,
    shown as "300+"), each a data table (UX.41) that keeps the search and
    every other panel's page. The API's results are handed to the tables
    so the page asks it once."""
    g.search_results = {panel: result for panel, result in results.items() if result is not TIMED_OUT}
    path = url_for("web.search")
    source = {}
    for key, value in state.params():
        source.setdefault(key, []).append(value)
    panels = []
    for panel, heading in RESULT_PANELS:
        result = results.get(panel)
        if result is TIMED_OUT:
            panels.append({"panel": panel, "heading": heading, "timed_out": True, "message": TIMED_OUT_MESSAGE,
                           "total": None, "total_capped": False, "table": None})
            continue
        if result is None or not result["total"]:
            continue  # a group with no matches is left out (UX.54)
        keep = tuple(key for key in _ALL_PARAMS if not key.startswith(f"{panel}_"))
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
    "tag_groups", "active_filters", "result_panels", "api_offsets", "value_form", "VALUE_FIELDS",
]
