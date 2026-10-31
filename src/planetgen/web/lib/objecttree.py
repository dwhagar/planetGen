# planetgen/web/lib/objecttree.py

"""
The summary tree (UX.172): the rows that list a sector's systems, a system's stars, a star's planets and a planet's moons,
each linked to its own page, with one line of what it is and the count of what lies below it.

A row holds the kind's mark, the name (the link to the object's page) and the summary line; a row with children also
holds a toggle button and a `<ul class="tree-children">` of them. Links never sit inside the toggle (a link inside a
button or `<summary>` is a nested control to assistive technology), so the toggle is its own button, hidden until
`static/objecttree.js` runs: without script every level is simply listed. A row whose children are not on the page yet
carries the route that lists them (`data-tree-url`, `GET /object/<ref>/tree`) and a link to its page for a reader
without script; the script replaces that link with the children when the row is opened, a page at a time (a level longer
than a page ends in a "Show more" row), and adds the arrow-key behaviour of the ARIA tree pattern. A row drawn open
(`data-open`) has its children on the page. A row with nothing below it is a plain row with no toggle.

The items are what `GET /api/objects/<ref>/children` returns (`db/children.py`): `ref`, `kind`, `name`, `summary`,
`counts` and `child_count`. Every value from the database is escaped here; the only raw markup is `fmt`'s unit
formatters (built from floats and fixed unit literals), the habitability symbols and the point-of-interest symbol of a body.
"""

import dataclasses

from planetgen.names.bodies import close_pair_label
from planetgen.web.lib.classref import PARABOLIC_COMET_CLASS, class_entry, star_type_classes, star_type_split
from planetgen.web.lib.fmt import distance_km_html, esc, measured_html, period_years_html
from planetgen.web.lib import poi_icons
from planetgen.web.lib.habitability_symbols import symbol_row

KIND_LABELS = {
    "sector": "Sector", "system": "Star system", "star": "Star", "planet": "Planet", "moon": "Moon",
    "belt": "Asteroid belt", "comet": "Comet", "nebula": "Nebula", "asteroid_field": "Asteroid field",
    "black_hole": "Black hole", "neutron_star": "Neutron star", "supernova_remnant": "Supernova remnant",
    "rogue_planet": "Rogue planet", "interstellar_comet": "Interstellar comet", "quasar": "Quasar",
}
"""dict: What a row's kind is called, for the reader of a screen reader and for a phenomenon's summary."""

COUNT_LABELS = (
    ("star", "star", "stars"), ("planet", "planet", "planets"), ("moon", "moon", "moons"),
    ("belt", "asteroid belt", "asteroid belts"), ("comet", "comet", "comets"),
)
"""tuple: The kinds a row counts below it, in the order it names them, with their singular and plural."""

CHILD_NOUNS = {
    "sector": "systems and phenomena", "system": "stars and bodies", "star": "planets, belts and comets",
    "planet": "moons", "rogue_planet": "moons",
}
"""dict: What a parent kind's children are called, for the list's label."""


@dataclasses.dataclass(frozen=True)
class TreeContext:
    """
    What a tree needs from the page that draws it.

    Attributes:
        item_url: A function from a reference (`planet:FE81...`) to the URL of its page.
        tree_url: A function from a reference to the URL of its children (`GET /object/<ref>/tree`).
        class_url: The class-page hook `(type_slug, code) -> url or None`, or `None` for plain text.
        sprite_url: The icon sprite, for the PHI-6 symbols of a scored planet or moon (none drawn without it).
        help_url: The habitability page the symbols link into, or `None`.
        static_url: A function from a path under `static/` to its URL, for the point-of-interest symbol of a body (UX.233; none
            drawn without it).
    """

    item_url: object
    tree_url: object
    class_url: object = None
    sprite_url: str = ""
    help_url: str = None
    static_url: object = None


def _link(url, text_html):
    return f'<a href="{esc(url)}">{text_html}</a>' if url else text_html


def star_type_html(star_type, class_url):
    """A star type as HTML (UX.182): with `class_url`, its spectral part ("G2") links to its spectral class page and its
    luminosity part ("V") to its luminosity class page, the words after them stay plain; unlinked where either page
    is missing or the type does not parse."""
    parts = star_type_split(star_type) if class_url and star_type else None
    classes = star_type_classes(star_type) if parts else None
    if parts is None or classes is None:
        return esc(star_type)
    spectral, luminosity, rest = parts
    spectral_url = class_url("star-spectral", classes[0])
    luminosity_url = class_url("star-luminosity", classes[1])
    return f"{_link(spectral_url, esc(spectral))}{_link(luminosity_url, esc(luminosity))}{esc(rest)}"


def counts_text(counts):
    """What lies below a row as words: "2 stars, 5 planets, 12 moons", or `""` for nothing."""
    parts = [f"{counts[kind]} {one if counts[kind] == 1 else many}" for kind, one, many in COUNT_LABELS
             if counts.get(kind)]
    return ", ".join(parts)


def _stat(text_html):
    return f'<span class="stat">{text_html}</span>' if text_html else ""


def _star_stats(summary, context):
    return [_stat(star_type_html(summary.get("star_type") or "", context.class_url)),
            _stat(esc(summary["role"].capitalize()) if summary.get("role") not in (None, "single", "member") else ""),
            _stat(f"{measured_html('kg', summary['mass_kg'])}" if summary.get("mass_kg") is not None else ""),
            _stat(f"{measured_html('W', summary['luminosity_w'])}" if summary.get("luminosity_w") is not None else "")]


def poi_mark_html(kind, summary, context):
    """The point-of-interest symbol of a row (UX.233) as a `.stat`, or `""` for a row with none or a context with no `static_url`."""
    mark = poi_icons.mark_of_body(kind, summary) if context.static_url else None
    if mark is None:
        return ""
    return f'<span class="stat tree-poi">{poi_icons.symbol_html(mark[0], context.static_url, mark[1])}</span>'


def has_poi_mark(items):
    """Whether any row of `items` (the items of a level) carries a point-of-interest symbol, whatever the context: the page shows its key
    when it does."""
    return any(poi_icons.mark_of_body(item["kind"], item.get("summary") or {}) for item in items)


def _body_stats(summary, context, page=None):
    code = summary.get("planet_class")
    class_html = _link(context.class_url("planet", code) if context.class_url and code else None, esc(code)) if code else ""
    symbols = symbol_row(summary, context.sprite_url, help_url=context.help_url, page=page) \
        if context.sprite_url and summary.get("phi4") is not None else ""
    return [_stat(f"Class {class_html}" if class_html else ""),
            _stat(distance_km_html(summary["distance_km"]) if summary.get("distance_km") is not None else ""),
            _stat(period_years_html(summary["period_years"]) if summary.get("period_years") is not None else ""),
            symbols]


def _belt_stats(summary, context):
    if summary.get("lower_limit_km") is not None and summary.get("upper_limit_km") is not None:
        where = f"{distance_km_html(summary['lower_limit_km'])} to {distance_km_html(summary['upper_limit_km'])}"
    else:
        where = distance_km_html(summary["distance_km"]) if summary.get("distance_km") is not None else ""
    return [_stat(esc(summary["density"]).capitalize() if summary.get("density") else ""), _stat(where)]


def _comet_stats(summary, context):
    elliptical = summary.get("orbit_type") == "elliptical"
    code = summary.get("period_class") if elliptical else PARABOLIC_COMET_CLASS
    url = context.class_url("comet", code) if context.class_url and code else None
    entry = class_entry("comet", code) if code else None
    kind = _link(url, esc(entry["name"])) if url and entry else f"{'Elliptical' if elliptical else 'Parabolic'} comet"
    period = period_years_html(summary["orbital_period_years"]) if summary.get("orbital_period_years") is not None else ""
    return [_stat(kind), _stat("Active" if summary.get("is_active") else "Dormant"), _stat(period)]


def _system_stats(summary, context):
    config = "Binary system" if summary.get("is_binary") else "Single star"
    types = ", ".join(star_type_html(star_type, context.class_url) for star_type in summary.get("star_types") or ())
    return [_stat(config), _stat(types)]


_STATS = {"star": _star_stats, "planet": _body_stats, "moon": _body_stats, "belt": _belt_stats, "comet": _comet_stats,
          "system": _system_stats}


def stats_html(item, context):
    """The row's summary line: its own facts as `.stat` chips, then what lies below it."""
    kind = item["kind"]
    if kind in ("planet", "moon"):    # a body's symbols link to the Habitability group of its own page
        stats = _body_stats(item.get("summary") or {}, context, context.item_url(item["ref"]))
        stats.insert(0, poi_mark_html(kind, item.get("summary") or {}, context))
    else:
        stats = _STATS[kind](item.get("summary") or {}, context) if kind in _STATS else [_stat(esc(KIND_LABELS.get(kind, kind)))]
        if kind == "belt":
            stats.insert(0, poi_mark_html(kind, {}, context))
    if item.get("orbits"):
        stats.insert(0, _stat(f"Orbits {esc(item['orbits'])}"))
    below = counts_text(item.get("counts") or {})
    return "".join(stats) + _stat(esc(below))


def _row_inner(item, context):
    name = esc(item["name"] or KIND_LABELS.get(item["kind"], item["kind"]))
    return (f'<span class="tree-kind tree-kind-{esc(item["kind"])}" aria-hidden="true"></span>'
            f'<a class="tree-name" href="{esc(context.item_url(item["ref"]))}">{name}</a>'
            f'<span class="sr-only"> ({esc(KIND_LABELS.get(item["kind"], item["kind"]).lower())})</span>'
            f'<span class="tree-stats">{stats_html(item, context)}</span>')


def more_html(ref, token):
    """The "Show more" row that ends a level longer than a page."""
    return (f'<li class="tree-more"><button type="button" class="btn btn-small" data-tree-more '
            f'data-next="{esc(token)}">Show more</button></li>')


def node_html(item, context, loaded=None):
    """
    One row. `loaded` is `{ref: (items, next_token)}` for the levels already read (a row listed there is drawn open with
    its children, recursively); any other row with children is drawn closed with the route that lists them.
    """
    loaded = loaded or {}
    ref = esc(item["ref"])
    if not item.get("child_count"):
        return f'<li class="tree-node tree-leaf" data-ref="{ref}"><div class="tree-row">{_row_inner(item, context)}</div></li>'
    noun = CHILD_NOUNS.get(item["kind"], "children")
    label = f'{esc(noun.capitalize())} of {esc(item["name"])}'
    if item["ref"] in loaded:
        children, token = loaded[item["ref"]]
        body = nodes_html(children, context, loaded) + (more_html(item["ref"], token) if token else "")
        state, expanded = " data-open", "true"
    else:
        link = _link(context.item_url(item["ref"]), f'Open its page to see its {esc(noun)}')
        body = f'<li class="tree-pending">{link}</li>'
        state, expanded = "", "false"
    toggle = f'<button type="button" class="tree-toggle" aria-expanded="{expanded}" aria-label="{label}" hidden></button>'
    return (f'<li class="tree-node" data-ref="{ref}" data-tree-url="{esc(context.tree_url(item["ref"]))}"{state}>'
            f'<div class="tree-row">{toggle}{_row_inner(item, context)}</div>'
            f'<ul class="tree-children" aria-label="{label}">{body}</ul></li>')


def _pair_label(items):
    """What the pair of a level that also holds a further star is called ("Abralinar A / B"), or `None` for a level that is not
    that: the stars beyond a close pair sit beside it, and the bodies that go round the pair rather than a star are listed
    after the stars (db/children.py)."""
    stars = {item["summary"].get("role"): item["name"] for item in items if item["kind"] == "star"}
    if "member" not in stars or "primary" not in stars or "secondary" not in stars:
        return None
    return close_pair_label(stars["primary"], stars["secondary"])


def nodes_html(items, context, loaded=None):
    """The rows of one level. On a level that holds a close pair and a star beyond it, the further stars and the bodies that go
    round the whole pair say so, so the rows read as a hierarchy."""
    pair = _pair_label(items)
    if pair:
        items = [{**item, "orbits": pair} if item["kind"] in ("planet", "belt", "comet")
                 or (item["kind"] == "star" and item["summary"].get("role") == "member") else item for item in items]
    return "".join(node_html(item, context, loaded) for item in items)


def tree_html(items, context, label, token=None, ref=None, loaded=None):
    """A whole tree: the list of one level's rows, ending in a "Show more" row when `token` says there are more (`ref` is
    the parent the token pages)."""
    more = more_html(ref, token) if token else ""
    url = f' data-tree-url="{esc(context.tree_url(ref))}"' if ref and token else ""
    return (f'<ul class="object-tree" data-tree-root aria-label="{esc(label)}"{url}>'
            f'{nodes_html(items, context, loaded)}{more}</ul>')
