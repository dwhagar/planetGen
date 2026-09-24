# html/lib/fmt.py

"""
Small, dependency-free formatting/escaping helpers shared by every script
in `html/` -- what's left of the old `dbutil.py` once its database-access
functions moved out: every page now fetches its data from the planetGen
API (`html/lib/apiclient.py`) instead of querying MySQL directly, so this
module only ever operates on plain values already handed back as JSON,
never a database row or connection.
"""

import html
import json
import os
import re
from urllib.parse import quote

try:
    from stellarObjects.physical_constants import LOCAL_STELLAR_DENSITY_LY3
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # density is still shown, just without the "% of local average"
    # comparison, rather than failing outright (matches every other
    # html/lib module's own fallback for an optional stellarObjects import).
    LOCAL_STELLAR_DENSITY_LY3 = None


def _read_package_version():
    """
    The planetGen package version (`src/stellarObjects/_version.py`'s
    `__version__`, which the post-merge stamp Action updates), read as
    text with a regex the way `setup.py` does -- importing
    `stellarObjects` would pull in `nltk` and friends just for a string.
    Falls back to `"dev"` if the file can't be found or parsed, so a page
    still renders (just without a meaningful cache-busting value).
    """
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                        "stellarObjects", "_version.py")
    try:
        with open(path, encoding="utf-8") as handle:
            match = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', handle.read(), re.MULTILINE)
    except OSError:
        return "dev"
    return match.group(1) if match else "dev"


STATIC_VERSION = _read_package_version()


def static_url(name):
    """
    The URL of a file under `html/static/`, with `?v=<package version>`
    appended -- the one place every `<link>`/`<script src>` in the shell
    and the pages gets its static URL from. The version changes with every
    release, so Apache can tell browsers to cache `static/` for a year
    (`Cache-Control: immutable`, see examples/apache/planetgen.conf.example)
    while an update still reaches everyone on their next page view.

    ES modules that import siblings (`sectormap.js` -> `bodyRendering.js`,
    `vendor/three.module.min.js`) copy their own `?v=` onto those imports
    (`import.meta.url`), so each module has exactly one URL per page.

    Args:
        name (str): Path relative to `static/`, e.g. `"style.css"` or
                    `"vendor/three.module.min.js"`.

    Returns:
        str: e.g. `"static/style.css?v=5.52.0"` (relative, like every
             other link in `html/`).
    """
    return f"static/{name}?v={quote(STATIC_VERSION, safe='')}"


def esc(value):
    """
    HTML-escapes any value for safe interpolation into a page -- API
    content (system/star/planet names, flavor text, generated wikitext)
    is user-influenced-adjacent (from `--name`, `--system-file`, etc.)
    and must never be trusted as pre-sanitized HTML.

    Args:
        value: Any value; `None` becomes `""`.

    Returns:
        str: The escaped string.
    """
    if value is None:
        return ""
    return html.escape(str(value), quote=True)


def post_link(action, params, label, css_class="", attrs=""):
    """
    Builds a same-effect replacement for `<a href="{action}?{a query
    string built from params}">{label}</a>` that carries `params` as
    hidden POST fields instead -- so following it never puts them in the
    browser's own address bar the way a plain link's target URL always
    does. The submit button is styled (`static/style.css`'s `.link-btn`,
    plus whatever `css_class` adds -- `.badge`/`.tag`/`.sidenav-item`/
    etc. all already style a plain class, not specifically an `<a>`, so
    they apply here unchanged) to be visually indistinguishable from the
    `<a>` it replaces. The wrapping `<form>` reuses `.table-form`'s own
    `display: inline` so it flows inline exactly like a link would, inside
    a table cell, breadcrumb, or badge.

    A page reads a followed link's `params` back via `page.nav_params()`/
    `page.nav_multi_params()`, not `page.query_params()` -- see those
    functions.

    Args:
        action (str): The bare target script, e.g. `"system.py"` -- never
                      a query string; that's what `params` replaces.
        params (dict or list[tuple]): `name: value` hidden fields -- every
                       value is escaped here, so callers pass raw values.
                       A `dict` covers every single-valued case (most
                       callers); `search.py`'s own repeated tag-facet
                       fields (several values under the same name, e.g.
                       two `type` filters at once) need a list of `(name,
                       value)` pairs instead, since a `dict` can't hold
                       more than one value per key.
        label (str): Pre-built label HTML (already escaped by the caller,
                     same convention as every other f-string in `html/` --
                     this may legitimately contain markup, e.g. a
                     `<span class="tag-count">`).
        css_class (str): Extra class(es) on the submit button, alongside
                         the shared `link-btn` reset.
        attrs (str): Extra raw HTML attributes on the submit button (e.g.
                     a `title="..."` tooltip).

    Returns:
        str: A complete `<form>` element.
    """
    items = params.items() if isinstance(params, dict) else params
    hidden = "".join(
        f'<input type="hidden" name="{esc(str(k))}" value="{esc(str(v))}">'
        for k, v in items
    )
    klass = esc(f"link-btn {css_class}".strip())
    extra = f" {attrs}" if attrs else ""
    return (
        f'<form method="post" action="{esc(action)}" class="table-form">{hidden}'
        f'<button type="submit" class="{klass}"{extra}>{label}</button></form>'
    )


def data_nav_params(params):
    """
    HTML-attribute-escaped JSON for a `data-nav-params` attribute --
    `static/navform.js`'s delegated click handler reads this (alongside a
    `data-nav-target`) off any element that still has to navigate via a
    real `<a>`/clickable marker rather than `page.post_link`'s own
    `<form>` (a `<form>` can't nest inside an SVG shape the way a Galaxy
    Map/Sector Map/NAV Map marker does -- see `lib/galaxymap.py`/
    `lib/starmap.py`/`lib/navmap.py`), and POSTs it instead of putting it
    in the address bar the way following a plain `href` would.

    Args:
        params (dict): `name: value` hidden fields to post on click.

    Returns:
        str: Escaped JSON, safe to interpolate directly into a
            double-quoted HTML attribute.
    """
    return esc(json.dumps(params, separators=(",", ":")))


_LOCATION_NEIGHBOR_MARKER = " -- nearest: "
_LOCATION_NEIGHBOR_RE = re.compile(r'^(.*) (\([\d.]+ ly\))$')


def linkify_location(db_name, location, name_to_id):
    """
    HTML-escapes a `star_systems.location` string and turns each nearest-
    neighbor name it lists into a link to that system's page.

    `location` is plain text baked in at generation time by
    `stellarObjects._db._format_location_string`, e.g.
    `"Voranthis Kelmoor -- nearest: Alpha Prime (4.2 ly), Beta (5.1 ly)"` --
    the sector name, then up to 3 "Name (distance ly)" entries
    comma-joined after a fixed `" -- nearest: "` marker (empty when the
    sector has no other systems, in which case this is just the sector
    name with nothing to link). This matches that exact format to pull the
    names back out; anything that doesn't fit it (older data predating the
    "-- nearest:" suffix, or a name not found in `name_to_id`) is left as
    plain escaped text rather than guessed at.

    Args:
        db_name (str): The current `?db=` value, for linking to each
                       neighbor's own `system.py` page.
        location (str): The raw `star_systems.location` value.
        name_to_id (dict[str, int]): Every `star_systems.name` -> `id` in
                                     the same sector, for resolving each
                                     neighbor name to a link target.

    Returns:
        str: HTML-safe markup, neighbor names linked where resolvable.
    """
    if not location:
        return ""
    if _LOCATION_NEIGHBOR_MARKER not in location:
        return esc(location)

    prefix, neighbors_part = location.split(_LOCATION_NEIGHBOR_MARKER, 1)
    linked_entries = []
    for entry in neighbors_part.split(", "):
        match = _LOCATION_NEIGHBOR_RE.match(entry)
        name = match.group(1) if match else None
        system_id = name_to_id.get(name) if name is not None else None
        if match and system_id is not None:
            distance = match.group(2)
            link = post_link("system.py", {"db": db_name, "id": system_id}, esc(name))
            linked_entries.append(f'{link} {esc(distance)}')
        else:
            linked_entries.append(esc(entry))

    return f"{esc(prefix)}{_LOCATION_NEIGHBOR_MARKER}" + ", ".join(linked_entries)


def nearest_neighbors_location(db_name, location, neighbors):
    """
    The "Location:" text for a system page, built from live
    `queryDb.system_detail` `nearest_neighbors` data: the sector name
    (the stored `location` string's own prefix) followed by each nearest
    neighbor as a link to its own `system.py` page with its distance.
    Every neighbor here has a real id, so every one is linked, unlike
    `linkify_location`, which can only link names that still match a row.

    Args:
        db_name (str): The current `?db=` value.
        location (str): The raw `star_systems.location` value (only its
                        sector-name prefix is used).
        neighbors (list[dict]): `{id, name, distance_ly}`, nearest first.

    Returns:
        str: HTML-safe markup.
    """
    prefix = (location or "").split(_LOCATION_NEIGHBOR_MARKER, 1)[0]
    entries = [
        f'{post_link("system.py", {"db": db_name, "id": n["id"]}, esc(n["name"]))} ({n["distance_ly"]:.1f} ly)'
        for n in neighbors
    ]
    return f"{esc(prefix)}{_LOCATION_NEIGHBOR_MARKER}" + ", ".join(entries)


def format_distance_ly(distance_ly):
    """Formats a distance in light-years for a table cell, e.g.
    `"26,012.4 ly"`, or an en dash when there is none (`None`, an unplaced
    sector or system)."""
    if distance_ly is None:
        return "&ndash;"
    return f"{distance_ly:,.1f} ly"


def format_density(edge_ly, system_count):
    """
    Formats a sector's star density as systems per cubic light-year, with a
    percentage relative to `LOCAL_STELLAR_DENSITY_LY3` (the real local
    stellar density sector generation targets -- see
    `stellarObjects.spaceSector`'s module docstring) when that comparison
    can be computed.

    Args:
        edge_ly (float): The sector's cube edge, in light-years (the
                         API's `edge_ly` field -- already converted from
                         `sectors.edge_mpc`).
        system_count (int): How many systems are placed in the sector.

    Returns:
        str: e.g. `"0.00329 systems/ly&sup3; (116% of local average)"`.
    """
    if not edge_ly:
        return "n/a"

    density_ly3 = system_count / (edge_ly ** 3)
    text = f"{density_ly3:.5f} systems/ly&sup3;"
    if LOCAL_STELLAR_DENSITY_LY3:
        relative_pct = (density_ly3 / LOCAL_STELLAR_DENSITY_LY3) * 100
        text += f" ({relative_pct:,.0f}% of local average)"
    return text
