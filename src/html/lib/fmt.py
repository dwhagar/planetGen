# html/lib/fmt.py

"""
Small, dependency-free formatting/escaping helpers shared by every page
in `html/` -- what's left of the old `dbutil.py` once its database-access
functions moved out: every page now fetches its data from the planetGen
API (`html/lib/apiclient.py`) instead of querying MySQL directly, so this
module only ever operates on plain values already handed back as JSON,
never a database row or connection.
"""

import datetime as _dt
import html
import math
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

try:
    from stellarObjects.utils import (
        format_distance_au as _ladder_au,
        format_distance_km as _ladder_km,
        format_distance_ly as _ladder_ly,
        format_distance_pc as _ladder_pc,
        format_duration_seconds,
        format_number,
        format_period_years,
        format_pressure_pa,
        format_speed_kms,
        format_temperature_k,
    )
except ImportError:
    # Without the planetGen package there is no ladder; plain units still
    # read correctly.
    def _ladder_km(km):
        return f"{km:,.0f} km"

    def _ladder_au(au):
        return f"{au:,.3g} AU"

    def _ladder_ly(ly):
        return f"{ly:,.1f} ly"

    def _ladder_pc(pc):
        return f"{pc:,.2f} pc"

    def format_number(value, spec=",.0f"):
        """`stellarObjects.utils.format_number`'s rule, without its module."""
        text = format(value, spec)
        if math.isfinite(value) and ("e" in text or len(text.lstrip("-+").split(".")[0].replace(",", "")) >= 5):
            mantissa, exponent = f"{value:.2e}".split("e")
            superscript = str.maketrans("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")
            return f"{mantissa} \u00d7 10{str(int(exponent)).translate(superscript)}"
        return text

    def format_speed_kms(kms):
        return "\u2013" if kms is None else f"{kms:,.3g} km/s"

    def format_duration_seconds(seconds):
        return "\u2013" if seconds is None else f"{seconds:,.3g} s"

    def format_period_years(years):
        return "\u2013" if years is None else f"{years:,.3g} years"

    def format_temperature_k(kelvin):
        return "\u2013" if kelvin is None else f"{kelvin:,.0f} K"

    def format_pressure_pa(pascals):
        return "\u2013" if pascals is None else f"{pascals:,.3g} Pa"


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


_LOCATION_NEIGHBOR_MARKER = " -- nearest: "
_LOCATION_NEIGHBOR_RE = re.compile(r'^(.*) (\([\d.]+ ly\))$')


def linkify_location(location, name_to_id, system_url):
    """
    HTML-escapes a `star_systems.location` string and turns each nearest-
    neighbor name it lists into a link to that system's page.

    `location` is plain text baked in at generation time by
    `stellarObjects._db._format_location_string`, e.g.
    `"Voranthis Kelmoor -- nearest: Alpha Vesta (4.2 ly), Beta (5.1 ly)"` --
    the sector name, then up to 3 "Name (distance ly)" entries
    comma-joined after a fixed `" -- nearest: "` marker (empty when the
    sector has no other systems, in which case this is just the sector
    name with nothing to link). This matches that exact format to pull the
    names back out; anything that doesn't fit it (older data predating the
    "-- nearest:" suffix, or a name not found in `name_to_id`) is left as
    plain escaped text rather than guessed at.

    Args:
        location (str): The raw `star_systems.location` value.
        name_to_id (dict[str, int]): Every `star_systems.name` -> `id` in
                                     the same sector, for resolving each
                                     neighbor name to a link target.
        system_url (callable): `system_url(system_id)` -> the URL each
                       neighbor links to.

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
            link = f'<a href="{esc(system_url(system_id))}">{esc(name)}</a>'
            linked_entries.append(f'{link} {esc(distance)}')
        else:
            linked_entries.append(esc(entry))

    return f"{esc(prefix)}{_LOCATION_NEIGHBOR_MARKER}" + ", ".join(linked_entries)


def nearest_neighbors_location(location, neighbors, system_url):
    """
    The "Location:" text for a system page, built from live
    `queryDb.system_detail` `nearest_neighbors` data: the sector name
    (the stored `location` string's own prefix) followed by each nearest
    neighbor as a link to its own system page with its distance.
    Every neighbor here has a real id, so every one is linked, unlike
    `linkify_location`, which can only link names that still match a row.

    Args:
        location (str): The raw `star_systems.location` value (only its
                        sector-name prefix is used).
        neighbors (list[dict]): `{id, name, distance_ly}`, nearest first.
        system_url (callable): As for `linkify_location`.

    Returns:
        str: HTML-safe markup.
    """
    prefix = (location or "").split(_LOCATION_NEIGHBOR_MARKER, 1)[0]
    return f"{esc(prefix)}{_LOCATION_NEIGHBOR_MARKER}" + nearest_systems_html(neighbors, system_url)


def nearest_systems_html(neighbors, system_url):
    """`{id, name, distance_ly}` neighbors (`queryDb.nearest_systems`) as
    comma-separated links with their distances. HTML-safe markup."""
    return ", ".join(f'<a href="{esc(system_url(n["id"]))}">{esc(n["name"])}</a> ({n["distance_ly"]:.1f} ly)'
                     for n in neighbors)


def format_distance_km(distance_km):
    """
    A distance (or a nebula, belt or field radius) for display, in the most
    meaningful unit on the ladder km < AU < mpc < cpc < ly < pc < kpc < Mpc
    < Gpc, with a parsec value's ly/AU/km parenthetical: the web face of
    `stellarObjects.utils.format_distance_m`. Every page passes its
    distances through this or its `_au`/`_ly`/`_pc` siblings. `None` (an
    unplaced sector or system) gives an en dash.

    Planet, moon and star radii are the exception: they are always km in
    scientific notation (`tabledisplay.format_star_radius`).
    """
    if distance_km is None:
        return "&ndash;"
    return _ladder_km(distance_km)


def format_distance_au(distance_au):
    """`format_distance_km` for a value in AU."""
    if distance_au is None:
        return "&ndash;"
    return _ladder_au(distance_au)


def format_distance_ly(distance_ly):
    """`format_distance_km` for a value in light-years, e.g. a sector's or
    system's distance from the galactic center ("8 kpc (26,093 ly)")."""
    if distance_ly is None:
        return "&ndash;"
    return _ladder_ly(distance_ly)


def format_distance_pc(distance_pc):
    """`format_distance_km` for a value in parsecs (galaxy geometry)."""
    if distance_pc is None:
        return "&ndash;"
    return _ladder_pc(distance_pc)


def utc_time_html(value):
    """
    A time for a page: `<time datetime="...Z" data-local-time>` whose text
    reads in UTC ("2026-09-30 21:26 UTC"), which `static/localtime.js`
    rewrites in the viewer's own time zone. Without script the page still
    reads correctly, labelled UTC.

    Args:
        value: A Unix time (int/float), a `datetime` (naive ones are UTC,
            the database connection's zone), or an ISO 8601 string (no
            offset means UTC). `None` or `""` gives `""`.

    Returns:
        str: HTML-safe markup, or `""`; an unreadable string comes back
             escaped as-is.
    """
    if value is None or value == "":
        return ""
    if isinstance(value, (int, float)):
        moment = _dt.datetime.fromtimestamp(value, _dt.timezone.utc)
    elif isinstance(value, _dt.datetime):
        moment = value
    else:
        text = str(value).strip()
        try:
            moment = _dt.datetime.fromisoformat(text[:-1] + "+00:00" if text.endswith("Z") else text)
        except ValueError:
            return esc(text)
    if moment.tzinfo is None:
        moment = moment.replace(tzinfo=_dt.timezone.utc)
    moment = moment.astimezone(_dt.timezone.utc)
    return (f'<time datetime="{moment.strftime("%Y-%m-%dT%H:%M:%SZ")}" data-local-time>'
            f'{moment.strftime("%Y-%m-%d %H:%M")} UTC</time>')


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

    # Divided one factor at a time rather than by `edge_ly ** 3`, which
    # raises OverflowError for a huge edge and underflows to 0.0
    # (ZeroDivisionError) for a tiny one; an infinite result is n/a too.
    density_ly3 = system_count / edge_ly / edge_ly / edge_ly
    if math.isinf(density_ly3):
        return "n/a"
    text = f"{format_number(density_ly3, '.5f')} systems/ly&sup3;"
    if LOCAL_STELLAR_DENSITY_LY3:
        relative_pct = (density_ly3 / LOCAL_STELLAR_DENSITY_LY3) * 100
        text += f" ({format_number(relative_pct)}% of local average)"
    return text


def inside_text(row):
    """`"Inside <name>"` for a system or phenomenon that sits in a nebula
    or supernova remnant (`queryDb.containing_cloud`, schema v39), else
    `None`."""
    inside = row.get("inside")
    return f"Inside {inside['name']}" if inside else None


def runaway_text(system):
    """`"Runaway star, 84.4 km/s"` / `"Hypervelocity star, 1.23 Mm/s"` (the
    speed on the shared ladder, `format_speed_kms`) for a system flagged
    fast (`star_systems.runaway_class`, schema v37), else `None`."""
    kind = system.get("runaway_class")
    if not kind:
        return None
    text = f"{kind.capitalize()} star"
    speed = system.get("runaway_speed_kms")
    return f"{text}, {format_speed_kms(speed)}" if speed is not None else text
