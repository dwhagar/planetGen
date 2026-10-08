# planetgen/web/lib/fmt.py

"""
The HTML side of formatting, shared by every page: escaping, links,
times, static URLs, and the shared number and unit formatters
(`planetgen.util.format`) with a dash for a missing value. Every page
fetches its data from the planetGen API (`planetgen/web/lib/apiclient.py`),
so this module only ever operates on plain values already handed back as
JSON, never a database row or connection.
"""

import datetime as _dt
import html
import math
import os
import re
from urllib.parse import quote

from planetgen.physics.constants import LOCAL_STELLAR_DENSITY_LY3
from planetgen.util import format as _format


_NON_FINITE_TEXT = re.compile(r"\b(?:nan|inf)\b", re.IGNORECASE)
"""re.Pattern: A "nan"/"inf" a formatter let through (TEST.53)."""


def _finite(value):
    """`value` as a finite float, or `None` for `None`, NaN, an infinity
    or anything `float()` can't take (an int past float range, a string
    that isn't a number)."""
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def dash_unless_finite(formatter, dash="\u2013"):
    """`formatter`, but `dash` for a value that isn't a finite number (or
    whose formatted text would read "nan"/"inf", e.g. 1e300 AU overflowing
    to infinite km) instead of raising or showing "nan km" (TEST.53)."""
    def guarded(value, *args, **kwargs):
        if _finite(value) is None:
            return dash
        text = formatter(value, *args, **kwargs)
        return dash if _NON_FINITE_TEXT.search(text) else text
    guarded.__name__ = formatter.__name__
    guarded.__doc__ = formatter.__doc__
    return guarded


# The shared formatters (`planetgen.util.format`), each giving a dash for
# a missing or non-finite value. A distance is shown on the ladder km < AU
# < mpc < cpc < ly < pc < kpc < Mpc < Gpc, with a parsec value's ly/AU/km
# parenthetical; every page passes its distances through these. Planet,
# moon and star radii are the exception: they are always km in scientific
# notation (`tabledisplay.format_star_radius`).
format_number = dash_unless_finite(_format.format_number)
format_speed_kms = dash_unless_finite(_format.format_speed_kms)
format_duration_seconds = dash_unless_finite(_format.format_duration_seconds)
format_period_years = dash_unless_finite(_format.format_period_years)
format_temperature_k = dash_unless_finite(_format.format_temperature_k)
format_pressure_pa = dash_unless_finite(_format.format_pressure_pa)
format_distance_km = dash_unless_finite(_format.format_distance_km, "&ndash;")
format_distance_au = dash_unless_finite(_format.format_distance_au, "&ndash;")
format_distance_ly = dash_unless_finite(_format.format_distance_ly, "&ndash;")
format_distance_pc = dash_unless_finite(_format.format_distance_pc, "&ndash;")


def _read_package_version():
    """
    The planetGen package version (`src/planetgen/_version.py`'s
    `__version__`, which the post-merge stamp Action updates), read as
    text with a regex the way `setup.py` does -- importing the generator
    modules would pull in `nltk` and friends just for a string.
    Falls back to `"dev"` if the file can't be found or parsed, so a page
    still renders (just without a meaningful cache-busting value).
    """
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                        "_version.py")
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

    ES modules that import siblings (`sectorscene.js` -> `bodyRendering.js`,
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


NEAREST_SHOWN = 3
"""int: How many neighbours a "Nearest:" list shows."""


def linkify_nearest(location, name_to_id, system_url):
    """
    The neighbours listed in a `star_systems.location` string, as HTML with
    each name linked to that system's page, or `""` when it lists none.
    The sector's own name is left out (the breadcrumb carries it, UX.52).

    `location` is plain text baked in at generation time by
    `planetgen.db.store._format_location_string`, e.g.
    `"Voranthis Kelmoor -- nearest: Alpha Vesta (4.2 ly), Beta (5.1 ly)"`:
    the sector name, then up to 3 "Name (distance ly)" entries comma-joined
    after a fixed `" -- nearest: "` marker. This matches that exact format
    to pull the names back out; an entry that doesn't fit it, or a name not
    found in `name_to_id`, is left as plain escaped text rather than
    guessed at. Used where live neighbour data (`nearest_systems_html`) is
    not at hand.

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
    if not location or _LOCATION_NEIGHBOR_MARKER not in location:
        return ""
    neighbors_part = location.split(_LOCATION_NEIGHBOR_MARKER, 1)[1]
    linked_entries = []
    for entry in neighbors_part.split(", ")[:NEAREST_SHOWN]:
        match = _LOCATION_NEIGHBOR_RE.match(entry)
        name = match.group(1) if match else None
        system_id = name_to_id.get(name) if name is not None else None
        if match and system_id is not None:
            distance = match.group(2)
            link = f'<a href="{esc(system_url(system_id))}">{esc(name)}</a>'
            linked_entries.append(f'{link} {esc(distance)}')
        elif entry:
            linked_entries.append(esc(entry))
    return ", ".join(linked_entries)


def nearest_systems_html(neighbors, system_url):
    """`{id, name, distance_ly}` neighbors (`queryDb.nearest_systems`) as
    comma-separated links with their distances (the nearest `NEAREST_SHOWN`).
    HTML-safe markup."""
    def _distance(neighbor):
        distance = _finite(neighbor.get("distance_ly"))
        return f" ({distance:.1f} ly)" if distance is not None else ""
    return ", ".join(f'<a href="{esc(system_url(n["id"]))}">{esc(n["name"])}</a>{_distance(n)}'
                     for n in neighbors[:NEAREST_SHOWN])


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
        try:
            moment = _dt.datetime.fromtimestamp(value, _dt.timezone.utc)
        except (OverflowError, ValueError, OSError):
            return ""  # NaN, an infinity or a time past what a datetime holds (TEST.53)
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
    `planetgen.galaxy.sector`'s module docstring) when that comparison
    can be computed.

    Args:
        edge_ly (float): The sector's cube edge, in light-years (the
                         API's `edge_ly` field -- already converted from
                         `sectors.edge_mpc`).
        system_count (int): How many systems are placed in the sector.

    Returns:
        str: e.g. `"0.00329 systems/ly&sup3; (116% of local average)"`.
    """
    edge_ly, system_count = _finite(edge_ly), _finite(system_count)
    if not edge_ly or edge_ly < 0 or system_count is None:
        return "n/a"

    # Divided one factor at a time rather than by `edge_ly ** 3`, which
    # raises OverflowError for a huge edge and underflows to 0.0
    # (ZeroDivisionError) for a tiny one; an infinite result is n/a too.
    density_ly3 = system_count / edge_ly / edge_ly / edge_ly
    if not math.isfinite(density_ly3):
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
