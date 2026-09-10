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
import re

try:
    from stellarObjects.physical_constants import LOCAL_STELLAR_DENSITY_LY3
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # density is still shown, just without the "% of local average"
    # comparison, rather than failing outright (matches every other
    # html/lib module's own fallback for an optional stellarObjects import).
    LOCAL_STELLAR_DENSITY_LY3 = None


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
        db_name (str): The current `?db=` value, for building system.py URLs.
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
            linked_entries.append(
                f'<a href="system.py?db={esc(db_name)}&id={system_id}">{esc(name)}</a> {esc(distance)}'
            )
        else:
            linked_entries.append(esc(entry))

    return f"{esc(prefix)}{_LOCATION_NEIGHBOR_MARKER}" + ", ".join(linked_entries)


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
