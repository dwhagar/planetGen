# planetgen/galaxy/objectref.py

"""
One reference form for every object (NAV.7): `<kind>:<id>`.

Kinds: `sector`, `system`, `star`, `planet`, `moon`, `belt`, `comet`, and the
standalone phenomenon types (`nebula`, `black_hole`, ...). `id` is the
object's own row id inside the code (`parse`, `format`) and its printed object
ID in everything a person or a client sees (`parse_public`, `format_public`, API.23):
`system:FE81000A2B-0000005-000`, `sector:FE81000A2B`. A bare number means a system
(row id), a bare three-part ID is a system's too. The galaxy itself is the root of
every parent chain and has the reference `galaxy`.

`static/objectref.js` is the JavaScript twin of `parse` and `format`.
"""

import re

GALAXY = "galaxy"
"""str: The reference of the galaxy, the root of every parent chain."""

PHENOMENON_KINDS = (
    "nebula",
    "asteroid_field",
    "black_hole",
    "neutron_star",
    "supernova_remnant",
    "rogue_planet",
    "interstellar_comet",
    "quasar",
)
"""tuple: The standalone phenomenon kinds, matching the `type` values
`queryDb.list_phenomena` returns."""

BODY_KINDS = ("star", "planet", "moon", "belt", "comet")
"""tuple: The kinds that live inside a system."""

KINDS = ("sector", "system") + BODY_KINDS + PHENOMENON_KINDS
"""tuple: Every kind a reference can name (the galaxy has no id)."""

TABLES = {
    "sector": "sectors",
    "system": "star_systems",
    "star": "stars",
    "planet": "planets",
    "moon": "moons",
    "belt": "asteroid_belts",
    "comet": "comets",
}
"""dict: Each non-phenomenon kind's table (phenomena use
`queryDb._PHENOMENON_TYPE_TO_TABLE`)."""

_REF_RE = re.compile(r"^(?:([a-z_]+):)?(\d+)$")
_PUBLIC_RE = re.compile(r"^(?:([a-z_]+):)?([0-9A-Fa-f]+(?:-[0-9A-Fa-f]+){0,2})$")


def format(kind, object_id):
    """
    The reference for an object: `format("moon", 7)` is `"moon:7"`.

    Raises:
        ValueError: For an unknown kind or a negative id.
    """
    if kind not in KINDS:
        raise ValueError(f"unknown object kind: {kind!r}")
    object_id = int(object_id)
    if object_id < 0:
        raise ValueError(f"object id must not be negative: {object_id}")
    return f"{kind}:{object_id}"


def parse(raw):
    """
    Parses a reference.

    Returns:
        tuple: `(kind, id)`; a bare number is `("system", id)`.

    Raises:
        ValueError: For anything that is not `<kind>:<id>` with a known kind.
    """
    match = _REF_RE.match(str(raw).strip())
    if not match or match.group(1) not in (None, *KINDS):
        raise ValueError(f"not an object reference: {raw!r}")
    return (match.group(1) or "system"), int(match.group(2))


def parse_public(raw):
    """
    Parses a reference as a person or a client writes it: `<kind>:<printed ID>` (a bare three-part object ID is a
    system's). The ID is returned as the upper-case text it was written in, never as a number: a sector's printed
    ID can be all decimal digits.

    Returns:
        tuple: `(kind, printed id)`.

    Raises:
        ValueError: For anything that is not `<kind>:<id>` with a known kind.
    """
    match = _PUBLIC_RE.match(str(raw).strip())
    if not match or match.group(1) not in (None, *KINDS) or (match.group(1) is None and "-" not in match.group(2)):
        raise ValueError(f"not an object reference: {raw!r}")
    return (match.group(1) or "system"), match.group(2).upper()


def format_public(kind, printed_id):
    """The reference for an object by its printed ID: `format_public("moon", "FE81000A2B-0000005-003")`."""
    if kind not in KINDS:
        raise ValueError(f"unknown object kind: {kind!r}")
    if not _PUBLIC_RE.match(str(printed_id)):
        raise ValueError(f"not an object ID: {printed_id!r}")
    return f"{kind}:{str(printed_id).upper()}"

