# planetgen/galaxy/objectref.py

"""
One reference form for every object (NAV.7): `<kind>:<id>`.

Kinds: `sector`, `system`, `star`, `planet`, `moon`, `belt`, `comet`, and the
standalone phenomenon types (`nebula`, `black_hole`, ...). `id` is the
object's own row id. A bare number means a system, so the `system:12` that
`/nav` and the Sector Map already write is the same reference. The galaxy
itself is the root of every parent chain and has the reference `galaxy`.

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
