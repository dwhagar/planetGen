# stellarObjects/objectId.py

"""
The 64-bit position ID every interstellar object is named by (GEN.64).

Every object in sector space that isn't a generated star system -- rogue
planets, standalone black holes and neutron stars, nebulae, supernova
remnants (whose core is named "<ID> Core"), quasars, interstellar
comets, asteroid fields -- and every star system built around a
bright-sweep star takes this ID, printed as 16 uppercase hex digits, as
its name. It is packed from where the object sits, seen from the
galactic center:

    bits 63-60  type       (4 bits, `KIND_CODES`)
    bits 59-57  unit       (3 bits, `UNITS`: mpc, cpc, pc, kpc, Mpc, Gpc)
    bits 56-40  distance   (17 bits, 0 to 131,071 in that unit)
    bits 39-20  bearing    (20 bits, 0-360 degrees)
    bits 19-0   mark       (20 bits, 0-360 degrees)

Bearing and mark are `navigation.course_between`'s galactic-frame course
from the core to the object, the same "bearing mark mark" the NAV page
uses. The distance takes the smallest unit it fits in, so it is kept to
the hundredth or thousandth of a parsec near the core and to one parsec
out in the disk. At 8 kpc one ID covers about 0.05 x 0.05 x 1 pc.

Two objects can sit closer than that, so an ID can clash. The caller
settles a clash by bumping the ID by one (`bump`) until it is free, in
the order the objects were generated, which keeps it independent of
which worker saved first.
"""

import math

from .navigation import course_between

KIND_CODES = {
    "rogue-planet": 1,
    "black-hole": 2,
    "neutron-star": 3,
    "nebula": 4,
    "supernova-remnant": 5,
    "quasar": 6,
    "comet": 7,
    "asteroid-field": 8,
    "bright-star": 9,
}
"""dict: Object kind (`program_constants.PHENOMENON_TYPE_CHOICES`'
spelling, plus `"bright-star"` for a bright-sweep system) -> the 4-bit
type code. 0 and 10-15 are unused."""

KINDS_BY_CODE = {code: kind for kind, code in KIND_CODES.items()}

UNITS = (("mpc", 1e-3), ("cpc", 1e-2), ("pc", 1.0), ("kpc", 1e3), ("Mpc", 1e6), ("Gpc", 1e9))
"""tuple: `(symbol, parsecs)` per 3-bit unit code, smallest first."""

TYPE_BITS = 4
UNIT_BITS = 3
DISTANCE_BITS = 17
BEARING_BITS = 20
MARK_BITS = 20

_MARK_SHIFT = 0
_BEARING_SHIFT = MARK_BITS
_DISTANCE_SHIFT = _BEARING_SHIFT + BEARING_BITS
_UNIT_SHIFT = _DISTANCE_SHIFT + DISTANCE_BITS
_TYPE_SHIFT = _UNIT_SHIFT + UNIT_BITS

MAX_DISTANCE = (1 << DISTANCE_BITS) - 1
ID_HEX_DIGITS = 16


def _angle_steps(degrees, bits):
    """`degrees` (any value) as a `bits`-bit step count around the circle."""
    steps = 1 << bits
    return int(math.floor((degrees % 360.0) / 360.0 * steps)) % steps


def _distance_fields(distance_pc):
    """`(unit_code, distance)`: the smallest unit `distance_pc` fits in."""
    for code, (_symbol, parsecs) in enumerate(UNITS):
        value = round(distance_pc / parsecs)
        if value <= MAX_DISTANCE:
            return code, value
    raise ValueError(f"distance {distance_pc} pc is too far for an object ID")


def pack(kind, position_pc):
    """
    The object ID for an object of `kind` at galaxy-frame `position_pc`.

    Args:
        kind (str): A `KIND_CODES` key.
        position_pc (tuple): `(x, y, z)` in parsecs, the galactic center
            at the origin.

    Returns:
        int: The 64-bit ID.

    Raises:
        ValueError: An unknown kind, or a position past `Gpc` range.
    """
    if kind not in KIND_CODES:
        raise ValueError(f"no object ID type for {kind!r}")
    course = course_between((0.0, 0.0, 0.0), tuple(float(v) for v in position_pc))
    unit, distance = _distance_fields(course.distance_ly)  # same unit in as out: parsecs here
    return (
        (KIND_CODES[kind] << _TYPE_SHIFT)
        | (unit << _UNIT_SHIFT)
        | (distance << _DISTANCE_SHIFT)
        | (_angle_steps(course.bearing_deg, BEARING_BITS) << _BEARING_SHIFT)
        | (_angle_steps(course.mark_deg, MARK_BITS) << _MARK_SHIFT)
    )


def bump(object_id):
    """The next ID to try after a clash: one mark step on, wrapping
    within the mark field so the type, distance and bearing stay put."""
    mark = (object_id + 1) & ((1 << MARK_BITS) - 1)
    return (object_id & ~((1 << MARK_BITS) - 1)) | mark


def format_id(object_id):
    """`object_id` as its name: 16 uppercase hex digits."""
    return f"{object_id:0{ID_HEX_DIGITS}X}"


def parse_id(name):
    """
    The fields an ID name encodes.

    Args:
        name (str): 16 hex digits, any case.

    Returns:
        dict: `kind`, `unit`, `distance` (in that unit), `distance_pc`,
            `bearing_deg`, `mark_deg` (each the low edge of its step).

    Raises:
        ValueError: Not 16 hex digits, or an unused type or unit code.
    """
    if len(name) != ID_HEX_DIGITS:
        raise ValueError(f"an object ID is {ID_HEX_DIGITS} hex digits, got {name!r}")
    try:
        value = int(name, 16)
    except ValueError:
        raise ValueError(f"not a hex object ID: {name!r}") from None
    code = value >> _TYPE_SHIFT
    unit = (value >> _UNIT_SHIFT) & ((1 << UNIT_BITS) - 1)
    if code not in KINDS_BY_CODE or unit >= len(UNITS):
        raise ValueError(f"not a valid object ID: {name!r}")
    distance = (value >> _DISTANCE_SHIFT) & MAX_DISTANCE
    bearing = (value >> _BEARING_SHIFT) & ((1 << BEARING_BITS) - 1)
    mark = value & ((1 << MARK_BITS) - 1)
    symbol, parsecs = UNITS[unit]
    return {
        "kind": KINDS_BY_CODE[code],
        "unit": symbol,
        "distance": distance,
        "distance_pc": distance * parsecs,
        "bearing_deg": bearing * 360.0 / (1 << BEARING_BITS),
        "mark_deg": mark * 360.0 / (1 << MARK_BITS),
    }


def is_object_id(name):
    """Whether `name` reads as an object ID (`parse_id` accepts it)."""
    try:
        parse_id(name)
    except (TypeError, ValueError):
        return False
    return True
