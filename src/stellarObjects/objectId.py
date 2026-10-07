# stellarObjects/objectId.py

"""
The 76-bit position ID every interstellar object is named by (GEN.64).

Every object in sector space that isn't a generated star system -- rogue
planets, standalone black holes and neutron stars, nebulae, supernova
remnants and their collapsed cores, quasars, interstellar comets,
asteroid fields -- and every star system built around a bright-sweep
star takes this ID, printed as 19 uppercase hex digits, as its name. It
is packed from where the object sits, seen from the galactic center:

    bits 75-70  type       (6 bits, `KIND_CODES`)
    bits 69-67  unit       (3 bits, `UNITS`: mpc, cpc, pc, kpc, Mpc, Gpc)
    bits 66-48  distance   (19 bits, 0 to 524,287 in that unit)
    bits 47-26  bearing    (22 bits, 0-360 degrees)
    bits 25-4   mark       (22 bits, 0-360 degrees)
    bits 3-0    collision  (4 bits, 0-15)

Bearing and mark are `navigation.course_between`'s galactic-frame course
from the core to the object, the same "bearing mark mark" the NAV page
uses. The distance takes the smallest unit it fits in, so it is kept to
the thousandth or hundredth of a parsec out to 5.2 kpc and to one parsec
beyond. At 8 kpc one position covers about 0.012 x 0.012 x 1 pc.

Two objects of one type can sit closer than that. The collision number
tells up to sixteen of them apart, in the order they were generated,
which keeps it independent of which worker saved first (`bump`). A
seventeenth moves on one mark step.
"""

import math

from planetgen.galaxy.navigation import course_between

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
    "black-hole-core": 10,
    "neutron-star-core": 11,
}
"""dict: Object kind (`program_constants.PHENOMENON_TYPE_CHOICES`'
spelling, plus `"bright-star"` for a bright-sweep system and the two
`-core` kinds for a supernova remnant's collapsed core) -> the 6-bit
type code. 0 and 12-63 are unused."""

KINDS_BY_CODE = {code: kind for kind, code in KIND_CODES.items()}

UNITS = (("mpc", 1e-3), ("cpc", 1e-2), ("pc", 1.0), ("kpc", 1e3), ("Mpc", 1e6), ("Gpc", 1e9))
"""tuple: `(symbol, parsecs)` per 3-bit unit code, smallest first."""

TYPE_BITS = 6
UNIT_BITS = 3
DISTANCE_BITS = 19
BEARING_BITS = 22
MARK_BITS = 22
COLLISION_BITS = 4

_COLLISION_SHIFT = 0
_MARK_SHIFT = COLLISION_BITS
_BEARING_SHIFT = _MARK_SHIFT + MARK_BITS
_DISTANCE_SHIFT = _BEARING_SHIFT + BEARING_BITS
_UNIT_SHIFT = _DISTANCE_SHIFT + DISTANCE_BITS
_TYPE_SHIFT = _UNIT_SHIFT + UNIT_BITS
ID_BITS = _TYPE_SHIFT + TYPE_BITS

MAX_DISTANCE = (1 << DISTANCE_BITS) - 1
MAX_COLLISION = (1 << COLLISION_BITS) - 1
ID_HEX_DIGITS = (ID_BITS + 3) // 4


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
        int: The 76-bit ID, collision number 0.

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


def collision(object_id):
    """The collision number (0-15) of `object_id`."""
    return (object_id >> _COLLISION_SHIFT) & MAX_COLLISION


def bump(object_id):
    """
    The next ID to try after a clash: the next collision number, or, past
    the sixteenth, one mark step on with the collision number back at 0
    (wrapping within the mark field, so the type, distance and bearing
    stay put).
    """
    if collision(object_id) < MAX_COLLISION:
        return object_id + 1
    mark_mask = ((1 << MARK_BITS) - 1) << _MARK_SHIFT
    mark = (((object_id & mark_mask) >> _MARK_SHIFT) + 1) & ((1 << MARK_BITS) - 1)
    return (object_id & ~mark_mask & ~MAX_COLLISION) | (mark << _MARK_SHIFT)


def format_id(object_id):
    """`object_id` as its name: 19 uppercase hex digits."""
    return f"{object_id:0{ID_HEX_DIGITS}X}"


def parse_id(name):
    """
    The fields an ID name encodes.

    Args:
        name (str): 19 hex digits, any case.

    Returns:
        dict: `kind`, `unit`, `distance` (in that unit), `distance_pc`,
            `bearing_deg`, `mark_deg` (each the low edge of its step),
            `collision`.

    Raises:
        ValueError: Not 19 hex digits, or an unused type or unit code.
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
    mark = (value >> _MARK_SHIFT) & ((1 << MARK_BITS) - 1)
    symbol, parsecs = UNITS[unit]
    return {
        "kind": KINDS_BY_CODE[code],
        "unit": symbol,
        "distance": distance,
        "distance_pc": distance * parsecs,
        "bearing_deg": bearing * 360.0 / (1 << BEARING_BITS),
        "mark_deg": mark * 360.0 / (1 << MARK_BITS),
        "collision": collision(value),
    }


def is_object_id(name):
    """Whether `name` reads as an object ID (`parse_id` accepts it)."""
    try:
        parse_id(name)
    except (TypeError, ValueError):
        return False
    return True
