# planetgen/galaxy/object_uid.py

"""
The object ID (GEN.170, docs/design/architecture/object-id-options.md section 0): one
fixed-length ID for every object in the galaxy, made of where it was born
and which birth it was.

    | birth sector address | serial in that sector | body number |
    | shell, band, wedge   |        28 bits        |   12 bits   |
    |       56 bits        |   (top 2 bits: kind)  |             |

The three sector fields are 19, 18 and 19 bits, so the whole ID is 96 bits:
24 hex digits, `BINARY(12)`, the same for every object in every galaxy
(GEN.230). That holds every cell of the grid out to shell 83,037 (332 kpc
of radius at the 4 pc edge): the outermost shell of the 262 kpc reach
(shell 65,535) needs a band field for bands up to about plus or minus
102,900 and 411,600 wedges in its widest band (`spherical.check_reach`). The 80-bit ID this replaced
was sized for the disk's outline, and its band field could not name a polar
cell beyond shell 1,303.

* The band is stored biased (`band + 2 ** (band_bits - 1)`) so negative
  bands pack as plain integers.
* The serial's top two bits say how it was given out, so no allocator needs
  another's state: `SERIAL_GENERATED` (the object's rank in generation order
  in its birth sector), `SERIAL_RUNTIME` (added after generation: an admin's
  addition, an ejection, a facility) and `SERIAL_FIELD` (an object stored by
  a sector other than the one holding its centre, such as a nebula) and `SERIAL_REGION` (a region of space stored on its
  own, a globular cluster's or a superbubble's container, DB.85: its rank among the regions whose centre lies in the
  birth sector).
* The body number is 0 for the top-level object (a system, a rogue planet, a
  nebula, a facility) and 1 and up for the stars, planets, moons, belts and
  comets born in that system, on one counter.
* The kind of object is not in the ID: a planet ejected from its star becomes
  a rogue planet and keeps its ID.

Printed as `FE81000A2B-0000005-000`: sector, serial and body in hex. Nothing
here reads a database; the allocators that hand out serials are DB.20 and
GEN.171/172.
"""

import re
from typing import NamedTuple

SERIAL_BITS = 28
"""int: Width of the serial field."""

BODY_BITS = 12
"""int: Width of the body-number field."""

SERIAL_KIND_BITS = 2
"""int: How many top bits of the serial say how it was given out."""

SERIAL_GENERATED = 0
SERIAL_RUNTIME = 1
SERIAL_FIELD = 2
SERIAL_REGION = 3
"""int: The serial kinds (the value of the serial's top two bits)."""

SERIAL_KINDS = (SERIAL_GENERATED, SERIAL_RUNTIME, SERIAL_FIELD, SERIAL_REGION)

_TEXT_RE = re.compile(r"^[0-9A-F]+(-[0-9A-F]+){2}$")


class Layout(NamedTuple):
    """The field widths of one galaxy's IDs, in bits."""

    shell_bits: int = 19
    band_bits: int = 18
    wedge_bits: int = 19
    serial_bits: int = SERIAL_BITS
    body_bits: int = BODY_BITS

    @property
    def sector_bits(self):
        return self.shell_bits + self.band_bits + self.wedge_bits

    @property
    def total_bits(self):
        return self.sector_bits + self.serial_bits + self.body_bits

    @property
    def length(self):
        """Bytes in the stored ID (`BINARY(length)`)."""
        return self.total_bits // 8

    @property
    def band_bias(self):
        return 1 << (self.band_bits - 1)

    @property
    def serial_count_bits(self):
        """Bits of the serial left for the number once its kind is taken."""
        return self.serial_bits - SERIAL_KIND_BITS


LAYOUT = Layout()
"""Layout: 19/18/19 + 28 + 12 = 96 bits, every galaxy's (the `uid` columns are `BINARY(12)`)."""


class ObjectId(NamedTuple):
    """An ID taken apart: the birth sector's address, the serial's kind and
    number, and the body number."""

    shell: int
    band: int
    wedge: int
    serial_kind: int
    serial: int
    body: int


def _check(name, value, bits):
    if not isinstance(value, int) or isinstance(value, bool) or not 0 <= value < (1 << bits):
        raise ValueError(f"{name} {value!r} does not fit in {bits} bits")


def check_layout(layout):
    """Raises `ValueError` unless every field is wide enough to use, the
    sector fields and total are whole hex digits and bytes."""
    if min(layout) < 1 or layout.serial_bits <= SERIAL_KIND_BITS or layout.band_bits < 2:
        raise ValueError(f"an object ID layout needs every field wide enough to use: {tuple(layout)}")
    if layout.sector_bits % 4 or layout.serial_bits % 4 or layout.body_bits % 4:
        raise ValueError(f"the sector, serial and body fields must be whole hex digits: {tuple(layout)}")
    if layout.total_bits % 8:
        raise ValueError(f"an object ID is a whole number of bytes: {layout.total_bits}")
    return layout


def pack(shell, band, wedge, serial_kind, serial, body=0, layout=LAYOUT):
    """
    The ID, as an integer, of body `body` of the object whose serial in
    birth sector (`shell`, `band`, `wedge`) is `serial` of kind `serial_kind`.

    Raises:
        ValueError: A field that does not fit its width, or an unknown kind.
    """
    check_layout(layout)
    if serial_kind not in SERIAL_KINDS:
        raise ValueError(f"serial_kind {serial_kind!r} is not one of {SERIAL_KINDS}")
    biased = band + layout.band_bias
    _check("shell", shell, layout.shell_bits)
    _check("band", biased, layout.band_bits)
    _check("wedge", wedge, layout.wedge_bits)
    _check("serial", serial, layout.serial_count_bits)
    _check("body", body, layout.body_bits)
    sector = (((shell << layout.band_bits) | biased) << layout.wedge_bits) | wedge
    number = (serial_kind << layout.serial_count_bits) | serial
    return (((sector << layout.serial_bits) | number) << layout.body_bits) | body


def sector_bytes(shell, band, wedge, layout=LAYOUT):
    """The birth sector's address on its own, as the big-endian bytes the
    sector part of an ID holds (the key of a sector's run-time counter)."""
    return (pack(shell, band, wedge, SERIAL_GENERATED, 0, 0, layout) >> (layout.serial_bits + layout.body_bits)
            ).to_bytes(layout.sector_bits // 8, "big")


def max_body(layout=LAYOUT):
    """The highest body number a system can hold."""
    return (1 << layout.body_bits) - 1


def body_of(value, layout=LAYOUT):
    """The body number of an ID (0 for the top-level object itself)."""
    return value & max_body(layout)


def with_body(value, body, layout=LAYOUT):
    """The ID of body `body` of the system whose ID is `value` (any body of it)."""
    _check("body", body, layout.body_bits)
    return (value & ~max_body(layout)) | body


def unpack(value, layout=LAYOUT):
    """The `ObjectId` an integer from `pack` holds.

    Raises:
        ValueError: A value too big for the layout, or with an unknown
            serial kind.
    """
    check_layout(layout)
    _check("object ID", value, layout.total_bits)
    body = value & ((1 << layout.body_bits) - 1)
    value >>= layout.body_bits
    number = value & ((1 << layout.serial_bits) - 1)
    value >>= layout.serial_bits
    wedge = value & ((1 << layout.wedge_bits) - 1)
    value >>= layout.wedge_bits
    biased = value & ((1 << layout.band_bits) - 1)
    shell = value >> layout.band_bits
    kind = number >> layout.serial_count_bits
    if kind not in SERIAL_KINDS:
        raise ValueError(f"unknown serial kind {kind} in object ID")
    return ObjectId(shell, biased - layout.band_bias, wedge, kind, number & ((1 << layout.serial_count_bits) - 1), body)


def to_bytes(value, layout=LAYOUT):
    """The big-endian bytes a `BINARY(layout.length)` column holds."""
    _check("object ID", value, layout.total_bits)
    return value.to_bytes(layout.length, "big")


def from_bytes(raw, layout=LAYOUT):
    """The integer a `BINARY` column's bytes hold.

    Raises:
        ValueError: Bytes that are not `layout.length` long.
    """
    raw = bytes(raw)
    if len(raw) != layout.length:
        raise ValueError(f"an object ID is {layout.length} bytes, got {len(raw)}")
    return int.from_bytes(raw, "big")


def format_id(value, layout=LAYOUT):
    """The printed form, uppercase hex in three groups: sector, serial and
    body, such as `FE81000A2B-0000005-000`."""
    _check("object ID", value, layout.total_bits)
    body_digits = layout.body_bits // 4
    serial_digits = layout.serial_bits // 4
    sector_digits = layout.sector_bits // 4
    text = f"{value:0{sector_digits + serial_digits + body_digits}X}"
    return f"{text[:sector_digits]}-{text[sector_digits:sector_digits + serial_digits]}-{text[-body_digits:]}"


def parse_id(text, layout=LAYOUT):
    """The integer an ID's printed form names. The hyphens are optional, the
    case is not significant, and the digit count must be the layout's.

    Raises:
        ValueError: Text that is not exactly the layout's hex digits.
    """
    cleaned = str(text).strip().upper()
    digits = layout.total_bits // 4
    if _TEXT_RE.match(cleaned):
        groups = cleaned.split("-")
        want = (layout.sector_bits // 4, layout.serial_bits // 4, layout.body_bits // 4)
        if tuple(len(group) for group in groups) != want:
            raise ValueError(f"an object ID is written {'-'.join('X' * n for n in want)}, got {text!r}")
        cleaned = "".join(groups)
    if len(cleaned) != digits or not re.fullmatch(r"[0-9A-F]+", cleaned):
        raise ValueError(f"an object ID is {digits} hex digits, got {text!r}")
    return int(cleaned, 16)


# --- members of a belt (GEN.1065, docs/design/systems/belt-evolution.md section 9) -----------------------------------------------

MEMBER_SEPARATOR = "."
"""str: Between a belt's printed ID and the rank of one of its members."""
MEMBER_RANK_DIGITS = 4
"""int: The fewest hex digits a printed rank has (zero-padded: `...-000.0042`); a rank past `FFFF` takes as many more as it needs."""
MAX_MEMBER_RANK = 16 ** 8 - 1
"""int: The largest rank a belt member has (eight hex digits)."""


def _check_member(belt_id, rank, layout):
    _check("object ID", belt_id, layout.total_bits)
    if body_of(belt_id, layout) == 0:
        raise ValueError(f"{format_id(belt_id, layout)} is a system's own ID, not a belt's (a belt is a body of the system)")
    if isinstance(rank, bool) or not isinstance(rank, int) or not 1 <= rank <= MAX_MEMBER_RANK:
        raise ValueError(f"a member rank is a whole number from 1 to {MAX_MEMBER_RANK}, got {rank!r}")


def designation(belt_id, rank, layout=LAYOUT):
    """
    The designation of the member of a belt that has `rank` by mass (1 is the heaviest), `<belt ID>.<rank hex>` such as `FE81000A2B-0000005-003.0042`.
    A pure function of the two numbers: no name, no salad name, no codec name, and no draw, so one member always prints one text. The rank is
    the member's place by the mass drawn once at birth, so it stays the same as the belt evolves (a body that is no longer there keeps its number).

    Args:
        belt_id (int): The belt's object ID (a body of its system).
        rank (int): From 1; a rank above the belt's individual count (`is_on_demand`) addresses a body drawn on demand and never stored.

    Raises:
        ValueError: A belt ID that is a system's own, or a rank outside 1 to `MAX_MEMBER_RANK`.
    """
    _check_member(belt_id, rank, layout)
    return f"{format_id(belt_id, layout)}{MEMBER_SEPARATOR}{rank:0{MEMBER_RANK_DIGITS}X}"


def parse_designation(text, layout=LAYOUT):
    """
    `(belt_id, rank)` of a member designation: the belt's ID read by `parse_id` and the rank in hex. The case is not significant; the rank has the digits
    `designation` prints (at least `MEMBER_RANK_DIGITS`, no leading zeros beyond them), so a member has one text.

    Raises:
        ValueError: Text that is not exactly a belt's ID, a point and a canonical rank.
    """
    head, point, tail = str(text).strip().upper().partition(MEMBER_SEPARATOR)
    if not point or MEMBER_SEPARATOR in tail:
        raise ValueError(f"a member designation is <belt ID>{MEMBER_SEPARATOR}<rank hex>, got {text!r}")
    belt_id = parse_id(head, layout)
    if not re.fullmatch(r"[0-9A-F]+", tail):
        raise ValueError(f"a member rank is hex digits, got {tail!r} in {text!r}")
    rank = int(tail, 16)
    _check_member(belt_id, rank, layout)
    if tail != f"{rank:0{MEMBER_RANK_DIGITS}X}":
        raise ValueError(f"the rank of {text!r} is written {rank:0{MEMBER_RANK_DIGITS}X}")
    return belt_id, rank


def is_on_demand(rank, individual_count):
    """Whether the member of this rank is drawn on demand and never stored: its rank is above the number of individually stored bodies of its belt."""
    return rank > individual_count


# --- captured bodies (TEST.495, docs/design/systems/captured-asteroid-moons.md section 4) -------------------------------------------

class Captured(NamedTuple):
    """A captured body's designation and the parent link it gains: `designation` is the text it prints and `captured_by` the printed ID of the planet
    that holds it (the object it is a moon of)."""

    designation: str
    captured_by: str


def same_system(first, second, layout=LAYOUT):
    """Whether two object IDs name bodies of one system: they differ at most in the body number."""
    _check("object ID", first, layout.total_bits)
    _check("object ID", second, layout.total_bits)
    return with_body(first, 0, layout) == with_body(second, 0, layout)


def captured_member(belt_id, rank, planet_id, layout=LAYOUT):
    """
    A member of a belt in the same system that a planet has captured keeps its belt designation, `<belt ID>.<rank hex>`, and gains a parent link to
    the planet (the capture is a flag and a link on the existing object, not a new number): `Captured(designation(belt_id, rank), format_id(planet_id))`.
    Pure; which members are captured is a function of the belt seed, the planet ID and the window, never of insertion order.

    Raises:
        ValueError: A belt ID or rank `designation` refuses, a planet ID that is a system's own ID, a planet of another system, or the belt itself.
    """
    text = designation(belt_id, rank, layout)
    _check_planet(planet_id, layout)
    if not same_system(belt_id, planet_id, layout):
        raise ValueError(f"{format_id(planet_id, layout)} is not in the system of {format_id(belt_id, layout)}")
    if planet_id == belt_id:
        raise ValueError("a belt does not capture its own members")
    return Captured(text, format_id(planet_id, layout))


def outside_reservoir_member(planet_id, rank, layout=LAYOUT):
    """
    A body from a reservoir that is not a belt of the system (the outer disc before it is a belt, a passing population) that a planet holds prints
    `<planet ID>.<rank hex>`, the same sub-number pattern with the planet's ID in the belt's place and the 12-bit body field untouched. Its
    parent link is the planet itself: `Captured(designation, format_id(planet_id))`. Which reservoir it came from is recorded apart.
    """
    _check_planet(planet_id, layout)
    return Captured(designation(planet_id, rank, layout), format_id(planet_id, layout))


def _check_planet(planet_id, layout):
    _check("object ID", planet_id, layout.total_bits)
    if body_of(planet_id, layout) == 0:
        raise ValueError(f"{format_id(planet_id, layout)} is a system's own ID, not a planet's")


def parse_captured(captured, layout=LAYOUT):
    """`(parent_id, rank, planet_id)` of a `Captured`: `parse_designation` of its designation and `parse_id` of its parent link. The inverse of
    `captured_member` and `outside_reservoir_member`; `parent_id == planet_id` for an outside-reservoir body.

    Raises:
        ValueError: A designation or a link that does not parse, or a planet that is not in the system of the parent.
    """
    parent_id, rank = parse_designation(captured.designation, layout)
    planet_id = parse_id(captured.captured_by, layout)
    _check_planet(planet_id, layout)
    if not same_system(parent_id, planet_id, layout):
        raise ValueError(f"{captured.captured_by} is not in the system of {captured.designation}")
    return parent_id, rank, planet_id

