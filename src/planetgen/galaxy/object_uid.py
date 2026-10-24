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
