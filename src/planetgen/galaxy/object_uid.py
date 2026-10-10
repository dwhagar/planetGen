# planetgen/galaxy/object_uid.py

"""
The object ID (GEN.170, docs/design/object-id-options.md section 0): one
fixed-length ID for every object in the galaxy, made of where it was born
and which birth it was.

    | birth sector address | serial in that sector | body number |
    |  ring, layer, slot   |        28 bits        |   12 bits   |
    |       40 bits        |   (top 2 bits: kind)  |             |

By default the three sector fields are 12, 12 and 16 bits, so the whole ID
is 80 bits: 20 hex digits, `BINARY(10)`. `fit_layout` picks wider sector
fields for a galaxy whose bounds do not fit, and the ID then grows to 96 or
at most 128 bits; the length is the same for every object in one galaxy.

* The layer is stored biased (`layer + 2 ** (layer_bits - 1)`) so negative
  layers pack as plain integers.
* The serial's top two bits say how it was given out, so no allocator needs
  another's state: `SERIAL_GENERATED` (the object's rank in generation order
  in its birth sector), `SERIAL_RUNTIME` (added after generation: an admin's
  addition, an ejection, a facility) and `SERIAL_FIELD` (an object stored by
  a sector other than the one holding its centre, such as a nebula).
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
"""int: The serial kinds (the value of the serial's top two bits)."""

SERIAL_KINDS = (SERIAL_GENERATED, SERIAL_RUNTIME, SERIAL_FIELD)

MAX_BITS = 128
"""int: The longest ID a galaxy may need."""

_TEXT_RE = re.compile(r"^[0-9A-F]+(-[0-9A-F]+){2}$")


class Layout(NamedTuple):
    """The field widths of one galaxy's IDs, in bits."""

    ring_bits: int = 12
    layer_bits: int = 12
    slot_bits: int = 16
    serial_bits: int = SERIAL_BITS
    body_bits: int = BODY_BITS

    @property
    def sector_bits(self):
        return self.ring_bits + self.layer_bits + self.slot_bits

    @property
    def total_bits(self):
        return self.sector_bits + self.serial_bits + self.body_bits

    @property
    def length(self):
        """Bytes in the stored ID (`BINARY(length)`)."""
        return self.total_bits // 8

    @property
    def layer_bias(self):
        return 1 << (self.layer_bits - 1)

    @property
    def serial_count_bits(self):
        """Bits of the serial left for the number once its kind is taken."""
        return self.serial_bits - SERIAL_KIND_BITS


DEFAULT_LAYOUT = Layout()
"""Layout: 12/12/16 + 28 + 12 = 80 bits, the default Milky Way's."""


class ObjectId(NamedTuple):
    """An ID taken apart: the birth sector's address, the serial's kind and
    number, and the body number."""

    ring: int
    layer: int
    slot: int
    serial_kind: int
    serial: int
    body: int


def _check(name, value, bits):
    if not isinstance(value, int) or isinstance(value, bool) or not 0 <= value < (1 << bits):
        raise ValueError(f"{name} {value!r} does not fit in {bits} bits")


def check_layout(layout):
    """Raises `ValueError` unless every field is wide enough to use, the
    sector fields and total are whole hex digits and bytes, and the whole is
    at most `MAX_BITS`."""
    if min(layout) < 1 or layout.serial_bits <= SERIAL_KIND_BITS or layout.layer_bits < 2:
        raise ValueError(f"an object ID layout needs every field wide enough to use: {tuple(layout)}")
    if layout.sector_bits % 4 or layout.serial_bits % 4 or layout.body_bits % 4:
        raise ValueError(f"the sector, serial and body fields must be whole hex digits: {tuple(layout)}")
    if layout.total_bits % 8 or layout.total_bits > MAX_BITS:
        raise ValueError(f"an object ID is a whole number of bytes, at most {MAX_BITS} bits: {layout.total_bits}")
    return layout


def pack(ring, layer, slot, serial_kind, serial, body=0, layout=DEFAULT_LAYOUT):
    """
    The ID, as an integer, of body `body` of the object whose serial in
    birth sector (`ring`, `layer`, `slot`) is `serial` of kind `serial_kind`.

    Raises:
        ValueError: A field that does not fit its width, or an unknown kind.
    """
    check_layout(layout)
    if serial_kind not in SERIAL_KINDS:
        raise ValueError(f"serial_kind {serial_kind!r} is not one of {SERIAL_KINDS}")
    biased = layer + layout.layer_bias
    _check("ring", ring, layout.ring_bits)
    _check("layer", biased, layout.layer_bits)
    _check("slot", slot, layout.slot_bits)
    _check("serial", serial, layout.serial_count_bits)
    _check("body", body, layout.body_bits)
    sector = (((ring << layout.layer_bits) | biased) << layout.slot_bits) | slot
    number = (serial_kind << layout.serial_count_bits) | serial
    return (((sector << layout.serial_bits) | number) << layout.body_bits) | body


def sector_bytes(ring, layer, slot, layout=DEFAULT_LAYOUT):
    """The birth sector's address on its own, as the big-endian bytes the
    sector part of an ID holds (the key of a sector's run-time counter)."""
    return (pack(ring, layer, slot, SERIAL_GENERATED, 0, 0, layout) >> (layout.serial_bits + layout.body_bits)
            ).to_bytes(layout.sector_bits // 8, "big")


def max_body(layout=DEFAULT_LAYOUT):
    """The highest body number a system can hold."""
    return (1 << layout.body_bits) - 1


def body_of(value, layout=DEFAULT_LAYOUT):
    """The body number of an ID (0 for the top-level object itself)."""
    return value & max_body(layout)


def with_body(value, body, layout=DEFAULT_LAYOUT):
    """The ID of body `body` of the system whose ID is `value` (any body of it)."""
    _check("body", body, layout.body_bits)
    return (value & ~max_body(layout)) | body


def unpack(value, layout=DEFAULT_LAYOUT):
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
    slot = value & ((1 << layout.slot_bits) - 1)
    value >>= layout.slot_bits
    biased = value & ((1 << layout.layer_bits) - 1)
    ring = value >> layout.layer_bits
    kind = number >> layout.serial_count_bits
    if kind not in SERIAL_KINDS:
        raise ValueError(f"unknown serial kind {kind} in object ID")
    return ObjectId(ring, biased - layout.layer_bias, slot, kind, number & ((1 << layout.serial_count_bits) - 1), body)


def to_bytes(value, layout=DEFAULT_LAYOUT):
    """The big-endian bytes a `BINARY(layout.length)` column holds."""
    _check("object ID", value, layout.total_bits)
    return value.to_bytes(layout.length, "big")


def from_bytes(raw, layout=DEFAULT_LAYOUT):
    """The integer a `BINARY` column's bytes hold.

    Raises:
        ValueError: Bytes that are not `layout.length` long.
    """
    raw = bytes(raw)
    if len(raw) != layout.length:
        raise ValueError(f"an object ID is {layout.length} bytes, got {len(raw)}")
    return int.from_bytes(raw, "big")


def format_id(value, layout=DEFAULT_LAYOUT):
    """The printed form, uppercase hex in three groups: sector, serial and
    body, such as `FE81000A2B-0000005-000`."""
    _check("object ID", value, layout.total_bits)
    body_digits = layout.body_bits // 4
    serial_digits = layout.serial_bits // 4
    sector_digits = layout.sector_bits // 4
    text = f"{value:0{sector_digits + serial_digits + body_digits}X}"
    return f"{text[:sector_digits]}-{text[sector_digits:sector_digits + serial_digits]}-{text[-body_digits:]}"


def parse_id(text, layout=DEFAULT_LAYOUT):
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


def _bits_for(count):
    """Bits to hold the values 0 .. count - 1 (at least 1)."""
    return max(1, (count - 1).bit_length())


def fit_layout(max_ring, min_layer, max_layer, max_slots):
    """
    The layout for a galaxy with rings 0..`max_ring`, layers `min_layer` to
    `max_layer` and at most `max_slots` slots in a ring: the default 80 bits
    when its bounds fit, else the sector fields widened (the layer field to
    whole hex digits, the rest to fill the nibbles) until the ID is 96 or 128
    bits.

    Raises:
        ValueError: Bounds that do not fit 128 bits.
    """
    ring_bits = _bits_for(max_ring + 1)
    slot_bits = _bits_for(max_slots)
    reach = max(abs(min_layer), abs(max_layer) + 1)
    layer_bits = max(2, reach.bit_length() + 1)
    if (ring_bits <= DEFAULT_LAYOUT.ring_bits and layer_bits <= DEFAULT_LAYOUT.layer_bits
            and slot_bits <= DEFAULT_LAYOUT.slot_bits):
        return DEFAULT_LAYOUT
    need = ring_bits + layer_bits + slot_bits
    for total in (96, 128):
        sector_bits = total - SERIAL_BITS - BODY_BITS
        if need <= sector_bits:
            spare = sector_bits - need
            # Spare bits go to the ring first, then the slot, keeping the sum.
            ring_bits += spare
            if (ring_bits + layer_bits + slot_bits) % 4:
                raise ValueError("internal: sector fields not nibble aligned")
            return check_layout(Layout(ring_bits, layer_bits, slot_bits))
    raise ValueError(f"a galaxy of {max_ring + 1} rings, layers {min_layer} to {max_layer} and {max_slots} slots "
                     f"does not fit a {MAX_BITS}-bit object ID")
