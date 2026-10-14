"""
The packed form of a Galaxy Map tile (MAP.159): what the tile cache stores and
what `/galaxy/tiles` sends.

A tile's stars are most of its bytes, and as JSON each costs about 150 bytes
(28 gzipped). Packed, a star is 19 bytes in columns (about 12 gzipped), which
is what the 78% cut in docs/design/galaxy-map-wire-format.md measures. Clouds,
points and the `incomplete` list stay JSON: they are few and carry text.

One tile (all numbers little-endian)::

    "PGT1"                        4 bytes, the format
    meta_length                   u32
    meta                          UTF-8 JSON: the sections that are not stars
                                  (`clouds`, `points`, `incomplete`), or `{}`
    then `stars`, then `generated`, each:
    count                         u32
    if count > 0:
      x, y, z                     3 x (f64 minimum, f64 step): a coordinate is
                                  minimum + code * step
      x, y, z codes               3 x u16[count]
      luminosity                  u16[count]: 0 is "none or zero", else 1..65535
                                  spread over log10 LUMINOSITY_LOG_RANGE
      temperature                 u8[count]: 0 is "none", else 1..255 over
                                  log10 TEMPERATURE_LOG_RANGE
      radius                      u8[count]: 0 is "none", else 1..255 over
                                  log10 RADIUS_LOG_RANGE
      class                       u8[count]: the first letter of the star
                                  type, 0 for none
      charted                     u8[count]: 1 when the star's system is
                                  known to the database (`system_id`)
      id kind                     u8: 0 whole numbers, 1 text of one length,
                                  2 anything else
        kind 0                    u32[count] low halves, then u32[count]
                                  high halves (ids are 64 bit)
        kind 1                    u8 length, then count x that many ASCII
                                  bytes (a printed object ID is 22)
        kind 2                    u32 length, then a JSON list of the ids

A star keeps its place in the list (the list order is its rank, MAP.154). Many
tiles in one response are joined by `pack_response`.
"""

import json
import math
import struct

MIMETYPE = "application/vnd.planetgen.tiles"
"""str: What `/galaxy/tiles` is sent as (the web server's compression lists name it)."""

MAGIC = b"PGT1"
RESPONSE_MAGIC = b"PGTB"

# What a code means. The ends are chosen far outside what the generator
# makes (a brown dwarf is 1e-5 L, the brightest star here 1e6; 300 K to
# 500,000 K; a hundredth of a solar radius to 10,000), so a clamp is a
# corrupt record, never a real one.
LUMINOSITY_LOG_RANGE = (-8.0, 8.0)
TEMPERATURE_LOG_RANGE = (2.0, 6.0)
RADIUS_LOG_RANGE = (-3.0, 4.0)

# The sections that stay JSON.
_META_SECTIONS = ("clouds", "points", "incomplete")
_STAR_SECTIONS = ("stars", "generated")
_U16 = 65535
_U8 = 255
_U32 = 0xFFFFFFFF


class TilePackError(ValueError):
    """Bytes that are not a packed tile."""


def _encode_log(value, lo, hi, top):
    """`value` as a code 1..`top` over log10 `lo..hi`, or 0 for none."""
    if value is None or value <= 0:
        return 0
    position = (math.log10(value) - lo) / (hi - lo)
    return 1 + min(top - 1, max(0, int(round(position * (top - 1)))))


def _decode_log(code, lo, hi, top):
    """The value a code stands for; None for code 0."""
    if code == 0:
        return None
    return 10.0 ** (lo + (code - 1) / (top - 1) * (hi - lo))


def _axis(values):
    """`(minimum, step)` that code `values` into u16 steps."""
    low, high = min(values), max(values)
    return low, (high - low) / _U16 if high > low else 0.0


def _codes(values, low, step):
    if step == 0.0:
        return [0] * len(values)
    return [min(_U16, max(0, int(round((value - low) / step)))) for value in values]


def _pack_stars(stars):
    count = len(stars)
    out = [struct.pack("<I", count)]
    if not count:
        return out
    axes = []
    planes = []
    for axis_name in ("x", "y", "z"):
        values = [float(star[axis_name]) for star in stars]
        low, step = _axis(values)
        axes.append(struct.pack("<dd", low, step))
        planes.append(struct.pack(f"<{count}H", *_codes(values, low, step)))
    out.extend(axes)
    out.extend(planes)

    luminosity_lo, luminosity_hi = LUMINOSITY_LOG_RANGE
    out.append(struct.pack(f"<{count}H", *(
        _encode_log(star.get("luminosity_sol"), luminosity_lo, luminosity_hi, _U16) for star in stars)))
    temperature_lo, temperature_hi = TEMPERATURE_LOG_RANGE
    out.append(bytes(
        _encode_log(star.get("temperature_k"), temperature_lo, temperature_hi, _U8) for star in stars))
    radius_lo, radius_hi = RADIUS_LOG_RANGE
    out.append(bytes(_encode_log(star.get("radius_sol"), radius_lo, radius_hi, _U8) for star in stars))
    out.append(bytes(ord(str(star["star_type"])[0]) if star.get("star_type") else 0 for star in stars))
    out.append(bytes(1 if star.get("system_id") is not None else 0 for star in stars))
    out.extend(_pack_ids([star["id"] for star in stars]))
    return out


def _pack_ids(ids):
    count = len(ids)
    if all(isinstance(i, int) and not isinstance(i, bool) and 0 <= i < 1 << 64 for i in ids):
        return [b"\x00", struct.pack(f"<{count}I", *(i & _U32 for i in ids)),
                struct.pack(f"<{count}I", *(i >> 32 for i in ids))]
    if all(isinstance(i, str) and i.isascii() for i in ids) and len({len(i) for i in ids}) == 1 and 0 < len(ids[0]) < 256:
        return [b"\x01", bytes([len(ids[0])]), "".join(ids).encode("ascii")]
    text = json.dumps(ids, separators=(",", ":")).encode("utf-8")
    return [b"\x02", struct.pack("<I", len(text)), text]


def pack_tile(tile):
    """
    A trimmed tile (`tilecache.trim_tile`'s dict) as bytes.

    Args:
        tile (dict): `stars` and `generated` (lists of star dicts with
            `id`, `x`, `y`, `z`, `luminosity_sol` and optionally
            `temperature_k`, `radius_sol`, `star_type`, `system_id`) and the
            JSON sections `clouds`, `points` and `incomplete`.

    Returns:
        bytes: The packed tile.
    """
    meta = {name: tile[name] for name in _META_SECTIONS if tile.get(name) is not None}
    meta_bytes = json.dumps(meta, separators=(",", ":")).encode("utf-8")
    out = [MAGIC, struct.pack("<I", len(meta_bytes)), meta_bytes]
    for section in _STAR_SECTIONS:
        out.extend(_pack_stars(tile.get(section) or []))
    return b"".join(out)


def is_incomplete(packed):
    """Whether a packed tile says some piece of it is missing (PERF.76): it is served but
    not kept."""
    try:
        length = struct.unpack_from("<I", packed, 4)[0]
        return b'"incomplete":[' in packed[8:8 + length]
    except struct.error:
        return False


_BYTES_PER_STAR = 6 + 2 + 1 + 1 + 1 + 1
_AXES_BYTES = 3 * 16


def _stars_end(raw, offset):
    """Where the star section that starts at `offset` ends, from its lengths alone."""
    (count,) = struct.unpack_from("<I", raw, offset)
    offset += 4
    if not count:
        return offset
    offset += _AXES_BYTES + count * _BYTES_PER_STAR
    kind = raw[offset]
    offset += 1
    if kind == 0:
        return offset + 8 * count
    if kind == 1:
        return offset + 1 + raw[offset] * count
    if kind == 2:
        (length,) = struct.unpack_from("<I", raw, offset)
        return offset + 4 + length
    raise struct.error("unknown id kind")


def is_whole(raw):
    """Whether `raw` is exactly one packed tile, found from its lengths alone (a torn or
    foreign file is not), without decoding a star."""
    if raw is None or raw[:4] != MAGIC or len(raw) < 8:
        return False
    try:
        offset = 8 + struct.unpack_from("<I", raw, 4)[0]
        for _ in _STAR_SECTIONS:
            offset = _stars_end(raw, offset)
    except (struct.error, IndexError):
        return False
    return offset == len(raw)


def _unpack_stars(buffer, offset):
    (count,) = struct.unpack_from("<I", buffer, offset)
    offset += 4
    if not count:
        return [], offset
    axes = []
    for _ in range(3):
        axes.append(struct.unpack_from("<dd", buffer, offset))
        offset += 16
    positions = []
    for _ in range(3):
        positions.append(struct.unpack_from(f"<{count}H", buffer, offset))
        offset += 2 * count
    luminosity = struct.unpack_from(f"<{count}H", buffer, offset)
    offset += 2 * count
    temperature = buffer[offset:offset + count]
    offset += count
    radius = buffer[offset:offset + count]
    offset += count
    classes = buffer[offset:offset + count]
    offset += count
    charted = buffer[offset:offset + count]
    offset += count
    kind = buffer[offset]
    offset += 1
    if kind == 0:
        low_ids = struct.unpack_from(f"<{count}I", buffer, offset)
        offset += 4 * count
        high_ids = struct.unpack_from(f"<{count}I", buffer, offset)
        offset += 4 * count
        ids = [low | (high << 32) for low, high in zip(low_ids, high_ids)]
    elif kind == 1:
        width = buffer[offset]
        offset += 1
        text = buffer[offset:offset + width * count].decode("ascii")
        offset += width * count
        ids = [text[i * width:(i + 1) * width] for i in range(count)]
    elif kind == 2:
        (length,) = struct.unpack_from("<I", buffer, offset)
        offset += 4
        ids = json.loads(buffer[offset:offset + length].decode("utf-8"))
        offset += length
    else:
        raise TilePackError("a packed tile has an unknown id kind")
    if offset > len(buffer) or len(ids) != count:
        raise TilePackError("a packed tile ends early")

    stars = []
    for i in range(count):
        star = {"id": ids[i]}
        for axis_index, axis_name in enumerate(("x", "y", "z")):
            low, step = axes[axis_index]
            star[axis_name] = low + positions[axis_index][i] * step
        star["luminosity_sol"] = _decode_log(luminosity[i], *LUMINOSITY_LOG_RANGE, _U16) or 0.0
        value = _decode_log(temperature[i], *TEMPERATURE_LOG_RANGE, _U8)
        if value is not None:
            star["temperature_k"] = value
        value = _decode_log(radius[i], *RADIUS_LOG_RANGE, _U8)
        if value is not None:
            star["radius_sol"] = value
        if classes[i]:
            star["star_type"] = chr(classes[i])
        star["charted"] = bool(charted[i])
        stars.append(star)
    return stars, offset


def unpack_tile(packed):
    """
    The tile dict a packed tile stands for (the page's script builds the same
    one): `stars` and `generated` as lists of dicts with `charted` (a bool)
    where the JSON form had `system_id`, and the JSON sections as stored.

    Raises:
        TilePackError: On bytes that are not a whole packed tile.
    """
    if packed[:4] != MAGIC:
        raise TilePackError("not a packed tile")
    try:
        (meta_length,) = struct.unpack_from("<I", packed, 4)
        meta = json.loads(packed[8:8 + meta_length].decode("utf-8"))
        offset = 8 + meta_length
        tile = dict(meta)
        for section in _STAR_SECTIONS:
            tile[section], offset = _unpack_stars(packed, offset)
    except (struct.error, ValueError) as exc:
        raise TilePackError(f"a packed tile is damaged: {exc}")
    if offset != len(packed):
        raise TilePackError("a packed tile has bytes left over")
    return tile


def pack_response(header, parts):
    """
    A `/galaxy/tiles` answer: `"PGTB"`, the header's JSON length (u32), the header (its
    own fields and `tiles`, a list of `[key, byte length]` in order) and
    then every tile's packed bytes one after the other.

    Args:
        header (dict): `stamp`, `generation`, ... as `tilecache` builds it.
        parts (list): `[(tile key, packed bytes), ...]`.
    """
    head = dict(header)
    head["tiles"] = [[key, len(raw)] for key, raw in parts]
    head_bytes = json.dumps(head, separators=(",", ":")).encode("utf-8")
    return b"".join([RESPONSE_MAGIC, struct.pack("<I", len(head_bytes)), head_bytes] + [raw for _, raw in parts])


def unpack_response(body):
    """
    `pack_response`'s inverse: `(header, {key: tile dict})` with the header's
    `tiles` list removed.
    """
    if body[:4] != RESPONSE_MAGIC:
        raise TilePackError("not a packed tile response")
    (head_length,) = struct.unpack_from("<I", body, 4)
    header = json.loads(body[8:8 + head_length].decode("utf-8"))
    offset = 8 + head_length
    tiles = {}
    for key, length in header.pop("tiles"):
        tiles[key] = unpack_tile(body[offset:offset + length])
        offset += length
    return header, tiles
