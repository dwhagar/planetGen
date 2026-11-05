"""
The packed form of a system's sky (VIEW.29, S11 of reports/system-sky-from-the-system-2026-11-04.md): what the sky route caches and
sends so that the browser can put every star at its light-time position on the graphics card, from the orbit clock.

The server works in doubles and sends offsets from the observer (`galaxy.sky.Sky`), so 32-bit floats on the card are enough: a star at 1 kpc keeps
a hundredth of an arcsecond. A star is 26 bytes, in columns (about 15 gzipped); 25,000 stars are 0.65 MB. static/skypack.js reads it and
`tests/test_skypack.py` with `tests/js/skypack.test.mjs` (fed by `test_js_unit.sky_pack`) hold the two to the same bytes.

One list (all numbers little-endian)::

    "PGK1"                        4 bytes, the format
    count                         u32
    epoch_unix                    f64: the time the offsets hold at (the sky's time)
    limit                         f32: the limiting magnitude the list was built to
    velocity_scale                f32: km/s of one velocity code
    observer_velocity             3 x f64, km/s: a star's own velocity is its relative one plus this
    meta_length                   u32
    meta                          UTF-8 JSON (the observer, the cache's version key, the tiers' report), or `{}`
    then, each a column of `count` values:
    offset x, y, z                3 x f32: parsecs, the star's place at `epoch_unix` as an offset from the observer, not seen back by its
                                  light time (the card does that: apparent = offset + relative velocity * t - velocity * tau)
    velocity x, y, z              3 x i16: km/s relative to the observer, code * velocity_scale
    magnitude                     i16: the star's absolute V magnitude plus the dust of its line of sight, in steps of 1/512 magnitude
                                  (the magnitude at distance d parsecs is that + 5 log10(d / 10))
    temperature                   u8: 0 is "none", else 1..255 over log10 TEMPERATURE_LOG_RANGE (the card colours it)
    kind                          u8: the tier, `KIND_CODES`
    index                         u32: the star's place in the list the caller built it from (to find a star again, to pick it)
"""

import hashlib
import json
import math
import struct

import numpy as np

MAGIC = b"PGK1"
MIMETYPE = "application/vnd.planetgen.sky"
"""str: What the sky list is sent as."""

TEMPERATURE_LOG_RANGE = (2.0, 6.0)
MAGNITUDE_STEPS_PER_MAG = 512
"""int: Codes of the magnitude column to one magnitude (range -64 to +64)."""

KIND_CODES = {"charted": 0, "bright_star": 1, "phenomenon": 2, "cluster": 3}
"""dict: The tier of each star, as the kind column numbers it."""

SKY_VERSION = 1
"""int: Raised when the sky an unchanged galaxy gives changes (a new field, a new constant), so every stored copy is built again (VIEW.31)."""

DEFAULT_LIMIT = 7.5
"""float: The limiting magnitude the sky route takes when `?limit=` is left out (`galaxy.sky.DEFAULT_LIMIT`; tests/test_system_sky_route.py holds the two equal)."""

MIN_LIMIT, MAX_LIMIT = 2.0, 9.0
"""float: The faintest and the brightest limiting magnitude the route takes (stand-in SKV1: the star count grows 2.5 to 3 times a magnitude, and the stored
stars and the sky tiers stop at 9; under 2 the sky holds a few dozen stars)."""


class SkyLimitError(ValueError):
    """A limiting magnitude that is not a number or is out of range."""


def check_limit(limit):
    """`limit` (text or a number, `None` or empty for `DEFAULT_LIMIT`) as a float from `MIN_LIMIT` to `MAX_LIMIT`, rounded to 0.01.

    Raises:
        SkyLimitError: Not a number, or out of range.
    """
    if limit is None or limit == "":
        return DEFAULT_LIMIT
    try:
        value = float(limit)
    except (TypeError, ValueError):
        raise SkyLimitError("limit must be a number") from None
    if not math.isfinite(value) or not MIN_LIMIT <= value <= MAX_LIMIT:
        raise SkyLimitError(f"limit must be from {MIN_LIMIT:g} to {MAX_LIMIT:g}")
    return round(value, 2)


def version_key(tile_stamp):
    """The key a stored sky is good for: 16 hex characters from the galaxy's tile stamp (`db.query.galaxy_tile_stamp`: sectors, systems, bright stars and
    the generator epoch) and `SKY_VERSION`."""
    return hashlib.sha256(f"{tile_stamp}.sky{SKY_VERSION}".encode("utf-8")).hexdigest()[:16]


BYTES_PER_STAR = 12 + 6 + 2 + 1 + 1 + 4
_HEADER_BYTES = 4 + 4 + 8 + 4 + 4 + 24 + 4
_I16 = 32767
_U8 = 255


class SkyPackError(ValueError):
    """Bytes that are not a packed sky, or a sky that cannot be packed."""


def _log_codes(values, lo, hi, top):
    values = np.asarray(values, dtype=float)
    codes = np.zeros(len(values), dtype=np.uint8)
    good = np.isfinite(values) & (values > 0)
    position = (np.log10(values[good]) - lo) / (hi - lo)
    codes[good] = 1 + np.clip(np.rint(position * (top - 1)), 0, top - 1).astype(np.uint8)
    return codes


def pack_sky(sky, epoch_unix, limit, observer_velocity_kms=(0.0, 0.0, 0.0), meta=None, indices=None):
    """
    The packed list of a `galaxy.sky.Sky`.

    Args:
        sky (Sky): The sky (its `offset_pc`, `velocity_kms`, `dimmed_absolute_magnitude`, `temperature_k` and `tier`).
        epoch_unix (float): The sky's time, which the offsets hold at.
        limit (float): The limiting magnitude it was built to.
        observer_velocity_kms (sequence): The observer's velocity, km/s.
        meta (dict): JSON-able extras to ride along (default none).
        indices (sequence): Each star's index (default its place in the sky).

    Returns:
        bytes: The packed list.

    Raises:
        SkyPackError: A tier with no code, a magnitude out of range, or an index past 32 bits.
    """
    count = len(sky.uid)
    observer = np.asarray(observer_velocity_kms, dtype=float).reshape(3)
    offsets = np.asarray(sky.offset_pc, dtype=float).reshape(count, 3)
    relative = np.asarray(sky.velocity_kms, dtype=float).reshape(count, 3) - observer[None, :]
    peak = float(np.max(np.abs(relative))) if count else 0.0
    scale = np.float32(peak / _I16 if peak > 0.0 else 1.0e-3)
    velocity_codes = np.rint(relative / float(scale)).astype(np.int16)
    magnitude = np.asarray(sky.dimmed_absolute_magnitude, dtype=float)
    steps = np.rint(magnitude * MAGNITUDE_STEPS_PER_MAG)
    if count and (not np.all(np.isfinite(steps)) or np.abs(steps).max() > _I16):
        raise SkyPackError("an absolute magnitude is outside -64 to +64")
    try:
        kinds = np.array([KIND_CODES[str(t)] for t in sky.tier], dtype=np.uint8)
    except KeyError as exc:
        raise SkyPackError(f"a star of tier {exc.args[0]!r} has no kind code") from exc
    index = np.arange(count, dtype=np.int64) if indices is None else np.asarray(indices, dtype=np.int64)
    if count and (index.min() < 0 or index.max() >= 1 << 32):
        raise SkyPackError("an index is outside 0 to 2^32 - 1")
    meta_bytes = json.dumps(meta or {}, separators=(",", ":")).encode("utf-8")
    low, high = TEMPERATURE_LOG_RANGE
    out = [MAGIC, struct.pack("<IdffdddI", count, float(epoch_unix), float(limit), float(scale), *observer, len(meta_bytes)), meta_bytes]
    for axis in range(3):
        out.append(offsets[:, axis].astype("<f4").tobytes())
    for axis in range(3):
        out.append(velocity_codes[:, axis].astype("<i2").tobytes())
    out.append(steps.astype("<i2").tobytes())
    out.append(_log_codes(sky.temperature_k, low, high, _U8).tobytes())
    out.append(kinds.tobytes())
    out.append(index.astype("<u4").tobytes())
    return b"".join(out)


def is_whole(raw):
    """Whether `raw` is exactly one packed list, found from its lengths alone."""
    if raw is None or raw[:4] != MAGIC or len(raw) < _HEADER_BYTES:
        return False
    count, meta_length = struct.unpack_from("<I", raw, 4)[0], struct.unpack_from("<I", raw, _HEADER_BYTES - 4)[0]
    return len(raw) == _HEADER_BYTES + meta_length + count * BYTES_PER_STAR


def unpack_sky(raw):
    """
    A packed list as a dict: `count`, `epoch_unix`, `limit`, `velocity_scale`, `observer_velocity_kms`, `meta` and the columns as numpy arrays
    decoded the way static/skypack.js decodes them: `offset_pc` and `velocity_kms` `(count, 3)` (the velocity relative to the observer),
    `magnitude` (the dimmed absolute magnitude), `temperature_k` (`nan` for none), `kind`, `index`.

    Raises:
        SkyPackError: Bytes that are not one whole packed list.
    """
    raw = bytes(raw)
    if not is_whole(raw):
        raise SkyPackError("not a packed sky")
    count, epoch, limit, scale, ox, oy, oz, meta_length = struct.unpack_from("<IdffdddI", raw, 4)
    at = _HEADER_BYTES
    meta = json.loads(raw[at:at + meta_length].decode("utf-8"))
    at += meta_length

    def column(dtype):
        nonlocal at
        values = np.frombuffer(raw, dtype=dtype, count=count, offset=at)
        at += values.nbytes
        return values

    offsets = np.stack([column("<f4") for _ in range(3)], axis=1).astype(float)
    velocity = np.stack([column("<i2") for _ in range(3)], axis=1).astype(float) * scale
    magnitude = column("<i2").astype(float) / MAGNITUDE_STEPS_PER_MAG
    code = column("u1")
    low, high = TEMPERATURE_LOG_RANGE
    temperature = np.where(code == 0, np.nan, 10.0 ** (low + (code.astype(float) - 1.0) / (_U8 - 1) * (high - low)))
    return {"count": count, "epoch_unix": epoch, "limit": limit, "velocity_scale": scale, "observer_velocity_kms": (ox, oy, oz), "meta": meta,
            "offset_pc": offsets, "velocity_kms": velocity, "magnitude": magnitude, "temperature_k": temperature, "kind": column("u1").copy(),
            "index": column("<u4").astype(np.int64)}
