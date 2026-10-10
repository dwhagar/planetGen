# planetgen/web/lib/tilecache.py

"""
On-disk cache of the 3D Galaxy Map's cube tiles, kept by the web layer so
the map doesn't hit the API (and through it, the database) every time the
camera moves.

A tile's contents (see `planetgen.galaxy.viewport`'s "Cube tiles"
section and `queryDb.galaxy_tiles`) depend only on its key and on the
database's contents, which `queryDb.galaxy_content_stamp` summarizes as a
short "stamp". Tiles are cached as `<cache dir>/<db>/<generation>/
<tile>.json` (the finished, trimmed JSON bytes, see `trim_tile`), where the generation is the stamp at the last time every
tile went stale. `fetch_tiles` asks the API only for the tiles it doesn't
already have, and when every tile is cached it doesn't call the API for
tiles at all.

Freshness is checked at most once per `STAMP_TTL_SECONDS` per database,
with one cheap API call (`GET /api/galaxy/changes`, see
`queryDb.galaxy_changes`), which reads the schema-v27 `modified_at`
columns and answers which tiles changed since the last check. Only those
tile files are deleted, so renaming a sector refetches the dozen tiles
holding it rather than the whole map. A change that can't be pinned to
tiles (a deleted sector, a re-planned galaxy, a new planetGen release)
answers `full`, and a new generation starts with nothing cached.

A request is answered by joining those stored bytes, never by parsing and
serialising the tiles again (`fetch_tiles_wire`); Apache compresses the
response (`mod_brotli`, else `mod_deflate`, see `examples/apache`). Sending
the tiles' own gzip copies as the members of one gzip stream was tried and
dropped: Chromium decodes only the first member.

The browser keeps its own copy of every tile too (`static/galaxymap3d.js`,
in IndexedDB, keyed by the generation, MAP.158). `stamp.json` remembers the
last few checks' changed tiles (`history`), and `fetch_tiles` hands them
to the browser whenever its stamp is out of date, so the browser drops
only those tiles too.

Everything here fails open: an unwritable or missing cache directory, a
corrupt file, or a race with another request just means the tile is
fetched from the API again, never an error page.
"""

import hashlib
import json
import os
import random
import re
import tempfile
import time

from planetgen.web.lib import apiclient
from planetgen.web.lib.apiclient import get_galaxy_changes, get_galaxy_stage, get_galaxy_tiles
from planetgen.web.lib.privatedir import ensure_private_dir
from planetgen.util.settings import get_settings
from planetgen.galaxy.drill import format_drill_key, parse_drill_key
from planetgen.galaxy.viewport import parse_tile_key, tile_key

DEFAULT_CACHE_DIR = "/var/cache/planetgen/tiles"
"""str: Used when neither `PLANETGEN_TILE_CACHE_DIR` nor `config.json`'s
`tile_cache.dir` names a directory. Falls back to a `planetgen-tiles`
folder in the system temp directory when Apache can't create this one."""

STAMP_TTL_SECONDS = 60

BUSY_KEEP_SECONDS = 600
"""int: While the galaxy changes faster than tiles can be refetched (a fill
adds sectors and stars every second, PERF.34), a cached generation is kept
and served at most this long past its start, instead of being thrown away at
every check, which made every map request recompute its tiles on a database
already under load. Then one full refresh brings it up to date."""

FAILED_CHECK_RETRY_SECONDS = 15
"""int: After the freshness check itself fails (a database too busy to answer
in time), the cache is served as it is and the check is not tried again for
this long."""

_failed_checks = {}
"""dict: `db -> time.monotonic()` of that database's last failed check."""
"""int: How long a database's stamp is trusted before asking the API again
-- the longest a newly generated or edited sector can take to appear on
an open map."""

HISTORY_MAX_KEYS = 1500
"""int: Most changed-tile keys `stamp.json`'s `history` keeps, across all
its entries. Oldest entries go first; a browser whose stamp is older than
what's left just drops all its tiles."""

HISTORY_MAX_ENTRIES = 32
"""int: Most checks `history` remembers."""

PRUNE_PROBABILITY = 0.05
"""float: Chance that a request which wrote new tiles also checks the
cache's total size (a directory walk, so not on every request)."""

MAX_TILES_PER_REQUEST = 128
"""int: Mirrors `queryDb.MAX_TILES_PER_REQUEST` -- rejected here before
anything is read or fetched."""

_STAMP_RE = re.compile(r"^[0-9a-f]{16}$")
_ROUND_DIGITS = 4

_TILE_STAR_KEYS = ("id", "x", "y", "z", "luminosity_sol", "temperature_k", "radius_sol", "star_type", "system_id")
"""tuple: The only fields of a tile's star the page's script reads (MAP.157;
`tileStars`, `setStars` and `starShown` in `static/galaxymap3d.js`)."""

_TILE_KEEP_SECTIONS = ("clouds", "stars", "generated", "points")
"""tuple: The sections of a tile the page reads; `placed`, `planned` and
`filled` (9 to 35% of every tile) are not served."""



class TileRequestError(ValueError):
    """A malformed tile request from the map's own client-side JS (bad key,
    too many keys) -- reported as a 400, not an API failure."""


def _config():
    return get_settings().tile_cache


def max_cache_bytes():
    """The cache's size budget, bytes. `0` disables the disk cache."""
    return max(0, int(_config().max_mb * 1024 * 1024))


def configured_cache_dir():
    """
    The directory the cache is meant to live in (`PLANETGEN_TILE_CACHE_DIR`,
    then `tile_cache.dir`, then `DEFAULT_CACHE_DIR`), or `None` when the
    disk cache is off. Doesn't create or check anything --
    `examples/apache/create-cache-dir.sh` asks this where to create it.
    """
    if max_cache_bytes() == 0:
        return None
    return _config().dir or DEFAULT_CACHE_DIR


def cache_dir():
    """
    The writable cache root, created if needed, or `None` when the disk
    cache is off or no candidate directory is writable.
    """
    configured = configured_cache_dir()
    if configured is None:
        return None
    candidates = [configured]
    if configured == DEFAULT_CACHE_DIR:
        # Private: one another local user created or can write to is
        # refused (they could plant tiles, or read the cache).
        candidates.append(os.path.join(tempfile.gettempdir(), "planetgen-tiles"))
    for index, candidate in enumerate(candidates):
        try:
            if index == 0:
                os.makedirs(candidate, mode=0o750, exist_ok=True)
            else:
                ensure_private_dir(candidate)
        except OSError:
            continue
        if os.access(candidate, os.W_OK | os.X_OK):
            return candidate
    return None


def _db_dir(root, db):
    # Hashed rather than the raw name, so no database name can ever form
    # an unexpected path.
    return os.path.join(root, hashlib.sha256(db.encode("utf-8")).hexdigest()[:16])


def _read_json(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def _write_json(path, payload):
    """Atomic write (temp file + rename), so a concurrent reader never sees
    half a file. Returns whether it worked."""
    directory = os.path.dirname(path)
    tmp_path = None
    try:
        os.makedirs(directory, mode=0o750, exist_ok=True)
        fd, tmp_path = tempfile.mkstemp(dir=directory, prefix=".tmp-", suffix=".json")
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(payload, f, separators=(",", ":"))
        os.replace(tmp_path, path)
        return True
    except OSError:
        if tmp_path:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
        return False


def _remove_tree(path):
    for dirpath, dirnames, filenames in os.walk(path, topdown=False):
        for name in filenames:
            try:
                os.unlink(os.path.join(dirpath, name))
            except OSError:
                pass
        for name in dirnames:
            try:
                os.rmdir(os.path.join(dirpath, name))
            except OSError:
                pass
    try:
        os.rmdir(path)
    except OSError:
        pass


def _read_remembered(db_dir):
    """`stamp.json`'s contents when well-formed, else `None` (including
    the pre-generation `{"stamp"}` format, which just starts afresh)."""
    remembered = _read_json(os.path.join(db_dir, "stamp.json"))
    if not (
        isinstance(remembered, dict)
        and _STAMP_RE.match(str(remembered.get("stamp", "")))
        and _STAMP_RE.match(str(remembered.get("generation", "")))
        and isinstance(remembered.get("state"), str)
        and isinstance(remembered.get("history"), list)
    ):
        return None
    return remembered


def _delete_tiles(generation_dir, tile_keys):
    for key in tile_keys:
        try:
            path = os.path.join(generation_dir, _tile_filename("t", key))
        except ValueError:
            continue
        try:
            os.unlink(path)
        except OSError:
            pass


def _stage_filename(key):
    """A drill-down stage's cache file: `sgalaxy.json` for the galaxy,
    else `s<m>_<ring>_<wedge>_<slab>.json`."""
    if key == "galaxy":
        return "sgalaxy.json"
    return "s" + format_drill_key(parse_drill_key(key)).replace(".", "_") + ".json"


def _delete_stages(generation_dir, stage_keys):
    for key in stage_keys:
        try:
            os.unlink(os.path.join(generation_dir, _stage_filename(key)))
        except (OSError, ValueError):
            pass


def _trim_history(history):
    """The newest `history` entries that fit `HISTORY_MAX_ENTRIES` and
    `HISTORY_MAX_KEYS`."""
    kept = []
    total = 0
    for entry in reversed(history[-HISTORY_MAX_ENTRIES:]):
        total += len(entry["tiles"])
        if total > HISTORY_MAX_KEYS:
            break
        kept.append(entry)
    return list(reversed(kept))


def current_stamp(db, root=None):
    """
    The database's current stamp, generation and change history: from
    `<db dir>/stamp.json` when that's younger than `STAMP_TTL_SECONDS`,
    else from `GET /api/galaxy/changes` (and then remembered). A check
    that finds changes deletes just the changed tiles, or on a `full`
    answer starts a new generation and deletes every other one.

    Returns:
        dict: `stamp`, `generation` (`None` when there's no disk cache),
            and `history` (a list of `{"from", "to", "tiles"}`, oldest
            first: each check that found changes, and the tile keys it
            found).
    """
    if root is None:
        return {"stamp": get_galaxy_changes(db)["stamp"], "generation": None, "history": []}

    db_dir = _db_dir(root, db)
    stamp_path = os.path.join(db_dir, "stamp.json")
    remembered = _read_remembered(db_dir)
    try:
        age = time.time() - os.path.getmtime(stamp_path)
    except OSError:
        age = None
    if remembered is not None and age is not None and 0 <= age < STAMP_TTL_SECONDS:
        return remembered

    if remembered is not None and time.monotonic() - _failed_checks.get(db, -1e9) < FAILED_CHECK_RETRY_SECONDS:
        return remembered
    try:
        changes = get_galaxy_changes(db, remembered["state"] if remembered else None)
    except apiclient.ApiError:
        if remembered is None:
            raise
        # A busy database that can't answer the check is no reason to fail a
        # request the cache can serve (PERF.34).
        _failed_checks[db] = time.monotonic()
        return remembered
    _failed_checks.pop(db, None)
    stamp = changes.get("stamp")
    if changes.get("busy") and remembered is not None:
        generation_age = time.time() - float(remembered.get("started") or 0)
        if 0 <= generation_age < BUSY_KEEP_SECONDS:
            # Keep the stored state, so the next check still sees everything
            # that changed since; just look again in a minute.
            _write_json(stamp_path, remembered)
            return remembered
    if not _STAMP_RE.match(str(stamp)):
        return {"stamp": stamp, "generation": None, "history": []}

    stale = []
    stale_stages = []
    started = time.time()
    if remembered is not None and not changes.get("full"):
        generation = remembered["generation"]
        started = remembered.get("started") or started
        history = remembered["history"]
        if stamp != remembered["stamp"]:
            stale = [str(key) for key in changes.get("tiles") or []]
            stale_stages = [str(key) for key in changes.get("stages") or []]
            history = _trim_history(history + [{"from": remembered["stamp"], "to": stamp, "tiles": stale}])
    else:
        generation = stamp
        history = []

    info = {"stamp": stamp, "state": str(changes.get("state") or ""), "generation": generation, "history": history,
            "started": started}
    generation_dir = os.path.join(db_dir, generation)
    # Deleted both before and after the new stamp is written: a request
    # still working under the old stamp checks it before writing a tile
    # (see `fetch_tiles`), so whichever side of the write it lands on, a
    # tile fetched before the change doesn't survive it.
    _delete_tiles(generation_dir, stale)
    _delete_stages(generation_dir, stale_stages)
    _write_json(stamp_path, info)
    _delete_tiles(generation_dir, stale)
    _delete_stages(generation_dir, stale_stages)
    try:
        for entry in os.scandir(db_dir):
            if entry.is_dir() and entry.name != generation:
                _remove_tree(entry.path)
    except OSError:
        pass
    return info


def _stale_since(known_stamp, info):
    """
    The `history` entries a browser holding `known_stamp` hasn't seen yet,
    so it drops just their tiles: `None` when its stamp is current
    (nothing to send), all of `history` when its stamp isn't known (the
    first page render, before the browser's own storage is read), and
    `[]` when its stamp is older than the history, which the browser
    reads as "drop everything".
    """
    if known_stamp == info["stamp"]:
        return None
    history = info.get("history") or []
    if known_stamp is None:
        return history
    for index in range(len(history) - 1, -1, -1):
        if history[index]["from"] == known_stamp:
            return history[index:]
    return []


def _round_floats(value):
    """Rounds every float in a JSON-able structure -- tiles are mostly
    coordinates, and full double precision roughly doubles their size for
    no visible difference."""
    if isinstance(value, float):
        return round(value, _ROUND_DIGITS)
    if isinstance(value, dict):
        return {k: _round_floats(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_round_floats(v) for v in value]
    return value


def _significant(value, digits):
    """`value` rounded to `digits` significant digits (floats only)."""
    if not isinstance(value, float) or value == 0.0:
        return value
    return float(f"{value:.{digits - 1}e}")


def _trim_star(star):
    """One star of a tile as the page reads it (MAP.157): the fields in
    `_TILE_STAR_KEYS`, `star_type` as its class letter, positions to 3
    decimals, luminosity to 4 and radius to 3 significant digits and
    temperature to 10 K. No `system_id` when there is none."""
    trimmed = {
        "id": star["id"],
        "x": round(star["x"], 3), "y": round(star["y"], 3), "z": round(star["z"], 3),
        "luminosity_sol": _significant(float(star["luminosity_sol"]), 4),
    }
    temperature = star.get("temperature_k")
    if temperature is not None:
        trimmed["temperature_k"] = int(round(temperature / 10.0)) * 10
    radius = star.get("radius_sol")
    if radius is not None:
        trimmed["radius_sol"] = _significant(float(radius), 3)
    star_type = star.get("star_type")
    if star_type:
        trimmed["star_type"] = str(star_type)[0]
    if star.get("system_id") is not None:
        trimmed["system_id"] = star["system_id"]
    return trimmed


def trim_tile(tile):
    """
    A tile as the page's script reads it (MAP.157): without `placed`,
    `planned` and `filled`, its stars cut to `_TILE_STAR_KEYS` and rounded
    (`_trim_star`), every other float rounded to `_ROUND_DIGITS` decimals.
    Measured on 400 real tiles this is 48% fewer gzipped bytes.
    """
    trimmed = {}
    for section in _TILE_KEEP_SECTIONS:
        value = tile.get(section)
        if value is None:
            continue
        if section in ("stars", "generated"):
            trimmed[section] = [_trim_star(star) for star in value]
        else:
            trimmed[section] = _round_floats(value)
    return trimmed


def _tile_bytes(tile):
    """A tile's finished JSON bytes."""
    return json.dumps(trim_tile(tile), separators=(",", ":")).encode("utf-8")


def _tile_filename(prefix, key):
    level, ix, iy, iz = parse_tile_key(key)
    return f"{prefix}{level}_{ix}_{iy}_{iz}.json"


def _read_bytes(path):
    try:
        with open(path, "rb") as f:
            return f.read()
    except OSError:
        return None


def _write_bytes(path, data):
    """Atomic write of raw bytes, like `_write_json`. Returns whether it
    worked."""
    directory = os.path.dirname(path)
    tmp_path = None
    try:
        os.makedirs(directory, mode=0o750, exist_ok=True)
        fd, tmp_path = tempfile.mkstemp(dir=directory, prefix=".tmp-", suffix=".part")
        with os.fdopen(fd, "wb") as f:
            f.write(data)
        os.replace(tmp_path, path)
        return True
    except OSError:
        if tmp_path:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
        return False


def validate_request(tile_keys):
    """
    Canonicalizes and validates a tile request.

    Returns:
        list[str]: The tile keys with duplicates removed, in canonical
            form.

    Raises:
        TileRequestError: On a malformed key or too many keys.
    """
    try:
        keys = list(dict.fromkeys(tile_key(*parse_tile_key(key)) for key in tile_keys))
    except ValueError as exc:
        raise TileRequestError(str(exc))
    if len(keys) > MAX_TILES_PER_REQUEST:
        raise TileRequestError(f"at most {MAX_TILES_PER_REQUEST} tiles per request")
    return keys


def _collect_tiles(db, tile_keys, known_stamp):
    """
    What `fetch_tiles` and `fetch_tiles_wire` share: the requested tiles'
    stored bytes, from the disk cache where possible and from
    `GET /api/galaxy/tiles` (trimmed, then cached) for the rest.

    Returns:
        tuple: `(header, parts)`. `header` is the response without its
            tiles (`stamp`, `generation`, `edge_pc`, `has_shape`, `cached`
            and, only when `known_stamp` isn't current, `history`);
            `parts` is `[(key, json bytes), ...]` in request order.
    """
    tile_keys = validate_request(tile_keys)
    root = cache_dir()
    info = current_stamp(db, root)
    stamp = info["stamp"]
    generation = info["generation"]
    generation_dir = os.path.join(_db_dir(root, db), generation) if root and generation else None

    stored = {}
    meta = None
    if generation_dir:
        meta = _read_json(os.path.join(generation_dir, "meta.json"))
        for key in tile_keys:
            path = os.path.join(generation_dir, _tile_filename("t", key))
            raw = _read_bytes(path)
            # A torn or foreign file is refetched, never sent.
            if raw is not None and raw[:1] == b"{" and raw[-1:] == b"}":
                stored[key] = raw

    cached_count = len(stored)
    missing = [key for key in tile_keys if key not in stored]

    if missing or not isinstance(meta, dict):
        fetched = get_galaxy_tiles(db, missing)
        meta = {"edge_pc": fetched.get("edge_pc"), "has_shape": bool(fetched.get("has_shape"))}
        made = {key: _tile_bytes(value) for key, value in (fetched.get("tiles") or {}).items()}

        # Another request may have found changes while this one was
        # fetching; what it fetched could be from before them.
        if generation_dir and (_read_remembered(_db_dir(root, db)) or {}).get("stamp") == stamp:
            _write_json(os.path.join(generation_dir, "meta.json"), meta)
            for key, raw in made.items():
                if key in missing:
                    _write_bytes(os.path.join(generation_dir, _tile_filename("t", key)), raw)
            if random.random() < PRUNE_PROBABILITY:
                prune(root)

        stored.update(made)

    header = {
        "stamp": stamp,
        "generation": generation,
        "edge_pc": meta.get("edge_pc"),
        "has_shape": bool(meta.get("has_shape")),
        "cached": cached_count,
    }
    history = _stale_since(known_stamp, info)
    if history is not None:
        header["history"] = history
    return header, [(key, stored[key]) for key in tile_keys if key in stored]


def fetch_tiles(db, tile_keys, known_stamp=None):
    """
    The requested tiles, from the disk cache
    where possible and from `GET /api/galaxy/tiles` for the rest, caching
    whatever the API returns (trimmed, see `trim_tile`).

    Args:
        db (str): The `?db=` value.
        tile_keys (list[str]): `level/ix/iy/iz` keys.
        known_stamp (str or None): The stamp the browser's own cache is
            at, when it has one.

    Returns:
        dict: `stamp` and `generation` (see `current_stamp` -- the browser
            keys its own cache by the generation), `history` (only when
            `known_stamp` isn't current: the changed tiles the browser
            hasn't seen, see `_stale_since`), `tiles` (`{key: {"clouds",
            "stars", "generated", "points"}}`), `edge_pc`, `has_shape`, and `cached` (how many of the
            requested parts came from disk -- for diagnostics).

    Raises:
        TileRequestError: On a malformed request.
        apiclient.ApiError / NotFoundError: If the API is needed and fails.
    """
    header, parts = _collect_tiles(db, tile_keys, known_stamp)
    return {**header, "tiles": {key: json.loads(raw) for key, raw in parts}}


def fetch_tiles_wire(db, tile_keys, known_stamp=None):
    """
    `fetch_tiles`' payload as the bytes to send (`/galaxy/tiles`), joined
    from the stored tile bytes without parsing them.
    """
    header, parts = _collect_tiles(db, tile_keys, known_stamp)
    rest = json.dumps(header, separators=(",", ":")).encode("utf-8")
    # `{"tiles":{"k":<tile>,"k2":<tile2>},"stamp":...}`: the header's object,
    # opened with "tiles" first.
    pieces = [b'{"tiles":{']
    for index, (key, raw) in enumerate(parts):
        pieces.append((b"," if index else b"") + json.dumps(key).encode("utf-8") + b":")
        pieces.append(raw)
    pieces.append(b"}," + rest[1:])
    return b"".join(pieces)


def fetch_stage(db, at=None):
    """
    One Galaxy Map drill-down stage (`queryDb.galaxy_stage`), from the disk
    cache when it's there, else from `GET /api/galaxy/stage` (and then
    cached). A sector change deletes just its chain's stages (`stages` in
    `/api/galaxy/changes`), the same way tiles are refreshed.

    Args:
        db (str): The `?db=` value.
        at (str or None): A block key `m.ring.wedge.slab`, or `None` for
            the galaxy.

    Returns:
        dict: The stage, plus `stamp`.

    Raises:
        TileRequestError: On a malformed key.
        apiclient.ApiError / NotFoundError: If the API is needed and fails.
    """
    try:
        key = format_drill_key(parse_drill_key(at)) if at else "galaxy"
    except ValueError as exc:
        raise TileRequestError(str(exc))
    root = cache_dir()
    info = current_stamp(db, root)
    stamp = info["stamp"]
    generation = info["generation"]
    generation_dir = os.path.join(_db_dir(root, db), generation) if root and generation else None
    path = os.path.join(generation_dir, _stage_filename(key)) if generation_dir else None

    stage = _read_json(path) if path else None
    if not isinstance(stage, dict):
        stage = get_galaxy_stage(db, None if key == "galaxy" else key)
        if path and (_read_remembered(_db_dir(root, db)) or {}).get("stamp") == stamp:
            _write_json(path, stage)
    return {**stage, "stamp": stamp}


def prune(root, max_bytes=None):
    """
    Deletes the least recently written cache files until the cache is
    under 80% of `max_bytes` (default `max_cache_bytes()`), then removes
    any directories that left empty. Stamp files are skipped.
    """
    if max_bytes is None:
        max_bytes = max_cache_bytes()
    files = []
    total = 0
    for dirpath, _dirnames, filenames in os.walk(root):
        for name in filenames:
            if name == "stamp.json":
                continue
            path = os.path.join(dirpath, name)
            try:
                info = os.stat(path)
            except OSError:
                continue
            files.append((info.st_mtime, info.st_size, path))
            total += info.st_size
    if total <= max_bytes:
        return
    target = max_bytes * 0.8
    for _mtime, size, path in sorted(files):
        if total <= target:
            break
        try:
            os.unlink(path)
            total -= size
        except OSError:
            pass
    for dirpath, dirnames, filenames in os.walk(root, topdown=False):
        if dirpath != root and not dirnames and not filenames:
            try:
                os.rmdir(dirpath)
            except OSError:
                pass
