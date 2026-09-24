# html/lib/tilecache.py

"""
On-disk cache of the 3D Galaxy Map's cube tiles, kept by the web layer so
the map doesn't hit the API (and through it, the database) every time the
camera moves.

A tile's contents (see `stellarObjects.galaxyViewport`'s "Cube tiles"
section and `queryDb.galaxy_tiles`) depend only on its key and on the
database's contents, which `queryDb.galaxy_content_stamp` summarizes as a
short "stamp". Tiles are cached as `<cache dir>/<db>/<generation>/
<tile>.json`, where the generation is the stamp at the last time every
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

The browser keeps its own copy of every tile too (`static/galaxymap3d.js`,
in `localStorage`, keyed by the generation). `stamp.json` remembers the
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

from apiclient import get_galaxy_changes, get_galaxy_tiles
from stellarObjects.appconfig import load_config
from stellarObjects.galaxyViewport import parse_tile_key, tile_key

DEFAULT_CACHE_DIR = "/var/cache/planetgen/tiles"
"""str: Used when neither `PLANETGEN_TILE_CACHE_DIR` nor `config.json`'s
`tile_cache.dir` names a directory. Falls back to a `planetgen-tiles`
folder in the system temp directory when Apache can't create this one."""

STAMP_TTL_SECONDS = 60
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


class TileRequestError(ValueError):
    """A malformed tile request from the map's own client-side JS (bad key,
    too many keys) -- reported as a 400, not an API failure."""


def _config():
    return load_config().get("tile_cache") or {}


def max_cache_bytes():
    """The cache's size budget, bytes. `0` disables the disk cache."""
    raw = os.environ.get("PLANETGEN_TILE_CACHE_MAX_MB")
    if raw is None:
        raw = _config().get("max_mb", 200)
    try:
        return max(0, int(float(raw) * 1024 * 1024))
    except (TypeError, ValueError):
        return 200 * 1024 * 1024


def configured_cache_dir():
    """
    The directory the cache is meant to live in (`PLANETGEN_TILE_CACHE_DIR`,
    then `tile_cache.dir`, then `DEFAULT_CACHE_DIR`), or `None` when the
    disk cache is off. Doesn't create or check anything --
    `examples/apache/create-cache-dir.sh` asks this where to create it.
    """
    if max_cache_bytes() == 0:
        return None
    return os.environ.get("PLANETGEN_TILE_CACHE_DIR") or _config().get("dir") or DEFAULT_CACHE_DIR


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
        candidates.append(os.path.join(tempfile.gettempdir(), "planetgen-tiles"))
    for candidate in candidates:
        try:
            os.makedirs(candidate, mode=0o750, exist_ok=True)
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
            os.unlink(os.path.join(generation_dir, _tile_filename("t", key)))
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

    changes = get_galaxy_changes(db, remembered["state"] if remembered else None)
    stamp = changes.get("stamp")
    if not _STAMP_RE.match(str(stamp)):
        return {"stamp": stamp, "generation": None, "history": []}

    stale = []
    if remembered is not None and not changes.get("full"):
        generation = remembered["generation"]
        history = remembered["history"]
        if stamp != remembered["stamp"]:
            stale = [str(key) for key in changes.get("tiles") or []]
            history = _trim_history(history + [{"from": remembered["stamp"], "to": stamp, "tiles": stale}])
    else:
        generation = stamp
        history = []

    info = {"stamp": stamp, "state": str(changes.get("state") or ""), "generation": generation, "history": history}
    generation_dir = os.path.join(db_dir, generation)
    # Deleted both before and after the new stamp is written: a request
    # still working under the old stamp checks it before writing a tile
    # (see `fetch_tiles`), so whichever side of the write it lands on, a
    # tile fetched before the change doesn't survive it.
    _delete_tiles(generation_dir, stale)
    _write_json(stamp_path, info)
    _delete_tiles(generation_dir, stale)
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


def _tile_filename(prefix, key):
    level, ix, iy, iz = parse_tile_key(key)
    return f"{prefix}{level}_{ix}_{iy}_{iz}.json"


def validate_request(tile_keys, density_key):
    """
    Canonicalizes and validates a tile request.

    Returns:
        tuple: `(tile_keys, density_key)` with duplicates removed and keys
            in canonical form.

    Raises:
        TileRequestError: On a malformed key or too many keys.
    """
    try:
        keys = list(dict.fromkeys(tile_key(*parse_tile_key(key)) for key in tile_keys))
        density = tile_key(*parse_tile_key(density_key)) if density_key else None
    except ValueError as exc:
        raise TileRequestError(str(exc))
    if len(keys) > MAX_TILES_PER_REQUEST:
        raise TileRequestError(f"at most {MAX_TILES_PER_REQUEST} tiles per request")
    return keys, density


def fetch_tiles(db, tile_keys, density_key=None, known_stamp=None):
    """
    The requested tiles (and optional density cloud), from the disk cache
    where possible and from `GET /api/galaxy/tiles` for the rest, caching
    whatever the API returns.

    Args:
        db (str): The `?db=` value.
        tile_keys (list[str]): `level/ix/iy/iz` keys.
        density_key (str or None): A tile key to anchor a density cloud on.
        known_stamp (str or None): The stamp the browser's own cache is
            at, when it has one.

    Returns:
        dict: `stamp` and `generation` (see `current_stamp` -- the browser
            keys its own cache by the generation), `history` (only when
            `known_stamp` isn't current: the changed tiles the browser
            hasn't seen, see `_stale_since`), `tiles` (`{key: {"placed",
            "planned"}}`), `density` (`{"key", "points"}` or `None`),
            `edge_pc`, `has_shape`, and `cached` (how many of the
            requested parts came from disk -- for diagnostics).

    Raises:
        TileRequestError: On a malformed request.
        apiclient.ApiError / NotFoundError: If the API is needed and fails.
    """
    tile_keys, density_key = validate_request(tile_keys, density_key)
    root = cache_dir()
    info = current_stamp(db, root)
    stamp = info["stamp"]
    generation = info["generation"]
    generation_dir = os.path.join(_db_dir(root, db), generation) if root and generation else None

    tiles = {}
    density = None
    meta = None
    if generation_dir:
        meta = _read_json(os.path.join(generation_dir, "meta.json"))
        for key in tile_keys:
            cached = _read_json(os.path.join(generation_dir, _tile_filename("t", key)))
            if isinstance(cached, dict):
                tiles[key] = cached
        if density_key:
            cached = _read_json(os.path.join(generation_dir, _tile_filename("d", density_key)))
            if isinstance(cached, list):
                density = {"key": density_key, "points": cached}

    cached_count = len(tiles) + (1 if density is not None else 0)
    missing = [key for key in tile_keys if key not in tiles]
    need_density = bool(density_key) and density is None

    if missing or need_density or not isinstance(meta, dict):
        fetched = get_galaxy_tiles(db, missing, density_key if need_density else None)
        meta = {"edge_pc": fetched.get("edge_pc"), "has_shape": bool(fetched.get("has_shape"))}
        fetched_tiles = {key: _round_floats(value) for key, value in (fetched.get("tiles") or {}).items()}
        fetched_density = fetched.get("density")
        if fetched_density is not None:
            fetched_density = {"key": density_key, "points": _round_floats(fetched_density.get("points") or [])}

        # Another request may have found changes while this one was
        # fetching; what it fetched could be from before them.
        if generation_dir and (_read_remembered(_db_dir(root, db)) or {}).get("stamp") == stamp:
            _write_json(os.path.join(generation_dir, "meta.json"), meta)
            for key, value in fetched_tiles.items():
                if key in missing:
                    _write_json(os.path.join(generation_dir, _tile_filename("t", key)), value)
            if fetched_density is not None:
                _write_json(os.path.join(generation_dir, _tile_filename("d", density_key)), fetched_density["points"])
            if random.random() < PRUNE_PROBABILITY:
                prune(root)

        tiles.update(fetched_tiles)
        if fetched_density is not None:
            density = fetched_density

    result = {
        "stamp": stamp,
        "generation": generation,
        "tiles": {key: tiles[key] for key in tile_keys if key in tiles},
        "density": density,
        "edge_pc": meta.get("edge_pc"),
        "has_shape": bool(meta.get("has_shape")),
        "cached": cached_count,
    }
    history = _stale_since(known_stamp, info)
    if history is not None:
        result["history"] = history
    return result


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
