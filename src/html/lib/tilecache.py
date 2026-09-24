# html/lib/tilecache.py

"""
On-disk cache of the 3D Galaxy Map's cube tiles, kept by the web layer so
the map doesn't hit the API (and through it, the database) every time the
camera moves.

A tile's contents (see `stellarObjects.galaxyViewport`'s "Cube tiles"
section and `queryDb.galaxy_tiles`) depend only on its key and on the
database's contents, which `queryDb.galaxy_content_stamp` summarizes as a
short "stamp". So every tile is cached under
`<cache dir>/<db>/<stamp>/<tile>.json`, and a cached tile is reused until
the stamp changes (new sectors generated, the galaxy re-planned, or a new
planetGen release). `fetch_tiles` asks the API only for the tiles it
doesn't already have, and when every tile is cached it doesn't call the
API for tiles at all.

The stamp itself comes from the API (`GET /api/galaxy/stamp`, one cheap
query) but is remembered on disk for `STAMP_TTL_SECONDS`, so a database
change shows up on the map within that long and busy traffic still costs
at most one stamp lookup per database per `STAMP_TTL_SECONDS`. When the
stamp changes, the old stamp's tiles are deleted.

The browser keeps its own copy of every tile too (`static/galaxymap3d.js`,
in `localStorage`, keyed by the same stamp), so a visitor panning back
over space they've already seen doesn't even reach this cache.

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

from apiclient import get_galaxy_stamp, get_galaxy_tiles
from stellarObjects.appconfig import load_config
from stellarObjects.galaxyViewport import parse_tile_key, tile_key

DEFAULT_CACHE_DIR = "/var/cache/planetgen/tiles"
"""str: Used when neither `PLANETGEN_TILE_CACHE_DIR` nor `config.json`'s
`tile_cache.dir` names a directory. Falls back to a `planetgen-tiles`
folder in the system temp directory when Apache can't create this one."""

STAMP_TTL_SECONDS = 60
"""int: How long a database's stamp is trusted before asking the API again
-- the longest a newly generated sector can take to appear on an open
map."""

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


def cache_dir():
    """
    The writable cache root, created if needed, or `None` when the disk
    cache is off or no candidate directory is writable.
    """
    if max_cache_bytes() == 0:
        return None
    configured = os.environ.get("PLANETGEN_TILE_CACHE_DIR") or _config().get("dir") or ""
    candidates = [configured] if configured else [
        DEFAULT_CACHE_DIR, os.path.join(tempfile.gettempdir(), "planetgen-tiles"),
    ]
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


def current_stamp(db, root=None):
    """
    The database's current content stamp: from `<db dir>/stamp.json` when
    that's younger than `STAMP_TTL_SECONDS`, else from the API (and then
    remembered). When the API reports a new stamp, every other stamp's
    cached tiles for this database are deleted.
    """
    if root is None:
        return get_galaxy_stamp(db)

    db_dir = _db_dir(root, db)
    stamp_path = os.path.join(db_dir, "stamp.json")
    remembered = _read_json(stamp_path)
    try:
        age = time.time() - os.path.getmtime(stamp_path)
    except OSError:
        age = None
    if (
        isinstance(remembered, dict)
        and _STAMP_RE.match(str(remembered.get("stamp", "")))
        and age is not None
        and 0 <= age < STAMP_TTL_SECONDS
    ):
        return remembered["stamp"]

    stamp = get_galaxy_stamp(db)
    if not _STAMP_RE.match(str(stamp)):
        return stamp
    _write_json(stamp_path, {"stamp": stamp})
    try:
        for entry in os.scandir(db_dir):
            if entry.is_dir() and entry.name != stamp:
                _remove_tree(entry.path)
    except OSError:
        pass
    return stamp


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


def fetch_tiles(db, tile_keys, density_key=None):
    """
    The requested tiles (and optional density cloud), from the disk cache
    where possible and from `GET /api/galaxy/tiles` for the rest, caching
    whatever the API returns.

    Args:
        db (str): The `?db=` value.
        tile_keys (list[str]): `level/ix/iy/iz` keys.
        density_key (str or None): A tile key to anchor a density cloud on.

    Returns:
        dict: `stamp` (see `current_stamp` -- the browser keys its own
            cache by it), `tiles` (`{key: {"placed", "planned"}}`),
            `density` (`{"key", "points"}` or `None`), `edge_pc`,
            `has_shape`, and `cached` (how many of the requested parts
            came from disk -- for diagnostics).

    Raises:
        TileRequestError: On a malformed request.
        apiclient.ApiError / NotFoundError: If the API is needed and fails.
    """
    tile_keys, density_key = validate_request(tile_keys, density_key)
    root = cache_dir()
    stamp = current_stamp(db, root)
    stamp_dir = os.path.join(_db_dir(root, db), stamp) if root and _STAMP_RE.match(str(stamp)) else None

    tiles = {}
    density = None
    meta = None
    if stamp_dir:
        meta = _read_json(os.path.join(stamp_dir, "meta.json"))
        for key in tile_keys:
            cached = _read_json(os.path.join(stamp_dir, _tile_filename("t", key)))
            if isinstance(cached, dict):
                tiles[key] = cached
        if density_key:
            cached = _read_json(os.path.join(stamp_dir, _tile_filename("d", density_key)))
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

        if stamp_dir:
            _write_json(os.path.join(stamp_dir, "meta.json"), meta)
            for key, value in fetched_tiles.items():
                if key in missing:
                    _write_json(os.path.join(stamp_dir, _tile_filename("t", key)), value)
            if fetched_density is not None:
                _write_json(os.path.join(stamp_dir, _tile_filename("d", density_key)), fetched_density["points"])
            if random.random() < PRUNE_PROBABILITY:
                prune(root)

        tiles.update(fetched_tiles)
        if fetched_density is not None:
            density = fetched_density

    return {
        "stamp": stamp,
        "tiles": {key: tiles[key] for key in tile_keys if key in tiles},
        "density": density,
        "edge_pc": meta.get("edge_pc"),
        "has_shape": bool(meta.get("has_shape")),
        "cached": cached_count,
    }


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
