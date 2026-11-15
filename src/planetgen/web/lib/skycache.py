# planetgen/web/lib/skycache.py

"""
The packed skies of the System View, kept on disk beside the Galaxy Map's nodes (VIEW.31; the cache is PERF.129's, its admin page PERF.134's).

A sky (`api/system_sky.py`, `skypack.py`) is the stars one system sees at a limiting magnitude, about 0.65 MB at 7.5, and the nebulae and
remnants beside them. It depends on the stored stars, so it is good while the galaxy's tile stamp stands (`tilecache.current_stamp`: the sectors, systems
and bright stars stored, and the generator epoch) and no longer. It does not depend on the orbit clock: the page moves the stars on from the list's
`epoch_unix`, so a sky stays right across orbital updates.

A sky is stored as `<cache dir>/<db>/<generation>/k<system>_<limit in hundredths>_<version key>.bin` (`skypack.version_key`: the stamp and
`SKY_VERSION`) in the generation directory of `tilecache`, so a new generation drops every sky with the nodes, `prune` counts and drops them with
the rest (least recently asked for first), and clearing the tile cache clears them. A sky built for a newer version key removes the older copies of
its own system and limit. The cache fails open like the node cache: an unwritable folder, a torn file or a race just means the API builds it again.

`fetch_sky(db, system_id, limit)` returns `(bytes, header)`: the packed list and `{"version_key", "source": "cache" | "built", "count", "sprites", "report"}`
where `report` is the tier state (`docs/design/interface/system-sky-wire-format.md` section 4; read from the list's own `meta`, so a hit is as full as a build).
"""

import base64
import glob
import os
import random
import re

from planetgen.util import log
from planetgen.web.lib import apiclient, singleflight, skypack, tilecache

NAME_PREFIX = "k"
"""str: What a sky's file name starts with (`n` is a node's, `s` a stage's)."""

BAND_PREFIX = "b"
"""str: What a glow image's file name starts with: `b<key>.bin`, the key being the band's cache key (a hash of its cell, the galaxy, the limit and every constant
it depends on), so the name says what the file holds and a changed galaxy is a new name."""

_BAND_KEY = re.compile(r"^[0-9a-f]{16,64}$")


_BAND_SIZE = (256, 128)
_BAND_BYTES = _BAND_SIZE[0] * _BAND_SIZE[1] * 6


def sky_filename(system_id, limit, key):
    """The cache file of the sky of `system_id` at `limit` (a float from `skypack.check_limit`) for the version key `key`."""
    return f"{NAME_PREFIX}{system_id}_{int(round(limit * 100))}_{key}.bin"


def _pattern(system_id, limit):
    return f"{NAME_PREFIX}{glob.escape(str(system_id))}_{int(round(limit * 100))}_*.bin"


def _current_directory(db):
    """`(cache root, generation directory, version key)` of the current generation of the tile cache for `db`; the directory and key are `None` when
    there is no cache, or the galaxy's stamp cannot be read."""
    root = tilecache.cache_dir()
    if not root:
        return root, None, None
    try:
        info = tilecache.current_stamp(db, root)
    except apiclient.ApiError:
        info = None
    if not (info and info.get("generation")):
        return root, None, None
    return root, os.path.join(tilecache._db_dir(root, db), info["generation"]), skypack.version_key(info["stamp"])


def fetch_sky(db, system_id, limit=None):
    """
    The packed sky of the system `system_id` (a printed ID) at `limit` (text or number; the default when `None`), from the disk cache when it is stored for
    the current version key, else from `GET /api/systems/<id>/sky` (and then stored).

    Returns:
        tuple: `(bytes, header)`; `header` is `{"version_key", "source", "count", "sprites", "report"}`.

    Raises:
        skypack.SkyLimitError: A limit out of range.
        apiclient.NotFoundError / ApiError: If the API is needed and fails.
    """
    limit = skypack.check_limit(limit)
    root, directory, key = _current_directory(db)
    if directory:
        path = os.path.join(directory, sky_filename(system_id, limit, key))
        raw = tilecache._read_bytes(path)
        if skypack.is_whole(raw):
            tilecache._visited(path)
            return raw, _describe(raw, key, "cache")
    name = f"{root}:{db}:sky:{system_id}:{limit}" if directory else None
    owned = []
    if name:
        owned, _held = singleflight.claim([name])
        if not owned:
            raw = singleflight.wait([name], lambda _name: _stored(directory, system_id, limit, key)).get(name)
            if skypack.is_whole(raw):
                return raw, _describe(raw, key, "cache")
    try:
        built = apiclient.get_system_sky(db, system_id, limit)
        raw = base64.b64decode(built["data"])
        if not skypack.is_whole(raw):
            raise apiclient.ApiError("the API sent a sky that is not a whole packed list")
        if directory and (tilecache._read_remembered(tilecache._db_dir(root, db)) or {}).get("stamp") is not None:
            _store(directory, system_id, limit, built["version_key"], raw, root)
        return raw, _describe(raw, built["version_key"], "built")
    finally:
        if owned:
            singleflight.release(owned)


def _stored(directory, system_id, limit, key):
    raw = tilecache._read_bytes(os.path.join(directory, sky_filename(system_id, limit, key)))
    return raw if skypack.is_whole(raw) else None


def _describe(raw, key, source):
    meta = skypack.unpack_sky(raw)
    info = meta["meta"]
    report = info.get("report") or {}
    return {"version_key": info.get("version_key") or key, "source": source, "count": int(meta["count"]), "sprites": len(info.get("nebulae") or []),
            "report": report}


def _store(directory, system_id, limit, key, raw, root):
    """Writes the sky and removes the older copies of the same system and limit; sometimes checks the cache's size."""
    path = os.path.join(directory, sky_filename(system_id, limit, key))
    if not tilecache._write_bytes(path, raw):
        return
    for other in glob.glob(os.path.join(glob.escape(directory), _pattern(system_id, limit))):
        if other != path:
            try:
                os.unlink(other)
            except OSError:
                pass
    if random.random() < tilecache.PRUNE_PROBABILITY:
        try:
            tilecache.prune(root)
        except OSError as exc:
            log.debug(f"Sky cache: not pruned ({exc}).")


def fetch_band(db, system_id, limit=None, key=None):
    """
    The glow image of the sky of `system_id` at `limit` (VIEW.13): the bytes (little-endian half floats, `rgb16f`) and `{"key", "width", "height", "format",
    "source"}`. A band is named by its key, which `meta.band.key` of the sky gives, so a stored `b<key>.bin` is served without the API; otherwise the
    API works it out (once for each 50 parsec cell and limit, in its process) and the file is stored under the key it answers.

    Raises:
        skypack.SkyLimitError: A limit out of range.
        ValueError: A key that is not a hash.
        apiclient.NotFoundError / ApiError: If the API is needed and fails.
    """
    limit = skypack.check_limit(limit)
    if key is not None and not _BAND_KEY.match(key):
        raise ValueError("key must be a hash")
    root, directory, _version = _current_directory(db)
    if directory and key:
        path = os.path.join(directory, f"{BAND_PREFIX}{key}.bin")
        raw = tilecache._read_bytes(path)
        if raw is not None and len(raw) == _BAND_BYTES:
            tilecache._visited(path)
            return raw, {"key": key, "width": _BAND_SIZE[0], "height": _BAND_SIZE[1], "format": "rgb16f", "source": "cache"}
    name = f"{root}:{db}:band:{system_id}:{limit}" if directory else None
    owned = []
    if name:
        owned, _held = singleflight.claim([name])
    try:
        built = apiclient.get_system_sky_band(db, system_id, limit)
        raw = base64.b64decode(built["data"])
        if len(raw) != built["width"] * built["height"] * 6:
            raise apiclient.ApiError("the API sent a glow image of the wrong size")
        if directory and _BAND_KEY.match(str(built["key"])):
            tilecache._write_bytes(os.path.join(directory, f"{BAND_PREFIX}{built['key']}.bin"), raw)
        return raw, {"key": built["key"], "width": built["width"], "height": built["height"], "format": built["format"], "source": "built"}
    finally:
        if owned:
            singleflight.release(owned)


def clear(db=None, root=None):
    """Removes every stored sky (of the database `db`, or of all when `None`) and returns how many files went: the clear the admin page of PERF.134 calls
    beside the nodes'."""
    root = root or tilecache.cache_dir()
    if not root:
        return 0
    removed = 0
    top = tilecache._db_dir(root, db) if db is not None else root
    for dirpath, _dirnames, filenames in os.walk(top):
        for name in filenames:
            if name.endswith(".bin") and ((name.startswith(NAME_PREFIX) and "_" in name) or name.startswith(BAND_PREFIX)):
                try:
                    os.unlink(os.path.join(dirpath, name))
                    removed += 1
                except OSError:
                    pass
    return removed
