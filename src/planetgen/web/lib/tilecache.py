# planetgen/web/lib/tilecache.py

"""
On-disk cache of the 3D Galaxy Map's sector-block nodes (and its drill-down stages), kept by the web layer so
the map doesn't hit the API (and through it, the database) every time the
camera moves.

A node's contents (see `planetgen.db.nodes.galaxy_nodes`) depend only on its key and on the
database's contents, which `planetgen.db.query.galaxy_tile_stamp` summarizes as a
short "stamp" (it leaves the release number out: a release that changes neither `query.TILE_FORMAT` nor the generator
epoch keeps every stored node, PERF.133). Nodes are cached as `<cache dir>/<db>/<generation>/
n<m>_<shell>_<band>_<wedge>.bin` (the finished, trimmed and packed bytes, see `trim_node` and `tilepack`), where the
generation is the stamp at the last time every node went stale. `fetch_nodes` asks the API only for the nodes it doesn't
already have, and when every node is cached it doesn't call the API for
nodes at all.

Freshness is checked at most once per `STAMP_TTL_SECONDS` per database,
with one cheap API call (`GET /api/galaxy/changes`, see
`planetgen.db.query.galaxy_changes`), which reads the schema-v27 `modified_at`
columns and answers which nodes changed since the last check. Only those
node files are deleted, so renaming a sector refetches the few dozen nodes
over it rather than the whole map. A change that can't be pinned to
nodes (a deleted sector, a re-planned galaxy, a new format or generator epoch)
answers `full`, and a new generation starts with nothing cached.

A request is answered by joining those stored bytes, never by unpacking and
packing the nodes again (`fetch_nodes_wire`, MAP.159); Apache compresses the
response (`mod_brotli`, else `mod_deflate`, see `examples/apache`). Sending
the nodes' own gzip copies as the members of one gzip stream was tried and
dropped: Chromium decodes only the first member.

The browser keeps its own copy of every node too (`static/galaxymap3d.js`,
in IndexedDB, keyed by the generation, MAP.158). `stamp.json` remembers the
last few checks' changed nodes (`history`), and `fetch_nodes` hands them
to the browser whenever its stamp is out of date, so the browser drops
only those nodes too.

The name `tilecache` is the module's from when the map was cut into cube tiles (MAP.262 retired them); the cache
directory and the `tile_cache` settings keep it, since a running site's folder and settings carry it.

Everything here fails open: an unwritable or missing cache directory, a
corrupt file, or a race with another request just means the node is
fetched from the API again, never an error page.
"""

import hashlib
import json
import os
import random
import re
import tempfile
import threading
import time

from planetgen.util import log
from planetgen.web.lib import apiclient, singleflight, tilepack
from planetgen.web.lib.apiclient import get_galaxy_changes, get_galaxy_nodes, get_galaxy_stage, get_galaxy_window
from planetgen.web.lib.privatedir import ensure_private_dir
from planetgen.util.settings import get_settings
from planetgen.galaxy.drill import format_drill_key, parse_node_key, parse_stage_key
from planetgen.galaxy.window import DEFAULT_WINDOW, check_width, crop_stages, parse_window_key, window_cells, window_parents

DEFAULT_CACHE_DIR = "/var/cache/planetgen/tiles"
"""str: Used when neither `PLANETGEN_TILE_CACHE_DIR` nor `config.json`'s
`tile_cache.dir` names a directory. Falls back to a `planetgen-tiles`
folder in the system temp directory when Apache can't create this one."""

STAMP_TTL_SECONDS = 60
"""int: How long a database's stamp is trusted before asking the API again
-- the longest a newly generated or edited sector can take to appear on
an open map."""

BUSY_KEEP_SECONDS = 600
"""int: While the galaxy changes faster than nodes can be refetched (a fill
adds sectors and stars every second, PERF.34), a cached generation is kept
and served at most this long past its start, instead of being thrown away at
every check, which made every map request recompute its nodes on a database
already under load. Then one full refresh brings it up to date."""

FAILED_CHECK_RETRY_SECONDS = 15
"""int: After the freshness check itself fails (a database too busy to answer
in time), the cache is served as it is and the check is not tried again for
this long."""

_failed_checks = {}
"""dict: `db -> time.monotonic()` of that database's last failed check."""

HISTORY_MAX_KEYS = 1500
"""int: Most changed-node keys `stamp.json`'s `history` keeps, across all
its entries. Oldest entries go first; a browser whose stamp is older than
what's left just drops all its nodes."""

HISTORY_MAX_ENTRIES = 32
"""int: Most checks `history` remembers."""

PRUNE_PROBABILITY = 0.05
"""float: Chance that a request which wrote new nodes also checks the
cache's total size (a directory walk, so not on every request)."""

MAX_NODES_PER_REQUEST = 128
"""int: Mirrors `planetgen.db.nodes.MAX_NODES_PER_REQUEST` -- rejected here before anything is read or fetched."""

_STAMP_RE = re.compile(r"^[0-9a-f]{16}$")
_ROUND_DIGITS = 4

_STAR_KEYS = ("id", "x", "y", "z", "luminosity_sol", "temperature_k", "radius_sol", "star_type", "system_id", "name")
"""tuple: The only fields of a node's star the page's script reads (MAP.157;
`nodeStars`, `setStars` and `starShown` in `static/galaxymap3d.js`)."""

_KEEP_SECTIONS = ("aggregate", "clouds", "stars", "generated", "points", "worlds", "placed", "incomplete", "name", "parent")
"""tuple: The sections of a node the page reads. `incomplete` (PERF.76) names the pieces the database was too busy
to give; such a node is served but never stored."""



class TileRequestError(ValueError):
    """A malformed node request from the map's own client-side JS (bad key,
    too many keys) -- reported as a 400, not an API failure."""


def _config():
    return get_settings().tile_cache


def budget_bytes():
    """The disk caches' budget in all, bytes (`cache.max_mb`). `0` disables the disk cache."""
    return max(0, int(get_settings().cache.max_mb * 1024 * 1024))


def max_cache_bytes():
    """The map cache's size budget, bytes: the budget less what the object index takes, which is kept first
    (PERF.140; at least one byte, so a budget the index fills leaves the cache full, not off). `0` disables the
    disk cache."""
    budget = budget_bytes()
    if budget == 0:
        return 0
    from planetgen.db import objectindex
    return max(1, budget - objectindex.reserved_bytes())


def configured_cache_dir():
    """
    The directory the cache is meant to live in (`PLANETGEN_TILE_CACHE_DIR`,
    then `tile_cache.dir`, then `DEFAULT_CACHE_DIR`), or `None` when the
    disk cache is off. Doesn't create or check anything --
    `examples/apache/create-cache-dir.sh` asks this where to create it.
    """
    if budget_bytes() == 0:
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
        # refused (they could plant nodes, or read the cache).
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


def _delete_nodes(generation_dir, node_keys):
    """Deletes the stored files of `node_keys`; returns the keys that were stored (and are gone now)."""
    removed = []
    for key in node_keys:
        try:
            path = os.path.join(generation_dir, _node_filename(key))
        except ValueError:
            continue
        try:
            os.unlink(path)
            removed.append(key)
        except OSError:
            pass
    return removed


_owed = {}
_owed_lock = threading.Lock()


def _owe_rebuild(db, generation, keys):
    """Remembers the stored nodes a change deleted, for `take_rebuilds` (PERF.133)."""
    with _owed_lock:
        had_generation, had = _owed.get(db, (generation, []))
        _owed[db] = (generation, (had if had_generation == generation else []) + list(keys))


def take_rebuilds(db):
    """
    The stored nodes of `db` that changes deleted since the last call, so they can be built again at once, in the
    background (PERF.133, `tileprefetch.rebuild`), not on the next visit that wants them.

    Returns:
        tuple: `(generation, keys)`, or `(None, [])` when no stored node was deleted.
    """
    with _owed_lock:
        generation, keys = _owed.pop(db, (None, []))
    return generation, list(dict.fromkeys(keys))


def _stage_filename(key):
    """A drill-down stage's cache file: `sgalaxy.json` for the galaxy,
    else `s<m>_<shell>_<band>_<wedge>.json`."""
    if key == "galaxy":
        return "sgalaxy.json"
    return "s" + format_drill_key(parse_stage_key(key)).replace(".", "_") + ".json"


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
        total += len(entry["nodes"])
        if total > HISTORY_MAX_KEYS:
            break
        kept.append(entry)
    return list(reversed(kept))


def current_stamp(db, root=None):
    """
    The database's current stamp, generation and change history: from
    `<db dir>/stamp.json` when that's younger than `STAMP_TTL_SECONDS`,
    else from `GET /api/galaxy/changes` (and then remembered). A check
    that finds changes deletes just the changed nodes, or on a `full`
    answer starts a new generation and deletes every other one.

    Returns:
        dict: `stamp`, `generation` (`None` when there's no disk cache),
            and `history` (a list of `{"from", "to", "nodes"}`, oldest
            first: each check that found changes, and the node keys it
            found).
    """
    if root is None:
        changes = get_galaxy_changes(db)
        return {"stamp": changes.get("tile_stamp") or changes["stamp"], "generation": None, "history": []}

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
    stamp = changes.get("tile_stamp") or changes.get("stamp")
    full = changes.get("tile_full", changes.get("full"))
    if changes.get("busy") and remembered is not None:
        generation_age = time.time() - float(remembered.get("started") or 0)
        if 0 <= generation_age < BUSY_KEEP_SECONDS:
            # Keep the stored state, so the next check still sees everything
            # that changed since; just look again in a minute.
            _write_json(stamp_path, remembered)
            return remembered
    if not _STAMP_RE.match(str(stamp)):
        return {"stamp": stamp, "generation": None, "history": []}

    stale_stages = []
    stale_nodes = []
    started = time.time()
    if remembered is not None and not full:
        generation = remembered["generation"]
        started = remembered.get("started") or started
        history = remembered["history"]
        if stamp != remembered["stamp"]:
            stale_stages = [str(key) for key in changes.get("stages") or []]
            stale_nodes = [str(key) for key in changes.get("nodes") or []]
            history = _trim_history(history + [
                {"from": remembered["stamp"], "to": stamp, "nodes": stale_nodes}])
    else:
        generation = stamp
        history = []

    info = {"stamp": stamp, "state": str(changes.get("state") or ""), "generation": generation, "history": history,
            "started": started}
    generation_dir = os.path.join(db_dir, generation)
    # Deleted both before and after the new stamp is written: a request
    # still working under the old stamp checks it before writing a node
    # (see `fetch_nodes`), so whichever side of the write it lands on, a
    # node fetched before the change doesn't survive it.
    removed_nodes = _delete_nodes(generation_dir, stale_nodes)
    _delete_stages(generation_dir, stale_stages)
    _write_json(stamp_path, info)
    removed_nodes = list(dict.fromkeys(removed_nodes + _delete_nodes(generation_dir, stale_nodes)))
    _delete_stages(generation_dir, stale_stages)
    if removed_nodes:
        _owe_rebuild(db, generation, removed_nodes)
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
    so it drops just their nodes: `None` when its stamp is current
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
    """Rounds every float in a JSON-able structure -- nodes are mostly
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
    """One star of a node as the page reads it (MAP.157): the fields in
    `_STAR_KEYS`, `star_type` as its class letter, positions to 3
    decimals, luminosity to 4 and radius to 3 significant digits and
    temperature to 10 K. No `system_id` when there is none, no `tier` for a star of the disc, the bulge and the rest, no `name`
    for a star that has none."""
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
    if star.get("tier"):
        trimmed["tier"] = star["tier"]  # 1 the halo, 2 a satellite galaxy: the first to dim (MAP.252)
    if star.get("name"):
        trimmed["name"] = str(star["name"])  # the map writes it beside the star (MAP.298); `tilepack` keeps it in the node's `names`
    return trimmed


def trim_node(node):
    """
    A node as the page's script reads it (MAP.260): its stars cut to `_STAR_KEYS` and rounded (`_trim_star`), every
    other float rounded to `_ROUND_DIGITS` decimals; what is kept is `_KEEP_SECTIONS`.
    """
    trimmed = {}
    for section in _KEEP_SECTIONS:
        value = node.get(section)
        if value is None:
            continue
        if section in ("stars", "generated"):
            trimmed[section] = [_trim_star(star) for star in value]
        else:
            trimmed[section] = _round_floats(value)
    return trimmed


def _node_bytes(node):
    """A node's finished packed bytes (`tilepack.pack_node` of `trim_node`)."""
    return tilepack.pack_node(trim_node(node))


def _incomplete(raw):
    """Whether a node's bytes say some piece of it is missing (PERF.76): it is served but not kept, so the next
    request builds it again."""
    return tilepack.is_incomplete(raw)


def _whole_node(raw):
    """Whether stored bytes are a whole packed node (a torn or foreign file is refetched,
    never sent)."""
    return tilepack.is_whole(raw)


def _node_filename(key):
    """A node's cache file: `n<m>_<shell>_<band>_<wedge>.bin`."""
    return "n" + format_drill_key(parse_node_key(key)).replace(".", "_") + ".bin"


def _flight_name(root, db, generation, kind, key):
    """The name a node or stage of `generation` goes by in `singleflight`
    (the cache directory is in it: two caches don't wait for each other)."""
    return f"{root}:{db}:{generation}:{kind}:{key}"


def _read_bytes(path):
    try:
        with open(path, "rb") as f:
            return f.read()
    except OSError:
        return None


AHEAD_MODE = 0o640
"""int: The mode of a node built ahead of any visit (the warm-up, the background build). The group-read bit is what
marks it: `prune` never deletes such a node, only the files a visit asked for (mode 0o600), least recently visited
first (PERF.133)."""


def _write_bytes(path, data, ahead=False):
    """Atomic write of raw bytes, like `_write_json`; `ahead` marks the file as built ahead (`AHEAD_MODE`).
    Returns whether it worked."""
    directory = os.path.dirname(path)
    tmp_path = None
    try:
        os.makedirs(directory, mode=0o750, exist_ok=True)
        fd, tmp_path = tempfile.mkstemp(dir=directory, prefix=".tmp-", suffix=".part")
        with os.fdopen(fd, "wb") as f:
            f.write(data)
        if ahead:
            os.chmod(tmp_path, AHEAD_MODE)
        os.replace(tmp_path, path)
        return True
    except OSError:
        if tmp_path:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
        return False


def _visited(path):
    """Records a visit to a stored file: its modification time is when it was last asked for, which is how `prune`
    orders what to drop."""
    try:
        os.utime(path)
    except OSError:
        pass


def validate_request(node_keys):
    """
    Canonicalizes and validates a node request.

    Returns:
        list[str]: The node keys with duplicates removed, in canonical form.

    Raises:
        TileRequestError: On a malformed key or too many keys.
    """
    try:
        keys = list(dict.fromkeys(format_drill_key(parse_node_key(key)) for key in node_keys))
    except ValueError as exc:
        raise TileRequestError(str(exc))
    if len(keys) > MAX_NODES_PER_REQUEST:
        raise TileRequestError(f"at most {MAX_NODES_PER_REQUEST} nodes per request")
    return keys


def _collect(db, keys, known_stamp, ahead=False):
    """
    What `fetch_nodes` and `fetch_nodes_wire` share: the requested nodes' stored bytes, from the disk cache where
    possible and from `GET /api/galaxy/nodes` (trimmed, then cached) for the rest.

    Returns:
        tuple: `(header, parts)`. `header` is the response without its
            nodes (`stamp`, `generation`, `edge_pc`, `has_shape`, `cached`
            and, only when `known_stamp` isn't current, `history`);
            `parts` is `[(key, packed bytes), ...]` in request order.

    `ahead` says the nodes are built for no visitor yet (the warm-up, the background build): what is built is marked
    to stay when the cache is full (`AHEAD_MODE`) and a node found stored is not counted as visited.
    """
    node_keys = validate_request(keys)
    root = cache_dir()
    info = current_stamp(db, root)
    stamp = info["stamp"]
    generation = info["generation"]
    generation_dir = os.path.join(_db_dir(root, db), generation) if root and generation else None

    stored = {}
    meta = None
    if generation_dir:
        meta = _read_json(os.path.join(generation_dir, "meta.json"))
        for key in node_keys:
            path = os.path.join(generation_dir, _node_filename(key))
            raw = _read_bytes(path)
            # A torn or foreign file is refetched, never sent.
            if _whole_node(raw):
                stored[key] = raw
                if not ahead:
                    _visited(path)

    cached_count = len(stored)
    missing = [key for key in node_keys if key not in stored]
    if log.debug_log_active():
        log.trace(f"node cache: {cached_count} of {len(node_keys)} node(s) found stored in database {db} "
                  f"(generation {generation}); {len(missing)} to build" + (" (ahead of any visitor)" if ahead else ""),
                  area="cache")

    if missing or not isinstance(meta, dict):
        # One build of each node at a time (PERF.38): what another request
        # is already fetching is waited for, not fetched again.
        owned, held = singleflight.claim(
            [_flight_name(root, db, generation, "n", key) for key in missing]) if generation_dir else (None, [])
        by_name = {_flight_name(root, db, generation, "n", key): key for key in missing}
        to_fetch = missing if owned is None else [by_name[name] for name in owned]
        made = {}
        try:
            if to_fetch or (not held and not isinstance(meta, dict)):
                fetched = get_galaxy_nodes(db, to_fetch)
                meta = {"edge_pc": fetched.get("edge_pc"), "has_shape": bool(fetched.get("has_shape"))}
                made = {key: _node_bytes(value) for key, value in (fetched.get("nodes") or {}).items()}

                # Another request may have found changes while this one was
                # fetching; what it fetched could be from before them.
                if generation_dir and (_read_remembered(_db_dir(root, db)) or {}).get("stamp") == stamp:
                    _write_json(os.path.join(generation_dir, "meta.json"), meta)
                    for key, raw in made.items():
                        if key in missing and not _incomplete(raw):
                            _write_bytes(os.path.join(generation_dir, _node_filename(key)), raw, ahead=ahead)
                    if random.random() < PRUNE_PROBABILITY:
                        prune(root)
        finally:
            if owned:
                singleflight.release(owned)

        stored.update(made)
        if held:
            def read_node(name):
                raw = _read_bytes(os.path.join(generation_dir, _node_filename(by_name[name])))
                return raw if _whole_node(raw) else None

            for name, raw in singleflight.wait(held, read_node).items():
                stored[by_name[name]] = raw
                cached_count += 1
            late = [by_name[name] for name in held if by_name[name] not in stored]
            if late:
                # The other request never delivered (or its node belongs to
                # a generation that has moved on): build them after all.
                fetched = get_galaxy_nodes(db, late)
                for key, value in (fetched.get("nodes") or {}).items():
                    stored[key] = _node_bytes(value)
                if not isinstance(meta, dict):
                    meta = {"edge_pc": fetched.get("edge_pc"), "has_shape": bool(fetched.get("has_shape"))}
            if not isinstance(meta, dict):
                meta = _read_json(os.path.join(generation_dir, "meta.json"))
                if not isinstance(meta, dict):
                    fetched = get_galaxy_nodes(db, [])
                    meta = {"edge_pc": fetched.get("edge_pc"), "has_shape": bool(fetched.get("has_shape"))}

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
    return header, [(key, stored[key]) for key in node_keys if key in stored]


def fetch_nodes_wire(db, node_keys, known_stamp=None, ahead=False):
    """
    The requested nodes (MAP.260) as the bytes to send (`/galaxy/nodes`): the header as JSON and the stored packed nodes
    one after the other (`tilepack.pack_response`), from the disk cache where possible and from `GET /api/galaxy/nodes`
    for the rest, caching whatever the API returns (trimmed, see `trim_node`). The header and the stamp are the
    nodes' (`current_stamp`), and `history` lists the changed node keys in `nodes`.

    Args:
        db (str): The `?db=` value.
        node_keys (list[str]): `m.shell.band.wedge` keys.
        known_stamp (str or None): The stamp the browser's own cache is at, when it has one.
        ahead (bool): The nodes are built for no visitor yet (see `_collect`).

    Raises:
        TileRequestError: On a malformed request.
        apiclient.ApiError / NotFoundError: If the API is needed and fails.
    """
    header, parts = _collect(db, node_keys, known_stamp, ahead)
    return tilepack.pack_response(header, parts)


def fetch_nodes(db, node_keys, known_stamp=None, ahead=False):
    """`fetch_nodes_wire`'s payload with every node unpacked: `nodes` (`{key: {...}}`) beside the header."""
    header, parts = _collect(db, node_keys, known_stamp, ahead)
    return {**header, "nodes": {key: tilepack.unpack_node(raw) for key, raw in parts}}


def fetch_stage(db, at=None):
    """
    One Galaxy Map drill-down stage (`planetgen.db.query.galaxy_stage`), from the disk
    cache when it's there, else from `GET /api/galaxy/stage` (and then
    cached). A sector change deletes just its chain's stages (`stages` in
    `/api/galaxy/changes`), the same way nodes are refreshed.

    Args:
        db (str): The `?db=` value.
        at (str or None): A block key `m.shell.band.wedge`, or `None` for
            the galaxy.

    Returns:
        dict: The stage, plus `stamp`.

    Raises:
        TileRequestError: On a malformed key.
        apiclient.ApiError / NotFoundError: If the API is needed and fails.
    """
    try:
        key = format_drill_key(parse_stage_key(at)) if at else "galaxy"
    except ValueError as exc:
        raise TileRequestError(str(exc))
    root = cache_dir()
    info = current_stamp(db, root)
    stamp = info["stamp"]
    generation = info["generation"]
    generation_dir = os.path.join(_db_dir(root, db), generation) if root and generation else None
    path = os.path.join(generation_dir, _stage_filename(key)) if generation_dir else None

    stage = _read_json(path) if path else None
    if isinstance(stage, dict):
        _visited(path)
    else:
        name = _flight_name(root, db, generation, "s", key)
        owned, _held = singleflight.claim([name]) if path else ([name], [])
        try:
            if not owned:
                # Another request is building this stage (PERF.38): wait for its file.
                stage = singleflight.wait([name], lambda _name: _read_json(path)).get(name)
            if not isinstance(stage, dict):
                stage = get_galaxy_stage(db, None if key == "galaxy" else key)
                if path and (_read_remembered(_db_dir(root, db)) or {}).get("stamp") == stamp:
                    _write_json(path, stage)
        finally:
            if owned and path:
                singleflight.release(owned)
    return {**stage, "stamp": stamp}


def fetch_window(db, at, g=DEFAULT_WINDOW):
    """
    One Galaxy Map data window (MAP.151, `galaxy.window`): `g` cells across
    around the cell `at`. Cells bigger than a sector come from the stored
    block totals in one query (`planetgen.db.query.galaxy_window`, MAP.151b); a window
    of sectors is cut from the stages of the parent blocks it touches,
    which come from `fetch_stage`, so a warm cache serves it without a
    query and it shares its files with every other window and stage in the
    neighbourhood.

    Args:
        db (str): The `?db=` value.
        at (str): A cell key `m.shell.band.wedge` of any level, `1.shell.band.wedge` for a sector.
        g (int): Cells across, one of `window.WINDOW_WIDTHS`.

    Returns:
        dict: `window.crop_stages`' result, plus `stamp`.

    Raises:
        TileRequestError: On a malformed key or width.
        apiclient.ApiError / NotFoundError: If the API is needed and fails.
    """
    try:
        centre = parse_window_key(at)
        width = check_width(g)
    except ValueError as exc:
        raise TileRequestError(str(exc))
    if centre.m > 1:
        # Cells bigger than a sector: the stored block totals, in one query.
        stamp = current_stamp(db, cache_dir())["stamp"]
        return {**get_galaxy_window(db, format_drill_key(centre), width), "stamp": stamp}
    # The sector edge comes with every stage; the galaxy's is always cached.
    galaxy = fetch_stage(db, None)
    edge_pc = galaxy["edge_pc"]
    stamp = galaxy.get("stamp")
    stages = {}
    for parent in window_parents(window_cells(centre, edge_pc, width)):
        key = format_drill_key(parent) if parent else "galaxy"
        stages[key] = galaxy if parent is None else fetch_stage(db, key)
    return {**crop_stages(centre, width, stages, edge_pc), "stamp": stamp}


KEPT_FILES = ("stamp.json", "meta.json", ".full")
"""tuple: Files in the cache that `prune` never counts or deletes."""


def cache_files(root):
    """
    Every cache file under `root` that counts against the budget, as `(path, size, mtime, ahead)`: the stamp files
    and the full flag are left out, and `ahead` is whether `AHEAD_MODE`'s bit marks the file as built ahead.
    """
    found = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for name in filenames:
            if name in KEPT_FILES:
                continue
            path = os.path.join(dirpath, name)
            try:
                info = os.stat(path)
            except OSError:
                continue
            found.append((path, info.st_size, info.st_mtime, bool(info.st_mode & 0o040)))
    return found


def cache_is_full(root=None, max_bytes=None):
    """Whether the files under `root` (default `cache_dir()`) already fill the budget (default `max_cache_bytes()`):
    the warm-up builds nothing more then (PERF.133). A cache that is off is full."""
    root = root or cache_dir()
    if max_bytes is None:
        max_bytes = max_cache_bytes()
    if root is None or max_bytes <= 0:
        return True
    return sum(size for _path, size, _mtime, _ahead in cache_files(root)) >= max_bytes


def prune(root, max_bytes=None):
    """
    Brings the cache under 80% of `max_bytes` (default `max_cache_bytes()`) when it is over `max_bytes`, and then
    removes any directories that left empty (PERF.133): only the files a visit asked for are deleted, the
    least recently visited first (a visit sets a file's modification time); the nodes built ahead
    (`AHEAD_MODE`) stay, so a cache that holds nothing but those stays over the budget and the warm-up stops
    (`cache_is_full`). The stamp files are skipped.
    """
    if max_bytes is None:
        max_bytes = max_cache_bytes()
    files = cache_files(root)
    total = sum(size for _path, size, _mtime, _ahead in files)
    if total <= max_bytes:
        return
    from planetgen.web.lib import tileprefetch
    tileprefetch.mark_full(root)     # the background build waits while the cache is this full (PERF.132)
    target = max_bytes * 0.8
    for path, size, _mtime, _ahead in sorted((entry for entry in files if not entry[3]), key=lambda entry: entry[2]):
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
