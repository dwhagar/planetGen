# planetgen/web/lib/tileprefetch.py

"""
Builds the Galaxy Map's tiles around a visit before the visitor asks for them (PERF.132; design in
docs/design/cache-warm-up.md).

A request for tiles (`GET /galaxy/tiles`, `schedule` runs when its response has gone) names the view the visitor is
in, and the page says where its camera and target are (`at`). The tiles the next move needs are (1) the next zoom
steps: the tiles of the cut for the camera `ZOOM_STEPS` fade bands closer along the same line to the target
(`zoom_candidates`, the same cut the map asks for, nearest step first), and (2) the neighbours of the requested tiles
at their own level and the tiles of the next level inside them (`candidates`, nearest the camera first: a request lists
its tiles nearest first), for a visitor who pans instead. Those that are not stored are queued as low-priority
background jobs on the bulk queue (never the interactive one, so an admin's work is never behind them), a chunk of
`CHUNK` tiles to a job, and a worker builds them into the same disk cache a visitor's request fills (`build_tiles`,
`tilecache.fetch_tiles`).

It stays small by rules, all named settings of `tile_cache`: at most `prefetch_per_visit` tiles queued for one request,
`prefetch_per_minute` in all, nothing while `prefetch_max_queued` jobs wait on the bulk queue, nothing while the
cache is full (`prune` has just deleted files, `full`), nothing during an update; a tile already stored or already
claimed by an earlier visit (a Redis key per tile and generation, `claim`) is not queued again, so a visit to an
area that is built queues nothing and a second visitor to the same area causes no build. Everything fails open: a
Redis that does not answer or a cache that cannot be read just means no prefetch, never an error for the visitor.
"""

import itertools
import os
import time

import redis

from planetgen.admin import updateguard
from planetgen.galaxy.viewport import DISTANCE_CUT_FADE_BAND, TILE_MAX_LEVEL, fetch_cut_tiles, parse_tile_key, tile_key
from planetgen.queue import redisclient, redisqueue, workers
from planetgen.util import log
from planetgen.util.settings import get_settings
from planetgen.web.lib import tilecache

CHUNK = 16
"""int: Tiles one queued job builds."""

ZOOM_STEPS = 4
"""int: How many zoom steps ahead (each one fade band closer, `viewport.DISTANCE_CUT_FADE_BAND`) are built."""

CLAIM_SECONDS = 900
"""int: How long a tile stays claimed once queued: a job that never ran frees the tile after this."""

MINUTE_KEY = "planetgen:tileprefetch:minute:"
"""str: The per-minute counter's key prefix (the minute's number follows)."""

CLAIM_KEY = "planetgen:tileprefetch:claim:"
"""str: A claimed tile's key prefix (`<database hash>:<generation>:<tile>` follows)."""

FULL_FILE = ".full"
"""str: The file `prune` leaves in the cache root when it had to delete files; its age says whether the cache is full."""

FULL_SECONDS = 600
"""int: How long after a prune that deleted files the cache counts as full."""

JOB = "planetgen.web.lib.tileprefetch.build_tiles"
"""str: The function a worker runs for a queued chunk (a string, the queue carries JSON)."""

SOCKET_SECONDS = 0.5
"""float: How long one Redis call may take before it counts as unreachable."""


def _redis():
    return redisclient.client(decode_responses=True, timeout=SOCKET_SECONDS, retries=1)


def candidates(keys):
    """
    The tiles worth building after a request for `keys`, in order: for each requested tile (nearest first), its eight
    children (the next level inside it) and then its neighbours at its own level (up to 26, inside the root cube).
    Neither the requested tiles nor a tile named twice appear.

    Yields:
        str: Tile keys.
    """
    seen = set(keys)
    for key in keys:
        level, ix, iy, iz = parse_tile_key(key)
        found = []
        if level < TILE_MAX_LEVEL:
            found += [tile_key(level + 1, 2 * ix + dx, 2 * iy + dy, 2 * iz + dz)
                      for dx, dy, dz in itertools.product((0, 1), repeat=3)]
        span = 2 ** level
        for dx, dy, dz in itertools.product((-1, 0, 1), repeat=3):
            if (dx, dy, dz) != (0, 0, 0) and all(0 <= index < span for index in (ix + dx, iy + dy, iz + dz)):
                found.append(tile_key(level, ix + dx, iy + dy, iz + dz))
        for tile in found:
            if tile not in seen:
                seen.add(tile)
                yield tile


def zoom_candidates(camera, look, bounds, steps=ZOOM_STEPS):
    """
    The tiles of the next `steps` zoom steps: the cut (`viewport.fetch_cut_tiles`) for the camera one fade band
    closer to `look` along its line, then two bands, and so on; each tile once, the nearest step first.

    Args:
        camera (tuple): The camera, parsecs.
        look (tuple): The target it looks at, parsecs.
        bounds (tuple): The box the cut covers (`galaxymap3d.tile_bounds_pc`).

    Yields:
        str: Tile keys.
    """
    seen = set()
    for step in range(1, steps + 1):
        scale = DISTANCE_CUT_FADE_BAND ** -step
        nearer = tuple(target + (position - target) * scale for position, target in zip(camera, look))
        for key in fetch_cut_tiles(nearer, look, bounds):
            if key not in seen:
                seen.add(key)
                yield key


def wanted(keys, at, bounds):
    """The tiles to build after a request for `keys` from the camera and target `at` (six numbers, or `None`):
    the next zoom steps first (when `at` is given), then the neighbours and the next level (`candidates`); never a
    tile the request itself named."""
    asked = set(keys)
    if at is not None:
        for key in zoom_candidates(tuple(at[:3]), tuple(at[3:]), bounds):
            if key not in asked:
                yield key
    yield from (key for key in candidates(keys) if key not in asked)


def mark_full(root):
    """Tells the prefetch the cache at `root` is full: `prune` calls this after it deleted files."""
    try:
        with open(os.path.join(root, FULL_FILE), "w", encoding="utf-8"):
            pass
    except OSError:
        pass


def is_full(root):
    """Whether `prune` deleted files from the cache at `root` within `FULL_SECONDS`."""
    try:
        return time.time() - os.stat(os.path.join(root, FULL_FILE)).st_mtime < FULL_SECONDS
    except OSError:
        return False


def _claim_name(db, generation, key):
    return f"{CLAIM_KEY}{os.path.basename(tilecache._db_dir('x', db))}:{generation}:{key}"


def _bounds(db):
    """The box the map's cut covers for this database's galaxy (`galaxymap3d.tile_bounds_pc`)."""
    from planetgen.physics.units import ly_to_pc
    from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY
    from planetgen.web.lib import apiclient
    from planetgen.web.maps.galaxymap3d import tile_bounds_pc
    shape = apiclient.get_galaxy_shape(db)
    return tile_bounds_pc(shape["edge_pc"] if shape else ly_to_pc(DEFAULT_SECTOR_EDGE_LY), shape)


def schedule(db, keys, at=None):
    """
    Queues the tiles around a request for `keys` that are not stored (see the module's text for the rules).

    Args:
        db (str): The database the site shows.
        keys (list[str]): The request's tile keys, nearest first.
        at (tuple, optional): The page's camera and target, six numbers in parsecs.

    Returns:
        int: How many tiles were queued (0 when the prefetch is off, nothing is missing or a rule held it back).
    """
    settings = get_settings().tile_cache
    if settings.prefetch_per_visit <= 0 or settings.prefetch_per_minute <= 0 or not keys:
        return 0
    if updateguard.is_updating():
        return 0
    root = tilecache.cache_dir()
    if root is None or is_full(root):
        return 0
    generation = tilecache.current_stamp(db, root)["generation"]
    if not generation:
        return 0
    generation_dir = os.path.join(tilecache._db_dir(root, db), generation)
    bounds = _bounds(db) if at is not None else None
    missing = []
    for key in wanted(keys, at, bounds):
        if key in missing:
            continue
        if not os.path.exists(os.path.join(generation_dir, tilecache._tile_filename("t", key))):
            missing.append(key)
            if len(missing) >= settings.prefetch_per_visit:
                break
    if not missing:
        return 0
    try:
        connection = _redis()
        queue = redisqueue.queue(redisqueue.BULK_QUEUE, connection)
        if queue.count >= settings.prefetch_max_queued:
            return 0
        names = [_claim_name(db, generation, key) for key in missing]
        pipe = connection.pipeline()
        for name in names:
            pipe.set(name, "1", nx=True, ex=CLAIM_SECONDS)
        claimed = [key for key, won in zip(missing, pipe.execute()) if won]
        if not claimed:
            return 0
        minute = f"{MINUTE_KEY}{int(time.time() // 60)}"
        used = connection.incrby(minute, len(claimed))
        connection.expire(minute, 120)
        room = max(0, settings.prefetch_per_minute - (used - len(claimed)))
        for key in claimed[room:]:
            connection.delete(_claim_name(db, generation, key))
        claimed = claimed[:room]
        for start in range(0, len(claimed), CHUNK):
            queue.enqueue(JOB, args=(db, claimed[start:start + CHUNK]), job_timeout=600, result_ttl=0,
                          failure_ttl=600, description=f"prefetch {len(claimed[start:start + CHUNK])} tiles")
        if claimed:
            workers.ensure_workers(redisqueue.BULK_QUEUE, connection, redisqueue.redis_url(), count=1)
        return len(claimed)
    except (redis.RedisError, OSError, redisqueue.Unavailable) as exc:
        log.debug(f"Tile prefetch: not queued ({exc}).")
        return 0


REBUILD_MAX = 2048
"""int: Most tiles one change queues for a rebuild (coarsest first); the rest are built by the next visit that wants them."""


def rebuild(db, generation, keys):
    """
    Queues the building of tiles a change deleted from the stored generation (`tilecache.current_stamp`): the changed
    tiles that were stored and their coarser parents, coarsest first, so the next visit finds them built again
    (PERF.133). The rules are the prefetch's except the per-visit and per-minute counts and the queue depth (a change
    touches a few hundred tiles at most, and it is the data the visitors are looking at): nothing while the
    prefetch is off (`prefetch_per_minute` 0), during an update, or while the cache is full; a tile already claimed
    is not queued twice.

    Args:
        db (str): The database the site shows.
        generation (str): The generation the tiles belong to.
        keys (list[str]): The tile keys that were stored and were deleted.

    Returns:
        int: How many tiles were queued.
    """
    settings = get_settings().tile_cache
    if settings.prefetch_per_minute <= 0 or not keys or updateguard.is_updating():
        return 0
    root = tilecache.cache_dir()
    if root is None or is_full(root):
        return 0
    ordered = sorted(keys, key=lambda key: (parse_tile_key(key)[0], key))[:REBUILD_MAX]
    try:
        connection = _redis()
        queue = redisqueue.queue(redisqueue.BULK_QUEUE, connection)
        pipe = connection.pipeline()
        for key in ordered:
            pipe.set(_claim_name(db, generation, key), "1", nx=True, ex=CLAIM_SECONDS)
        claimed = [key for key, won in zip(ordered, pipe.execute()) if won]
        for start in range(0, len(claimed), CHUNK):
            queue.enqueue(JOB, args=(db, claimed[start:start + CHUNK]), job_timeout=600, result_ttl=0,
                          failure_ttl=600, description=f"rebuild {len(claimed[start:start + CHUNK])} tiles")
        if claimed:
            workers.ensure_workers(redisqueue.BULK_QUEUE, connection, redisqueue.redis_url(), count=1)
        return len(claimed)
    except (redis.RedisError, OSError, redisqueue.Unavailable) as exc:
        log.debug(f"Tile rebuild: not queued ({exc}).")
        return 0


def build_tiles(db, keys):
    """
    The queued job: builds the tiles `keys` into the disk cache the way a visitor's request does
    (`tilecache.fetch_tiles`), so what is stored is what a visitor is served. A tile stored since it was queued is
    not built again.

    Returns:
        dict: `tiles` asked and `built` (those that were not stored).
    """
    from planetgen.web import warmup
    app = warmup.build_app()
    with app.test_request_context("/"):
        result = tilecache.fetch_tiles(db, list(keys), ahead=True)
    return {"tiles": len(keys), "built": len(keys) - result["cached"]}
