# planetgen/web/lib/tileprefetch.py

"""
Builds the Galaxy Map's nodes around a visit before the visitor asks for them (PERF.132, MAP.260; design in
docs/design/operations/cache-warm-up.md).

A request for nodes (`GET /galaxy/nodes`, `schedule` runs when its response has gone) names the blocks the visitor's
view holds. The nodes the next move needs are, for each requested node (the request lists its nodes nearest first), its
children (the next level inside it) and then its parent and the parent's other children (`candidates`), for a visitor
who zooms in or out or pans. Those that are not stored are queued as low-priority background jobs on the bulk queue
(never the interactive one, so an admin's work is never behind them), a chunk of `CHUNK` nodes to a job, and a worker
builds them into the same disk cache a visitor's request fills (`build_nodes`, `tilecache.fetch_nodes`).

It stays small by rules, all named settings of `tile_cache`: at most `prefetch_per_visit` nodes queued for one request,
`prefetch_per_minute` in all, nothing while `prefetch_max_queued` jobs wait on the bulk queue, nothing while the
cache is full (`prune` has just deleted files, `full`), nothing during an update; a node already stored or already
claimed by an earlier visit (a Redis key per node and generation, `claim`) is not queued again, so a visit to an
area that is built queues nothing and a second visitor to the same area causes no build. Everything fails open: a
Redis that does not answer or a cache that cannot be read just means no prefetch, never an error for the visitor.
"""

import os
import time

import redis

from planetgen.admin import updateguard
from planetgen.galaxy.drill import drill_children, drill_parent, format_drill_key, parse_node_key
from planetgen.queue import redisclient, redisqueue, stopall, workers
from planetgen.util import log
from planetgen.util.settings import get_settings
from planetgen.web.lib import tilecache

CHUNK = 16
"""int: Nodes one queued job builds."""

NODE_SECONDS = 2.0
"""float: What the size of a chunk's work is counted at, for its time limit (PERF.144): the seconds a cold node is
expected to take to build, an estimate of the size and not a recorded time (none is recorded for nodes)."""

JOB_START_SECONDS = 30.0
"""float: What a chunk's start-up (the app a worker builds) is counted at, with `NODE_SECONDS` a node."""

MIN_JOB_SECONDS = 600
"""int: The shortest limit of a chunk, and the ceiling is `redisqueue.MAX_TIMEOUT_SECONDS`."""

CLAIM_SECONDS = 900
"""int: How long a node stays claimed once queued: a job that never ran frees the node after this."""

MINUTE_KEY = "planetgen:tileprefetch:minute:"
"""str: The per-minute counter's key prefix (the minute's number follows)."""

CLAIM_KEY = "planetgen:tileprefetch:claim:"
"""str: A claimed node's key prefix (`<database hash>:<generation>:<node>` follows)."""

FULL_FILE = ".full"
"""str: The file `prune` leaves in the cache root when it had to delete files; its age says whether the cache is full."""

FULL_SECONDS = 600
"""int: How long after a prune that deleted files the cache counts as full."""

JOB = "planetgen.web.lib.tileprefetch.build_nodes"
"""str: The function a worker runs for a queued chunk of nodes (a string, the queue carries JSON)."""

SOCKET_SECONDS = 0.5
"""float: How long one Redis call may take before it counts as unreachable."""


def _redis():
    return redisclient.client(decode_responses=True, timeout=SOCKET_SECONDS, retries=1)


def candidates(keys):
    """
    The nodes worth building after a request for `keys`, in order: for each requested node (nearest first), its
    children (the next level inside it: three shells of rows and wedges, a sector's own node for a level-3 block), then
    its parent and the parent's other children. Neither the requested nodes nor one named twice appear.

    Yields:
        str: Node keys.
    """
    seen = set(keys)
    for key in keys:
        try:
            block = parse_node_key(key)
        except ValueError:
            continue
        found = []
        if block.m > 1:
            found += [format_drill_key(child) for _shell, group in drill_children(block) for child in group]
        parent = drill_parent(block)
        if parent is not None:
            found.append(format_drill_key(parent))
            found += [format_drill_key(child) for _shell, group in drill_children(parent) for child in group]
        for node in found:
            if node not in seen:
                seen.add(node)
                yield node


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


def chunk_timeout(count):
    """The seconds RQ lets a job of `count` nodes run (PERF.144): `redisqueue.job_timeout` of the size of its work,
    `JOB_START_SECONDS` and `NODE_SECONDS` a node; at least `MIN_JOB_SECONDS`."""
    return redisqueue.job_timeout(JOB_START_SECONDS + NODE_SECONDS * count, floor=MIN_JOB_SECONDS)


def _enqueue_chunks(queue, db, keys, what):
    """Queues `keys` as jobs of `CHUNK` nodes each on `queue`, `what` ("prefetch" or "rebuild") in their description.
    Nothing is queued while everything is stopped (PERF.186)."""
    if stopall.is_stopped():
        return
    for start in range(0, len(keys), CHUNK):
        chunk = keys[start:start + CHUNK]
        queue.enqueue(JOB, args=(db, chunk), job_timeout=chunk_timeout(len(chunk)), result_ttl=0, failure_ttl=600,
                      description=f"{what} {len(chunk)} nodes")


def schedule(db, keys):
    """
    Queues the nodes around a request for `keys` that are not stored (`candidates`; see the module's text for the
    rules).

    Args:
        db (str): The database the site shows.
        keys (list[str]): The request's node keys, nearest first.

    Returns:
        int: How many nodes were queued (0 when the prefetch is off, nothing is missing or a rule held it back).
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
    missing = []
    for key in candidates(keys):
        if key in missing:
            continue
        if not os.path.exists(os.path.join(generation_dir, tilecache._node_filename(key))):
            missing.append(key)
            if len(missing) >= settings.prefetch_per_visit:
                break
    if not missing:
        return 0
    try:
        connection = _redis()
        queue = redisqueue.queue(redisqueue.MAINTENANCE_QUEUE, connection)
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
        _enqueue_chunks(queue, db, claimed, "prefetch")
        if claimed:
            workers.ensure_workers(redisqueue.MAINTENANCE_QUEUE, connection, redisqueue.redis_url())
        return len(claimed)
    except (redis.RedisError, OSError, redisqueue.Unavailable) as exc:
        log.debug(f"Node prefetch: not queued ({exc}).")
        return 0


REBUILD_MAX = 2048
"""int: Most nodes one change queues for a rebuild (coarsest first); the rest are built by the next visit that wants them."""


def rebuild(db, generation, keys):
    """
    Queues the building of nodes a change deleted from the stored generation (`tilecache.current_stamp`): the changed
    nodes that were stored, coarsest first (the changed block chains are the nodes over a changed sector, so the root
    nodes come first), so the next visit finds them built again (PERF.133, MAP.260). The rules are the prefetch's except
    the per-visit and per-minute counts and the queue depth (a change touches a few hundred nodes at most, and it is
    the data the visitors are looking at): nothing while the prefetch is off (`prefetch_per_minute` 0), during an
    update, or while the cache is full; a node already claimed is not queued twice.

    Args:
        db (str): The database the site shows.
        generation (str): The generation the nodes belong to.
        keys (list[str]): The node keys that were stored and were deleted.

    Returns:
        int: How many nodes were queued.
    """
    settings = get_settings().tile_cache
    if settings.prefetch_per_minute <= 0 or not keys or updateguard.is_updating():
        return 0
    root = tilecache.cache_dir()
    if root is None or is_full(root):
        return 0
    ordered = sorted(keys, key=lambda key: (-parse_node_key(key).m, key))[:REBUILD_MAX]
    try:
        connection = _redis()
        queue = redisqueue.queue(redisqueue.MAINTENANCE_QUEUE, connection)
        pipe = connection.pipeline()
        for key in ordered:
            pipe.set(_claim_name(db, generation, key), "1", nx=True, ex=CLAIM_SECONDS)
        claimed = [key for key, won in zip(ordered, pipe.execute()) if won]
        _enqueue_chunks(queue, db, claimed, "rebuild")
        if claimed:
            workers.ensure_workers(redisqueue.MAINTENANCE_QUEUE, connection, redisqueue.redis_url())
        return len(claimed)
    except (redis.RedisError, OSError, redisqueue.Unavailable) as exc:
        log.debug(f"Node rebuild: not queued ({exc}).")
        return 0


def build_nodes(db, keys):
    """
    The queued job: builds the nodes `keys` into the disk cache the way a visitor's request does
    (`tilecache.fetch_nodes`), so what is stored is what a visitor is served. A node stored since it was queued is not
    built again.

    Returns:
        dict: `nodes` asked and `built` (those that were not stored).
    """
    from planetgen.web import warmup
    app = warmup.build_app()
    with app.test_request_context("/"):
        result = tilecache.fetch_nodes(db, list(keys), ahead=True)
    return {"nodes": len(keys), "built": len(keys) - result["cached"]}
