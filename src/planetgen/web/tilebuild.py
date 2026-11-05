# planetgen/web/tilebuild.py

"""
The whole-galaxy build of the Galaxy Map's node cache (PERF.233): every node of every level the database holds, built into a
new cache beside the live one and swapped in when it is complete, so visitors never meet an empty cache and a killed run
leaves the old cache serving. The nightly maintenance timer runs it after the orbital update (`planetgen.cli.tiles`).

A "tile" is a node (`planetgen.db.nodes`, `web/lib/tilecache.py`): one per block of each level of the ladder that holds a
generated sector or a pre-placed object, and one per generated sector. The map's address space is 23.9 billion cells, so
"the entire map" is the populated part: `nodes.all_node_keys`.

How it runs:

1. The nodes are counted per level, and the projected size (`SECTOR_NODE_BYTES`, `BLOCK_NODE_BYTES`, measured in
   reports/perf233-tile-cache-cost-2026-11-04.md) is compared with the free disk: the run refuses to start unless twice
   the projected size is free (`DISK_MARGIN`), saying the figures.
2. The nodes are built into `<cache root>/.build/` (`tilecache.building_into`: the same layout as the live cache, no pruning
   while it runs and the budget not asked), a level at a time, coarsest first, `tile_cache.build_workers` processes at
   once, one progress line per level. A node already stored in the staging folder is not built again, so a run that was
   killed picks up where it was (the staging folder is only removed by a swap or by hand).
3. The staging cache is checked against the database's stamp (a change meanwhile deletes just the changed nodes, which the
   next pass builds), up to `MAX_PASSES` passes, and then swapped in: the live folder is renamed away and the staging one
   takes its name. A node request in the instant between is answered from the API and cached afresh (the cache fails
   open); nothing is ever served from half a cache.
"""

import concurrent.futures
import multiprocessing
import os
import shutil
import time

from planetgen.db import nodes, store
from planetgen.util import log
from planetgen.util.settings import get_settings
from planetgen.web.lib import tilecache

STAGING_NAME = ".build"
"""str: The folder under the cache root the new cache is built in."""

NODE_CHUNK = 64
"""int: Nodes one request builds at a time (the most `/galaxy/nodes` takes)."""

SECTOR_NODE_BYTES = 110_000
"""int: The most a sector's node takes on disk (p95 of reports/perf233-tile-cache-cost-2026-11-04.md: 107 KB)."""

BLOCK_NODE_BYTES = 45_000
"""int: The most a block's node takes on disk (p95: 44.6 KB at the finest block levels)."""

DISK_MARGIN = 2.0
"""float: Free disk the run needs, as a multiple of the projected size: the live cache stays while the new one is built."""

MAX_PASSES = 3
"""int: How many times the build looks again for nodes a change in the database took away while it ran."""


class BuildRefused(Exception):
    """The build did not start, or could not finish; the message says why. The live cache was not touched."""


def levels_of(keys):
    """`{m: count}` of node keys, the biggest blocks first."""
    counts = {}
    for key in keys:
        m = int(key.split(".", 1)[0])
        counts[m] = counts.get(m, 0) + 1
    return dict(sorted(counts.items(), reverse=True))


def projected_bytes(counts):
    """The disk the new cache is expected to take, bytes: `SECTOR_NODE_BYTES` a sector's node, `BLOCK_NODE_BYTES` a block's."""
    return sum(count * (SECTOR_NODE_BYTES if m == 1 else BLOCK_NODE_BYTES) for m, count in counts.items())


def check_disk(root, counts, free_bytes=None):
    """Raises `BuildRefused` unless `DISK_MARGIN` times the projected size is free under `root`; returns the projected bytes."""
    need = projected_bytes(counts)
    os.makedirs(root, mode=0o750, exist_ok=True)
    free = shutil.disk_usage(root).free if free_bytes is None else free_bytes
    if free < DISK_MARGIN * need:
        raise BuildRefused(
            f"the whole map is projected at {need / 1073741824:.1f} GB ({sum(counts.values())} nodes) and the build needs "
            f"{DISK_MARGIN:.0f} times that free while the live cache stays, but {free / 1073741824:.1f} GB is free under {root}")
    return need


def staging_root(root):
    """The folder the new cache is built in."""
    return os.path.join(root, STAGING_NAME)


def workers_setting():
    """How many processes build at once (`tile_cache.build_workers`; 0 is one per core)."""
    configured = get_settings().tile_cache.build_workers
    return configured if configured > 0 else max(1, os.cpu_count() or 1)


_app = None


def _get_app():
    """The site's app the nodes are fetched through, made once per process."""
    global _app
    if _app is None:
        from planetgen.web import warmup
        _app = warmup.build_app()
    return _app


def _init_worker(root):
    """A pool process: builds into the staging cache."""
    tilecache.building_into(root).__enter__()
    _get_app()


def _build_chunk(database, chunk):
    """Builds one chunk of nodes through the map's own request (`tilecache.fetch_nodes`, marked built ahead); returns how many
    were already stored."""
    with _get_app().test_request_context("/"):
        return tilecache.fetch_nodes(database, chunk, ahead=True)["cached"]


def _build_level(database, keys, pool):
    """Builds `keys` (one level); returns `(already stored, built)`."""
    chunks = [keys[start:start + NODE_CHUNK] for start in range(0, len(keys), NODE_CHUNK)]
    cached = 0
    if pool is None:
        for chunk in chunks:
            cached += _build_chunk(database, chunk)
    else:
        futures = [pool.submit(_build_chunk, database, chunk) for chunk in chunks]
        for future in concurrent.futures.as_completed(futures):
            cached += future.result()
    return cached, len(keys) - cached


def missing_keys(database, keys, staging):
    """The `keys` with no stored node in the staging cache's current generation (a stat each, no node is read)."""
    info = tilecache.current_stamp(database, staging)
    if not info.get("generation"):
        return list(keys)
    generation_dir = os.path.join(tilecache._db_dir(staging, database), info["generation"])
    return [key for key in keys if not os.path.isfile(os.path.join(generation_dir, tilecache._node_filename(key)))]


def swap_in(root, database, staging=None):
    """
    Puts the staging cache in the live one's place: the live folder of `database` is renamed away, the staging one takes its
    name and the old one is removed. Returns the old folder's name when there was one.
    """
    staging = staging or staging_root(root)
    live = tilecache._db_dir(root, database)
    new = tilecache._db_dir(staging, database)
    if not os.path.isdir(new):
        raise BuildRefused("nothing was built to swap in")
    old = None
    if os.path.isdir(live):
        old = os.path.join(root, f".old-{os.getpid()}-{int(time.time())}")
        os.rename(live, old)
    os.rename(new, live)
    if old:
        tilecache._remove_tree(old)
    return old


def build_all(database, workers=None, force=False, swap=True, progress=None, keys=None, free_bytes=None):
    """
    Builds every node of the whole galaxy into a new cache and swaps it in (see the module's note).

    Args:
        database (str): The database whose map is built.
        workers (int, optional): Processes to build in; `tile_cache.build_workers` when not given.
        force (bool): Build although the disk check refuses.
        swap (bool): Swap the finished cache in; `False` leaves it in the staging folder.
        progress (callable, optional): Called with each progress line.
        keys (list[str], optional): The nodes to build, instead of every node the database holds (for tests).
        free_bytes (int, optional): The free disk to check against, instead of the disk's own (for tests).

    Returns:
        dict: `nodes`, `built`, `already`, `levels` (`{m: count}`), `seconds`, `swapped`.

    Raises:
        BuildRefused: The cache is off, the disk is too small, or the data changed under every pass.
    """
    say = progress or (lambda line: log.normal(line))
    started = time.monotonic()
    root = tilecache.cache_dir()
    if root is None or tilecache.budget_bytes() == 0:
        raise BuildRefused("the disk cache is off (cache.max_mb is 0 or no cache folder is writable)")
    if keys is None:
        conn = store.get_connection(store.MySQLConfig(database=database), ensure_schema=False)
        try:
            keys = nodes.all_node_keys(conn)
        finally:
            conn.close()
    counts = levels_of(keys)
    if not keys:
        raise BuildRefused("the database holds no generated sector and no placed object: there is no node to build")
    projected = projected_bytes(counts)
    try:
        check_disk(root, counts, free_bytes)
    except BuildRefused as exc:
        if not force:
            raise
        say(f"warning: {exc}; building anyway (--force)")
    workers = workers or workers_setting()
    say(f"whole-map build: {len(keys)} nodes in {len(counts)} levels, about {projected / 1073741824:.2f} GB, {workers} process(es)")
    staging = staging_root(root)
    by_level = {}
    for key in keys:
        by_level.setdefault(int(key.split(".", 1)[0]), []).append(key)
    total_built = total_already = 0
    pool = None
    try:
        if workers > 1:
            pool = concurrent.futures.ProcessPoolExecutor(
                max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
                initializer=_init_worker, initargs=(staging,))
        with tilecache.building_into(staging):
            for m, level_keys in by_level.items():
                level_started = time.monotonic()
                cached, built = _build_level(database, level_keys, pool)
                total_built += built
                total_already += cached
                say(f"level {m}: {len(level_keys)} nodes ({built} built, {cached} already stored) in "
                    f"{time.monotonic() - level_started:.1f} s")
            # A change in the database meanwhile deletes just the nodes it touched (`current_stamp`): they are built again.
            for attempt in range(1, MAX_PASSES + 1):
                lost = missing_keys(database, keys, staging)
                if not lost:
                    break
                if attempt == MAX_PASSES:
                    raise BuildRefused(f"the data kept changing: {len(lost)} nodes are still missing after {MAX_PASSES} passes")
                say(f"{len(lost)} nodes were taken away by a change in the data while the build ran: building them again")
                _cached, built = _build_level(database, lost, pool)
                total_built += built
    finally:
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
    swapped = False
    if swap:
        swap_in(root, database, staging)
        swapped = True
        say("swapped in: the new cache is live")
    seconds = time.monotonic() - started
    say(f"whole-map build done: {len(keys)} nodes ({total_built} built, {total_already} already stored) in {seconds:.0f} s")
    return {"nodes": len(keys), "built": total_built, "already": total_already, "levels": counts,
            "seconds": seconds, "swapped": swapped}
