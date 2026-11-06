#!/usr/bin/env python3
"""
Measures what a full build of the Galaxy Map's tile cache costs (PERF.233): how many nodes the whole galaxy has at each zoom
level, how many bytes one takes on disk, how many seconds one takes to build, and from those the disk and the hours of a build
of all of them on this machine.

A "tile" of the single-space Galaxy Map is a sector-block node (`src/planetgen/db/nodes.py`): a block of the ladder (`galaxy/drill.py`, 177,147
shells down to 3) that holds a generated sector or a pre-placed object, and a node for each generated sector (the ladder's
last step). The blocks are the stored rows of `galaxy_block_stats` and `galaxy_block_lights`; the sectors are `sectors`.

It only reads the database. It builds a random sample of each level through `tilecache.fetch_nodes`, the call the map's own
request and `planetgen warm` make, into a private folder it makes and removes (`--scratch`), never the site's cache folder, and
its numbers come from that sample: the count of a level is exact, the bytes and seconds are the sample's mean. The sample is built
once on one process for the seconds and the bytes of each node, then again on `--workers` processes in the chunks `planetgen warm`
builds (`warmup.NODE_CHUNK`) for the throughput of a parallel build.

    python scripts/measure_tile_cache.py --database planetgen --sample 40 --workers 4
"""

import argparse
import json
import math
import multiprocessing
import os
import random
import shutil
import statistics
import sys
import tempfile
import time

LEVEL_ONE = 1
"""int: The level of a sector's own node."""


def _parse(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--database", required=True, help="the galaxy's database (the site's `?db=` value)")
    parser.add_argument("--sample", type=int, default=40, help="nodes built per level for the timings (default 40)")
    parser.add_argument("--workers", type=int, default=4, help="processes of the parallel pass (default 4)")
    parser.add_argument("--seed", type=int, default=1, help="seed of the sample (default 1)")
    parser.add_argument("--scratch", default=None, help="a folder to build the sample in (default: a new one in the temp folder)")
    parser.add_argument("--json", default=None, help="also write the numbers to this file")
    parser.add_argument("--keep", action="store_true", help="leave the scratch folder in place")
    return parser.parse_args(argv)


def level_counts(conn):
    """`{m: populated blocks}` for the levels above a sector (the stored block rows) and `{1: generated sectors}`."""
    rows = conn.execute(
        "SELECT m, COUNT(*) AS n FROM (SELECT m, shell_index, band_index, wedge_index FROM galaxy_block_stats WHERE m >= 3 UNION "
        "SELECT m, shell_index, band_index, wedge_index FROM galaxy_block_lights WHERE m >= 3) u GROUP BY m").fetchall()
    counts = {int(row["m"]): int(row["n"]) for row in rows}
    counts[LEVEL_ONE] = int(conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"])
    return counts


def sample_keys(conn, counts, per_level, rng):
    """`{m: [node keys]}`: up to `per_level` populated nodes of each level, chosen at random."""
    from planetgen.galaxy.drill import DrillBlock, format_drill_key
    found = {}
    for m, total in counts.items():
        if not total:
            continue
        if m == LEVEL_ONE:
            rows = [conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors ORDER BY id LIMIT 1 OFFSET ?",
                                 (offset,)).fetchone() for offset in sorted(rng.sample(range(total), min(per_level, total)))]
        else:
            rows = conn.execute(
                "SELECT shell_index, band_index, wedge_index FROM (SELECT shell_index, band_index, wedge_index FROM galaxy_block_stats "
                "WHERE m = ? UNION SELECT shell_index, band_index, wedge_index FROM galaxy_block_lights WHERE m = ?) u "
                "ORDER BY RAND(?) LIMIT ?", (m, m, rng.randrange(1 << 30), per_level)).fetchall()
        found[m] = [format_drill_key(DrillBlock(m, r["shell_index"], r["band_index"], r["wedge_index"])) for r in rows]
    return found


def _node_files(root):
    """`{file name: size}` of the node files under `root`."""
    found = {}
    for dirpath, _dirs, names in os.walk(root):
        for name in names:
            if name.startswith("n") and name.endswith(".bin"):
                found[name] = os.path.getsize(os.path.join(dirpath, name))
    return found


def _build(database, keys):
    """Builds `keys` through `tilecache.fetch_nodes`, one request per call, the way the warm-up does; returns the seconds."""
    from planetgen.web import warmup
    from planetgen.web.lib import tilecache
    app = warmup.build_app()
    started = time.monotonic()
    for start in range(0, len(keys), warmup.NODE_CHUNK):
        with app.test_request_context("/"):
            tilecache.fetch_nodes(database, keys[start:start + warmup.NODE_CHUNK], ahead=True)
    return time.monotonic() - started


def _build_chunk(job):
    return _build(*job)


def time_each(database, keys_by_level, root):
    """`{m: [(seconds, bytes), ...]}`: every sampled node built alone, with the size of the file it left."""
    from planetgen.web import warmup
    from planetgen.web.lib import tilecache
    app = warmup.build_app()
    with app.test_request_context("/"):
        tilecache.fetch_nodes(database, [], ahead=True)         # the first call reads the stamp and makes the generation folder
    measured = {}
    for m, keys in keys_by_level.items():
        rows = []
        for key in keys:
            before = _node_files(root)
            started = time.monotonic()
            with app.test_request_context("/"):
                tilecache.fetch_nodes(database, [key], ahead=True)
            seconds = time.monotonic() - started
            after = _node_files(root)
            size = sum(after[name] for name in after if name not in before) or after.get(
                "n" + key.replace(".", "_") + ".bin", 0)
            rows.append((seconds, size))
        measured[m] = rows
    return measured


def parallel_pass(database, keys_by_level, workers, root):
    """Builds every sampled node once more on `workers` processes in chunks of `warmup.NODE_CHUNK`; returns (nodes, seconds)."""
    from planetgen.web import warmup
    shutil.rmtree(root, ignore_errors=True)
    os.makedirs(root, exist_ok=True)
    keys = [key for level in keys_by_level.values() for key in level]
    chunk = max(1, min(warmup.NODE_CHUNK, math.ceil(len(keys) / (workers * 4))))
    jobs = [(database, keys[start:start + chunk]) for start in range(0, len(keys), chunk)]
    started = time.monotonic()
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        pool.map(_build_chunk, jobs, chunksize=1)
    return len(keys), time.monotonic() - started


def summarize(counts, measured, nodes_parallel, seconds_parallel, workers, budget_mb):
    """The report: per level the count, the sample's mean bytes and seconds and the projection; and the totals."""
    levels = []
    total_bytes = total_seconds = 0.0
    for m in sorted(counts, reverse=True):
        rows = measured.get(m, [])
        n = counts[m]
        if rows:
            seconds = [row[0] for row in rows]
            sizes = [row[1] for row in rows]
            mean_seconds, mean_bytes = statistics.fmean(seconds), statistics.fmean(sizes)
            high = sorted(sizes)[min(len(sizes) - 1, int(0.95 * len(sizes)))]
        else:
            mean_seconds = mean_bytes = high = 0.0
        levels.append({"m": m, "nodes": n, "sampled": len(rows), "mean_bytes": mean_bytes, "p95_bytes": high,
                       "mean_seconds": mean_seconds, "total_bytes": mean_bytes * n, "total_seconds": mean_seconds * n})
        total_bytes += mean_bytes * n
        total_seconds += mean_seconds * n
    rate_one = sum(len(rows) for rows in measured.values()) / max(sum(row[0] for rows in measured.values() for row in rows), 1e-9)
    rate_parallel = nodes_parallel / max(seconds_parallel, 1e-9)
    total_nodes = sum(counts.values())
    speedup = rate_parallel / rate_one if rate_one else 0.0
    return {"levels": levels, "nodes": total_nodes, "bytes": total_bytes, "serial_seconds": total_seconds,
            "workers": workers, "parallel_nodes_per_second": rate_parallel, "serial_nodes_per_second": rate_one,
            "speedup": speedup, "parallel_seconds": total_nodes / rate_parallel if rate_parallel else None,
            "budget_mb": budget_mb}


def _hours(seconds):
    return f"{seconds / 3600:.2f} h" if seconds >= 3600 else f"{seconds / 60:.1f} min"


def render(report):
    lines = [f"{'level (shells)':>14} {'nodes':>12} {'sampled':>8} {'bytes/node':>11} {'p95 bytes':>10} {'s/node':>8} "
             f"{'disk':>11} {'serial time':>12}"]
    for level in report["levels"]:
        lines.append(f"{level['m']:>14,} {level['nodes']:>12,} {level['sampled']:>8} {level['mean_bytes']:>11,.0f} "
                     f"{level['p95_bytes']:>10,.0f} {level['mean_seconds']:>8.3f} {level['total_bytes'] / 1048576:>8,.1f} MB "
                     f"{_hours(level['total_seconds']):>12}")
    lines.append(f"{'all':>14} {report['nodes']:>12,} {'':>8} {'':>11} {'':>10} {'':>8} "
                 f"{report['bytes'] / 1048576:>8,.1f} MB {_hours(report['serial_seconds']):>12}")
    lines.append("")
    lines.append(f"One process builds {report['serial_nodes_per_second']:.2f} nodes/s; {report['workers']} processes build "
                 f"{report['parallel_nodes_per_second']:.2f} nodes/s ({report['speedup']:.2f}x).")
    if report["parallel_seconds"]:
        lines.append(f"A full build of {report['nodes']:,} nodes on {report['workers']} processes: about "
                     f"{_hours(report['parallel_seconds'])}, {report['bytes'] / 1048576:,.0f} MB on disk "
                     f"(cache.max_mb is {report['budget_mb']:,.0f}).")
    return "\n".join(lines)


def main(argv=None):
    args = _parse(argv)
    from planetgen.util.settings import get_settings
    budget_mb = float(get_settings().cache.max_mb)       # the site's own, before this run sets its own
    scratch = args.scratch or tempfile.mkdtemp(prefix="planetgen-tile-measure-")
    os.makedirs(scratch, exist_ok=True)
    root = os.path.join(scratch, "tiles")
    os.makedirs(root, exist_ok=True)
    os.environ["PLANETGEN_TILE_CACHE_DIR"] = root
    os.environ["PLANETGEN_CACHE_MAX_MB"] = str(10 ** 7)          # nothing in the sample is pruned or left out
    try:
        from planetgen.db import store
        conn = store.get_connection(store.MySQLConfig(database=args.database), ensure_schema=False)
        try:
            counts = level_counts(conn)
            keys = sample_keys(conn, counts, args.sample, random.Random(args.seed))
        finally:
            conn.close()
        print(f"{sum(counts.values()):,} nodes in {args.database}; building {sum(len(k) for k in keys.values())} of them in {root}",
              flush=True)
        measured = time_each(args.database, keys, root)
        nodes, seconds = parallel_pass(args.database, keys, args.workers, root)
        report = summarize(counts, measured, nodes, seconds, args.workers, budget_mb)
        print(render(report))
        if args.json:
            with open(args.json, "w", encoding="utf-8") as handle:
                json.dump(report, handle, indent=1)
        return 0
    finally:
        if not args.keep:
            shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
