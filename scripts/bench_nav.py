#!/usr/bin/env python3
"""
Times `query.nav_between` on a large synthetic galaxy database (NAV.10).

Builds a throwaway database of `--sectors` generated sectors on a flat grid
(4 pc apart), `--per-sector` bare star systems in each, written with bulk SQL
(no generation), then times a route between opposite corners with the
corridor search and, for comparison, what rebuilding the graph over every
placed system would cost. The database is dropped at the end.

Usage (a MySQL account that can create databases, via the PLANETGEN_MYSQL_*
variables):

    python scripts/bench_nav.py [--sectors 100000] [--per-sector 4] [--keep]
"""

import argparse
import math
import random
import sys
import time

sys.path.insert(0, "src")

from planetgen import tuning  # noqa: E402
from planetgen.db import query, store  # noqa: E402
from planetgen.galaxy.nav_graph import build_route_graph  # noqa: E402
from planetgen.physics.units import pc_to_ly  # noqa: E402

BATCH = 5000


def build(conn, sectors, per_sector, edge_pc):
    side = math.ceil(math.sqrt(sectors))
    rng = random.Random(1)
    config_id = conn.execute("INSERT INTO system_configs (markdown) VALUES (0)").lastrowid
    count = 0
    rows = []
    for index in range(sectors):
        i, j = divmod(index, side)
        rows.append((f"S{index}", edge_pc * 1000.0, i * edge_pc, j * edge_pc, 0.0,
                     math.hypot(i * edge_pc, j * edge_pc), i, 0, j))
        if len(rows) == BATCH or index == sectors - 1:
            conn.executemany(
                "INSERT INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc,"
                " ring_index, layer_index, ring_slot_index) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", rows)
            conn.commit()
            count += len(rows)
            rows = []
            print(f"  {count:,} sectors", end="\r", flush=True)
    print()
    ids = [row["id"] for row in conn.execute("SELECT id FROM sectors ORDER BY id").fetchall()]
    half = edge_pc * 500.0
    rows = []
    total = 0
    for sector_id in ids:
        for _ in range(per_sector):
            rows.append((sector_id, config_id, f"Star {total}", rng.uniform(-half, half), rng.uniform(-half, half),
                         rng.uniform(-half, half)))
            total += 1
        if len(rows) >= BATCH:
            _flush(conn, rows)
            rows = []
            print(f"  {total:,} systems", end="\r", flush=True)
    _flush(conn, rows)
    print()
    return total


def _flush(conn, rows):
    if rows:
        conn.executemany(
            "INSERT INTO star_systems (sector_id, system_config_id, name, position_x_mpc, position_y_mpc,"
            " position_z_mpc) VALUES (?, ?, ?, ?, ?, ?)", rows)
        conn.commit()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sectors", type=int, default=100000)
    parser.add_argument("--per-sector", type=int, default=4)
    parser.add_argument("--keep", action="store_true", help="Keep the database afterwards.")
    args = parser.parse_args()

    base = store.DEFAULT_MYSQL_CONFIG
    name = "planetgen_bench_nav"
    admin = store.get_connection(store.MySQLConfig(host=base.host, port=base.port, user=base.user,
                                                   password=base.password, database=""), ensure_schema=False)
    admin.execute(f"DROP DATABASE IF EXISTS {name}")
    admin.execute(f"CREATE DATABASE {name}")
    admin.commit()
    config = store.MySQLConfig(host=base.host, port=base.port, user=base.user, password=base.password, database=name)
    conn = store.get_connection(config)
    try:
        edge_pc = float(tuning.DEFAULT_SECTOR_EDGE_PC)
        started = time.monotonic()
        systems = build(conn, args.sectors, args.per_sector, edge_pc)
        print(f"Built {args.sectors:,} sectors and {systems:,} systems in {time.monotonic() - started:.0f} s.")
        first = conn.execute("SELECT id FROM star_systems ORDER BY id LIMIT 1").fetchone()["id"]
        last = conn.execute("SELECT id FROM star_systems ORDER BY id DESC LIMIT 1").fetchone()["id"]

        started = time.monotonic()
        result = query.nav_between(conn, first, last)
        elapsed = time.monotonic() - started
        route = result["route"]
        print(f"Route over {result['direct'].distance_ly:,.0f} ly: {len(route['path'])} stops, "
              f"{route['distance_ly']:,.0f} ly, corridor search {elapsed:.2f} s.")

        started = time.monotonic()
        rows = conn.execute(
            "SELECT ss.id, ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc, sec.center_x_pc,"
            " sec.center_y_pc, sec.center_z_pc FROM star_systems ss JOIN sectors sec ON sec.id = ss.sector_id"
        ).fetchall()
        positions = {row["id"]: (pc_to_ly(row["center_x_pc"] + row["position_x_mpc"] / 1000.0),
                                 pc_to_ly(row["center_y_pc"] + row["position_y_mpc"] / 1000.0),
                                 pc_to_ly(row["center_z_pc"] + row["position_z_mpc"] / 1000.0)) for row in rows}
        loaded = time.monotonic() - started
        started = time.monotonic()
        build_route_graph(positions, tuning.NAV_ADJACENCY_K, tuning.NAV_ISLAND_LINKS)
        print(f"For comparison, loading every system takes {loaded:.1f} s and building the graph over all "
              f"{len(positions):,} another {time.monotonic() - started:.1f} s.")
    finally:
        conn.close()
        if not args.keep:
            admin.execute(f"DROP DATABASE IF EXISTS {name}")
            admin.commit()


if __name__ == "__main__":
    main()
