#!/usr/bin/env python3
# planetgen.cli.warm_shape

"""
Works out, during the update, what building a galaxy's density shape needs for every stored cluster: the mean mass of a
system of its age (`generation.system_mass`). A galaxy of tens of thousands of clusters needs tens of thousands of them, a
minute or more, which a web request may not spend (the web daemon is stopped at 60 seconds), so the results are written
to the data cache folder (`system_mass.save_means`) and the web processes read them from there.

`update.sh` runs it after the migration. It never fails the update: a database without clusters, or without the
isochrone data, is skipped with a word.

Usage:
    python3 -m planetgen.cli.warm_shape [--mysql-host HOST] ...
"""

import argparse
import sys
import time

from planetgen.db import store
from planetgen.generation import system_mass
from planetgen.util import log


def cluster_ages(config):
    """The distinct ages (Gy, to the millionth) of the stored clusters and satellites of every galaxy database."""
    ages = set()
    for entry in store.list_databases(config):
        conn = store.get_connection(store.MySQLConfig(
            host=config.host, port=config.port, user=config.user, password=config.password,
            database=entry["name"]), ensure_schema=False)
        try:
            if not conn.execute("SELECT 1 FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE() "
                                "AND TABLE_NAME = 'galaxy_components'").fetchone():
                continue
            ages.update(round(float(row["age_gy"]), 6) for row in conn.execute("SELECT DISTINCT age_gy FROM galaxy_components").fetchall()
                        if row["age_gy"] is not None)
        finally:
            conn.close()
    return sorted(ages)


def _save_with_bar(ages):
    """`system_mass.save_means` with a bar on a terminal (the age count and the time left), a line a minute without one."""
    if not sys.stdout.isatty():
        last = [time.monotonic()]

        def line(done, total):
            if time.monotonic() - last[0] >= 60 or done == total:
                last[0] = time.monotonic()
                print(f"  {done} of {total} ages done")
        return system_mass.save_means(ages, progress=line)
    from rich.progress import BarColumn, MofNCompleteColumn, Progress, TextColumn, TimeRemainingColumn
    with Progress(TextColumn("  cluster ages"), BarColumn(), MofNCompleteColumn(), TimeRemainingColumn()) as bar:
        task = bar.add_task("ages", total=len(ages))
        return system_mass.save_means(ages, progress=lambda done, total: bar.update(task, completed=done))


def main():
    parser = argparse.ArgumentParser(description="Work out the mean system mass of every stored cluster age, once, for the web app.")
    store.add_mysql_connection_args(parser)
    log.add_debug_arguments(parser)
    args = parser.parse_args()
    log.apply_debug_arguments(args, parser)
    config = store.mysql_config_from_args(args)
    try:
        ages = cluster_ages(config)
        if not ages:
            print("No stored clusters: nothing to work out.")
            return
        started = time.perf_counter()
        print(f"Working out the mean system mass of {len(ages)} cluster ages (a big galaxy takes a few minutes):")
        held = _save_with_bar(ages)
        print(f"Mean system mass worked out for {len(ages)} cluster ages in {time.perf_counter() - started:.1f} s; "
              f"{held} ages are stored in {system_mass.STORE_FILE}.")
    except Exception as exc:  # noqa: BLE001 -- the update carries on
        print(f"warning: the mean system masses could not be worked out ({exc}); the web app works them out when asked.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
