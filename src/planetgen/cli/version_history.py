#!/usr/bin/env python3
# planetgen.cli.version_history

"""
Records the version key every galaxy database was updated under (OPS.13).

`update.sh` runs it after the code and databases are current. For each
galaxy database that has a seed it adds a history row when the version key,
the generator epoch, the generator battery's digest, the lock file or the
environment differs from the newest row, and otherwise only sets that row's
last-seen time (`planetgen.galaxy.version_history`), and says whether the generated output changed (the battery's digest
against the newest row's, naming the first case that differs, OPS.15); `planetgen versions`
lists what was recorded. The battery (OPS.57) is worked out once, about five
seconds, however many galaxies there are.

Usage:
    python3 -m planetgen.cli.version_history [--mysql-host HOST] ...
"""

import argparse
import sys

from planetgen.db import store
from planetgen.galaxy import version_history
from planetgen.generation import generator_battery
from planetgen.util import log


def _galaxy_databases(config):
    """Every galaxy database on the server, with its seed when it has one."""
    found = []
    for entry in store.list_databases(config):
        conn = store.get_connection(store.MySQLConfig(
            host=config.host, port=config.port, user=config.user, password=config.password,
            database=entry["name"]), ensure_schema=False)
        try:
            try:
                seed = store.get_galaxy_seed(conn)
            except Exception:  # noqa: BLE001 -- not a planetGen database, or older than v51
                seed = None
        finally:
            conn.close()
        found.append((entry["name"], seed))
    return found


def record_all(config, progress=None):
    """Adds a history row for every planned galaxy; returns `[(database, Recorded)]` (OPS.15: the key and whether the
    generated output changed). `progress(done, total)` is called as the generator battery's cases are worked out."""
    control = store.get_control_connection(store.control_mysql_config(config), ensure_schema=False)
    try:
        planned = [(name, seed) for name, seed in _galaxy_databases(config) if seed is not None]
        battery = generator_battery.running_digest(progress) if planned else None
        return [(name, version_history.record(control, name, seed, battery=battery)) for name, seed in planned]
    finally:
        control.close()


def _record_with_bar(config):
    """`record_all` with a bar on a terminal (the cases done and the time left), a line per case without one."""
    if not sys.stdout.isatty():
        return record_all(config, progress=lambda done, total: print(f"  {done} of {total} generator test cases done") if done else None)
    from rich.progress import BarColumn, MofNCompleteColumn, Progress, TextColumn, TimeRemainingColumn
    with Progress(TextColumn("  generator test cases"), BarColumn(), MofNCompleteColumn(), TimeRemainingColumn(), transient=True) as bar:
        task = bar.add_task("cases", total=None)

        def advance(done, total):
            bar.update(task, completed=done, total=total)
        return record_all(config, progress=advance)


def format_history(config):
    """The listing `planetgen versions` prints."""
    control = store.get_control_connection(store.control_mysql_config(config), ensure_schema=False)
    lines = []
    try:
        for name, seed in _galaxy_databases(config):
            if seed is None:
                continue
            lines.append(f"{name} (seed {seed.hex().upper()})")
            rows = version_history.history(control, name)
            if not rows:
                lines.append("  no updates recorded yet")
            for position, row in enumerate(rows):
                lock = (row["requirements_sha256"] or "no lock file")[:12]
                epoch = "epoch ?" if row["generator_epoch"] is None else f"epoch {row['generator_epoch']}"
                battery = (row["battery_digest"] or "no battery")[:12]
                lines.append(f"  {row['recorded_at']:%Y-%m-%d %H:%M}  {row['version_key']}  "
                             f"PlanetGen {row['planetgen_version']}  {epoch}  battery {battery}  lock {lock}  "
                             f"last seen {row['last_seen_at']:%Y-%m-%d %H:%M}")
                changed = version_history.changes(row, rows[position + 1] if position + 1 < len(rows) else None)
                if changed:
                    lines.append("    changed: " + ", ".join(changed))
    finally:
        control.close()
    return "\n".join(lines) if lines else "No planned galaxy databases."


def main():
    parser = argparse.ArgumentParser(description="Record the version key each galaxy was updated under.")
    store.add_mysql_connection_args(parser)
    log.add_debug_arguments(parser)
    args = parser.parse_args()
    log.apply_debug_arguments(args, parser)
    config = store.mysql_config_from_args(args)
    try:
        recorded = _record_with_bar(config)
    except Exception as exc:  # noqa: BLE001
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
    for name, done in recorded:
        print(f"Recorded version key {done.key} for {name}. {done.output}")
    if not recorded:
        print("No planned galaxy database to record a version key for.")


if __name__ == "__main__":
    main()
