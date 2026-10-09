#!/usr/bin/env python3
# planetgen.cli.version_history

"""
Records the version key every galaxy database was updated under (OPS.13),
or lists what was recorded.

`update.sh` and `update.ps1` run it after the code and databases are
current. With no option it adds one history row per galaxy database that has
a seed (`planetgen.galaxy.version_history`); `--list` prints each database's
rows, newest first, and records nothing. `planetgen versions` is the same
listing.

Usage:
    python3 -m planetgen.cli.version_history [--list] [--mysql-host HOST] ...
"""

import argparse
import sys

from planetgen.db import store
from planetgen.galaxy import version_history


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


def record_all(config):
    """Adds a history row for every planned galaxy; returns `[(database, key)]`."""
    control = store.get_control_connection(store.control_mysql_config(config), ensure_schema=False)
    try:
        return [(name, version_history.record(control, name, seed))
                for name, seed in _galaxy_databases(config) if seed is not None]
    finally:
        control.close()


def format_history(config):
    """The listing `--list` and `planetgen versions` print."""
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
            for row in rows:
                lock = (row["requirements_sha256"] or "no lock file")[:12]
                lines.append(f"  {row['recorded_at']:%Y-%m-%d %H:%M}  {row['version_key']}  "
                             f"PlanetGen {row['planetgen_version']}  lock {lock}")
    finally:
        control.close()
    return "\n".join(lines) if lines else "No planned galaxy databases."


def main():
    parser = argparse.ArgumentParser(description="Record, or list, the version key each galaxy was updated under.")
    store.add_mysql_connection_args(parser)
    parser.add_argument("--list", action="store_true", help="Print the recorded history and record nothing.")
    args = parser.parse_args()
    config = store.mysql_config_from_args(args)
    try:
        if args.list:
            print(format_history(config))
            return
        recorded = record_all(config)
    except Exception as exc:  # noqa: BLE001
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
    for name, key in recorded:
        print(f"Recorded version key {key} for {name}.")
    if not recorded:
        print("No planned galaxy database to record a version key for.")


if __name__ == "__main__":
    main()
