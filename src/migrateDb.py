#!/usr/bin/env python3
# src/migrateDb.py

"""
Brings the configured planetGen MySQL database's `schema_migrations`
bookkeeping up to the current schema (`stellarObjects/schema.sql`),
applying any migration step in between -- a no-op for a database that's
already current. See `stellarObjects/_db.py`'s `migrate_database` for how
a single database's version is checked/advanced.

Run automatically by `install.sh` (and so by `update.sh`, which calls it)
on every deploy, so a database created under an older schema keeps
working after a `git pull` brings in a newer one. Also runnable directly
for a one-off check/migration outside of a deployment.

This file lives alongside `stellarObjects/` under `src/`, so Python's own
sys.path[0] (the running script's directory) already makes
`stellarObjects` importable -- no sys.path shim needed.

Usage:
    python3 src/migrateDb.py [--mysql-host HOST] [--mysql-port PORT]
                             [--mysql-user USER] [--mysql-password PASSWORD]
                             [--mysql-database DATABASE]

    Every flag defaults to the same $PLANETGEN_MYSQL_* environment
    variable every other entry point in this project reads (see
    `stellarObjects._db.MySQLConfig`) -- unlike the pre-MySQL-port version
    of this script, there is exactly one database to migrate (a MySQL
    server, not a directory of `*.db` files), so this needs a connection
    to point at rather than a directory to scan.
"""

import argparse
import sys

from stellarObjects._db import SCHEMA_VERSION, add_mysql_connection_args, migrate_database, mysql_config_from_args


def main():
    parser = argparse.ArgumentParser(
        description="Bring the configured MySQL database's schema_migrations bookkeeping up to date.",
    )
    add_mysql_connection_args(parser)
    args = parser.parse_args()

    try:
        version = migrate_database(mysql_config_from_args(args))
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)

    if version == SCHEMA_VERSION:
        print(f"Database is at schema v{SCHEMA_VERSION} (current).")
    else:
        print(f"Database is at schema v{version}, expected v{SCHEMA_VERSION} -- no migration path available yet.")
        sys.exit(1)


if __name__ == "__main__":
    main()
