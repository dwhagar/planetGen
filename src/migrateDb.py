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

from stellarObjects import adminAuth
from stellarObjects._db import (
    SCHEMA_VERSION,
    add_mysql_connection_args,
    control_mysql_config,
    migrate_database,
    mysql_config_from_args,
)


def main():
    parser = argparse.ArgumentParser(
        description="Bring the configured MySQL database's schema_migrations bookkeeping up to date, "
                     "and the control schema (admin logins/sessions/API keys -- see control_schema.sql) "
                     "alongside it.",
    )
    add_mysql_connection_args(parser)
    args = parser.parse_args()

    config = mysql_config_from_args(args)

    try:
        version = migrate_database(config)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)

    if version == SCHEMA_VERSION:
        print(f"Database is at schema v{SCHEMA_VERSION} (current).")
    else:
        print(f"Database is at schema v{version}, expected v{SCHEMA_VERSION} -- no migration path available yet.")
        sys.exit(1)

    # The control schema (admin_users/admin_sessions/admin_api_keys/
    # admin_audit_log) is a separate schema from the content database just
    # migrated above (see control_schema.sql's header comment) -- ensured/
    # seeded here too so a fresh deploy's default admin/password login
    # exists without a separate manual step. Reuses this same account's
    # host/user/password (the full-access account this script already
    # runs as, same as install.sh's existing step 2) against the control
    # schema's own name instead of the content database's.
    try:
        adminAuth.bootstrap_control_schema(control_mysql_config(config))
    except Exception as exc:
        print(f"error: could not set up the control schema ({exc}).", file=sys.stderr)
        sys.exit(1)
    print("Control schema (admin logins) is up to date.")


if __name__ == "__main__":
    main()
