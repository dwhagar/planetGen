#!/usr/bin/env python3
# src/migrateDb.py

"""
Migrates every planetGen SQLite database in a directory to the current
schema (`stellarObjects/schema.sql`'s `PRAGMA user_version`), backing
each one up first -- a no-op for a database that's already current. See
`stellarObjects/_db.py`'s `migrate_database` for how a single database is
converted and backed up.

Run automatically by `install.sh` (and so by `update.sh`, which calls it)
on every deploy, so a database generated under an older schema keeps
working after a `git pull` brings in a newer one. Also runnable directly
for a one-off migration outside of a deployment.

This file lives alongside `stellarObjects/` under `src/`, so Python's own
sys.path[0] (the running script's directory) already makes
`stellarObjects` importable -- no sys.path shim needed.

Usage:
    python3 src/migrateDb.py [db_dir]

    db_dir defaults to `db/` at the repo root (two directories up from
    this file's location under src/), not alongside this script.
"""

import argparse
import glob
import os
import sys

from stellarObjects._db import BACKUP_MARKER, SCHEMA_VERSION, UnsupportedSchemaVersionError, migrate_database

DEFAULT_DB_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "db")


def main():
    parser = argparse.ArgumentParser(description="Migrate every *.db file in a directory to the current schema.")
    parser.add_argument(
        "db_dir", nargs="?", default=DEFAULT_DB_DIR,
        help=f"Directory containing *.db files (default: {DEFAULT_DB_DIR}).",
    )
    args = parser.parse_args()

    if not os.path.isdir(args.db_dir):
        print(f"No database directory at {args.db_dir} -- nothing to migrate.")
        return

    # Backups made by a previous run of this same migration are gzip
    # files (`.db.gz`, from `migrate_database`) so they wouldn't match
    # this `*.db` glob anyway -- but they're also excluded by name here,
    # explicitly, so a backup is never handed back into
    # `migrate_database` (which would re-migrate it and stamp out a
    # nested backup) even if the naming scheme changes again later.
    paths = [
        path for path in sorted(glob.glob(os.path.join(args.db_dir, "*.db")))
        if BACKUP_MARKER not in os.path.basename(path)
    ]
    if not paths:
        print(f"No *.db files in {args.db_dir} -- nothing to migrate.")
        return

    for path in paths:
        try:
            backup_path = migrate_database(path)
        except UnsupportedSchemaVersionError as exc:
            print(f"error: {exc}", file=sys.stderr)
            sys.exit(1)
        if backup_path:
            print(f"Migrated {path} to schema v{SCHEMA_VERSION} (backup: {backup_path})")
        else:
            print(f"{path}: already at schema v{SCHEMA_VERSION}, skipped")


if __name__ == "__main__":
    main()
