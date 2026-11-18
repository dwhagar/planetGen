#!/usr/bin/env python3
# legacy/bump_database_to_8_0.py

"""
Brings a planetGen content database made by a release before 8.0 up to the
oldest structure 8.0 and above accept (DB.126, Boss's hard cut of 2026-11-12:
no database structure before 8.0 is accepted).

Run it by hand, on purpose, once per database; `update.sh`, the tests and the
site never run it. It does nothing else: it runs the Alembic revisions
v62 to v130 in `legacy/db_revisions/` (copies of the steps `planetgen/db/migrations/`
held before the re-baseline, DB.98) and writes a `schema_migrations` row for
each. From v130 on, `update.sh` (or `python3 -m planetgen.cli.migrate`)
carries the database to the current structure as usual.

A database that is at v130 or later has nothing to do and is left alone. One
below v61 (the last structure the old `schema.sql` steps made, before Alembic) or
with no recorded version cannot be converted here: restore it from a backup or
make a new database.

Back the database up first (`mysqldump`); a step that stops half way is run again
by running this script again (each step is written to be repeated, DB.18), but a
backup is the only way back.

Usage:
    python3 legacy/bump_database_to_8_0.py --status [--mysql-database NAME]
    python3 legacy/bump_database_to_8_0.py --yes    [--mysql-host H --mysql-port P
                                                     --mysql-user U --mysql-password W
                                                     --mysql-database NAME]

    --status  prints the database's version and the steps that would run; changes nothing.
    --yes     runs the steps without asking.

The revisions import a few helpers from `planetgen` (`planetgen.db.migrations.helpers`,
`planetgen.galaxy.cellkey`, `spherical`, `skeleton`, `density`, `planetgen.db.block_stats`,
`nebula_blocks`, `planetgen.tuning`). They were written against the code of their day;
if one of them has changed shape since, run this script from a checkout of the last
release before the re-baseline (git history, commit b3f60a487^) instead. No test covers
this script (Boss: it only has to exist).
"""

import argparse
import os
import sys

import sqlalchemy
from alembic import command
from alembic.config import Config
from alembic.runtime.migration import MigrationContext
from alembic.script import ScriptDirectory

from planetgen.db import login, store

HERE = os.path.dirname(os.path.abspath(__file__))
REVISIONS_DIR = os.path.join(HERE, "db_revisions")
"""str: The Alembic directory (`env.py`, `versions/`) of the legacy steps v61 to v130."""

OLDEST_VERSION = 61
"""int: The oldest structure these steps start from (the Alembic baseline of the old history)."""

TARGET_VERSION = 130
"""int: The structure 8.0 and above accept as their oldest (`alembic_runner.BASELINE_VERSION`)."""


def _engine(config):
    return sqlalchemy.create_engine(
        "mysql+pymysql://",
        creator=lambda: login.connect(
            host=config.host, port=config.port, user=config.user, password=config.password,
            database=config.database, charset="utf8mb4", autocommit=False,
            init_command="SET time_zone = '+00:00', SESSION lock_wait_timeout = 3600"),
        poolclass=sqlalchemy.pool.NullPool)


def stored_version(connection):
    """The newest `schema_migrations` version, or `None` when the table is missing or empty."""
    exists = connection.execute(sqlalchemy.text(
        "SELECT COUNT(*) FROM information_schema.tables WHERE table_schema = DATABASE() AND table_name = 'schema_migrations'")).scalar()
    if not exists:
        return None
    return connection.execute(sqlalchemy.text("SELECT MAX(version) FROM schema_migrations")).scalar()


def steps(from_version):
    """The `(version, revision id)` pairs from just after `from_version` up to `TARGET_VERSION`."""
    script = ScriptDirectory.from_config(_config())
    ordered = list(reversed(list(script.walk_revisions())))
    return [(int(rev.revision), rev.revision) for rev in ordered if from_version < int(rev.revision) <= TARGET_VERSION]


def _config(connection=None):
    cfg = Config()
    cfg.set_main_option("script_location", REVISIONS_DIR)
    if connection is not None:
        cfg.attributes["connection"] = connection
    return cfg


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    store.add_mysql_connection_args(parser)
    parser.add_argument("--status", action="store_true", help="print the version and the steps to run; change nothing")
    parser.add_argument("--yes", action="store_true", help="run the steps without asking")
    args = parser.parse_args(argv)
    config = store.mysql_config_from_args(args)
    engine = _engine(config)
    try:
        with engine.connect() as connection:
            version = stored_version(connection)
            if version is None:
                print(f"{config.database}: no schema version recorded; this script cannot convert it. Restore a backup or make a new database.")
                return 2
            if version >= TARGET_VERSION:
                print(f"{config.database}: v{version} is already v{TARGET_VERSION} or later; nothing to do (run update.sh for the rest).")
                return 0
            if version < OLDEST_VERSION:
                print(f"{config.database}: v{version} is older than v{OLDEST_VERSION}, the oldest structure these steps start from. Restore a backup or make a new database.")
                return 2
            todo = steps(version)
            print(f"{config.database}: v{version} -> v{TARGET_VERSION} in {len(todo)} steps")
            if args.status:
                return 0
            if not args.yes and input("Have you backed this database up, and run it now? [y/N] ").strip().lower() != "y":
                print("Stopped; nothing changed.")
                return 1
            cfg = _config(connection)
            current = MigrationContext.configure(connection).get_current_revision()
            if current != f"{version:04d}":
                command.stamp(cfg, f"{version:04d}", purge=current is not None)
                connection.commit()
            for number, (target, revision) in enumerate(todo, start=1):
                print(f"[{number}/{len(todo)}] v{target} ({revision})", flush=True)
                command.upgrade(cfg, revision)
                connection.commit()
                connection.execute(sqlalchemy.text("INSERT INTO schema_migrations (version) VALUES (:v)"), {"v": target})
                connection.commit()
            print(f"{config.database}: now v{TARGET_VERSION}. Run update.sh (or python3 -m planetgen.cli.migrate) to bring it to the current structure.")
            return 0
    finally:
        engine.dispose()


if __name__ == "__main__":
    sys.exit(main())
