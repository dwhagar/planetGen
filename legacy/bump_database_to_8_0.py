#!/usr/bin/env python3
"""
Brings a planetGen database from before 8.0 up to the oldest structure 8.0 accepts.

Run only by hand, once, on a database you want to keep (back it up first). Nothing in `install.sh`, `update.sh`, the migrations or the program calls this. Boss's hard cut (2026-11-12):
no database structure before 8.0 is accepted, so a database below v130 is refused (`store.SchemaTooOldError`) and the program never converts it. This script is the one way to keep one:
it runs the old Alembic revisions of `legacy/db_revisions/versions/` (0062 to 0130, as they were before DB.98 replaced them with an empty baseline) on the database, one after another, writing each
version it reaches to `schema_migrations`, and stops at v130. Then run `update.sh`: the normal migrations take it from v130 to the current schema.

It does that and nothing else. It does not generate, delete or change any galaxy data beyond what those revisions did. Run it from the root of a checkout of this version::

    python3 legacy/bump_database_to_8_0.py --database planetgen --host 127.0.0.1 --user planetgen

(the password is asked for, or `PLANETGEN_MYSQL_PASSWORD`). `--yes` skips the question. A database already at v130 or newer, or one with no version recorded, is left alone.
"""

import argparse
import getpass
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

FIRST_VERSION = 61
"""int: The oldest version these revisions start from (revision 0061 was the baseline of the history before DB.98)."""

LAST_VERSION = 130
"""int: The version the old history ends at: `alembic_runner.BASELINE_VERSION`, the oldest structure 8.0 upgrades from."""


def _stored_version(connection):
    import sqlalchemy as sa
    tables = {row[0] for row in connection.execute(sa.text("SHOW TABLES"))}
    if "schema_migrations" not in tables:
        return None
    return connection.execute(sa.text("SELECT MAX(version) FROM schema_migrations")).scalar()


def bump(config, assume_yes=False, say=print):
    """
    Runs the old revisions on `config`'s database up to v130.

    Returns:
        int: The version the database is at after the run (unchanged when nothing was done).
    """
    import sqlalchemy as sa
    from alembic import command
    from alembic.config import Config
    from alembic.runtime.migration import MigrationContext
    from planetgen.db import alembic_runner

    engine = alembic_runner._engine(config)
    try:
        with engine.connect() as connection:
            version = _stored_version(connection)
            if version is None:
                say(f"Database '{config.database}' has no schema version recorded: nothing to bump. Restore it from a backup, or drop it and run update.sh.")
                return None
            if version >= LAST_VERSION:
                say(f"Database '{config.database}' is at v{version}, which this planetGen accepts: nothing to do. Run update.sh.")
                return version
            if version < FIRST_VERSION:
                say(f"Database '{config.database}' is at v{version}, older than v{FIRST_VERSION}, the oldest this script starts from: nothing was done.")
                return version
            if not assume_yes:
                answer = input(f"Bring database '{config.database}' from v{version} to v{LAST_VERSION}? Back it up first. Type 'bump' to go on: ")
                if answer.strip() != "bump":
                    say("Nothing was done.")
                    return version
            cfg = Config()
            cfg.set_main_option("script_location", os.path.join(HERE, "db_revisions"))
            cfg.attributes["connection"] = connection
            current = MigrationContext.configure(connection).get_current_revision()
            if current != alembic_runner.revision_id(version):
                command.stamp(cfg, alembic_runner.revision_id(version), purge=current is not None)
                connection.commit()
            for target in range(version + 1, LAST_VERSION + 1):
                say(f"v{version} to v{target} ...")
                command.upgrade(cfg, alembic_runner.revision_id(target))
                connection.commit()
                connection.execute(sa.text("INSERT INTO schema_migrations (version) VALUES (:version)"), {"version": target})
                connection.commit()
                version = target
            say(f"Database '{config.database}' is at v{version}. Now run update.sh.")
            return version
    finally:
        engine.dispose()


def main(argv=None):
    from planetgen.db import store
    defaults = store.DEFAULT_MYSQL_CONFIG
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--host", default=defaults.host)
    parser.add_argument("--port", type=int, default=defaults.port)
    parser.add_argument("--user", default=defaults.user)
    parser.add_argument("--database", default=defaults.database)
    parser.add_argument("--yes", action="store_true", help="do not ask first")
    args = parser.parse_args(argv)
    password = os.environ.get("PLANETGEN_MYSQL_PASSWORD") or defaults.password or getpass.getpass("MySQL password: ")
    config = store.MySQLConfig(host=args.host, port=args.port, user=args.user, password=password, database=args.database)
    bump(config, assume_yes=args.yes)


if __name__ == "__main__":
    main()
