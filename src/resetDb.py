#!/usr/bin/env python3
# src/resetDb.py

"""
Wipes every generated sector/system/galaxy row out of the configured
planetGen database, leaving it empty and ready for a fresh galaxy --
`TRUNCATE TABLE` (not `DROP DATABASE`/`DROP TABLE`) against every real
table `schema.sql` defines except `schema_migrations`, so:

  - The schema itself (tables, indexes, foreign keys, the `sector_objects`
    view) is left completely intact -- nothing needs re-creating
    afterward, and `schema_migrations`'s own bookkeeping (this database's
    current DDL version) is untouched, so a later `migrateDb.py` run still
    correctly sees "already current" instead of re-running migrations
    against a wiped-but-still-v22 database.
  - `TRUNCATE` (unlike `DELETE FROM`) resets every table's own
    `AUTO_INCREMENT` counter back to 1 -- the new galaxy's first sector/
    system/star/... row gets id 1 again, not some arbitrarily high number
    left over from the wiped one.
  - The table list is discovered live via `SHOW FULL TABLES ... WHERE
    Table_type = 'BASE TABLE'` (excluding `schema_migrations`), not
    hand-maintained here -- a future schema change that adds a new content
    table is picked up automatically, and a view (`sector_objects`) is
    never a `BASE TABLE` so it's excluded for free, no special-casing
    needed.

Deliberately does NOT touch the separate control schema
(`control_schema.sql` -- admin logins/sessions/API keys/audit log,
`stellarObjects.control_mysql_config`'s own database): resetting the
galaxy's content has no reason to also sign every admin out or forget
their accounts, and this script never even connects to that schema.

This is destructive and cannot be undone -- by default, prompts for the
exact database name to be typed back before doing anything (`--yes` skips
the prompt, for scripted/automated use only); `--dry-run` lists what would
be wiped without touching anything.

Usage:
    python3 src/resetDb.py [--mysql-host HOST] [--mysql-port PORT]
                           [--mysql-user USER] [--mysql-password PASSWORD]
                           [--mysql-database DATABASE]
                           [--yes] [--dry-run]

    Every `--mysql-*` flag defaults to the same `$PLANETGEN_MYSQL_*`
    environment variable (or `config.json`'s `mysql` section) every other
    entry point in this project reads -- see
    `stellarObjects._db.MySQLConfig`.
"""

import argparse
import sys

from stellarObjects._db import add_mysql_connection_args, get_connection, mysql_config_from_args

_EXCLUDED_TABLES = {"schema_migrations"}
"""set: Real tables that exist in every fresh database but hold DDL
bookkeeping, not galaxy content -- never truncated. `sector_objects` (a
VIEW) needs no equivalent entry here; `SHOW FULL TABLES ... BASE TABLE`
already excludes it."""


def _content_tables(conn, database):
    """
    Returns the live list of tables to truncate: every `BASE TABLE` in
    `database` except `_EXCLUDED_TABLES` -- see this module's own
    docstring for why this is discovered rather than hand-maintained.
    Reads `information_schema.tables` (a clean, explicitly-named
    `table_name`/`table_type` result) rather than `SHOW FULL TABLES`
    (whose table-name column is instead named after the database itself,
    `Tables_in_<db>`, awkward to read generically).

    Args:
        conn (stellarObjects._db.Connection): An open connection.
        database (str): The database name to list tables from --
            `information_schema` spans every database on the server, so
            this must be explicit rather than implied by the connection.

    Returns:
        list[str]: Table names, in whatever order the query returned them
            (irrelevant here -- `FOREIGN_KEY_CHECKS` is disabled for the
            whole truncate pass below, so no table's own foreign keys can
            block truncating it regardless of order).
    """
    rows = conn.execute(
        "SELECT table_name AS name FROM information_schema.tables "
        "WHERE table_schema = ? AND table_type = 'BASE TABLE'",
        (database,),
    ).fetchall()
    return [row["name"] for row in rows if row["name"] not in _EXCLUDED_TABLES]


def _row_counts(conn, tables):
    """Returns `{table: row_count}` for every name in `tables` -- used to
    show what a `--dry-run`/confirmation prompt is actually about to wipe,
    not just which tables exist."""
    counts = {}
    for table in tables:
        counts[table] = conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
    return counts


def _confirm(database, host, total_rows):
    """
    Prompts for the exact database name before a real (non-`--dry-run`)
    reset -- the standard "type the name back" confirmation for a
    destructive, unrecoverable operation, rather than a bare y/n a
    fat-fingered Enter could pass.

    Returns:
        bool: Whether the operator confirmed.
    """
    print(f"About to permanently wipe {total_rows:,} row(s) of galaxy content from:")
    print(f"  database: {database}")
    print(f"  host:     {host}")
    print()
    print("This cannot be undone. Admin logins/sessions (the separate control")
    print("schema) are not affected.")
    print()
    try:
        typed = input(f"Type the database name ({database}) to confirm, or anything else to cancel: ")
    except EOFError:
        return False
    return typed.strip() == database


def reset_database(config, dry_run=False, assume_yes=False):
    """
    Truncates every content table in `config.database`, after confirming
    unless `assume_yes`.

    Args:
        config (MySQLConfig): Connection parameters -- `config.database`
            is the one database this wipes; nothing else on the server is
            touched.
        dry_run (bool): If `True`, only prints what would be truncated.
        assume_yes (bool): If `True`, skips the interactive confirmation
            prompt (for scripted/automated use -- the caller is asserting
            they already know exactly what this will do).

    Returns:
        bool: Whether a reset actually happened (`False` for `--dry-run`
            or a declined confirmation).
    """
    conn = get_connection(config)
    try:
        tables = _content_tables(conn, config.database)
        if not tables:
            print(f"'{config.database}' has no content tables to wipe (already empty schema).")
            return False

        counts = _row_counts(conn, tables)
        total_rows = sum(counts.values())

        if dry_run:
            print(f"Would truncate {len(tables)} table(s) in '{config.database}' ({total_rows:,} row(s) total):")
            for table in sorted(tables):
                print(f"  {table}: {counts[table]:,} row(s)")
            return False

        if not assume_yes and not _confirm(config.database, config.host, total_rows):
            print("Cancelled -- no changes made.")
            return False

        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        try:
            for table in tables:
                conn.execute(f"TRUNCATE TABLE {table}")
        finally:
            conn.execute("SET FOREIGN_KEY_CHECKS = 1")
        conn.commit()

        print(f"Wiped {len(tables)} table(s) ({total_rows:,} row(s)) from '{config.database}'.")
        print("Ready for a new galaxy -- e.g.:")
        print("  python3 generate.py plan ...   # rebuild the density skeleton")
        print("  python3 generate.py galaxy ... # generate sectors into it")
        return True
    finally:
        conn.close()


def main():
    parser = argparse.ArgumentParser(
        description="Wipe every generated sector/system/galaxy row from the configured planetGen "
                     "database, leaving its schema and control (admin) schema intact.",
    )
    add_mysql_connection_args(parser)
    parser.add_argument('--yes', '-y', action='store_true',
                         help="Skip the interactive confirmation prompt. For scripted/automated use "
                              "only -- this is a destructive, unrecoverable operation.")
    parser.add_argument('--dry-run', action='store_true',
                         help="List every table (and its row count) that would be wiped, without "
                              "changing anything.")
    args = parser.parse_args()

    config = mysql_config_from_args(args)
    did_reset = reset_database(config, dry_run=args.dry_run, assume_yes=args.yes)

    if not did_reset and not args.dry_run:
        sys.exit(1)


if __name__ == "__main__":
    main()
