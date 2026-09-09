#!/usr/bin/env python3
# src/migrateSqliteToMysql.py

"""
One-time import of an existing planetGen SQLite database (from before the
MySQL port, TODO.md Phase 5) into a MySQL database.

Only migrates a SQLite database already at schema v5 (`PRAGMA
user_version`) -- the version every schema-version-N-to-N+1 migration
this project ever shipped converted up to, and the same version this
project's MySQL schema (`stellarObjects/schema.sql`) starts at. A SQLite
database still on an older schema needs to go through a pre-MySQL-port
release of this project first (any version through `migrateDb.py`'s old
SQLite-to-SQLite `stellarObjects._db.migrate_database`, since removed --
see git history/CHANGELOG.md for the last release that still had it) to
reach v5, then through this script.

The copy is column-name-preserving and generic across all 13 tables
(`_TABLES_IN_FK_ORDER` below): every table in the current schema has the
exact same column names in both the SQLite source and the MySQL
destination (only column *types* changed during the port -- see
`schema.sql`'s "MySQL port" header note), including primary key `id`
columns, which are copied verbatim rather than remapped. This is what
lets every foreign key (`star_system_id`, `planet_id`, `belt_id`, ...)
keep pointing at the same logical row across the copy with no id-mapping
table of its own, and lets this script use one small generic
`_copy_table` helper instead of 13 hand-written per-table copies. MySQL's
`AUTO_INCREMENT` counter self-adjusts to stay past the highest id it's
ever seen inserted (standard InnoDB behavior), so generation continues
correctly afterward without any manual counter reset.

Usage:
    python3 src/migrateSqliteToMysql.py path/to/old.db [--mysql-host ...]

    Every --mysql-* flag is the same as every other entry point in this
    project (see `stellarObjects._db.MySQLConfig`) -- this typically
    wants the read-write account (`sectorGen.py`/`systemGen.py` use),
    not a read-only one, since it writes every migrated row.
"""

import argparse
import sqlite3
import sys

from stellarObjects._db import (
    SCHEMA_VERSION,
    add_mysql_connection_args,
    get_connection,
    mysql_config_from_args,
)

_TABLES_IN_FK_ORDER = (
    "sectors",
    "system_configs",
    "system_config_slots",
    "star_systems",
    "stars",
    "planets",
    "planet_evolutionary_paragraphs",
    "planet_reflection_spectrum",
    "moons",
    "moon_evolutionary_paragraphs",
    "moon_reflection_spectrum",
    "asteroid_belts",
    "asteroid_belt_composition",
)
"""tuple: Every data table, in an order that never inserts a row before a
table it foreign-keys to -- see the module docstring. `schema_migrations`
and the `sector_objects` view are deliberately excluded: the former is
schema-version bookkeeping the destination's own `get_connection` already
seeds, not source data to copy, and the latter has no rows of its own to
copy at all."""

_BATCH_SIZE = 500
"""int: Rows read from SQLite and committed to MySQL per batch, so a
large table doesn't hold its entire row set in memory at once or leave
one unbroken multi-million-row transaction open."""


def _sqlite_columns(sqlite_conn, table):
    """Returns `table`'s column names, in schema order, via SQLite's
    `PRAGMA table_info` -- used instead of a hardcoded per-table column
    list so this script keeps working even if a future schema version
    adds/removes columns, as long as both databases still share the same
    names for the columns that exist in both."""
    rows = sqlite_conn.execute(f"PRAGMA table_info({table})").fetchall()
    return [row[1] for row in rows]  # row[1] is the column name


def _copy_table(sqlite_conn, mysql_conn, table):
    """
    Copies every row of `table` from `sqlite_conn` to `mysql_conn`,
    column-for-column (see the module docstring), in batches of
    `_BATCH_SIZE`.

    Returns:
        int: Number of rows copied.
    """
    columns = _sqlite_columns(sqlite_conn, table)
    column_list = ", ".join(columns)
    placeholders = ", ".join("?" for _ in columns)
    insert_sql = f"INSERT INTO {table} ({column_list}) VALUES ({placeholders})"

    cursor = sqlite_conn.execute(f"SELECT {column_list} FROM {table}")
    total = 0
    while True:
        batch = cursor.fetchmany(_BATCH_SIZE)
        if not batch:
            break
        for row in batch:
            mysql_conn.execute(insert_sql, tuple(row))
        mysql_conn.commit()
        total += len(batch)
    return total


def migrate(sqlite_path, mysql_config):
    """
    Migrates one SQLite database file into the MySQL database `mysql_config`
    points at.

    Args:
        sqlite_path (str): Path to the source `.db` file.
        mysql_config (MySQLConfig): Destination connection parameters.

    Raises:
        SystemExit: If the source isn't schema v5, or can't be opened.
    """
    try:
        sqlite_conn = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
    except sqlite3.OperationalError as exc:
        raise SystemExit(f"Error: could not open {sqlite_path!r} ({exc}).")

    try:
        version = sqlite_conn.execute("PRAGMA user_version").fetchone()[0]
        if version != SCHEMA_VERSION:
            raise SystemExit(
                f"Error: {sqlite_path!r} is at schema v{version}, expected v{SCHEMA_VERSION}. "
                f"Migrate it to v{SCHEMA_VERSION} with a pre-MySQL-port release of this project first "
                f"(see this script's module docstring)."
            )

        mysql_conn = get_connection(mysql_config)
        try:
            for table in _TABLES_IN_FK_ORDER:
                count = _copy_table(sqlite_conn, mysql_conn, table)
                print(f"{table}: {count} rows copied")
        finally:
            mysql_conn.close()
    finally:
        sqlite_conn.close()


def main():
    parser = argparse.ArgumentParser(
        description="One-time import of an existing schema-v5 SQLite database into MySQL.",
    )
    parser.add_argument("sqlite_path", help="Path to the source SQLite .db file.")
    add_mysql_connection_args(parser)
    args = parser.parse_args()

    migrate(args.sqlite_path, mysql_config_from_args(args))
    print("Done.")


if __name__ == "__main__":
    main()
