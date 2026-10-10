#!/usr/bin/env python3
# planetgen.cli.dedupe

"""
One-off backfill: scans an existing database for sector/system names
that duplicate each other (within or across those two levels) and
resolves them with the same Greek/Roman and diminutive decoration
`planetgen/db/store.py`'s `insert_sector`/`insert_star_system` already
apply automatically to every *new* row (v24, see
`planetgen/names/uniqueness.py`'s own module docstring for the
sector > system hierarchy) -- for a database that predates that feature,
or already holds duplicate names from before it existed. Planets and
moons are named from their system (`planetgen/names/bodies.py`, v34), so
a renamed system's derived names follow it (`store.rename_star_system`).

Live generation already guarantees no new duplicate ever lands in the
database going forward; this script is only for cleaning up whatever's
already there. Processes one level at a time, top-down (sectors, then
systems) -- exactly the order `nameUniqueness`'s own hierarchy needs,
since a system's cross-level check against `sector_name_registry` has to
see the sectors' *final*, already-resolved state. Reuses `store.py`'s own `reserve_*_name`/`confirm_*_name` functions
directly -- the same two-phase reservation `insert_sector`/etc. call --
rather than a second, parallel implementation of the same rules.

Idempotent: a base name whose registry `occurrence_count` already matches
how many rows currently share it is left untouched, so re-running this
script against a database it's already cleaned (with no new duplicates
added by some other means in between) is a no-op.

Run it as `python3 -m planetgen.cli.dedupe` from anywhere (the editable
install makes `planetgen` importable).

Usage:
    python3 -m planetgen.cli.dedupe [--mysql-host HOST] [--mysql-port PORT]
                               [--mysql-user USER] [--mysql-password PASSWORD]
                               [--mysql-database DATABASE]

    Every flag defaults to the same $PLANETGEN_MYSQL_* environment
    variable every other entry point in this project reads (see
    `planetgen.db.store.MySQLConfig`). The configured account needs
    ordinary `SELECT`/`INSERT`/`UPDATE` grants on the content database --
    no `CREATE`/`ALTER` (this never touches DDL, unlike `planetgen.cli.migrate`).
"""

import argparse
import sys
from collections import defaultdict

import pymysql

from planetgen.db import store
from planetgen.generation import steps
from planetgen.names.uniqueness import strip_decoration


def _group_by_base(rows):
    """`rows`: an iterable of dicts with `id`/`name`. Returns `{base_name:
    [rows sharing it, sorted by id ascending]}` -- id order is this
    script's stand-in for creation order (see the module docstring's note
    on why that's not available directly)."""
    groups = defaultdict(list)
    for row in rows:
        groups[strip_decoration(row["name"])].append(row)
    for base in groups:
        groups[base].sort(key=lambda r: r["id"])
    return groups


def _already_resolved_count(conn, registry_table, base_name):
    """How many rows this base name's registry row already accounts for
    -- `0` if no registry row exists yet (nothing resolved so far)."""
    row = conn.execute(
        f"SELECT occurrence_count FROM {registry_table} WHERE base_name = ?", (base_name,),
    ).fetchone()
    return row["occurrence_count"] if row else 0


def _display(show):
    """`steps.UNSET` (the open display, if any) for a command-line run, `None` (draw nothing) otherwise."""
    return steps.UNSET if show else None


def _dedupe_sectors(conn, show=False):
    """Resolves every sector-vs-sector duplicate, and (via `reserve_sector_name`'s
    own cross-level check) retroactively decorates any already-existing
    system that happens to share a sector's base name. Must run before
    `_dedupe_systems` -- see the module docstring.

    Returns:
        int: How many `sectors` rows this call actually renamed.
    """
    rows = conn.execute("SELECT id, name FROM sectors ORDER BY id").fetchall()
    renamed = 0
    groups = _group_by_base(rows)
    with steps.Step("Sector names", "dedupe-sectors", max(len(groups), 1), own=show, progress=_display(show)) as bar:
        for base, group in groups.items():
            already_done = _already_resolved_count(conn, "sector_name_registry", base)
            for row in group[already_done:]:
                new_name, _name_base = store.reserve_sector_name(conn, base, row["id"])
                if new_name != row["name"]:
                    conn.execute("UPDATE sectors SET name = ? WHERE id = ?", (new_name, row["id"]))
                    renamed += 1
            bar.advance()
    return renamed


def _dedupe_systems(conn, show=False):
    """Resolves every system-vs-system duplicate and every system-vs-sector
    collision (diminutive prefix, system side only). Must run after
    `_dedupe_sectors`.

    Returns:
        int: How many `star_systems` rows this call actually renamed.
    """
    rows = conn.execute("SELECT id, name FROM star_systems ORDER BY id").fetchall()
    renamed = 0
    groups = _group_by_base(rows)
    with steps.Step("System names", "dedupe-systems", max(len(groups), 1), own=show, progress=_display(show)) as bar:
        for base, group in groups.items():
            already_done = _already_resolved_count(conn, "system_name_registry", base)
            for row in group[already_done:]:
                new_name, name_base, diminutive_index = store.reserve_system_name(conn, base)
                if new_name != row["name"]:
                    store.rename_star_system(conn, row["id"], new_name)
                    renamed += 1
                store.confirm_system_name(conn, name_base, row["id"], diminutive_index)
            bar.advance()
    return renamed


def dedupe_names(config=None, show=False):
    """
    Runs the full sector -> system dedup pass, in one
    transaction (commits only if every step succeeds).

    Args:
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.
        show (bool): Draw a bar per level when it is long (the command
            line; `planetgen.generation.steps`, UX.84).

    Returns:
        dict: `sectors`/`star_systems`, each the
            count of rows this run actually renamed (`0` for every key
            means nothing needed fixing).
    """
    conn = store.get_connection(config)
    try:
        with conn:
            sectors_renamed = _dedupe_sectors(conn, show)
            systems_renamed = _dedupe_systems(conn, show)
        return {
            "sectors": sectors_renamed,
            "star_systems": systems_renamed,
        }
    finally:
        conn.close()


def main():
    parser = argparse.ArgumentParser(
        description="Scan the configured database for sector/system names that duplicate "
                     "each other (within or across those two levels) and resolve them with the same "
                     "Greek/Roman and diminutive decoration new generation runs already apply "
                     "automatically (see planetgen/names/uniqueness.py).",
    )
    store.add_mysql_connection_args(parser)
    args = parser.parse_args()

    config = store.mysql_config_from_args(args)
    try:
        counts = dedupe_names(config, show=True)
    except (pymysql.MySQLError, store.SchemaTooNewError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)

    total = sum(counts.values())
    if total == 0:
        print("No duplicate names found -- nothing to do.")
    else:
        print(
            f"Renamed {counts['sectors']} sector(s) and {counts['star_systems']} system(s) "
            f"to resolve name collisions."
        )


if __name__ == "__main__":
    main()
