#!/usr/bin/env python3
# planetgen.cli.migrate

"""
Brings the configured planetGen MySQL database's `schema_migrations`
bookkeeping up to the current schema (`planetgen/db/schema.sql`),
applying any migration step in between -- a no-op for a database that's
already current. See `planetgen/db/store.py`'s `migrate_database` for how
a single database's version is checked/advanced.

Run automatically by `install.sh` and `update.sh`
on every deploy, so a database created under an older schema keeps
working after a `git pull` brings in a newer one. Also runnable directly
for a one-off check/migration outside of a deployment.

Run it as `python3 -m planetgen.cli.migrate` from anywhere (the editable
install makes `planetgen` importable).

Shows a progress bar while it migrates: one tick per migration step,
with the elapsed time and an estimate of the time left, and under it a
bar for a step that works in batches (`alembic_runner.report_progress`),
both on `planetgen.generation.steps` (DB.15). When output is not a terminal
it prints a line per step instead, and for a step in batches at most one
every 30 seconds and at each 10 percent.

Usage:
    python3 -m planetgen.cli.migrate [--mysql-host HOST] [--mysql-port PORT]
                             [--mysql-user USER] [--mysql-password PASSWORD]
                             [--mysql-database DATABASE] [--status]

    --status prints "<current version> <target version> <pending steps>
    <database>" and changes nothing (update.sh reads it to decide whether
    to ask about the database).

    Every flag defaults to the same $PLANETGEN_MYSQL_* environment
    variable every other entry point in this project reads (see
    `planetgen.db.store.MySQLConfig`) -- unlike the pre-MySQL-port version
    of this script, there is exactly one database to migrate (a MySQL
    server, not a directory of `*.db` files), so this needs a connection
    to point at rather than a directory to scan.
"""

import argparse
import contextlib
import os
import sys
import time

from planetgen.admin import auth, maintenance, updateguard
from planetgen.cli import stage_progress
from planetgen.generation import run_common
from planetgen.db.store import (
    SCHEMA_VERSION,
    add_mysql_connection_args,
    control_mysql_config,
    migrate_database,
    mysql_config_from_args,
    schema_status,
)


LINE_SECONDS = 30.0
"""float: At most one progress line this often when output is not a terminal."""


class _LineReporter:
    """What stands in for the bars when output is not a terminal (a log, `update.sh`
    piped): one line per step, and for a step that works in batches at most one line
    every `LINE_SECONDS` and at each 10 percent."""

    def __init__(self, write=print, clock=time.monotonic):
        self.write, self.clock = write, clock
        self.label = None
        self.last_time = self.last_decile = None

    def on_step(self, number, total, from_version, to_version):
        self.write(f"Migrating v{from_version} -> v{to_version} (step {number} of {total})")
        self.label = None

    def on_progress(self, label, done, total):
        now = self.clock()
        decile = int(10 * done / total) if total else 0
        if label != self.label:
            self.label, self.last_time, self.last_decile = label, None, -1
        if done < total and self.last_time is not None and decile == self.last_decile \
                and now - self.last_time < LINE_SECONDS:
            return
        if done >= total and self.last_decile == 10:
            return
        self.last_time, self.last_decile = now, decile
        self.write(f"  {label}: {done:,} of {total:,} ({100 * done // total if total else 100}%)")


def _migrate_with_progress(config):
    """Runs `migrate_database` with a bar over the steps pending and, under
    it, a bar for a step that works in batches (`alembic_runner.report_progress`),
    both with the time left; plain lines instead when output is not a
    terminal. Shows nothing when there is nothing to migrate."""
    if not sys.stdout.isatty():
        lines = _LineReporter()
        return migrate_database(config, on_step=lines.on_step, on_progress=lines.on_progress)

    stage = None

    def on_step(number, total, from_version, to_version):
        nonlocal stage
        if stage is None:
            stage = stage_progress.StageProgress(total, kind="migrate", stats=_stats_or_none(config))
            stage.__enter__()
        stage.stage(f"Migrating v{from_version} -> v{to_version} (step {number} of {total})")

    def on_progress(label, done, total):
        if stage is not None:
            stage.detail(label, done, total)

    failed = True
    try:
        version = migrate_database(config, on_step=on_step, on_progress=on_progress)
        failed = False
    finally:
        if stage is not None:
            stage.__exit__(*((RuntimeError, None, None) if failed else (None, None, None)))
    return version


def _stats_or_none(config):
    """The stored speeds the bars' time left starts from, or `None` when the control database is
    not there yet (a first install)."""
    try:
        return run_common._stats_for_config(config)
    except Exception:  # noqa: BLE001 -- statistics never stop a migration
        return None


def print_initial_admin_login(username, password):
    """
    Shows the first admin login `auth.bootstrap_control_schema` just
    seeded. This is the only time the password exists in plain text (only
    its hash is stored), so it goes to the console once, on the run that
    created it, and never again.
    """
    print()
    print("=" * 72)
    print("New admin login password (shown only this once):")
    print(f"    username: {username}")
    print(f"    password: {password}")
    print("Log in at /login now and change both; the admin pages stay locked")
    print("until you do. If you lose this password, see docs/reference/api.md")
    print("(\"Resetting the admin login\").")
    print("=" * 72)
    print()


SECTOR_SWITCH_VERSION = 91
"""int: The schema version that converts the sector grid from the cylinder to shells, bands and wedges."""


def close_site_for_the_sector_switch(config, write=print):
    """
    Turns maintenance mode on (UX.112) before a database that still has the
    cylinder sector grid is migrated past `SECTOR_SWITCH_VERSION`, so nobody
    but an admin sees the site while the conversion and the map rebuild run.
    The map rebuild turns it off with `maintenance.turn_off_for(maintenance.SECTOR_SWITCH)`.
    A fresh database (no version yet), one already past the switch, and one
    an admin already closed are left alone; a file the web user's folder
    won't take only warns, because the migration is what matters.
    """
    try:
        version, _pending = schema_status(config)
    except Exception:  # noqa: BLE001 -- the migration itself reports what is wrong
        return False
    state = maintenance.read()
    if not 0 < version < SECTOR_SWITCH_VERSION <= SCHEMA_VERSION:
        return False
    if state["enabled"] and state["cause"] != updateguard.UPDATE_CAUSE:
        return False        # an admin's own closing stays theirs; an update's is taken over (the warm-up opens it)
    try:
        maintenance.turn_on(by="migration", reason=maintenance.SECTOR_SWITCH_REASON,
                            cause=maintenance.SECTOR_SWITCH)
    except OSError as exc:
        write(f"warning: could not turn maintenance mode on ({exc}); turn it on by hand "
              f"(python3 -m planetgen.cli.maintenance on) if visitors should not see the site mid-conversion.")
        return False
    write("Maintenance mode is on: only admins can use the site until the new Galaxy Map is built.")
    return True


@contextlib.contextmanager
def _update_guard(config, cancel_running):
    """
    OPS.41 and OPS.53 for a migration run by hand: when there is something to migrate and `update.sh` is not
    running it (it holds the guard already), list the running work and ask before cancelling it, close the site
    for the migration and open it again when it has finished. A failure leaves the site closed (the message says
    how to open it). Exits with status 1, having changed nothing, when the work is not cancelled.
    """
    if updateguard.inside_update():
        yield
        return
    try:
        _version, pending = schema_status(config)
    except Exception:  # noqa: BLE001 -- the migration reports what is wrong
        pending = 0
    if not pending:
        yield
        return
    pid = os.getpid()
    if not updateguard.begin(pid, by="migrate", config=config, interactive=sys.stdin.isatty(),
                             cancel_running=cancel_running):
        sys.exit(1)
    ok = False
    try:
        yield
        ok = True
    except SystemExit as exc:
        ok = exc.code in (None, 0)
        raise
    finally:
        updateguard.end(pid, ok=ok, step="the migration")


def main():
    parser = argparse.ArgumentParser(
        description="Bring the configured MySQL database's schema_migrations bookkeeping up to date, "
                     "and the control schema (admin logins/sessions/API keys -- see control_schema.sql) "
                     "alongside it.",
    )
    add_mysql_connection_args(parser)
    parser.add_argument("--cancel-running", action="store_true",
                        help="Cancel running and queued work without asking, when a migration is pending "
                             "and this is not run by update.sh.")
    parser.add_argument("--status", action="store_true",
                        help="Print '<current version> <target version> <pending steps> <database>' "
                             "and change nothing.")
    args = parser.parse_args()

    config = mysql_config_from_args(args)

    if args.status:
        try:
            version, pending = schema_status(config)
        except Exception as exc:
            print(f"error: {exc}", file=sys.stderr)
            sys.exit(1)
        print(version, SCHEMA_VERSION, pending, config.database)
        return

    with _update_guard(config, args.cancel_running):
        close_site_for_the_sector_switch(config)
        try:
            version = _migrate_with_progress(config)
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
        # seeded here too so a fresh deploy's first admin login exists without
        # a separate manual step. Reuses this same account's host/user/password
        # (the full-access account this script already runs as, same as
        # install.sh's existing step 2) against the control schema's own name
        # instead of the content database's.
        try:
            seeded = auth.bootstrap_control_schema(control_mysql_config(config))
        except Exception as exc:
            print(f"error: could not set up the control schema ({exc}).", file=sys.stderr)
            sys.exit(1)
        if seeded is not None:
            print_initial_admin_login(*seeded)
        print("Control schema (admin logins) is up to date.")


if __name__ == "__main__":
    main()
