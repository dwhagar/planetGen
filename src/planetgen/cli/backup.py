#!/usr/bin/env python3
# planetgen.cli.backup

"""
Takes, lists, verifies and restores backup archives (ADM.87).

    python3 -m planetgen.cli.backup take [--part NAME ...]   # one archive of every part, or of the named ones
    python3 -m planetgen.cli.backup list                     # the archives, newest first, with their parts
    python3 -m planetgen.cli.backup verify ARCHIVE           # check every hash
    python3 -m planetgen.cli.backup restore ARCHIVE [--part NAME ...]
                                                             # the whole archive, or only the named parts
    python3 -m planetgen.cli.backup parts                    # the parts this build can back up

An archive is one `backups/<UTC time>.tar.bz2` in the state folder (`/var/lib/planetGen/state`), holding a part
per kind of data: `settings` (your `config.json`, secrets included) now, the control database and the DB repair parity
data later. `settings.backups_kept` of them are kept; the newest one holding every part is never deleted. A
restore verifies the archive first, backs up what it is about to replace in a safety archive, and skips (and
reports) a part this build does not know. The archive files are mode 0600, owned by the web user: run the
command as that user or as root. Each take and restore is written to the activity log, and a restore to the
audit table when the control database is reachable.
"""

import argparse
import sys

from planetgen import backup
from planetgen.cli import config as config_cli
from planetgen.util import log


def _size(value):
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024 or unit == "GB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024


def run_parts():
    for name, part in backup.registered().items():
        print(f"{name}  {part.label}")
    return 0


def run_list():
    rows = backup.list_archives()
    if not rows:
        print("There are no backups yet.")
    for row in rows:
        if not row.readable:
            print(f"{row.name}  (manifest unreadable, {_size(row.size)})")
            continue
        parts = ", ".join(f"{p.name} {_size(p.size)}{'' if p.known else ' (unknown here)'}" for p in row.parts)
        print(f"{row.name}  {row.taken_at}  {row.taken_by or '-'}  {row.reason or '-'}  [{parts}]"
              f"{'' if row.complete else '  (not every part)'}")
    return 0


def run_take(args):
    try:
        archive = backup.take(args.part or None, who=config_cli._actor(), reason="manual")
    except backup.BackupError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Backed up {', '.join(archive.part_names)}: {archive.name}")
    return 0


def run_verify(args):
    try:
        checked = backup.verify(args.archive)
    except backup.BackupError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    for problem in checked.problems:
        print(f"problem: {problem}", file=sys.stderr)
    print(f"{args.archive}: " + ("sound." if checked.ok else "NOT sound."))
    return 0 if checked.ok else 1


def run_restore(args):
    audit, close = config_cli._open_audit(args)
    try:
        report = backup.restore(args.archive, args.part or None,
                                backup.RestoreContext(who=config_cli._actor(), role="host", audit=audit))
    except backup.BackupError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    finally:
        close()
    if report.restored:
        print(f"Restored {', '.join(report.restored)} from {args.archive}. What they replaced is in {report.safety}.")
    for name in report.skipped:
        print(f"skipped: this build has no provider for the {name!r} part; the archive keeps it.", file=sys.stderr)
    if not report.restored and not report.skipped:
        print("Nothing to restore.")
    return 0 if report.restored or not report.skipped else 1


def main(argv=None):
    from planetgen.db.store import add_mysql_connection_args
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    take = commands.add_parser("take", help="take a backup now")
    take.add_argument("--part", action="append", help="back up only this part (repeatable)")
    commands.add_parser("list", help="list the archives")
    commands.add_parser("parts", help="list the parts this build can back up")
    verify = commands.add_parser("verify", help="check an archive's hashes")
    verify.add_argument("archive")
    restore = commands.add_parser("restore", help="restore an archive or some of its parts")
    restore.add_argument("archive")
    restore.add_argument("--part", action="append", help="restore only this part (repeatable)")
    add_mysql_connection_args(restore)
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.apply_debug_arguments(args, parser)
    return {"take": lambda: run_take(args), "list": run_list, "parts": run_parts,
            "verify": lambda: run_verify(args), "restore": lambda: run_restore(args)}[args.command]()


if __name__ == "__main__":
    sys.exit(main())
