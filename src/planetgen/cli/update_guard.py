#!/usr/bin/env python3
# planetgen.cli.update_guard

"""
What `update.sh` runs to close the site and clear the machine before it changes anything, and to open the site
again when it has finished (OPS.41, OPS.53; `planetgen.admin.updateguard` says how).

    python3 -m planetgen.cli.update_guard begin --pid PID [--cancel-running]
    python3 -m planetgen.cli.update_guard end --pid PID [--failed [--step NAME]]
    python3 -m planetgen.cli.update_guard status

`begin` lists every running or queued piece of work and goes on only if there is none or the operator agrees to
cancel it (a question on the terminal; where there is none, `--cancel-running`). It exits 0 to go on and 1 when
the update must not (nothing is changed then). `end` opens the site again, or, with `--failed`, leaves it closed
and says how to open it. `status` prints whether an update runs and what runs.
"""

import argparse
import sys

from planetgen.admin import maintenance, updateguard
from planetgen.db import store


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=("begin", "end", "status"))
    parser.add_argument("--pid", type=int, help="The update's process id (it holds the lock while it lives).")
    parser.add_argument("--cancel-running", action="store_true",
                        help="Cancel the running and queued work without asking.")
    parser.add_argument("--failed", action="store_true", help="With end: the update failed.")
    parser.add_argument("--step", default=None, help="With end --failed: the step that failed.")
    store.add_mysql_connection_args(parser)
    args = parser.parse_args(argv)
    config = store.mysql_config_from_args(args)

    if args.command == "status":
        lock = updateguard.read_lock()
        state = maintenance.read()
        if lock:
            print(f"An update is running (process {lock['pid']}, since {lock['since']}).")
        elif updateguard.is_updating():
            print(f"An update failed and left the site closed ({updateguard.TURN_OFF_HINT} opens it).")
        else:
            print("No update is running.")
        items = updateguard.list_running(updateguard.default_sources(config))
        print("Running or queued work:" if items else "Nothing is running or queued.")
        for line in updateguard.describe(items):
            print(line)
        print("Maintenance mode is " + ("ON" if state["enabled"] else "off") + ".")
        return 0
    if args.pid is None:
        parser.error(f"{args.command} needs --pid")
    if args.command == "begin":
        ok = updateguard.begin(args.pid, config=config, interactive=sys.stdin.isatty(),
                               cancel_running=args.cancel_running)
        return 0 if ok else 1
    updateguard.end(args.pid, ok=not args.failed, step=args.step)
    return 0


if __name__ == "__main__":
    sys.exit(main())
