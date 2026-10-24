# planetgen/cli/queue.py

"""
`planetgen queue`: stop everything, hold the host until told otherwise, and look at the hold (PERF.186).

    planetgen queue stop-all [--no-hold]
    planetgen queue resume
    planetgen queue status

`stop-all` is the emergency stop (`planetgen.queue.stopall`): it writes the hold, cancels the Generate page's job and every
live run, empties every queue, and ends every planetgen worker, job runner and command process on this host with what
they started. While the hold is on nothing starts: not a worker, not a run, not a job, not the orbital update, not the
site's tile prefetch. `resume` lifts the hold. `status` says whether it is on, and what planetgen processes run. Run it as
root (or the web user) so every process can be ended; what could not be is listed.

`--no-hold` is for `update.sh --kill-all`: the same stop and clear, without leaving the hold on afterwards (the update's
own lock keeps new work out while it runs).

It also runs on its own, `python -m planetgen.cli.queue stop-all`, with nothing but the standard library and Redis
imported, so it works when the rest of the install is in trouble.
"""

import argparse
import sys
from datetime import datetime, timezone

from planetgen.queue import stopall


def add_queue_arguments(parser):
    """Adds the `queue` command's action and options."""
    parser.add_argument("action", choices=("stop-all", "resume", "status"),
                        help="stop-all: end everything and hold; resume: lift the hold; status: show it.")
    parser.add_argument("--no-hold", action="store_true",
                        help="With stop-all: end everything but do not leave the hold on (update.sh --kill-all).")


def _report(summary):
    lines = []
    if summary["already"]:
        lines.append("Everything was already on hold; stopped again.")
    lines.append(f"Cancelled {summary['runs']} live run(s)" + (" and the Generate page's job" if summary["job"] else "")
                 + f"; removed {summary['queued']} queued job(s); ended {summary['killed']} process(es).")
    if summary["survivors"]:
        lines.append(f"These processes could not be ended by this user (run it as root): "
                     f"{', '.join(str(pid) for pid in summary['survivors'])}")
    return lines


def run_queue(args):
    """Runs the action; returns the exit status (1 when some process could not be ended)."""
    if args.action == "stop-all":
        summary = stopall.stop_all("planetgen queue", hold=not args.no_hold)
        for line in _report(summary):
            print(line)
        if not args.no_hold:
            print("Nothing will start until you resume: planetgen queue resume (or Resume on the admin Queue page).")
        return 1 if summary["survivors"] else 0
    if args.action == "resume":
        was = stopall.resume_all()
        print("Resumed: work can start again." if was else "It was not stopped.")
        return 0
    held = stopall.state()
    if held["stopped"]:
        when = datetime.fromtimestamp(held["at"], timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ") if held["at"] else "an unknown time"
        print(f"Everything is on hold since {when}" + (f", stopped by {held['by']}" if held["by"] else "") + ".")
    else:
        print("Nothing is on hold.")
    table = stopall.process_table()
    running = stopall.find_targets(table)
    print(f"{len(running)} planetgen worker or command process(es) running.")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog="planetgen.cli.queue", description=__doc__.strip().splitlines()[0])
    add_queue_arguments(parser)
    return run_queue(parser.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
