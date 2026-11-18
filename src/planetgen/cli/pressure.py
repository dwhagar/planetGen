# planetgen.cli.pressure

"""
Reads the record the job engine's memory controller keeps of the host's pressure (OPS.249):

    python3 -m planetgen.cli.pressure report                    # what the record holds and the edges it gives
    python3 -m planetgen.cli.pressure report --csv FILE         # also write every recorded minute to FILE and the kills to FILE-kills

The record is `pressure_minutes` of the control database, one row a minute while a supervisor runs (`planetgen.queue.pressurelearn`), with the ledger's memory kills beside it. The report says how
much record there is, how many of the kills the hold did not see coming and how many holds were not needed (what each rule of the learning needs before it moves an edge), how many minutes carry a probe of
the database (the record the factors of `planetgen.queue.dbsignal` are to be learned from), and the edges the record gives now. It reads only: nothing is changed, whatever `PLANETGEN_PRESSURE_LEARN` says.
A gate or a trial run ends with `report --csv` so its record is kept beside its other results.
"""

import argparse
import csv
import sys
import time

from planetgen.db import store
from planetgen.queue import pressurelearn
from planetgen.util import log
from planetgen.util.format import format_duration_seconds


def lines(found):
    """The sentences of a `pressurelearn.Report`, one a line."""
    if not found.minutes:
        return [f"The record is empty: no minute of the last {pressurelearn.WINDOW_DAYS} days was recorded (no supervisor ran, or {pressurelearn.ENV_VAR} was off), so the typed edges stand."]
    edges = found.edges
    return [
        f"Record: {found.minutes} minutes over {format_duration_seconds(found.span_seconds)} of the last {pressurelearn.WINDOW_DAYS} days, {found.psi_minutes} of them read by PSI.",
        f"Memory kills: {found.kills}, of which {found.unforeseen_kills} no pressed reading came before (the edges move down at {pressurelearn.MIN_KILLS}).",
        f"Holds: {found.holds}, of which {found.idle_holds} were not needed (the edges move up at {pressurelearn.MIN_EPISODES} holds with {pressurelearn.IDLE_SHARE:.0%} of them not needed and no kill).",
        f"Database probes: {found.probe_minutes} minutes carry one" + (f", the slowest took {found.slowest_probe_ms:.0f} ms." if found.probe_minutes else "."),
        f"Edges now: sensitivity {edges.sensitivity:.2f}, {edges.words()}.",
    ]


def main(argv=None):
    parser = argparse.ArgumentParser(description="Read the memory controller's record of the host's pressure.")
    commands = parser.add_subparsers(dest="command", required=True)
    shown = commands.add_parser("report", help="say what the record holds and the edges it gives")
    shown.add_argument("--csv", metavar="FILE", help="write every recorded minute of the window to FILE, and the times of the memory kills to FILE with -kills before its extension")
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.apply_debug_arguments(args, parser, errors_to_stderr=True)
    try:
        conn = store.get_control_connection()
        try:
            minutes, kills = pressurelearn.read_record(conn)
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- the command says why it could not read and exits 1
        log.error(f"The control database could not be read ({exc}).")
        return 1
    for line in lines(pressurelearn.report(minutes, kills)):
        print(line)
    if args.csv:
        with open(args.csv, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(pressurelearn.CSV_COLUMNS)
            writer.writerows(pressurelearn.csv_rows(minutes))
        stem, dot, extension = args.csv.rpartition(".")
        kills_path = f"{stem}-kills{dot}{extension}" if dot and stem else f"{args.csv}-kills"
        with open(kills_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("killed_at_utc",))
            writer.writerows((time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(killed)),) for killed in kills)
        print(f"Wrote {len(minutes)} minutes to {args.csv} and {len(kills)} kills to {kills_path}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
