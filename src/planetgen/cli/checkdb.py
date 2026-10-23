# planetgen.cli.checkdb

"""
`planetgen check-db`: checks the database for damage and changes nothing (DB.8, PERF.126).

    planetgen check-db [CHECK ...] [--shell N] [--sector SHELL BAND WEDGE] [--wide-binaries] [--yes]

Each check (`planetgen.db.check.CHECKS`) is a command of its own (`check-revision`, `check-tables`, `check-orphans`,
... listed by `planetgen commands`) that does that check and nothing else, reads the database and writes nothing, and
can run alone (`planetgen command-run`). This command only builds the list of them (`queue/wrapper.py`,
`commands.build_steps`): the checks named, or `commands.DEFAULT_CHECKS` (every check but the two slow ones) when none
is named; `tables-extended` and `systems` are the slow ones, and naming them is what a deep check is. The wrapper
queues them one after the other as a check run on the Redis queue (it does not own the machine, so it neither waits
for nor holds up a generation run and shows no banner to the site's visitors), and the answers are collected into the
report: one pass or fail line per check, the offenders under it, and the verdict.

It exits 1 when it found damage and 2 when a check could not run, a Redis that does not answer or a run already
going included. The systems check takes long: its time is estimated first and a yes asked for (`--yes` skips the
question, which is needed off a terminal; `--estimate-only` prints the estimate and stops, which the Generate page
uses).
"""

import json
import signal
import sys

from planetgen.db import check as db_check, store
from planetgen.generation import run_common
from planetgen.queue import commands, redisqueue, runrecord, units, wrapper
from planetgen.util import log


def add_check_db_arguments(parser):
    """
    Adds the `check-db` subcommand's options: the checks to run, the sectors to scope the value, count, version-key
    and systems checks to (`--shell`, `--sector`), `--wide-binaries`, `--yes`, `--estimate-only`, and the MySQL
    connection args (the logging args are `generate.py`'s).

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    parser.add_argument("checks", nargs="*", metavar="CHECK",
                        help=f"the checks to run (default: {' '.join(commands.DEFAULT_CHECKS)}; "
                             f"every one: {' '.join(commands.KIND_OF_CHECK)}; see 'planetgen commands' for the "
                             f"check-* commands)")
    parser.add_argument("--shell", type=int, action="append", metavar="SHELL",
                        help="Check only this shell's sectors for values, counts, version keys and systems "
                             "(repeatable).")
    parser.add_argument("--sector", type=int, nargs=3, action="append", metavar=("SHELL", "BAND", "WEDGE"),
                        help="Check only the sector at this address for values, counts, version keys and systems "
                             "(repeatable).")
    parser.add_argument("--wide-binaries", action="store_true",
                        help="With the systems check: validate only the wide binaries (S-type pairs), for the "
                             "count of pairs whose orbits overlap.")
    parser.add_argument("--yes", action="store_true",
                        help="With the systems check: go ahead after the estimate without asking (needed off a "
                             "terminal).")
    parser.add_argument("--estimate-only", action="store_true",
                        help="With the systems check: print the time estimate and stop (the Generate page uses it).")
    store.add_mysql_connection_args(parser)


def _stop_on_sigterm():
    """Makes SIGTERM end the run as cancelled in Redis, as Ctrl+C does."""
    def on_term(_signum, _frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, on_term)


def _confirm_systems(args, config, scope):
    """The estimate of the systems check: printed, and whether to go on: not with `--estimate-only`, and only with a yes
    from the terminal unless `--yes`."""
    conn = store.get_connection(config, ensure_schema=False)
    try:
        estimate = db_check.deep_estimate(conn, config, scope)
    finally:
        conn.close()
    log.normal(estimate.summary())
    if args.estimate_only:
        sys.stdout.write(run_common.ESTIMATE_PREFIX + json.dumps(estimate.as_dict()) + "\n")
        sys.stdout.flush()
        return False
    if args.yes:
        return True
    if not run_common._interactive():
        log.normal("Nothing was checked: the systems check off a terminal needs --yes.")
        raise SystemExit(db_check.EXIT_UNCHECKED)
    if run_common._ask(None, "Run the checks? [y/N] ") not in ("y", "yes"):
        log.normal("Nothing was checked.")
        return False
    return True


def run_check_db(args):
    """
    `planetgen check-db`: builds the list of check commands the arguments ask for, runs it through the wrapper and
    prints the report; exits 1 on damage, 2 when a check could not run.
    """
    config = store.mysql_config_from_args(args)
    args.database = wrapper.database_args(config)
    try:
        steps = wrapper.steps_for("check-db", args)
    except ValueError as exc:
        raise SystemExit(f"check-db: {exc}")
    asked = {step.key for step in steps}
    if args.estimate_only and "systems" not in asked:
        raise SystemExit("check-db: --estimate-only is for the systems check.")
    if "systems" in asked:
        scope = db_check.scope_of({"shells": args.shell, "addresses": args.sector,
                                   "wide_binaries": args.wide_binaries})
        if not _confirm_systems(args, config, scope):
            return
    _stop_on_sigterm()
    results = []

    def on_result(index, _unit, result):
        results.append(db_check.CheckResult.from_dict(result))

    try:
        wrapper.run("Database check", steps, lambda index, plan: None, on_result, database=config.database,
                    argv=sys.argv[1:], kind=runrecord.CHECK)
    except (redisqueue.Unavailable, runrecord.Busy, units.RunStopped) as exc:
        log.error(f"error: the database check could not run: {exc}")
        raise SystemExit(db_check.EXIT_UNCHECKED)
    report = db_check.Report(results)
    for line in report.lines():
        log.normal(line)
    if report.exit_code():
        raise SystemExit(report.exit_code())
