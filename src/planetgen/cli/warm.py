# planetgen.cli.warm

"""
`planetgen warm`: builds the site's caches ahead of the first visitor (MAP.134, PERF.104, PERF.127): the Galaxy Map's
opening view, every list, every sector, system and phenomenon page, and the search, so the first visit after a release
(or after a cache was cleared) isn't the one who waits for it. update.sh runs it in two stages (`planetgen warm opening --progress` before the Apache
reload, so the reload and the health check find what the site opens with, then `planetgen warm --progress` after the
health check, which update.sh waits for); run it by hand after clearing a cache. Run it as the web server's user so the files it writes are readable by Apache.

    planetgen warm [FAMILY|PHASE ...] [--cap N]

The warm-up also starts by itself (`planetgen.queue.warmstart`, PERF.141): after a generation run, `planetgen orbits`,
an edit or anything else that moves a data version tag, a background `planetgen warm --on-change DATABASE` waits
`cache.warm_settle_seconds` for more changes, then warms once (and again if more came meanwhile). Nothing starts during
an update (it warms at its end) or with `cache.warm_on_change` off.

The warm-up runs in five phases (`planetgen.web.warmup.PHASES`, PERF.141), each stopping at the cache budget
(`cache.max_mb`): (1) `opening`: the opening views and the first page of each list, (2) `built-levels`: the Galaxy
Map's nodes to the built levels, (3) `deep-nodes`: the deep nodes along the zoom run, (4) `objects`: the sector, system
and phenomenon pages, (5) `search`: the search page. After each phase it prints one line: what it built, what was
already stored, what the budget or a cap left out, or `unchanged, nothing built`. A phase name stands for all of its
families.

Each family (`planetgen.web.warmup.FAMILIES`) is a command of its own (`warm-galaxy-map`, `warm-lists`, ... listed by
`planetgen commands`) that warms that family and nothing else, and can run alone (`planetgen command-run`). This
command only builds the list of them (`queue/wrapper.py`, `commands.build_steps`), every family in order when none is
named, and the wrapper queues them one after the other as a maintenance run on the Redis queue (it does not own the
machine, so it neither waits for nor holds up a generation run and shows no banner to the site's visitors). It prints
one line per family, with how many pages it built, how many were already cached and how many a cap left out
(`page_cache.warm_max_objects`, or `--cap`). It exits non-zero when a family failed, saying which; a Redis that does
not answer, or an update still running when the run would start, is an error printed, never a silent skip. It waits
for a running update to finish before it queues anything.

The answers go to Redis and to the disk cache the web server reads, so they are the ones a visitor gets; each is kept
only while the data it read is unchanged (docs/design/operations/cache-warm-up.md). A list of every family ends with the command
that opens the site after the sector-grid switch closed it (`maintenance.SECTOR_SWITCH`, UX.112): the site is closed
only while the maps are being rebuilt.
"""

import signal
import sys
import time

from planetgen.admin import updateguard
from planetgen.cli import warm_progress
from planetgen.queue import redisqueue, runrecord, units, warmphases, warmstart, wrapper
from planetgen.util import log

UPDATE_WAIT_SECONDS = 1200
"""int: How long it waits for a running update to finish before it starts the run."""


class Failed(Exception):
    """The warm-up could not run or a command of it failed; the message says what."""


def add_warm_arguments(parser):
    """Adds the `warm` subcommand's options: the families to warm and the cap."""
    parser.add_argument("families", nargs="*", metavar="FAMILY",
                        help="the cache families or phases to warm, in order (default: every family; a phase is all "
                             "its families; see 'planetgen commands' for the warm-* commands)")
    parser.add_argument("--cap", type=int, default=None,
                        help="most sector, system and phenomenon pages to warm (default page_cache.warm_max_objects)")
    parser.add_argument("--progress", action="store_true",
                        help="show progress: rich bars with a terminal (the stage's chunk N of M, elapsed, time left), "
                             "otherwise one line per chunk (update.sh shows it while it waits for the warm-up)")
    parser.add_argument("--on-change", metavar="DATABASE", default=None,
                        help="what a change starts (warmstart.request): wait cache.warm_settle_seconds, then warm once "
                             "for each time DATABASE was changed meanwhile; does nothing when DATABASE is not the one "
                             "the site shows")


def _wait_for_update():
    """Waits for a running update to finish: the queue takes no job during one (OPS.53)."""
    if not updateguard.is_updating() or updateguard.inside_update():
        return
    log.normal("The update is still running: waiting for it to finish before starting the warm-up.")
    deadline = time.monotonic() + UPDATE_WAIT_SECONDS
    while updateguard.is_updating():
        if time.monotonic() > deadline:
            raise Failed(f"the update was still running after {UPDATE_WAIT_SECONDS} s")
        time.sleep(2.0)


def _site_database():
    """The database the site shows (`web.helpers.db_name`), which only an app context can say."""
    from planetgen.web import warmup
    from planetgen.web.helpers import db_name
    with warmup.build_app().app_context():
        return db_name()


def _stop_on_sigterm():
    """Makes SIGTERM (the update guard's cancel, update.sh's timeout) end the run as cancelled in Redis, as Ctrl+C does."""
    def on_term(_signum, _frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, on_term)


def _warm_once(args, reason=None):
    """
    Builds the list of cache-warming commands the arguments ask for and runs it through the wrapper, printing each
    command's line.

    Raises:
        Failed: A command failed, Redis did not answer or another warm-up runs.
        SystemExit: The arguments name no such family or phase.
    """
    from planetgen.web import warmup
    try:
        steps = wrapper.steps_for("warm", args)
    except ValueError as exc:
        raise SystemExit(f"warm: {exc}")
    # Registered for the update guard (OPS.53), which lists it and stops it, the caches being safe to build again.
    with updateguard.running("warm-map", "Site cache warm-up", refuse=reason is None):
        try:
            _wait_for_update()
            tally = warmup.PhaseTally([step.key for step in steps])
            phases = [name for name, _text, families in warmup.PHASES if any(step.key in families for step in steps)]
            warmphases.begin(args.database, phases, reason)
            bars = warm_progress.WarmProgress(
                [step.key for step in steps], log.normal, show=getattr(args, "progress", False),
                rich_bars=warm_progress.stdout_is_terminal())

            def report_ended():
                for name in tally.ended:
                    warmphases.finish_phase(args.database, name, tally.result(name))

            def on_plan(index, plan):
                for line in tally.start(steps[index].key):
                    bars.say(line)
                report_ended()
                if steps[index].key in warmup.FAMILIES:
                    warmphases.start_phase(args.database, warmup.phase_of(steps[index].key))
                bars.plan(index, plan)

            def on_result(index, _unit, result):
                if steps[index].key in tally.totals:
                    tally.add(steps[index].key, result)
                bars.result(result)

            with bars:
                wrapper.run("Cache warm-up", steps, on_plan, on_result, database=args.database,
                            argv=sys.argv[1:], kind=runrecord.MAINTENANCE)
                for line in tally.finish():
                    bars.say(line)
            report_ended()
            warmphases.end(args.database)
        except (redisqueue.Unavailable, runrecord.Busy, units.RunStopped) as exc:
            warmphases.end(args.database, error=str(exc))
            raise Failed(str(exc)) from exc
        except BaseException as exc:
            warmphases.end(args.database, error=f"{type(exc).__name__}: {exc}" if str(exc) else type(exc).__name__)
            raise
    warmed = [step.key for step in steps if step.key in warmup.FAMILIES]
    log.normal(f"Cache warm-up done: {len(warmed)} families.")


def _warm_on_change(args):
    """
    `planetgen warm --on-change DATABASE` (what `warmstart.request` starts): waits for more changes, then warms once
    for each batch of changes (`warmstart.drain`). Does nothing when `DATABASE` is not the one the site shows. A
    failed run stops it, said in the log; the next change starts it again.
    """
    if args.on_change != args.database:
        log.normal(f"The change was on {args.on_change}, the site shows {args.database}: no warm-up.")
        warmstart.release(args.on_change)
        return

    def run(kinds):
        log.normal(f"Warming after a change ({', '.join(kinds)}).")
        args.index_mode = "build" if set(kinds) & set(warmstart.FULL_INDEX_TAGS) else "refresh"
        _warm_once(args, reason=', '.join(kinds))

    try:
        rounds = warmstart.drain(args.on_change, run)
    except Failed as exc:
        log.error(f"error: the cache warm-up after a change stopped: {exc}")
        raise SystemExit(1)
    log.normal(f"The warm-up after a change ran {rounds} time{'' if rounds == 1 else 's'}.")


def run_warm(args):
    """
    `planetgen warm`: runs the cache-warming commands the arguments ask for (`_warm_once`), or with `--on-change`
    what a change asks for (`_warm_on_change`). Exits 1 on any failure.
    """
    args.database = _site_database()
    _stop_on_sigterm()
    if args.on_change is not None:
        _warm_on_change(args)
        return
    try:
        _warm_once(args)
    except Failed as exc:
        log.error(f"error: the cache warm-up stopped: {exc}")
        raise SystemExit(1)
