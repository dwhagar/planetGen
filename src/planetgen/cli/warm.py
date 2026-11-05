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

A warm-up goes on without its terminal (PERF.235): the process that started it is only a client of a run whose units queue workers do, and
when its terminal is closed (an SSH session lost, a hangup) it stops drawing and writing to the terminal, says in the run's log that the
terminal was closed, and carries on to the end of the run; the admin page's phase list keeps moving. `planetgen warm --follow` watches a
running warm-up from any other session, drawing the same bars from the run record (the stages already done show as done); it starts nothing,
and ending it, with Ctrl+C or a hangup, leaves the warm-up as it was. A SIGTERM (the update guard's cancel, a stop-all, update.sh's time limit)
and Ctrl+C in the process that started the run still end it.

The answers go to Redis and to the disk cache the web server reads, so they are the ones a visitor gets; each is kept
only while the data it read is unchanged (docs/design/operations/cache-warm-up.md). A list of every family ends with the command
that opens the site after the sector-grid switch closed it (`maintenance.SECTOR_SWITCH`, UX.112): the site is closed
only while the maps are being rebuilt.
"""

import os
import signal
import sys
import time

from planetgen.admin import stop_reason, updateguard
from planetgen.cli import warm_progress
from planetgen.queue import redisqueue, runrecord, units, warmphases, warmstart, wrapper
from planetgen.util import log

RUN_TITLE = "Cache warm-up"
"""str: The title of the warm-up's run in the run record: what `--follow` looks for among the maintenance runs."""

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
    parser.add_argument("--follow", action="store_true",
                        help="watch the warm-up that is running (started anywhere, by update.sh or by hand) from its run record, "
                             "drawing the same progress; starts nothing, and ending it leaves the warm-up running")
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


TIME_LIMIT_ENV = "PLANETGEN_WARM_TIME_LIMIT_SECONDS"
"""str: update.sh's time limit for the stage (`WARM_STAGE_SECONDS`), told to the run so it can say the limit ended it."""

_started = time.monotonic()


class Stopped(KeyboardInterrupt):
    """The run was told to stop by a signal; `reason` says by what. A `KeyboardInterrupt`, so the run is marked canceled
    in Redis and its workers stop, as Ctrl+C does."""

    def __init__(self, reason, status):
        super().__init__(reason)
        self.reason = reason
        self.status = status


def _who_sent_sigterm():
    """What ended the run on SIGTERM: the reason a canceller left (`stop_reason`), else update.sh's time limit when the limit it
    told this run has passed, else that nothing said."""
    left = stop_reason.take(os.getpid())
    if left:
        return left
    try:
        limit = float(os.environ.get(TIME_LIMIT_ENV, ""))
    except ValueError:
        limit = 0.0
    if limit > 0 and time.monotonic() - _started >= limit - 5:
        return f"update.sh's time limit for this stage ({limit:g} s, WARM_STAGE_SECONDS)"
    return "a SIGTERM that left no reason (a kill, or a stop-all or cancel that could not write one)"


_hung_up = [False]
"""list[bool]: Whether the terminal this process ran in has been closed (`_detach_from_terminal`)."""


def _silence(fds=(1, 2)):
    """Points file descriptors `fds` (standard output and error) at /dev/null, so writing to a terminal that is gone (the bars, the log) fails
    nowhere."""
    try:
        null = os.open(os.devnull, os.O_RDWR)
    except OSError:
        return
    try:
        for fd in fds:
            try:
                os.dup2(null, fd)
            except OSError:
                pass
    finally:
        os.close(null)


def _detach_from_terminal():
    """What a hangup (SIGHUP: the terminal closed, an SSH session lost) does to the warm-up (PERF.235): it stops writing to the terminal and goes
    on. The run is the queue workers' and the process that started it only follows it, so ending the run for a lost terminal would waste what
    was built; the line goes to the run's log, which the admin page and `planetgen warm --follow` show."""
    if _hung_up[0]:
        return
    _hung_up[0] = True
    _silence()
    log.normal("The terminal this warm-up ran in was closed: it goes on without it. Watch it with `planetgen warm --follow`.")


def _stop_on_signals():
    """Makes SIGTERM (the update guard's cancel, a stop-all, update.sh's time limit) end the run as canceled in Redis, as Ctrl+C does, with the
    reason in a `Stopped` that `run_warm` prints as one line. A hangup (SIGHUP: the terminal closed) does not end it: the run goes on without the
    terminal (`_detach_from_terminal`)."""
    global _started
    _started = time.monotonic()
    _hung_up[0] = False
    stop_reason.clear(os.getpid())

    def on_term(_signum, _frame):
        raise Stopped(_who_sent_sigterm(), 143)

    signal.signal(signal.SIGTERM, on_term)
    signal.signal(signal.SIGHUP, lambda _signum, _frame: _detach_from_terminal())


def _say_stopped(reason, status):
    """The one line a stopped warm-up prints, and the exit status."""
    log.error(f"The cache warm-up was ended by {reason}. What it stored is kept: run `planetgen warm --progress` again and it "
              f"continues where it stopped.")
    raise SystemExit(status)


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
                rich_bars=warm_progress.stdout_is_terminal(), expected=warmup.EXPECTED_SECONDS)

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
                wrapper.run(RUN_TITLE, steps, on_plan, on_result, database=args.database,
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


def _follow(args):
    """
    `planetgen warm --follow` (PERF.235): draws the progress of the warm-up that is running, from its run record. It starts nothing and says
    nothing to the run, so ending it leaves the warm-up as it was.

    Raises:
        Failed: Redis did not answer, or the warm-up stopped before its end.
    """
    from planetgen.web import warmup
    try:
        conn = redisqueue.connect(redisqueue.redis_url())
        conn.ping()
    except Exception as exc:  # noqa: BLE001 -- any connection failure
        raise Failed(f"no Redis server answers: {exc}") from exc
    run_id = runrecord.active_run(conn, runrecord.MAINTENANCE)
    run = runrecord.read_run(conn, run_id) if run_id else None
    if run is None or run.get("title") != RUN_TITLE:
        log.normal("No cache warm-up is running.")
        return
    keys = [item["key"] for item in runrecord.read_steps(conn, run_id)]
    tally = warmup.PhaseTally(keys)
    bars = warm_progress.WarmProgress(keys, log.normal, show=True, rich_bars=warm_progress.stdout_is_terminal(),
                                      expected=warmup.EXPECTED_SECONDS)
    log.normal(f"Following the cache warm-up {run_id}: {len(keys)} steps. Ctrl+C stops watching; the warm-up goes on.")

    def on_plan(index, plan):
        for line in tally.start(keys[index]):
            bars.say(line)
        bars.plan(index, plan)

    def on_result(index, _unit, result):
        if keys[index] in tally.totals:
            tally.add(keys[index], result)
        bars.result(result)

    try:
        with bars:
            wrapper.follow(run_id, on_plan, on_result)
            for line in tally.finish():
                bars.say(line)
    except (redisqueue.Unavailable, units.RunStopped) as exc:
        raise Failed(str(exc)) from exc
    log.normal("Cache warm-up done.")


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
    what a change asks for (`_warm_on_change`), or with `--follow` only watches the one that is running (`_follow`). Exits 1 on any failure.
    A run told to stop (SIGTERM, Ctrl+C) prints one line saying what ended it and that a rerun continues, and exits 143 or 130, never with
    a traceback; a hangup does not stop it (`_detach_from_terminal`).
    """
    args.database = _site_database()
    if getattr(args, "follow", False):
        try:
            _follow(args)
        except Failed as exc:
            log.error(f"error: the cache warm-up stopped: {exc}")
            raise SystemExit(1)
        except KeyboardInterrupt:
            log.normal("Stopped watching. The warm-up goes on; `planetgen warm --follow` watches it again.")
            raise SystemExit(130)
        return
    _stop_on_signals()
    try:
        if args.on_change is not None:
            _warm_on_change(args)
            return
        try:
            _warm_once(args)
        except Failed as exc:
            log.error(f"error: the cache warm-up stopped: {exc}")
            raise SystemExit(1)
    except Stopped as stop:
        _say_stopped(stop.reason, stop.status)
    except KeyboardInterrupt:
        _say_stopped("Ctrl+C", 130)
