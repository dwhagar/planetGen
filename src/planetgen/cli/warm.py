# planetgen.cli.warm

"""
`planetgen warm`: builds the site's caches ahead of the first visitor (MAP.134, PERF.104, PERF.127): the Galaxy Map's
opening view, every list, every sector, system and phenomenon page, and the search, so the first visit after a release
(or after a cache was cleared) isn't the one who waits for it. update.sh runs it in the background after every update;
run it by hand after clearing a cache. Run it as the web server's user so the files it writes are readable by Apache.

    planetgen warm [FAMILY ...] [--cap N]

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
only while the data it read is unchanged (docs/design/cache-warm-up.md). A list of every family ends with the command
that opens the site after the sector-grid switch closed it (`maintenance.SECTOR_SWITCH`, UX.112): the site is closed
only while the maps are being rebuilt.
"""

import signal
import sys
import time

from planetgen.admin import updateguard
from planetgen.db import dataversion
from planetgen.queue import redisqueue, runrecord, units, wrapper
from planetgen.util import log
from planetgen.util.settings import get_settings

UPDATE_WAIT_SECONDS = 1200
"""int: How long it waits for a running update to finish before it starts the run."""


MAX_PASSES = 3
"""int: Most times one run builds everything: once, and again while the data changed under it (PERF.141)."""

KIND = "warm-map"
"""str: What the warm-up registers as for the update guard (`updateguard.running`)."""


class Failed(Exception):
    """The warm-up could not run or a command of it failed; the message says what."""


def add_warm_arguments(parser):
    """Adds the `warm` subcommand's options: the families to warm and the cap."""
    parser.add_argument("families", nargs="*", metavar="FAMILY",
                        help="the cache families to warm, in order (default: every family; see 'planetgen commands' "
                             "for the warm-* commands)")
    parser.add_argument("--cap", type=int, default=None,
                        help="most sector, system and phenomenon pages to warm (default page_cache.warm_max_objects)")


def _wait_for_update():
    """Waits for a running update to finish: the queue takes no job during one (OPS.53)."""
    if not updateguard.is_updating():
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


def run_warm(args):
    """
    `planetgen warm`: builds the list of cache-warming commands the arguments ask for and runs it through the wrapper,
    printing each command's line. Exits 1 on any failure.
    """
    from planetgen.web import warmup
    args.database = _site_database()
    try:
        steps = wrapper.steps_for("warm", args)
    except ValueError as exc:
        raise SystemExit(f"warm: {exc}")
    _stop_on_sigterm()
    # Registered for the update guard (OPS.53), which lists it and stops it, the caches being safe to build again.
    with updateguard.running(KIND, "Site cache warm-up"):
        try:
            _wait_for_update()
            results = {}

            def on_result(index, _unit, result):
                results[steps[index].key] = result
                log.normal(result["line"])

            for number in range(1, MAX_PASSES + 1):
                changes = _change_count(args.database)
                wrapper.run("Cache warm-up", steps, lambda index, plan: None, on_result, database=args.database,
                            argv=sys.argv[1:], kind=runrecord.MAINTENANCE, workers_count=1)
                if number == MAX_PASSES or _change_count(args.database) == changes:
                    break
                log.normal("The data changed while the warm-up ran: building once more.")
        except (Failed, redisqueue.Unavailable, runrecord.Busy, units.RunStopped) as exc:
            log.error(f"error: the cache warm-up stopped: {exc}")
            raise SystemExit(1)
    warmed = [step.key for step in steps if step.key in warmup.FAMILIES]
    log.normal(f"Cache warm-up done: {len(warmed)} families.")


def _change_count(database):
    """How many data changes are recorded for `database`, or `None` when Redis does not say."""
    try:
        return dataversion.change_count(database)
    except dataversion.Unavailable:
        return None


def start_after_change(reason):
    """
    Starts `planetgen warm` in the background once a run has changed the galaxy (PERF.141): the end of a generation
    run and of the orbital update, so the caches are built again without anyone running a command. Does nothing when
    `cache.warm_after_change` is off, an update is running (it starts its own warm-up at its end), or no Redis
    server answers (the queue the warm-up runs on); a warm-up already running is left alone, it builds once more
    when it finds the data changed under it (`run_warm`).

    Args:
        reason (str): What changed, for the log line (`"generation"`, `"orbital update"`).

    Returns:
        bool: Whether a warm-up was started.
    """
    if not get_settings().cache.warm_after_change or updateguard.is_updating():
        return False
    if any(item["kind"] == KIND for item in updateguard.MarkerSource().items()):
        log.normal(f"A cache warm-up is already running: it builds once more for the {reason}.")
        return False
    try:
        redisqueue.connect().ping()
    except Exception as exc:    # no Redis, or one that refuses: the first visit builds what the warm-up would have
        log.normal(f"The cache warm-up did not start after the {reason}: {exc}")
        return False
    try:
        redisqueue.start_detached([sys.executable, "-m", "planetgen.cli.generate", "warm"], "/")
    except OSError as exc:
        log.normal(f"The cache warm-up did not start after the {reason}: {exc}")
        return False
    log.normal(f"The cache warm-up started in the background after the {reason}.")
    return True
