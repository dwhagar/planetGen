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
from planetgen.queue import redisqueue, runrecord, units, wrapper
from planetgen.util import log

UPDATE_WAIT_SECONDS = 1200
"""int: How long it waits for a running update to finish before it starts the run."""


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
    with updateguard.running("warm-map", "Site cache warm-up"):
        try:
            _wait_for_update()
            results = {}

            def on_result(index, _unit, result):
                results[steps[index].key] = result
                log.normal(result["line"])

            wrapper.run("Cache warm-up", steps, lambda index, plan: None, on_result, database=args.database,
                        argv=sys.argv[1:], kind=runrecord.MAINTENANCE, workers_count=1)
        except (Failed, redisqueue.Unavailable, runrecord.Busy, units.RunStopped) as exc:
            log.error(f"error: the cache warm-up stopped: {exc}")
            raise SystemExit(1)
    warmed = [step.key for step in steps if step.key in warmup.FAMILIES]
    log.normal(f"Cache warm-up done: {len(warmed)} families.")
