# planetgen/cli/supervisor.py

"""
The supervisor process (OPS.207): one per host, started on demand by a submit, a status poll or the cache warm-up, and ending after an idle spell.

Usage:
    python -m planetgen.cli.supervisor [--url REDIS_URL] [--idle-seconds N]

It takes the host's lock first (`planetgen.queue.supervisor.take_lock`) and exits at once with status 0 when another supervisor holds it. For now it only
publishes that it is alive and logs its turns; the stages that follow give it the watch, the worker starts and the reconciler.
"""

import argparse
import os
import signal

from planetgen import _version
from planetgen.admin import updateguard
from planetgen.queue import redisqueue, supervisor
from planetgen.util import log


def main(argv=None):
    parser = argparse.ArgumentParser(prog="planetgen.cli.supervisor", description=__doc__.strip().splitlines()[0])
    parser.add_argument("--url", help="the Redis server (default: PLANETGEN_REDIS_URL, else config.json's redis.url)")
    parser.add_argument("--idle-seconds", type=float, default=supervisor.IDLE_SECONDS,
                        help=f"exit after this long with nothing to watch (default {supervisor.IDLE_SECONDS:.0f})")
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.set_component("supervisor")
    log.apply_debug_arguments(args, parser, errors_to_stderr=True)
    try:
        updateguard.refuse_while_updating("The supervisor")
    except updateguard.UpdateRunning as exc:
        log.debug("The supervisor was not started: %s", exc, area="queue")
        return updateguard.EX_TEMPFAIL
    lock = supervisor.take_lock()
    if lock is None:
        log.debug("The supervisor was not started: another one holds the host's lock.", area="queue")
        return 0
    stopping = []
    for number in (signal.SIGTERM, signal.SIGINT):
        signal.signal(number, lambda *_: stopping.append(True))
    log.debug("The supervisor %s starts with release %s.", os.getpid(), _version.__version__, area="queue")
    connection = redisqueue.connect(args.url)
    reason = supervisor.Supervisor(supervisor.SystemClock(), supervisor.RedisStore(connection), supervisor.ProcessTable(), duties=(supervisor.Reconciler(connection), supervisor.RunWatches(connection), supervisor.WorkerStarts(connection)),
                                   stop=lambda: bool(stopping) or updateguard.is_updating(), idle_seconds=args.idle_seconds).run()
    log.debug("The supervisor exits: %s.", reason, area="queue")
    lock.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
