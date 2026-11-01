# planetgen/cli/worker.py

"""
One RQ worker for planetGen's queues, at the lowered priority
generation workers always run at. `planetgen.queue.workers.ensure_workers`
starts these, each serving one queue and exiting after `--idle-seconds`
without work (a run's queue is empty for a moment at every stage barrier,
so a burst worker would leave too soon).

Usage:
    python -m planetgen.cli.worker [--burst | --idle-seconds N] [--name NAME] [--url REDIS_URL] [QUEUE ...]
"""

import argparse
import gc
import os
import re
import sys
import time

import rq

from planetgen import _version
from planetgen.admin import updateguard
from planetgen.queue import pause, progress_file, redisqueue, stopall, units, work, workers
from planetgen.util import log

WARM_MODULES = (
    "planetgen.generation.run_plan",
    "planetgen.generation.run_galaxy",
    "planetgen.generation.run_common",
    "planetgen.generation.bright_stars",
    "planetgen.queue.api_jobs",
    "planetgen.queue.units",
    "planetgen.web.job_runner",
)
"""tuple: The modules the queue's tasks run in. Imported before the worker
forks, so a work horse starts with them loaded instead of importing them
again for every job (PERF.42)."""


def warm(queues=None):
    """
    PERF.42: loads in the worker, before RQ forks a horse per job, what every
    job would otherwise rebuild in its own process: the generation modules
    (about 2.4 s each job) and, for the generation queue, the bright-star
    sampling tables (`star_population._bright_table`, one per luminosity
    floor and population) at the scatter's floor. Nothing
    here opens a database or Redis connection, so the horses share it safely.

    Args:
        queues (iterable, optional): The queue names served; the tables are
            built when the generation queue is among them (or none are
            given).
    """
    import importlib
    for name in WARM_MODULES:
        importlib.import_module(name)
    if queues and not {redisqueue.GENERATION_QUEUE, redisqueue.BULK_QUEUE} & set(queues):
        return
    from planetgen import tuning
    from planetgen.generation import star_population
    for population in star_population.POPULATIONS:
        star_population.warm_bins(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL, population)


def code_version_on_disk():
    """The release in `_version.py` as it is on disk now (the loaded module
    keeps the one this process started with)."""
    path = os.path.splitext(_version.__file__)[0] + ".py"
    try:
        with open(path, "r", encoding="utf-8") as handle:
            match = re.search(r'^__version__\s*=\s*"([^"]+)"', handle.read(), re.M)
    except OSError:
        return None
    return match.group(1) if match else None


def stale(loaded=None, on_disk=None):
    """Whether the code on disk is another release than this worker loaded."""
    on_disk = on_disk or code_version_on_disk()
    return on_disk is not None and on_disk != (loaded or _version.__version__)


def worker_class():
    """
    RQ's worker class with PERF.42's guard: a worker that pre-imports its
    code would keep serving the old release after an update (a horse used to
    import the new files on every job), so before it waits for the next job
    it starts itself over when the release on disk has changed. PERF.146: before it forks a horse for a unit, it
    works out what the unit would work out from the galaxy alone (`units.warm_unit`), so the horse starts with it
    instead of building it again. That is done before RQ creates the unit's execution record, not after: the record
    expires 90 s after its last heartbeat and nothing beats while the worker warms, so a warm of over a minute on a
    loaded host let it expire before the horse's first heartbeat, and the unit, finished, failed its run on
    `ValueError: Execution ... not found in Redis` when RQ cleaned up after it.

    PERF.152: a worker of a queue that `pause.holds` takes no job while the Queue page's pause flag is set. It waits
    before it asks for one and, if a job came to it in the moment the flag was set, holds that job (still queued,
    untouched) until the flag is cleared, so a paused queue starts no new unit. A unit already running finishes.

    TEST.371: when it registers, it stops counting among the workers `ensure_workers` is waiting for
    (`workers.registered`), so one killed soon after it started is replaced at once.
    """
    class WarmWorker(rq.Worker):
        def register_birth(self):
            super().register_birth()
            workers.registered(self.connection, self.queue_names())

        def prepare_execution(self, job):
            units.warm_unit(job)
            gc.collect()
            gc.freeze()      # PERF.226: the collector would write the header of every object the worker holds in each horse, copying
            return super().prepare_execution(job)           # the pages of its tables (a worker and its horse took 2.2 GB for 1.5 of tables)

        def _wait_while_paused(self):
            if not any(pause.holds(name) for name in self.queue_names()):
                return
            log.debug("worker holds: the queue is paused", area="queue")
            while pause.is_paused(self.connection) and not self._stop_requested:
                self.heartbeat()
                time.sleep(pause.WAIT_SECONDS)

        def dequeue_job_and_maintain_ttl(self, *args, **kwargs):
            if stale():
                log.debug("worker restarts itself: the release on disk is not the one it loaded", area="queue")
                self.register_death()
                os.execv(sys.executable, [sys.executable, "-m", "planetgen.cli.worker", *sys.argv[1:]])
            self._exit_if_stopped()
            self._wait_while_paused()
            result = super().dequeue_job_and_maintain_ttl(*args, **kwargs)
            if result is not None:
                self._exit_if_stopped()
                self._wait_while_paused()
            return result

        def _exit_if_stopped(self):
            """PERF.186: a worker that is still alive when everything is stopped ends before it takes or starts a job."""
            if stopall.is_stopped():
                log.debug("worker exits: everything was stopped", area="queue")
                self.register_death()
                os._exit(0)

    return WarmWorker


def main(argv=None):
    parser = argparse.ArgumentParser(prog="planetgen.cli.worker", description=__doc__.strip().splitlines()[0])
    parser.add_argument("queues", nargs="*", default=[redisqueue.INTERACTIVE_QUEUE],
                        help="the queues to serve, first one first (default: the interactive queue; a worker "
                             "serves one kind)")
    parser.add_argument("--burst", action="store_true", help="exit once the queues are empty")
    parser.add_argument("--idle-seconds", type=int, default=workers.IDLE_SECONDS,
                        help=f"exit after this long without work (default {workers.IDLE_SECONDS}; ignored with "
                             "--burst)")
    parser.add_argument("--name", help="the worker's name (default: RQ picks one)")
    parser.add_argument("--url", help="the Redis server (default: PLANETGEN_REDIS_URL, else config.json's redis.url)")
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.set_component("worker")
    os.environ.pop(progress_file.ENV_VAR, None)     # a worker serves many runs; the progress file is its starter's
    log.apply_debug_arguments(args, parser)
    try:
        updateguard.refuse_while_updating("A worker")     # OPS.53: the workers take no job during an update
    except updateguard.UpdateRunning as exc:
        log.debug("worker not started: %s", exc, area="queue")
        print(exc, file=sys.stderr)
        return updateguard.EX_TEMPFAIL
    log.debug("worker starts: queues %s, %s, release %s", args.queues, "burst" if args.burst else
              f"idle exit after {args.idle_seconds} s", _version.__version__, area="queue")
    work.lower_priority(workers.niceness(args.queues))
    warm(args.queues)
    serializers = {redisqueue.serializer_for(name) for name in args.queues}
    if len(serializers) > 1:
        parser.error("one worker serves queues of one kind: the web queues (JSON) or a run's own queue (pickle), "
                     "not both.")
    connection = redisqueue.connect(args.url, blocking=True)   # RQ blocks on the server; keepalive finds it gone
    queues = [redisqueue.queue(name, connection) for name in args.queues]
    workers.reap_dead_jobs(connection)      # the jobs a killed worker left "running"
    worker = worker_class()(queues, connection=connection, name=args.name, serializer=serializers.pop(),
                            work_horse_killed_handler=units.horse_killed)     # PERF.106: a killed unit fails its run now
    worker.work(burst=args.burst, with_scheduler=False,
                max_idle_time=None if args.burst else args.idle_seconds)
    log.debug("worker exits: no work for %s s", "its idle time" if not args.burst else "the queue was empty",
              area="queue")


if __name__ == "__main__":
    sys.exit(main())
