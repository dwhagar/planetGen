# planetgen/cli/worker.py

"""
One RQ worker for planetGen's queues, at the lowered priority generation workers always run at.
Only the Generate page's job still starts one (a burst worker for its own queue) until OPS.244 moves that job to the ledger (OPS.235); `--idle-seconds N` is a worker that leaves after N seconds without work.

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
from rq.intermediate_queue import IntermediateQueue

from planetgen import _version
from planetgen.admin import updateguard
from planetgen.queue import allowlist, claims, memorygate, pause, progress_file, redisqueue, scheduler, stopall, work, workers
from planetgen.util import log
from planetgen.util.format import format_duration_seconds, format_megabytes

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
    (about 2.4 s each job). Nothing here opens a database or Redis connection, so the horses share it safely.
    OPS.271: a run's units and openings are jobs the supervisor's launcher forks from a template (`jobhost.WARM_SETS`, which holds the bright-star sampling tables of `units.warm_star_tables`), not jobs of a worker.

    OPS.229: the site's app, which a worker of the maintenance lane built once here (PERF.245), is the warm set of the maintenance lane's template now (`jobhost.WARM_SETS`, `warmup.warm_app`).

    Args:
        queues (iterable, optional): The queue names served (nothing here depends on them now).
    """
    import importlib
    for name in WARM_MODULES:
        importlib.import_module(name)



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
    it starts itself over when the release on disk has changed.

    PERF.152: a worker of a queue that `pause.holds` takes no job while the Queue page's pause flag is set. It waits
    before it asks for one and, if a job came to it in the moment the flag was set, holds that job (still queued,
    untouched) until the flag is cleared, so a paused queue starts no new unit. A unit already running finishes.

    PERF.254: every job passes the memory gate (`queue/memorygate.py`) before it is prepared: the worker takes the job's claim
    (`claims.claim_for`) and waits, with its heartbeat going, until the host has room for it; it gives the claim up when the job ends,
    renews it at every heartbeat while the job runs, and keeps the most the job held (the gate learns the size of each kind from it). A
    job waiting for room is still queued, and goes back to the front of its queue if the worker is asked to stop.

    OPS.229: a worker serves one queue. The maintenance lane's units (the cache warm-up and the checks) are forked jobs of the supervisor's launcher now, and the one-queue-per-worker rule is RQ's reliable dequeue
    back (a worker killed between taking a job and starting it no longer loses the job).
    """
    base = rq.Worker

    class WarmWorker(base):
        def execute_job(self, job, queue):
            try:
                super().execute_job(job, queue)
            finally:
                self._release_claim()

        def maintain_heartbeats(self, job):
            super().maintain_heartbeats(job)
            try:
                memorygate.refresh(self.connection, self.name)
            except Exception as exc:  # noqa: BLE001 -- a missed sample is a claim that expires a little later, never a failed job
                log.debug("The memory-gate claim of %s was not renewed: %s.", self.name, exc, area="queue")

        def _release_claim(self):
            try:
                memorygate.release(self.connection, self.name)
            except Exception as exc:  # noqa: BLE001 -- the claim expires by itself
                log.debug("The memory-gate claim of %s was not released: %s.", self.name, exc, area="queue")

        def _take_claim(self, job):
            """PERF.254: waits at the memory gate for the room this job claims, then holds the claim until the job ends (`execute_job`)."""
            try:
                claim = claims.claim_for(job, self.connection)
            except claims.NoClaim as exc:
                log.error(f"The claim of job {job.id} fell back to the default {format_megabytes(claims.DEFAULT_MB)}: {exc}")
                claim = claims.default_claim(job)
            admission = memorygate.acquire(self.connection, claim, self.name, job=job.id, cap_of=claims.cap_of, scheduler=scheduler.Scheduler(),
                                           on_wait=lambda decision, waited: self._while_gated(job, decision))

        def _while_gated(self, job, decision=None):
            """What the worker does at every poll of a wait at the gate: keeps its heartbeat, notes in the ledger why a unit waits (OPS.224: the gate's own words, `decision`), and gives the job back to the front of its queue
            and ends if everything was stopped or it was asked to stop, so a stopped worker never starts a job it was holding."""
            self.heartbeat()
            origin = next(queue for queue in self.queues if queue.name == job.origin)
            if stopall.is_stopped() or self._stop_requested:
                log.debug("The worker exits, asked to stop while %s waited at the memory gate.", job.id, area="queue")
                origin.push_job_id(job.id, at_front=True)
                self.register_death()
                os._exit(0)
            # A job a worker has taken and not started sits in the queue's intermediate list, and another worker's maintenance (every 10 minutes
            # and at its start) fails a job it has seen there for over a minute ("stuck in the intermediate queue"). A wait at the gate can
            # last longer, so each poll makes the job new to that check.
            self.connection.delete(IntermediateQueue(origin.key, self.connection).get_first_seen_key(job.id))

        def prepare_execution(self, job):
            self._take_claim(job)
            gc.collect()
            return super().prepare_execution(job)

        def _wait_while_paused(self):
            if not any(pause.holds(name) for name in self.queue_names()):
                return
            if not pause.is_paused(self.connection):
                return
            log.debug("The worker holds, because the queue is paused.", area="queue")
            while pause.is_paused(self.connection) and not self._stop_requested:
                self.heartbeat()
                time.sleep(pause.WAIT_SECONDS)

        def dequeue_job_and_maintain_ttl(self, *args, **kwargs):
            while True:
                if stale():
                    log.debug("The worker restarts itself, because the release on disk is not the one it loaded.", area="queue")
                    self.register_death()
                    os.execv(sys.executable, [sys.executable, "-m", "planetgen.cli.worker", *sys.argv[1:]])
                self._exit_if_stopped()
                self._wait_while_paused()
                result = super().dequeue_job_and_maintain_ttl(*args, **kwargs)
                if result is None:
                    return result
                self._exit_if_stopped()
                self._wait_while_paused()
                return result

        def _exit_if_stopped(self):
            """PERF.186: a worker that is still alive when everything is stopped ends before it takes or starts a job."""
            if stopall.is_stopped():
                log.debug("The worker exits, because everything was stopped.", area="queue")
                self.register_death()
                os._exit(0)

    return WarmWorker


def main(argv=None):
    parser = argparse.ArgumentParser(prog="planetgen.cli.worker", description=__doc__.strip().splitlines()[0])
    parser.add_argument("queues", nargs="*", default=[redisqueue.BULK_QUEUE],
                        help="the queues to serve, first one first (default: the bulk queue; a worker "
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
    log.apply_debug_arguments(args, parser, errors_to_stderr=True)
    try:
        updateguard.refuse_while_updating("A worker")     # OPS.53: the workers take no job during an update
    except updateguard.UpdateRunning as exc:
        log.debug("The worker was not started: %s", exc, area="queue")
        log.error(exc)
        return updateguard.EX_TEMPFAIL
    log.debug("The worker starts on queues %s (%s) with release %s.", args.queues, "burst" if args.burst else
              f"idle exit after {format_duration_seconds(args.idle_seconds)}", _version.__version__, area="queue")
    work.lower_priority(work.WORKER_NICENESS)
    warm(args.queues)
    connection = redisqueue.connect(args.url, blocking=True)   # RQ blocks on the server; keepalive finds it gone
    queues = [redisqueue.queue(name, connection) for name in args.queues]
    workers.reap_dead_jobs(connection)      # the jobs a killed worker left "running"
    worker = worker_class()(queues, connection=connection, name=args.name, serializer=redisqueue.JsonSerializer,
                            job_class=allowlist.GuardedJob,       # SEC.45: a job naming a function off the list is refused, not run
                            job_monitoring_interval=int(memorygate.POLL_SECONDS))     # PERF.254: the claim is renewed and sampled this often
    worker.work(burst=args.burst, with_scheduler=False,
                max_idle_time=None if args.burst else args.idle_seconds)
    log.debug("The worker exits because %s.", "its idle time ran out" if not args.burst else "the queue was empty", area="queue")


if __name__ == "__main__":
    sys.exit(main())
