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
from rq.intermediate_queue import IntermediateQueue
from rq.worker import RoundRobinWorker

from planetgen import _version
from planetgen.admin import updateguard
from planetgen.queue import allowlist, claims, memorygate, pause, progress_file, redisqueue, runrecord, scheduler, stopall, units, work, workers
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
    (about 2.4 s each job). Nothing here opens a database or Redis connection, so the horses share it safely.
    PERF.226: the bright-star sampling tables (`star_population._bright_table`, about 1.4 GB with the datasets they read) are not
    built here any more: a worker built them at its start for every stage, and a worker of a stage that draws no star
    (a link, a settle, a re-score) took 1.5 GB of the host for nothing. A unit that draws stars builds them in the worker before
    its horse forks (`units.warm_unit`, `units.warm_star_tables`), and a worker that holds them starts itself over before a
    unit that does not (`WarmWorker.prepare_execution`).

    Args:
        queues (iterable, optional): The queue names served; a worker of `APP_QUEUES` also builds the site's app once (PERF.245).
    """
    import importlib
    for name in WARM_MODULES:
        importlib.import_module(name)
    if APP_QUEUES.intersection(queues or ()):
        from planetgen.web import warmup
        warmup.build_app()      # PERF.245: made once here and dropped; what it imports and sets up stays, so a horse's build_app takes 0.07 s, not 6.8 s


GATE_REPORT_SECONDS = 5.0
"""float: How long a unit's or an opening's wait at the memory gate must last for the run's log to say so."""

APP_QUEUES = frozenset({*redisqueue.maintenance_lane(), redisqueue.INTERACTIVE_QUEUE, redisqueue.LONG_QUEUE})
"""frozenset: The queues whose jobs build the site's app (`web/warmup.build_app`): a warm-up unit, a node prefetch, a slow page's answer
(`queue/pagejobs.build_answer`). Their workers build it once before they fork a horse for a job (PERF.245); building it opens no
connection or thread, so the horses share it safely. On a 4-core test box it cost a horse 6.8 s a job and costs a worker 75 MB."""


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


def worker_class(round_robin=False):
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

    OPS.145: the flag may be a solo hold (`pause.hold_solo`), for the one restarted unit that is to run alone. A worker then takes jobs as usual to reach that unit, which is at the
    front of its queue, and turns every other job it takes back to the queue and waits (`_turn_away_during_solo_hold`).

    PERF.254: every job passes the memory gate (`queue/memorygate.py`) before it is prepared: the worker takes the job's claim
    (`claims.claim_for`) and waits, with its heartbeat going, until the host has room for it; it gives the claim up when the job ends,
    renews it at every heartbeat while the job runs, and keeps the most the job held (the gate learns the size of each kind from it). A
    job waiting for room is still queued, and goes back to the front of its queue if the worker is asked to stop.

    TEST.371: when it registers, it stops counting among the workers `ensure_workers` is waiting for
    (`workers.registered`), so one killed soon after it started is replaced at once.

    PERF.245: a worker of several queues (`round_robin`, the maintenance lane's `redisqueue.maintenance_lane()`) is RQ's
    `RoundRobinWorker`: it takes one job from each queue in turn, so a cache warm-up that queued thousands of units never makes a
    database check or a visitor's node wait behind them; with one queue it is the plain `rq.Worker`.
    """
    base = RoundRobinWorker if round_robin else rq.Worker

    class WarmWorker(base):
        def register_birth(self):
            super().register_birth()
            workers.registered(self.connection, self.queue_names())

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
                log.error(f"The claim of job {job.id} fell back to the default {claims.DEFAULT_MB} MB: {exc}")
                claim = claims.default_claim(job)
            admission = memorygate.acquire(self.connection, claim, self.name, job=job.id, cap_of=claims.cap_of, scheduler=scheduler.Scheduler(),
                                           on_wait=lambda decision, waited: self._while_gated(job))
            if admission.waited >= GATE_REPORT_SECONDS:
                self._report_wait(job, claim, admission)

        def _while_gated(self, job):
            """What the worker does at every poll of a wait at the gate: keeps its heartbeat, and gives the job back to the front of its queue
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

        def _report_wait(self, job, claim, admission):
            """A job of a run that waited at the gate says so in the run's log: how long, for what, and why it was let in."""
            if job.func_name not in (units.UNIT_FUNCTION, units.OPEN_FUNCTION):
                return
            try:
                runrecord.log(self.connection, job.args[0], f"Memory gate: {job.description or job.id} waited {admission.waited:,.0f} s for "
                              f"{claim.mb:,.0f} MB ({admission.decision.reason}).", kind="gate.wait", level="warning" if admission.timed_out else "info")
            except Exception as exc:  # noqa: BLE001 -- a log line is never a reason to fail the job
                log.swallowed("cli/worker", "worker_class.WarmWorker._report_wait", exc)

        def prepare_execution(self, job):
            if units.sheds_tables_for(job):
                self._start_over_with(job)
            self._take_claim(job)
            units.warm_unit(job)
            gc.collect()
            gc.freeze()      # PERF.226: the collector would write the header of every object the worker holds in each horse, copying
            return super().prepare_execution(job)           # the pages of its tables (a worker and its horse took 2.2 GB for 1.5 of tables)

        def _start_over_with(self, job):
            """PERF.226: a worker that holds the star tables and is given a unit that needs none puts the unit back at the front of its
            queue and starts itself over, so the unit runs in a worker of its own size. The unit is only queued still (RQ starts it
            after this), so it is not lost."""
            log.debug("The worker restarts itself, because it holds the star tables and %s needs none.", job.id, area="queue")
            next(queue for queue in self.queues if queue.name == job.origin).push_job_id(job.id, at_front=True)
            self.register_death()
            os.execv(sys.executable, [sys.executable, "-m", "planetgen.cli.worker", *sys.argv[1:]])

        def _wait_while_paused(self, job_id=None):
            if not any(pause.holds(name) for name in self.queue_names()):
                return
            if not pause.blocks(self.connection, job_id):
                return
            log.debug("The worker holds, because the queue is paused.", area="queue")
            while pause.blocks(self.connection, job_id) and not self._stop_requested:
                self.heartbeat()
                time.sleep(pause.WAIT_SECONDS)

        def _turn_away_during_solo_hold(self, job):
            """
            OPS.145: a job taken while the queues are held for a restarted unit that runs alone, and that is not that unit, goes back to the back of its queue (the unit is at the front, and a job put in front of it
            could keep every worker from reaching it) and this worker holds off until the hold lifts, so nothing starts beside the unit. Returns whether it was turned away. (A pause by hand keeps the job in hand, as before.)
            """
            if not any(pause.holds(name) for name in self.queue_names()) or pause.solo_of(self.connection) is None or not pause.blocks(self.connection, job.id):
                return False
            log.debug("The worker turns %s away, because the queue is held for a unit that runs alone.", job.id, area="queue")
            next(queue for queue in self.queues if queue.name == job.origin).push_job_id(job.id)
            self._wait_while_paused(job.id)
            return True

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
                if self._turn_away_during_solo_hold(result[0]):
                    continue
                self._wait_while_paused(result[0].id)
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
    log.apply_debug_arguments(args, parser, errors_to_stderr=True)
    try:
        updateguard.refuse_while_updating("A worker")     # OPS.53: the workers take no job during an update
    except updateguard.UpdateRunning as exc:
        log.debug("The worker was not started: %s", exc, area="queue")
        log.error(exc)
        return updateguard.EX_TEMPFAIL
    log.debug("The worker starts on queues %s (%s) with release %s.", args.queues, "burst" if args.burst else
              f"idle exit after {args.idle_seconds} s", _version.__version__, area="queue")
    work.lower_priority(workers.niceness(args.queues))
    warm(args.queues)
    connection = redisqueue.connect(args.url, blocking=True)   # RQ blocks on the server; keepalive finds it gone
    queues = [redisqueue.queue(name, connection) for name in args.queues]
    workers.reap_dead_jobs(connection)      # the jobs a killed worker left "running"
    worker = worker_class(round_robin=len(queues) > 1)(queues, connection=connection, name=args.name, serializer=redisqueue.JsonSerializer,
                            job_class=allowlist.GuardedJob,       # SEC.45: a job naming a function off the list is refused, not run
                            work_horse_killed_handler=units.horse_killed,     # PERF.106: a killed unit fails its run now (PERF.254: after one retry)
                            job_monitoring_interval=int(memorygate.POLL_SECONDS))     # PERF.254: the claim is renewed and sampled this often
    worker.work(burst=args.burst, with_scheduler=False,
                max_idle_time=None if args.burst else args.idle_seconds)
    log.debug("The worker exits because %s.", "its idle time ran out" if not args.burst else "the queue was empty", area="queue")


if __name__ == "__main__":
    sys.exit(main())
