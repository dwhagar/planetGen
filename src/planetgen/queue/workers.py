# planetgen/queue/workers.py

"""
The workers of the shared queues (PERF.95, folding in PERF.40, PERF.78 and ADM.15): who starts them, how
many there are, and how a running job is stopped. Section 5.4 of docs/design/operations/redis-and-the-work-queue.md.

* Lanes (PERF.151): `redisqueue.BULK_QUEUE` (a run's units) and the maintenance lane; a worker serves one of them. The site's short and long jobs are
  light jobs the supervisor forks (`queue/lightjobs.py`, OPS.230), not RQ jobs.
* Workers are not burst workers: a run's queue is briefly empty at every stage barrier, so a worker waits for
  `IDLE_SECONDS` before it exits. Nothing runs when there is no work, and no service has to be installed.
* OPS.235: nothing starts RQ workers for a run or the maintenance lane any more (`ensure_workers`, `start_workers` and the supervisor's `WorkerStarts` are gone: the launcher forks each unit); the
  Generate page's job still starts its own worker until OPS.244 moves it to the ledger.
* OPS.235 slice 1c: the admin's worker count (`wanted`, `set_wanted`, "Ludicrous Speed", ADM.15) is gone with the Queue page's Workers panel; the launcher sizes every stage by memory (`stage_worker_cap`) and the gate.
* `cancel_now` stops a job at once: a waiting one is canceled, a running one gets RQ's stop command.
* `reap_dead_jobs` clears the jobs a killed worker left "running" in Redis: the web queues' jobs have no time
  limit, so RQ never gives them up.
"""

import os
import socket

import rq
from rq import command, worker_registration

from planetgen.queue import redisqueue
from planetgen.util import log
from planetgen.util import process

WORKER_MEMORY_MB = 300
"""int: What one worker takes (its libraries, sampling tables and the work horse it forks), for sizing the host cap."""

HOST_MEMORY_SHARE = 0.6
"""float: The share of the host's memory the workers may take at most; the rest is for MySQL, Redis and the web server."""

BULK_WORKER_MEMORY_MB = 2100
"""int: What one worker takes at most in a stage that draws stars from the sampling tables but is not a sector stage (PERF.226): the
scatters. Raised from 1,700 to 2,100 MB (PERF.252): the bright stars of the full-galaxy run of 2026-11-08 held 1,879 MB, and 10 percent more is the margin. Such a worker holds the tables of every star population (about 1.4 GB with the datasets they read), and a worker and its
horse measured 1,456 to 1,564 MB together (`units.peak_resident_mb`, shared pages counted once) before the tables moved out of the
workers of the other stages. The kernel killed the run of 2026-10-31 (36 workers, 16,285 MB on a 15,992 MB host, MySQL first at
3:22 AM EDT) with the scatters sized by a guess of 300 MB."""

LEAN_WORKER_MEMORY_MB = 700
"""int: What one worker takes at most in a stage that draws no star (PERF.226): a worker starts without the sampling tables
(`cli/worker.warm`) and starts over to shed them before such a unit. Measured on a populated database, the worker and its horse
together: 601 MB for a link, 588 for a settle, 552 for a habitability re-score and 540 for a map-luminosity fill (1,456 to 1,564 MB
while every worker held the tables), 700 MB with the margin."""

SCATTER_WORKER_MEMORY_MB = BULK_WORKER_MEMORY_MB
"""int: What one worker takes at most while it scatters stars (PERF.220, PERF.226): a table worker's (`BULK_WORKER_MEMORY_MB`)."""

SECTOR_WORKER_MEMORY_MB = 3050
"""int: What one worker takes at most while it generates sectors (PERF.226). Measured on a populated database (a galaxy of one
shell, 4 workers), the worker's and its horse's memory together: 2,083 MB for the starting sector and 2,370 MB for a neighbourhood
sector, 10 percent more with the margin. The horse alone showed 1,487 MB (resident, which counts the worker's tables but not the
worker beside it), and a count by that figure let 10 workers in where 5 fit a 16 GB host."""

STAGE_MEMORY_SHARE = 0.5
"""float: The share of the host's memory a stage's workers may take (less than `HOST_MEMORY_SHARE`): MySQL's buffer pool, Redis,
the web server and the kernel keep the rest, so a stage that grows past its figure is slowed, not killed. 50 percent of a 16 GB
host in 700 MB workers is 11, in 2,100 MB ones 3, in 3,050 MB ones 2."""

STAGE_WORKER_MEMORY_MB = {
    "scatter-mass": SCATTER_WORKER_MEMORY_MB, "scatter-luminosity": SCATTER_WORKER_MEMORY_MB,
    "scatter-backfill": SCATTER_WORKER_MEMORY_MB, "scatter-band": SCATTER_WORKER_MEMORY_MB,
    "generate-start-sector": SECTOR_WORKER_MEMORY_MB, "generate-neighborhood": SECTOR_WORKER_MEMORY_MB,
    "generate-sectors": SECTOR_WORKER_MEMORY_MB, "generate-unplaced-sector": SECTOR_WORKER_MEMORY_MB,
    "scatter-phenomena": LEAN_WORKER_MEMORY_MB, "fill-map-luminosity": LEAN_WORKER_MEMORY_MB,
    "rescore-habitability": LEAN_WORKER_MEMORY_MB, "refill-effect-ledger": LEAN_WORKER_MEMORY_MB, "draw-nebula-field": LEAN_WORKER_MEMORY_MB, "empty-nebula-field": LEAN_WORKER_MEMORY_MB, "clear-bright-stars": LEAN_WORKER_MEMORY_MB, "store-scatter-limits": LEAN_WORKER_MEMORY_MB,
    "empty-phenomenon-scatter": LEAN_WORKER_MEMORY_MB, "drop-phenomenon-indexes": LEAN_WORKER_MEMORY_MB,
    "link-neighbors": LEAN_WORKER_MEMORY_MB, "settle-sector-paths": LEAN_WORKER_MEMORY_MB,
    "warm-galaxy-nodes": LEAN_WORKER_MEMORY_MB, "warm-galaxy-deep-nodes": LEAN_WORKER_MEMORY_MB,
    "warm-sectors": LEAN_WORKER_MEMORY_MB, "warm-systems": LEAN_WORKER_MEMORY_MB, "warm-phenomena": LEAN_WORKER_MEMORY_MB,
    "refresh-block-lights": LEAN_WORKER_MEMORY_MB,      # `one_unit` (no unit-length ceiling): the roll-up is the one unit (OPS.154 ended the uncharted index units of OPS.153)
    "warm-cluster-field": BULK_WORKER_MEMORY_MB,        # not measured yet (PERF.225 cuts it into parts): the larger figure
}
"""dict: The command kinds of the stages with units, and what one worker takes at most in MB (`stage_worker_cap`). Every other
command is one unit by nature (`CommandSpec.one_unit`) or inline, so one worker and nothing to limit. `one_unit` only exempts a command
from the unit-length ceiling: a one-unit command that is cut into several units has a figure too (the block lights: the run of
2026-11-08 held 15 workers at 15,128 MB of 15,992 MB with no figure and no cap)."""

IDLE_SECONDS = 60
"""int: How long a worker waits for work before it exits (about 60 s covers every stage barrier)."""

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
"""str: The planetGen checkout, where workers start."""


process_exists = process.pid_alive
"""Whether a process with this id is running here (`util.process.pid_alive`; a test patches this name)."""


def live_workers(queue):
    """
    The workers serving `queue` that are alive (PERF.93). `rq.Worker.count` is the size of a Redis set that
    only a worker's own start-up cleans, so a worker killed with `kill -9` (or by the out-of-memory killer)
    stayed in it, and nothing replaced it. This first drops the registry entries whose worker key is gone, then
    keeps a worker only if its process still exists on this host (a worker on another host is taken as alive: its
    key's time-to-live is all there is to go by), and removes the others from both registries, with their key.

    A job a killed worker was running is not touched here (`reap_dead_jobs` clears it).
    """
    worker_registration.clean_worker_registry(queue)
    here = socket.gethostname()
    alive = []
    for worker in rq.Worker.all(queue=queue):
        pid = getattr(worker, "pid", None)
        if worker.hostname and worker.hostname != here or not isinstance(pid, int) or process_exists(pid):
            alive.append(worker)
        else:
            queue.connection.srem(worker_registration.WORKERS_BY_QUEUE_KEY % queue.name, worker.key)
            queue.connection.srem(worker_registration.REDIS_WORKER_KEYS, worker.key)
            queue.connection.delete(worker.key)
    return alive


def _physical_memory_bytes():
    """The host's memory in bytes, or `None` where the OS does not say."""
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None


def stage_worker_mb(kind, conn=None, run_id=None):
    """What one worker of a stage of command `kind` takes at most, MB (PERF.226): its figure in `STAGE_WORKER_MEMORY_MB`, or `None`
    for a command that is one unit or inline (one worker, nothing to limit). With `conn` (OPS.157) the figure is the claim the gate makes for a unit of the kind
    (`profileclaim.claim`: from the job profile when this machine has history, held for the run `run_id`), which is the typed figure with no history."""
    typed = STAGE_WORKER_MEMORY_MB.get(kind)
    if conn is None or typed is None:
        return typed
    from planetgen.queue import profileclaim
    return profileclaim.claim(conn, f"unit:{kind}", typed, token=run_id) or typed


def stage_worker_cap(kind, conn=None, run_id=None):
    """
    The most workers a stage of command `kind` may use at once (PERF.220, PERF.226): `STAGE_MEMORY_SHARE` of the host's memory in
    workers of what that stage's worker takes (`stage_worker_mb`), at least 1. `None` for a command that is one unit or inline or
    where the OS does not say how much memory there is. What a stage generates does not depend on how many workers share it
    (`work.task_seed`), so the cap only changes the time. With `conn` the workers are sized by the claim from the job profile (OPS.157, see `stage_worker_mb`).
    """
    per_worker_mb = stage_worker_mb(kind, conn, run_id)
    total = _physical_memory_bytes()
    if per_worker_mb is None or total is None:
        return None
    return max(1, int(total * STAGE_MEMORY_SHARE / (per_worker_mb * 2 ** 20)))


def cancel_now(job_id, connection=None, serializer=None):
    """
    Stops the job `job_id` at once (PERF.40): a waiting job is canceled, a running one is sent RQ's stop command
    (`cancel()` alone leaves a running job running; safe because a unit's save is one transaction). A finished or
    unknown job is left alone.

    Returns:
        str | None: What was done, `"canceled"` or `"stopped"`, or `None`.
    """
    connection = connection or redisqueue.connect()
    job = redisqueue.fetch_job(job_id, connection, serializer)
    if job is None:
        return None
    state = job.get_status(refresh=True)
    if state in ("queued", "deferred", "scheduled"):
        job.cancel()
        return "canceled"
    if state == "started":
        command.send_stop_job_command(connection, job.id)
        return "stopped"
    return None


def _worker_gone(connection, name, here):
    """Whether the RQ worker `name` is dead: its key has expired, or it ran on this host and its process is gone
    (a worker on another host is taken as alive while its key lives)."""
    worker = rq.Worker.find_by_key(rq.Worker.redis_worker_namespace_prefix + name, connection=connection)
    if worker is None:
        return True
    pid = getattr(worker, "pid", None)
    return bool(worker.hostname == here and isinstance(pid, int) and not process_exists(pid))


def reap_dead_jobs(connection=None):
    """
    Clears the jobs that a worker that is gone left "running" in Redis. A web job (the Generate page's, the
    API's) is queued with no time limit, so RQ gives its entry in the started registry no expiry, and a worker
    killed with `kill -9`, by the out-of-memory killer or by a restart never gives it up: the job stayed "running"
    for ever, the update guard listed it, could not stop it and refused to go on. A job whose record has expired is
    dropped from the registry. A job whose worker is gone (`_worker_gone`) is retried when it has a retry left (the
    generation queue's units) and otherwise fails with a message saying why.

    Args:
        connection: The Redis client (default the configured server).

    Returns:
        list: `(queue name, job id)` of each job cleared.
    """
    connection = connection or redisqueue.connect()
    here = socket.gethostname()
    cleared = []
    for name in sorted(queue.name for queue in rq.Queue.all(connection=connection)):
        queue = redisqueue.queue(name, connection)
        registry = queue.started_job_registry
        entries = {}
        for raw in connection.zrange(registry.key, 0, -1):
            entry = raw.decode() if isinstance(raw, bytes) else raw
            entries.setdefault(registry.parse_job_id(entry), []).append(entry)
        for job_id, keys in entries.items():
            job = redisqueue.fetch_job(job_id, connection, redisqueue.serializer_for(name))
            if job is None:
                connection.zrem(registry.key, *keys)
                cleared.append((name, job_id))
                continue
            owner = job.worker_name
            if job.get_status(refresh=True) != rq.job.JobStatus.STARTED or not owner \
                    or not _worker_gone(connection, owner, here):
                continue
            with connection.pipeline() as pipeline:
                pipeline.zrem(registry.key, *keys)
                log.get("queue").error(f"Worker {owner} was gone (killed or restarted) before job {job_id} on queue {name} ended; "
                                       + ("the job is queued again." if job.retries_left and job.retries_left > 0 else "the job failed."),
                                       worker=owner)
                if job.retries_left and job.retries_left > 0:
                    job.retry(queue, pipeline)
                else:
                    job.set_status(rq.job.JobStatus.FAILED, pipeline=pipeline)
                    job._handle_failure(f"Its worker {owner} was gone (killed or restarted) before the job ended.",
                                        pipeline, worker_name=owner)
                pipeline.execute()
            cleared.append((name, job_id))
    return cleared
