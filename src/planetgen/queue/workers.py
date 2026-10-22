# planetgen/queue/workers.py

"""
The workers of the shared queues (PERF.95, folding in PERF.40, PERF.78 and ADM.15): who starts them, how
many there are, and how a running job is stopped. Section 5.4 of docs/design/redis-and-the-work-queue.md.

* Lanes (PERF.151): `redisqueue.BULK_QUEUE` (a run's units), `redisqueue.INTERACTIVE_QUEUE` (short API jobs) and
  `redisqueue.LONG_QUEUE` (API jobs expected to run long); a worker serves one of them, so a long run never holds up
  an admin's estimate and the interactive worker (PERF.40, PERF.78) is a separate, mostly idle, warm process, not
  a priority order. Each lane has its own worker cap: `wanted` for bulk, `queue.interactive_workers` (outside
  `wanted`, at `INTERACTIVE_NICENESS`) and `queue.long_workers` for the other two.
* Workers are not burst workers: a run's queue is briefly empty at every stage barrier, so a worker waits for
  `IDLE_SECONDS` before it exits. Nothing runs when there is no work, and no service has to be installed.
* `ensure_workers` starts workers up to `wanted` (the admin's count in `WANTED_KEY`, else `work.worker_count()`),
  counts only live workers (a worker killed with `kill -9` leaves its registry entry behind, PERF.93), and is
  called by whoever queues work and by every status poll, so a dead worker is replaced.
* `set_wanted` (ADM.15) changes the count; lowering it tells the extra bulk workers to finish their unit and
  exit. "Ludicrous Speed" is a number (`ludicrous_count`), not a mode.
* `cancel_now` stops a job at once: a waiting one is cancelled, a running one gets RQ's stop command.
* `reap_dead_jobs` clears the jobs a killed worker left "running" in Redis: the web queues' jobs have no time
  limit, so RQ never gives them up.
"""

import os
import socket
import sys

import rq
from rq import command, worker_registration

from planetgen.admin import updateguard
from planetgen.queue import redisqueue, runrecord, work

WANTED_KEY = "planetgen:workers:wanted"
"""str: The admin's worker count for the bulk queue, a number in Redis (no expiry)."""

STARTING_PREFIX = "planetgen:workers:starting:"
STARTING_SECONDS = 15
"""The Redis counter, per queue, of workers started but not yet registered with RQ, which expires this long after
the last start: counted with the live ones so a burst of jobs doesn't start more than wanted."""

LENT_PREFIX = "planetgen:workers:lent:"
LENT_SECONDS = 3600
"""The Redis key, per bulk worker, of a worker asked to finish its unit and exit so a maintenance worker can have its
place (the lend rule, PERF.150); it expires this long after the request, so a worker that never got the command
is not skipped for ever."""

IDLE_SECONDS = 60
"""int: How long a worker waits for work before it exits (about 60 s covers every stage barrier)."""

INTERACTIVE_NICENESS = 5
"""int: How much `os.nice` lowers an interactive worker (the bulk and long workers use `work.WORKER_NICENESS`): the
short jobs get the CPU before a generation run's units."""

LUDICROUS_SHARE = 0.95
"""float: The share of the host's cores "Ludicrous Speed" aims at (ADM.15)."""

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
"""str: The planetGen checkout, where workers start."""


def ludicrous_count():
    """The worker count that aims at `LUDICROUS_SHARE` of the host's cores, whether or not MySQL runs here."""
    return max(1, int(LUDICROUS_SHARE * work.cpu_count()))


def wanted(connection, mysql_host=None):
    """The bulk worker count: the admin's (`set_wanted`), else `work.worker_count()`."""
    value = connection.get(WANTED_KEY)
    if value is not None:
        try:
            return max(1, int(value))
        except ValueError:
            pass
    return work.worker_count(mysql_host=mysql_host)


def lane_cap(queue_name):
    """
    The worker cap of a shared queue that has a lane of its own (PERF.151), from the settings: the interactive
    queue's reserve and the long queue's; `None` for any other queue, whose cap is `wanted`.
    """
    from planetgen.util.settings import get_settings
    section = get_settings().queue
    if queue_name == redisqueue.INTERACTIVE_QUEUE:
        return section.interactive_workers
    if queue_name == redisqueue.LONG_QUEUE:
        return section.long_workers
    return None


def niceness(queue_names):
    """How much a worker for `queue_names` lowers its priority: `INTERACTIVE_NICENESS` for the interactive queue,
    `work.WORKER_NICENESS` for every other."""
    return INTERACTIVE_NICENESS if redisqueue.INTERACTIVE_QUEUE in queue_names else work.WORKER_NICENESS


def process_exists(pid):
    """Whether a process with this id is running here (a process of another user counts as running)."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    return True


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


def live_worker_count(queue):
    """How many workers serve `queue` and are alive (`live_workers`)."""
    return len(live_workers(queue))


def ensure_workers(queue_name, connection=None, url=None, mysql_host=None, count=None):
    """
    Starts idle-timeout workers for the waiting jobs of the queue `queue_name`, never more than the queue's cap in
    all (those serving the queue and alive, and those started a moment ago and not registered yet). The cap is
    `count` when given (a run's own `--workers`, PERF.97, which leaves the admin's count alone), else the lane's
    (`lane_cap`, `maintenance_cap`), else `wanted`; the bulk queue's is never more than `wanted` less the maintenance
    lane's workers (`maintenance_share`), and an ask for the maintenance lane's workers sheds the bulk workers over
    that (`shed_bulk`).

    Returns:
        int: How many were started.
    """
    if updateguard.is_updating() and not updateguard.inside_update():
        return 0        # OPS.53: no worker starts during an update
    url = url or redisqueue.redis_url()
    connection = connection or redisqueue.connect(url)
    queue = redisqueue.queue(queue_name, connection)
    waiting = queue.count
    if not waiting:
        return 0
    starting_key = STARTING_PREFIX + queue_name
    starting = int(connection.get(starting_key) or 0)
    if queue_name == redisqueue.MAINTENANCE_QUEUE:
        cap = count or maintenance_cap(connection, mysql_host)
    else:
        cap = count or lane_cap(queue_name) or wanted(connection, mysql_host)
        if queue_name == redisqueue.BULK_QUEUE:
            cap = min(cap, max(1, wanted(connection, mysql_host) - maintenance_share(connection, mysql_host)))
    live = live_worker_count(queue)
    need = min(waiting, cap) - live - starting
    for _ in range(max(0, need)):
        connection.incr(starting_key)
        connection.expire(starting_key, STARTING_SECONDS)
        redisqueue.start_detached(
            redisqueue.worker_argv([queue_name], url, python=sys.executable, burst=False), REPO_DIR)
    if queue_name == redisqueue.MAINTENANCE_QUEUE:
        shed_bulk(connection, max(1, wanted(connection, mysql_host) - maintenance_share(connection, mysql_host)))
    return max(0, need)


def maintenance_cap(connection, mysql_host=None):
    """
    The most workers the maintenance lane may have (PERF.150): `queue.maintenance_workers` while a generation run is
    live, so upkeep never takes more than its guaranteed share from the run, and the whole of `wanted` when none is.
    """
    if runrecord.active_run(connection, runrecord.GENERATION):
        from planetgen.util.settings import get_settings
        return get_settings().queue.maintenance_workers
    return wanted(connection, mysql_host)


def maintenance_share(connection, mysql_host=None):
    """How many of the host's workers the maintenance lane has or needs: those busy, plus those its waiting jobs ask
    for, never more than `maintenance_cap`. (An idle maintenance worker waits for its 60 seconds to pass without
    using a core, so the bulk lane is back to `wanted` as soon as the lane has drained.) The bulk lane gives up that
    many places."""
    queue = redisqueue.queue(redisqueue.MAINTENANCE_QUEUE, connection)
    busy = sum(1 for worker in live_workers(queue) if worker.get_state() != "idle")
    return min(maintenance_cap(connection, mysql_host), busy + queue.count)


def shed_bulk(connection, cap):
    """
    The lend rule (PERF.150): asks the bulk workers over `cap` to finish their unit and exit (RQ's shutdown command),
    idle ones first and the newest first, so a maintenance worker can have a place without the host's total growing.
    A worker already asked is not asked again (`LENT_PREFIX`), and counts as gone. The bulk lane starts workers
    again by itself (`ensure_workers`) when the maintenance lane no longer needs the places.

    Returns:
        int: How many workers were asked.
    """
    live = live_workers(redisqueue.queue(redisqueue.BULK_QUEUE, connection))
    staying = [worker for worker in live if not connection.exists(LENT_PREFIX + worker.name)]
    over = len(staying) - cap
    if over <= 0:
        return 0
    staying.sort(key=lambda worker: (worker.get_state() != "idle", -(worker.birth_date.timestamp()
                                                                        if worker.birth_date else 0.0)))
    for worker in staying[:over]:
        connection.set(LENT_PREFIX + worker.name, "1", ex=LENT_SECONDS)
        command.send_shutdown_command(connection, worker.name)
    return over


def set_wanted(connection, count):
    """
    Sets the bulk worker count (ADM.15) and, when it went down, asks the extra bulk workers to finish their
    current unit and exit (RQ's shutdown command). Returns the count stored.

    Raises:
        ValueError: `count` is not a whole number of at least 1.
    """
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError("the worker count is a whole number of at least 1")
    connection.set(WANTED_KEY, count)
    extra = live_workers(redisqueue.queue(redisqueue.BULK_QUEUE, connection))[count:]
    for worker in extra:
        command.send_shutdown_command(connection, worker.name)
    return count


def cancel_now(job_id, connection=None, serializer=None):
    """
    Stops the job `job_id` at once (PERF.40): a waiting job is cancelled, a running one is sent RQ's stop command
    (`cancel()` alone leaves a running job running; safe because a unit's save is one transaction). A finished or
    unknown job is left alone.

    Returns:
        str | None: What was done, `"cancelled"` or `"stopped"`, or `None`.
    """
    connection = connection or redisqueue.connect()
    job = redisqueue.fetch_job(job_id, connection, serializer)
    if job is None:
        return None
    state = job.get_status(refresh=True)
    if state in ("queued", "deferred", "scheduled"):
        job.cancel()
        return "cancelled"
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
                if job.retries_left and job.retries_left > 0:
                    job.retry(queue, pipeline)
                else:
                    job.set_status(rq.job.JobStatus.FAILED, pipeline=pipeline)
                    job._handle_failure(f"Its worker {owner} was gone (killed or restarted) before the job ended.",
                                        pipeline, worker_name=owner)
                pipeline.execute()
            cleared.append((name, job_id))
    return cleared
