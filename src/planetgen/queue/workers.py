# planetgen/queue/workers.py

"""
The workers of the two shared queues (PERF.95, folding in PERF.40, PERF.78 and ADM.15): who starts them, how
many there are, and how a running job is stopped. Section 5.4 of docs/design/redis-and-the-work-queue.md.

* Two queues, `redisqueue.BULK_QUEUE` and `redisqueue.INTERACTIVE_QUEUE`; a worker serves one of them, so a long
  run never holds up an admin's estimate and the interactive worker (PERF.40, PERF.78) is a separate, mostly idle,
  warm process, not a priority order.
* Workers are not burst workers: a run's queue is briefly empty at every stage barrier, so a worker waits for
  `IDLE_SECONDS` before it exits. Nothing runs when there is no work, and no service has to be installed.
* `ensure_workers` starts workers up to `wanted` (the admin's count in `WANTED_KEY`, else `work.worker_count()`),
  counts only live workers (a worker killed with `kill -9` leaves its registry entry behind, PERF.93), and is
  called by whoever queues work and by every status poll, so a dead worker is replaced.
* `set_wanted` (ADM.15) changes the count; lowering it tells the extra bulk workers to finish their unit and
  exit. "Ludicrous Speed" is a number (`ludicrous_count`), not a mode.
* `cancel_now` stops a job at once: a waiting one is cancelled, a running one gets RQ's stop command.
"""

import os
import socket
import sys

import rq
from rq import command, worker_registration

from planetgen.admin import updateguard
from planetgen.queue import redisqueue, work

WANTED_KEY = "planetgen:workers:wanted"
"""str: The admin's worker count for the bulk queue, a number in Redis (no expiry)."""

STARTING_PREFIX = "planetgen:workers:starting:"
STARTING_SECONDS = 15
"""The Redis counter, per queue, of workers started but not yet registered with RQ, which expires this long after
the last start: counted with the live ones so a burst of jobs doesn't start more than wanted."""

IDLE_SECONDS = 60
"""int: How long a worker waits for work before it exits (about 60 s covers every stage barrier)."""

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

    A job a killed worker was running is not touched here: RQ notices it when the worker's heartbeat runs out,
    about 90 s later, and runs it again (`Retry(max=1)`).
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


def ensure_workers(queue_name, connection=None, url=None, mysql_host=None):
    """
    Starts idle-timeout workers for the waiting jobs of the queue `queue_name`, never more than `wanted` in all
    (those serving the queue and alive, and those started a moment ago and not registered yet).

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
    need = min(waiting, wanted(connection, mysql_host)) - live_worker_count(queue) - starting
    for _ in range(max(0, need)):
        connection.incr(starting_key)
        connection.expire(starting_key, STARTING_SECONDS)
        redisqueue.start_detached(
            redisqueue.worker_argv([queue_name], url, python=sys.executable, burst=False), REPO_DIR)
    return max(0, need)


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
