# planetgen/queue/workers.py

"""
The workers of the shared queues (PERF.95, folding in PERF.40, PERF.78 and ADM.15): who starts them, how
many there are, and how a running job is stopped. Section 5.4 of docs/design/operations/redis-and-the-work-queue.md.

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
import subprocess
import sys

import rq
from rq import command, worker_registration

from planetgen.admin import updateguard
from planetgen.queue import redisqueue, runrecord, stopall, work
from planetgen.util import log

WANTED_KEY = "planetgen:workers:wanted"
"""str: The admin's worker count for the bulk queue, a number in Redis (no expiry)."""

STARTING_PREFIX = "planetgen:workers:starting:"
STARTING_SECONDS = 90
"""The Redis counter, per queue, of workers started but not yet registered with RQ, which expires this long after
the last start: counted with the live ones so a burst of jobs doesn't start more than wanted. A worker on a host that
is short of memory or busy can take far longer than a few seconds to import and register, and a window shorter than
that starts a second worker for the same need, and a third (PERF.181)."""

DECIDING_KEY = "planetgen:workers:deciding"
DECIDING_SECONDS = 20
"""The Redis lock, one for the host, held by the caller that counts the workers and starts the missing ones. Callers come from
every thread of the web app (each map request that queues a tile), so several counted "none starting" at the same
moment and each started the whole cap: about 130 workers of 125 MB took a 16 GB host on 2026-10-24 and the kernel
killed MySQL (PERF.181). A caller that does not get the lock leaves the starting to the one that has it. The lock is
the host's, not the queue's, so the cap on all the host's workers (`host_worker_cap`) is checked and used up in one place."""

WORKER_COMMAND = "planetgen.cli.worker"
"""str: What every planetgen worker process has in its command line."""

WORKER_MEMORY_MB = 300
"""int: What one worker takes (its libraries, sampling tables and the work horse it forks), for sizing the host cap."""

HOST_MEMORY_SHARE = 0.6
"""float: The share of the host's memory the workers may take at most; the rest is for MySQL, Redis and the web server."""

BULK_WORKER_MEMORY_MB = 1700
"""int: What one worker takes at most in a stage that draws stars from the sampling tables but is not a sector stage (PERF.226): the
scatters. Such a worker holds the tables of every star population (about 1.4 GB with the datasets they read), and a worker and its
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

SECTOR_WORKER_MEMORY_MB = 2600
"""int: What one worker takes at most while it generates sectors (PERF.226). Measured on a populated database (a galaxy of one
shell, 4 workers), the worker's and its horse's memory together: 2,083 MB for the starting sector and 2,370 MB for a neighbourhood
sector, 10 percent more with the margin. The horse alone showed 1,487 MB (resident, which counts the worker's tables but not the
worker beside it), and a count by that figure let 10 workers in where 5 fit a 16 GB host."""

STAGE_MEMORY_SHARE = 0.5
"""float: The share of the host's memory a stage's workers may take (less than `HOST_MEMORY_SHARE`): MySQL's buffer pool, Redis,
the web server and the kernel keep the rest, so a stage that grows past its figure is slowed, not killed. 50 percent of a 16 GB
host in 700 MB workers is 11, in 1,700 MB ones 4, in 2,600 MB ones 3."""

STAGE_WORKER_MEMORY_MB = {
    "scatter-mass": SCATTER_WORKER_MEMORY_MB, "scatter-luminosity": SCATTER_WORKER_MEMORY_MB,
    "scatter-backfill": SCATTER_WORKER_MEMORY_MB, "scatter-band": SCATTER_WORKER_MEMORY_MB,
    "generate-start-sector": SECTOR_WORKER_MEMORY_MB, "generate-neighborhood": SECTOR_WORKER_MEMORY_MB,
    "generate-sectors": SECTOR_WORKER_MEMORY_MB, "generate-unplaced-sector": SECTOR_WORKER_MEMORY_MB,
    "scatter-phenomena": LEAN_WORKER_MEMORY_MB, "fill-map-luminosity": LEAN_WORKER_MEMORY_MB,
    "rescore-habitability": LEAN_WORKER_MEMORY_MB, "draw-nebula-field": LEAN_WORKER_MEMORY_MB,
    "link-neighbours": LEAN_WORKER_MEMORY_MB, "settle-sector-paths": LEAN_WORKER_MEMORY_MB,
    "warm-galaxy-nodes": LEAN_WORKER_MEMORY_MB, "warm-galaxy-deep-nodes": LEAN_WORKER_MEMORY_MB,
    "warm-sectors": LEAN_WORKER_MEMORY_MB, "warm-systems": LEAN_WORKER_MEMORY_MB, "warm-phenomena": LEAN_WORKER_MEMORY_MB,
    "warm-cluster-field": BULK_WORKER_MEMORY_MB,        # not measured yet (PERF.225 cuts it into parts): the larger figure
}
"""dict: The command kinds of the stages with units, and what one worker takes at most in MB (`stage_worker_cap`). Every other
command is one unit by nature (`CommandSpec.whole`) or inline, so one worker and nothing to limit."""

REFUSED_KEY = "planetgen:workers:refused"
REFUSED_NOTE_KEY = "planetgen:workers:refused-note"
REFUSED_NOTE_SECONDS = 60
"""The Redis counter of starts the host cap refused (kept an hour after the last), and the key that lets the log line
about it out once a minute instead of at every poll."""

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
    if stopall.is_stopped():
        return 0        # PERF.186: no worker starts while everything is stopped
    url = url or redisqueue.redis_url()
    connection = connection or redisqueue.connect(url)
    queue = redisqueue.queue(queue_name, connection)
    waiting = queue.count
    if not waiting:
        return 0
    if not connection.set(DECIDING_KEY, str(os.getpid()), nx=True, ex=DECIDING_SECONDS):
        return 0            # another caller is counting and starting on this host (PERF.181)
    try:
        need = _start_missing(queue, queue_name, connection, url, mysql_host, count, waiting)
    finally:
        connection.delete(DECIDING_KEY)
    if queue_name == redisqueue.MAINTENANCE_QUEUE:
        shed_bulk(connection, max(1, wanted(connection, mysql_host) - maintenance_share(connection, mysql_host)))
    return max(0, need)


def registered(connection, queue_names):
    """
    A worker has registered with RQ (`cli.worker`'s `register_birth`): it is live now, so it stops counting among the
    ones `ensure_workers` started and is waiting for. Without this the counter held it for the whole of
    `STARTING_SECONDS` on top of its place among the live ones, and a worker killed in that time was not replaced
    until the counter ran out (about 90 s). A worker not started by `ensure_workers` leaves no count to take it from.
    """
    for name in queue_names:
        key = STARTING_PREFIX + name
        if connection.exists(key) and connection.decr(key) <= 0:
            connection.delete(key)


def _start_missing(queue, queue_name, connection, url, mysql_host, count, waiting):
    """The counting and starting of `ensure_workers`, run by the one caller that holds the queue's lock. Returns how many
    workers it started."""
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
    log.trace(f"workers for {queue_name}: {waiting} job(s) waiting, cap {cap}, {live} live, {starting} starting; "
              f"{max(0, need)} to start", area="queue")
    started = 0
    for _ in range(max(0, need)):
        running, host_cap = host_worker_processes(), host_worker_cap()
        if running >= host_cap:
            _note_refused(connection, queue_name, running, host_cap)
            break
        connection.incr(starting_key)
        connection.expire(starting_key, STARTING_SECONDS)
        proc = redisqueue.start_detached(
            redisqueue.worker_argv([queue_name], url, python=sys.executable, burst=False), REPO_DIR)
        started += 1
        log.trace(f"workers for {queue_name}: started worker process {getattr(proc, 'pid', '?')}", area="queue")
    return started


def _physical_memory_bytes():
    """The host's memory in bytes, or `None` where the OS does not say."""
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None


def host_worker_cap():
    """
    The most planetgen worker processes this host may run at once, all lanes and commands together (PERF.185):
    `queue.max_host_workers` when set, else `HOST_MEMORY_SHARE` of the host's memory in workers of `WORKER_MEMORY_MB`
    each, no more than four per core, and never fewer than a core per bulk worker plus the other lanes (`cores + 4`).
    A 16 GB host with 8 cores gets 31.
    """
    from planetgen.util.settings import get_settings
    configured = get_settings().queue.max_host_workers
    if configured:
        return configured
    cores = work.cpu_count()
    floor = cores + 4
    total = _physical_memory_bytes()
    if total is None:
        return max(floor, 4 * cores)
    by_memory = int(total * HOST_MEMORY_SHARE / (WORKER_MEMORY_MB * 2 ** 20))
    return max(floor, min(by_memory, 4 * cores))


def stage_worker_mb(kind):
    """What one worker of a stage of command `kind` takes at most, MB (PERF.226): its figure in `STAGE_WORKER_MEMORY_MB`, or `None`
    for a command that is one unit or inline (one worker, nothing to limit)."""
    return STAGE_WORKER_MEMORY_MB.get(kind)


def stage_worker_cap(kind):
    """
    The most workers a stage of command `kind` may use at once (PERF.220, PERF.226): `STAGE_MEMORY_SHARE` of the host's memory in
    workers of what that stage's worker takes (`stage_worker_mb`), at least 1. `None` for a command that is one unit or inline or
    where the OS does not say how much memory there is. What a stage generates does not depend on how many workers share it
    (`work.task_seed`), so the cap only changes the time.
    """
    per_worker_mb = stage_worker_mb(kind)
    total = _physical_memory_bytes()
    if per_worker_mb is None or total is None:
        return None
    return max(1, int(total * STAGE_MEMORY_SHARE / (per_worker_mb * 2 ** 20)))


def count_worker_processes(listing):
    """
    How many worker processes `listing` (the output of `ps -axo pid=,ppid=,args=`) holds. A work horse is a fork of
    its worker with the same command line, so a worker's children do not count again.
    """
    rows = []
    for line in listing.splitlines():
        parts = line.split(None, 2)
        if len(parts) == 3 and parts[0].isdigit() and parts[1].isdigit() and WORKER_COMMAND in parts[2]:
            rows.append((int(parts[0]), int(parts[1])))
    pids = {pid for pid, _parent in rows}
    return sum(1 for _pid, parent in rows if parent not in pids)


def host_worker_processes():
    """How many planetgen worker processes run on this host now, started by anyone and registered or not (read from
    the process list, so a worker still importing its libraries counts too). 0 when the list cannot be read."""
    try:
        listing = subprocess.run(["ps", "-axww", "-o", "pid=,ppid=,args="], capture_output=True, text=True, timeout=10,
                                 check=False).stdout
    except (OSError, subprocess.SubprocessError):
        return 0
    return count_worker_processes(listing)


def _note_refused(connection, queue_name, running, host_cap):
    """Counts a start the host cap refused and says so in the log, once a minute."""
    connection.incr(REFUSED_KEY)
    connection.expire(REFUSED_KEY, 3600)
    if connection.set(REFUSED_NOTE_KEY, "1", nx=True, ex=REFUSED_NOTE_SECONDS):
        log.normal(f"Worker start refused for {queue_name}: {running} planetgen worker processes already run on this "
                   f"host, and the cap is {host_cap} (queue.max_host_workers, 0 sizes it from memory and cores).")


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
