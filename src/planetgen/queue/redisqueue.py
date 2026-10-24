# planetgen/queue/redisqueue.py

"""
The work queue on Redis with RQ (PERF.24, step 1): the connection, the
queues, and the workers that run them.

Workers are started on demand, not installed as a service
(docs/design/operations/work-queue-audit.md): a run queues its tasks on a queue of
its own and starts `planetgen.cli.worker --burst` processes for it
(`RQExecutor`), up to the run's worker count. A burst worker exits once
the queue is empty, so nothing runs while there's no work. The control
database's lease still keeps two runs' workers off the machine at once
(`planetgen.queue.work`).
"""

import datetime
import decimal
import json
import math
import queue as queue_module
import subprocess
import sys
import threading
import time

import redis
import rq

from planetgen.queue import redisclient
from planetgen.util.settings import get_settings

GENERATION_QUEUE = "planetgen-generation"
"""str: Sector fills, shell scatters and backfill blocks."""

BULK_QUEUE = "planetgen-bulk"
"""str: A run's units of work (PERF.95): the scatters, sector fills and backfill blocks of the work queue
rebuild. Its workers serve this queue and no other."""

INTERACTIVE_QUEUE = "planetgen-interactive"
"""str: Work an admin or the API starts and waits for (PERF.40, PERF.78): neighborhoods, estimates, exports, the
deep database check. Its workers serve this queue and no other, so a long run never holds one up."""

LONG_QUEUE = "planetgen-long"
"""str: The API jobs expected to take longer than `queue.long_job_seconds` (PERF.151): a neighbourhood or a large
settle, which can run for hours. Its own workers serve it, so a long job never holds the workers kept for the short
ones on `INTERACTIVE_QUEUE`."""

MAINTENANCE_QUEUE = "planetgen-maintenance"
"""str: The units of a `MAINTENANCE` or `CHECK` run (the cache warm-up, the database check) and the tile prefetch and
rebuild jobs (PERF.150). Its workers have a cap of their own (`workers.maintenance_cap`) and share the host's worker
count with the bulk workers, so a generation stage that queued thousands of units on `BULK_QUEUE` never makes this
work wait for them."""

JSON_QUEUE_PREFIXES = (BULK_QUEUE, MAINTENANCE_QUEUE, INTERACTIVE_QUEUE, LONG_QUEUE, "planetgen-web-", "planetgen-api-")
"""tuple: A queue whose name starts with one of these carries JSON (`serializer_for`): the shared queues, each
Generate page job's own (`planetgen-web-<id>`), and the API queues the tests give a name of their own."""

QUEUES = (GENERATION_QUEUE, BULK_QUEUE, MAINTENANCE_QUEUE, INTERACTIVE_QUEUE, LONG_QUEUE)
"""tuple: Every shared queue."""

def redis_url():
    """The Redis server: `PLANETGEN_REDIS_URL`, else `config.json`'s
    `redis.url` (the settings model applies both)."""
    return get_settings().redis.url


def connect(url=None, **options):
    """The shared client for `url` (default `redis_url()`); `options` are `redisclient.client`'s (PERF.94)."""
    return redisclient.client(url or redis_url(), **options)


def _json_default(value):
    """What `JsonSerializer` writes for a value `json` has no form for: a `Decimal` as the number it holds, a
    time as ISO text."""
    if isinstance(value, decimal.Decimal):
        return int(value) if value == value.to_integral_value() else float(value)
    if isinstance(value, (datetime.datetime, datetime.date)):
        return value.isoformat()
    raise TypeError(f"{type(value).__name__} can't be sent to the queue as JSON")


class JsonSerializer:
    """RQ's serializer for the queues that carry plain data (PERF.94, finding F8): a job's arguments and its
    result travel as JSON, not pickle, so a process that can reach Redis can queue data but not code to run."""

    @staticmethod
    def dumps(value, *args, **kwargs):
        return json.dumps(value, default=_json_default, *args, **kwargs).encode("utf-8")

    @staticmethod
    def loads(raw, *args, **kwargs):
        return json.loads(raw.decode("utf-8"), *args, **kwargs)


def serializer_for(name):
    """
    The serializer the queue `name` uses: JSON for the web queues (the Generate page's jobs and the API's,
    whose arguments are paths, ids and numbers); RQ's default, pickle, for a run's own queue, whose tasks carry
    Python objects until the work queue rebuild (PERF.95 to 98) turns them into ids. One worker serves one kind.
    """
    return JsonSerializer if name.startswith(JSON_QUEUE_PREFIXES) else None


TIMEOUT_MULTIPLE = 10
"""int: How many times its expected seconds a job may run before RQ ends it (PERF.106, PERF.144): work is expected
to take that long on all its workers, so a job that takes ten times as long is hung or starved, not slow."""

MIN_TIMEOUT_SECONDS = 15 * 60
"""int: The shortest limit a job gets, so work that is usually quick is not ended by one slow start-up."""

MAX_TIMEOUT_SECONDS = 6 * 3600
"""int: The longest limit a job gets from `job_timeout`'s defaults, and the one of a job whose work has no
estimate at all (a fixed ceiling, not an estimate); a site names a different floor, ceiling or `unestimated`
limit where that fits its work."""


def job_timeout(expected_seconds=None, ceiling=MAX_TIMEOUT_SECONDS, floor=MIN_TIMEOUT_SECONDS, unestimated=None):
    """
    How long a queued job may run before RQ ends it (PERF.144): `TIMEOUT_MULTIPLE` times the seconds the system
    expects it to take (recorded timings, then the stage ETAs, then the size of the work: each site's own best
    estimate), but at least `floor` and at most `ceiling`. Never -1 (no limit).

    Args:
        expected_seconds (float, optional): The best estimate of the job's run time; `None`, zero, negative or not
            finite counts as no estimate.
        ceiling (int): The most seconds.
        floor (int): The fewest seconds.
        unestimated (int, optional): The limit of a job with no estimate, a fixed ceiling and not an estimate
            (default `ceiling`).

    Returns:
        int: Whole seconds, at least 1.
    """
    if expected_seconds is None or not math.isfinite(expected_seconds) or expected_seconds <= 0:
        return max(int(ceiling if unestimated is None else unestimated), 1)
    return max(int(min(max(TIMEOUT_MULTIPLE * expected_seconds, floor), ceiling)), 1)


class TimedQueue(rq.Queue):
    """
    The RQ queue every site uses (`queue`): it refuses a job that has no finite time limit (PERF.144), so RQ always
    has a deadline to expire the entry of a job whose worker was killed, and a new kind of job cannot skip
    `job_timeout`. Every way RQ queues a job (`enqueue`, `enqueue_call`, `enqueue_many`, `enqueue_at`) goes through
    `create_job`.
    """

    def create_job(self, func, *args, timeout=None, **kwargs):
        """
        Raises:
            ValueError: `timeout` is missing (RQ would then use its default), zero or negative (-1 is RQ's "no
                limit").
        """
        if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or timeout <= 0:
            raise ValueError(f"a job queued to {self.name} needs a time limit in seconds (redisqueue.job_timeout), "
                             f"not {timeout!r}")
        return super().create_job(func, *args, timeout=timeout, **kwargs)


def queue(name, connection):
    """The RQ queue `name` on `connection`, with its serializer (`serializer_for`); jobs on it need a finite
    time limit (`TimedQueue`)."""
    return TimedQueue(name, connection=connection, serializer=serializer_for(name))


def fetch_job(job_id, connection=None, serializer=None):
    """The RQ job `job_id` (of a queue using `serializer`, default pickle), or `None` when there is none (never
    queued, or expired).

    Raises:
        Unavailable: No Redis server answers.
    """
    connection = connection or connect()
    try:
        return rq.job.Job.fetch(job_id, connection=connection, serializer=serializer)
    except rq.exceptions.NoSuchJobError:
        return None
    except redis.RedisError as exc:
        raise Unavailable(f"no Redis server answers: {exc}") from exc


def worker_argv(names, url, name=None, python=None, burst=True):
    """The command line of one worker for queues `names`, called `name` when given, run by `python` (default
    this interpreter). A burst worker exits once its queues are empty; with `burst=False` it waits for work
    and exits after its idle time (`planetgen.queue.workers.IDLE_SECONDS`)."""
    return [python or sys.executable, "-m", "planetgen.cli.worker", *(["--burst"] if burst else []), "--url", url,
            *(["--name", name] if name else []), *names]


def start_detached(argv, cwd):
    """
    Starts `argv` apart from this process, in its own session, and reaps
    it when it exits, so it doesn't linger as a zombie of the web server
    process; if this process goes first, init adopts it instead.

    Returns:
        subprocess.Popen: The started process.
    """
    proc = subprocess.Popen(
        argv, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        cwd=cwd, close_fds=True, start_new_session=True,
    )
    threading.Thread(target=proc.wait, name=f"planetgen-worker-{proc.pid}", daemon=True).start()
    return proc


class Unavailable(Exception):
    """No Redis server answers at the configured URL."""


class WorkerDied(RuntimeError):
    """A task's worker died under it (killed, out of memory) twice: RQ
    retried it once on a fresh worker, and that one died too."""


TASK_RETRIES = 1
"""int: How many times RQ runs a task again after its worker died."""

POLL_SECONDS = 0.05
"""float: How often `RQExecutor.wait` reads the tasks' states."""

WORKER_CRASH_LIMIT = 10
"""int: How many of a run's workers may exit with a failure in a row, with no task finished in between, before the run
fails. A worker that cannot start (an import error, no memory, a Redis it cannot reach) is started again at once while
tasks wait, and without a limit that goes on for ever, a process every few milliseconds (TEST.295)."""

STOP_GRACE_SECONDS = 1.5
"""float: How long a stopped worker gets between the two stop signals. RQ's first SIGTERM is a warm shutdown (the
work horse finishes its task), the second a cold one (the horse is killed); RQ ignores a second signal that comes
within a second of the first."""

STOP_WAIT_SECONDS = 30
"""float: How long `RQExecutor.stop` waits for the workers to exit after the second signal before it kills them."""

_DONE_STATES = ("finished", "failed", "stopped", "canceled")


class RQFuture:
    """One queued task, as `RQExecutor.submit` returns it: the RQ job, or
    the error that kept it from being queued (`error`)."""

    def __init__(self, job=None, error=None):
        self.job = job
        self.error = error
        self.status = "failed" if error is not None else "queued"

    def done(self):
        return self.status in _DONE_STATES

    def result(self):
        """`(result, seconds)` from the worker; raises the task's own
        exception, or `WorkerDied`."""
        if self.error is not None:
            raise self.error
        if self.status == "finished":
            outcome = self.job.return_value()
            if outcome[0] == "error":
                raise outcome[1]
            return outcome[1], outcome[2]
        latest = self.job.latest_result()
        text = (latest.exc_string if latest is not None else None) or f"the task's job ended {self.status}"
        raise WorkerDied(text.strip().splitlines()[-1])


class Channel:
    """
    A one-way message channel from a run's workers back to the run, on a
    Redis list: the scatter's per-shell progress (PERF.4). Items are plain
    data, sent as JSON (PERF.94). It pickles as just its URL and key, so a
    task's payload can carry it to a burst worker (a multiprocessing
    manager's queue can't: the worker isn't the manager's child, so it
    can't authenticate to it).

    Args:
        url (str): The Redis server.
        key (str): The list's key; `close` deletes it.
    """

    EXPIRE_SECONDS = 86400
    """int: A list a crashed run leaves behind goes after a day."""

    def __init__(self, url, key):
        self.url = url
        self.key = key
        self._connection = None

    def __getstate__(self):
        return {"url": self.url, "key": self.key}

    def __setstate__(self, state):
        self.__init__(state["url"], state["key"])

    def _redis(self, blocking=False):
        if self._connection is None:
            self._connection = connect(self.url)
        return connect(self.url, blocking=True) if blocking else self._connection

    def put(self, item):
        """Appends `item` (numbers, strings, booleans, lists; a tuple arrives as a list). A lost Redis is
        `OSError`."""
        try:
            pipe = self._redis().pipeline()
            pipe.rpush(self.key, json.dumps(item).encode("utf-8"))
            pipe.expire(self.key, self.EXPIRE_SECONDS)
            pipe.execute()
        except redis.RedisError as exc:
            raise OSError(str(exc)) from exc

    def get(self, timeout):
        """The oldest item, waiting up to `timeout` seconds; `queue.Empty`
        when none comes, `OSError` when Redis is gone."""
        try:
            popped = self._redis(blocking=True).blpop([self.key], timeout=timeout)
        except redis.RedisError as exc:
            raise OSError(str(exc)) from exc
        if popped is None:
            raise queue_module.Empty
        return json.loads(popped[1].decode("utf-8"))

    def close(self):
        """Deletes the list."""
        try:
            self._redis().delete(self.key)
        except redis.RedisError:
            pass


class RQExecutor:
    """
    Runs one run's tasks on RQ (PERF.24): a queue of the run's own and
    up to `workers` burst workers it starts for it, which inherit this
    process's environment and exit once the queue is empty. Stands in for
    the process pool `work.WorkQueue` used before, with the same
    `submit`/`wait`/`shutdown` shape.

    Args:
        workers (int): Worker processes.
        task (callable): The module-level function every job runs,
            `task(*args)`.
        url (str, optional): The Redis server (default `redis_url()`).
        name (str): The queue's name, unique to the run.

    Raises:
        Unavailable: No Redis server answers.
    """

    def __init__(self, workers, task, name, url=None):
        self.workers = max(1, int(workers))
        self.task = task
        self.url = url or redis_url()
        self.connection = connect(self.url)
        try:
            self.connection.ping()
        except redis.RedisError as exc:
            raise Unavailable(f"no Redis server answers at {self.url}: {exc}") from exc
        self.queue = queue(name, self.connection)
        self._procs = {}  # worker name -> its process
        self._started = 0
        self._counted = set()  # names of the workers whose failed exit is in `_crash_streak`
        self._crash_streak = 0
        self._broken = None  # the `WorkerDied` that ends the run once `WORKER_CRASH_LIMIT` workers failed in a row

    def _fill_workers(self):
        """Starts workers until `workers` of the run's own are alive, unless `WORKER_CRASH_LIMIT` of them in a row
        exited with a failure (then `_broken` says so and none is started)."""
        for name, proc in self._procs.items():
            code = proc.poll()
            if code not in (None, 0) and name not in self._counted:
                self._counted.add(name)
                self._crash_streak += 1
                last = (name, code)
                if self._crash_streak >= WORKER_CRASH_LIMIT and self._broken is None:
                    self._broken = WorkerDied(
                        f"{self._crash_streak} workers in a row exited with a failure and no task finished in "
                        f"between (the last, {last[0]}, with code {last[1]}); their output is discarded, so run "
                        f"`{' '.join(worker_argv([self.queue.name], self.url, last[0]))}` to see why")
        if self._broken is not None:
            return
        alive = sum(1 for proc in self._procs.values() if proc.poll() is None)
        for _ in range(self.workers - alive):
            self._started += 1
            name = f"{self.queue.name}-{self._started}"
            self._procs[name] = subprocess.Popen(
                worker_argv([self.queue.name], self.url, name),
                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def _worker_gone(self, name):
        """Whether `name` is one of the run's workers and has exited."""
        proc = self._procs.get(name)
        return proc is not None and proc.poll() is not None

    def submit(self, *args, expected_seconds=None):
        """Queues `task(*args)`, to be ended by RQ after `job_timeout(expected_seconds)`; a task that can't be
        queued (an argument that won't pickle) comes back as an already failed future."""
        try:
            job = self.queue.enqueue(self.task, *args, job_timeout=job_timeout(expected_seconds), result_ttl=3600,
                                     failure_ttl=86400, retry=rq.Retry(max=TASK_RETRIES))
        except Exception as exc:  # noqa: BLE001 -- the task's failure, reported like any other
            return RQFuture(error=exc)
        self._fill_workers()
        return RQFuture(job)

    def _refresh(self, futures):
        pending = [future for future in futures if not future.done()]
        if not pending:
            return
        jobs = rq.job.Job.fetch_many([future.job.id for future in pending], connection=self.connection)
        for future, job in zip(pending, jobs):
            if job is None:
                future.status = "canceled"
                continue
            status = job.get_status(refresh=False)
            future.status = getattr(status, "value", status)
            if future.status == "finished":
                self._crash_streak = 0
            if future.status == "failed" and job.retries_left:
                future.status = "queued"  # RQ puts it back for another worker
            elif future.status == "started" and self._worker_gone(job.worker_name):
                # The whole worker died, not just its work horse: nothing is
                # left to report the job, so it fails here.
                future.status = "failed"
                future.error = WorkerDied(f"worker {job.worker_name} died while running the task")
            future.job = job

    def wait(self, futures, block):
        """The futures in `futures` that are done, waiting for at least
        one when `block`."""
        futures = list(futures)
        while True:
            self._refresh(futures)
            done = [future for future in futures if future.done()]
            if done or not block or not futures:
                return done
            if self.queue.count:
                self._fill_workers()
            if self._broken is not None:
                raise self._broken
            time.sleep(POLL_SECONDS)

    def channel(self, name):
        """A `Channel` named `name` on this run's Redis, for its tasks to
        report through."""
        return Channel(self.url, f"{self.queue.name}:{name}")

    def stop(self):
        """
        Stops the run's workers and the tasks they are running (OPS.59): SIGTERM to each worker, then, after the
        grace RQ needs to tell a second signal from a repeat, SIGTERM again, which makes RQ take its work horse
        down with it; a worker still alive after `STOP_WAIT_SECONDS` is killed. Every process is reaped. A task
        that was running ends with its database transaction rolled back, as when a worker dies.
        """
        live = [proc for proc in self._procs.values() if proc.poll() is None]
        try:
            for proc in live:
                proc.terminate()
            grace = time.monotonic() + STOP_GRACE_SECONDS
            while time.monotonic() < grace and any(proc.poll() is None for proc in live):
                time.sleep(POLL_SECONDS)
            for proc in live:
                if proc.poll() is None:
                    proc.terminate()
            for proc in live:
                try:
                    proc.wait(timeout=STOP_WAIT_SECONDS)
                except subprocess.TimeoutExpired:
                    pass
        finally:
            for proc in live:
                if proc.poll() is None:
                    proc.kill()
                    proc.wait()

    def shutdown(self, futures=(), stop=False):
        """Cancels the tasks no worker has started, waits for the ones
        running (or, with `stop`, stops the workers and the tasks with
        them), then for the workers to exit and drops the queue."""
        futures = list(futures)
        for future in futures:
            if future.job is not None and not future.done():
                try:
                    if future.job.get_status(refresh=True) in ("queued", "deferred", "scheduled"):
                        future.job.cancel()
                except Exception:  # noqa: BLE001 -- best effort; the queue is dropped below
                    pass
        if stop:
            self.stop()
        while not stop:
            self._refresh(futures)
            if all(future.done() for future in futures):
                break
            if self.queue.count:
                self._fill_workers()
            if self._broken is not None:
                break
            time.sleep(POLL_SECONDS)
        for proc in self._procs.values():
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                proc.terminate()
        self._procs = {}
        try:
            self.queue.delete(delete_jobs=True)
        except Exception:  # noqa: BLE001
            pass
