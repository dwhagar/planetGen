# planetgen/queue/redisqueue.py

"""
The work queue on Redis with RQ (PERF.24, step 1): the connection, the
queues, and the workers that run them.

Workers are started on demand, not installed as a service
(docs/design/work-queue-audit.md): a run queues its tasks on a queue of
its own and starts `planetgen.cli.worker --burst` processes for it
(`RQExecutor`), up to the run's worker count. A burst worker exits once
the queue is empty, so nothing runs while there's no work. The control
database's lease still keeps two runs' workers off the machine at once
(`planetgen.queue.work`).

Windows can't fork, so its workers are RQ's `SpawnWorker`
(`worker_class`), with Redis itself in WSL2 (OPS.27).
"""

import os
import pickle
import queue as queue_module
import subprocess
import sys
import time

import redis
import rq

from planetgen.util.appconfig import load_config

GENERATION_QUEUE = "planetgen-generation"
"""str: Sector fills, layer scatters and backfill blocks."""

WEB_QUEUE = "planetgen-web"
"""str: Generate page runs and work the API hands off."""

QUEUES = (GENERATION_QUEUE, WEB_QUEUE)
"""tuple: Every queue, in the order a worker serves them."""

URL_ENV_VAR = "PLANETGEN_REDIS_URL"
"""str: Overrides `config.json`'s `redis.url`."""


def redis_url(config=None):
    """The Redis server: `PLANETGEN_REDIS_URL`, else `config.json`'s
    `redis.url`."""
    url = os.environ.get(URL_ENV_VAR)
    if url:
        return url
    config = load_config() if config is None else config
    return config["redis"]["url"]


def connect(url=None):
    """A Redis connection to `url` (default `redis_url()`)."""
    return redis.Redis.from_url(url or redis_url())


def queue(name, connection):
    """The RQ queue `name` on `connection`."""
    return rq.Queue(name, connection=connection)


def worker_class():
    """RQ's forking `Worker`, or `SpawnWorker` where there's no fork
    (Windows)."""
    return rq.Worker if hasattr(os, "fork") else rq.SpawnWorker


def worker_argv(names, url, name=None, python=None):
    """The command line of one burst worker for queues `names`, called
    `name` when given, run by `python` (default this interpreter)."""
    return [python or sys.executable, "-m", "planetgen.cli.worker", "--burst", "--url", url,
            *(["--name", name] if name else []), *names]


class Unavailable(Exception):
    """No Redis server answers at the configured URL."""


class WorkerDied(RuntimeError):
    """A task's worker died under it (killed, out of memory) twice: RQ
    retried it once on a fresh worker, and that one died too."""


TASK_RETRIES = 1
"""int: How many times RQ runs a task again after its worker died."""

POLL_SECONDS = 0.05
"""float: How often `RQExecutor.wait` reads the tasks' states."""

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
    Redis list: the scatter's per-layer progress (PERF.4). It pickles as
    just its URL and key, so a task's payload can carry it to a burst
    worker (a multiprocessing manager's queue can't: the worker isn't
    the manager's child, so it can't authenticate to it).

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

    def _redis(self):
        if self._connection is None:
            self._connection = connect(self.url)
        return self._connection

    def put(self, item):
        """Appends `item` (anything picklable). A lost Redis is `OSError`."""
        try:
            pipe = self._redis().pipeline()
            pipe.rpush(self.key, pickle.dumps(item))
            pipe.expire(self.key, self.EXPIRE_SECONDS)
            pipe.execute()
        except redis.RedisError as exc:
            raise OSError(str(exc)) from exc

    def get(self, timeout):
        """The oldest item, waiting up to `timeout` seconds; `queue.Empty`
        when none comes, `OSError` when Redis is gone."""
        try:
            popped = self._redis().blpop([self.key], timeout=timeout)
        except redis.RedisError as exc:
            raise OSError(str(exc)) from exc
        if popped is None:
            raise queue_module.Empty
        return pickle.loads(popped[1])

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

    def _fill_workers(self):
        """Starts workers until `workers` of the run's own are alive."""
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

    def submit(self, *args):
        """Queues `task(*args)`; a task that can't be queued (an argument
        that won't pickle) comes back as an already failed future."""
        try:
            job = self.queue.enqueue(self.task, *args, job_timeout=-1, result_ttl=3600, failure_ttl=86400,
                                     retry=rq.Retry(max=TASK_RETRIES))
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
            time.sleep(POLL_SECONDS)

    def channel(self, name):
        """A `Channel` named `name` on this run's Redis, for its tasks to
        report through."""
        return Channel(self.url, f"{self.queue.name}:{name}")

    def shutdown(self, futures=()):
        """Cancels the tasks no worker has started, waits for the ones
        running, then for the workers to exit and drops the queue."""
        futures = list(futures)
        for future in futures:
            if future.job is not None and not future.done():
                try:
                    if future.job.get_status(refresh=True) in ("queued", "deferred", "scheduled"):
                        future.job.cancel()
                except Exception:  # noqa: BLE001 -- best effort; the queue is dropped below
                    pass
        while True:
            self._refresh(futures)
            if all(future.done() for future in futures):
                break
            if self.queue.count:
                self._fill_workers()
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
