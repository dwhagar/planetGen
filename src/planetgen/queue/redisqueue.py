# planetgen/queue/redisqueue.py

"""
The work queue on Redis with RQ (PERF.24): the connection, the queues, the
serializers and timeouts of their jobs, and the detached workers that run them.

Workers are started on demand, not installed as a service
(docs/design/operations/work-queue-audit.md): a run's units are RQ jobs
(`planetgen.queue.units`) and a run starts `planetgen.cli.worker` processes
for them (`planetgen.queue.workers`), up to the run's worker count.
"""

import datetime
import decimal
import json
import math
import subprocess
import sys
import threading

import redis
import rq

from planetgen.queue import redisclient
from planetgen.util.settings import get_settings

GENERATION_QUEUE = "planetgen-generation"
"""str: The queue of the first design (PERF.24), where the supervising process of a run put its tasks. Nothing is queued on
it since the rebuild: a run's units go to `BULK_QUEUE`. The name stays in `QUEUES` and in the worker's warm-up
(`cli/worker.warm`), which treat it as a generation queue."""

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
"""str: The units of a `MAINTENANCE` or `CHECK` run (the cache warm-up, the database check) and the node prefetch and
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

