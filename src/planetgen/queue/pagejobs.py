# planetgen/queue/pagejobs.py

"""
Slow page requests as jobs (PERF.105; design in docs/design/operations/slow-page-answers.md).

A page that expects a slow answer does not build it inside the request. `ensure` queues one job
(`build_answer`, on the interactive queue like the cache warm-up's) that makes the request the way the page
would, which stores the answer in the tagged page cache (`web/lib/pagecache.py`) under the tags it read, so an
edit that touches the answer invalidates it exactly as it does any cached answer. The page gets the job's id at
once and polls `status` (`GET /api/page-jobs/<id>`) until the answer is there, then reads it from the cache.

Everything a job needs lives in Redis (PERF.86, Option A), all under `planetgen:slow...` and all kept
`tuning.SLOW_PAGE_JOB_SECONDS`:

- `planetgen:slowjob:<db>:<digest of the request>` -> the id of the job building that request. Every visitor
  who asks for the same request attaches to it, so one request is one job (it is claimed with `SET NX`).
- `planetgen:slowmeta:<id>` -> the job's record: what it builds, the estimated row count and seconds, when it
  was queued, started and finished, whether a cancel was asked, and the id of its queue job.
- `planetgen:slowtime:<db>:<digest>` -> how long the last build of the request took (`record_duration`), which
  the cost rule trusts over the row estimate.
- `planetgen:slowpath:<db>:<digest of the route>` -> the slowest build of any request on the route in the last
  `tuning.SLOW_PAGE_PATH_SECONDS`, which the cost rule uses for a request not timed itself (a search for a new
  name), before it falls back on the row estimate.

A job's progress is an estimate, not a count: the answer is one query, so there is no step to count. While the
job runs, `done` of `total` rows is the elapsed share of the expected time (never past 95 %), `eta_seconds` the
time left (`None` once the expected time has passed), and the status says so (`estimated`). It reaches `done ==
total` only when the answer is stored.

A job is also a node on the admin queue page (kind `slow-page`). An admin cancels it there
(`work.request_control`) or with `POST /api/page-jobs/<id>/cancel` (`cancel`): a job that has not started is
removed and one that has not begun its query stops; a query already running finishes, and its answer is kept in the
cache (it is a correct answer the next visitor gets), but the job reports `cancelled`. Every visitor waiting on a
request shares its job, so a visitor's own "Cancel" only stops that page waiting; it never cancels the job.
"""

import hashlib
import secrets
import time

import redis

from planetgen import tuning
from planetgen.queue import api_jobs, redisclient, work as workQueue
from planetgen.util import log

SOCKET_SECONDS = 0.5
"""float: How long one Redis call may take before it counts as unreachable."""

ID_RE = api_jobs.JOB_ID_RE
"""re.Pattern: A page job's id (16 hex digits). Checked before an id from a URL reaches Redis."""

STATES = ("queued", "running", "succeeded", "failed", "cancelled")
"""tuple: The states `status` reports."""

PROGRESS_CEILING = 0.95
"""float: The share of its expected time a running job's progress never passes, so it never reads as finished
before the answer is stored."""

CLAIM_WAIT_SECONDS = 2.0
"""float: How long a visitor who lost the race to queue a job waits for the winner's job to be recorded."""

MAX_TARGET_CHARS = 2000
"""int: The longest request a job is made for."""

NoQueue = api_jobs.NoQueue


def _redis():
    return redisclient.client(decode_responses=True, timeout=SOCKET_SECONDS, retries=1)


def _digest(db, target):
    return f"{db}:{hashlib.sha256(target.encode('utf-8')).hexdigest()[:32]}"


def job_key(db, target):
    """The Redis key holding the id of the job building `target` in `db`."""
    return f"planetgen:slowjob:{_digest(db, target)}"


def meta_key(job_id):
    """The Redis key of job `job_id`'s record."""
    return f"planetgen:slowmeta:{job_id}"


def time_key(db, target):
    """The Redis key holding how long `target` took to build last time."""
    return f"planetgen:slowtime:{_digest(db, target)}"


def path_key(db, path):
    """The Redis key holding the slowest recent build of any request on route `path` in `db`."""
    return f"planetgen:slowpath:{_digest(db, path)}"


def record_duration(db, target, seconds):
    """Remembers that building `target` in `db` took `seconds`, for the request and for its route (the slowest
    build within `SLOW_PAGE_PATH_SECONDS`). Never raises: the time is only a hint."""
    path = target.partition("?")[0]
    try:
        connection = _redis()
        connection.set(time_key(db, target), f"{seconds:.3f}", ex=tuning.SLOW_PAGE_TIMING_SECONDS)
        slowest = connection.get(path_key(db, path))
        if slowest is None or seconds >= float(slowest):
            connection.set(path_key(db, path), f"{seconds:.3f}", ex=tuning.SLOW_PAGE_PATH_SECONDS)
    except (redis.RedisError, OSError, ValueError):
        pass


def _duration(key):
    try:
        found = _redis().get(key)
    except (redis.RedisError, OSError):
        return None
    try:
        return None if found is None else float(found)
    except ValueError:
        return None


def last_duration(db, target):
    """How long the last build of `target` took, in seconds, or `None` when none was timed or Redis can't say."""
    return _duration(time_key(db, target))


def slowest_recent_duration(db, path):
    """The slowest build of any request on route `path` within `SLOW_PAGE_PATH_SECONDS`, or `None`."""
    return _duration(path_key(db, path))


def _read_meta(job_id):
    if not isinstance(job_id, str) or not ID_RE.match(job_id):
        return None
    try:
        return _redis().hgetall(meta_key(job_id)) or None
    except (redis.RedisError, OSError) as exc:
        raise NoQueue(f"no Redis server answers: {exc}") from exc


def _float(meta, field):
    try:
        return float(meta[field])
    except (KeyError, ValueError):
        return None


def status(job_id):
    """
    Where page job `job_id` stands.

    Returns:
        dict | None: `id`, `state` (one of `STATES`), `done` and `total` (rows, estimated; `None` when the request's
            size is not known), `unit` (`"rows"`), `estimated` (always true: see the module's text), `eta_seconds`
            (the expected time left, or `None`), `seconds` (how long it has run), `cancelling` (a cancel was asked
            and has not taken effect), `error` (once failed), `cache_key` (where the finished answer is stored)
            and `target` (the request). `None` for an unknown, expired or bad id.

    Raises:
        NoQueue: No Redis server answers.
    """
    meta = _read_meta(job_id)
    if meta is None:
        return None
    queued = api_jobs.status(meta["queue_job"]) if meta.get("queue_job") else None
    state, error, result = "queued", None, {}
    if queued is not None:
        state = queued["state"]
        if state == "succeeded":
            result = queued["result"] or {}
            state = "cancelled" if result.get("cancelled") else "succeeded"
        elif state == "failed":
            error = queued["error"]
    elif meta.get("queue_job") or time.time() - (_float(meta, "created") or 0) > api_jobs.RESULT_TTL_SECONDS:
        state, error = "failed", "the job is gone from the queue"
    if meta.get("cancelled") == "1" and state == "failed" and queued is not None and "ended without" in (error or ""):
        state, error = "cancelled", None        # removed from the queue before it started
    total = int(meta["rows"]) if meta.get("rows") else None
    expected = _float(meta, "expected") or 0.0
    started, finished = _float(meta, "started"), _float(meta, "finished")
    now = time.time()
    elapsed = 0.0 if started is None else (finished or now) - started
    done, eta = 0, None
    if state == "succeeded":
        done, eta = total, 0.0
    elif state == "running" and started is not None:
        share = min(PROGRESS_CEILING, elapsed / expected) if expected > 0 else 0.0
        done = None if total is None else int(total * share)
        eta = max(0.0, expected - elapsed) if elapsed < expected else None
    elif state == "queued":
        done = 0 if total is not None else None
    elif total is not None:
        done = 0
    from planetgen.web.lib import pagecache   # not at import: the queue modules load without the web package
    return {"id": job_id, "state": state, "done": done, "total": total, "unit": "rows", "estimated": True,
            "eta_seconds": eta, "seconds": round(elapsed, 3), "error": error,
            "cancelling": meta.get("cancelled") == "1" and state in ("queued", "running"),
            "cache_key": pagecache.shared_key(meta["db"], meta["target"]), "target": meta["target"]}


def ensure(db, target, rows=None, expected=None, config=None, retry=False):
    """
    The job that builds `target` in `db`: the one already queued or running for it, else a new one. A job that
    finished (its answer has since been invalidated or evicted) is replaced; one that failed or was cancelled is
    shown for `SLOW_PAGE_FAILED_SECONDS` and then replaced, or at once with `retry`.

    Args:
        db (str): The database the page shows.
        target (str): The request under `/api` with its query, as `apiclient` makes it.
        rows (int, optional): The estimated rows it reads, for the progress.
        expected (float, optional): The seconds it is expected to take.
        config (MySQLConfig, optional): The database the job connects to; the site's own when not given.
        retry (bool): Replace a failed or cancelled job at once.

    Returns:
        dict: `status` of the job.

    Raises:
        NoQueue: No Redis server answers, or an update is running (the caller builds the answer itself).
    """
    if len(target) > MAX_TARGET_CHARS:
        raise NoQueue("the request is too long to queue")
    try:
        connection = _redis()
        key = job_key(db, target)
        deadline = time.monotonic() + CLAIM_WAIT_SECONDS
        while True:
            current = connection.get(key)
            if current is not None:
                found = status(current)
                if found is not None and (found["state"] in ("queued", "running")
                                          or (found["state"] in ("failed", "cancelled") and not retry)):
                    if found["state"] in ("failed", "cancelled"):
                        connection.expire(key, tuning.SLOW_PAGE_FAILED_SECONDS)
                    return found
                connection.delete(key)       # finished, expired or retried: start afresh
            token = secrets.token_hex(8)
            now = time.time()
            connection.hset(meta_key(token), mapping={
                "db": db, "target": target, "created": f"{now:.3f}",
                "rows": "" if rows is None else int(rows), "expected": f"{float(expected or 0):.3f}"})
            connection.expire(meta_key(token), tuning.SLOW_PAGE_JOB_SECONDS)
            if connection.set(key, token, nx=True, ex=tuning.SLOW_PAGE_JOB_SECONDS):
                break
            connection.delete(meta_key(token))   # another visitor queued it between our read and our claim
            if time.monotonic() >= deadline:
                raise NoQueue("another request is queuing this one and did not finish")
        try:
            queue_job = api_jobs.submit(build_answer, db, target, token, config, expected_seconds=expected)
        except Exception:
            connection.delete(key, meta_key(token))
            raise
        connection.hset(meta_key(token), "queue_job", queue_job)
    except (redis.RedisError, OSError) as exc:
        raise NoQueue(f"no Redis server answers: {exc}") from exc
    return status(token)


def cancel(job_id):
    """
    Asks page job `job_id` to stop (see the module's text), and frees its request so the next visitor starts a new
    job.

    Returns:
        dict | None: `status` after the request, or `None` for an unknown id.

    Raises:
        NoQueue: No Redis server answers.
    """
    meta = _read_meta(job_id)
    if meta is None:
        return None
    try:
        connection = _redis()
        connection.hset(meta_key(job_id), "cancelled", "1")
        connection.expire(job_key(meta["db"], meta["target"]), tuning.SLOW_PAGE_FAILED_SECONDS)
    except (redis.RedisError, OSError) as exc:
        raise NoQueue(f"no Redis server answers: {exc}") from exc
    if meta.get("queue_job"):
        api_jobs.cancel(meta["queue_job"])
    return status(job_id)


def _cancel_asked(job_id):
    try:
        return _redis().hget(meta_key(job_id), "cancelled") == "1"
    except (redis.RedisError, OSError):
        return False


def _mark(job_id, **fields):
    try:
        _redis().hset(meta_key(job_id), mapping={name: f"{value:.3f}" for name, value in fields.items()})
    except (redis.RedisError, OSError) as exc:
        log.debug(f"Slow page job {job_id}: can't record {', '.join(fields)} ({exc}).")


def build_answer(db, target, job_id, config=None):
    """
    The queued body of a page job (a module-level function so a worker can import it): makes the request `target`
    through the site's own app, the way a visitor's page does, which stores the answer in the page cache.

    Args:
        db (str): The database the page shows.
        target (str): The request under `/api`, with its query.
        job_id (str): The page job's id (`ensure`).
        config (MySQLConfig, optional): The database to build the app on; the site's own when not given.

    Returns:
        dict: `seconds` and `bytes` of the answer, or `{"cancelled": True}`.
    """
    from planetgen.web import warmup
    from planetgen.web.lib import apiclient

    title = f"Slow page: /api{target[:150]}"
    node = workQueue.open_node("slow-page", title, config, database=db)
    state, started = "failed", time.monotonic()
    try:
        _mark(job_id, started=time.time())
        asked = _cancel_asked(job_id)
        if not asked and node.store.available:
            try:
                asked = bool(node.store.controls(node.chain).get("cancel"))
            except Exception as exc:  # noqa: BLE001 -- the admin's request is lost, the job goes on
                log.debug(f"Slow page job {job_id}: can't read the admin's requests ({exc}).")
        if asked:
            state = "cancelled"
            workQueue.record_task(node, "slow-page-line", f"{target[:200]}: cancelled before it started",
                                  state="done")
            return {"cancelled": True}
        app = warmup.build_app(config)
        with app.test_request_context("/"):
            size = len(apiclient.fetch_target(target))
        seconds = time.monotonic() - started
        workQueue.record_task(node, "slow-page-line", f"{target[:200]}: {size} bytes in {seconds:.1f} s",
                              seconds=seconds)
        state = "done"
        # A cancel that came while the query ran does not undo it: the answer is stored and correct.
        return {"cancelled": True, "seconds": seconds} if _cancel_asked(job_id) else {"seconds": seconds, "bytes": size}
    except Exception as exc:
        workQueue.record_task(node, "slow-page-line", f"{target[:200]}: {exc}", state="failed", error=str(exc))
        raise
    finally:
        _mark(job_id, finished=time.time())
        workQueue.close_node(node, state)
