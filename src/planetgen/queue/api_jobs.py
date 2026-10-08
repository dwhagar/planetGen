# planetgen/queue/api_jobs.py

"""
Long API work on the Redis queue (PERF.24 step 4): a route that would
hold a web worker for minutes queues the work with `submit`, answers
`202 Accepted` with the job's id at once, and the client reads the
outcome from `GET /api/jobs/<id>` (`status`).

Each job has a queue of its own and one burst worker started for it,
detached from the web server (`redisqueue.start_detached`), as the
admin Generate page's jobs do (`planetgen.web.jobs`): a request timeout
or a graceful reload of the web server doesn't stop it, and the worker
exits when the job is done. The result is kept for `RESULT_TTL_SECONDS`.
"""

import os
import re
import secrets
import sys
import time

from planetgen.queue import redisqueue

QUEUE_PREFIX = "planetgen-api-"
"""str: A job's RQ queue and its worker are this plus its id."""

RESULT_TTL_SECONDS = 86400
"""int: How long a finished job's result or error is kept."""

JOB_ID_RE = re.compile(r"^[0-9a-f]{16}$")
"""re.Pattern: A job id. Checked before any id from a URL reaches Redis."""

_STATES = {
    "queued": "queued", "deferred": "queued", "scheduled": "queued",
    "started": "running", "finished": "succeeded",
    "failed": "failed", "stopped": "failed", "canceled": "failed",
}
"""dict: RQ's job states -> the API's."""

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
"""str: The planetGen checkout, where workers start."""


class NoQueue(Exception):
    """No Redis server answers, so nothing can be queued."""


def new_job_id():
    return secrets.token_hex(8)


SHORT_WAIT_SECONDS = 8.0
"""float: How long a route that queues quick work waits for the answer
before it gives up and answers `202` with the job's id."""


def execute(function, args):
    """
    What a worker runs for a queued job: `function(*args)`, with a refusal
    the work raises on purpose (an `ApiError`: 404, 409 ...) returned as
    data, so `status` can report the HTTP status it carries.
    """
    from planetgen.api.common import ApiError
    try:
        return {"result": function(*args)}
    except ApiError as exc:
        return {"refused": str(exc), "status": exc.status_code}


def submit(function, *args):
    """
    Queues `function(*args)` (a module-level function; its arguments and
    return value must pickle) and starts a burst worker for it.

    Returns:
        str: The job's id, for `status`.

    Raises:
        NoQueue: No Redis server answers at `redis.url`.
    """
    url = redisqueue.redis_url()
    connection = redisqueue.connect(url)
    try:
        connection.ping()
    except Exception as exc:  # noqa: BLE001 -- any connection failure
        raise NoQueue(f"no Redis server answers at {url} (config.json's redis.url); "
                      f"start Redis to run this: {exc}") from exc
    job_id = new_job_id()
    name = QUEUE_PREFIX + job_id
    redisqueue.queue(name, connection).enqueue(
        execute, args=(function, args), job_id=name, job_timeout=-1,
        result_ttl=RESULT_TTL_SECONDS, failure_ttl=RESULT_TTL_SECONDS,
    )
    redisqueue.start_detached(redisqueue.worker_argv([name], url, name=name, python=sys.executable), REPO_DIR)
    return job_id


def wait(job_id, seconds=SHORT_WAIT_SECONDS):
    """
    `status(job_id)` once the job has finished, or as it stands after
    `seconds`.
    """
    deadline = time.monotonic() + seconds
    while True:
        job = status(job_id)
        if job is None or job["state"] in ("succeeded", "failed") or time.monotonic() >= deadline:
            return job
        time.sleep(0.05)


def status(job_id):
    """
    One job's state.

    Returns:
        dict | None: `id`, `state` (`queued`, `running`, `succeeded` or
            `failed`), `result` (what the function returned, once it
            succeeded), `error` (the refusal, or the last line of the
            error text, once it failed) and `error_status` (the HTTP
            status of a refusal, else `None`). `None` for an unknown or expired id, or a bad one.

    Raises:
        NoQueue: No Redis server answers.
    """
    if not isinstance(job_id, str) or not JOB_ID_RE.match(job_id):
        return None
    try:
        connection = redisqueue.connect()
        job = redisqueue.fetch_job(QUEUE_PREFIX + job_id, connection)
    except redisqueue.Unavailable as exc:
        raise NoQueue(str(exc)) from exc
    if job is None:
        return None
    state = _STATES.get(job.get_status(refresh=True), "failed")
    body = {"id": job_id, "state": state, "result": None, "error": None, "error_status": None}
    if state == "succeeded":
        outcome = job.return_value()
        if "refused" in outcome:
            body.update(state="failed", error=outcome["refused"], error_status=outcome["status"])
        else:
            body["result"] = outcome["result"]
    elif state == "failed":
        latest = job.latest_result()
        text = (latest.exc_string if latest is not None else None) or "the job ended without finishing"
        body["error"] = text.strip().splitlines()[-1]
    return body


def generate_neighborhood(sector_id, radius_ly, config):
    """The queued body of `POST /api/sectors/<id>/generate-neighborhood`
    (a module-level function so a worker can import it by name)."""
    from planetgen.generation import run_galaxy
    return run_galaxy.generate_sector_neighborhood(sector_id, radius_ly=radius_ly, config=config)


def regenerate_sector(sector_id, config):
    """The queued body of `POST /api/sectors/<id>/regenerate`: deletes the
    galaxy-placed sector with everything in it, then generates its slot
    again.

    Returns:
        dict: `deleted` (the row counts), `sector_id` and `sector_name` of
            the new sector (`None` when the slot is outside the outline).
    """
    import random
    from planetgen.db import edits as editStore, store
    from planetgen.generation import run_galaxy
    conn = store.get_connection(config or store.DEFAULT_MYSQL_CONFIG)
    try:
        with conn:
            address = editStore.sector_address(conn, sector_id)
            counts = editStore.delete_sector_with_contents(conn, sector_id)
    finally:
        conn.close()
    random.seed()
    result = run_galaxy.ensure_sector_generated(*address, config=config)
    return {"deleted": counts, "sector_id": result["sector_id"], "sector_name": result["sector_name"]}
