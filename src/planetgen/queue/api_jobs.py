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


def submit(function, *args, **kwargs):
    """
    Queues `function(*args, **kwargs)` (a module-level function, given
    as its dotted name or the function itself; its arguments and return
    value must pickle) and starts a burst worker for it.

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
        function, args=args, kwargs=kwargs, job_id=name, job_timeout=-1,
        result_ttl=RESULT_TTL_SECONDS, failure_ttl=RESULT_TTL_SECONDS,
    )
    redisqueue.start_detached(redisqueue.worker_argv([name], url, name=name, python=sys.executable), REPO_DIR)
    return job_id


def status(job_id):
    """
    One job's state.

    Returns:
        dict | None: `id`, `state` (`queued`, `running`, `succeeded` or
            `failed`), `result` (what the function returned, once it
            succeeded) and `error` (its last line of error text, once it
            failed). `None` for an unknown or expired id, or a bad one.

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
    body = {"id": job_id, "state": state, "result": None, "error": None}
    if state == "succeeded":
        body["result"] = job.return_value()
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
