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

import calendar
import os
import re
import secrets
import sys
import time

from planetgen.queue import redisqueue
from planetgen.util import draw

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
            error text, once it failed), `error_status` (the HTTP
            status of a refusal, else `None`) and `made` (ADM.31: for a
            finished job that makes sectors, `{"since", "until"}`, the
            Unix seconds window the Galaxy Map's `made` view takes, else
            `None`). `None` for an unknown or expired id, or a bad one.

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
    body = {"id": job_id, "state": state, "result": None, "error": None, "error_status": None, "made": None}
    if state == "succeeded":
        outcome = job.return_value()
        if "refused" in outcome:
            body.update(state="failed", error=outcome["refused"], error_status=outcome["status"])
        else:
            body["result"] = outcome["result"]
            body["made"] = _made_window(job)
    elif state == "failed":
        latest = job.latest_result()
        text = (latest.exc_string if latest is not None else None) or "the job ended without finishing"
        body["error"] = text.strip().splitlines()[-1]
    return body


MAKES_SECTORS = ("generate_neighborhood", "regenerate_sector")
"""tuple[str]: The queued functions that generate sectors (ADM.31)."""


def _made_window(job):
    """The `made` window of a finished job that generates sectors: from
    when it started to a second past when it ended; else `None`."""
    function = job.args[0] if job.args else None
    if getattr(function, "__name__", None) not in MAKES_SECTORS or job.started_at is None or job.ended_at is None:
        return None
    return {"since": calendar.timegm(job.started_at.timetuple()), "until": calendar.timegm(job.ended_at.timetuple()) + 1}


def settle_sectors(sector_ids, config):
    """The queued body of a settle (GEN.126): saves the sector paths of
    `sector_ids` and the sectors around them. Returns `{"saved": n}`."""
    from planetgen.db import sector_paths, store
    conn = store.get_connection(config or store.DEFAULT_MYSQL_CONFIG)
    try:
        return {"saved": sector_paths.settle_sectors(conn, sector_ids)}
    finally:
        conn.close()


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
    from planetgen.db import edits as editStore, sector_paths, store
    from planetgen.generation import run_galaxy
    conn = store.get_connection(config or store.DEFAULT_MYSQL_CONFIG)
    try:
        with conn:
            address = editStore.sector_address(conn, sector_id)
            around = [other for other in sector_paths.sectors_to_settle(conn, [sector_id]) if other != sector_id]
            counts = editStore.delete_sector_with_contents(conn, sector_id)
    finally:
        conn.close()
    with draw.bound(secrets.randbits(128)):
        result = run_galaxy.ensure_sector_generated(*address, config=config, settle=False)
    # GEN.126: the new sector's masses and the old neighbours' paths, now the neighbour set is final.
    conn = store.get_connection(config or store.DEFAULT_MYSQL_CONFIG)
    try:
        wanted = set(around)
        if result["sector_id"] is not None:
            wanted.add(result["sector_id"])
        sector_paths.settle_sectors(conn, wanted, expand=False)
    finally:
        conn.close()
    return {"deleted": counts, "sector_id": result["sector_id"], "sector_name": result["sector_name"]}


def run_command(argv, env, cwd, timeout, merge_stderr=False):
    """
    Runs a command to completion (the queued body of `command_and_wait`).

    Returns:
        dict: `returncode`, `stdout`, `stderr` (text; empty when merged
            into `stdout`), `timed_out`, and `error` (why it couldn't be
            started, else `None`).
    """
    import subprocess
    out = {"returncode": None, "stdout": "", "stderr": "", "timed_out": False, "error": None}
    try:
        done = subprocess.run(argv, env=env, cwd=cwd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT if merge_stderr else subprocess.PIPE, timeout=timeout,
                              check=False)
    except subprocess.TimeoutExpired as exc:
        out["timed_out"] = True
        out["stdout"] = (exc.stdout or b"").decode("utf-8", errors="replace")
    except OSError as exc:
        out["error"] = str(exc)
    else:
        out["returncode"] = done.returncode
        out["stdout"] = (done.stdout or b"").decode("utf-8", errors="replace")
        out["stderr"] = (done.stderr or b"").decode("utf-8", errors="replace")
    return out


def command_and_wait(argv, env, cwd, timeout, merge_stderr=False):
    """
    Runs a command on the queue and waits for it (the one-off system page
    and the Generate page's estimate, PERF.24). Without a Redis server
    (Windows without WSL's Redis) the command runs here instead, as the
    Generate page's jobs do.

    Returns:
        dict: `run_command`'s.
    """
    try:
        job_id = submit(run_command, argv, env, cwd, timeout, merge_stderr)
    except NoQueue:
        return run_command(argv, env, cwd, timeout, merge_stderr)
    job = wait(job_id, timeout + 15)
    if job is None or job["state"] != "succeeded":
        reason = job["error"] if job and job["state"] == "failed" else "the queue did not answer in time"
        return {"returncode": None, "stdout": "", "stderr": "", "timed_out": False,
                "error": f"the work queue could not run it: {reason}"}
    return job["result"]
