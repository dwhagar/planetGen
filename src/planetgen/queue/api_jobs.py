# planetgen/queue/api_jobs.py

"""
Long API work on the Redis queue (PERF.24 step 4): a route that would
hold a web worker for minutes queues the work with `submit`, answers
`202 Accepted` with the job's id at once, and the client reads the
outcome from `GET /api/jobs/<id>` (`status`).

The jobs share two queues, the interactive one (`QUEUE_NAME`) for short jobs and
the long one (`LONG_QUEUE_NAME`) for jobs expected to run longer than
`queue.long_job_seconds` (PERF.151), each served by workers started
detached from the web server (`workers.ensure_workers`), as the admin
Generate page's jobs do (`planetgen.web.jobs`): a request timeout or a
graceful reload of the web server doesn't stop them, and a worker exits
after `workers.IDLE_SECONDS` without work. A job starts a worker only while
fewer than its lane's cap serve the queue (`queue.interactive_workers`,
`queue.long_workers`), so a batch of uploads doesn't start a process (and its
imports) per job (PERF.39) and a long job never holds a worker kept for short
ones. The result is kept for `RESULT_TTL_SECONDS`.
"""

import calendar
import importlib
import math
import os
import re
import secrets
import sys
import time


from planetgen.admin import updateguard
from planetgen.queue import redisqueue, workers

QUEUE_PREFIX = "planetgen-api-"
"""str: A job's id in RQ is this plus its own id."""

QUEUE_NAME = redisqueue.INTERACTIVE_QUEUE
"""str: The queue a short job goes on, and its workers serve."""

LONG_QUEUE_NAME = redisqueue.LONG_QUEUE
"""str: The queue a job expected to run longer than `queue.long_job_seconds` goes on (PERF.151)."""

RESULT_TTL_SECONDS = 86400
"""int: How long a finished job's result or error is kept."""

UNESTIMATED_SECONDS = 3600
"""int: The time limit of a job with no estimate (PERF.144): a one-system edit, a
regenerate of one sector. A fixed ceiling, not an estimate."""

ESTIMATED_CEILING_SECONDS = 7 * 86400
"""int: The longest limit of a job with an estimate (a neighbourhood of a great many sectors can run for days)."""

JOB_ID_RE = re.compile(r"^[0-9a-f]{16}$")
"""re.Pattern: A job id. Checked before any id from a URL reaches Redis."""

_STATES = {
    "queued": "queued", "deferred": "queued", "scheduled": "queued",
    "started": "running", "finished": "succeeded",
    "failed": "failed", "stopped": "failed", "canceled": "failed",
}
"""dict: RQ's job states -> the API's."""

class NoQueue(Exception):
    """No Redis server answers, so nothing can be queued."""


def new_job_id():
    return secrets.token_hex(8)


SHORT_WAIT_SECONDS = 8.0
"""float: How long a route that queues quick work waits for the answer
before it gives up and answers `202` with the job's id."""


CONFIG_TAG = "__mysql_config__"
"""str: Marks a MySQL connection setting among a job's JSON arguments (`encode_argument`)."""


def function_reference(function):
    """`module:name` of a module-level function, the form a job names it in (JSON can't carry the function).

    Raises:
        ValueError: A function a worker could not import by that name.
    """
    module = getattr(function, "__module__", None)
    name = getattr(function, "__qualname__", None)
    if not module or not name or "<" in name or "." in name:
        raise ValueError(f"{function!r} is not a module-level function, so a worker cannot import it")
    return f"{module}:{name}"


def encode_argument(value):
    """A job argument as plain JSON data: numbers, strings, lists and dicts as they are, the MySQL connection
    settings tagged so `decode_argument` makes them again."""
    if type(value).__name__ == "MySQLConfig":
        return {CONFIG_TAG: dict(vars(value))}
    return value


def decode_argument(value):
    """`encode_argument` undone."""
    if isinstance(value, dict) and set(value) == {CONFIG_TAG}:
        from planetgen.db import store
        return store.MySQLConfig(**value[CONFIG_TAG])
    return value


def execute(reference, args):
    """
    What a worker runs for a queued job: the function named by `reference`
    (`function_reference`) called with `args` (`decode_argument`), with a refusal
    the work raises on purpose (an `ApiError`: 404, 409 ...) returned as
    data, so `status` can report the HTTP status it carries.

    `ApiError` is looked up among the loaded modules, not imported: a
    function that raises one has imported `planetgen.api.common` already,
    and importing it here would cost every job 2 s (PERF.39).
    """
    module, _, name = reference.partition(":")
    function = getattr(importlib.import_module(module), name)
    try:
        return {"result": function(*[decode_argument(arg) for arg in args])}
    except Exception as exc:  # noqa: BLE001 -- anything but a refusal is re-raised
        common = sys.modules.get("planetgen.api.common")
        if common is not None and isinstance(exc, common.ApiError):
            return {"refused": str(exc), "status": exc.status_code}
        raise


def _connection(url=None):
    """The shared client (`redisclient`) for `submit`, `status` and `wait`."""
    return redisqueue.connect(url or redisqueue.redis_url())


def time_limit(expected_seconds=None):
    """The seconds RQ lets an API job run (PERF.144): `redisqueue.job_timeout` of `expected_seconds`, at most
    `ESTIMATED_CEILING_SECONDS`, and `UNESTIMATED_SECONDS` for a job with no estimate."""
    return redisqueue.job_timeout(expected_seconds, ceiling=ESTIMATED_CEILING_SECONDS, unestimated=UNESTIMATED_SECONDS)


def recorded_task_seconds(kind, config=None, workers=1):
    """
    How long a task of `kind` (a generation stats kind: `sector`, ...) took a worker in earlier runs
    (`generation.stats`), the estimate of a job that does one; `None` with nothing recorded or no control
    database to read (the test suite turns the stats off, as the generation does).
    """
    from planetgen.db import store
    from planetgen.generation import stats
    if os.environ.get(stats.STATS_ENV_VAR, "1").strip().lower() in ("0", "off", "no", "false"):
        return None
    try:
        return stats.GenerationStats(store.control_mysql_config(config or store.DEFAULT_MYSQL_CONFIG)) \
            .seconds_per_task(kind, workers)
    except Exception:  # noqa: BLE001 -- no estimate then: the job gets the fixed limit
        return None


def is_long(expected_seconds):
    """Whether a job expected to take `expected_seconds` goes to the long lane (PERF.151): more than
    `queue.long_job_seconds`. A job with no estimate counts as short."""
    from planetgen.util.settings import get_settings
    if expected_seconds is None or not math.isfinite(expected_seconds):
        return False
    return expected_seconds > get_settings().queue.long_job_seconds


def lane(expected_seconds=None):
    """
    The queue and the worker cap of a job expected to take `expected_seconds`: `(LONG_QUEUE_NAME, queue.long_workers)`
    when `is_long`, else `(QUEUE_NAME, queue.interactive_workers)`.
    """
    from planetgen.util.settings import get_settings
    section = get_settings().queue
    if is_long(expected_seconds):
        return LONG_QUEUE_NAME, section.long_workers
    return QUEUE_NAME, section.interactive_workers


def submit(function, *args, expected_seconds=None):
    """
    Queues `function(*args)` (a module-level function; its arguments and
    return value travel as JSON, `encode_argument`) on its lane (`lane`: a job expected to run
    longer than `queue.long_job_seconds` goes to the long queue) and starts a worker for it unless
    enough are serving the queue already (`workers.ensure_workers`). RQ ends the job after
    `time_limit(expected_seconds)`, so a job whose worker was killed expires.

    Returns:
        str: The job's id, for `status`.

    Raises:
        NoQueue: No Redis server answers at `redis.url`, or an update is running (OPS.53).
    """
    try:
        updateguard.refuse_while_updating("This job")      # OPS.53: the queue takes nothing during an update
    except updateguard.UpdateRunning as exc:
        raise NoQueue(str(exc)) from exc
    url = redisqueue.redis_url()
    connection = _connection(url)
    try:
        connection.ping()
    except Exception as exc:  # noqa: BLE001 -- any connection failure
        raise NoQueue(f"no Redis server answers at {url} (config.json's redis.url); "
                      f"start Redis to run this: {exc}") from exc
    job_id = new_job_id()
    queue_name, cap = lane(expected_seconds)
    redisqueue.queue(queue_name, connection).enqueue(
        execute, args=(function_reference(function), [encode_argument(arg) for arg in args]), job_id=QUEUE_PREFIX + job_id,
        job_timeout=time_limit(expected_seconds),
        result_ttl=RESULT_TTL_SECONDS, failure_ttl=RESULT_TTL_SECONDS,
    )
    workers.ensure_workers(queue_name, connection, url, count=cap)
    return job_id


def _cap_of(queue_name):
    """The worker cap of the lane whose queue is `queue_name`."""
    from planetgen.util.settings import get_settings
    section = get_settings().queue
    return section.long_workers if queue_name == LONG_QUEUE_NAME else section.interactive_workers


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


def cancel(job_id):
    """
    Removes job `job_id` from the queue if no worker has taken it yet (PERF.105).

    Returns:
        bool: Whether it was removed; a job already running, finished or unknown is left alone.

    Raises:
        NoQueue: No Redis server answers.
    """
    if not isinstance(job_id, str) or not JOB_ID_RE.match(job_id):
        return False
    try:
        job = redisqueue.fetch_job(QUEUE_PREFIX + job_id, _connection(), redisqueue.serializer_for(QUEUE_NAME))
    except redisqueue.Unavailable as exc:
        raise NoQueue(str(exc)) from exc
    if job is None or _STATES.get(job.get_status(refresh=True)) != "queued":
        return False
    job.cancel()
    return True


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
            `None`), and `waiting` (PERF.153: for a running job that waits in line for the generation run in front of
            it, `{"run", "position", "text"}`, else `None`). `None` for an unknown or expired id, or a bad one.

    Raises:
        NoQueue: No Redis server answers.
    """
    if not isinstance(job_id, str) or not JOB_ID_RE.match(job_id):
        return None
    try:
        connection = _connection()
        job = redisqueue.fetch_job(QUEUE_PREFIX + job_id, connection, redisqueue.serializer_for(QUEUE_NAME))
    except redisqueue.Unavailable as exc:
        raise NoQueue(str(exc)) from exc
    if job is None:
        return None
    state = _STATES.get(job.get_status(refresh=True), "failed")
    if state == "queued":
        # A worker that was about to exit when the job arrived leaves it waiting: whoever asks starts another.
        try:
            workers.ensure_workers(job.origin, connection, count=_cap_of(job.origin))
        except Exception:  # noqa: BLE001 -- a status read does not fail for want of a worker
            pass
    body = {"id": job_id, "state": state, "result": None, "error": None, "error_status": None, "made": None,
            "waiting": job.meta.get("waiting") if state == "running" else None}
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


MAKES_SECTORS = ("generate_neighborhood", "generate_block", "regenerate_sector")
"""tuple[str]: The queued functions that generate sectors (ADM.31)."""


def _made_window(job):
    """The `made` window of a finished job that generates sectors: from
    when it started to a second past when it ended; else `None`."""
    function = job.args[0] if job.args else None
    if getattr(function, "__name__", None) not in MAKES_SECTORS or job.started_at is None or job.ended_at is None:
        return None
    return {"since": calendar.timegm(job.started_at.timetuple()), "until": calendar.timegm(job.ended_at.timetuple()) + 1}


def generate_neighborhood(sector_id, radius_ly, config, disregard_missing_data=False):
    """The queued body of `POST /api/sectors/<id>/generate-neighborhood`
    (a module-level function so a worker can import it by name): the sectors, the backfill and
    the settle as runs through the wrapper (PERF.119), the job waiting in line for the machine (PERF.153).
    `disregard_missing_data` is the request's GEN.599 switch: a missing
    downloaded data file then takes its labelled fallback."""
    from planetgen.datasets import fetch
    from planetgen.generation import run_galaxy
    fetch.set_disregard_missing_data(disregard_missing_data)
    return run_galaxy.generate_sector_neighborhood(sector_id, radius_ly=radius_ly, config=config)


def generate_block(block_key, config, disregard_missing_data=False):
    """The queued body of `POST /api/galaxy/blocks/<key>/generate` (GEN.310; a module-level function so a worker can
    import it by name): the sectors the outline holds in the block, the backfill and the settle as runs through the
    wrapper, the job waiting in line for the machine (PERF.153). `disregard_missing_data` is the request's GEN.599
    switch."""
    from planetgen.datasets import fetch
    from planetgen.generation import run_galaxy
    fetch.set_disregard_missing_data(disregard_missing_data)
    return run_galaxy.generate_block(block_key, config=config)


def empty_sector(sector_id, config):
    """The first half of `regenerate_sector` (GEN.845): deletes the galaxy-placed sector with everything in it and says what
    to make again.

    Returns:
        dict: `address` of its cell, `around` (the ids of the sectors whose paths it touched, not itself) and `deleted`
            (the row counts).
    """
    from planetgen.db import edits as editStore, sector_paths, store
    conn = store.get_connection(config or store.DEFAULT_MYSQL_CONFIG)
    try:
        with conn:
            address = editStore.sector_address(conn, sector_id)
            around = [other for other in sector_paths.sectors_to_settle(conn, [sector_id]) if other != sector_id]
            counts = editStore.delete_sector_with_contents(conn, sector_id)
    finally:
        conn.close()
    return {"address": address, "around": around, "deleted": counts}


def fill_emptied_sector(emptied, config):
    """The second half of `regenerate_sector` (GEN.845): generates the cell `empty_sector` emptied again, as runs through the
    wrapper (PERF.119: the sector and its backfill, then the settle of the new sector and the old neighbours that still
    exist).

    Returns:
        dict: `sector_id` and `sector_name` of the new sector (`None` only when the cell is beyond the grid's reach).
    """
    from planetgen.db import store
    from planetgen.generation import run_galaxy
    config = config or store.DEFAULT_MYSQL_CONFIG
    result = run_galaxy.regenerate_sector_cell(emptied["address"], config)
    # GEN.126: the new sector's masses and the old neighbours' paths, now the neighbour set is final.
    wanted = set(emptied["around"])
    if result["sector_id"] is not None:
        wanted.add(result["sector_id"])
    conn = store.get_connection(config)
    try:
        if wanted:
            marks = ", ".join("?" for _ in wanted)
            wanted = {row["id"] for row in conn.execute(f"SELECT id FROM sectors WHERE id IN ({marks})", sorted(wanted)).fetchall()}
    finally:
        conn.close()
    run_galaxy.settle_named_sectors(wanted, config, expand=False)
    return result


def regenerate_sector(sector_id, config):
    """The queued body of `POST /api/sectors/<id>/regenerate`: deletes the
    galaxy-placed sector with everything in it (`empty_sector`), then generates its cell
    again (`fill_emptied_sector`), the job waiting in line for the machine (PERF.153).

    Returns:
        dict: `deleted` (the row counts), `sector_id` and `sector_name` of
            the new sector (`None` only when the cell is beyond the grid's reach).
    """
    emptied = empty_sector(sector_id, config)
    result = fill_emptied_sector(emptied, config)
    return {"deleted": emptied["deleted"], "sector_id": result["sector_id"], "sector_name": result["sector_name"]}


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
    the command runs here instead.

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
