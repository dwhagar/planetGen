# planetgen/queue/api_jobs.py

"""
Long API work as light jobs (PERF.24 step 4, OPS.230): a route that would hold a web worker for minutes queues the work with `submit`, answers
`202 Accepted` with the job's id at once, and the client reads the outcome from `GET /api/jobs/<id>` (`status`).

The jobs are light jobs (`queue/lightjobs.py`): a record in Redis that the supervisor hands out, as a process forked from a warm template, when the
host has room for the job's claim (a short job on the interactive lane, a job expected to run longer than `queue.long_job_seconds` on the long lane,
PERF.151). The supervisor is started by the submit and by every status poll that finds the job waiting. A restart or reload of the web server does not
stop a job, and a job that outlives its time limit (`lightjobs.time_limit`) is ended. The result is kept for `RESULT_TTL_SECONDS`.
"""

import importlib
import os
import sys
import time

import redis

from planetgen.admin import updateguard
from planetgen.queue import jobbeat, lightjobs, redisqueue
from planetgen.util import log, secrets_mask
from planetgen.util.format import format_duration_seconds

RESULT_TTL_SECONDS = lightjobs.RESULT_TTL_SECONDS
"""int: How long a finished job's result or error is kept."""

JOB_ID_RE = lightjobs.JOB_ID_RE
"""re.Pattern: A job id. Checked before any id from a URL reaches Redis."""

NoQueue = lightjobs.NoQueue
"""Exception: No Redis server answers, or an update is running, so nothing can be queued."""


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


OWN_BEAT = frozenset({"planetgen.web.lib.tileprefetch:build_nodes"})
"""frozenset: The API functions that write a beat of their own (`jobbeat.beating` of a kind and label that say more than `api`): `execute` does not beat for them."""


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
    from planetgen.queue import allowlist
    module, _, name = allowlist.check_reference(reference).partition(":")      # SEC.45: only the API's own functions
    function = getattr(importlib.import_module(module), name)
    try:
        if reference in OWN_BEAT:                                               # the function beats under its own kind and label
            return {"result": function(*[decode_argument(arg) for arg in args])}
        with jobbeat.beating("api", lightjobs.current_id(), label=name):        # OPS.142: the Queue page sees this job is alive
            return {"result": function(*[decode_argument(arg) for arg in args])}
    except Exception as exc:  # noqa: BLE001 -- anything but a refusal is re-raised
        common = sys.modules.get("planetgen.api.common")
        if common is not None and isinstance(exc, common.ApiError):
            return {"refused": str(exc), "status": exc.status_code}
        raise


def _connection(url=None):
    """The shared client (`redisclient`) for `submit`, `status` and `wait`."""
    return redisqueue.connect(url or redisqueue.redis_url())


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
    except Exception as exc:  # noqa: BLE001 -- no estimate then: the job gets the fixed limit
        log.swallowed("queue/api_jobs", "recorded_task_seconds", exc)
        return None


def submit(function, *args, expected_seconds=None, wanted_seconds=None, drop_keys=(), lane=None, limit_seconds=None):
    """
    Queues `function(*args)` (a module-level function on `allowlist.API_FUNCTIONS`; its arguments and return value travel as JSON, `encode_argument`) as a light job
    (`lightjobs.add_job`) and asks for the supervisor that hands it out. A job expected to run longer than `queue.long_job_seconds` is on the long lane; it is ended after
    `lightjobs.time_limit(expected_seconds)`, so a job that hangs does not hold its claim.

    Args:
        expected_seconds (float, optional): How long it is expected to take.
        wanted_seconds (float, optional): How long it stays wanted once nobody asks after it: a job that has not started by then is dropped (W12 of the audit); without it the job is never dropped.
        drop_keys (iterable): Redis keys deleted when the job is dropped.
        lane (str, optional): The lane when it is not the one the expected time gives (`lightjobs.LANE_MAINTENANCE`, the node prefetch).
        limit_seconds (int, optional): The time limit when it is not the one the expected time gives.

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
        job_id = lightjobs.add_job(connection, function_reference(function), [encode_argument(arg) for arg in args], expected_seconds=expected_seconds,
                                   wanted_seconds=wanted_seconds, drop_keys=drop_keys, lane=lane, limit_seconds=limit_seconds)
    except (redis.exceptions.RedisError, OSError) as exc:
        raise NoQueue(f"no Redis server answers at {secrets_mask.mask_text(url)} (config.json's redis.url); start Redis to run this: {exc}") from exc
    lightjobs.ask_supervisor(connection, url)
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


def cancel(job_id):
    """
    Removes job `job_id` from the pending list if the supervisor has not taken it yet (PERF.105).

    Returns:
        bool: Whether it was removed; a job already running, finished or unknown is left alone.

    Raises:
        NoQueue: No Redis server answers.
    """
    try:
        return lightjobs.cancel(_connection(), job_id)
    except (redis.exceptions.RedisError, OSError) as exc:
        raise NoQueue(str(exc)) from exc


def status(job_id):
    """
    One job's state.

    Returns:
        dict | None: `id`, `state` (`queued`, `running`, `succeeded` or `failed`), `result` (what the function returned, once it succeeded), `error` (the refusal, or the last line of
            the error text, once it failed), `error_status` (the HTTP status of a refusal, else `None`) and `waiting` (PERF.153: for a running job that waits in line for the generation run in
            front of it, `{"run", "position", "text"}`, else `None`). `None` for an unknown or expired id, or a bad one. Asking after a queued job that has a deadline keeps it wanted.

    Raises:
        NoQueue: No Redis server answers.
    """
    try:
        return lightjobs.status(_connection(), job_id)
    except (redis.exceptions.RedisError, OSError) as exc:
        raise NoQueue(str(exc)) from exc


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


def build_gravity_grid(sector_id, config, folder):
    """The queued body of `GET /api/sectors/<id>/gravity` when the grid is not cached (API.24; a module-level function so a worker can import it by name):
    works the sector's local field out (`db.gravity_grid.grid`, which imports numba, so only here) and stores it in `folder`.

    Returns:
        dict: `sector_id`, `stamp` (the key of the file) and `built` (false when another job had stored it already).
    """
    from planetgen.db import gravity_grid, store
    from planetgen.galaxy import potential
    conn = store.get_connection(config or store.DEFAULT_MYSQL_CONFIG)
    try:
        info = store.get_galaxy_shape(conn, light=True)
        pot = potential.current() if info is None else potential.for_galaxy(info.shape, info.galaxy_seed)
        result = gravity_grid.grid(conn, sector_id, folder, pot=pot)
    finally:
        conn.close()
    return {"sector_id": sector_id, "stamp": result.stamp, "built": not result.hit}


def empty_sector(sector_id, config):
    """The first half of `regenerate_sector` (GEN.845): deletes the galaxy-placed sector with what was born in it
    (`db.edits.empty_sector`: the row stays when objects born elsewhere are filed in it) and says what to make again.

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
            counts = editStore.empty_sector(conn, sector_id)
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
        log.get("queue").error(f"The command {_command_words(argv)} timed out after {format_duration_seconds(timeout)}.")
    except OSError as exc:
        out["error"] = str(exc)
        log.get("queue").error(f"The command {_command_words(argv)} could not be started: {exc}")
    else:
        out["returncode"] = done.returncode
        out["stdout"] = (done.stdout or b"").decode("utf-8", errors="replace")
        out["stderr"] = (done.stderr or b"").decode("utf-8", errors="replace")
        log.process_ended("queue", f"The command {_command_words(argv)}", done.returncode, output=out["stderr"] or out["stdout"])
    return out


def _command_words(argv):
    """The first words of a command line, for a log line (a long argument list is not one)."""
    words = " ".join(str(part) for part in argv[:4])
    return words if len(argv) <= 4 else words + " ..."


def command_and_wait(argv, env, cwd, timeout, merge_stderr=False):
    """
    Runs a command on the queue and waits for it (the Generate page's estimate, PERF.24). Without a Redis server
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


def call_and_wait(function, *args, timeout):
    """
    Runs `function(*args)` as a light job in the supervisor's warm process (the one-off system page, OPS.244: no process of its own to start, no imports to pay for) and waits up to `timeout`
    seconds for what it returns.

    Raises:
        NoQueue: No Redis server answers, or an update is running: nothing was queued.
        RuntimeError: The job failed, was stopped for taking longer than `timeout`, or was not answered in time; the message says which.
    """
    job_id = submit(function, *args, expected_seconds=timeout, limit_seconds=int(timeout))
    job = wait(job_id, timeout + 15)
    if job is not None and job["state"] == "succeeded":
        return job["result"]
    if job is not None and job["state"] == "failed":
        raise RuntimeError(job["error"] or "the job failed")
    raise RuntimeError("the work queue did not answer in time")
