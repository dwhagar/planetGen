# planetgen/web/jobs.py

"""
Background jobs for the admin Generate page (`web/generate_page.py`):
resetting the galaxy, building its density skeleton (`planetgen plan`)
and generating sectors (`planetgen galaxy`), started from the browser
and left running after the request that started them has returned.

A job is a directory under the jobs directory (`jobs_dir()`):

    <jobs dir>/
        queued-<ticket>.json    the recipe of a job waiting in the line (0600; it holds the database environment)
        20260930-124433-1a2b/
            job.json            what to run: id, title, steps (argv lists), env
            runner.pid          the job's worker's pid, written when it is started
            cancel              written by Cancel; the runner stops when it sees it
            state.json          written by planetgen.web.job_runner as it goes
            progress.json       written by planetgen (planetgen.queue.progress_file)
            output.log          every step's stdout and stderr

`start_job` writes `job.json`, takes the lock (a hash in Redis, `<prefix>:job:active`) and queues
`planetgen.web.job_runner.run(<job dir>)` as an RQ job on Redis, on a
queue of the job's own (PERF.24 step 3). It then starts one burst worker
for that queue (`planetgen.cli.worker`), named after the job, in its own
session, detached from the web server, so a mod_wsgi request timeout or a
graceful Apache reload doesn't stop it. The worker exits once the job is
done. The page then only ever reads these files, so it works the same
whichever server process answers the next request. Without a Redis
server at `redis.url`, no job starts.

Only one job runs at a time: generation, planning and reset all write the
same database. A lock whose runner has died (the server was rebooted mid
job) is noticed and cleared by the next `start_job`, and that job shows
as "interrupted". A job that finds another running does not fail: `submit`
puts it in the one waiting list (`planetgen.queue.runline`, PERF.154), up to
`queue.max_waiting` entries, and it starts by itself when the job before it
succeeded or was cancelled (`promote_queued`).

The database a job writes is the one this site shows (`helpers.db_name`),
passed to the child through the same `PLANETGEN_MYSQL_*` environment
variables every command-line tool reads, never on the command line.
"""

import contextlib
import dataclasses
import json
import os
import re
import secrets
import shutil
import signal
import sys
import tempfile
import time

import redis

from planetgen.admin import updateguard
from planetgen.queue import redisqueue, runline, runrecord
from planetgen.web.lib.privatedir import ensure_private_dir
from planetgen.util import log, stagemark
from planetgen.util.settings import get_settings

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
"""str: The planetGen checkout."""

GENERATE_COMMAND = ["-m", "planetgen.cli.generate"]
"""list: Runs the generator, `planetgen.cli.generate`."""
RESET_COMMAND = ["-m", "planetgen.cli.reset"]
"""list: Runs `planetgen.cli.reset`."""
DATASETS_COMMAND = ["-m", "planetgen.cli.datasets"]
"""list: Runs `planetgen.cli.datasets` (GEN.548): `status` is the data check in front of every generating job (GEN.590)."""

DEFAULT_JOBS_DIR = "/var/lib/planetGen/jobs"
"""str: Used when neither `PLANETGEN_JOBS_DIR` nor `config.json`'s
`jobs.dir` names a directory. Falls back to a `planetgen-jobs` folder in
the system temp directory when Apache can't create this one."""

DEFAULT_KEEP = 20
"""int: How many finished jobs (and their logs) are kept; older ones are
deleted when a new job starts."""

QUEUED_PREFIX = "queued-"
"""str: A waiting job's recipe is the file `queued-<ticket>.json` in the jobs directory (ADM.123): the same spec as a job's
`job.json`. Mode 0600 like `job.json`: it holds the same database environment, which does not go into Redis. The
waiting list (`planetgen.queue.runline`) holds the rest."""

CANCEL_NAME = "cancel"
"""str: Must match `planetgen.web.job_runner.CANCEL_NAME`."""

RECORDS_NAME = "records.log"
"""str: The job's log as records (UX.122), in the job's folder. Must match `planetgen.web.job_runner.RECORDS_FILE`."""

QUEUE_PREFIX = "planetgen-web-"
"""str: A job's RQ queue and its worker are this plus the job's id."""

STARTING_GRACE_SECONDS = 15
"""int: How long a just-spawned job may go without a `state.json` before
it counts as interrupted (the runner writes one as its first act)."""

RUNNER_START_LIMIT_SECONDS = 300
"""int: How long a job whose runner is still alive may go without a
`state.json` (a slow start on a busy server) before it counts as
interrupted anyway: where only the pid can be checked (no `/proc`), a reused pid can't block jobs for longer than this."""

JOB_ID_RE = re.compile(r"^\d{8}-\d{6}-[0-9a-f]{4}$")
"""re.Pattern: A job id: `YYYYMMDD-HHMMSS-xxxx`. Checked before any id
from a URL becomes part of a path."""

FINISHED = frozenset({"succeeded", "failed", "cancelled", "interrupted"})


class JobBusy(Exception):
    """Another job is still running; `job` is its `get_job` dict."""

    def __init__(self, job):
        super().__init__(f"job {job['id'] if job else '?'} is still running")
        self.job = job


def _config():
    return get_settings().jobs


def configured_jobs_dir():
    """`PLANETGEN_JOBS_DIR`, then `jobs.dir`, then `DEFAULT_JOBS_DIR`.
    `examples/apache/create-cache-dir.sh` asks this where to create it."""
    return _config().dir or DEFAULT_JOBS_DIR


def jobs_dir():
    """
    The writable jobs directory, created if needed. When the default
    directory can't be created, falls back to a private `planetgen-jobs`
    in the system temp directory (see `privatedir.ensure_private_dir`:
    one another local user created or can write to is refused).

    Raises:
        OSError: When no candidate directory is writable.
    """
    configured = configured_jobs_dir()
    candidates = [configured]
    if configured == DEFAULT_JOBS_DIR:
        candidates.append(os.path.join(tempfile.gettempdir(), "planetgen-jobs"))
    for index, candidate in enumerate(candidates):
        try:
            if index == 0:
                os.makedirs(candidate, mode=0o750, exist_ok=True)
            else:
                ensure_private_dir(candidate)
        except OSError as exc:
            if index > 0:
                log.debug("jobs: not using %s: %s", candidate, exc)
            continue
        if os.access(candidate, os.W_OK | os.X_OK):
            return candidate
    raise OSError(f"no writable jobs directory (tried {', '.join(candidates)}); "
                  f"set jobs.dir in config.json to a directory Apache's user can write")


def keep_count():
    return _config().keep


def python_executable():
    """
    The Python that runs jobs: `PLANETGEN_PYTHON`, then `jobs.python`,
    then this process's own interpreter. Under mod_wsgi `sys.executable`
    can be Apache itself rather than Python, so that case falls back to
    `<sys.prefix>/bin/python3`, then `python3` on the PATH.
    """
    configured = _config().python
    if configured:
        return configured
    if sys.executable and os.path.basename(sys.executable).startswith("python"):
        return sys.executable
    candidate = os.path.join(sys.prefix, "bin", "python3")
    if os.path.exists(candidate):
        return candidate
    return shutil.which("python3") or "python3"


def mysql_env(mysql_config, database):
    """The `PLANETGEN_MYSQL_*` variables that point a child process at
    `database` on the server `mysql_config` describes."""
    return {
        "PLANETGEN_MYSQL_HOST": str(mysql_config.host),
        "PLANETGEN_MYSQL_PORT": str(mysql_config.port),
        "PLANETGEN_MYSQL_USER": str(mysql_config.user),
        "PLANETGEN_MYSQL_PASSWORD": str(mysql_config.password or ""),
        "PLANETGEN_MYSQL_DATABASE": str(database),
    }


# ---------------------------------------------------------------------
# Reading jobs
# ---------------------------------------------------------------------

def _read_json(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def _job_dir(root, job_id):
    if not isinstance(job_id, str) or not JOB_ID_RE.match(job_id):
        return None
    return os.path.join(root, job_id)


def _runner_alive(pid, job_id):
    """Whether `pid` is still this job's runner. Reads the process's
    command line where `/proc` exists."""
    if not pid:
        return False
    cmdline = f"/proc/{int(pid)}/cmdline"
    if os.path.isdir("/proc/self"):
        try:
            with open(cmdline, "rb") as f:
                return job_id.encode() in f.read()
        except OSError:
            return False
    try:
        os.kill(int(pid), 0)
    except (OSError, ValueError):
        return False
    return True


def _runner_pid(path, state):
    pid = (state or {}).get("pid")
    if pid:
        return pid
    try:
        with open(os.path.join(path, "runner.pid"), "r", encoding="utf-8") as f:
            return int(f.read().strip() or 0)
    except (OSError, ValueError):
        return None


def get_job(job_id, root=None):
    """
    One job as the page shows it, or `None` for an unknown id.

    Returns:
        dict: `job.json`'s fields (`id`, `title`, `kind`, `steps` as
            labels, `created_at`, `admin`, `database`) plus `status`
            (`starting`, `running`, `succeeded`, `failed`, `cancelled`,
            `interrupted`), `step` (1-based, 0 before the first),
            `step_label`, `step_started_at`, `step_estimates` (seconds per step in earlier
            runs, `None` where unknown; PERF.55), `started_at`, `finished_at`, `elapsed_s`,
            `error`, and `progress` (`progressFile`'s dict or `None`).
    """
    root = root or jobs_dir()
    path = _job_dir(root, job_id)
    if path is None:
        return None
    job = _read_json(os.path.join(path, "job.json"))
    if not job:
        return None
    state = _read_json(os.path.join(path, "state.json")) or {}
    status = state.get("status")
    now = time.time()
    if status is None:
        # A job is starting for the grace period whatever its pid says
        # (TEST.92: a runner just started may not have exec'd yet, so its
        # command line doesn't name the job), then while its runner is
        # still alive, until `RUNNER_START_LIMIT_SECONDS` (TEST.90: calling
        # a slow start interrupted let a finished-looking job come back as
        # running).
        pid = _runner_pid(path, state)
        age = now - job.get("created_at", 0)
        if age < STARTING_GRACE_SECONDS:
            status = "starting"
        elif pid is not None and age < RUNNER_START_LIMIT_SECONDS and _runner_alive(pid, job_id):
            status = "starting"
        else:
            status = "interrupted"
    elif status == "running" and not _runner_alive(state.get("pid"), job_id):
        status = "interrupted"

    labels = [step["label"] for step in job.get("steps", [])]
    step = state.get("step") or 0
    stages = job_stages(job.get("steps", []))
    started = state.get("started_at")
    finished = state.get("finished_at")
    progress = _read_json(os.path.join(path, "progress.json")) if step else None
    error = state.get("error")
    if status == "interrupted" and not error:
        error = "The job stopped without finishing (the server may have restarted)."
    in_step = sum(1 for stage in stages if stage["step"] == step)
    step_marker = (stagemark.marker(next(stage["n"] for stage in stages if stage["step"] == step), len(stages), in_step)
                   if in_step else (stagemark.marker(step, len(labels)) if step else ""))
    return {
        "id": job["id"],
        "kind": job.get("kind"),
        "title": job.get("title"),
        "admin": job.get("admin"),
        "database": job.get("database"),
        "draft": job.get("draft") or [],
        "created_at": job.get("created_at"),
        "steps": labels,
        "status": status,
        "finished": status in FINISHED,
        "step": step,
        "step_label": labels[step - 1] if 0 < step <= len(labels) else None,
        "step_marker": step_marker,
        "stages": stages,
        "step_started_at": state.get("step_started_at"),
        "step_estimates": state.get("step_estimates") or [],
        "stage_estimates": state.get("stage_estimates") or [],
        "started_at": started,
        "finished_at": finished,
        "elapsed_s": ((finished or now) - started) if started else None,
        "error": error,
        "progress": progress,
    }


def job_stages(steps):
    """
    Every stage of a job, in order (UX.89): `{"n", "step", "label", "skipped"}`
    for each stage of each step (a step with no stage list is one stage of its
    own label), `skipped` the reason a stage will not run or `None`.
    """
    stages = []
    for index, step in enumerate(steps, start=1):
        for stage in step.get("stages") or [{"label": step["label"], "skipped": None}]:
            stages.append({"n": len(stages) + 1, "step": index, "label": stage["label"],
                           "skipped": stage.get("skipped")})
    return stages


def _job_names(root):
    """Every job's id under `root`, newest first (`[]` when unreadable)."""
    try:
        return sorted((name for name in os.listdir(root) if JOB_ID_RE.match(name)), reverse=True)
    except OSError:
        return []


def count_jobs(root=None):
    """How many jobs are kept under `root`."""
    return len(_job_names(root or jobs_dir()))


def list_jobs(limit=10, root=None, offset=0):
    """The newest `limit` jobs from `offset`, newest first (`get_job` dicts)."""
    root = root or jobs_dir()
    names = _job_names(root)[offset:]
    jobs = []
    for name in names:
        job = get_job(name, root)
        if job is not None:
            jobs.append(job)
        if len(jobs) >= limit:
            break
    return jobs


def active_job(root=None):
    """The running job (`get_job` dict), or `None`.

    Raises:
        OSError: No writable jobs directory, or no Redis server answers (the lock is kept there).
    """
    root = root or jobs_dir()
    with _redis_errors():
        holder = _read_lock(_connection())
    if not holder:
        return None
    job = get_job(holder, root)
    if job is None or job["finished"]:
        return None
    return job


def remaining_steps(job_id, root=None):
    """
    What Retry on the admin queue page (ADM.10) runs again for a finished
    job that didn't succeed: its steps from the one that failed, was
    cancelled or was cut off, onward (a step that already succeeded
    isn't run twice).

    Returns:
        tuple | None: `(job, steps)` -- the `get_job` dict and the
            `{"label", "argv"}` steps -- or `None` for an unknown,
            running or succeeded job.
    """
    root = root or jobs_dir()
    job = get_job(job_id, root)
    if job is None or not job["finished"] or job["status"] == "succeeded":
        return None
    spec = _read_json(os.path.join(_job_dir(root, job_id), "job.json")) or {}
    steps = spec.get("steps") or []
    first = max(int(job.get("step") or 1), 1) - 1
    remaining = steps[first:]
    return (job, remaining) if remaining else None


def log_tail(job_id, max_bytes=64 * 1024, root=None):
    """
    The last `max_bytes` of a job's output, as text, with rich's carriage
    return redraws collapsed to their final state. `None` for an unknown
    job.
    """
    root = root or jobs_dir()
    path = _job_dir(root, job_id)
    if path is None or not os.path.isfile(os.path.join(path, "job.json")):
        return None
    try:
        with open(os.path.join(path, "output.log"), "rb") as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            f.seek(max(0, size - max_bytes))
            data = f.read()
    except OSError:
        return ""
    text = data.decode("utf-8", errors="replace")
    if size > max_bytes:
        text = text.split("\n", 1)[-1] if "\n" in text else text
    lines = [line.rstrip("\r").rsplit("\r", 1)[-1] for line in text.split("\n")]
    return "\n".join(lines)


def log_size(job_id, root=None):
    """The size in bytes of a job's `output.log` (0 before the first output); `None` for an unknown job."""
    root = root or jobs_dir()
    path = _job_dir(root, job_id)
    if path is None or not os.path.isfile(os.path.join(path, "job.json")):
        return None
    try:
        return os.path.getsize(os.path.join(path, "output.log"))
    except OSError:
        return 0


def read_log(job_id, offset=0, max_bytes=64 * 1024, root=None):
    """
    The raw output a job has written from byte `offset` on (at most
    `max_bytes`), for the live log stream (ADM.22): carriage returns and
    colour codes are left in, since the browser's terminal draws them.

    Returns:
        tuple: `(text, next_offset)`. A multi-byte character cut by the read
            is left for the next call, so `next_offset` always falls between
            characters. An offset past the end of the file (the log was
            replaced) starts again from 0. `("", offset)` when there is
            nothing new; `None` for an unknown job.
    """
    root = root or jobs_dir()
    path = _job_dir(root, job_id)
    if path is None or not os.path.isfile(os.path.join(path, "job.json")):
        return None
    try:
        with open(os.path.join(path, "output.log"), "rb") as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            if offset > size or offset < 0:
                offset = 0
            f.seek(offset)
            data = f.read(max_bytes)
    except OSError:
        return "", max(offset, 0)
    for cut in range(4):
        try:
            text = data[:len(data) - cut].decode("utf-8")
        except UnicodeDecodeError:
            continue
        return text, offset + len(data) - cut
    return data.decode("utf-8", errors="replace"), offset + len(data)


# ---------------------------------------------------------------------
# Starting and stopping jobs
# ---------------------------------------------------------------------

class _LockHeld(Exception):
    """The job lock is taken."""


@contextlib.contextmanager
def _redis_errors():
    """A Redis failure as the `OSError` that the pages and `start_job`'s callers already handle."""
    try:
        yield
    except redis.exceptions.RedisError as exc:
        raise OSError(f"no Redis server answers (config.json's redis.url); start Redis to run jobs: {exc}") from exc


def _connection():
    return redisqueue.connect()


def _lock_key():
    return runrecord._k("job", "active")


_TAKE = """
if redis.call('HSETNX', KEYS[1], 'id', ARGV[1]) == 1 then
  redis.call('HSET', KEYS[1], 'at', ARGV[2])
  return 1
end
return 0
"""

_RELEASE = """
if redis.call('HGET', KEYS[1], 'id') == ARGV[1] then return redis.call('DEL', KEYS[1]) end
return 0
"""


def _write_lock(connection, job_id):
    """Creates the job lock holding `job_id` and the time, in one step, so nobody ever reads it half made.

    Raises:
        _LockHeld: The lock exists.
    """
    if not connection.register_script(_TAKE)(keys=[_lock_key()], args=[job_id, time.time()]):
        raise _LockHeld(job_id)


def _read_lock(connection):
    holder = connection.hget(_lock_key(), "id")
    return None if holder is None else runrecord._s(holder).strip()


def _lock_age(connection):
    taken = connection.hget(_lock_key(), "at")
    try:
        return time.time() - float(runrecord._s(taken))
    except (TypeError, ValueError):
        return None


def _path_age(path):
    try:
        return time.time() - os.path.getmtime(path)
    except OSError:
        return None


def _take_lock(root, job_id):
    """
    Creates the job lock for `job_id`, clearing a stale one. A lock counts as stale when the job it names has finished,
    or is unknown (garbage, an id whose `job.json` is missing or unreadable) and the lock is older than
    `STARTING_GRACE_SECONDS` -- a younger one may belong to a job still being written.

    Raises:
        JobBusy: A live job (or one still starting) holds it.
    """
    with _redis_errors():
        connection = _connection()
        for _attempt in range(3):
            try:
                _write_lock(connection, job_id)
                return
            except _LockHeld:
                pass
            holder = _read_lock(connection)
            if holder is None:
                continue  # gone meanwhile
            job = get_job(holder, root) if JOB_ID_RE.match(holder) else None
            if job is not None and not job["finished"]:
                raise JobBusy(job)
            if job is None:
                age = _lock_age(connection)
                if age is None:
                    continue
                if age < STARTING_GRACE_SECONDS:
                    raise JobBusy(None)
            connection.register_script(_RELEASE)(keys=[_lock_key()], args=[holder])        # only if it still names the holder
    raise JobBusy(active_job(root))


def release_lock(job_id):
    """Frees the job lock if it names `job_id`."""
    try:
        with _redis_errors():
            _connection().register_script(_RELEASE)(keys=[_lock_key()], args=[job_id])
    except OSError as exc:
        log.warning("jobs: could not free the lock for job %s: %s", job_id, exc)


def free_lock():
    """Frees the job lock whoever holds it (Stop everything)."""
    with _redis_errors():
        _connection().delete(_lock_key())


def _prune(root, keep):
    """Deletes all but the newest `keep` jobs' directories. Never the
    one the lock names, one still running or starting, or one whose
    `job.json` can't be read yet less than `STARTING_GRACE_SECONDS` old."""
    names = sorted((name for name in os.listdir(root) if JOB_ID_RE.match(name)), reverse=True)
    with _redis_errors():
        holder = _read_lock(_connection())
    for name in names[keep:]:
        if name == holder:
            continue
        job = get_job(name, root)
        if job is not None and not job["finished"]:
            continue
        if job is None:
            age = _path_age(os.path.join(root, name))
            if age is None or age < STARTING_GRACE_SECONDS:
                continue
        shutil.rmtree(os.path.join(root, name), ignore_errors=True)


def new_job_id():
    return time.strftime("%Y%m%d-%H%M%S") + "-" + secrets.token_hex(2)


def start_job(kind, title, steps, env=None, admin=None, database=None, root=None, spawn=True, draft=None):
    """
    Starts a job in the background and returns its id straight away.

    Args:
        kind (str): `new_galaxy`, `new_plan`, `plan`, `galaxy`, `bring_to_epoch`, `finish_run`, `check_db` or `reset`.
        title (str): What the page calls it ("Generate sectors").
        steps (list[dict]): `{"label", "argv"}` per step, run in order; the generation commands among them run as the one
            step `planetgen run-job` (`job_steps.collapse`, PERF.203).
        env (dict, optional): Extra environment for every step (the
            `mysql_env` variables).
        admin (str, optional): Who started it, for the page.
        database (str, optional): Which database it writes, for the page.
        spawn (bool): Tests pass `False` to only write the files.
        draft (list, optional): The form fields (`[name, value]` pairs) the job was made from, kept so History can offer
            Run again with the same settings (ADM.115).

    Raises:
        JobBusy: Another job is still running.
        OSError: No writable jobs directory, no Redis server, or the
            worker can't start.
    """
    updateguard.refuse_while_updating("A Generate job")     # OPS.53
    from planetgen.web import job_steps
    steps = job_steps.collapse(steps, title)
    root = root or jobs_dir()
    # The directory first, under a fresh id if two jobs drew the same one
    # in the same second (TEST.41); then the lock.
    for _attempt in range(20):
        job_id = new_job_id()
        path = os.path.join(root, job_id)
        try:
            os.makedirs(path, mode=0o750)
            break
        except FileExistsError:
            continue
    else:
        raise OSError(f"could not create a job directory in {root}")
    try:
        _take_lock(root, job_id)
    except BaseException:
        shutil.rmtree(path, ignore_errors=True)
        raise
    try:
        job = {
            "id": job_id, "kind": kind, "title": title, "admin": admin, "database": database,
            "created_at": time.time(), "cwd": REPO_DIR, "steps": steps, "env": env or {}, "draft": draft or [],
        }
        fd = os.open(os.path.join(path, "job.json"), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(job, f)
        if spawn:
            _spawn(path)
    except BaseException:
        release_lock(job_id)
        raise
    log.debug("Started job %s (%s) for %s on %s: %s", job_id, kind, admin, database,
              [step["label"] for step in steps])
    try:
        _prune(root, keep_count())
    except OSError:
        pass
    return job_id


JOB_FLOOR_SECONDS = 3600
"""int: The shortest time limit of a job's RQ entry (PERF.144): a job of steps that usually take seconds (a reset, a
download) is not ended by one slow start."""

JOB_CEILING_SECONDS = 30 * 86400
"""int: The longest time limit of a job's RQ entry, and the limit of a job whose steps are not all recorded (a fixed
ceiling, not an estimate): a whole galaxy run takes hours to days."""


def job_time_limit(path):
    """
    The seconds RQ lets the job in `path` run (PERF.144): `redisqueue.job_timeout` of the time its steps are expected
    to take (`job_runner.expected_seconds`: the recorded stage times, else what each step took in earlier runs),
    `JOB_CEILING_SECONDS` when a step has no record.
    """
    from planetgen.web import job_runner

    with open(os.path.join(path, "job.json"), "r", encoding="utf-8") as f:
        job = json.load(f)
    return redisqueue.job_timeout(job_runner.expected_seconds(job), ceiling=JOB_CEILING_SECONDS,
                                  floor=JOB_FLOOR_SECONDS)


def _spawn(path):
    """Queues the job in `path` on Redis and starts its burst worker.

    Raises:
        OSError: No Redis server, or the process can't start.
    """
    job_id = os.path.basename(path)
    url = redisqueue.redis_url()
    connection = redisqueue.connect(url)
    try:
        connection.ping()
    except Exception as exc:  # noqa: BLE001 -- any connection failure
        raise OSError(f"no Redis server answers at {url} (config.json's redis.url); "
                      f"start Redis to run jobs: {exc}") from exc
    name = QUEUE_PREFIX + job_id
    redisqueue.queue(name, connection).enqueue(
        "planetgen.web.job_runner.run", path, job_id=name, job_timeout=job_time_limit(path),
        result_ttl=0, failure_ttl=86400,
    )
    _start_detached(path, redisqueue.worker_argv([name], url, name=name, python=python_executable()))


def _start_detached(path, argv):
    """Starts `argv` (whose command line names the job) detached from the
    web server and records its pid as the job's `runner.pid`."""
    proc = redisqueue.start_detached(argv, REPO_DIR)
    with open(os.path.join(path, "runner.pid"), "w", encoding="utf-8") as f:
        f.write(str(proc.pid))


def cancel_runs(job_id, who=None):
    """
    Cancels the live runs the job started, now (PERF.99, `planetgen.queue.runcontrol.cancel_run`): their waiting units
    are taken off the queue, their started units are sent the stop command and the machine is free, without waiting
    for the step that follows them to notice. Best effort: with no Redis the step's own stop is all there is.

    Returns:
        int: How many runs were cancelled.
    """
    try:
        from planetgen.queue import runcontrol, runrecord

        connection = redisqueue.connect()
        return sum(runcontrol.cancel_run(connection, run_id, who=who) is not None
                   for run_id in runrecord.job_runs(connection, job_id))
    except Exception as exc:  # noqa: BLE001 -- the runner's own stop still follows
        log.debug("Cancel of job %s: its runs were not cancelled (%s)", job_id, exc)
        return 0


def _end_orphan(root, path, job_id, who):
    """
    Cancel of a job whose runner is gone (killed by the out-of-memory killer, a restart or a crash; PERF.186): its step
    and the processes of that step's group were started apart from the runner and run on, so they are ended here
    (SIGTERM, then SIGKILL), its live runs are cancelled, and the job is marked cancelled and its lock freed. Before
    this, Cancel on such a job did nothing, and the page showed it as interrupted while its processes kept the host busy.

    Returns:
        bool: Whether there was something to end or the job was not marked yet.
    """
    state_path = os.path.join(path, "state.json")
    state = _read_json(state_path) or {}
    did = False
    step_pid = state.get("step_pid")
    if isinstance(step_pid, int) and step_pid > 1:
        for number, wait in ((signal.SIGTERM, 5.0), (signal.SIGKILL, 2.0)):
            try:
                os.killpg(step_pid, number)
            except OSError:
                break
            did = True
            deadline = time.time() + wait
            while time.time() < deadline:
                try:
                    os.killpg(step_pid, 0)
                except OSError:
                    break
                time.sleep(0.1)
            else:
                continue
            break
    did = bool(cancel_runs(job_id, who)) or did
    if state.get("status") in (None, "running"):
        state.update(status="cancelled", finished_at=time.time(), error=f"Cancelled by {who or 'an admin'} after its runner "
                                                                         "had stopped.")
        state.pop("step_pid", None)
        try:
            with open(state_path, "w", encoding="utf-8") as f:
                json.dump(state, f)
            did = True
        except OSError:
            pass
    release_lock(job_id)
    return did


def cancel_job(job_id, root=None, who=None):
    """
    Asks a running job to stop: writes a `cancel` file into its directory, which the runner checks for while a step
    runs, and cancels the job's live runs through the stop command (`cancel_runs`). The runner then
    stops the step and everything the step started, and marks the job
    cancelled. (A file rather than a signal, so the runner stops the step
    and releases the lock itself.)

    Returns:
        bool: Whether a running job was asked to stop.
    """
    root = root or jobs_dir()
    job = get_job(job_id, root)
    if job is None or (job["finished"] and job["status"] != "interrupted"):
        return False
    path = _job_dir(root, job_id)
    pid = _runner_pid(path, _read_json(os.path.join(path, "state.json")))
    if not _runner_alive(pid, job_id):
        return _end_orphan(root, path, job_id, who)
    try:
        with open(os.path.join(path, CANCEL_NAME), "w", encoding="utf-8") as f:
            f.write(str(time.time()))
    except OSError:
        return False
    cancel_runs(job_id, who)        # after the file: the step that finds its run stopped is then a cancelled one, not a failed one
    log.debug("Cancel requested for job %s (runner pid %s)", job_id, pid)
    return True


# ---------------------------------------------------------------------
# Jobs waiting in the line (ADM.123, PERF.154)
# ---------------------------------------------------------------------

@dataclasses.dataclass
class Submitted:
    """What `submit` did: `job_id` the job it started now, or `ticket` and `position` (1 is next) of the place it took in
    the waiting list, `running` the title of the job it waits for."""

    job_id: str = None
    ticket: str = None
    position: int = None
    running: str = None


def _spec_path(root, ticket):
    return os.path.join(root, f"{QUEUED_PREFIX}{ticket}.json")


def _write_spec(path, spec):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        json.dump(spec, f)


def _remove(path):
    with contextlib.suppress(OSError):
        os.remove(path)


def submit(kind, title, steps, env=None, admin=None, database=None, root=None, draft=None):
    """
    Starts the job, or, when another job is running, puts it in the waiting list (PERF.154): it starts by itself when the
    job before it succeeded or was cancelled and is held for an admin when that job failed.

    Returns:
        Submitted: The job started, or the place taken.

    Raises:
        runline.LineFull: The list holds `queue.max_waiting` entries; the message names them.
        updateguard.UpdateRunning, OSError: As `start_job`.
    """
    root = root or jobs_dir()
    try:
        return Submitted(job_id=start_job(kind, title, steps, env=env, admin=admin, database=database, root=root,
                                          draft=draft))
    except JobBusy as busy:
        running = busy.job["title"] if busy.job else "Another job"
    ticket = queue_next(kind, title, steps, env=env, admin=admin, database=database, root=root, draft=draft)
    started = promote_queued(root)      # the running job may have ended since the lock was refused
    waiting = [entry["ticket"] for entry in queued_jobs(root)]
    if started and ticket not in waiting:
        return Submitted(job_id=started)
    return Submitted(ticket=ticket, position=waiting.index(ticket) + 1 if ticket in waiting else None, running=running)


def queue_next(kind, title, steps, env=None, admin=None, database=None, root=None, draft=None):
    """
    Puts a job at the end of the waiting list: it starts by itself when the job before it succeeded or was cancelled, and
    is held for the admin to start or discard when that job failed.

    Returns:
        str: The entry's ticket.

    Raises:
        runline.LineFull: The list holds `queue.max_waiting` entries.
        OSError: No writable jobs directory, or no Redis server.
    """
    root = root or jobs_dir()
    running = active_job(root)
    ticket = secrets.token_hex(8)
    spec = {"kind": kind, "title": title, "admin": admin, "database": database, "steps": steps, "env": env or {},
            "draft": draft or []}
    path = _spec_path(root, ticket)
    _write_spec(path, spec)
    try:
        with _redis_errors():
            connection = _connection()
            ahead = runline.of_type(connection, runline.JOB, root)
            runline.join(connection, runline.JOB, title, who=admin, database=database, kind=kind, ticket=ticket,
                         fields={"root": root, "spec": os.path.basename(path),
                                 "after_job": running["id"] if running and not ahead else None})
    except BaseException:
        _remove(path)
        raise
    log.debug("Queued %s (%s) for %s behind %s", title, kind, admin, running["id"] if running else None)
    return ticket


def queued_jobs(root=None):
    """The jobs waiting in the list, next first, as `{"ticket", "kind", "title", "admin", "database", "queued_at",
    "after_job", "held", "position"}` (no steps or environment)."""
    root = root or jobs_dir()
    with _redis_errors():
        found = runline.of_type(_connection(), runline.JOB, root)
    return [{"ticket": entry["ticket"], "kind": entry["kind"], "title": entry["title"], "admin": entry["who"],
             "database": entry["database"], "queued_at": entry["queued_at"], "after_job": entry["after_job"],
             "held": entry["held"], "position": entry["position"]} for entry in found]


def line_view(root=None):
    """The whole waiting list as the strip shows it, in the order it arrived: `{"ticket", "type", "title", "who",
    "held", "position"}` for the jobs waiting for the running job and the runs waiting for the machine."""
    with _redis_errors():
        found = runline.entries(_connection())
    root = root or jobs_dir()
    return [{"ticket": entry["ticket"], "type": entry["type"], "title": entry["title"], "who": entry["who"],
             "held": entry["held"], "position": entry["position"]}
            for entry in found if entry["type"] == runline.RUN or entry["root"] == root]


def discard_queued(root=None, ticket=None):
    """Takes a waiting job out of the list (the next one when no `ticket` is given); True when there was one."""
    root = root or jobs_dir()
    with _redis_errors():
        connection = _connection()
        found = runline.of_type(connection, runline.JOB, root)
        entry = next((e for e in found if ticket is None or e["ticket"] == ticket), None)
        if entry is None:
            return False
        gone = runline.leave(connection, entry["ticket"])
    _remove(_spec_path(root, entry["ticket"]))
    return gone


def discard_all_queued(root=None):
    """Takes every waiting job out of the list; how many there were."""
    root = root or jobs_dir()
    count = 0
    while discard_queued(root):
        count += 1
    return count


def promote_queued(root=None, force=False):
    """
    Starts the job at the head of the waiting list if the job it waited for is over. Called when a job ends
    (`job_runner.run`), when an admin page asks for the running job, and by Start now (`force`, for a held job).

    Returns:
        str | None: The new job's id, or `None` (nothing waiting, still running, held, or the start failed).
    """
    root = root or jobs_dir()
    with _redis_errors():
        connection = _connection()
        waiting = runline.of_type(connection, runline.JOB, root)
        if not waiting or active_job(root) is not None:
            return None
        entry = waiting[0]
        after = get_job(entry["after_job"], root) if entry["after_job"] else None
        if not force and after is not None:
            if not after["finished"]:
                return None
            if after["status"] not in ("succeeded", "cancelled"):
                if not entry["held"]:
                    runline.set_field(connection, entry["ticket"], "held",
                                      f"The job before it did not succeed ({after['status']}).")
                return None
        if not runline.claim(connection, entry["ticket"]):        # one process wins this
            return None
    spec_path = _spec_path(root, entry["ticket"])
    spec = _read_json(spec_path)
    try:
        if not spec:
            raise OSError("its recipe is gone")
        job_id = start_job(spec["kind"], spec["title"], spec["steps"], env=spec.get("env"), admin=spec.get("admin"),
                           database=spec.get("database"), root=root, draft=spec.get("draft"))
    except JobBusy:
        with _redis_errors():
            runline.restore(connection, entry["ticket"])
        return None
    except Exception as exc:  # noqa: BLE001 -- keep the job for the admin to start or discard
        log.error("Could not start the waiting job %s: %s", entry["title"], exc)
        with _redis_errors():
            runline.set_field(connection, entry["ticket"], "held", f"It could not be started: {exc}")
            runline.restore(connection, entry["ticket"])
        return None
    with _redis_errors():
        runline.leave(connection, entry["ticket"])
        later = runline.of_type(connection, runline.JOB, root)
        if later:
            runline.set_field(connection, later[0]["ticket"], "after_job", job_id)      # it waits for the job just started
    _remove(spec_path)
    return job_id
