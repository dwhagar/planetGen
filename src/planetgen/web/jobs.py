# planetgen/web/jobs.py

"""
Background jobs for the admin Generate page (`web/generate_page.py`):
resetting the galaxy, building its density skeleton (`planetgen plan`)
and generating sectors (`planetgen galaxy`), started from the browser
and left running after the request that started them has returned.

A job is a run of the job ledger (`planetgen.queue.pagerun`, OPS.244): `start_job` writes a row for the job and a row for each of its steps and asks the supervisor
(`lightjobs.ask_supervisor`), which forks one job process for each step in turn (`planetgen.web.job_runner.run_step`), so a mod_wsgi request timeout or a graceful Apache
reload does not stop it. What the page shows (the status, the step, the times, why it failed) is read from the ledger's rows (`get_job`); Cancel is a request to the supervisor
(`controls.stop_engine_run`). The job's folder under the jobs directory (`jobs_dir()`) holds only what a process writes as it goes:

    <jobs dir>/
        queued-<ticket>.json    the recipe of a job waiting in the line (0600)
        20260930-124433-1a2b/
            output.log          every step's stdout and stderr
            records.log         the same as records (UX.122)
            progress.json       written by planetgen (planetgen.queue.progress_file)
            memory.json         the memory peak of each step
            step.pid            the process group of the step that runs, for Cancel when its runner is gone

The folder's name is the job's id; the jobs directory's listing is the list of jobs (`list_jobs`, newest first). Without a Redis server at `redis.url` or a control database, no job starts.

Only one job runs at a time: generation, planning and reset all write the
same database. A lock whose job has ended (the server was rebooted mid
job) is noticed and cleared by the next `start_job`, and that job shows
as "interrupted". A job that finds another running does not fail: `submit`
puts it in the one waiting list (`planetgen.queue.runline`, PERF.154), up to
`queue.max_waiting` entries, and it starts by itself when the job before it
succeeded or was canceled (`promote_queued`).

The database a job writes is the one this site shows (`helpers.db_name`),
named on the job; its connection variables (`PLANETGEN_MYSQL_*`, `mysql_env`) are built from the
site's settings when a step starts and given to the step's process, never stored or put on a command line.
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
from planetgen.queue import controls, lightjobs, pagerun, redisqueue, runline, runrecord, supervisor
from planetgen.web.lib.privatedir import ensure_private_dir
from planetgen.util import log, process, secrets_mask, stagemark
from planetgen.util.settings import get_settings
from planetgen.web import job_runner

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
`start_job` takes (`kind`, `title`, `steps`, `admin`, `database`, `draft`), mode 0600. The waiting list (`planetgen.queue.runline`) holds the rest."""

RECORDS_NAME = "records.log"
"""str: The job's log as records (UX.122), in the job's folder. Must match `planetgen.web.job_runner.RECORDS_FILE`."""

STARTING_GRACE_SECONDS = 15
"""int: How long a job nothing has begun may be younger than before a lock naming a job the ledger does not know counts as stale."""

JOB_ID_RE = re.compile(r"^\d{8}-\d{6}-[0-9a-f]{4}$")
"""re.Pattern: A job id: `YYYYMMDD-HHMMSS-xxxx`. Checked before any id
from a URL becomes part of a path."""

FINISHED = frozenset({"succeeded", "failed", "canceled", "interrupted"})


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
                log.debug("Not using %s as the jobs folder: %s.", candidate, exc)
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


def get_job(job_id, root=None):
    """
    One job as the page shows it, or `None` for an unknown id (one with no folder under `root` or no row in the ledger, `pagerun.read`).

    Returns:
        dict: The job's fields (`id`, `title`, `kind`, `steps` as
            labels, `created_at`, `admin`, `database`) plus `status`
            (`starting`, `running`, `succeeded`, `failed`, `canceled`,
            `interrupted`), `step` (1-based, 0 before the first),
            `step_label`, `step_started_at`, `step_estimates` (seconds per step in earlier
            runs, `None` where unknown; PERF.55), `started_at`, `finished_at`, `elapsed_s`,
            `error`, and `progress` (`progressFile`'s dict or `None`).
    """
    root = root or jobs_dir()
    path = _job_dir(root, job_id)
    if path is None or not os.path.isdir(path):
        return None
    with _ledger() as control:
        job = pagerun.read(control, job_id)
    if job is None:
        return None
    status = job["status"]
    now = time.time()
    labels = [step["label"] for step in job["steps"]]
    step = job["step"]
    stages = job_stages(job["steps"])
    started = job["started_at"]
    finished = job["finished_at"] if status in FINISHED else None
    progress = _read_json(os.path.join(path, "progress.json")) if step else None
    error = job["error"]
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
        "step_started_at": job["step_started_at"],
        "step_estimates": job["step_estimates"],
        "stage_estimates": job["stage_estimates"],
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
        connection = _connection()
        holder = _read_lock(connection)
    if not holder:
        return None
    job = get_job(holder, root)
    if job is None or job["finished"]:
        return None
    lightjobs.ask_supervisor(connection)        # a job that is not over has a supervisor to hand its steps out: one that ended (or never started) is started by the poll
    return job


def remaining_steps(job_id, root=None):
    """
    What Retry on the admin queue page (ADM.10) runs again for a finished
    job that didn't succeed: its steps from the one that failed, was
    canceled or was cut off, onward (a step that already succeeded
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
    with _ledger() as control:
        found = pagerun.read(control, job_id)
    steps = (found or {}).get("steps") or []
    first = max(int(job.get("step") or 1), 1) - 1
    remaining = steps[first:]
    return (job, remaining) if remaining else None


def _known_job_dir(root, job_id):
    """The folder of a job that exists under `root` (the jobs folder by default), else `None`."""
    path = _job_dir(root or jobs_dir(), job_id)
    if path is None or not os.path.isdir(path):
        return None
    return path


def log_tail(job_id, max_bytes=64 * 1024, root=None):
    """
    The last `max_bytes` of a job's output, as text, with rich's carriage
    return redraws collapsed to their final state. `None` for an unknown
    job.
    """
    path = _known_job_dir(root, job_id)
    if path is None:
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
    path = _known_job_dir(root, job_id)
    if path is None:
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
    path = _known_job_dir(root, job_id)
    if path is None:
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


@contextlib.contextmanager
def _ledger():
    """A connection to the control database (`pagerun.connect`), closed after use.

    Raises:
        OSError: There is no control database.
    """
    control = pagerun.connect()
    try:
        yield control
    finally:
        control.close()


def hand_over(connection, job_id):
    """Tells the supervisor a job is waiting (`lightjobs.ask_supervisor`: one is started when none is alive); it forks the process of each step. A test replaces this to run the steps itself."""
    lightjobs.ask_supervisor(connection)


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
    or is unknown (garbage, an id the ledger does not know) and the lock is older than
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
        log.warning("The lock for job %s could not be freed: %s.", job_id, exc)


def free_lock():
    """Frees the job lock whoever holds it (Stop everything)."""
    with _redis_errors():
        _connection().delete(_lock_key())


def _prune(root, keep):
    """Deletes all but the newest `keep` jobs' directories. Never the
    one the lock names, one still running or starting, or one the ledger does not know yet less than `STARTING_GRACE_SECONDS` old. The ledger's rows of a job
    that is pruned go with its folder."""
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
        with contextlib.suppress(OSError), _ledger() as control:
            pagerun.delete(control, name)


def new_job_id():
    return time.strftime("%Y%m%d-%H%M%S") + "-" + secrets.token_hex(2)


def start_job(kind, title, steps, admin=None, database=None, root=None, spawn=True, draft=None):
    """
    Starts a job in the background and returns its id straight away.

    Args:
        kind (str): `new_galaxy`, `new_plan`, `plan`, `galaxy`, `bring_to_epoch`, `finish_run`, `check_db` or `reset`.
        title (str): What the page calls it ("Generate sectors").
        steps (list[dict]): `{"label", "argv"}` per step, run in order; the generation commands among them run as the one
            step `planetgen run-job` (`job_steps.collapse`, PERF.203).
        admin (str, optional): Who started it, for the page.
        database (str, optional): Which database it writes, for the page.
        spawn (bool): Tests pass `False` to only write the rows, with no supervisor asked.
        draft (list, optional): The form fields (`[name, value]` pairs) the job was made from, kept so History can offer
            Run again with the same settings (ADM.115).

    Raises:
        JobBusy: Another job is still running.
        OSError: No writable jobs directory, no Redis server or control database.
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
            "created_at": time.time(), "cwd": REPO_DIR, "root": root, "draft": draft or [],
        }
        with _ledger() as control:
            pagerun.start(control, job, steps)
        with _redis_errors():
            connection = _connection()
            connection.hset(supervisor.pages_key(), job_id, time.time())       # a supervisor with a job to watch does not end (`supervisor.Store.busy`)
            if spawn:
                hand_over(connection, job_id)
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


def cancel_runs(job_id, who=None):
    """
    Cancels the live runs the job started, now (PERF.99, `planetgen.queue.runcontrol.cancel_run`): their waiting units
    are taken off the queue, their started units are sent the stop command and the machine is free, without waiting
    for the step that follows them to notice. Best effort: with no Redis the step's own stop is all there is.

    Returns:
        int: How many runs were canceled.
    """
    try:
        from planetgen.queue import held, runcontrol, runrecord

        connection = redisqueue.connect()
        return sum(runcontrol.cancel_run(connection, run_id, who=who) is not None
                   for run_id in runrecord.job_runs(connection, job_id) if not held.is_drained(connection, run_id))      # OPS.222: a run the update drained is not canceled with the job that followed it
    except Exception as exc:  # noqa: BLE001 -- the runner's own stop still follows
        log.debug("Cancel of job %s: its runs were not canceled (%s).", job_id, exc)
        return 0


def _end_step_group(path):
    """Ends the process group of the step `step.pid` names (SIGTERM, then SIGKILL): a step whose runner is gone was started apart from it and runs on (PERF.186). Returns whether there was one to end."""
    try:
        with open(os.path.join(path, job_runner.STEP_PID_NAME), "r", encoding="utf-8") as handle:
            step_pid = int(handle.read().strip() or 0)
    except (OSError, ValueError):
        return False
    did = False
    if step_pid > 1:
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
    with contextlib.suppress(OSError):
        os.remove(os.path.join(path, job_runner.STEP_PID_NAME))
    return did


def cancel_job(job_id, root=None, who=None):
    """
    Cancels a job that has not ended: a request to the supervisor (`controls.stop_engine_run`) that ends the job's rows and asks the process that runs its step to stop it, and stops
    the step and everything the step started (SIGTERM to its process group, then SIGKILL). With no supervisor, the rows are ended here and the step's group is ended from `step.pid`; the same
    when the job's process is gone (killed by the out-of-memory killer, a restart or a crash; PERF.186). The job's live runs are canceled through the stop command (`cancel_runs`), and its lock is freed.

    Returns:
        bool: Whether a running job was asked to stop.
    """
    root = root or jobs_dir()
    job = get_job(job_id, root)
    if job is None or (job["finished"] and job["status"] != "interrupted"):
        return False
    path = _job_dir(root, job_id)
    asked, runner = None, None
    try:
        with _ledger() as control:
            run = pagerun.run_of(control, job_id)
            runner = (pagerun.read(control, job_id) or {}).get("runner_pid")
        with _redis_errors():
            connection = _connection()
            asked = controls.stop_engine_run(connection, run, who=who)
            connection.hdel(supervisor.pages_key(), job_id)
    except OSError as exc:
        log.warning("Cancel of job %s: the ledger or the supervisor could not be reached: %s.", job_id, exc)
    ended = False if process.pid_alive(runner) else _end_step_group(path)          # a step whose runner is gone; a runner that is alive stops its step itself, between two looks at its row
    canceled = cancel_runs(job_id, who)    # after the stop: the step that finds its run stopped is then a canceled one, not a failed one
    release_lock(job_id)
    log.debug("Cancel requested for job %s.", job_id)
    return bool(asked) or ended or bool(canceled)


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


def submit(kind, title, steps, admin=None, database=None, root=None, draft=None):
    """
    Starts the job, or, when another job is running, puts it in the waiting list (PERF.154): it starts by itself when the
    job before it succeeded or was canceled and is held for an admin when that job failed.

    Returns:
        Submitted: The job started, or the place taken.

    Raises:
        runline.LineFull: The list holds `queue.max_waiting` entries; the message names them.
        updateguard.UpdateRunning, OSError: As `start_job`.
    """
    root = root or jobs_dir()
    try:
        return Submitted(job_id=start_job(kind, title, steps, admin=admin, database=database, root=root, draft=draft))
    except JobBusy as busy:
        running = busy.job["title"] if busy.job else "Another job"
    ticket = queue_next(kind, title, steps, admin=admin, database=database, root=root, draft=draft)
    started = promote_queued(root)      # the running job may have ended since the lock was refused
    waiting = [entry["ticket"] for entry in queued_jobs(root)]
    if started and ticket not in waiting:
        return Submitted(job_id=started)
    return Submitted(ticket=ticket, position=waiting.index(ticket) + 1 if ticket in waiting else None, running=running)


def queue_next(kind, title, steps, admin=None, database=None, root=None, draft=None):
    """
    Puts a job at the end of the waiting list: it starts by itself when the job before it succeeded or was canceled, and
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
    spec = {"kind": kind, "title": title, "admin": admin, "database": database, "steps": steps, "draft": draft or []}
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
            if after["status"] not in ("succeeded", "canceled"):
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
        job_id = start_job(spec["kind"], spec["title"], spec["steps"], admin=spec.get("admin"),
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
