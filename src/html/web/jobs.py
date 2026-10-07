# html/web/jobs.py

"""
Background jobs for the admin Generate page (`web/generate_page.py`):
resetting the galaxy, building its density skeleton (`generate.py plan`)
and generating sectors (`generate.py galaxy`), started from the browser
and left running after the request that started them has returned.

A job is a directory under the jobs directory (`jobs_dir()`):

    <jobs dir>/
        active                  the running job's id (the one-job lock)
        20260930-124433-1a2b/
            job.json            what to run: id, title, steps (argv lists), env
            runner.pid          the runner's pid, written when it is spawned
            cancel              written by Cancel; the runner stops when it sees it
            state.json          written by planetgen.cli.job as it goes
            progress.json       written by generate.py (planetgen.queue.progress_file)
            output.log          every step's stdout and stderr

`start_job` writes `job.json`, takes the lock and spawns
`python3 -m planetgen.cli.job <job dir>` in its own session (on Windows, a
detached process in its own process group), detached from the web server (a mod_wsgi request timeout or a graceful Apache reload
doesn't stop it). The page then only ever reads these files, so it works
the same whichever server process answers the next request.

Only one job runs at a time: generation, planning and reset all write the
same database. A lock whose runner has died (the server was rebooted mid
job) is noticed and cleared by the next `start_job`, and that job shows
as "interrupted".

The database a job writes is the one this site shows (`helpers.db_name`),
passed to the child through the same `PLANETGEN_MYSQL_*` environment
variables every command-line tool reads, never on the command line.
"""

import json
import os
import re
import secrets
import shutil
import subprocess
import sys
import tempfile
import threading
import time

from planetgen.web.lib.privatedir import ensure_private_dir
from planetgen.util import log
from planetgen.util.appconfig import load_config

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
"""str: The planetGen checkout (`generate.py` lives here)."""

GENERATE_SCRIPT = os.path.join(REPO_DIR, "generate.py")
RESET_COMMAND = ["-m", "planetgen.cli.reset"]
"""list: Runs `planetgen.cli.reset` (the job runner puts the checkout's
`src/` on the steps' `PYTHONPATH`)."""
SRC_DIR = os.path.join(REPO_DIR, "src")
"""str: The checkout's `src/`, where the job runner starts so that
`python -m planetgen.cli.job` finds the planetgen package."""

DEFAULT_JOBS_DIR = "/var/lib/planetgen/jobs"
"""str: Used when neither `PLANETGEN_JOBS_DIR` nor `config.json`'s
`jobs.dir` names a directory. Falls back to a `planetgen-jobs` folder in
the system temp directory when Apache can't create this one."""

DEFAULT_KEEP = 20
"""int: How many finished jobs (and their logs) are kept; older ones are
deleted when a new job starts."""

LOCK_NAME = "active"
"""str: Must match `planetgen.cli.job.LOCK_NAME`."""

CANCEL_NAME = "cancel"
"""str: Must match `planetgen.cli.job.CANCEL_NAME`."""

WINDOWS = os.name == "nt"

STARTING_GRACE_SECONDS = 15
"""int: How long a just-spawned job may go without a `state.json` before
it counts as interrupted (the runner writes one as its first act)."""

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
    return load_config().get("jobs") or {}


def configured_jobs_dir():
    """`PLANETGEN_JOBS_DIR`, then `jobs.dir`, then `DEFAULT_JOBS_DIR`.
    `examples/apache/create-cache-dir.sh` asks this where to create it."""
    return os.environ.get("PLANETGEN_JOBS_DIR") or _config().get("dir") or DEFAULT_JOBS_DIR


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
    try:
        return max(1, int(_config().get("keep") or DEFAULT_KEEP))
    except (TypeError, ValueError):
        return DEFAULT_KEEP


def python_executable():
    """
    The Python that runs jobs: `PLANETGEN_PYTHON`, then `jobs.python`,
    then this process's own interpreter. Under mod_wsgi `sys.executable`
    can be Apache itself rather than Python, so that case falls back to
    `<sys.prefix>/bin/python3`, then `python3` on the PATH.
    """
    configured = os.environ.get("PLANETGEN_PYTHON") or _config().get("python")
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


def _windows_process_alive(pid, created_by=None):
    """
    Whether process `pid` is running, on Windows (`os.kill(pid, 0)` there
    sends CTRL_C_EVENT instead of checking). With `created_by` (a Unix
    time), a process created after it doesn't count: that pid was reused.
    """
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel32.GetExitCodeProcess.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
    kernel32.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(wintypes.FILETIME)] * 4
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    process_query_limited_information = 0x1000
    still_active = 259
    handle = kernel32.OpenProcess(process_query_limited_information, False, int(pid))
    if not handle:
        return False
    try:
        code = wintypes.DWORD()
        if not kernel32.GetExitCodeProcess(handle, ctypes.byref(code)) or code.value != still_active:
            return False
        if created_by is not None:
            times = [wintypes.FILETIME() for _ in range(4)]
            if kernel32.GetProcessTimes(handle, *(ctypes.byref(t) for t in times)):
                ticks = (times[0].dwHighDateTime << 32) | times[0].dwLowDateTime
                created = ticks / 1e7 - 11644473600  # 100 ns ticks since 1601 -> Unix time
                if created > created_by:
                    return False
        return True
    finally:
        kernel32.CloseHandle(handle)


def _runner_alive(pid, job_id, created_at=None):
    """Whether `pid` is still this job's runner. Reads the process's
    command line where `/proc` exists, and on Windows checks it was
    created within `STARTING_GRACE_SECONDS` of the job (`created_at`),
    so a reused pid doesn't count."""
    if not pid:
        return False
    if WINDOWS:
        try:
            created_by = created_at + STARTING_GRACE_SECONDS if created_at else None
            return _windows_process_alive(int(pid), created_by)
        except (OSError, ValueError):
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
            `step_label`, `started_at`, `finished_at`, `elapsed_s`,
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
        pid = _runner_pid(path, state)
        if now - job.get("created_at", 0) < STARTING_GRACE_SECONDS and (
                pid is None or _runner_alive(pid, job_id, job.get("created_at"))):
            status = "starting"
        else:
            status = "interrupted"
    elif status == "running" and not _runner_alive(state.get("pid"), job_id, job.get("created_at")):
        status = "interrupted"

    labels = [step["label"] for step in job.get("steps", [])]
    step = state.get("step") or 0
    started = state.get("started_at")
    finished = state.get("finished_at")
    progress = _read_json(os.path.join(path, "progress.json")) if step else None
    error = state.get("error")
    if status == "interrupted" and not error:
        error = "The job stopped without finishing (the server may have restarted)."
    return {
        "id": job["id"],
        "kind": job.get("kind"),
        "title": job.get("title"),
        "admin": job.get("admin"),
        "database": job.get("database"),
        "created_at": job.get("created_at"),
        "steps": labels,
        "status": status,
        "finished": status in FINISHED,
        "step": step,
        "step_label": labels[step - 1] if 0 < step <= len(labels) else None,
        "started_at": started,
        "finished_at": finished,
        "elapsed_s": ((finished or now) - started) if started else None,
        "error": error,
        "progress": progress,
    }


def list_jobs(limit=10, root=None):
    """The newest `limit` jobs, newest first (`get_job` dicts)."""
    root = root or jobs_dir()
    try:
        names = sorted((name for name in os.listdir(root) if JOB_ID_RE.match(name)), reverse=True)
    except OSError:
        return []
    jobs = []
    for name in names:
        job = get_job(name, root)
        if job is not None:
            jobs.append(job)
        if len(jobs) >= limit:
            break
    return jobs


def active_job(root=None):
    """The running job (`get_job` dict), or `None`."""
    root = root or jobs_dir()
    try:
        with open(os.path.join(root, LOCK_NAME), "r", encoding="utf-8") as f:
            holder = f.read().strip()
    except OSError:
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


# ---------------------------------------------------------------------
# Starting and stopping jobs
# ---------------------------------------------------------------------

CLEARING_NAME = "active.clearing"
"""str: Held (created exclusively) by whoever is clearing a stale lock,
so two admins starting a job at once can't both clear it, the second
deleting the lock the first has just taken (TEST.40)."""


def _lock_age(path):
    try:
        return time.time() - os.path.getmtime(path)
    except OSError:
        return None


def _write_lock(lock, job_id):
    """Creates `lock` holding `job_id` in one step, so nobody ever reads
    it empty: the id goes into a temporary file that is then hard-linked
    to the lock's name (which fails if the lock exists). Where hard links
    aren't available, creates it exclusively and then writes it.

    Raises:
        FileExistsError: The lock exists.
    """
    directory = os.path.dirname(lock)
    fd, tmp = tempfile.mkstemp(prefix=".active-", dir=directory)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(job_id)
        try:
            os.link(tmp, lock)
            return
        except FileExistsError:
            raise
        except (OSError, NotImplementedError):
            pass
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass
    fd = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        f.write(job_id)


def _read_lock(lock):
    try:
        with open(lock, "r", encoding="utf-8") as f:
            return f.read(256).strip()
    except OSError:
        return None


def _clear_stale_lock(root, lock, holder):
    """Removes `lock` if it still holds `holder` (judged stale), under
    `CLEARING_NAME`. Returns whether the caller may try again now."""
    clearing = os.path.join(root, CLEARING_NAME)
    try:
        fd = os.open(clearing, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    except FileExistsError:
        age = _lock_age(clearing)
        if age is not None and age > STARTING_GRACE_SECONDS:
            # Left by a process that died while clearing.
            try:
                os.remove(clearing)
            except OSError:
                pass
        return False
    os.close(fd)
    try:
        if _read_lock(lock) == holder:
            try:
                os.remove(lock)
            except OSError:
                pass
        return True
    finally:
        try:
            os.remove(clearing)
        except OSError:
            pass


def _take_lock(root, job_id):
    """
    Creates the `active` lock for `job_id`, clearing a stale one. A lock
    counts as stale when the job it names has finished, or is unknown
    (garbage, an id whose `job.json` is missing or unreadable) and the
    lock is older than `STARTING_GRACE_SECONDS` -- a younger one may
    belong to a job still being written.

    Raises:
        JobBusy: A live job (or one still starting) holds it.
    """
    lock = os.path.join(root, LOCK_NAME)
    for _attempt in range(3):
        try:
            _write_lock(lock, job_id)
            return
        except FileExistsError:
            pass
        holder = _read_lock(lock)
        if holder is None:
            continue  # gone meanwhile
        job = get_job(holder, root) if JOB_ID_RE.match(holder) else None
        if job is not None and not job["finished"]:
            raise JobBusy(job)
        if job is None:
            age = _lock_age(lock)
            if age is None:
                continue
            if age < STARTING_GRACE_SECONDS:
                raise JobBusy(None)
        if not _clear_stale_lock(root, lock, holder):
            raise JobBusy(active_job(root))
    raise JobBusy(active_job(root))


def _release_lock(root, job_id):
    lock = os.path.join(root, LOCK_NAME)
    try:
        with open(lock, "r", encoding="utf-8") as f:
            if f.read().strip() == job_id:
                os.remove(lock)
    except OSError:
        pass


def _prune(root, keep):
    """Deletes all but the newest `keep` jobs' directories. Never the
    one the lock names, one still running or starting, or one whose
    `job.json` can't be read yet less than `STARTING_GRACE_SECONDS` old."""
    names = sorted((name for name in os.listdir(root) if JOB_ID_RE.match(name)), reverse=True)
    holder = _read_lock(os.path.join(root, LOCK_NAME))
    for name in names[keep:]:
        if name == holder:
            continue
        job = get_job(name, root)
        if job is not None and not job["finished"]:
            continue
        if job is None:
            age = _lock_age(os.path.join(root, name))
            if age is None or age < STARTING_GRACE_SECONDS:
                continue
        shutil.rmtree(os.path.join(root, name), ignore_errors=True)


def new_job_id():
    return time.strftime("%Y%m%d-%H%M%S") + "-" + secrets.token_hex(2)


def start_job(kind, title, steps, env=None, admin=None, database=None, root=None, spawn=True):
    """
    Starts a job in the background and returns its id straight away.

    Args:
        kind (str): `new_galaxy`, `plan`, `galaxy` or `reset`.
        title (str): What the page calls it ("Generate sectors").
        steps (list[dict]): `{"label", "argv"}` per step, run in order.
        env (dict, optional): Extra environment for every step (the
            `mysql_env` variables).
        admin (str, optional): Who started it, for the page.
        database (str, optional): Which database it writes, for the page.
        spawn (bool): Tests pass `False` to only write the files.

    Raises:
        JobBusy: Another job is still running.
        OSError: No writable jobs directory, or the runner can't start.
    """
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
            "created_at": time.time(), "cwd": REPO_DIR, "steps": steps, "env": env or {},
        }
        fd = os.open(os.path.join(path, "job.json"), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(job, f)
        if spawn:
            _spawn(path)
    except BaseException:
        _release_lock(root, job_id)
        raise
    log.debug("Started job %s (%s) for %s on %s: %s", job_id, kind, admin, database,
              [step["label"] for step in steps])
    try:
        _prune(root, keep_count())
    except OSError:
        pass
    return job_id


def _detached_options():
    """
    Popen options that detach the runner from the web server, most
    detached first: its own session on POSIX; on Windows (where
    start_new_session is ignored) its own process group with no console,
    first also broken away from the server's job object (IIS and some
    service wrappers kill a job object's processes on a recycle), then
    without that when the job object doesn't allow it.
    """
    if not WINDOWS:
        return [{"start_new_session": True}]
    flags = subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS | subprocess.CREATE_NO_WINDOW
    return [{"creationflags": flags | subprocess.CREATE_BREAKAWAY_FROM_JOB}, {"creationflags": flags}]


def _spawn(path):
    options = _detached_options()
    for index, extra in enumerate(options):
        try:
            proc = subprocess.Popen(
                [python_executable(), "-m", "planetgen.cli.job", path],
                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                cwd=SRC_DIR, close_fds=True, **extra,
            )
            break
        except PermissionError:  # breakaway refused by the job object
            if index == len(options) - 1:
                raise
    with open(os.path.join(path, "runner.pid"), "w", encoding="utf-8") as f:
        f.write(str(proc.pid))
    # Reap it when it exits, so it doesn't linger as a zombie of the web
    # server process. If this process goes first, init adopts it instead.
    threading.Thread(target=proc.wait, name=f"planetgen-job-{proc.pid}", daemon=True).start()


def cancel_job(job_id, root=None):
    """
    Asks a running job to stop: writes a `cancel` file into its
    directory, which the runner checks for while a step runs. It then
    stops the step and everything the step started, and marks the job
    cancelled. (A file rather than a signal: on Windows `os.kill` is
    TerminateProcess, which would kill the runner and leave its step
    running and the lock taken.)

    Returns:
        bool: Whether a running job was asked to stop.
    """
    root = root or jobs_dir()
    job = get_job(job_id, root)
    if job is None or job["finished"]:
        return False
    path = _job_dir(root, job_id)
    pid = _runner_pid(path, _read_json(os.path.join(path, "state.json")))
    if not _runner_alive(pid, job_id, job.get("created_at")):
        return False
    try:
        with open(os.path.join(path, CANCEL_NAME), "w", encoding="utf-8") as f:
            f.write(str(time.time()))
    except OSError:
        return False
    log.debug("Cancel requested for job %s (runner pid %s)", job_id, pid)
    return True
