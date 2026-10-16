# planetgen/admin/updateguard.py

"""
The update guard (OPS.41, OPS.53, OPS.54)
=========================================

While `update.sh` (or a manual `planetgen.cli.migrate` that has something to
migrate) changes the checkout and the databases, nothing else may run:

* **The site is closed** (OPS.41): maintenance mode (`admin/maintenance.py`) is turned on with the cause
  `UPDATE_CAUSE` before anything changes and turned off only when the update has finished. An update that fails
  leaves it on, with a message naming the step and the command that turns it off, because the database may be
  half migrated.
* **Running work is listed first and stopped only with the operator's word** (OPS.53): every Generate job, work
  queue run, queued job on the shared queues and long command-line run (the orbital update, the Galaxy Map
  warm-up) is printed; the update goes on only when the operator agrees to cancel them (a question on the
  terminal, or `--cancel-running` where there is none), each is stopped through its own cancel path, and the
  update waits until none is left. A job that does not stop in time is reported and the update stops before it has
  changed anything.
* **Nothing new starts** (OPS.53, OPS.54): the update takes a lock, `update.lock` in the state folder, and
  `refuse_while_updating` (called by the Generate page's job starter, the API's job queue, the work queue's runs,
  the queue's workers and the orbital update) refuses for as long as the lock is held by a live process or the
  maintenance state still has the cause `UPDATE_CAUSE`.

A run that is not a queue job registers itself in the state folder's `running` folder (`running`), which is how
the guard sees it. This module imports no Flask and, at import time, nothing outside the standard library and
`admin/maintenance.py`; the sources of running work import what they need when they are asked.
"""

import contextlib
import json
import os
import signal
import tempfile
import time
from datetime import datetime, timezone

from planetgen.admin import maintenance

UPDATE_CAUSE = "update"
"""str: The maintenance `cause` of an update."""

UPDATE_REASON = "The site is being updated and will be back in a moment."
"""str: What visitors are told while an update runs."""

EX_TEMPFAIL = 75
"""int: The exit code of a run that did nothing because an update is running (BSD `EX_TEMPFAIL`): a timer simply
tries again at its next firing."""

LOCK_ENV = "PLANETGEN_UPDATE_LOCK_FILE"
"""str: The environment variable that names the lock file."""

OWNER_ENV = "PLANETGEN_UPDATE_PID"
"""str: The process id of the update that runs a command (set by `update.sh`), so a `migrate` it starts knows the
update already holds the lock."""

LOCK_NAME = "update.lock"
RUNNING_DIR = "running"

STOP_SECONDS = 600
"""int: How long the update waits for the work it cancelled to end."""

POLL_SECONDS = 2.0
"""float: How often it looks again while it waits."""

TURN_OFF_HINT = "python -m planetgen.cli.maintenance off"
"""str: The command that opens the site again after a failed update."""


class UpdateRunning(OSError):
    """Work was refused because an update is running; the message says what and why."""


def _now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def pid_alive(pid):
    """Whether a process with this id is running here (one of another user counts)."""
    if not isinstance(pid, int) or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    try:
        with open(f"/proc/{pid}/stat", "rb") as handle:
            return handle.read().rpartition(b")")[2].split()[0] != b"Z"      # a zombie has ended
    except (OSError, IndexError):
        return True


def lock_path():
    """The lock file: `PLANETGEN_UPDATE_LOCK_FILE`, else `update.lock` in the state folder."""
    return os.environ.get(LOCK_ENV) or os.path.join(maintenance.state_dir(), LOCK_NAME)


def read_lock():
    """The update's lock (`{"pid", "since", "by", "kept"}`) while its process is alive, else `None`. The lock of a
    process that died is ignored (and replaced by the next update)."""
    try:
        with open(lock_path(), "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict) or not pid_alive(data.get("pid")):
        return None
    return data


def is_updating():
    """Whether an update runs: its lock is held by a live process, or maintenance mode still has the update's cause
    (an update that failed leaves the site closed, and nothing may run until an admin opens it)."""
    state = maintenance.read()
    return read_lock() is not None or (state["enabled"] and state["cause"] == UPDATE_CAUSE)


def inside_update():
    """Whether this process was started by the update that holds the lock (`OWNER_ENV`)."""
    lock = read_lock()
    try:
        return lock is not None and lock["pid"] == int(os.environ.get(OWNER_ENV, ""))
    except ValueError:
        return False


def refuse_while_updating(what):
    """
    Raises `UpdateRunning` naming `what` when an update runs. The update's own commands (`inside_update`) pass.

    Args:
        what (str): What was about to start ("A Generate job", "The orbital update").
    """
    if inside_update() or not is_updating():
        return
    raise UpdateRunning(f"{what} can't start: the site is being updated. Try again when the update has finished "
                        f"(if it failed, an admin opens the site again with `{TURN_OFF_HINT}`).")


def _write_lock(pid, by, kept):
    path = lock_path()
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    temp = f"{path}.{os.getpid()}.tmp"
    with open(temp, "w", encoding="utf-8") as handle:
        json.dump({"pid": pid, "since": _now(), "by": by, "kept": kept}, handle)
        handle.write("\n")
    os.chmod(temp, 0o644)
    os.replace(temp, path)


def _remove_lock():
    with contextlib.suppress(OSError):
        os.unlink(lock_path())


# ---------------------------------------------------------------------
# What is running
# ---------------------------------------------------------------------

def running_dir():
    """The folder the long command-line runs register in."""
    return os.path.join(maintenance.state_dir(), RUNNING_DIR)


@contextlib.contextmanager
def running(kind, title):
    """
    Registers this process as running `kind` (`"orbits"`, `"warm-map"`) for the guard to list and cancel, for the
    length of the block. A folder that can't be written only means the guard does not see this run.
    """
    path = None
    try:
        folder = running_dir()
        os.makedirs(folder, exist_ok=True)
        with contextlib.suppress(OSError):
            os.chmod(folder, 0o1777)    # every user a run starts as may register
        path = os.path.join(folder, f"{kind}-{os.getpid()}.json")
        with open(path, "w", encoding="utf-8") as handle:
            json.dump({"kind": kind, "title": title, "pid": os.getpid(), "since": _now()}, handle)
        os.chmod(path, 0o644)
    except OSError:
        path = None
    try:
        yield
    finally:
        if path:
            with contextlib.suppress(OSError):
                os.unlink(path)


class MarkerSource:
    """The registered command-line runs (`running`). Cancelled with SIGTERM."""

    name = "run"

    def items(self):
        found = []
        try:
            names = sorted(os.listdir(running_dir()))
        except OSError:
            return found
        for name in names:
            path = os.path.join(running_dir(), name)
            try:
                with open(path, "r", encoding="utf-8") as handle:
                    data = json.load(handle)
            except (OSError, ValueError):
                continue
            if not isinstance(data, dict):
                continue
            if not pid_alive(data.get("pid")):
                with contextlib.suppress(OSError):
                    os.unlink(path)
                continue
            found.append({"source": self.name, "kind": str(data.get("kind")), "id": str(data["pid"]),
                          "title": str(data.get("title") or data.get("kind")), "started": data.get("since"),
                          "holds": "a command-line run", "queued": False, "pid": data["pid"]})
        return found

    def cancel(self, item):
        with contextlib.suppress(OSError):
            os.kill(item["pid"], signal.SIGTERM)
        return True


class GenerateJobSource:
    """The Generate page's job (`web/jobs.py`), running or just started. Cancelled by its own cancel file. Looks
    where the job directory is without creating it: an update on a new machine must not make it as root."""

    name = "generate-job"

    def items(self):
        from planetgen.web import jobs
        found = []
        for root in dict.fromkeys((jobs.configured_jobs_dir(), os.path.join(tempfile.gettempdir(), "planetgen-jobs"))):
            job = jobs.active_job(root) if os.path.isdir(root) else None
            if job is None:
                continue
            started = (datetime.fromtimestamp(job["created_at"], timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
                       if job.get("created_at") else None)
            found.append({"source": self.name, "kind": job.get("kind") or "job", "id": job["id"],
                          "title": job.get("title") or "Generate job", "started": started,
                          "holds": "the Generate page's job lock", "queued": job["status"] == "starting",
                          "root": root})
        return found

    def cancel(self, item):
        from planetgen.web import jobs
        return jobs.cancel_job(item["id"], item["root"])


class WorkRunSource:
    """The work queue's live runs (`queue/work.py`) that no Generate job owns: a `planetgen` run from the command
    line. Cancelled through the run's control request, which it obeys between tasks."""

    name = "work-run"

    def __init__(self, config):
        self.config = config        # the control database's (`store.control_mysql_config`), or None

    def _connect(self):
        from planetgen.db import store
        return store.get_control_connection(self.config, ensure_schema=False)

    def items(self):
        if self.config is None:
            return []
        from planetgen.queue import work
        try:
            conn = self._connect()
        except Exception:  # noqa: BLE001 -- no control database yet (first install): no runs
            return []
        try:
            try:
                status = work.queue_status(conn)
                rows = conn.execute(
                    "SELECT id, kind, title, state, web_job_id, UNIX_TIMESTAMP(COALESCE(started_at, created_at)) AS t"
                    " FROM work_jobs WHERE parent_id IS NULL AND state IN ('waiting', 'running', 'paused')"
                    " AND heartbeat_at >= NOW(6) - INTERVAL ? SECOND", (work.STALE_SECONDS,)).fetchall()
            except Exception:  # noqa: BLE001 -- a control schema older than the job tree: no runs
                return []
            conn.rollback()
        finally:
            conn.close()
        found = []
        for row in rows:
            if row.get("web_job_id"):
                continue    # the Generate job that started it is listed (and cancelled) instead
            started = datetime.fromtimestamp(float(row["t"]), timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ") \
                if row.get("t") else None
            holds = "the work lease" if status.get("job_id") == row["id"] else "no lease yet"
            found.append({"source": self.name, "kind": row["kind"], "id": row["id"], "title": row["title"],
                          "started": started, "holds": holds, "queued": row["state"] == "waiting"})
        return found

    def cancel(self, item):
        from planetgen.queue import work
        conn = self._connect()
        try:
            return work.request_control(conn, item["id"], "cancel")
        finally:
            conn.close()


class QueuedJobSource:
    """The jobs waiting or running on the shared queues (`queue/redisqueue.py`: bulk, interactive and
    generation). A run's own queue belongs to the run and is cancelled with it. Cancelled with
    `workers.cancel_now`."""

    name = "queue-job"

    def __init__(self, url=None):
        self.url = url

    def _connection(self):
        from planetgen.queue import redisqueue
        connection = redisqueue.connect(self.url or redisqueue.redis_url())
        connection.ping()
        return connection

    def items(self):
        from planetgen.queue import redisqueue
        try:
            connection = self._connection()
        except Exception:  # noqa: BLE001 -- no Redis server: nothing is queued
            return []
        found = []
        for name in redisqueue.QUEUES:
            queue = redisqueue.queue(name, connection)
            for state, ids in (("queued", queue.get_job_ids()),
                               ("running", queue.started_job_registry.get_job_ids())):
                for job_id in ids:
                    job = redisqueue.fetch_job(job_id, connection, redisqueue.serializer_for(name))
                    started = job.started_at or job.enqueued_at if job is not None else None
                    found.append({"source": self.name, "kind": name, "id": job_id,
                                  "title": (job.description if job is not None else None) or job_id,
                                  "started": started.strftime("%Y-%m-%dT%H:%M:%SZ") if started else None,
                                  "holds": f"the {name} queue", "queued": state == "queued", "queue": name})
        return found

    def cancel(self, item):
        from planetgen.queue import redisqueue, workers
        workers.cancel_now(item["id"], self._connection(), redisqueue.serializer_for(item["queue"]))
        return True


def default_sources(config=None, url=None):
    """The sources the guard looks at: Generate jobs, work queue runs, queued jobs, registered runs."""
    control = None
    if config is not None:
        from planetgen.db import store
        control = store.control_mysql_config(config)
    return [GenerateJobSource(), WorkRunSource(control), QueuedJobSource(url), MarkerSource()]


def list_running(sources):
    """Every running or queued piece of work from `sources`, each as a dict (`source`, `kind`, `id`, `title`,
    `started`, `holds`, `queued`)."""
    found = []
    for source in sources:
        found.extend(source.items())
    return found


def describe(items):
    """The lines that list `items` for the operator."""
    lines = []
    for item in items:
        state = "queued" if item.get("queued") else "running"
        started = f", started {item['started']}" if item.get("started") else ""
        lines.append(f"  - {item['title']} ({item['kind']} {item['id']}, {state}{started}; holds {item['holds']})")
    return lines


def cancel_all(sources, items, write=print):
    """Asks each source to stop `items` through its own cancel path."""
    by_name = {source.name: source for source in sources}
    for item in items:
        try:
            by_name[item["source"]].cancel(item)
        except Exception as exc:  # noqa: BLE001 -- reported; the wait below decides whether it stopped
            write(f"  could not ask {item['title']} ({item['kind']} {item['id']}) to stop: {exc}")


def wait_until_idle(sources, timeout=STOP_SECONDS, poll=POLL_SECONDS, write=print, clock=time.monotonic,
                    sleep=time.sleep):
    """
    Cancels everything running and waits until nothing is, looking again every `poll` seconds (work that appeared
    meanwhile is cancelled too: only the first look asked the operator).

    Returns:
        list: What is still there after `timeout` seconds (empty when all stopped).
    """
    deadline = clock() + timeout
    asked = set()
    while True:
        items = list_running(sources)
        if not items:
            return []
        fresh = [item for item in items if (item["source"], item["id"]) not in asked]
        if fresh:
            cancel_all(sources, fresh, write)
            asked.update((item["source"], item["id"]) for item in fresh)
        if clock() >= deadline:
            return items
        sleep(poll)


# ---------------------------------------------------------------------
# Beginning and ending an update
# ---------------------------------------------------------------------

def begin(pid, by="update.sh", sources=None, config=None, interactive=False, cancel_running=False,
          answer=input, write=print, timeout=STOP_SECONDS, poll=POLL_SECONDS):
    """
    Starts an update: lists what runs and gets the operator's word to cancel it, takes the lock, closes the site and
    stops the work, before the first change.

    Args:
        pid (int): The process that holds the lock for as long as it lives (`update.sh`'s shell).
        sources: What to look at (`default_sources` for `config`).
        interactive (bool): Whether to ask on the terminal.
        cancel_running (bool): The operator already agreed (`--cancel-running`).
        answer (callable): Asks a question and returns the reply.

    Returns:
        bool: Whether the update may go on. `False` leaves everything as it was (nothing cancelled that was not
        agreed, the site open, no lock).
    """
    sources = sources if sources is not None else default_sources(config)
    existing = read_lock()
    if existing is not None and existing["pid"] != pid:
        write(f"error: an update is already running (process {existing['pid']}, since {existing['since']}).")
        return False
    items = list_running(sources)
    agreed = bool(cancel_running)
    if items:
        write("The update needs the machine to itself. Running or queued work:")
        for line in describe(items):
            write(line)
        if agreed:
            pass
        elif interactive:
            reply = answer("Cancel these and go on with the update? [y/N] ")
            agreed = str(reply).strip().lower() in ("y", "yes")
        if not agreed:
            write("Not updating: nothing was changed and nothing was cancelled. "
                  + ("Run it again when the work has finished, or answer y to cancel it."
                     if interactive else
                     "Run it again when the work has finished, or add --cancel-running to cancel it."))
            return False
    state = maintenance.read()
    already_closed = state["enabled"] and state["cause"] != UPDATE_CAUSE    # an admin's own closing stays theirs
    try:
        _write_lock(pid, by, kept=already_closed)
        if not already_closed:
            maintenance.turn_on(by=by, reason=UPDATE_REASON, cause=UPDATE_CAUSE)
    except OSError as exc:
        _remove_lock()
        write(f"error: could not close the site for the update ({exc}); nothing was changed.")
        return False
    # Whatever appeared between the first look and the lock is cancelled only with the same word.
    leftover = list_running(sources)
    if leftover and not agreed:
        _abort(already_closed)
        write("Not updating: work started while the update began:")
        for line in describe(leftover):
            write(line)
        write("Run it again when it has finished.")
        return False
    if leftover:
        write("Cancelling the work above and waiting for it to stop...")
        remaining = wait_until_idle(sources, timeout, poll, write)
        if remaining:
            _abort(already_closed)
            write(f"Not updating: this did not stop within {int(timeout)} seconds, and nothing was changed:")
            for line in describe(remaining):
                write(line)
            return False
        write("Everything stopped; the update has the machine to itself.")
    return True


def _abort(already_closed):
    """Undoes `begin`'s lock and, when the update closed the site, opens it again."""
    _remove_lock()
    if not already_closed:
        maintenance.turn_off_for(UPDATE_CAUSE, by="update")


def end(pid, ok, step=None, write=print):
    """
    Ends the update `pid` holds the lock for. A good run opens the site again (unless an admin had closed it
    before, or the sector-grid switch took the closing over: the Galaxy Map warm-up opens it then). A failed one
    leaves it closed, says which step failed and how to open it, and still releases the lock so the operator can run
    the update again.

    Returns:
        bool: Whether the site is open afterwards.
    """
    lock = read_lock()
    if lock is not None and lock["pid"] != pid:
        write(f"error: the update lock belongs to process {lock['pid']}, not {pid}; left alone.")
        return not maintenance.is_on()
    _remove_lock()
    state = maintenance.read()
    if ok:
        if lock is not None and not lock.get("kept"):
            maintenance.turn_off_for(UPDATE_CAUSE, by="update")
        return not maintenance.is_on()
    if state["enabled"] and state["cause"] == UPDATE_CAUSE:
        write(f"error: the update failed{f' at {step}' if step else ''}. The site stays closed (maintenance mode) "
              f"because the database may be half migrated. Fix the problem and run ./update.sh again; to open the "
              f"site without updating: {TURN_OFF_HINT}")
    else:
        write(f"error: the update failed{f' at {step}' if step else ''}.")
    return False
