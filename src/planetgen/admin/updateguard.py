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
* **Running work is listed first and stopped only with the operator's word** (OPS.53, DB.29): every Generate job,
  run in the run record (PERF.96), queued job on the shared queues and long command-line run (the orbital update,
  the Galaxy Map warm-up) is printed; the update goes on only when the operator agrees to cancel them (`--yes`, or
  a typed y on the terminal), each is stopped through its own cancel path, and the update waits until none is
  left. A job that does not stop in time is reported and the update stops before it has changed anything. Nobody
  to ask is a no: a run started by systemd or cron (`started_unattended`) and one with no terminal refuse unless
  `--yes` was given, and a question left unanswered for `PROMPT_SECONDS` is a no too.
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
import re
import select
import signal
import sys
import tempfile
import time
from datetime import datetime, timezone

from planetgen.admin import maintenance
from planetgen.queue import stopall

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

PROMPT_SECONDS = 30
"""int: How long the question to cancel the running work waits for an answer; silence is a no (DB.29)."""

TIMER_PARENTS = ("cron", "crond", "anacron", "atd")
"""tuple: The names of the processes whose children have no person at the keyboard."""

UNIT_JOB = re.compile(r"^\d{8}-\d{6}-[0-9a-f]{8}-s\d+-u\d+$")
"""re.Pattern: The id of a queue job that is a unit of a run (`queue/units.py`'s `job_id`): the run lists it."""

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
    Raises `UpdateRunning` too when everything was stopped (`stopall`, PERF.186), saying so.

    Args:
        what (str): What was about to start ("A Generate job", "The orbital update").
    """
    try:
        stopall.refuse_while_stopped(what)      # PERF.186: stopped, then not even an update's own commands start work
    except stopall.Stopped as exc:
        raise UpdateRunning(str(exc)) from exc
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
def running(kind, title, **extra):
    """
    Registers this process as running `kind` (`"orbits"`, `"warm-map"`) for the guard to list and cancel, for the
    length of the block. A folder that can't be written only means the guard does not see this run. `extra` is kept in
    the marker (`admin/cachebuild.py` reads `refuse`).
    """
    path = None
    try:
        folder = running_dir()
        os.makedirs(folder, exist_ok=True)
        with contextlib.suppress(OSError):
            os.chmod(folder, 0o1777)    # every user a run starts as may register
        path = os.path.join(folder, f"{kind}-{os.getpid()}.json")
        with open(path, "w", encoding="utf-8") as handle:
            json.dump({"kind": kind, "title": title, "pid": os.getpid(), "since": _now(), **extra}, handle)
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
            try:
                job = jobs.active_job(root) if os.path.isdir(root) else None
            except OSError:         # no Redis server answers: the lock is kept there, so no job holds it
                job = None
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


class QueuedJobSource:
    """The jobs waiting or running on the shared queues (`queue/redisqueue.py`: bulk, interactive and
    generation). A run's own queue belongs to the run and is cancelled with it. Cancelled with
    `workers.cancel_now`. A job whose worker is gone is cleared before the listing (`workers.reap_dead_jobs`), so
    it does not hold the update for ever."""

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
        from planetgen.queue import workers
        workers.reap_dead_jobs(connection)      # a job whose worker was killed is not running
        found = []
        for name in redisqueue.QUEUES:
            queue = redisqueue.queue(name, connection)
            for state, ids in (("queued", queue.get_job_ids()),
                               ("running", queue.started_job_registry.get_job_ids())):
                for job_id in ids:
                    if UNIT_JOB.match(job_id):
                        continue        # a unit of a run: the run is listed (and cancelled) instead
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


class RunRecordSource:
    """The runs whose record in Redis says they are live (`queue/runrecord.py`, PERF.96): a `planetgen` run from
    the command line, one a Generate job started, a check or a warm-up. Cancelled through `runcontrol.cancel_run`
    (PERF.99), which takes the run's waiting units off the queue and stops the started ones."""

    name = "run-record"

    def __init__(self, url=None):
        self.url = url

    def _connection(self):
        from planetgen.queue import redisqueue
        connection = redisqueue.connect(self.url or redisqueue.redis_url())
        connection.ping()
        return connection

    def items(self):
        from planetgen.queue import runrecord
        try:
            connection = self._connection()
        except Exception:  # noqa: BLE001 -- no Redis server: no run records
            return []
        found = []
        for run_id in runrecord.recent_runs(connection, limit=100_000):
            run = runrecord.read_run(connection, run_id)
            if run is None or run["state"] not in runrecord.LIVE_STATES:
                continue
            started = run["started"] or run["created"]
            found.append({"source": self.name, "kind": run["kind"], "id": run_id, "title": run["title"],
                          "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started)) if started else None,
                          "holds": "the generation slot" if run["kind"] == runrecord.GENERATION else "a run slot",
                          "queued": False})
        return found

    def cancel(self, item):
        from planetgen.queue import runcontrol
        return runcontrol.cancel_run(self._connection(), item["id"], who="update") is not None


def default_sources(url=None):
    """The sources the guard looks at: Generate jobs, runs in the run record, queued jobs, registered runs."""
    return [GenerateJobSource(), RunRecordSource(url), QueuedJobSource(url), MarkerSource()]


# ---------------------------------------------------------------------
# Whether there is anybody to ask
# ---------------------------------------------------------------------

def _parent_names(pid=None):
    """The names of this process's ancestors, nearest first (Linux `/proc`; empty where there is none)."""
    names = []
    pid = os.getppid() if pid is None else pid
    while pid > 1 and len(names) < 32:
        try:
            with open(f"/proc/{pid}/stat", "rb") as handle:
                text = handle.read().decode("utf-8", "replace")
        except OSError:
            break
        names.append(text.partition("(")[2].rpartition(")")[0])
        fields = text.rpartition(")")[2].split()
        pid = int(fields[1]) if len(fields) > 1 else 0
    return names


def started_unattended(environ=None, parents=None, stdin_is_terminal=None):
    """
    Whether systemd or cron started this run, so nobody is at the keyboard (DB.29): it takes the answer to "cancel
    the running work?" as no unless `--yes` was given.

    systemd sets `JOURNAL_STREAM` for a service whose output goes to the journal, and `INVOCATION_ID` for every
    process of a unit (a terminal opened in a desktop session has that one too, and a terminal on its stdin);
    cron is the name of an ancestor process.

    Returns:
        str | None: `"systemd"`, `"cron"`, or `None` when a person may be there.
    """
    environ = os.environ if environ is None else environ
    if stdin_is_terminal is None:
        stdin_is_terminal = sys.stdin is not None and sys.stdin.isatty()
    if environ.get("JOURNAL_STREAM") or (environ.get("INVOCATION_ID") and not stdin_is_terminal):
        return "systemd"
    names = _parent_names() if parents is None else parents
    return "cron" if any(name in TIMER_PARENTS for name in names) else None


def answer_within(question, seconds=PROMPT_SECONDS, stdin=None):
    """
    Asks `question` on the terminal and waits `seconds` for a line.

    Returns:
        str | None: The reply, or `None` when none came in time (no answer is a no).
    """
    stdin = sys.stdin if stdin is None else stdin
    sys.stdout.write(question)
    sys.stdout.flush()
    ready, _, _ = select.select([stdin], [], [], seconds)
    if not ready:
        sys.stdout.write("\n")
        return None
    return stdin.readline()


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

def begin(pid, by="update.sh", sources=None, interactive=False, yes=False, unattended=None,
          answer=None, write=print, timeout=STOP_SECONDS, poll=POLL_SECONDS):
    """
    Starts an update: lists what runs and gets the operator's word to cancel it, takes the lock, closes the site and
    stops the work, before the first change.

    Args:
        pid (int): The process that holds the lock for as long as it lives (`update.sh`'s shell).
        sources: What to look at (`default_sources`).
        interactive (bool): Whether there is a terminal to ask on.
        yes (bool): The operator already agreed (`--yes`).
        unattended (str | None): Who started this run with nobody at the keyboard (`started_unattended`): it asks
            nothing, and with running work and no `yes` the answer is no.
        answer (callable): Asks a question and returns the reply, `None` for none in time (default
            `answer_within`: `PROMPT_SECONDS` of silence is a no).

    Returns:
        bool: Whether the update may go on. `False` leaves everything as it was (nothing cancelled that was not
        agreed, the site open, no lock).
    """
    sources = sources if sources is not None else default_sources()
    answer = answer_within if answer is None else answer
    existing = read_lock()
    if existing is not None and existing["pid"] != pid:
        write(f"error: an update is already running (process {existing['pid']}, since {existing['since']}).")
        return False
    items = list_running(sources)
    agreed = bool(yes)
    if items:
        write("The update needs the machine to itself. Running or queued work:")
        for line in describe(items):
            write(line)
        asks = interactive and not unattended
        if agreed:
            pass
        elif asks:
            reply = answer("Cancel these and go on with the update? [y/N] ")
            if reply is None:
                write(f"No answer in {PROMPT_SECONDS} seconds: taking that as no.")
            agreed = str(reply or "").strip().lower() in ("y", "yes")
        elif unattended:
            write(f"This run was started by {unattended}, with nobody to ask: taking that as no.")
        if not agreed:
            write("Not updating: nothing was changed and nothing was cancelled. "
                  + ("Run it again when the work has finished, or answer y to cancel it."
                     if asks else
                     "Run it again when the work has finished, or add --yes to cancel it."))
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
