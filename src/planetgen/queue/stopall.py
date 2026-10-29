# planetgen/queue/stopall.py

"""
Stop everything (PERF.186): the one way to end all work on the host and keep it from starting again until an admin says
so. `stop_all` writes the hold, cancels every live run and the Generate page's job and its waiting one, empties every
queue, and kills every planetgen worker, job runner and command process on this host together with the processes they
started. While the hold is there nothing starts: `refuse_while_stopped` (called by `updateguard.refuse_while_updating`,
which every start of a run, a job, a worker and an orbital update already calls) raises `Stopped`, `workers.ensure_workers`
starts no worker, the tile prefetch queues nothing, and a worker that is still alive exits before it takes a job.
`resume_all` removes the hold and the queue pause (`pause`).

The hold is a file in the state folder (`maintenance.state_dir`, `stopped.json`), not a flag in Redis, so it survives a
restart of Redis, the web server and the host, and is seen by the command line, the web server and the workers alike. The
Queue page and `planetgen queue stop-all` / `resume` are the controls.

Only processes this user may signal can be killed: run it as root (or as the web user when the work was started from the
site); what could not be signalled is reported, and the hold still keeps the rest from starting anything.
"""

import contextlib
import json
import os
import signal
import subprocess
import tempfile
import time
from datetime import datetime, timezone

from planetgen.admin import maintenance

FILE_NAME = "stopped.json"
"""str: The hold's file in the state folder: `{"by": name, "at": Unix seconds}` as JSON."""

RESUME_HINT = "planetgen queue resume"
"""str: The command that lifts the hold (the Queue page's Resume button does the same)."""

TERM_SECONDS = 5.0
"""float: How long processes asked to end (SIGTERM) get before they are killed (SIGKILL)."""

ROUNDS = 4
"""int: How many times the process list is looked at and what is found killed: a process that a dying parent started
in its last moment is found by the next look."""

ROOT_MODULES = ("planetgen.cli.worker", "planetgen.cli.generate", "planetgen.web.job_runner", "planetgen.cli.warm",
                "planetgen.cli.orbits", "planetgen.cli.smoke", "planetgen.cli.benchmark", "planetgen.cli.reset",
                "planetgen.cli.checkdb", "planetgen.cli.stage_progress")
"""tuple: The modules whose processes `stop_all` ends, with every process they started (their command line has `-m`
and the module). Not the maintenance, migrate and config commands, nor this one."""

ORPHAN_MARKERS = ("multiprocessing.spawn", "multiprocessing.forkserver")
"""tuple: The command lines of a process pool's children; one whose parent is gone (adopted by init) is a left-over of a
job that was killed, and is ended too."""


class Stopped(Exception):
    """Work was refused because everything was stopped; the message says what and how to resume."""


def _now():
    return time.time()


def path():
    """The hold's file."""
    return os.path.join(maintenance.state_dir(), FILE_NAME)


def state():
    """The hold as it stands: `{"stopped": bool, "by": name or None, "at": Unix seconds or None}`."""
    try:
        with open(path(), "r", encoding="utf-8") as handle:
            body = json.load(handle)
    except OSError:
        return {"stopped": False, "by": None, "at": None}
    except ValueError:
        body = {}       # a hold whose file is damaged still holds
    if not isinstance(body, dict):
        body = {}
    return {"stopped": True, "by": body.get("by"), "at": body.get("at")}


def is_stopped():
    """Whether the hold is on."""
    return state()["stopped"]


def _since(held):
    if not held.get("at"):
        return ""
    return " since " + datetime.fromtimestamp(held["at"], timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")


def refuse_while_stopped(what):
    """
    Raises `Stopped` naming `what` when the hold is on.

    Args:
        what (str): What was about to start ("A Generate job", "A worker").
    """
    held = state()
    if not held["stopped"]:
        return
    by = f" by {held['by']}" if held.get("by") else ""
    raise Stopped(f"{what} can't start: everything was stopped{by}{_since(held)}. Nothing starts until an admin resumes "
                  f"it (`{RESUME_HINT}`, or Resume on the admin Queue page).")


def _write_hold(who):
    target = path()
    os.makedirs(os.path.dirname(target) or ".", exist_ok=True)
    if os.path.exists(target):
        return False        # the first hold is the one kept
    fd, temp = tempfile.mkstemp(prefix=".stopped-", dir=os.path.dirname(target) or ".")
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        json.dump({"by": str(who or "")[:64] or None, "at": _now()}, handle)
        handle.write("\n")
    os.chmod(temp, 0o644)
    os.replace(temp, target)
    return True


# ---------------------------------------------------------------------
# The processes
# ---------------------------------------------------------------------

def process_table():
    """`{pid: (ppid, args)}` of every process on the host, `{}` where `ps` fails."""
    try:
        listing = subprocess.run(["ps", "-axww", "-o", "pid=,ppid=,args="], capture_output=True, text=True, timeout=20,
                                 check=False).stdout
    except (OSError, subprocess.SubprocessError):
        return {}
    table = {}
    for line in listing.splitlines():
        parts = line.split(None, 2)
        if len(parts) == 3 and parts[0].isdigit() and parts[1].isdigit():
            table[int(parts[0])] = (int(parts[1]), parts[2])
    return table


def _ancestors(table, pid):
    found = set()
    while pid and pid not in found:
        found.add(pid)
        pid = table.get(pid, (0, ""))[0]
    return found


def _is_root(args):
    words = args.split()
    for module in ROOT_MODULES:
        if module in words or f"-m{module}" in words:
            return True
    return False


def find_targets(table, protected=()):
    """
    The pids `stop_all` ends: every process running one of `ROOT_MODULES` or an orphaned pool child, and every
    descendant of those, less `protected` (this process and the processes that started it).

    Args:
        table (dict): `process_table`.
        protected (iterable): Pids never returned.

    Returns:
        list: Pids, parents before children.
    """
    roots = [pid for pid, (ppid, args) in table.items()
             if _is_root(args) or (ppid == 1 and any(marker in args for marker in ORPHAN_MARKERS))]
    children = {}
    for pid, (ppid, _args) in table.items():
        children.setdefault(ppid, []).append(pid)
    ordered, seen = [], set(protected)
    queue = sorted(roots)
    while queue:
        pid = queue.pop(0)
        if pid in seen:
            continue
        seen.add(pid)
        ordered.append(pid)
        queue.extend(sorted(children.get(pid, [])))
    return ordered


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    try:
        with open(f"/proc/{pid}/stat", "rb") as handle:
            return handle.read().rpartition(b")")[2].split()[0] != b"Z"
    except (OSError, IndexError):
        return True


def _signal(pids, number):
    """Sends `number` to each pid; returns the pids that could not be signalled for want of permission."""
    refused = []
    for pid in pids:
        try:
            os.kill(pid, number)
        except PermissionError:
            refused.append(pid)
        except OSError:
            pass
    return refused


def kill_processes(protected=None, term_seconds=None):
    """
    Ends every planetgen worker and command process on this host and what they started: SIGTERM, then SIGKILL for what is
    left after `term_seconds`, looking again up to `ROUNDS` times for what the dying ones started.

    Returns:
        dict: `{"killed": count, "survivors": [pid, ...]}`; a survivor is a process this user may not signal.
    """
    term_seconds = TERM_SECONDS if term_seconds is None else term_seconds
    mine = {os.getpid()}
    protected = set(protected or ()) | mine
    killed, survivors = set(), []
    for _round in range(ROUNDS):
        table = process_table()
        protected |= _ancestors(table, os.getpid())
        targets = find_targets(table, protected)
        if not targets:
            survivors = []
            break
        killed.update(targets)
        refused = set(_signal(targets, signal.SIGTERM))
        deadline = time.monotonic() + term_seconds
        while time.monotonic() < deadline and any(_alive(pid) for pid in targets if pid not in refused):
            time.sleep(0.1)
        left = [pid for pid in targets if pid not in refused and _alive(pid)]
        refused |= set(_signal(left, signal.SIGKILL))
        time.sleep(0.2)
        survivors = sorted(pid for pid in targets if _alive(pid))
        protected |= set(refused)       # a process this user may not signal is not tried again
    return {"killed": len(killed), "survivors": survivors}


# ---------------------------------------------------------------------
# Runs, queues and jobs
# ---------------------------------------------------------------------

def _cancel_runs(connection, who):
    from planetgen.queue import runcontrol, runrecord

    cancelled = 0
    for run_id in runrecord.recent_runs(connection, limit=100_000):
        try:
            if runcontrol.cancel_run(connection, run_id, who=who) is not None:
                cancelled += 1
        except Exception:       # noqa: BLE001 -- the next run, and the kill that follows, still go ahead
            pass
    return cancelled


def _empty_queues(connection):
    import rq

    from planetgen.queue import redisqueue

    removed = 0
    names = {queue.name for queue in rq.Queue.all(connection=connection)}
    names |= {redisqueue.BULK_QUEUE, redisqueue.MAINTENANCE_QUEUE, redisqueue.INTERACTIVE_QUEUE, redisqueue.LONG_QUEUE}
    for name in sorted(names):
        try:
            queue = redisqueue.queue(name, connection)
            removed += queue.count
            queue.empty()
        except Exception:       # noqa: BLE001 -- the other queues are still emptied
            pass
    return removed


def _clear_dead_workers(connection):
    """Takes the workers that were killed, and the jobs they left started, out of Redis."""
    import rq

    from planetgen.queue import workers

    with contextlib.suppress(Exception):
        for worker in rq.Worker.all(connection=connection):
            with contextlib.suppress(Exception):
                if not _alive(getattr(worker, "pid", 0) or 0):
                    worker.register_death()
    with contextlib.suppress(Exception):
        workers.reap_dead_jobs(connection)


def _job_roots():
    """The jobs folders that exist (looked for, never made: an update on a new machine must not create it as root)."""
    from planetgen.web import jobs

    roots = dict.fromkeys((jobs.configured_jobs_dir(), os.path.join(tempfile.gettempdir(), "planetgen-jobs")))
    return [root for root in roots if os.path.isdir(root)]


def _stop_jobs():
    """Drops the Generate page's waiting jobs and tells its running one to stop (the cancel file); the processes
    themselves are ended by `kill_processes`. Returns `(jobs stopped, folders)`."""
    from planetgen.web import jobs

    stopped, roots = 0, _job_roots()
    for root in roots:
        with contextlib.suppress(Exception):
            jobs.discard_all_queued(root)
        with contextlib.suppress(Exception):
            active = jobs.active_job(root)
            if active:
                with contextlib.suppress(OSError):
                    with open(os.path.join(root, active["id"], jobs.CANCEL_NAME), "w", encoding="utf-8") as handle:
                        handle.write(str(_now()))
                stopped += 1
    return stopped, roots


def _free_job_locks(roots):
    from planetgen.web import jobs

    if roots:
        with contextlib.suppress(OSError):
            jobs.free_lock()


def stop_all(who=None, connection=None, kill=True, hold=True):
    """
    Stops everything: the hold first (so nothing starts while the rest is done), then the Generate page's job, every live
    run, every queue, and last every planetgen process on the host.

    Args:
        who (str, optional): Who stopped it, for the hold, the run records and the audit.
        connection: The Redis client (default the configured server; `None` when none answers: the hold and the kill
            still happen).
        kill (bool): Tests pass `False` to leave the process table alone.
        hold (bool): Whether to leave the hold on. `update.sh` passes `False`: the update's own lock keeps new work
            out while it runs, and the site works again after it, but nothing the update left running survives it.

    Returns:
        dict: `{"already": the hold was on, "runs": cancelled, "queued": jobs removed from the queues, "job": a page
            job was stopped, "killed": processes ended, "survivors": pids this user could not signal}`.
    """
    from planetgen.queue import pause, redisqueue

    already = (not _write_hold(who)) if hold else is_stopped()
    if connection is None:
        with contextlib.suppress(Exception):
            connection = redisqueue.connect()
            connection.ping()
    summary = {"already": already, "runs": 0, "queued": 0, "job": 0, "killed": 0, "survivors": []}
    job_roots = []
    try:
        summary["job"], job_roots = _stop_jobs()
    except Exception:       # noqa: BLE001 -- the kill below stops the job's processes whatever this did
        pass
    if connection is not None:
        if hold:
            with contextlib.suppress(Exception):
                pause.pause(who or "stop everything", connection)
        summary["runs"] = _cancel_runs(connection, who)
        with contextlib.suppress(Exception):
            summary["queued"] = _empty_queues(connection)
    if kill:
        summary.update(kill_processes())
    if connection is not None:
        _clear_dead_workers(connection)
        with contextlib.suppress(Exception):
            summary["queued"] += _empty_queues(connection)      # what a dying process queued in its last moment
    _free_job_locks(job_roots)
    return summary


def resume_all(connection=None):
    """
    Lifts the hold and the queue pause, so work can start again; nothing is started by this.

    Returns:
        bool: Whether the hold was on.
    """
    from planetgen.queue import pause, redisqueue

    was = is_stopped()
    with contextlib.suppress(OSError):
        os.remove(path())
    try:
        pause.resume(connection or redisqueue.connect())
    except Exception:       # noqa: BLE001 -- no Redis, no pause to lift
        pass
    return was
