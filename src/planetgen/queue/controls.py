# planetgen/queue/controls.py

"""
Controls through the supervisor (OPS.234, D9 of `/mnt/project-files/reports/job-engine-items-stage-3-6-2026-11-13.md`): pause, resume, Stop All, cancel and clear are requests to the
supervisor, which knows every job process it forked, instead of a scan of the host's process list for the processes that look like a worker.

**The request** is a JSON line on a Redis list (`supervisor.requests_key`): `{"id", "op", "who", "at", "arguments"}`. `ask` posts one and, when the caller wants the answer, waits for the reply the
supervisor pushes on a list of its own (`reply_key`). A supervisor that is not alive is never started for a request: with none alive no job process of the ledger is running (a job whose
supervisor is gone finishes its piece and exits, `template.EXIT_SUPERVISOR_GONE`), so `ask` gives `None` and the caller does what the control always did, and `end_rows` settles the ledger itself.

**The duty** (`Controls`) takes the requests at each turn of the supervisor, before the launcher's own turn:

* `pause` and `resume` write and clear the flag the engine reads between two pieces (`pause`). While it is set the launcher's hand-outs are held (`supervisor.Launches`): jobs that run finish their
  piece, nothing new is forked.
* `stop_jobs` (a canceled or cleared run) and `stop_all` end the ledger's jobs: the rows still waiting or ready are canceled at once, each job process that runs is asked to stop between two pieces
  through the pipe the supervisor holds to it (`{"op": "stop"}`, the pieces done so far kept), and one still alive after `GRACE_SECONDS` is killed. The run's row and the rows of its jobs end
  `canceled`; an attempt that was running is closed `lost` with the kind `canceled`. Stop All's reply comes when the last process is gone, so the process scan that follows it (`stopall`, kept as the
  second line) has nothing of the engine's to find.

* `stop_unit` (the stuck-unit follower, `unitstuck`) kills the one job process that runs a unit the judge found stuck, and closes its attempt `stuck` with the reason and the stack the follower had it
  write (a file the attempt's row names); the launcher settles the process's end as it settles any other, the row going back to `ready` for one more try.

Nothing here decides which run to cancel or how its Redis record ends: that stays with `runcontrol`, which asks this for the ledger's half.
"""

import contextlib
import json
import time
import uuid

from planetgen.queue import attemptlog, ledger, planrows, runrecord, supervisor
from planetgen.util import log

OPS = ("pause", "resume", "stop_jobs", "stop_all", "stop_unit")
"""tuple: What a request can ask."""

ASK_SECONDS = 30.0
"""float: How long a caller that wants the answer waits for it (a Stop All waits for every job to be gone, up to `GRACE_SECONDS` and a turn or two)."""

GRACE_SECONDS = 10.0
"""float: How long a job process asked to stop between two pieces has before it is killed."""

REPLY_SECONDS = 120
"""int: How long a reply is kept for a caller that has gone."""

PER_TURN = 20
"""int: The most requests taken at one turn."""

STOPPED = "the run was canceled"
"""str: The text of the attempt a stop closes."""

TERMINAL = ("done", "failed", "canceled")
"""tuple: The states of a job that is over."""


def reply_key(request_id):
    """The list the supervisor pushes the reply to request `request_id` on."""
    return runrecord._k("supervisor", "reply", request_id)


def ask(connection, op, who=None, wait=ASK_SECONDS, **arguments):
    """
    Leaves a request for the supervisor.

    Args:
        connection (redis.Redis): The server the supervisor watches.
        op (str): One of `OPS`.
        who (str, optional): Who asked, for the supervisor's list of acts and the audit.
        wait (float): How long to wait for the reply; `0` leaves the request and returns at once.
        **arguments: What the op needs (`run_id` for `stop_jobs`, `state` for `stop_all`).

    Returns:
        dict | None: The reply (`{"ok": bool, ...}`); `{"ok": True, "pending": True}` when the request was left and not answered in `wait`; `None` when it was not left (no supervisor is alive, or
        Redis did not answer), so the caller does the work itself.
    """
    if op not in OPS:
        raise ValueError(f"No request is {op!r}: {OPS}.")
    try:
        if not supervisor.is_alive(connection):
            return None
        request = {"id": uuid.uuid4().hex, "op": op, "who": who, "at": time.time(), "arguments": arguments}
        raw = json.dumps(request)
        connection.rpush(supervisor.requests_key(), raw)
        connection.expire(supervisor.requests_key(), supervisor.RECENT_SECONDS)
        if not wait:
            return {"ok": True, "pending": True, "id": request["id"]}
        got = connection.blpop(reply_key(request["id"]), timeout=max(1, int(wait)))
        if got is None:
            if connection.lrem(supervisor.requests_key(), 1, raw):
                return None             # nobody took it: the caller's own path is the one that runs
            return {"ok": True, "pending": True, "id": request["id"]}
        return json.loads(runrecord._s(got[1]))
    except Exception as exc:  # noqa: BLE001 -- a control that cannot reach the supervisor does what it always did
        log.swallowed("queue/controls", "ask", exc)
        return None


# -- the ledger's half ---------------------------------------------------------------------------------------------------------

def engine_run_of(connection, run_id):
    """The `engine_runs` id of the Redis run `run_id`, or `None` for a run with no rows in the ledger."""
    found = connection.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)
    return int(runrecord._s(found)) if found else None


def end_rows(control, engine_run=None, run_state="canceled"):
    """
    Ends the rows of the ledger run `engine_run` (every run still going, for `None`): each job not over becomes `canceled`, each attempt still `running` is closed `lost`, and the run's row takes
    `run_state`. Does not stop a process (the supervisor's duty asks them to; with none alive there is none). The caller commits.

    Returns:
        dict: `{"rows": jobs canceled, "attempts": attempts closed, "runs": run rows ended}`.
    """
    where, params = ("run_id = ?", (engine_run,)) if engine_run is not None else ("1 = 1", ())
    running = [row["id"] for row in control.execute(f"SELECT id FROM engine_jobs WHERE {where} AND state = 'running'", params).fetchall()]
    attempts = sum(len(ledger.end_running(control, job_id, "lost", error_kind="canceled", error_text=STOPPED)) for job_id in running)
    rows = control.execute(f"UPDATE engine_jobs SET state = 'canceled', waited = NULL, updated_at = NOW(6) WHERE {where} AND state NOT IN ('done', 'failed', 'canceled')", params).rowcount
    runs = control.execute(f"UPDATE engine_runs SET state = ?, held_reason = NULL, finished_at = NOW(6) WHERE {'id = ?' if engine_run is not None else '1 = 1'} AND state NOT IN ('done', 'canceled')",
                           (run_state, *((engine_run,) if engine_run is not None else ()))).rowcount
    return {"rows": rows, "attempts": attempts, "runs": runs}


def stop_ledger_run(connection, run_id):
    """
    The ledger's half of a cancel or a clear done without a supervisor: the rows of the Redis run `run_id` ended (`end_rows`). Returns what it ended, or `None` for a run with no rows.
    """
    engine = engine_run_of(connection, run_id)
    if engine is None:
        return None
    control = planrows.control_connection(connection, run_id)
    try:
        ended = end_rows(control, engine)
        control.commit()
        return ended
    finally:
        control.close()


def stop_run(connection, run_id, who=None):
    """
    The ledger's half of a cancel or a clear (`runcontrol`): a request to the supervisor when one is alive (left, not waited for: the supervisor ends the rows and stops the jobs at its next turn),
    else the rows ended here. A run with no rows in the ledger is nothing to stop. Never raises: the run's Redis record is what the cancel has ended, the rows are its history.

    Returns:
        dict | None: What was ended or asked, or `None` for a run that is not in the ledger.
    """
    try:
        if engine_run_of(connection, run_id) is None:
            return None
        asked = ask(connection, "stop_jobs", who, wait=0, run_id=run_id)
        return asked if asked is not None else stop_ledger_run(connection, run_id)
    except Exception as exc:  # noqa: BLE001 -- see above
        log.swallowed("queue/controls", "stop_run", exc)
        return None


def running_jobs(connection, run_id):
    """How many jobs of the Redis run `run_id` the ledger still has `running` (a cancel's wait for the processes to end); 0 for a run with no rows there or when the ledger cannot be read."""
    try:
        engine = engine_run_of(connection, run_id)
        if engine is None:
            return 0
        control = planrows.control_connection(connection, run_id)
        try:
            return int(control.execute("SELECT COUNT(*) AS n FROM engine_jobs WHERE run_id = ? AND state = 'running'", (engine,)).fetchone()["n"])
        finally:
            control.close()
    except Exception as exc:  # noqa: BLE001 -- the wait is a courtesy
        log.swallowed("queue/controls", "running_jobs", exc)
        return 0


def stop_everything(connection, who=None):
    """
    Stop All's half for the ledger: the supervisor ends every job it forked and waits until they are gone (the reply says how many), else the rows are ended here.

    Returns:
        dict: `{"jobs": job processes stopped, "rows": rows canceled}` (zeros when there was nothing).
    """
    asked = ask(connection, "stop_all", who)
    if asked is not None and not asked.get("pending"):
        return {"jobs": int(asked.get("jobs", 0)), "rows": int(asked.get("rows", 0))}
    if asked is not None:
        return {"jobs": 0, "rows": 0}           # the supervisor has it and is still at it; the scan that follows is the second line
    done = {"jobs": 0, "rows": 0}
    with contextlib.suppress(Exception):
        from planetgen.db import store
        control = store.get_control_connection(store.control_mysql_config())
        try:
            done["rows"] = end_rows(control)["rows"]
            control.commit()
        finally:
            control.close()
    return done


# -- the supervisor's duty -----------------------------------------------------------------------------------------------------

class Pending:
    """Job processes asked to stop that have not been seen end: `ids` are their job-row ids, `engine` the ledger run they belong to (`None` for all), `reply` the request to answer when they are
    gone (or `None`), `summary` what the reply says."""

    def __init__(self, ids, deadline, engine, reply, summary):
        self.ids, self.deadline, self.engine, self.reply, self.summary = set(ids), deadline, engine, reply, summary
        self.killed = False


class Controls:
    """
    The supervisor's duty for the requests of `OPS` (see the module). `launches` is the `supervisor.Launches` duty whose launcher holds the jobs forked; `control` returns a connection to the control
    database for the moments there is no launcher (the setting that forks jobs off, or no job yet), which this duty closes after use.
    """

    def __init__(self, connection, launches=None, clock=time.time, grace=GRACE_SECONDS, control=None):
        self.connection, self.launches, self.clock, self.grace = connection, launches, clock, grace
        self._control = control
        self.pending = []
        self.supervisor = None

    def attach(self, running):
        self.supervisor = running

    def _act(self, what, reason=""):
        if self.supervisor is not None:
            self.supervisor.act(what, reason)

    def _launcher(self):
        return self.launches.launcher if self.launches is not None else None

    @contextlib.contextmanager
    def _control_connection(self, run_id=None):
        """The launcher's control connection, else one of its own (closed after)."""
        launcher = self._launcher()
        if launcher is not None:
            yield launcher.control
            return
        if self._control is not None:
            made = self._control()
        elif run_id is not None:
            made = planrows.control_connection(self.connection, run_id)
        else:
            from planetgen.db import store
            made = store.get_control_connection(store.control_mysql_config())
        try:
            yield made
        finally:
            made.close()

    # -- one turn ----------------------------------------------------------------------------------------------------------------

    def turn(self):
        """Takes the requests left, then looks at the jobs it asked to stop. Returns how many requests it took."""
        taken = 0
        for _ in range(PER_TURN):
            raw = self.connection.lpop(supervisor.requests_key())
            if raw is None:
                break
            taken += 1
            self._take(runrecord._s(raw))
        self._settle()
        return taken

    def _answer(self, request, body):
        if request is None or not request.get("id"):
            return
        key = reply_key(request["id"])
        pipe = self.connection.pipeline()
        pipe.rpush(key, json.dumps(body))
        pipe.expire(key, REPLY_SECONDS)
        pipe.execute()

    def _take(self, raw):
        try:
            request = json.loads(raw)
            op = request["op"]
        except (ValueError, KeyError, TypeError) as exc:
            log.swallowed("queue/controls", "_take", exc)
            return
        try:
            handler = {"pause": self._pause, "resume": self._resume, "stop_jobs": self._stop_jobs, "stop_all": self._stop_all,
                       "stop_unit": self._stop_unit}[op]
            body = handler(request, request.get("arguments") or {}, request.get("who"))
        except Exception as exc:  # noqa: BLE001 -- a request that failed is answered, and the next one is taken
            log.swallowed("queue/controls", f"the request {op!r}", exc)
            body = {"ok": False, "error": f"{type(exc).__name__}: {exc}"[:300]}
        if body is not None:
            self._answer(request, body)

    # -- the ops ------------------------------------------------------------------------------------------------------------------

    def _pause(self, request, _arguments, who):
        from planetgen.queue import pause
        ok = pause.pause(who or "the supervisor", self.connection)
        if ok:
            self._act("paused the queues", f"asked by {who}" if who else "asked")
        return {"ok": bool(ok)}

    def _resume(self, request, _arguments, who):
        from planetgen.queue import pause
        ok = pause.resume(self.connection)
        if ok:
            self._act("resumed the queues", f"asked by {who}" if who else "asked")
        return {"ok": bool(ok)}

    def _stop_jobs(self, request, arguments, who):
        run_id = arguments["run_id"]
        engine = engine_run_of(self.connection, run_id)
        if engine is None:
            return {"ok": True, "jobs": 0, "rows": 0}
        with self._control_connection(run_id) as control:
            return self._end(request, control, engine, who, f"run {run_id}", wait_for_reply=False)

    def _stop_all(self, request, _arguments, who):
        with self._control_connection() as control:
            return self._end(request, control, None, who, "everything", wait_for_reply=True)

    def _stop_unit(self, request, arguments, who):
        """Kills the job process that runs unit `unit` of stage `stage` of run `run_id` and closes its attempt `stuck` (`reason`, and the `stack` kept in a file). The job is not asked politely: a unit that
        shows no sign of life will not stop between two pieces. The row stays `running` for the launcher, which sees the process end and hands the unit out again while it has a try left."""
        run_id, stage, unit = arguments["run_id"], int(arguments["stage"]), str(arguments["unit"])
        engine = engine_run_of(self.connection, run_id)
        if engine is None:
            return {"ok": True, "jobs": 0}
        key = planrows._stage_key(self.connection, run_id, stage)
        with self._control_connection(run_id) as control:
            row = control.execute("SELECT id, attempt FROM engine_jobs WHERE run_id = ? AND stage = ? AND item_key = ? AND state = 'running'", (engine, key, unit)).fetchone()
            if row is None:
                return {"ok": True, "jobs": 0}
            launcher = self._launcher()
            started = launcher.running.get(row["id"]) if launcher is not None else None
            if started is not None:
                with contextlib.suppress(Exception):
                    started.job.kill()
            ledger.end_running(control, row["id"], "stuck", error_kind="stuck", error_text=arguments.get("reason") or "stuck",
                               stack_file=attemptlog.keep_stack(row["id"], row["attempt"], arguments.get("stack")))
            control.commit()
        self._act(f"stopped unit {unit} of run {run_id}", f"found stuck: {arguments.get('reason') or 'no reason given'}")
        return {"ok": True, "jobs": 1 if started is not None else 0}

    def _end(self, request, control, engine, who, what, wait_for_reply):
        """Asks every job process of the ledger run `engine` (all of them, for `None`) to stop and ends its rows. Returns the reply, or `None` when it is sent once the processes are gone."""
        launcher = self._launcher()
        mine = [started for started in (launcher.running.values() if launcher is not None else ()) if engine is None or started.row["run_id"] == engine]
        for started in mine:
            with contextlib.suppress(Exception):
                started.job.send({"op": "stop"})
        ended = end_rows(control, engine)
        control.commit()
        self._act(f"stopped {what}", f"{len(mine)} job process(es) asked to stop, {ended['rows']} row(s) canceled" + (f"; asked by {who}" if who else ""))
        summary = {"ok": True, "jobs": len(mine), "rows": ended["rows"]}
        if not mine:
            return summary
        self.pending.append(Pending([started.row["id"] for started in mine], self.clock() + self.grace, engine, request if wait_for_reply else None, summary))
        return None

    def _settle(self):
        """Looks at the processes asked to stop: the ones the launcher no longer holds are gone; what is still there after the grace is killed (and given one more grace for the launcher to see it end); the
        reply of a Stop All goes when none is left. The rows are ended once more then, for a unit that put its row back as it stopped."""
        launcher = self._launcher()
        keep = []
        for item in self.pending:
            held = [] if launcher is None else [started for job_id, started in launcher.running.items() if job_id in item.ids]
            now = self.clock()
            if held and not item.killed and now >= item.deadline:
                for started in held:
                    with contextlib.suppress(Exception):
                        started.job.kill()
                item.killed = True
                item.summary["killed"] = len(held)
                item.deadline = now + self.grace
                self._act("killed job processes that did not stop", f"{len(held)} still alive {self.grace:g} s after they were asked to stop")
            if held and now < item.deadline:
                keep.append(item)
                continue
            self._finish(item)
        self.pending = keep

    def _finish(self, item):
        with contextlib.suppress(Exception):
            with self._control_connection() as control:
                end_rows(control, item.engine)
                control.commit()
        self._answer(item.reply, item.summary)

    def clear(self):
        """The supervisor is ending: whoever waits is answered with what is known, and the jobs' pipes close with the launcher."""
        for item in self.pending:
            self._answer(item.reply, item.summary)
        self.pending = []
