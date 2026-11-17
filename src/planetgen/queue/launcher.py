# planetgen/queue/launcher.py

"""
The launcher (OPS.228, slice A; D3 of `/mnt/project-files/reports/job-engine-items-stage-3-6-2026-11-13.md`): the supervisor's turn that hands out the jobs of the ledger. Built and tested, used by
nothing yet: the bulk lane moves to it in the cutover slice, a run at a time.

Each turn (`Launcher.turn`):

1. **Ends.** Every job process that ended (`TemplateSet.poll`) is looked at: a unit writes its own attempt and its row's end (`units.run_unit_on`, as in an RQ worker), so a row that is `done` or `failed`
   is left as the unit wrote it. A row still `running` was left by a process that died without ending it: its attempt is closed `lost` with the signal named and the unit is handed out once more
   (`RETRIES`), else it is `failed`. What a job says on its pipe is read before its end is judged.
2. **Tickets.** Every `ready` job row whose `not_before` has passed is a `scheduler.Ticket` (its claim from the row: `claim_mb`, the lane of its `job_class`), oldest first.
3. **Dispatch.** `dispatch.dispatch` over the tickets, the host and the claims of the jobs this launcher has running says which go now. A job that does not go has its reason written to the row's `waited`.
4. **Starts.** Each admitted job is marked `running` (so no turn hands it out twice), then forked from the template of its lane (`TemplateSet.fork`) running `forkedjob.run`, which takes the attempt in the
   ledger and runs the unit as an RQ worker does. (Taking the attempt moves to this turn when the units' own bookkeeping is deleted, OPS.228's last slice.)

The launcher holds no decision of its own: the order and the reasons are `dispatch`'s, the fork is `jobhost`'s, the rows are `ledger`'s. Its state is the jobs it has started and not yet seen end; the ledger
holds the rest, so a launcher that starts over reads `running` rows it does not hold (recovery is the supervisor's, a later slice).
"""

import json
import os
import signal
import time

from planetgen.queue import dispatch, jobhost, ledger, memorygate, planrows, scheduler
from planetgen.queue.memorygate import Claim, Holder

RETRIES = 1
"""int: How many times a unit whose process died without ending its attempt is handed out again (what RQ's `Retry(max=1)` did)."""

RUNNER = "planetgen.queue.forkedjob.run"
"""str: What a handed-out job runs in its process."""

READY_WINDOW = 500
"""int: How many `ready` rows a turn reads (the oldest; `dispatch` judges the top 64 of them)."""

LANES = {"bulk": memorygate.LANE_BULK, "maintenance": memorygate.LANE_MAINTENANCE, "interactive": memorygate.LANE_INTERACTIVE, "long": memorygate.LANE_LONG}
"""dict: The lane of each `job_class` a row can have; any other class is the `other` lane."""

RECOVER_SECONDS = 30.0
"""float: How often the rows left `running` by a supervisor that ended are looked for."""

HOLDER_NAME = "job-{}"


def claim_of_row(row):
    """The claim of a job row: the kind is `unit:<command>`, the memory its `claim_mb` (the minimum when the plan gave none), one core, the lane of its class."""
    mb = float(row["claim_mb"]) if row["claim_mb"] else float(memorygate.MIN_CLAIM_MB)
    return Claim(f"unit:{row['command']}", mb, 1.0, LANES.get(row["job_class"], memorygate.LANE_OTHER))


def _arguments(row):
    try:
        loaded = json.loads(row["arguments"]) if row["arguments"] else {}
    except ValueError:
        return {}
    return loaded if isinstance(loaded, dict) else {"value": loaded}


class Started:
    """A job this launcher handed out and has not seen end."""

    def __init__(self, row, job, key, claim):
        self.row, self.job, self.key, self.claim = row, job, key, claim
        self.said = []

    def drain(self):
        self.said += self.job.messages()

    def skipped(self):
        """Whether the unit said it did nothing: its run is no longer live, its unit already done or held (`units.run_unit_on` returns `skipped` and writes no attempt)."""
        return any(message.get("event") == "unit" and (message.get("result") or {}).get("skipped") for message in self.said)

    def error(self):
        """The last thing the job said about its failure (`error`), or `None`."""
        for message in reversed(self.said):
            if message.get("event") in ("failed", "error"):
                return message.get("error")
        return None


class Launcher:
    """One supervisor's hand-outs. `control` is the control connection; `forker` is a `jobhost.TemplateSet` (`fork(kind, runner, job)`, `poll()`); `galaxy` is `run_id -> connection` for the
    galaxy the run writes to (or `None`)."""

    def __init__(self, control, forker, galaxy=None, clock=time.time, read_host=memorygate.read_host, used=None, runner=RUNNER):
        self.control, self.forker, self.galaxy = control, forker, galaxy
        self.clock, self.read_host, self.runner = clock, read_host, runner
        self.used = used if used is not None else (lambda: {})
        self.running = {}           # job id -> Started
        self.recovered_at = None

    # -- ends ------------------------------------------------------------------------------------------------------------------

    def _outcome(self, started, exit_):
        """`(outcome, error kind, error text)` of an end the unit did not write itself."""
        error = started.error()
        if exit_.signal is not None:
            return "lost", "killed", error or f"the job process ended on {signal.Signals(exit_.signal).name}"
        if exit_.code is None:
            return "lost", "lost", error or "the job process ended and nothing said how"
        return "failed", "failed", error or f"the job process exited with code {exit_.code}"

    def ended(self):
        """Settles every process that ended (see the module). Returns the ids of the jobs that ended."""
        done = []
        for _key, exit_ in self.forker.poll():
            found = next((job_id for job_id, started in self.running.items() if started.job.id == exit_.id and started.key == _key), None)
            if found is None:
                continue
            started = self.running.pop(found)
            started.drain()
            row = self.control.execute("SELECT state, attempt FROM engine_jobs WHERE id = ?", (found,)).fetchone()
            opening = started.row["command"] == planrows.OPENING
            skipped = started.skipped()
            if row is not None and row["state"] == "running" and exit_.code == 0 and (skipped or opening):
                ledger.end_running(self.control, found, "done")             # an opening has the launcher's attempt; a unit's attempt is its own (none when it skipped)
                finished = opening and not skipped
                ledger.set_job_state(self.control, found, "done" if finished else "canceled", waited=None if finished else "skipped: the run is no longer live or the unit was already done")
                self.control.commit()
            elif row is not None and row["state"] == "running":
                outcome, kind, text = self._outcome(started, exit_)
                ledger.end_running(self.control, found, outcome, error_kind=kind, error_text=text)
                ledger.set_job_state(self.control, found, "ready" if row["attempt"] <= RETRIES else "failed", waited=None if row["attempt"] <= RETRIES else text[:ledger.WAITED_CHARS])
                self.control.commit()
            started.job.close()
            done.append(found)
        return done

    def _gone(self, pid):
        """Whether the process `pid` is not there (a job whose attempt names no process is gone too)."""
        if not pid:
            return True
        try:
            os.kill(int(pid), 0)
        except ProcessLookupError:
            return True
        except PermissionError:
            return False
        return False

    def recover(self):
        """The rows `running` that this launcher did not start (a supervisor that ended took its launcher's memory with it): one whose attempt's process is gone, or that has no attempt (the supervisor
        ended between marking it and forking it), is settled like a job that died (its attempt closed `lost`, the job handed out once more, else `failed`). A row whose process is alive is the unit's own to end.
        Returns the ids settled."""
        now = self.clock()
        if self.recovered_at is not None and now - self.recovered_at < RECOVER_SECONDS:
            return []
        self.recovered_at = now
        settled = []
        for row in self.control.execute("SELECT j.id, j.attempt, a.pid AS pid FROM engine_jobs j LEFT JOIN engine_attempts a ON a.job_id = j.id AND a.attempt = j.attempt WHERE j.state = 'running'").fetchall():
            if row["id"] in self.running or not self._gone(row["pid"]):
                continue
            text = "the supervisor that handed the job out ended and its process is gone"
            ledger.end_running(self.control, row["id"], "lost", error_kind="lost", error_text=text)
            ledger.set_job_state(self.control, row["id"], "ready" if row["attempt"] <= RETRIES else "failed", waited=None if row["attempt"] <= RETRIES else text)
            settled.append(row["id"])
        self.control.commit()
        return settled

    # -- tickets and dispatch --------------------------------------------------------------------------------------------------

    def tickets(self):
        """`(tickets, rows)` of the ready jobs: the tickets for `dispatch`, the rows by ticket name."""
        rows = self.control.execute(
            "SELECT id, run_id, stage, command, item_key, arguments, claim_mb, job_class, waited, attempt, UNIX_TIMESTAMP(updated_at) AS since FROM engine_jobs"
            " WHERE state = 'ready' AND (not_before IS NULL OR not_before <= NOW(6)) ORDER BY id LIMIT ?", (READY_WINDOW,)).fetchall()
        named = {HOLDER_NAME.format(row["id"]): row for row in rows}
        return [scheduler.Ticket(name, claim_of_row(row), float(row["since"]), 0, "launcher") for name, row in named.items()], named

    def holders(self):
        return [Holder(HOLDER_NAME.format(job_id), started.claim, started.job.pid, "launcher", started.claim.mb, started.claim.mb) for job_id, started in self.running.items()]

    def _note(self, rows, turn):
        for name, reason in turn.reasons.items():
            row = rows.get(name)
            if row is None or row["waited"] == reason[:ledger.WAITED_CHARS]:
                continue
            self.control.execute("UPDATE engine_jobs SET waited = ? WHERE id = ? AND state = 'ready'", (reason[:ledger.WAITED_CHARS], row["id"]))
        self.control.commit()

    def _job_of(self, row):
        """What `forkedjob.run` takes for the unit `row`: the Redis run id and the stage's index (from the stage's `open` row, `planrows.record_cut`), the unit's number, command, arguments and cost, and its
        row. A unit is given a retry when it has not had one."""
        if row["command"] == planrows.OPENING:         # the cut of a stage (`planrows.queue_opening`)
            named = _arguments(row)
            return {"opening": True, "run_id": named["run_id"], "stage": named["index"], "kind": row["command"], "row": row["id"], "retries_left": max(0, RETRIES - row["attempt"])}
        opened = self.control.execute("SELECT arguments FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = 'open' AND item_key = ''", (row["run_id"], row["stage"])).fetchone()
        plan = json.loads(opened["arguments"]) if opened and opened["arguments"] else {}
        if "run_id" not in plan or "index" not in plan:
            raise LookupError(f"the stage {row['stage']!r} of ledger run {row['run_id']} has no open row that names its run and index")
        unit = _arguments(row)
        return {"run_id": plan["run_id"], "stage": plan["index"], "unit": int(row["item_key"]), "kind": row["command"], "args": unit.get("args", {}), "cost": unit.get("cost"), "row": row["id"],
                "retries_left": max(0, RETRIES - row["attempt"])}

    def start(self, row):
        """Hands the job out: marked `running`, then one fork. Returns the `Started`, or `None` for a row another turn took."""
        if self.control.execute("UPDATE engine_jobs SET state = 'running', waited = NULL WHERE id = ? AND state = 'ready'", (row["id"],)).rowcount == 0:
            self.control.rollback()
            return None
        self.control.commit()
        try:
            key, forked = self.forker.fork(row["command"], self.runner, self._job_of(row))
        except jobhost.NotReady as exc:                 # its template is still being built: the row is handed out at a later turn
            ledger.set_job_state(self.control, row["id"], "ready", waited=str(exc)[:ledger.WAITED_CHARS])
            self.control.commit()
            return None
        except Exception as exc:
            ledger.set_job_state(self.control, row["id"], "failed", waited=f"{type(exc).__name__}: {exc}"[:ledger.WAITED_CHARS])
            self.control.commit()
            raise
        if row["command"] == planrows.OPENING:          # an opening has no unit to begin its attempt: the launcher does, with the process that runs it
            ledger.begin_attempt(self.control, row["id"], pid=forked.pid)
        started = Started(row, forked, key, claim_of_row(row))
        self.running[row["id"]] = started
        return started

    def turn(self):
        """One turn: ends, tickets, dispatch, starts. Returns the ids of the jobs started."""
        self.ended()
        self.recover()
        waiting, rows = self.tickets()
        if not waiting:
            return []
        now = self.clock()
        decided = dispatch.dispatch(waiting, self.read_host(), self.holders(), now, self.used())
        self._note(rows, decided)
        started = [self.start(rows[name]) for name in decided.admitted]
        return [item.row["id"] for item in started if item is not None]

    def stop(self):
        """Lets go of every job (their pipes close, so each finishes and exits at its next safe point), stops the templates and closes the control connection."""
        for started in self.running.values():
            started.job.close()
        self.running = {}
        for held in (getattr(self.forker, "stop", None), getattr(self.control, "close", None)):
            if held is not None:
                held()
