# planetgen/queue/launcher.py

"""
The launcher (OPS.228, slice A; D3 of `/mnt/project-files/reports/job-engine-items-stage-3-6-2026-11-13.md`): the supervisor's turn that hands out the jobs of the ledger. Built and tested, used by
nothing yet: the bulk lane moves to it in the cutover slice, a run at a time.

Each turn (`Launcher.turn`):

1. **Ends.** Every job process that ended (`TemplateSet.poll`) gets its attempt row closed with how it ended and its job row moved on: `done` for exit 0, else `failed` with the error the job said last
   (a signal or an exit with nothing said is `lost`). What a job says on its pipe is read before its end is judged.
2. **Tickets.** Every `ready` job row whose `not_before` has passed is a `scheduler.Ticket` (its claim from the row: `claim_mb`, the lane of its `job_class`), oldest first.
3. **Dispatch.** `dispatch.dispatch` over the tickets, the host and the claims of the jobs this launcher has running says which go now. A job that does not go has its reason written to the row's `waited`.
4. **Starts.** Each admitted job is handed out: `ledger.begin_attempt` (the one write the turn waits for, which also writes the attempt number to the galaxy's `engine_items` so the fence has a number), then
   one fork from the template of its lane (`TemplateSet.fork`) running `forkedjob.run`.

The launcher holds no decision of its own: the order and the reasons are `dispatch`'s, the fork is `jobhost`'s, the rows are `ledger`'s. Its state is the jobs it has started and not yet seen end; the ledger
holds the rest, so a launcher that starts over reads `running` rows it does not hold (recovery is the supervisor's, a later slice).
"""

import json
import signal
import time

from planetgen.queue import dispatch, ledger, memorygate, scheduler
from planetgen.queue.memorygate import Claim, Holder

RUNNER = "planetgen.queue.forkedjob.run"
"""str: What a handed-out job runs in its process."""

READY_WINDOW = 500
"""int: How many `ready` rows a turn reads (the oldest; `dispatch` judges the top 64 of them)."""

LANES = {"bulk": memorygate.LANE_BULK, "maintenance": memorygate.LANE_MAINTENANCE, "interactive": memorygate.LANE_INTERACTIVE, "long": memorygate.LANE_LONG}
"""dict: The lane of each `job_class` a row can have; any other class is the `other` lane."""

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

    def __init__(self, row, attempt, key, job, claim):
        self.row, self.attempt, self.key, self.job, self.claim = row, attempt, key, job, claim
        self.said = []

    def drain(self):
        self.said += self.job.messages()

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

    # -- ends ------------------------------------------------------------------------------------------------------------------

    def _outcome(self, started, exit_):
        if exit_.code == 0:
            return "done", None, None
        error = started.error()
        if exit_.signal is not None:
            return "lost", "killed", error or f"the job process ended on {signal.Signals(exit_.signal).name}"
        if exit_.code is None:
            return "lost", "lost", error or "the job process ended and nothing said how"
        return "failed", "failed", error or f"the job process exited with code {exit_.code}"

    def ended(self):
        """Closes the attempt and moves the job of every process that ended. Returns the ids of the jobs that ended."""
        done = []
        for _key, exit_ in self.forker.poll():
            found = next((job_id for job_id, started in self.running.items() if started.job.id == exit_.id and started.key == _key), None)
            if found is None:
                continue
            started = self.running.pop(found)
            started.drain()
            outcome, kind, text = self._outcome(started, exit_)
            ledger.end_attempt(self.control, found, started.attempt, outcome, error_kind=kind, error_text=text)
            ledger.set_job_state(self.control, found, "done" if outcome == "done" else "failed")
            self.control.commit()
            started.job.close()
            done.append(found)
        return done

    # -- tickets and dispatch --------------------------------------------------------------------------------------------------

    def tickets(self):
        """`(tickets, rows)` of the ready jobs: the tickets for `dispatch`, the rows by ticket name."""
        rows = self.control.execute(
            "SELECT id, run_id, stage, command, item_key, arguments, claim_mb, job_class, waited, UNIX_TIMESTAMP(updated_at) AS since FROM engine_jobs"
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

    def start(self, row):
        """Hands the job out: the attempt is written (and its number written to the galaxy), then one fork. Returns the `Started`."""
        galaxy = self.galaxy(row["run_id"]) if self.galaxy is not None else None
        attempt = ledger.begin_attempt(self.control, row["id"], galaxy_conn=galaxy)
        job = dict(_arguments(row), row=row["id"], attempt=attempt, run=row["run_id"], stage=row["stage"], command=row["command"], item_key=row["item_key"])
        try:
            key, forked = self.forker.fork(row["command"], self.runner, job)
        except Exception as exc:
            ledger.end_attempt(self.control, row["id"], attempt, "failed", error_kind="fork", error_text=f"{type(exc).__name__}: {exc}")
            ledger.set_job_state(self.control, row["id"], "failed")
            self.control.commit()
            raise
        self.control.execute("UPDATE engine_attempts SET pid = ? WHERE job_id = ? AND attempt = ?", (forked.pid, row["id"], attempt))
        self.control.commit()
        started = Started(row, attempt, key, forked, claim_of_row(row))
        self.running[row["id"]] = started
        return started

    def turn(self):
        """One turn: ends, tickets, dispatch, starts. Returns the ids of the jobs started."""
        self.ended()
        waiting, rows = self.tickets()
        if not waiting:
            return []
        now = self.clock()
        decided = dispatch.dispatch(waiting, self.read_host(), self.holders(), now, self.used())
        self._note(rows, decided)
        return [self.start(rows[name]).row["id"] for name in decided.admitted]

    def stop(self):
        """Lets go of every job (their pipes close, so each finishes and exits at its next safe point)."""
        for started in self.running.values():
            started.job.close()
        self.running = {}
