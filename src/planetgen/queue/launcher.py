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
import types

from planetgen.queue import claims, dispatch, jobhost, ledger, memorygate, planrows, runrecord, scheduler, workers
from planetgen.queue.memorygate import Claim, Holder
from planetgen.util import log
from planetgen.util.settings import get_settings

RETRIES = 1
"""int: How many times a unit whose process died without ending its attempt is handed out again (what RQ's `Retry(max=1)` did)."""

RUNNER = "planetgen.queue.forkedjob.run"
"""str: What a handed-out job runs in its process."""

READY_WINDOW = 500
"""int: How many `ready` rows a turn reads (the oldest; `dispatch` judges the top 64 of them)."""

LANES = {"bulk": memorygate.LANE_BULK, "maintenance": memorygate.LANE_MAINTENANCE, "interactive": memorygate.LANE_INTERACTIVE, "long": memorygate.LANE_LONG}
"""dict: The lane of each `job_class` a row can have; any other class is the `other` lane."""

TEMPLATES = {"maintenance": "app"}
"""dict: The template (`jobhost.WARM_SETS`) a job of a class is forked from when it is not the one its command names (`jobhost.template_key`): the maintenance lane's units build the site's app, which its template holds built (OPS.229)."""

YIELDED_TEXT = "stopped to give memory back; it starts again when the host has room"
"""str: What the row of a job the memory controller asked to stop says (`engine_jobs.waited`)."""

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

    def __init__(self, row, job, key, claim, began=0.0):
        self.row, self.job, self.key, self.claim, self.began = row, job, key, claim, began
        self.said = []
        self.asked = False          # the controller asked it to stop at its next piece (`Launcher.yield_youngest`)

    def drain(self):
        self.said += self.job.messages()

    def skipped(self):
        """Whether the unit said it did nothing: its run is no longer live, its unit already done or held (`units.run_unit_on` returns `skipped` and writes no attempt)."""
        return any(message.get("event") == "unit" and (message.get("result") or {}).get("skipped") for message in self.said)

    def stopped(self):
        """Whether the unit said it stopped between two pieces (`units._work_of_unit` returns `stopped`): asked to by the supervisor, or its run was no longer live."""
        return any(message.get("event") == "unit" and (message.get("result") or {}).get("stopped") for message in self.said)

    def error(self):
        """The last thing the job said about its failure (`error`), or `None`."""
        for message in reversed(self.said):
            if message.get("event") in ("failed", "error"):
                return message.get("error")
        return None


class Launcher:
    """One supervisor's hand-outs. `control` is the control connection; `forker` is a `jobhost.TemplateSet` (`fork(kind, runner, job)`, `poll()`); `galaxy` is `run_id -> connection` for the
    galaxy the run writes to (or `None`); `also_held` returns the claims of the jobs another duty of the supervisor holds (OPS.230: the site's light jobs)."""

    def __init__(self, control, forker, galaxy=None, clock=time.time, read_host=memorygate.read_host, used=None, runner=RUNNER, redis=None, also_held=None, hold=None):
        self.control, self.forker, self.galaxy, self.redis = control, forker, galaxy, redis
        self.hold = hold or (lambda: None)      # the memory controller's `holding` (OPS.231): why nothing is admitted now, or `None`
        self.also_held = also_held or (lambda: [])      # the claims the site's light jobs hold (`lightjobs.LightJobs.holders`): the bulk lane is judged against them
        self.clock, self.read_host, self.runner = clock, read_host, runner
        self.used = used if used is not None else (lambda: {})
        self.running = {}           # job id -> Started
        self.recovered_at = None
        self._stages = {}           # (ledger run, stage) -> its `open` row's arguments (`_stage_of`)

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
            elif row is not None and row["state"] in ("running", "waiting") and exit_.code == 0 and started.asked and started.stopped():
                ledger.end_running(self.control, found, "lost", error_kind=ledger.STOPPED, error_text="stopped between two pieces to give memory back")       # the unit wrote its own end as `lost`, kind `stopped`
                ledger.set_job_state(self.control, found, "ready", waited=YIELDED_TEXT)         # it goes again when there is room, and the stop does not count as one of its attempts
                self.control.commit()
            elif row is not None and row["state"] == "running":
                outcome, kind, text = self._outcome(started, exit_)
                ledger.end_running(self.control, found, outcome, error_kind=kind, error_text=text)
                again = row["attempt"] - ledger.stopped_attempts(self.control, found) <= RETRIES
                ledger.set_job_state(self.control, found, "ready" if again else "failed", waited=None if again else text[:ledger.WAITED_CHARS])
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
        rows = self._within_caps(rows)
        named = {HOLDER_NAME.format(row["id"]): row for row in rows}
        return [scheduler.Ticket(name, self.claim(row), float(row["since"]), 0, "launcher") for name, row in named.items()], named

    def _stage_of(self, row):
        """What the `open` row of the stage of `row` says (`{"run_id", "index", "cap"}`; `{}` while the stage is not cut). Kept once read, for the cut never changes."""
        key = (row["run_id"], row["stage"])
        if key not in self._stages:
            found = self.control.execute("SELECT arguments FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = 'open' AND item_key = '' AND state = 'done'", key).fetchone()
            if found is None:
                return {}
            self._stages[key] = json.loads(found["arguments"] or "{}") or {}
        return self._stages[key]

    def _cap_of(self, row):
        """The most units of the stage of `row` that run at once: the run's own limit (`--workers`, the `cap` of the stage's `open` row) and, with the run's Redis in hand, what a share of the host's memory holds
        of the stage's units (`workers.stage_worker_cap`, PERF.220), the lesser of the two; `None` when neither limits it."""
        opened = self._stage_of(row)
        asked = opened.get("cap")
        if self.redis is None or row["command"] == planrows.OPENING or "run_id" not in opened:
            return asked
        try:
            held = workers.stage_worker_cap(row["command"], self.redis, opened["run_id"])
        except Exception as exc:  # noqa: BLE001 -- a figure that cannot be read leaves the memory gate to decide
            log.swallowed("queue/launcher", "_cap_of", exc)
            held = None
        return min(asked, held) if asked and held else (asked or held)

    def claim(self, row):
        """The claim of the job row: what the gate has learned of its kind (`memorygate.claim_mb`, the job profile's figure once this machine has history) raised from the typed figure, the cores and named
        resources of its command, the cut's whole lane for an opening; the typed figure of the row alone where the launcher has no Redis to read the learned ones from (the tests)."""
        if self.redis is None:
            return claim_of_row(row)
        opened = self._stage_of(row)
        token = opened.get("run_id")
        lane = LANES.get(row["job_class"], memorygate.LANE_OTHER)
        try:
            if row["command"] == planrows.OPENING:
                named = _arguments(row)
                command = claims.cut_command(self.redis, named["run_id"], named["index"])
                kind = f"cut:{command or 'unknown'}"
                return Claim(kind, memorygate.claim_mb(self.redis, kind, claims.CUT_MB, token=named["run_id"]), 1.0, lane, True)
            kind = f"unit:{row['command']}"
            return Claim(kind, memorygate.claim_mb(self.redis, kind, claims.unit_table_mb(row["command"]), token=token), claims.UNIT_CPU, lane, False, claims.unit_tokens(row["command"]))
        except claims.NoClaim:
            return claim_of_row(row)

    def _upkeep_limit(self):
        """The most units of the maintenance class (the cache warm-up and the checks, `planrows.JOB_CLASSES`) that run at once, or `None` for no limit of their own (OPS.229): `queue.maintenance_workers` while a generation run
        is live, so upkeep never takes more than its guaranteed share from the run (PERF.150); none when no generation run is (the memory gate and the stage's own cap still hold)."""
        if self.redis is None:
            return None
        try:
            live = runrecord.active_run(self.redis, runrecord.GENERATION)
        except Exception as exc:  # noqa: BLE001 -- a run that cannot be read leaves the upkeep to the memory gate
            log.swallowed("queue/launcher", "_upkeep_limit", exc)
            return None
        return get_settings().queue.maintenance_workers if live else None

    def _within_caps(self, rows):
        """`rows` less the units the stage's limit holds back (OPS.228): of each stage, no more units run at once than `_cap_of` says; and, of the maintenance class, no more at once than `_upkeep_limit` (OPS.229)."""
        taken, upkeep = {}, 0
        for started in self.running.values():
            if started.row["command"] != planrows.OPENING:
                key = (started.row["run_id"], started.row["stage"])
                taken[key] = taken.get(key, 0) + 1
            upkeep += started.row["job_class"] == "maintenance"
        limit = self._upkeep_limit() if any(row["job_class"] == "maintenance" for row in rows) else None
        kept = []
        for row in rows:
            if limit and row["job_class"] == "maintenance" and row["command"] != planrows.OPENING:
                if upkeep >= limit:
                    continue
                upkeep += 1
            if row["command"] == planrows.OPENING:
                kept.append(row)
                continue
            key, cap = (row["run_id"], row["stage"]), self._cap_of(row)
            if cap and taken.get(key, 0) >= cap:
                if limit and row["job_class"] == "maintenance":
                    upkeep -= 1
                continue
            taken[key] = taken.get(key, 0) + 1
            kept.append(row)
        return kept

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
            return {"opening": True, "run_id": named["run_id"], "stage": named["index"], "kind": row["command"], "row": row["id"], "retries_left": max(0, RETRIES - row["attempt"] + ledger.stopped_attempts(self.control, row["id"]))}
        opened = self.control.execute("SELECT arguments FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = 'open' AND item_key = ''", (row["run_id"], row["stage"])).fetchone()
        plan = json.loads(opened["arguments"]) if opened and opened["arguments"] else {}
        if "run_id" not in plan or "index" not in plan:
            raise LookupError(f"the stage {row['stage']!r} of ledger run {row['run_id']} has no open row that names its run and index")
        unit = _arguments(row)
        return {"run_id": plan["run_id"], "stage": plan["index"], "unit": int(row["item_key"]), "kind": row["command"], "args": unit.get("args", {}), "cost": unit.get("cost"), "row": row["id"],
                "retries_left": max(0, RETRIES - row["attempt"] + ledger.stopped_attempts(self.control, row["id"]))}

    def start(self, row):
        """Hands the job out: marked `running`, then one fork. Returns the `Started`, or `None` for a row another turn took."""
        if self.control.execute("UPDATE engine_jobs SET state = 'running', waited = NULL WHERE id = ? AND state = 'ready'", (row["id"],)).rowcount == 0:
            self.control.rollback()
            return None
        self.control.commit()
        try:
            key, forked = self.forker.fork(row["command"], self.runner, self._job_of(row), **({"template": TEMPLATES[row["job_class"]]} if row["job_class"] in TEMPLATES else {}))
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
        started = Started(row, forked, key, self.claim(row), self.clock())
        self.running[row["id"]] = started
        return started

    def turn(self):
        """One turn: ends, tickets, dispatch, starts. Returns the ids of the jobs started."""
        self.ended()
        self.recover()
        waiting, rows = self.tickets()
        if not waiting:
            return []
        held = self.hold()
        if held:                    # OPS.231: the host is short of memory; the bulk and maintenance lanes wait, each row saying why
            self._note(rows, types.SimpleNamespace(reasons={name: held for name in rows}))
            return []
        now = self.clock()
        decided = dispatch.dispatch(waiting, self.read_host(), [*self.holders(), *self.also_held()], now, self.used())
        self._note(rows, decided)
        started = [self.start(rows[name]) for name in decided.admitted]
        return [item.row["id"] for item in started if item is not None]

    def yield_youngest(self):
        """
        The memory controller's ask (OPS.231): asks the youngest unit that runs on the bulk lane (the one that began last, so the one that has done the least) to stop at its next piece, once; of the units
        whose command resumes after the pieces it finished (`commands.spec(kind).resumes`) when there are any, since a stop costs the others their work. Returns `(row id, command)` of the job asked, or `None`
        when none can be asked (an opening is never asked: its cut is whole or nothing).
        """
        from planetgen.queue import commands
        candidates = [started for started in self.running.values()
                      if not started.asked and started.row["command"] != planrows.OPENING and started.row["job_class"] == "bulk"]
        resumable = [started for started in candidates if started.row["command"] in commands.REGISTRY and commands.spec(started.row["command"]).resumes]
        pool = resumable or candidates
        if not pool:
            return None
        chosen = max(pool, key=lambda started: (started.began, started.row["id"]))
        try:
            chosen.job.send({"op": "stop"})
        except (OSError, ValueError) as exc:
            log.swallowed("queue/launcher", "yield_youngest", exc)
            return None
        chosen.asked = True
        return chosen.row["id"], chosen.row["command"]

    def stop(self):
        """Lets go of every job (their pipes close, so each finishes and exits at its next safe point), stops the templates and closes the control connection."""
        for started in self.running.values():
            started.job.close()
        self.running = {}
        for held in (getattr(self.forker, "stop", None), getattr(self.control, "close", None)):
            if held is not None:
                held()
