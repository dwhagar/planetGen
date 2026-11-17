# planetgen/queue/unitstuck.py

"""
What the follower does with a unit the judge says is stuck (OPS.140; Boss, 5:16 PM EST on 11-09): a unit whose beat stopped, or whose beat is fresh but whose counter has stood
still for far longer than this command has ever stood still (`beatjudge`: dead and quiet, in two checks in a row, `unitwatch.UnitWatch.agreed`), is

1. **dumped**: the follower sends its horse `unitbeat.DUMP_SIGNAL`, which has `faulthandler` write the stack of every thread to the unit's stack file (from the signal handler in C, so
   it answers even when the work loop holds the interpreter), and reads what was written;
2. **stopped**: for a run in the ledger, a request to the supervisor (`controls.ask`, `stop_unit`), which kills the job process it forked and closes the attempt `stuck` in the ledger with the
   reason and the stack (OPS.234); for a run with no rows in the ledger (until OPS.271), RQ's stop command (`rq.command.send_stop_job_command`), the way a cancel stops one;
3. **queued again, once**: for a ledgered run the launcher hands the unit out once more when it settles the process's end (the row goes back to `ready`), and this module only clears the unit's
   stale beat when it sees that; RQ does not retry a job that was stopped on purpose, so for the others this module does it: it takes the job out of the registry it sits in, clears the beat
   and queues the same job (`redisqueue.requeue_job`) once the horse is gone;
4. **run alone (OPS.145; RQ runs only, the launcher's hand-outs do not read a solo hold)**: the queue is held (`pause.hold_solo`, the one pause flag, naming the reason and the run) so that every other task waits, and the memory gate and the
   scheduler give the restarted unit the host (`memorygate.Claim.solo`); the hold lifts when the unit completes or the run ends, here or, if this follower is gone, when the unit's
   job has ended (`pause.state`), and an admin can lift it from the Queue page;
5. **failed, the second time**: a unit that is stuck again on its second attempt fails the run, with both stacks in the message.

A unit whose counter moves is never here, at any size over its estimate: `beatjudge` only calls dead and quiet actionable. Every stop or failure sets the run's killed flag
(`runrecord.mark_killed`, OPS.144; `finish_run(killed=True)` for the failure) with the reason. What was done to each unit is kept in `runrecord.stuck_key` (one JSON field to a unit), so a follower that is restarted
neither forgets a stop it began nor stops the same unit twice. The unit's time limit is only the backstop (`units.unit_timeout`, 6 hours).
"""

import json
import os
import socket
import time

import rq
from rq import command

from planetgen.queue import beatjudge, pause, redisqueue, runrecord, unitbeat
from planetgen.util import process
from planetgen.util.format import format_duration_seconds

DUMP_WAIT_SECONDS = 3.0
"""float: How long the follower waits for a horse to write its stack after the signal."""

STOP_WAIT_SECONDS = 120.0
"""float: How long a stopped unit may take to leave its horse before the stop is taken to have failed (and the run with it)."""

MESSAGE_STACK_CHARS = 6000
"""int: How much of a stack goes into the run's error message (the whole of it is in the run log)."""


def alive(pid, host):
    """Whether process `pid` on this machine exists; `True` for one on another host (it cannot be told, and a stop is then waited for)."""
    if host != socket.gethostname() or not isinstance(pid, int):
        return True
    return process.pid_alive(pid)


def request_dump(raw, signal_process=os.kill, wait=DUMP_WAIT_SECONDS, sleep=time.sleep):
    """
    Asks the horse that wrote the beat record `raw` for its stack and returns it as text (`""` with a note when it cannot be had: another host, no signal, a process that is gone).
    """
    path, pid, host = raw.get("stack_file"), raw.get("pid"), raw.get("host")
    if not path or not isinstance(pid, int) or host != socket.gethostname() or unitbeat.DUMP_SIGNAL is None:
        return "(no stack: the horse is not on this machine)"
    try:
        before = os.path.getsize(path) if os.path.exists(path) else 0
        signal_process(pid, unitbeat.DUMP_SIGNAL)
    except ProcessLookupError:
        return "(no stack: the process is gone)"
    except OSError as exc:
        return f"(no stack: {exc})"
    deadline = wait
    while deadline > 0:
        sleep(0.1)
        deadline -= 0.1
        try:
            if os.path.getsize(path) > before:
                sleep(0.1)          # the dump is written in one go; this lets the last lines land
                with open(path, "r", encoding="utf-8", errors="replace") as handle:
                    handle.seek(before)
                    return handle.read().strip()
        except OSError:
            break
    return "(no stack: the horse did not answer)"


class Stuck:
    """
    The action on the stuck units of one stage.

    Args:
        conn (redis.Redis): The run's Redis.
        run_id (str): The run.
        stage (int): The stage followed.
        dump (callable, optional): `dump(raw_record) -> str`; `request_dump` where left out.
        stop (callable, optional): `stop(job_id)`, how a started job is stopped; RQ's stop command where left out.
        record_error (callable, optional): `record_error(message, **fields)`, how a failure reaches `error.log`.
        forget (callable, optional): `forget(unit)`, told when a unit is queued again, so the watch starts afresh on it.
        solo (bool, optional): Hold the queue so the restarted unit runs alone (OPS.145; on; the tests of the stop alone turn it off).
    """

    def __init__(self, conn, run_id, stage, dump=None, stop=None, record_error=None, forget=None, solo=True):
        self.conn, self.run_id, self.stage = conn, run_id, stage
        self.key = runrecord.stuck_key(run_id, stage)
        self._dump = dump or request_dump
        self._stop = stop or self._stop_with_rq
        self._record_error = record_error
        self._forget = forget
        self._solo = solo
        self._ledger = None

    # --- one pass, made at every poll -----------------------------------------------------------------------

    def step(self, agreed, verdicts):
        """
        Follows the stops under way, then acts on the units in `agreed` (`unitwatch.UnitWatch.agreed()`), whose verdicts are in `verdicts`. Does nothing for a run that is not live.
        """
        live = self._live()
        done = {int(unit) for unit in self.conn.smembers(runrecord.done_key(self.run_id, self.stage))}
        self._lift_when_over(done, live)
        if not live:
            return
        now = unitbeat.server_time(self.conn)
        records = self.records()
        for unit, record in sorted(records.items()):
            if record.get("phase") == "stopping" and unit not in done:
                self._follow_stop(unit, record, now)
        records = self.records()
        for unit in agreed:
            record = records.get(unit)
            if unit in done or (record is not None and record.get("phase") == "stopping") or unit not in verdicts:
                continue        # a unit queued again by this very pass was forgotten by the watch: its verdict belongs to the attempt that was stopped
            if not self._live():
                return
            if verdicts[unit].state == beatjudge.DEAD and self._waiting_to_run(unit):
                continue        # OPS.268: its horse is gone because the attempt failed and the retry is queued (RQ's, or the ledger's `ready`): it is waiting for a worker, not dead, and queuing it again would run it twice
            self._act(unit, record, now, verdicts[unit])

    def _lift_when_over(self, done, live):
        """Lifts the solo hold of this stage once its unit is complete or the run has ended."""
        solo = pause.solo_of(self.conn)
        if solo is None or (solo.get("run"), solo.get("stage")) != (self.run_id, self.stage):
            return
        if not live or int(solo["unit"]) in done:
            if pause.lift_solo(self.run_id, self.stage, solo["unit"], self.conn):
                runrecord.log_said(self.conn, self.run_id, f"unit {solo['unit']} is {'complete' if live else 'no longer running'}: the queue is released", "info", "unit-solo-over",
                              stage=self.stage, unit=solo["unit"])

    def _ledgered(self):
        """Whether the run has rows in the ledger (its units are forked jobs the supervisor holds), as against one whose units are RQ jobs."""
        if self._ledger is None:
            from planetgen.queue import controls
            self._ledger = controls.engine_run_of(self.conn, self.run_id) is not None
        return self._ledger

    def _ledger_row(self, unit):
        """`{"id", "state", "attempt"}` of the unit's newest job row in the ledger, or `None`."""
        from planetgen.queue import planrows
        engine = self.conn.hget(runrecord.run_key(self.run_id), planrows.ENGINE_RUN_FIELD)
        if not engine:
            return None
        control = planrows.control_connection(self.conn, self.run_id)
        try:
            return control.execute("SELECT id, state, attempt FROM engine_jobs WHERE run_id = ? AND stage = ? AND item_key = ? ORDER BY id DESC LIMIT 1",
                                   (int(runrecord._s(engine)), planrows._stage_key(self.conn, self.run_id, self.stage), str(unit))).fetchone()
        finally:
            control.close()

    def _hold_for_solo_run(self, unit, job):
        """Holds the queue for the restarted `unit` (unless it is turned off or the queue is paused already, by hand or for another unit) and says so in the run log."""
        if not self._solo:
            return
        if pause.hold_solo(self.run_id, self.stage, unit, job.id, self.conn):
            runrecord.log_said(self.conn, self.run_id, f"the queue is held so that unit {unit} runs alone; it lifts when the unit completes or fails, or when you resume the queue", "warning",
                          "unit-solo-hold", stage=self.stage, unit=unit)
        else:
            runrecord.log_said(self.conn, self.run_id, f"the queue is paused already, so unit {unit} will not run alone", "warning", "unit-solo-hold", stage=self.stage, unit=unit)

    def _live(self):
        return runrecord.run_state(self.conn, self.run_id) in runrecord.LIVE_STATES

    def _beat_record(self, unit):
        """The beat record of `unit` as stored (`{}` where it has none)."""
        return unitbeat.parse(self.conn.hget(runrecord.progress_key(self.run_id, self.stage), unit) or b"{}") or {}

    def records(self):
        """`{unit: {attempts, phase, at, reason, dump}}` of what has been done in this stage."""
        return runrecord.unit_json_hash(self.conn, self.key)

    def stack_of(self, unit):
        """The stack of the unit's horse, asked for now (a busy unit that has run past its learned time: OPS.194); nothing is stopped."""
        raw = unitbeat.parse(self.conn.hget(runrecord.progress_key(self.run_id, self.stage), unit) or b"{}") or {}
        return self._dump(raw)

    # --- the action -------------------------------------------------------------------------------------------

    def _act(self, unit, record, now, verdict):
        raw = self._beat_record(unit)
        attempts = int(record.get("attempts", 0)) if record else 0
        reason = f"unit {unit} is {verdict.state}: {verdict.reason}"
        stack = self._dump(raw) if verdict.state != "dead" or alive(raw.get("pid"), raw.get("host")) else "(no stack: the process is gone)"
        if self._ledgered():
            self._stop_in_ledger(unit, reason, stack)
        else:
            self._stop_quietly(unit)
            self._end_attempt(unit, reason, stack)
        if attempts >= 1:
            self._fail(unit, reason, record.get("dump", ""), stack, now)
            return
        found = self._ledger_row(unit) if self._ledgered() else None
        self.conn.hset(self.key, unit, json.dumps({"attempts": 1, "phase": "stopping", "at": now, "reason": reason, "dump": stack,
                                                   **({"attempt": found["attempt"]} if found is not None else {})}))
        runrecord.mark_killed(self.conn, self.run_id, reason)
        runrecord.log_said(self.conn, self.run_id, f"{reason}; its stack is below, it is being stopped and will be queued again once\n{stack}", "warning", "unit-stuck",
                      stage=self.stage, unit=unit)
        self._to_error_log(f"run {self.run_id} stage {self.stage}: {reason}; stopped, to be queued again once\n{stack}", unit)

    def _stop_in_ledger(self, unit, reason, stack):
        """Asks the supervisor to kill the unit's job process and close its attempt `stuck` (`controls`, `stop_unit`); the launcher then hands the unit out again. Not waited for: `_follow_stop` looks at the row."""
        from planetgen.queue import controls
        controls.ask(self.conn, "stop_unit", "the stuck-unit follower", wait=0, run_id=self.run_id, stage=self.stage, unit=unit, reason=reason, stack=stack)

    def _end_attempt(self, unit, reason, stack):
        """OPS.223: the stopped attempt's row in the ledger ends `stuck`, with why and the stack the horse wrote (kept in a file the row names). A horse that is stopped cannot write its own end."""
        from planetgen.queue import attempts as ledger_attempts, units
        try:
            job = self._job(unit)
            if job is not None:
                ledger_attempts.close_lost(self.conn, self.run_id, units.row_of(job), "stuck", error_kind="stuck", error_text=reason, stack_text=stack)
        except Exception as exc:  # noqa: BLE001 -- the ledger's history never stops the follower
            from planetgen.util import log
            log.swallowed("queue/unitstuck", "_end_attempt", exc)

    def _follow_stop(self, unit, record, now):
        if self._ledgered():
            self._follow_ledger_stop(unit, record, now)
            return
        job = self._job(unit)
        status = None if job is None else job.get_status(refresh=True)
        raw = self._beat_record(unit)
        gone = not alive(raw.get("pid"), raw.get("host"))
        if job is not None and (status in (rq.job.JobStatus.STOPPED, rq.job.JobStatus.FAILED) or (status == rq.job.JobStatus.STARTED and gone)):
            self._requeue(unit, record, job)
        elif job is None or now - float(record.get("at", now)) > STOP_WAIT_SECONDS:
            reason = f"unit {unit} could not be stopped (job {'missing' if job is None else status}) within {format_duration_seconds(STOP_WAIT_SECONDS)} of being found stuck: {record.get('reason', '')}"
            self._fail(unit, reason, record.get("dump", ""), "", now)

    def _follow_ledger_stop(self, unit, record, now):
        """A ledgered unit is stopped when its row is no longer the `running` one the supervisor was asked to kill: back to `ready` (or already handed out again: a later attempt), or over."""
        row = self._ledger_row(unit)
        if row is not None and (row["state"] != "running" or row["attempt"] > int(record.get("attempt", row["attempt"]))):
            self.conn.hdel(runrecord.progress_key(self.run_id, self.stage), unit)       # the beat of the stopped attempt would read as dead beside the new one
            self.conn.hset(self.key, unit, json.dumps({**record, "phase": "requeued"}))
            if self._forget is not None:
                self._forget(unit)
            runrecord.log_said(self.conn, self.run_id, f"unit {unit} was stopped and is handed out again", "warning", "unit-requeued", stage=self.stage, unit=unit)
        elif row is None or now - float(record.get("at", now)) > STOP_WAIT_SECONDS:
            reason = f"unit {unit} could not be stopped (job {'missing' if row is None else row['state']}) within {format_duration_seconds(STOP_WAIT_SECONDS)} of being found stuck: {record.get('reason', '')}"
            self._fail(unit, reason, record.get("dump", ""), "", now)

    def _requeue(self, unit, record, job):
        self.conn.hdel(runrecord.progress_key(self.run_id, self.stage), unit)       # the beat of the stopped attempt would read as dead beside the new one
        redisqueue.requeue_job(job, self.conn, at_front=True)           # at the front: it is the one job the held workers are to reach
        self.conn.hset(self.key, unit, json.dumps({**record, "phase": "requeued"}))
        self._hold_for_solo_run(unit, job)
        if self._forget is not None:
            self._forget(unit)
        runrecord.log_said(self.conn, self.run_id, f"unit {unit} was stopped and is queued again", "warning", "unit-requeued", stage=self.stage, unit=unit)

    def _fail(self, unit, reason, first, second, now):
        from planetgen.queue import wrapper         # the wrapper queues units, so it imports this file
        name = wrapper.stage_name(self.conn, self.run_id, self.stage)
        pause.lift_solo(self.run_id, self.stage, unit, self.conn)
        runrecord.mark_killed(self.conn, self.run_id, reason)
        shown = (f"\nfirst time:\n{first[:MESSAGE_STACK_CHARS]}" if first else "") + (f"\nthis time:\n{second[:MESSAGE_STACK_CHARS]}" if second else "")
        runrecord.finish_run(self.conn, self.run_id, "failed", error=f"unit {unit} of stage {name} was stuck twice and the run stops: {reason}{shown}", now=now,
                             killed=True)
        self._to_error_log(f"run {self.run_id}: unit {unit} of stage {name} was stuck twice; {reason}{shown}", unit)

    # --- RQ ---------------------------------------------------------------------------------------------------

    def _job(self, unit):
        from planetgen.queue import units
        return redisqueue.fetch_job(units.job_id(self.run_id, self.stage, unit), self.conn, redisqueue.JsonSerializer)

    def _waiting_to_run(self, unit):
        """Whether the unit's RQ job is waiting for a worker (queued, scheduled or deferred), as a retry of a failed attempt is, so the old attempt's silent beat is not a unit that died."""
        if self._ledgered():
            row = self._ledger_row(unit)
            return row is not None and row["state"] in ("ready", "waiting")
        job = self._job(unit)
        if job is None:
            return False
        return job.get_status(refresh=True) in (rq.job.JobStatus.QUEUED, rq.job.JobStatus.SCHEDULED, rq.job.JobStatus.DEFERRED)

    def _stop_quietly(self, unit):
        from planetgen.queue import units
        try:
            self._stop(units.job_id(self.run_id, self.stage, unit))
        except (rq.exceptions.InvalidJobOperation, rq.exceptions.NoSuchJobError):
            pass        # it ended, or is not running, between the two reads

    def _stop_with_rq(self, job_id):
        command.send_stop_job_command(self.conn, job_id, redisqueue.JsonSerializer)

    def _to_error_log(self, message, unit):
        to_error_log(self._record_error, message, unit)


def to_error_log(record_error, message, unit):
    """`record_error(message, unit=unit)` where a follower was given one, else `log.record_error("queue", ...)` (`error.log`)."""
    if record_error is not None:
        record_error(message, unit=unit)
        return
    from planetgen.util import log
    log.record_error("queue", message, unit=unit)
