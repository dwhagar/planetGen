# planetgen/queue/unitstuck.py

"""
What the follower does with a unit the judge says is stuck (OPS.140; Boss, 5:16 PM EST on 11-09): a unit whose beat stopped, or whose beat is fresh but whose counter has stood
still for far longer than this command has ever stood still (`beatjudge`: dead and quiet, in two checks in a row, `unitwatch.UnitWatch.agreed`), is

1. **dumped**: the follower sends its horse `unitbeat.DUMP_SIGNAL`, which has `faulthandler` write the stack of every thread to the unit's stack file (from the signal handler in C, so
   it answers even when the work loop holds the interpreter), and reads what was written;
2. **stopped**: RQ's stop command (`rq.command.send_stop_job_command`), the way a cancel stops one;
3. **queued again, once**: RQ does not retry a job that was stopped on purpose, so this module does it: it takes the job out of the registry it sits in, clears the unit's stale beat
   and queues the same job (`redisqueue.requeue_job`) once the horse is gone;
4. **run alone (OPS.145)**: the queue is held (`pause.hold_solo`, the one pause flag, naming the reason and the run) so that every other task waits, and the memory gate and the
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

from planetgen.queue import pause, redisqueue, runrecord, unitbeat

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
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


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

    # --- one pass, made at every poll -----------------------------------------------------------------------

    def step(self, agreed, verdicts):
        """
        Follows the stops under way, then acts on the units in `agreed` (`unitwatch.UnitWatch.agreed()`), whose verdicts are in `verdicts`. Does nothing for a run that is not live.
        """
        live = runrecord._s(self.conn.hget(runrecord.run_key(self.run_id), "state")) in runrecord.LIVE_STATES
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
            if runrecord._s(self.conn.hget(runrecord.run_key(self.run_id), "state")) not in runrecord.LIVE_STATES:
                return
            self._act(unit, record, now, verdicts[unit])

    def _lift_when_over(self, done, live):
        """Lifts the solo hold of this stage once its unit is complete or the run has ended."""
        solo = pause.solo_of(self.conn)
        if solo is None or (solo.get("run"), solo.get("stage")) != (self.run_id, self.stage):
            return
        if not live or int(solo["unit"]) in done:
            if pause.lift_solo(self.run_id, self.stage, solo["unit"], self.conn):
                runrecord.log(self.conn, self.run_id, f"unit {solo['unit']} is {'complete' if live else 'no longer running'}: the queue is released", "info", "unit-solo-over",
                              stage=self.stage, unit=solo["unit"])

    def _hold_for_solo_run(self, unit, job):
        """Holds the queue for the restarted `unit` (unless it is turned off or the queue is paused already, by hand or for another unit) and says so in the run log."""
        if not self._solo:
            return
        if pause.hold_solo(self.run_id, self.stage, unit, job.id, self.conn):
            runrecord.log(self.conn, self.run_id, f"the queue is held so that unit {unit} runs alone; it lifts when the unit completes or fails, or when you resume the queue", "warning",
                          "unit-solo-hold", stage=self.stage, unit=unit)
        else:
            runrecord.log(self.conn, self.run_id, f"the queue is paused already, so unit {unit} will not run alone", "warning", "unit-solo-hold", stage=self.stage, unit=unit)

    def records(self):
        """`{unit: {attempts, phase, at, reason, dump}}` of what has been done in this stage."""
        out = {}
        for unit, raw in self.conn.hgetall(self.key).items():
            try:
                out[int(unit)] = json.loads(raw if isinstance(raw, str) else raw.decode())
            except (ValueError, TypeError):
                continue
        return out

    # --- the action -------------------------------------------------------------------------------------------

    def _act(self, unit, record, now, verdict):
        raw = unitbeat.parse(self.conn.hget(runrecord.progress_key(self.run_id, self.stage), unit) or b"{}") or {}
        attempts = int(record.get("attempts", 0)) if record else 0
        reason = f"unit {unit} is {verdict.state}: {verdict.reason}"
        stack = self._dump(raw) if verdict.state != "dead" or alive(raw.get("pid"), raw.get("host")) else "(no stack: the process is gone)"
        self._stop_quietly(unit)
        if attempts >= 1:
            self._fail(unit, reason, record.get("dump", ""), stack, now)
            return
        self.conn.hset(self.key, unit, json.dumps({"attempts": 1, "phase": "stopping", "at": now, "reason": reason, "dump": stack}))
        runrecord.mark_killed(self.conn, self.run_id, reason)
        runrecord.log(self.conn, self.run_id, f"{reason}; its stack is below, it is being stopped and will be queued again once\n{stack}", "warning", "unit-stuck",
                      stage=self.stage, unit=unit)
        self._to_error_log(f"run {self.run_id} stage {self.stage}: {reason}; stopped, to be queued again once\n{stack}", unit)

    def _follow_stop(self, unit, record, now):
        job = self._job(unit)
        status = None if job is None else job.get_status(refresh=True)
        raw = unitbeat.parse(self.conn.hget(runrecord.progress_key(self.run_id, self.stage), unit) or b"{}") or {}
        gone = not alive(raw.get("pid"), raw.get("host"))
        if job is not None and (status in (rq.job.JobStatus.STOPPED, rq.job.JobStatus.FAILED) or (status == rq.job.JobStatus.STARTED and gone)):
            self._requeue(unit, record, job)
        elif job is None or now - float(record.get("at", now)) > STOP_WAIT_SECONDS:
            reason = f"unit {unit} could not be stopped (job {'missing' if job is None else status}) within {STOP_WAIT_SECONDS:.0f} s of being found stuck: {record.get('reason', '')}"
            self._fail(unit, reason, record.get("dump", ""), "", now)

    def _requeue(self, unit, record, job):
        self.conn.hdel(runrecord.progress_key(self.run_id, self.stage), unit)       # the beat of the stopped attempt would read as dead beside the new one
        redisqueue.requeue_job(job, self.conn, at_front=True)           # at the front: it is the one job the held workers are to reach
        self.conn.hset(self.key, unit, json.dumps({**record, "phase": "requeued"}))
        self._hold_for_solo_run(unit, job)
        if self._forget is not None:
            self._forget(unit)
        runrecord.log(self.conn, self.run_id, f"unit {unit} was stopped and is queued again", "warning", "unit-requeued", stage=self.stage, unit=unit)

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

    def _stop_quietly(self, unit):
        from planetgen.queue import units
        try:
            self._stop(units.job_id(self.run_id, self.stage, unit))
        except (rq.exceptions.InvalidJobOperation, rq.exceptions.NoSuchJobError):
            pass        # it ended, or is not running, between the two reads

    def _stop_with_rq(self, job_id):
        command.send_stop_job_command(self.conn, job_id, redisqueue.JsonSerializer)

    def _to_error_log(self, message, unit):
        if self._record_error is not None:
            self._record_error(message, unit=unit)
            return
        from planetgen.util import log
        log.record_error("queue", message, unit=unit)
