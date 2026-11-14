# planetgen/queue/unitwatch.py

"""
The follower's side of the unit heartbeat (OPS.138): at every poll the process that started a run reads the beat of each unit running in the stage it follows
(`unitbeat.read`, written by the horse, OPS.136), has `beatjudge.judge` say what each is doing, logs **each change of state once**, and keeps the verdicts where a page
can read them (`runrecord.unitstate_key`, served by `GET /api/admin/work/<id>/units`).

* **A change is logged once.** A unit that goes quiet gets one line in the run log (and, for dead, quiet and suspect, one line in `error.log`); it is not repeated at every poll,
  and a unit that works again gets one line saying so. A unit that starts working is not announced. A reading that fails (Redis down for a moment) is skipped, never a verdict.
* **Acting needs two checks in a row** (`beatjudge.agreed`): `UnitWatch.agreed()` lists the units whose last two checks both say dead or quiet, which is what the stuck action
  (OPS.140) takes. A unit whose counter moves is never in it.
* **Stuck units are acted on** by `unitstuck.Stuck` (OPS.140): dumped, stopped, queued again once, then the run fails.
* **Only this process judges.** The process that only watches a run (`wrapper.follow`) does not make a watch, so two processes never log the same change twice.

    watch = unitwatch.UnitWatch(conn, run_id, stage, kind)
    watch.check()                  # at every poll: {unit: Verdict}
"""

import json
import time

from planetgen.queue import beatjudge, runrecord, taskprogress, unitbeat, unitstuck
from planetgen.util import log


class UnitWatch:
    """
    Judges the running units of one stage at every `check` and says what changed.

    Args:
        conn (redis.Redis): The run's Redis.
        run_id (str): The run.
        stage (int): The stage followed.
        kind (str): The command of the stage, for its history (`history_of`).
        history_of (callable, optional): `history_of(kind) -> beatjudge.History`; a plain `History()` where left out (OPS.135 supplies the recorded gaps).
        starved (callable, optional): `starved() -> bool`, whether the host is starved (the quiet limit doubles).
        record_error (callable, optional): `record_error(message, **fields)`, how an error-level change reaches `error.log` (`log.record_error("queue", ...)` where left out).
        rate_of (callable, optional): `rate_of(kind) -> items a second or None`, the rate recorded for the command (PERF.147: `unitrates.rate_of`); with it each unit's bar has an ETA before its own rate is known.
        stuck (unitstuck.Stuck, optional): What is done with the units that two checks in a row call dead or quiet (OPS.140); a `unitstuck.Stuck` of this stage where left out, `False` for none.
    """

    def __init__(self, conn, run_id, stage, kind, history_of=None, starved=None, record_error=None, stuck=None, rate_of=None):
        self.conn = conn
        self.run_id = run_id
        self.stage = stage
        self.kind = kind
        self.history_of = history_of
        self.starved = starved
        self.record_error = record_error
        self.rate_of = rate_of
        self._beats = {}            # the beat each unit had at the check before, for `judge`'s `before`
        self._previous = {}         # the verdict of the check before
        self._current = {}          # the verdict now
        self._announced = {}        # (state, severity) each unit was last logged in; a unit not in it has only worked
        self.stuck = unitstuck.Stuck(conn, run_id, stage, record_error=record_error, forget=self.forget) if stuck is None else (stuck or None)

    def history(self):
        return beatjudge.History() if self.history_of is None else self.history_of(self.kind)

    def check(self):
        """
        Reads the beats, judges each unit and logs the changes; returns `{unit: Verdict}` of the units that have a record now. A unit that finished (its record was removed)
        is forgotten.
        """
        try:
            found = unitbeat.read(self.conn, self.run_id, self.stage)
            now = unitbeat.server_time(self.conn)
        except Exception as exc:   # noqa: BLE001 -- a Redis that cannot answer for one poll says nothing about the units
            log.swallowed("queue/unitwatch", "UnitWatch.check", exc)
            return dict(self._current)
        history = self.history()
        starved = bool(self.starved()) if self.starved is not None else False
        verdicts = {unit: beatjudge.judge(beat, history, now, before=self._beats.get(unit), starved=starved) for unit, beat in found.items()}
        for unit in set(self._announced) - set(found):
            self._announced.pop(unit, None)
        self._previous, self._current, self._beats = self._current, verdicts, found
        for unit, verdict in sorted(verdicts.items()):
            self._announce(unit, verdict, found[unit])
        self._publish(verdicts, found, now, self._stored_rate())
        if self.stuck is not None:
            self.stuck.step(self.agreed(), self._current)
        return dict(verdicts)

    def forget(self, unit):
        """A unit that was queued again starts afresh: nothing of its last attempt is remembered."""
        for held in (self._beats, self._previous, self._current, self._announced):
            held.pop(unit, None)

    def clear(self):
        """Removes the published verdicts (the stage is over)."""
        self.conn.delete(runrecord.unitstate_key(self.run_id, self.stage))

    def agreed(self):
        """The units that dead or quiet in two checks in a row (`beatjudge.agreed`): what the stuck action may take."""
        return sorted(unit for unit, verdict in self._current.items() if beatjudge.agreed(self._previous.get(unit), verdict))

    def _announce(self, unit, verdict, beat):
        shown = (verdict.state, verdict.severity)
        before = self._announced.get(unit)
        if shown == (beatjudge.WORKING, "") and before is None:
            return
        if shown == before:
            return
        self._announced[unit] = shown
        if verdict.state == beatjudge.WORKING:
            runrecord.log(self.conn, self.run_id, f"unit {unit} works again: {verdict.reason or 'the count moves'}", "info", "unit-state", stage=self.stage, unit=unit)
            return
        where = f" in {beat.phase}" if beat.phase else ""
        message = f"unit {unit} is {verdict.state}{where}: {verdict.reason}"
        if verdict.state == beatjudge.BUSY and verdict.severity and self.stuck is not None:
            message += f"; nothing is stopped, its stack is below\n{self.stuck.stack_of(unit)}"
        error = verdict.state in (beatjudge.DEAD, beatjudge.QUIET) or verdict.severity == beatjudge.SUSPECT
        runrecord.log(self.conn, self.run_id, message, "error" if error else "warning" if verdict.severity or verdict.state == beatjudge.LOOPING else "info", "unit-state",
                      stage=self.stage, unit=unit)
        if error:
            self._to_error_log(f"run {self.run_id} stage {self.stage}: {message}", unit)

    def _to_error_log(self, message, unit):
        if self.record_error is not None:
            self.record_error(message, unit=unit)
            return
        from planetgen.util import log
        log.record_error("queue", message, unit=unit)

    def _stored_rate(self):
        """The items a second recorded for this stage's command (`rate_of`), `None` where there is none."""
        if self.rate_of is None:
            return None
        try:
            return self.rate_of(self.kind)
        except Exception as exc:       # noqa: BLE001 -- a bar without a recorded rate still counts
            log.swallowed("queue/unitwatch", "UnitWatch._stored_rate", exc)
            return None

    def _publish(self, verdicts, found, now, stored_rate=None):
        """The verdicts where the page reads them: one hash field to a unit, the units that ended removed."""
        key = runrecord.unitstate_key(self.run_id, self.stage)
        pipe = self.conn.pipeline()
        stale = {int(unit) for unit in self.conn.hkeys(key)} - set(verdicts)
        if stale:
            pipe.hdel(key, *sorted(stale))
        for unit, verdict in verdicts.items():
            pipe.hset(key, unit, json.dumps(view(verdict, found[unit], now, unit in self.agreed(), stored_rate)))
        pipe.execute()


def view(verdict, beat, now, agreed=False, stored_rate=None):
    """A verdict with the beat it was made from, as the page and `/units` show it (plain numbers and text). `progress` is the unit's own bar (PERF.147: `taskprogress.as_dict`)."""
    return {"state": verdict.state, "severity": verdict.severity, "reason": verdict.reason, "actionable": verdict.actionable, "agreed": bool(agreed),
            "counter": verdict.counter, "rate": verdict.rate, "projected_finish_seconds": verdict.projected_finish_seconds, "still_seconds": verdict.still_seconds,
            "over_limit": verdict.over_limit, "phase": beat.phase, "waiting_on": beat.waiting_on, "done": beat.done, "estimate": beat.estimate,
            "restarts": beat.restarts, "resident_mb": beat.resident_mb, "sql": beat.sql, "sql_seconds": None if beat.sql_since is None or not beat.sql else max(now - beat.sql_since, 0.0), "beat_age_seconds": max(now - beat.beat, 0.0), "running_seconds": max(now - beat.started, 0.0),
            "progress": taskprogress.as_dict(taskprogress.progress(beat.done, beat.estimate, beat.started, now, stored_rate)), "step": step_view(beat, now)}


def step_view(beat, now):
    """The unit's second-level count (OPS.192) as the pages show it, `None` where the unit has none: `name`, `done`, `total` and `progress` (its own bar and time left from its own rate)."""
    if not beat.step:
        return None
    started = beat.step_started if beat.step_started is not None else beat.started
    return {"name": beat.step, "done": beat.step_done, "total": beat.step_total,
            "progress": taskprogress.as_dict(taskprogress.progress(beat.step_done, beat.step_total, started, now))}


def read(conn, run_id, stage):
    """`{unit: view}` as the follower last published it (empty where it has not)."""
    out = {}
    for unit, raw in conn.hgetall(runrecord.unitstate_key(run_id, stage)).items():
        try:
            out[int(unit)] = json.loads(raw if isinstance(raw, str) else raw.decode())
        except (ValueError, TypeError):
            continue
    return out
