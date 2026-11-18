# planetgen/queue/held.py

"""
Held items and held runs (OPS.221, L7 of `docs/design/operations/job-engine.md`): a unit that fails for good no longer ends the run.

* **The unit is held.** A unit that fails for good (the error reached `units.unit_failed`: nothing in the engine retries it, a kill with no retry left comes here too) is written to the run's held hash
  of its stage (`runrecord.held_key`: unit to ledger row), its job row is `failed` (`engine_jobs`; its attempt rows carry the error, the traceback and the stack file, OPS.223), and the run's log
  says so. Its queued retry finds it held and does nothing (`is_held`).
* **The rest goes on.** Every other unit of the stage runs as it did. Nothing else in a run depends on a stage that did not end (the next opening is queued by the unit that completes the stage), so
  the units not held are all there is to do.
* **The run is held when nothing more can run** (`settle`, called when a unit is held and after every unit that ends while one is): every unit of the stage is done or held. The run takes the state
  `held` (`runrecord.hold_run`): it frees the machine (it does not count as the live run, so the line moves), no worker takes a unit of it, and the starter that follows it ends with the reason. The
  ledger's run row says `held` with why.
* **Retry held** (`retry_held`, from the Queue page, `planetgen queue retry-held` and `POST /api/admin/work/<id>/control`) takes the machine back for the run, puts each held job back to `waiting`
  and queues it again with the arguments its row holds (`units.requeue_unit`); the supervisor carries the run on from there (a run with no starter is its business already, OPS.210). Cancel ends a held run as it
  ends a live one.

A unit that is only a *piece* of an opening, a conclusion or an end of a stage is not held: those still fail the run (`units.open_failed`) until OPS.218's jobs carry the same state.
"""

import json
import os
import traceback

from planetgen.queue import ledger, planrows, runrecord
from planetgen.util import log, secrets_mask

WHERE_FRAMES = 3
"""int: How many of the last frames of a failure's traceback the held error text names (OPS.275)."""


def where(tb, frames=WHERE_FRAMES):
    """
    Where a failure happened, for the held error text (OPS.275): ` (at file:line in function, ...)` for the last `frames` frames of the traceback `tb`, the innermost first, files named from
    the `planetgen` package (or by their own name outside it). Empty when there is no traceback. A one-off failure under load can then be traced from the held text and the log without the stack file.
    """
    if tb is None:
        return ""
    named = []
    for frame in reversed(traceback.extract_tb(tb)[-frames:]):
        path = frame.filename.replace(os.sep, "/")
        named.append(f"{path[path.rfind('planetgen/'):] if 'planetgen/' in path else os.path.basename(path)}:{frame.lineno} in {frame.name}")
    return f" (at {', '.join(named)})" if named else ""


def is_held(conn, run_id, stage, unit):
    """Whether the unit is held (the retry RQ queued after its failure finds it so and does nothing)."""
    return bool(conn.hexists(runrecord.held_key(run_id, stage), unit))


def held_units(conn, run_id):
    """`[(stage, unit, row)]` of the units the run holds, in stage and unit order."""
    found = []
    for stage in range(len(runrecord.declared_stages(conn, run_id))):
        for unit, value in conn.hgetall(runrecord.held_key(run_id, stage)).items():
            found.append((stage, int(runrecord._s(unit)), int(runrecord._s(value).split("|", 1)[0])))
    return sorted(found)


def _errors(conn, run_id, stage):
    """`{unit: error}` of the stage's held units (the held hash's value is `row|error`)."""
    return {int(runrecord._s(unit)): runrecord._s(value).split("|", 1)[1] if "|" in runrecord._s(value) else ""
            for unit, value in conn.hgetall(runrecord.held_key(run_id, stage)).items()}


def hold_unit(conn, run_id, stage, unit, row, error, killed=False):
    """
    Holds the unit that failed for good: it joins its stage's held hash, its job row is `failed`, the log says what failed, and the run is held when this was the last thing that could run. `error` is the
    exception in words (`ValueError: ...`); `killed` says the process was killed, as the run record's flag does for a run that ended so.
    """
    from planetgen.queue import wrapper         # the wrapper queues units and imports this file's callers
    error = secrets_mask.mask_text(error)
    conn.hset(runrecord.held_key(run_id, stage), unit, f"{row}|{error}"[:1000])
    runrecord.log_said(conn, run_id, f"unit {unit} of stage {wrapper.stage_name(conn, run_id, stage)} is held: {error}", level="error", kind="unit.held", stage=stage, unit=unit)
    _mark_job(conn, run_id, row, "failed")
    if killed:
        runrecord.mark_killed(conn, run_id, error)
    settle(conn, run_id)


def settle(conn, run_id):
    """
    Holds the run if it has a held unit and every unit of that stage is done or held, so nothing more can run. A run with nothing held, or that is not live, is left alone.

    Returns:
        bool: Whether this call held the run.
    """
    if runrecord.run_state(conn, run_id) != "running":
        return False
    for stage in range(len(runrecord.declared_stages(conn, run_id))):
        held = [int(runrecord._s(unit)) for unit in conn.hkeys(runrecord.held_key(run_id, stage))]
        if not held:
            continue
        done = conn.scard(runrecord.done_key(run_id, stage))
        total = int(float(runrecord._s(conn.hget(runrecord.stage_key(run_id, stage), "units")) or 0))
        if done + sum(1 for unit in held if not conn.sismember(runrecord.done_key(run_id, stage), unit)) < total:
            return False                # units of the stage are still queued or running: they go on
        from planetgen.queue import wrapper
        name = wrapper.stage_name(conn, run_id, stage)
        errors = _errors(conn, run_id, stage)
        first = min(errors)
        reason = (f"unit {first} of stage {name} failed: {errors[first]}" if len(held) == 1 else
                  f"{len(held)} units of stage {name} failed for good (the first, unit {first}: {errors[first]})")
        reason += "; everything that did not depend on it is done. The run is held: Retry held continues it after a fix, cancel ends it"
        if runrecord.hold_run(conn, run_id, reason):
            mark_run(conn, run_id, "held", reason)
            return True
        return False
    return False


DRAIN_REASON = "stopped for update"
"""str: Why a run is held when an update drains it (OPS.222, part A)."""


def drain_run(conn, run_id, reason=DRAIN_REASON, now=None):
    """
    Drains the live run for an update (OPS.222, part A of L8): the run is held with `reason`, which frees the machine and takes it off every worker; each running unit stops at its next piece (the
    engine finds the run no longer live, `engine.UnitStopped`), is not failed, and goes back to `waiting` in the ledger with what it committed kept. The retry of an unfinished unit is `retry_held`, which
    the end of the update runs.

    Returns:
        bool: Whether this call held the run (`False` for a run that is not live).
    """
    if not runrecord.hold_run(conn, run_id, reason, now):
        return False
    conn.hset(runrecord.run_key(run_id), mapping={HELD_FOR: "update", HELD_EPOCH: _epoch()})
    mark_run(conn, run_id, "held", reason)
    return True


HELD_FOR, HELD_EPOCH = "held_for", "held_epoch"
"""str: The run-record fields of a drained run: what it is held for (`update`) and the generator epoch it was drained under."""


def is_drained(conn, run_id):
    """Whether the run is held for an update (`drain_run`): the update's stop of the job that follows it must not cancel it."""
    return runrecord.run_state(conn, run_id) == runrecord.HELD and runrecord._s(conn.hget(runrecord.run_key(run_id), HELD_FOR)) == "update"


def _epoch():
    """The generator epoch the code makes galaxies by now (`util.generator_epoch`), or `""` where it cannot be read."""
    try:
        from planetgen.util import generator_epoch
        return str(generator_epoch.generator_epoch())
    except Exception as exc:  # noqa: BLE001 -- a drain must not fail for it; a resume then finds no epoch and leaves the run held
        log.swallowed("queue/held", "_epoch", exc)
        return ""


def resume_after_update(conn, now=None):
    """
    The end of an update's drain (OPS.222, part A): every run held "for update" is carried on (`retry_held`) when the generator epoch is the one it was drained under, so the galaxy it finishes is the
    one an undisturbed run makes; under another epoch (or one that cannot be read) it stays held and says why. Another run owning the machine leaves it held too.

    Returns:
        list[str]: One line for each held run, for the update's output (`updateguard.end`).
    """
    lines = []
    here = _epoch()
    for run_id in runrecord.recent_runs(conn, limit=100_000):
        if runrecord.run_state(conn, run_id) != runrecord.HELD or runrecord._s(conn.hget(runrecord.run_key(run_id), HELD_FOR)) != "update":
            continue
        title = runrecord._s(conn.hget(runrecord.run_key(run_id), "title")) or run_id
        then = runrecord._s(conn.hget(runrecord.run_key(run_id), HELD_EPOCH))
        if not here or then != here:
            lines.append(f"{title} ({run_id}) stays held: the generator changed with the update (epoch {then or '?'} then, {here or '?'} now), so it is run again, not carried on.")
            continue
        try:
            carried = retry_held(conn, run_id, who="the update", now=now)
        except runrecord.Busy as exc:
            lines.append(f"{title} ({run_id}) stays held: {exc}")
            continue
        conn.hdel(runrecord.run_key(run_id), HELD_FOR, HELD_EPOCH)
        lines.append(f"{title} ({run_id}) goes on: {carried['requeued']} unfinished {'unit' if carried['requeued'] == 1 else 'units'} queued again." if carried else f"{title} ({run_id}) was not held any more.")
    return lines


def retry_held(conn, run_id, who=None, now=None):
    """
    Carries a held run on (OPS.221, and the end of an update's drain, OPS.222): takes the machine back for it, puts each unit that is not done (the held ones, and the ones a drain stopped) back to
    `waiting` and queues it again as it was cut (the arguments and the cost its row holds) unless the queue still has it, queues the opening of a stage whose opening the hold skipped, and makes the run live.
    The supervisor carries it from there.

    Raises:
        runrecord.Busy: Another run of that kind owns the machine; the run stays held.

    Returns:
        dict | None: `{"requeued": units}`, or `None` when the run is unknown or not held.
    """
    from planetgen.queue import supervisor, units
    if runrecord.run_state(conn, run_id) != runrecord.HELD:
        return None
    runrecord.take_held_back(conn, run_id, now)             # Busy leaves the run held and everything as it was
    control = planrows.control_connection(conn, run_id)
    try:
        engine = conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)
        items = _unfinished(conn, control, run_id, int(runrecord._s(engine))) if engine else []
        requeued = 0
        for _stage, _unit, row, _command, _planned in items:
            # the launcher hands a failed (held) or waiting unit out again as it was cut (the arguments and the cost its row holds); one that is ready or running already goes on as it is
            requeued += control.execute("UPDATE engine_jobs SET state = 'ready', waited = NULL, updated_at = NOW(6) WHERE id = ? AND state IN ('failed', 'waiting')", (row,)).rowcount
        if engine:
            control.execute("UPDATE engine_jobs SET state = 'waiting', waited = NULL, updated_at = NOW(6) WHERE run_id = ? AND command = ? AND item_key = '' AND state = 'failed'",
                            (int(runrecord._s(engine)), planrows.SPAWN))        # the spawner of a stage that could not write its next window (OPS.240) tries again
            ledger.set_run_state(control, int(runrecord._s(engine)), "running")
        control.commit()
    finally:
        control.close()
    for stage, unit, _row, _command, _planned in items:
        conn.hdel(runrecord.held_key(run_id, stage), unit)
    reopened = units.reopen_pending(conn, run_id)
    runrecord.log(conn, run_id, f"{requeued} unfinished {'unit' if requeued == 1 else 'units'} queued again" + (f" and the opening of stage {reopened + 1}" if reopened is not None else "")
                  + (f" by {who}" if who else ""), kind="run.retry_held")
    try:
        supervisor.ensure_started(conn)
    except Exception as exc:  # noqa: BLE001 -- the next poll of the page or the next supervisor start carries the run on
        log.swallowed("queue/held", "retry_held", exc)
    return {"requeued": requeued}


def _unfinished(conn, control, run_id, engine_run):
    """`[(stage index, unit, row, command, arguments)]` of the run's units the ledger does not say are done or canceled, in plan order."""
    index_of = {planrows._stage_key(conn, run_id, index): index for index in range(len(runrecord.declared_stages(conn, run_id)))}
    found = []
    for row in control.execute("SELECT id, stage, command, item_key, arguments FROM engine_jobs WHERE run_id = ? AND item_key <> '' AND state NOT IN ('done', 'canceled') ORDER BY id", (engine_run,)).fetchall():
        if row["stage"] in index_of:
            found.append((index_of[row["stage"]], int(row["item_key"]), row["id"], row["command"], json.loads(row["arguments"] or "{}")))
    return found


def _mark_job(conn, run_id, row, state):
    """Sets the job row's state in the ledger; best effort (the held hash is what the run goes by, the row is its history)."""
    try:
        control = planrows.control_connection(conn, run_id)
    except Exception as exc:  # noqa: BLE001
        log.swallowed("queue/held", "_mark_job", exc)
        return
    try:
        ledger.set_job_state(control, row, state)
        control.commit()
    except Exception as exc:  # noqa: BLE001
        log.swallowed("queue/held", "_mark_job", exc)
    finally:
        control.close()


def mark_run(conn, run_id, state, reason):
    """Sets the ledger's run row to `state` with its reason (the cancel of a held run uses it too); best effort, as `_mark_job`."""
    engine = conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)
    if not engine:
        return
    try:
        control = planrows.control_connection(conn, run_id)
    except Exception as exc:  # noqa: BLE001
        log.swallowed("queue/held", "mark_run", exc)
        return
    try:
        ledger.set_run_state(control, int(runrecord._s(engine)), state, reason[:255])
        control.commit()
    except Exception as exc:  # noqa: BLE001
        log.swallowed("queue/held", "mark_run", exc)
    finally:
        control.close()
