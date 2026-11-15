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

from planetgen.queue import ledger, planrows, runrecord
from planetgen.util import log


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
    conn.hset(runrecord.held_key(run_id, stage), unit, f"{row}|{error}"[:1000])
    runrecord.log(conn, run_id, f"unit {unit} of stage {wrapper.stage_name(conn, run_id, stage)} is held: {error}", level="error", kind="unit.held", stage=stage, unit=unit)
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


def retry_held(conn, run_id, who=None, now=None):
    """
    Carries a held run on (OPS.221): takes the machine back for it, puts each held job back to `waiting` and queues it again as it was cut (the arguments and the cost its row holds), and makes the
    run live. The supervisor carries it from there.

    Raises:
        runrecord.Busy: Another run of that kind owns the machine; the run stays held.

    Returns:
        dict | None: `{"requeued": units}`, or `None` when the run is unknown or not held.
    """
    from planetgen.queue import supervisor, units, wrapper
    if runrecord.run_state(conn, run_id) != runrecord.HELD:
        return None
    items = held_units(conn, run_id)
    runrecord.take_held_back(conn, run_id, now)             # Busy leaves the run held and everything as it was
    kind = runrecord._s(conn.hget(runrecord.run_key(run_id), "kind")) or runrecord.GENERATION
    control = planrows.control_connection(conn, run_id)
    try:
        rows = {row: control.execute("SELECT command, arguments FROM engine_jobs WHERE id = ?", (row,)).fetchone() for _stage, _unit, row in items}
        for row in rows:
            ledger.set_job_state(control, row, "waiting")
        engine = conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)
        if engine:
            ledger.set_run_state(control, int(runrecord._s(engine)), "running")
        control.commit()
    finally:
        control.close()
    for stage, unit, row in items:
        planned = json.loads(rows[row]["arguments"])
        conn.hdel(runrecord.held_key(run_id, stage), unit)
        units.requeue_unit(conn, run_id, stage, unit, rows[row]["command"], planned.get("args") or {}, planned.get("cost"), row, queue_name=wrapper.queue_of(kind))
    runrecord.log(conn, run_id, f"{len(items)} held {'unit' if len(items) == 1 else 'units'} queued again" + (f" by {who}" if who else ""), kind="run.retry_held")
    try:
        supervisor.nudge(conn, wrapper.queue_of(kind))
    except Exception as exc:  # noqa: BLE001 -- the next poll of the page or the next supervisor start carries the run on
        log.swallowed("queue/held", "retry_held", exc)
    return {"requeued": len(items)}


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
