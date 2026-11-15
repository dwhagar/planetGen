# planetgen/queue/planrows.py

"""
The plan is rows (OPS.216, L2 of `docs/design/operations/job-engine.md`): the cut of a stage writes its units into the job ledger (`engine_jobs`, OPS.215) and
marks itself done, in one transaction, and everything that wants to know what a stage is made of reads the rows.

* **The cut.** `record_cut` inserts one job row per unit (stage = the stage's key, command = the unit's kind, item key = the unit's index, arguments = the unit's arguments
  and predicted cost) and the row of the cut itself (command `open`, state `done`, `arguments` = the plan summary the stage's command gave; `makes` is a TEXT column and a plan with the starts of a big galaxy is longer), and commits once. A cut that
  stopped half way left nothing; a cut that finished is found by its `open` row and is not made again, so the retry of a killed opening reuses the same rows and the same
  ids, and a stage is never cut into two sets of units.
* **The ids.** `enqueue_units` hands RQ each unit's row id with the unit (`units.run_unit`'s `row`), so what runs is what the plan says.
* **The reading.** `read_plan` is the stage's plan summary (what the Redis string `<prefix>:run:<id>:s<k>:plan` held until now, and no longer does); `unit_rows` the
  stage's units. A plan never changes once it is cut, so a plan that was read is kept for the life of the process.
* **The run.** A run's rows hang on an `engine_runs` row, made by the first cut and named in the run record's `engine_run` field.

The ledger lives in the control database, found from the database the run's steps name (`wrapper.expected_seconds` reads the same one).
"""

import json

from planetgen.queue import ledger, runrecord

OPEN = "open"
"""str: The command of the row that is the cut of a stage."""

CONCLUDE, FINISH, RELEASE = "conclude", "finish", "release"
END_COMMANDS = (CONCLUDE, FINISH, RELEASE)
"""tuple[str]: The jobs that end a stage (OPS.218), each with an edge to the one before: `conclude` needs every unit, `finish` needs `conclude`, and `release` (the next opening queued) needs `finish`.
They are rows of the stage with no item key, as the opening is."""

ENGINE_RUN_FIELD = "engine_run"
"""str: The field of the run record that holds the run's `engine_runs` id."""

_PLANS = {}
"""dict: `(key prefix, run id, stage index)` to the plan read (it never changes once the stage is cut)."""


def control_connection(conn, run_id):
    """An open connection to the control database: on the server of the database the run's first step with one names, else on the configured one."""
    from planetgen.db import store
    for item in runrecord.read_steps(conn, run_id):
        named = (item.get("params") or {}).get("database")
        if named:
            return store.get_control_connection(store.control_mysql_config(store.MySQLConfig.from_args(named)))
    return store.get_control_connection(store.control_mysql_config())


def _galaxy_name(conn, run_id):
    for item in runrecord.read_steps(conn, run_id):
        named = (item.get("params") or {}).get("database")
        if named:
            return named.get("database")
    return None


def engine_run(conn, run_id, control):
    """
    The `engine_runs` id of the run: the one its record names, else a new row, committed and then named in the record (`HSETNX`: when another process got there first its
    id stands and this row is removed).
    """
    key = runrecord.run_key(run_id)
    found = conn.hget(key, ENGINE_RUN_FIELD)
    if found:
        return int(runrecord._s(found))
    made = ledger.add_run(control, runrecord._s(conn.hget(key, "kind")) or "generate", _galaxy_name(conn, run_id))
    control.commit()
    if conn.hsetnx(key, ENGINE_RUN_FIELD, made):
        return made
    control.execute("DELETE FROM engine_runs WHERE id = ?", (made,))
    control.commit()
    return int(runrecord._s(conn.hget(key, ENGINE_RUN_FIELD)))


def _stage_key(conn, run_id, index):
    return runrecord._s(runrecord.declared_stages(conn, run_id)[index]["key"])


def record_cut(conn, run_id, index, kind, units, cost, plan_detail):
    """
    Writes the cut of stage `index` as rows, in one transaction, and returns the unit rows' ids in unit order.

    Args:
        kind (str): The command of the stage (the kind of each unit).
        units (list[tuple]): `(args, cost)` per unit, as the stage's expansion gave them.
        cost (float): The units' total predicted cost.
        plan_detail: What the expansion said of the stage as a whole (`commands.Expansion.plan`).

    A stage already cut (its `open` row is done) is not cut again: the rows it has are returned, whatever this call was given.
    """
    stage = _stage_key(conn, run_id, index)
    control = control_connection(conn, run_id)
    try:
        run = engine_run(conn, run_id, control)
        if _cut_row(control, run, stage) is not None:
            return unit_rows(control, run, stage)
        summary = {"kind": kind, "units": len(units), "cost": cost, "detail": plan_detail}
        try:
            ledger.add_jobs(control, run, [{"stage": stage, "command": OPEN, "arguments": {"index": index, "plan": summary}}]
                            + [{"stage": stage, "command": kind, "item_key": str(number), "arguments": {"args": args, "cost": unit_cost}}
                               for number, (args, unit_cost) in enumerate(units)]
                            + [{"stage": stage, "command": command, "needs": [before]}
                               for command, before in zip(END_COMMANDS, ("units", *END_COMMANDS))])
            control.execute("UPDATE engine_jobs SET state = 'done', updated_at = NOW(6) WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''", (run, stage, OPEN))
            control.commit()
        except Exception:
            control.rollback()
            raise
        return unit_rows(control, run, stage)
    finally:
        control.close()


def _cut_row(control, run, stage):
    return control.execute("SELECT arguments FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = ? AND item_key = '' AND state = 'done'", (run, stage, OPEN)).fetchone()


def unit_rows(control, run, stage):
    """The ids of the unit rows of the stage, in unit order."""
    rows = control.execute("SELECT id, item_key FROM engine_jobs WHERE run_id = ? AND stage = ? AND item_key <> '' ORDER BY CAST(item_key AS UNSIGNED)", (run, stage)).fetchall()
    return [row["id"] for row in rows]


def read_plan(conn, run_id, index):
    """Stage `index`'s plan summary (`{"kind", "units", "cost", "detail"}`), or `None` while the stage has not been cut (or the run has no ledger run)."""
    keep = (runrecord.PREFIX, run_id, index)
    if keep in _PLANS:
        return _PLANS[keep]
    found = conn.hget(runrecord.run_key(run_id), ENGINE_RUN_FIELD)
    if not found or index >= len(runrecord.declared_stages(conn, run_id)):
        return None
    control = control_connection(conn, run_id)
    try:
        row = _cut_row(control, int(runrecord._s(found)), _stage_key(conn, run_id, index))
    finally:
        control.close()
    if row is None:
        return None
    _PLANS[keep] = json.loads(row["arguments"])["plan"]
    return _PLANS[keep]


def ended(conn, run_id, index):
    """The end jobs of stage `index` that are done (`END_COMMANDS`, a set): what the end of a stage has got through. Empty for a stage with no rows."""
    found = conn.hget(runrecord.run_key(run_id), ENGINE_RUN_FIELD)
    if not found or index >= len(runrecord.declared_stages(conn, run_id)):
        return set()
    control = control_connection(conn, run_id)
    try:
        rows = control.execute("SELECT command FROM engine_jobs WHERE run_id = ? AND stage = ? AND item_key = '' AND state = 'done' AND command IN (?, ?, ?)",
                               (int(runrecord._s(found)), _stage_key(conn, run_id, index), *END_COMMANDS)).fetchall()
    finally:
        control.close()
    return {row["command"] for row in rows}


def mark_ended(conn, run_id, index, command):
    """
    Marks end job `command` of stage `index` done (idempotent; it commits). The stage's end rows are made first when its cut did not (a run begun by hand, whose stages nobody cut: the rows the
    cut writes, added the same way, so a cut that comes later finds them).
    """
    control = control_connection(conn, run_id)
    try:
        run, stage = engine_run(conn, run_id, control), _stage_key(conn, run_id, index)
        ledger.add_jobs(control, run, [{"stage": stage, "command": name, "needs": [before]} for name, before in zip(END_COMMANDS, ("units", *END_COMMANDS))])
        control.execute("UPDATE engine_jobs SET state = 'done', updated_at = NOW(6) WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''", (run, stage, command))
        control.commit()
    finally:
        control.close()
