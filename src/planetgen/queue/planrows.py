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
END_ROLES = (CONCLUDE, FINISH, RELEASE)
"""tuple[str]: The jobs that end a stage (OPS.218), each with an edge to the one before: `conclude` needs every unit, `finish` needs `conclude`, and `release` (the next opening queued) needs `finish`.
They are rows of the stage with no item key, as the opening is. The command of the `conclude` and the `finish` row is the stage's own end command where its kind has one (`commands.end_name`,
OPS.260), else the role's name; `release` is always `release`."""


def end_rows(kind):
    """`{role: the command of its row}` for the end jobs of a stage of `kind` (`None`: a run whose kind is not known has the roles' own names)."""
    from planetgen.queue import commands
    return {role: (role if kind is None or role == RELEASE else commands.end_name(kind, role)) for role in END_ROLES}


def _kind_of(conn, run_id, index):
    """The command of stage `index` of the run, from the steps its record holds (`None` for a run begun by hand, which holds none)."""
    key = _stage_key(conn, run_id, index)
    return next((item["kind"] for item in runrecord.read_steps(conn, run_id) if item["key"] == key), None)

AFTER = "after:"
"""str: What a `needs` entry of a stage's `open` row starts with when it names a stage: `after:<stage key>` is met when that stage's `release` row is done."""

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


def record_plan(conn, run_id, steps):
    """
    Writes the whole run as rows when it starts (OPS.238, the stored plan): one `open` row for each step, waiting, which `needs` the stages it comes after (`after:<stage key>`: the step before it, and
    the earlier steps that make what it needs, `needs.required_stage_keys`) and `makes` the artifacts its command makes (`needs.makes_of`). The cut of a stage (`record_cut`) later fills the
    row in and marks it done, so a run's order can be read from the rows before its first stage opens. A plan that was written already is left as it is.

    Args:
        steps (list[commands.Step]): The run's steps, in order, as `wrapper.run` has them.
    """
    from planetgen.queue import needs
    control = control_connection(conn, run_id)
    try:
        run = engine_run(conn, run_id, control)
        jobs = []
        for index, step in enumerate(steps):
            earlier = ([steps[index - 1].key] if index else []) + needs.required_stage_keys(steps, index)
            jobs.append({"stage": step.key, "command": OPEN, "arguments": {"index": index, "kind": step.kind, "label": step.label, "planned": True},
                         "needs": [AFTER + key for key in dict.fromkeys(earlier)], "makes": list(needs.makes_of(step.kind))})
        ledger.add_jobs(control, run, jobs)
        control.commit()
    except Exception:
        control.rollback()
        raise
    finally:
        control.close()


def record_cut(conn, run_id, index, kind, units, cost, plan_detail, split=None):
    """
    Writes the cut of stage `index` as rows, in one transaction, and returns the unit rows' ids in unit order.

    Args:
        kind (str): The command of the stage (the kind of each unit).
        units (list[tuple]): `(args, cost)` per unit, as the stage's expansion gave them.
        cost (float): The units' total predicted cost.
        plan_detail: What the expansion said of the stage as a whole (`commands.Expansion.plan`).
        split (dict): The cut the queue engine chose (OPS.161, `splitter.Decision.as_dict`), kept in the plan summary so a restarted unit has the same one; `None` when none was chosen.

    A stage already cut (its `open` row is done) is not cut again: the rows it has are returned, whatever this call was given.
    """
    stage = _stage_key(conn, run_id, index)
    control = control_connection(conn, run_id)
    try:
        run = engine_run(conn, run_id, control)
        if _cut_row(control, run, stage) is not None:
            return unit_rows(control, run, stage)
        summary = {"kind": kind, "units": len(units), "cost": cost, "detail": plan_detail}
        if split is not None:
            summary["split"] = split
        try:
            ledger.add_jobs(control, run, [{"stage": stage, "command": OPEN, "arguments": {"index": index, "plan": summary}}]
                            + [{"stage": stage, "command": kind, "item_key": str(number), "arguments": {"args": args, "cost": unit_cost}}
                               for number, (args, unit_cost) in enumerate(units)]
                            + _end_jobs(stage, kind))
            control.execute("UPDATE engine_jobs SET state = 'done', arguments = ?, updated_at = NOW(6) WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''",
                            (ledger._text({"index": index, "plan": summary}), run, stage, OPEN))        # a stage the plan wrote at the run's start (`record_plan`) has its row filled in here
            control.commit()
        except Exception:
            control.rollback()
            raise
        return unit_rows(control, run, stage)
    finally:
        control.close()


def _end_jobs(stage, kind):
    """The rows of the end of a stage of `kind`, each needing the one before (the first needs the units)."""
    names = [end_rows(kind)[role] for role in END_ROLES]
    return [{"stage": stage, "command": command, "needs": [before]} for command, before in zip(names, ("units", *names))]


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
    """The end jobs of stage `index` that are done (roles of `END_ROLES`, a set): what the end of a stage has got through. Empty for a stage with no rows."""
    found = conn.hget(runrecord.run_key(run_id), ENGINE_RUN_FIELD)
    if not found or index >= len(runrecord.declared_stages(conn, run_id)):
        return set()
    names = end_rows(_kind_of(conn, run_id, index))
    control = control_connection(conn, run_id)
    try:
        rows = control.execute("SELECT command FROM engine_jobs WHERE run_id = ? AND stage = ? AND item_key = '' AND state = 'done' AND command IN (?, ?, ?)",
                               (int(runrecord._s(found)), _stage_key(conn, run_id, index), *names.values())).fetchall()
    finally:
        control.close()
    done = {row["command"] for row in rows}
    return {role for role, command in names.items() if command in done}


def mark_ended(conn, run_id, index, role):
    """
    Marks end job `role` (`CONCLUDE`, `FINISH` or `RELEASE`) of stage `index` done (idempotent; it commits). The stage's end rows are made first when its cut did not (a run begun by hand, whose stages nobody cut: the rows the
    cut writes, added the same way, so a cut that comes later finds them).
    """
    control = control_connection(conn, run_id)
    try:
        run, stage = engine_run(conn, run_id, control), _stage_key(conn, run_id, index)
        kind = _kind_of(conn, run_id, index)
        ledger.add_jobs(control, run, _end_jobs(stage, kind))
        control.execute("UPDATE engine_jobs SET state = 'done', updated_at = NOW(6) WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''", (run, stage, end_rows(kind)[role]))
        control.commit()
    finally:
        control.close()
