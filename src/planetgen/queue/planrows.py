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

* **A run with no usable ledger** (no control database, or one with no ledger tables yet: `generate system --output` touches none, TEST.769, and the update's math check runs before the migration, OPS.263) has none to write to. Its plan is kept in the run's own Redis record, which every process of the run
  shares, behind the same functions: `record_plan`, `record_cut`, `read_plan`, `ended` and `mark_ended` choose the store once, by `ledgered`, and give the same answers; its units have no row
  (`record_cut` returns `None` for each), so no attempt of theirs is written (`attempts.begin` takes `None`).
"""

import json

from planetgen.queue import ledger, runrecord
from planetgen.util import log

OPEN = "open"
"""str: The command of the row that is the cut of a stage."""

CONCLUDE, FINISH, RELEASE = "conclude", "finish", "release"
END_ROLES = (CONCLUDE, FINISH, RELEASE)
"""tuple[str]: The jobs that end a stage (OPS.218), each with an edge to the one before: `conclude` needs every unit, `finish` needs `conclude`, and `release` (the next opening queued) needs `finish`.
They are rows of the stage with no item key, as the opening is. The command of the `conclude` and the `finish` row is the stage's own end command where its kind has one (`commands.end_name`,
OPS.260), else the role's name; `release` is always `release`."""


LEDGER_FIELD = "ledger"
"""str: The field of the run record that holds which store the run's plan is in (`1` the job ledger, `0` the run's own record), chosen once, by the first call."""


def _ledger_usable(conn, run_id):
    """Whether the control database answers and has the ledger's tables: not so on a machine with no database, or during `update.sh` before its migration step (the math check runs first)."""
    try:
        control = control_connection(conn, run_id)
    except Exception as exc:  # noqa: BLE001 -- no way to reach the control database is "no ledger"
        log.swallowed("queue/planrows", "open the control database", exc)
        return False
    try:
        control.execute("SELECT 1 FROM engine_runs LIMIT 1").fetchall()
        return True
    except Exception as exc:  # noqa: BLE001 -- a missing table (1146) or any other error is "no ledger"
        log.swallowed("queue/planrows", "probe the ledger", exc)
        return False
    finally:
        control.close()


def _names_no_database(conn, run_id):
    """Whether the run says it touches no database: every step has `database` in its parameters, as `None` (`generate system --output`); a step with no `database` at all is a run begun by hand or a test."""
    params = [item.get("params") or {} for item in runrecord.read_steps(conn, run_id)]
    return bool(params) and all("database" in each and not each["database"] for each in params)


def ledgered(conn, run_id):
    """
    Whether the run's plan is rows of the job ledger. It is when the control database is there and has the ledger's tables; otherwise (OPS.263, TEST.769: the update's math check runs before the migration that
    makes the tables, and `generate system --output` has no database) the plan is kept in the run's own Redis record. Chosen by the first call and kept in the record, so a run never changes store.
    """
    key = runrecord.run_key(run_id)
    found = conn.hget(key, LEDGER_FIELD)
    if found is None:
        conn.hsetnx(key, LEDGER_FIELD, "1" if not _names_no_database(conn, run_id) and _ledger_usable(conn, run_id) else "0")
        found = conn.hget(key, LEDGER_FIELD)
    return runrecord._s(found) == "1"


def _kept(conn, run_id, field, value=None):
    """Reads (`value` is `None`) or writes the Redis field `field` of the run's record that holds a ledgerless run's plan."""
    key = runrecord.run_key(run_id)
    if value is None:
        found = conn.hget(key, field)
        return None if found is None else json.loads(runrecord._s(found))
    conn.hset(key, field, json.dumps(value))
    return value


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
            from planetgen.db import store
            return store.MySQLConfig.from_args(named).database
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


def made_artifacts(control, galaxy, run=None):
    """
    The artifacts the ledger's rows say the galaxy `galaxy` has (`needs.exists_after`): what the stages that ended (their `release` row is done) made, in the order they ran, forgetting what came before a reset or a new skeleton.
    The rows of run `run` (the run being planned) are left out. A galaxy with no rows has none, as does one whose runs began before the ledger kept plans.
    """
    from planetgen.queue import needs
    if galaxy is None:
        return set()
    rows = control.execute(
        "SELECT o.makes, o.arguments FROM engine_jobs o JOIN engine_runs r ON r.id = o.run_id"
        " WHERE r.galaxy = ? AND o.run_id <> ? AND o.command = ? AND o.item_key = '' AND o.state = 'done'"
        " AND EXISTS (SELECT 1 FROM engine_jobs e WHERE e.run_id = o.run_id AND e.stage = o.stage AND e.command = ? AND e.item_key = '' AND e.state = 'done')"
        " ORDER BY o.run_id, o.id", (galaxy, run if run is not None else 0, OPEN, RELEASE)).fetchall()
    stages = []
    for row in rows:
        kind = (json.loads(row["arguments"] or "{}").get("plan") or {}).get("kind")
        stages.append((kind, json.loads(row["makes"] or "[]")))
    return needs.exists_after(stages)


def record_plan(conn, run_id, steps):
    """
    Writes the whole run as rows when it starts (OPS.238, the stored plan): one `open` row for each step, waiting, which `needs` the stages it comes after (`after:<stage key>`: the earlier steps it may not overlap with, `touches.plan_edges` (the step before it, until `independence.json`
    lists a pair), and the earlier steps that make what it needs, `needs.required_stage_keys`) and `makes` the artifacts its command makes (`needs.makes_of`). The cut of a stage (`record_cut`) later fills the
    row in and marks it done, so a run's order can be read from the rows before its first stage opens. A plan that was written already is left as it is.

    Args:
        steps (list[commands.Step]): The run's steps, in order, as `wrapper.run` has them.

    Returns:
        dict: `{step key: [artifacts]}` for the steps that need what no earlier step makes and the galaxy's rows do not say exists (`needs.unmet` against `made_artifacts`); each such row also says so in its `unmet`.
    """
    from planetgen.queue import needs, touches
    if not ledgered(conn, run_id):
        return {}       # a run with no ledger keeps no rows (and so no unmet needs to say); its plan is read from its steps
    control = control_connection(conn, run_id)
    try:
        run = engine_run(conn, run_id, control)
        lacking = {}
        for step, artifact in needs.unmet(steps, made_artifacts(control, _galaxy_name(conn, run_id), run)):
            lacking.setdefault(step.key, []).append(artifact)
        jobs = []
        waits = touches.plan_edges([step.kind for step in steps])
        for index, step in enumerate(steps):
            earlier = [steps[place].key for place in waits[index]] + needs.required_stage_keys(steps, index)
            extra = {"unmet": lacking[step.key]} if step.key in lacking else {}
            jobs.append({"stage": step.key, "command": OPEN, "arguments": {"index": index, "kind": step.kind, "label": step.label, "planned": True, **extra},
                         "needs": [AFTER + key for key in dict.fromkeys(earlier)], "makes": list(needs.makes_of(step.kind))})
        ledger.add_jobs(control, run, jobs)
        control.commit()
    except Exception:
        control.rollback()
        raise
    finally:
        control.close()
    return lacking


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
    if not ledgered(conn, run_id):
        return _cut_kept(conn, run_id, index, kind, units, cost, plan_detail, split)
    control = control_connection(conn, run_id)
    try:
        run = engine_run(conn, run_id, control)
        if _cut_row(control, run, stage) is not None:
            return unit_rows(control, run, stage)
        summary = {"kind": kind, "units": len(units), "cost": cost, "detail": plan_detail}
        if split is not None:
            summary["split"] = split
        try:
            ledger.add_jobs(control, run, [{"stage": stage, "command": OPEN, "arguments": {"index": index, "plan": summary, "run_id": run_id}}]
                            + [{"stage": stage, "command": kind, "item_key": str(number), "arguments": {"args": args, "cost": unit_cost}}
                               for number, (args, unit_cost) in enumerate(units)]
                            + _end_jobs(stage, kind))
            control.execute("UPDATE engine_jobs SET state = 'done', arguments = ?, updated_at = NOW(6) WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''",
                            (ledger._text({"index": index, "plan": summary, "run_id": run_id}), run, stage, OPEN))        # a stage the plan wrote at the run's start (`record_plan`) has its row filled in here
            control.commit()
        except Exception:
            control.rollback()
            raise
        return unit_rows(control, run, stage)
    finally:
        control.close()


def _cut_kept(conn, run_id, index, kind, units, cost, plan_detail, split):
    """`record_cut` of a run with no ledger: the plan summary is kept in the run's record, once (a stage already cut is not cut again), and the units have no rows."""
    field = f"plan:{index}"
    kept = _kept(conn, run_id, field)
    if kept is None:
        kept = {"kind": kind, "units": len(units), "cost": cost, "detail": plan_detail}
        if split is not None:
            kept["split"] = split
        _kept(conn, run_id, field, kept)
    return [None] * kept["units"]


def _end_jobs(stage, kind):
    """The rows of the end of a stage of `kind`, each needing the one before (the first needs the units)."""
    names = [end_rows(kind)[role] for role in END_ROLES]
    return [{"stage": stage, "command": command, "needs": [before]} for command, before in zip(names, ("units", *names))]


def release_units(conn, run_id, index, skip=()):
    """
    Makes the unit rows of stage `index` that are still `waiting` `ready` for the launcher to hand out (OPS.228), except the units numbered in `skip` (the ones a run this one carries on from did). Returns how many it released.
    The caller has cut the stage (`record_cut`); a run with no ledger has no rows and releases none.
    """
    if not ledgered(conn, run_id):
        return 0
    control = control_connection(conn, run_id)
    try:
        run, stage = engine_run(conn, run_id, control), _stage_key(conn, run_id, index)
        marks = ", ".join(["?"] * len(skip))
        released = control.execute("UPDATE engine_jobs SET state = 'ready', updated_at = NOW(6) WHERE run_id = ? AND stage = ? AND item_key <> '' AND state = 'waiting'"
                                   + (f" AND item_key NOT IN ({marks})" if skip else ""), (run, stage, *map(str, skip))).rowcount
        control.commit()
        return released
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
    if not ledgered(conn, run_id):
        kept = _kept(conn, run_id, f"plan:{index}")
        if kept is not None:
            _PLANS[keep] = kept
        return kept
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
    if not ledgered(conn, run_id):
        return set(_kept(conn, run_id, f"ended:{index}") or ())
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
    if not ledgered(conn, run_id):
        _kept(conn, run_id, f"ended:{index}", sorted(set(_kept(conn, run_id, f"ended:{index}") or ()) | {role}))
        return
    control = control_connection(conn, run_id)
    try:
        run, stage = engine_run(conn, run_id, control), _stage_key(conn, run_id, index)
        kind = _kind_of(conn, run_id, index)
        ledger.add_jobs(control, run, _end_jobs(stage, kind))
        control.execute("UPDATE engine_jobs SET state = 'done', updated_at = NOW(6) WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''", (run, stage, end_rows(kind)[role]))
        control.commit()
    finally:
        control.close()
