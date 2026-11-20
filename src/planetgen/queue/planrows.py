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
* **The windows (OPS.240).** A stage whose command has `window` (a sector for each unit, whatever the galaxy's size) is cut into the rows of its first window only, and a `spawn` row (`SPAWN`, `waiting`; its arguments
  are the cursor `next` and the `total`) says there are more. The launcher's turn calls `spawn_window` for each such row: when the units of the stage that are not done fall under `LOW_WATER` of a window it writes the
  next window of unit rows, in one transaction with the cursor's move, so a restart finds the cursor where the rows are. The stage's `units` (the plan summary, the stage's hash) is the total, so the follower, the
  bar and the end of the stage count the units not yet written too.

The ledger lives in the control database, found from the database the run's steps name (`wrapper.expected_seconds` reads the same one).
"""

import json

from planetgen.queue import itemfile, ledger, runrecord
from planetgen.util import log

OPEN = "open"
"""str: The command of the row that is the cut of a stage."""

SPAWN = "spawn"
"""str: The command of the row that writes the later windows of a stage's unit rows (OPS.240). It stays `waiting` while there are windows to write, and is `done` when the last is written, `canceled` when the run is over first."""

WINDOW_PER_WORKER = 4
"""int: How many unit rows a window holds for each worker the stage may run at once."""

WINDOW_MIN, WINDOW_MAX = 200, 2000
"""tuple[int]: The least and the most unit rows of a window, whatever the workers."""

LOW_WATER = 0.5
"""float: The next window is written when the stage's units not yet done are fewer than this share of a window."""


FOLD_ABOVE = 10_000
"""int: A stage with more unit rows than this has the history of its finished rows folded when it ends (`fold_finished`, OPS.240)."""

FOLD_KEEP = 100
"""int: The slowest finished units of a folded stage whose rows and attempts stay whole."""

FOLD_BATCH = 20_000
"""int: How many job ids one statement of the fold reaches."""


def stage_window(conn, run_id, index):
    """The window of the stage `index` of the run (`window_size` of the workers it was declared with: the run's own limit, else its worker count)."""
    declared = runrecord.declared_stages(conn, run_id)[index]
    return window_size(declared.get("cap") or declared.get("workers"))


def window_size(workers):
    """How many unit rows a window of a stage that may run `workers` units at once holds: `WINDOW_PER_WORKER` of each, between `WINDOW_MIN` and `WINDOW_MAX`."""
    return max(WINDOW_MIN, min(WINDOW_MAX, WINDOW_PER_WORKER * max(1, int(workers or 1))))

CONCLUDE, FINISH, RELEASE = "conclude", "finish", "release"
END_ROLES = (CONCLUDE, FINISH, RELEASE)
"""tuple[str]: The jobs that end a stage (OPS.218), each with an edge to the one before: `conclude` needs every unit, `finish` needs `conclude`, and `release` (the next opening queued) needs `finish`.
They are rows of the stage with no item key, as the opening is. The command of the `conclude` and the `finish` row is the stage's own end command where its kind has one (`commands.end_name`,
OPS.260), else the role's name; `release` is always `release`."""


def end_rows(kind):
    """`{role: the command of its row}` for the end jobs of a stage of `kind` (`None`: a run whose kind is not known has the roles' own names)."""
    from planetgen.queue import commands
    return {role: (role if kind is None or role == RELEASE else commands.end_name(kind, role)) for role in END_ROLES}


JOB_CLASSES = {runrecord.GENERATION: "bulk", runrecord.MAINTENANCE: "maintenance", runrecord.CHECK: "maintenance"}
"""dict: The job class of the rows of a run of each kind (OPS.229): the warm-up and the checks are the maintenance lane, which the launcher claims on its own lane and holds to its own cap."""


def job_class(conn, run_id):
    """The class of the rows of the run: `JOB_CLASSES` of its kind (`bulk` for a kind it does not know)."""
    return JOB_CLASSES.get(runrecord._s(conn.hget(runrecord.run_key(run_id), "kind")), "bulk")


def _classed(jobs, name):
    return [{**job, "job_class": name} for job in jobs]


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
        from planetgen.queue import runstore
        runstore.persist(conn, run_id)          # OPS.284: from here the run has a durable copy
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
        ledger.add_jobs(control, run, _classed(jobs, job_class(conn, run_id)))
        control.commit()
    except Exception:
        control.rollback()
        raise
    finally:
        control.close()
    return lacking


def record_cut(conn, run_id, index, kind, units, cost, plan_detail, split=None, total=None, later=None, window=None, prior=None):
    """
    Writes the cut of stage `index` as rows, in one transaction, and returns the unit rows' ids in unit order.

    Args:
        kind (str): The command of the stage (the kind of each unit).
        units (list[tuple]): `(args, cost)` per unit, as the stage's expansion gave them.
        cost (float): The units' total predicted cost.
        plan_detail: What the expansion said of the stage as a whole (`commands.Expansion.plan`).
        split (dict): The cut the queue engine chose (OPS.161, `splitter.Decision.as_dict`), kept in the plan summary so a restarted unit has the same one; `None` when none was chosen.
        total (int): How many units the stage has in all when `units` is only its first window (OPS.240); `None` for `len(units)`.
        later (dict): What the command's `window` needs to say the later units (`commands.Expansion.later`); with it a `spawn` row is written that holds the cursor.
        window (int): How many unit rows each later window holds (`window_size`); `units` was cut to the same size.
        prior (str): The run this one carries on from (`resumed.RESUMES_FROM`), whose done units the later windows subtract as the first one did.

    A stage already cut (its `open` row is done) is not cut again: the rows it has are returned, whatever this call was given.
    """
    stage = _stage_key(conn, run_id, index)
    control = control_connection(conn, run_id)
    try:
        run = engine_run(conn, run_id, control)
        if _cut_row(control, run, stage) is not None:
            return unit_rows(control, run, stage)
        summary = {"kind": kind, "units": len(units) if total is None else int(total), "cost": cost, "detail": plan_detail}
        if split is not None:
            summary["split"] = split
        cap = runrecord.declared_stages(conn, run_id)[index].get("cap")      # the most units of it that run at once, when the run was asked
        opened = {"index": index, "plan": summary, "run_id": run_id, "cap": cap}
        spawn = []
        if later is not None and summary["units"] > len(units):
            spawn = [{"stage": stage, "command": SPAWN, "arguments": {"run_id": run_id, "index": index, "kind": kind, "next": len(units), "total": summary["units"], "window": int(window), "later": later,
                                                                           "released": len(units), "prior": prior}}]
        try:
            ledger.add_jobs(control, run, _classed([{"stage": stage, "command": OPEN, "arguments": opened}]
                                                    + [{"stage": stage, "command": kind, "item_key": str(number), "arguments": {"args": args, "cost": unit_cost}}
                                                       for number, (args, unit_cost) in enumerate(units)]
                                                    + spawn + _end_jobs(stage, kind), job_class(conn, run_id)))
            control.execute("UPDATE engine_jobs SET state = 'done', arguments = ?, updated_at = NOW(6) WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''",
                            (ledger._text(opened), run, stage, OPEN))        # a stage the plan wrote at the run's start (`record_plan`) has its row filled in here
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


OPENING = "opening"
"""str: The command of the row that opens a stage (OPS.228): the job that cuts stage `index` into its unit rows and releases them. Its arguments are `{run_id, index}`; the stage is the stage's key."""


def queue_opening(conn, run_id, index):
    """
    Makes the opening of stage `index` a job the supervisor's launcher hands out (OPS.228), in place of the RQ job `units.enqueue_open` queues: the stage's `opening` row, added if it is not there and made `ready`
    (a row of an opening that ended or failed is handed out again, as a queued job of it would be). Claims the cut's memory (`claims.CUT_MB`). Returns the row's id.
    """
    from planetgen.queue import claims
    control = control_connection(conn, run_id)
    try:
        run, stage = engine_run(conn, run_id, control), _stage_key(conn, run_id, index)
        ledger.add_jobs(control, run, [{"stage": stage, "command": OPENING, "arguments": {"run_id": run_id, "index": index}, "claim_mb": claims.CUT_MB, "job_class": job_class(conn, run_id)}])
        control.execute("UPDATE engine_jobs SET state = 'ready', waited = NULL, updated_at = NOW(6) WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''", (run, stage, OPENING))
        found = control.execute("SELECT id FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''", (run, stage, OPENING)).fetchone()["id"]
        control.commit()
        return found
    finally:
        control.close()


def opening_pending(conn, run_id, index):
    """Whether the opening of stage `index` is a row the launcher has not finished with (`ready` for it to hand out, or `running`)."""
    control = control_connection(conn, run_id)
    try:
        found = conn.hget(runrecord.run_key(run_id), ENGINE_RUN_FIELD)
        if not found:
            return False
        row = control.execute("SELECT state FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''", (int(runrecord._s(found)), _stage_key(conn, run_id, index), OPENING)).fetchone()
        return row is not None and row["state"] in ("ready", "running")
    finally:
        control.close()


def release_units(conn, run_id, index, skip=()):
    """
    Makes the unit rows of stage `index` that are still `waiting` `ready` for the launcher to hand out (OPS.228), except the units numbered in `skip` (the ones a run this one carries on from did). Returns how many it released.
    The caller has cut the stage (`record_cut`).
    """
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


def add_to_unit_args(conn, run_id, index, extra):
    """
    Adds the keys of `extra` to the `args` of the unit rows of stage `index` that are still `waiting` (OPS.228): what a run that carries on from an earlier one tells a unit that cannot resume (`resumed.AGAIN`),
    which the RQ path puts on the arguments it queues. Returns how many rows it changed.
    """
    control = control_connection(conn, run_id)
    try:
        run, stage = engine_run(conn, run_id, control), _stage_key(conn, run_id, index)
        changed = 0
        for row in control.execute("SELECT id, arguments FROM engine_jobs WHERE run_id = ? AND stage = ? AND item_key <> '' AND state = 'waiting'", (run, stage)).fetchall():
            body = json.loads(row["arguments"]) if row["arguments"] else {}
            body["args"] = {**body.get("args", {}), **extra}
            control.execute("UPDATE engine_jobs SET arguments = ? WHERE id = ?", (json.dumps(body), row["id"]))
            changed += 1
        control.commit()
        return changed
    finally:
        control.close()


def spawn_rows(control):
    """The `spawn` rows still `waiting` (OPS.240), oldest first: `(row id, ledger run, stage, arguments)` of each, the arguments as the dict the row holds."""
    found = control.execute("SELECT id, run_id, stage, arguments FROM engine_jobs WHERE command = ? AND state = 'waiting' AND item_key = '' ORDER BY id", (SPAWN,)).fetchall()
    return [(row["id"], row["run_id"], row["stage"], json.loads(row["arguments"] or "{}")) for row in found]


def _window_units(arguments, first, count):
    """The `(args, cost)` of units `first` up to `first + count` of the stage a `spawn` row's arguments describe, from its command's `window`; fewer than asked is a stage that lost its items (an error)."""
    from planetgen.queue import commands
    units = commands.window_function(arguments["kind"])(arguments["later"], first, count)
    if len(units) != count:
        raise RuntimeError(f"stage {arguments['index']} of run {arguments['run_id']}: units {first:,} to {first + count:,} were asked of its items and {len(units):,} came back")
    return units


def _release_window(conn, control, run, stage, row_id, arguments):
    """
    Carries (when the run carries on from another), tells and releases the unit rows written but not released yet (`released` up to `next` in the `spawn` row's arguments), then moves `released`. Run again after a
    stop, it does the same to the same rows: a row already `ready` or `done` is not touched.
    """
    from planetgen.queue import resumed
    run_id, index = arguments["run_id"], arguments["index"]
    first, last = arguments["released"], arguments["next"]
    if first >= last:
        return
    found = resumed.Carried()
    prior = arguments.get("prior")
    if prior:
        units = _window_units(arguments, first, last - first)
        rows = {int(item["item_key"]): item["id"] for item in control.execute(
            "SELECT id, item_key FROM engine_jobs WHERE run_id = ? AND stage = ? AND item_key IN (" + ", ".join(["?"] * (last - first)) + ")", (run, stage, *map(str, range(first, last)))).fetchall()}
        found = resumed.carry(conn, run_id, index, prior, arguments["kind"], [rows[unit] for unit in range(first, last)], units, first=first)
        resumed.settle(conn, run_id, index, found)
    from planetgen.queue import commands
    if prior and not commands.spec(arguments["kind"]).resumes:
        add_to_unit_args(conn, run_id, index, {resumed.AGAIN: True})
    release_units(conn, run_id, index, skip=sorted(found.done))
    arguments["released"] = last
    control.execute("UPDATE engine_jobs SET arguments = ?, updated_at = NOW(6) WHERE id = ?", (ledger._text(arguments), row_id))
    control.commit()


def spawn_window(conn, run_id, index):
    """
    Writes the next window of the unit rows of stage `index` when the stage wants one (OPS.240): the units written and not done fewer than `LOW_WATER` of a window and units left to write. The rows and the cursor's
    move are one transaction; the carry and the release of the rows that follow are run again by the next call when this one was stopped between them. A run that is not live gets none (a run that is over has its
    `spawn` row canceled). Returns how many unit rows it wrote.
    """
    control = control_connection(conn, run_id)
    try:
        run, stage = engine_run(conn, run_id, control), _stage_key(conn, run_id, index)
        found = control.execute("SELECT id, arguments FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = ? AND item_key = '' AND state = 'waiting'", (run, stage, SPAWN)).fetchone()
        if found is None:
            return 0
        arguments = json.loads(found["arguments"])
        state = runrecord.run_state(conn, run_id)
        if state not in runrecord.LIVE_STATES:
            if state != runrecord.HELD:
                ledger.set_job_state(control, found["id"], "canceled", waited="the run is over")
                control.commit()
                itemfile.remove(itemfile.path_of(run_id, stage))
            return 0
        _release_window(conn, control, run, stage, found["id"], arguments)
        total, written, window = arguments["total"], arguments["next"], arguments["window"]
        if written < total and written - conn.scard(runrecord.done_key(run_id, index)) - conn.hlen(runrecord.held_key(run_id, index)) < LOW_WATER * window:      # a held unit has run as far as it will: it does not keep a window from being written (the run is held when the last is done or held, `held.settle`)
            count = min(window, total - written)
            units = _window_units(arguments, written, count)
            ledger.add_jobs(control, run, _classed([{"stage": stage, "command": arguments["kind"], "item_key": str(written + number), "arguments": {"args": args, "cost": cost}}
                                                    for number, (args, cost) in enumerate(units)], job_class(conn, run_id)))
            arguments["next"] = written + count
            control.execute("UPDATE engine_jobs SET arguments = ?, updated_at = NOW(6) WHERE id = ?", (ledger._text(arguments), found["id"]))
            control.commit()
            _release_window(conn, control, run, stage, found["id"], arguments)
            written = count
        else:
            written = 0
        if arguments["released"] >= total:
            ledger.set_job_state(control, found["id"], "done")
            control.commit()
            if conn.scard(runrecord.done_key(run_id, index)) >= total:       # every unit of the last window was carried done: no unit is left to end the stage
                from planetgen.queue import units, wrapper
                units._end_stage(conn, run_id, index, wrapper)
        return written
    finally:
        control.close()


def fold_finished(conn, run_id, index):
    """
    Folds the history of a stage that ended with more than `FOLD_ABOVE` unit rows (OPS.240): a unit row that finished at its first attempt loses its `arguments` (about a kilobyte of what it was told to do) and its
    attempt row, except the `FOLD_KEEP` slowest. Every row stays, so the stage's count, the status lines and the Queue page's totals are as real as before; a unit that was retried or failed is whole, and so is
    every row of a stage that did not end. What was folded (`folded`, `attempts`, the seconds the folded attempts took, and the ids kept) is written to the stage's `release` row. Returns the number of rows folded, 0
    for a stage that is too small. The stage must have ended: a row without arguments cannot be run again or compared with a new cut (`resumed._same` takes it as different).
    """
    control = control_connection(conn, run_id)
    try:
        run, stage = engine_run(conn, run_id, control), _stage_key(conn, run_id, index)
        where = "run_id = ? AND stage = ? AND item_key <> ''"
        found = control.execute(f"SELECT COUNT(*) AS n, MIN(id) AS low, MAX(id) AS high FROM engine_jobs WHERE {where}", (run, stage)).fetchone()
        if found["n"] <= FOLD_ABOVE or control.execute(f"SELECT 1 AS x FROM engine_jobs WHERE {where} AND state <> 'done' LIMIT 1", (run, stage)).fetchone():
            return 0
        slow = [row["job_id"] for row in control.execute(
            "SELECT a.job_id FROM engine_attempts a JOIN engine_jobs j ON j.id = a.job_id WHERE j.run_id = ? AND j.stage = ? AND j.item_key <> '' AND a.ended_at IS NOT NULL"
            " ORDER BY TIMESTAMPDIFF(MICROSECOND, a.started_at, a.ended_at) DESC LIMIT ?", (run, stage, FOLD_KEEP)).fetchall()]
        keep = ", ".join(["?"] * len(slow)) or "NULL"
        folded = attempts = 0
        seconds = 0.0
        for low in range(found["low"], found["high"] + 1, FOLD_BATCH):
            span = (run, stage, low, low + FOLD_BATCH - 1, *slow)
            rule = f"j.run_id = ? AND j.stage = ? AND j.item_key <> '' AND j.id BETWEEN ? AND ? AND j.state = 'done' AND j.attempt <= 1 AND j.arguments IS NOT NULL AND j.id NOT IN ({keep})"
            seconds += float(control.execute(f"SELECT COALESCE(SUM(TIMESTAMPDIFF(MICROSECOND, a.started_at, a.ended_at)), 0) / 1000000 AS s FROM engine_attempts a JOIN engine_jobs j ON j.id = a.job_id"
                                              f" WHERE {rule} AND a.ended_at IS NOT NULL", span).fetchone()["s"])
            attempts += control.execute(f"DELETE a FROM engine_attempts a JOIN engine_jobs j ON j.id = a.job_id WHERE {rule}", span).rowcount
            folded += control.execute(f"UPDATE engine_jobs j SET j.arguments = NULL WHERE {rule}", span).rowcount
        summary = {"folded": folded, "attempts": attempts, "seconds": round(seconds, 3), "kept": slow}
        control.execute("UPDATE engine_jobs SET arguments = ?, updated_at = NOW(6) WHERE run_id = ? AND stage = ? AND command = ? AND item_key = ''",
                        (ledger._text({"fold": summary}), run, stage, RELEASE))
        control.commit()
        return folded
    except Exception:
        control.rollback()
        raise
    finally:
        control.close()


def _cut_row(control, run, stage):
    return control.execute("SELECT arguments FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = ? AND item_key = '' AND state = 'done'", (run, stage, OPEN)).fetchone()


def unit_rows(control, run, stage):
    """The ids of the unit rows of the stage, in unit order."""
    rows = control.execute("SELECT id, item_key FROM engine_jobs WHERE run_id = ? AND stage = ? AND item_key <> '' ORDER BY CAST(item_key AS UNSIGNED)", (run, stage)).fetchall()
    return [row["id"] for row in rows]


def read_plan(conn, run_id, index):
    """Stage `index`'s plan summary (`{"kind", "units", "cost", "detail"}`), or `None` while the stage has not been cut."""
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
        if role == RELEASE:
            itemfile.remove(itemfile.path_of(run_id, stage))       # the stage's item file (OPS.240) has no more windows to give once the stage has ended
    finally:
        control.close()
    if role == RELEASE:
        try:
            fold_finished(conn, run_id, index)
        except Exception as exc:  # noqa: BLE001 -- the fold is the history's housekeeping: a stage that ended stays ended whatever it did
            log.swallowed("queue/planrows", "fold_finished", exc)
