# planetgen/queue/wrapper.py

"""
The wrapper of a `planetgen galaxy` or `planetgen plan` run (PERF.118, docs/design/operations/redis-and-the-work-queue.md 5.12):
the one place that knows the order of the run's commands, and the only code that starts one.

* **It selects.** `steps` is `commands.build_steps` for the parsed arguments: the commands that will run, in order,
  and nothing else; a command the request does not ask for is not in the list.
* **It writes the run record** (`runrecord.start_run`) with that list, declares a stage for each step that is queued
  (`queued`), and queues the first stage's units. Each stage is cut into its units by its command's `expand`, from the
  step's parameters and the database alone (`commands.Expansion`); the cut and a JSON summary of it (`plan`) are
  recorded before any unit is queued.
* **The done set is the barrier.** The unit that completes a stage calls `advance`, which expands and queues the next
  stage from the run record (`runrecord.read_steps`), so the order is read from data and no command knows what comes
  next.
* **The starter follows** (`run`): it reads each stage's plan from Redis, tells the caller (`on_plan`), and follows
  the units' results and progress (`units.follow_stage`). It is a client: a second process can follow the same run
  from the run record alone.

A step whose command is inline (`commands.CommandSpec.inline`) is not queued and has no stage in the run record; the
process that started the run does it, between the queued stages. No registered command is inline now (PERF.106 to
PERF.108 and PERF.121 to PERF.128 gave each its units), so every step of a run is a stage.
"""

import importlib
import json
import logging
import os
import time

import redis

from planetgen.datasets import fetch
from planetgen.queue import commands, progress_file, redisqueue, runline, runrecord, stagetimes, units, workers
from planetgen.util import log
from planetgen.util.settings import get_settings


JOB_ENV = "PLANETGEN_JOB_ID"
"""str: The variable a Generate page job's runner sets for its steps (`web.job_runner.JOB_ENV`): the run records the
job it belongs to, so the page finds the job's runs (PERF.99)."""


def database_args(mysql_config):
    """The database a unit opens, as plain JSON: its name, host, port and user. Never the password; a worker reads
    that from its own settings."""
    return {"host": mysql_config.host, "port": mysql_config.port, "user": mysql_config.user,
            "database": mysql_config.database}


def steps_for(command, args):
    """The steps of a `galaxy`, `plan`, `scatter` or one-stage run for the parsed `args` (`commands.build_steps`): the same list for the
    same arguments."""
    return commands.build_steps(command, args)


def step_for(args, key, label, params):
    """
    The step for stage `key` of the run `args` asks for, with `params` (what the stage needs besides the request's own
    settings: the database, and values worked out from the stored galaxy), and the run's whole list with it.

    The list is `steps_for` when `args` came from `planetgen galaxy` or `planetgen plan`. A stage that runs although
    the request does not select it (a luminosity redo with no scatter stored yet needs the massive stars too) is added
    to the list in its place, before the luminosity pass; arguments from neither command give a list of that step.
    """
    from planetgen.generation import stages
    command = getattr(args, "command", None)
    listing = steps_for(command, args) if command in commands.STAGED_COMMANDS else []
    found = next((step for step in listing if step.key == key), None)
    step = commands.Step(key, commands.KIND_OF_STAGE[key], found.label if found else label,
                         {**(found.params if found else stages.stage_settings(key, args)), **params})
    if found is not None:
        return step, [step if item.key == key else item for item in listing]
    later = next((index for index, item in enumerate(listing) if item.key == "luminosity"), len(listing))
    return step, listing[:later] + [step] + listing[later:]


def plan_steps(args, parts, settings=None):
    """
    The steps of a run that follows several stages (PERF.172), and the run's whole list with them: `parts` is
    `(key, label, params)` for each stage the run queues, in order, `params` what the stage needs besides the request's
    own settings (`step_for` for one). A stage the request's list holds takes that step's label and settings; one it
    does not (a luminosity redo with no scatter stored yet needs the massive stars too) is added to the list before the
    luminosity pass, where `step_for` puts it. `settings` (`{key: settings}`, PERF.203) are the settings of the stages of a run of
    several commands, each the ones its own command was given, in place of `stages.stage_settings(key, args)` of the run's.

    Returns:
        tuple: `(steps, listing)`.
    """
    from planetgen.generation import stages
    command = getattr(args, "command", None)
    listing = steps_for(command, args) if command in commands.STAGED_COMMANDS else []
    steps = []
    for key, label, params in parts:
        found = next((step for step in listing if step.key == key), None)
        own = (settings or {}).get(key)
        steps.append(commands.Step(key, commands.KIND_OF_STAGE[key], found.label if found else label,
                                   {**(found.params if found else own if own is not None
                                       else stages.stage_settings(key, args)), **params}))
    by_key = {step.key: step for step in steps}
    listing = [by_key.pop(item.key, item) if item.key in by_key else item for item in listing]
    if by_key:
        later = next((index for index, item in enumerate(listing) if item.key == "luminosity"), len(listing))
        listing = listing[:later] + list(by_key.values()) + listing[later:]
    return steps, listing


def queued(steps):
    """The steps of `steps` that are queued as stages of units (the rest are inline), in order."""
    return [step for step in steps if not commands.spec(step.kind).inline]


def _expand(conn, run_id, index):
    """
    Cuts stage `index` of the run into its units, records the cut, begins the stage and queues the units. The
    step comes from the run record, so any process can do this.

    Returns:
        int: How many units the stage queued.
    """
    stages = json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), "stages")))
    key = stages[index]["key"]
    step = next(commands.Step.from_dict(item) for item in runrecord.read_steps(conn, run_id) if item["key"] == key)
    started = time.time()       # the stage's time counts from before it was cut: reading the galaxy is part of it
    if step.params.get(FINISHED_BY):
        _carry_finished_stage(conn, run_id, index, step, started)
        return 0
    if _is_referred_to(runrecord.read_steps(conn, run_id), key):
        runrecord.set_opened_clock(conn, run_id, index, _database_clock(step.params))
    params = resolve_params(conn, run_id, step.params)
    prepare = commands.prepare_function(step.kind)
    if prepare is not None:
        prepare(params)
    expansion = commands.expand_function(step.kind)(params)
    cost = sum(unit_cost for _args, unit_cost in expansion.units)
    runrecord.set_stage_plan(conn, run_id, index, len(expansion.units), cost,
                             {"kind": step.kind, "units": len(expansion.units), "cost": cost,
                              "detail": expansion.plan})
    runrecord.begin_stage(conn, run_id, index, now=started)
    if not expansion.units and commands.conclude_function(step.kind) is not None:
        _conclude_empty_stage(conn, run_id, index, stage_label(conn, run_id, index))
        return 0
    units.enqueue_units(conn, run_id, index, step.kind, expansion.units,
                        queue_name=queue_of(runrecord._s(conn.hget(runrecord.run_key(run_id), "kind"))),
                        expected_seconds=json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), "stages")))[index]
                        .get("expected_seconds"))
    return len(expansion.units)


DEFERRED = "deferred"
"""str: The key of a deferred parameter (PERF.213): a step's parameter that is a dict with this key is not a value but the name of
a resolver in `RESOLVERS`, and the stage is given its value when it opens (`resolve_params`), not when the run is built."""


def opened_by(key):
    """A deferred parameter (PERF.213): the database's clock when stage `key` of the run opened, as text, so a stage that works on
    "the sectors the run saved" takes the clock of the sectors stage, not of the moment its followers were built."""
    return {DEFERRED: "stage-clock", "stage": key}


def is_deferred(value):
    """Whether a parameter value is a deferred one."""
    return isinstance(value, dict) and DEFERRED in value


def _is_referred_to(recorded, key):
    """Whether any step of the run (`runrecord.read_steps`) has a deferred parameter that needs the clock of stage `key`."""
    return any(is_deferred(value) and value.get("stage") == key and value[DEFERRED] == "stage-clock"
               for item in recorded for value in (item.get("params") or {}).values())


def _database_clock(params):
    """The database's own clock now, as text (`store.database_now`), for the stage whose `params` name the database."""
    from planetgen.db import store
    conn = store.get_connection(_mysql_config(params["database"]))
    try:
        return store.database_now(conn).isoformat(sep=" ")
    finally:
        conn.close()


class Unresolved(LookupError):
    """A deferred parameter whose value is not there (yet)."""


def _stage_clock(conn, run_id, marker):
    """The clock stage `marker["stage"]` recorded when it opened (`runrecord.set_opened_clock`)."""
    declared = json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), "stages")))
    keys = [runrecord._s(item["key"]) for item in declared]
    if marker["stage"] not in keys:
        raise Unresolved(f"the run has no stage {marker['stage']!r} to take the clock of")
    clock = runrecord.opened_clock(conn, run_id, keys.index(marker["stage"]))
    if clock is None:
        raise Unresolved(f"stage {marker['stage']!r} has not opened, so it has no clock")
    return clock


def plan_detail(key, field):
    """A deferred parameter (PERF.214): `field` of the summary stage `key` made of itself when it was cut (`commands.Expansion.plan`,
    what a worker worked out when the stage opened), so a stage after it reads what the stage planned."""
    return {DEFERRED: "plan-detail", "stage": key, "field": field}


def _plan_field(conn, run_id, marker):
    """The `marker["field"]` of the plan summary stage `marker["stage"]` recorded (`runrecord.set_stage_plan`)."""
    declared = json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), "stages")))
    keys = [runrecord._s(item["key"]) for item in declared]
    if marker["stage"] not in keys:
        raise Unresolved(f"the run has no stage {marker['stage']!r} to read the plan of")
    plan = runrecord.read_plan(conn, run_id, keys.index(marker["stage"]))
    detail = (plan or {}).get("detail") or {}
    if marker["field"] not in detail:
        raise Unresolved(f"stage {marker['stage']!r} has not planned {marker['field']!r}")
    return detail[marker["field"]]


RESOLVERS = {"stage-clock": _stage_clock, "plan-detail": _plan_field,
             "scatter-settings": "planetgen.generation.run_plan:resolve_scatter_settings"}
"""dict: Each deferred parameter's name to `resolver(conn, run_id, marker)`, which gives its value (or raises `Unresolved`), or to
`"module:function"` naming it, for a resolver in a module that imports this one. A command that needs a value worked out when its
stage opens adds its resolver here. A marker that has `MERGE` set stands for several parameters at once: the resolver gives a dict
and the stage's parameters take it (PERF.212)."""

MERGE = "merge"
"""str: The key of a deferred parameter whose value is a dict of parameters, which join the stage's (the marker's own name is then
only a label)."""


def deferred(name, merge=False, **fields):
    """A deferred parameter (PERF.212): the resolver `name` of `RESOLVERS`, given `fields` (plain JSON) to work from. With `merge`
    it stands for the dict of parameters the resolver gives."""
    return {DEFERRED: name, **({MERGE: True} if merge else {}), **fields}


def cached(conn, run_id, name, compute):
    """
    The value `compute()` gave the first time a stage of the run asked for `name`, kept in the run record (PERF.212). A value read from
    the database when a stage opens would not be the same later: the stage may empty what it read (the scatter's own record), and a
    resolver runs again at the stage's end and when a stopped run is carried on. The first stage to open computes it, before its
    own `prepare`, and the later ones, the conclude hooks and `finish_run` take what was kept.
    """
    field = f"deferred:{name}"
    kept = conn.hget(runrecord.run_key(run_id), field)
    if kept is not None:
        return json.loads(runrecord._s(kept))
    value = compute()
    conn.hsetnx(runrecord.run_key(run_id), field, json.dumps(value))
    return json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), field)))


def _resolver(name):
    """The function that resolves deferred parameter `name`."""
    found = RESOLVERS[name]
    if isinstance(found, str):
        module, function = found.split(":")
        found = getattr(importlib.import_module(module), function)
        RESOLVERS[name] = found
    return found


def scatter_settings_name(marker):
    """The name a `scatter-settings` marker's value is kept under in the run record (`cached`): its scatter's scope and tag, so the
    scatters of one run (a job holds several, PERF.203) keep theirs apart."""
    return marker["scope"] + marker.get("tag", "")


def resolve_params(conn, run_id, params, strict=True):
    """
    `params` with each deferred parameter (`opened_by`) replaced by its value as the run record has it now. With `strict` false a
    parameter that cannot be resolved yet stays as it is (`generation.finish_run` copies the steps of a run that stopped before
    some stage opened), and so does a `scatter-settings` parameter no stage has worked out yet: reading the database for a stage that
    has not opened would keep a value the stages before it may still change (a clear, a mass pass that stores the limit).

    Raises:
        Unresolved: A deferred parameter has no value and `strict` is true.
        RuntimeError: A resolver refuses (and `strict` is true).
    """
    out = {}
    for name, value in params.items():
        if is_deferred(value):
            if (not strict and value[DEFERRED] == "scatter-settings"
                    and conn.hget(runrecord.run_key(run_id), f"deferred:{scatter_settings_name(value)}") is None):
                out[name] = value
                continue
            try:
                found = _resolver(value[DEFERRED])(conn, run_id, value)
            except (Unresolved, RuntimeError):      # a RuntimeError is a resolver's refusal (the galaxy has no skeleton, say)
                if strict:
                    raise
                found = value
            if value.get(MERGE) and found is not value:
                out.update(found)
                continue
            value = found
        out[name] = value
    return out


FINISHED_BY = "finished_by"
"""str: The step parameter that marks a stage a run carries on from (PERF.182, `generation/finish_run.py`): the id of the earlier
run that finished it. The stage is not cut, prepared, run or concluded again: it opens and ends at once with no unit, and
carries the conclusion the earlier run stored, which the process following the run reads as it would a stage that just ended."""


def _carry_finished_stage(conn, run_id, index, step, now):
    """Opens and ends stage `index` at once as one the run before it finished (`FINISHED_BY`): a plan with no unit that says
    so, and the earlier run's conclusion (`runrecord.stored_conclusion`), if the command has a `conclude`."""
    earlier = step.params[FINISHED_BY]
    runrecord.set_stage_plan(conn, run_id, index, 0, 0.0,
                             {"kind": step.kind, "units": 0, "cost": 0.0,
                              "detail": {"skip": f"finished by run {earlier}", FINISHED_BY: earlier}})
    runrecord.begin_stage(conn, run_id, index, now=now)
    runrecord.end_stage(conn, run_id, index, now=now)
    if commands.conclude_function(step.kind) is not None:
        stored = runrecord.stored_conclusion(conn, earlier, step.key)
        runrecord.set_conclusion(conn, run_id, index, stored if stored is not None else {})


def _conclude_empty_stage(conn, run_id, index, label):
    """A stage whose cut gave no unit (a plan with no shells, say) has no last unit to conclude it, so whichever process cut it
    does: it ended as it began, its `conclude` runs, and its row is stored, as the unit that completes a stage would."""
    runrecord.end_stage(conn, run_id, index)
    try:
        took = conclude_stage(conn, run_id, index)
    except Exception as exc:
        runrecord.finish_run(conn, run_id, "failed", error=f"the end of stage {label} failed: {type(exc).__name__}: {exc}")
        raise
    finish_stage(conn, run_id, index, took)


def stage_label(conn, run_id, index):
    """The label stage `index` of the run was declared with, `None` past the last."""
    declared = json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), "stages")))
    return runrecord._s(declared[index]["label"]) if 0 <= index < len(declared) else None


def stage_name(conn, run_id, index):
    """What a message calls stage `index` (PERF.221): the label the page shows for it, not the internal number; the
    number only where the run record cannot say (a run already gone from Redis)."""
    try:
        label = stage_label(conn, run_id, index)
    except Exception:  # noqa: BLE001 -- naming the stage must never hide the failure being reported
        label = None
    return label or str(index)


def _stage_step(conn, run_id, index):
    """`(declared, step)` of stage `index`: how the run declared it and the step the run record holds for it."""
    declared = json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), "stages")))[index]
    step = next(commands.Step.from_dict(item) for item in runrecord.read_steps(conn, run_id)
                if item["key"] == declared["key"])
    return declared, commands.Step(step.key, step.kind, step.label, resolve_params(conn, run_id, step.params))


def _stage_results(conn, run_id, index):
    """Every unit's result of stage `index`, in unit order."""
    return [json.loads(runrecord._s(raw)) for _unit, raw in sorted(
        conn.hgetall(runrecord.results_key(run_id, index)).items(), key=lambda pair: int(pair[0]))]


def conclude_stage(conn, run_id, index):
    """
    What the unit that completes stage `index` does first (PERF.172): the command's `conclude`, with the step's
    parameters, the stage's plan, every unit's result and the seconds the stage took. What it returns is stored
    with the stage (`read_conclusion`), and the stage counts as concluded (`wait_concluded`). A command with no
    `conclude` records nothing. Unlike `finish_stage` a failure here is not swallowed: it is the stage's own work,
    so it fails the run.

    Returns:
        float: The seconds the hook took (0.0 for a command with none).
    """
    try:
        _declared, step = _stage_step(conn, run_id, index)
    except StopIteration:        # a run that holds no steps (it was started by hand) has no command to conclude
        return 0.0
    hook = commands.conclude_function(step.kind)
    if hook is None:
        return 0.0
    record = {runrecord._s(name): runrecord._s(value)
              for name, value in conn.hgetall(runrecord.stage_key(run_id, index)).items()}
    plan = runrecord.read_plan(conn, run_id, index)
    began = time.time()
    conclusion = hook(step.params, plan["detail"], _stage_results(conn, run_id, index),
                      max(float(record["ended"]) - float(record["started"]), 0.0)) or {}
    took = time.time() - began
    runrecord.set_conclusion(conn, run_id, index, conclusion)
    return took


def finish_stage(conn, run_id, index, concluded_seconds=0.0):
    """
    What the unit that completes stage `index` does before the run moves on (PERF.106): the command's `finish`
    records its units' speeds and gives the stage's metrics (`CommandSpec.finish`), and the stage's row
    (`runrecord.Stage.row`) is stored in `generation_stage_runs` with the seconds from the stage's start to its last
    unit's end, and those the command's `conclude` took (`concluded_seconds`). A stage whose expansion said it had
    nothing to do (`detail["skip"]`, the reason) is stored as skipped. Statistics are a nicety: a failure here is
    logged and the run goes on.
    """
    from planetgen.generation import stages
    try:
        declared, step = _stage_step(conn, run_id, index)
        record = {runrecord._s(name): runrecord._s(value)
                  for name, value in conn.hgetall(runrecord.stage_key(run_id, index)).items()}
        _log_long_unit(conn, run_id, index, step, declared, float(record.get("longest_unit") or 0.0))
        seconds = max(float(record["ended"]) - float(record["started"]), 0.0) + concluded_seconds
        hook = commands.finish_function(step.kind)
        plan = runrecord.read_plan(conn, run_id, index)
        metrics = {}
        if hook is not None:
            metrics = hook(step.params, plan["detail"], _stage_results(conn, run_id, index), seconds) or {}
        if declared.get("row"):
            skipped = (plan["detail"] or {}).get("skip")
            if skipped:
                stages.write_skipped_row(_mysql_config(step.params["database"]), declared["key"], declared["row"],
                                         skipped)
            else:
                stages.write_queued_row(_mysql_config(step.params["database"]), declared["key"], declared["row"],
                                        seconds, metrics)
    except Exception as exc:  # noqa: BLE001 -- statistics are a nicety, never a reason to fail a run
        from planetgen.util import log
        log.debug(f"Stage statistics: not recorded ({exc}).")


def _log_long_unit(conn, run_id, index, step, declared, longest):
    """
    PERF.155: writes a warning to the run's log when the longest unit of the stage that just ended ran longer than
    `queue.unit_ceiling_seconds`, naming the stage, the unit's seconds and how many units the stage was cut into. A
    command that is one unit by nature (`CommandSpec.whole`) is not held to the ceiling.
    """
    ceiling = get_settings().queue.unit_ceiling_seconds
    if commands.spec(step.kind).whole or longest <= ceiling:
        return
    runrecord.log(conn, run_id, f"stage {declared['label']}: its longest unit ran {longest:,.0f} s, over the "
                  f"{ceiling:,.0f} s unit ceiling; the stage was cut into {declared.get('units', 1):,} units",
                  level="warning", kind="unit.long", stage=index)


def _mysql_config(named):
    """The `MySQLConfig` a step's `database` parameters name (`database_args`); the password is this process's
    settings'."""
    from planetgen.db import store
    return store.MySQLConfig(host=named["host"], port=named["port"], user=named["user"], database=named["database"])


def expected_seconds(step, row):
    """
    How long earlier runs took this step's stage with the settings it runs with (`generation.stats.stage_seconds`,
    PERF.56), or `None` when none is recorded, statistics are off (`row` is `None`) or the database cannot be read.
    The stage's ETA starts from it, and its units' time limit is a multiple of it (`units.unit_timeout`).
    """
    if row is None:
        return None
    try:
        from planetgen.db import store
        from planetgen.generation import stats
        conn = store.get_control_connection(store.control_mysql_config(_mysql_config(step.params["database"])))
        try:
            return stats.stage_seconds(conn, step.key, row["settings"])
        finally:
            conn.close()
    except Exception:  # noqa: BLE001 -- a time that cannot be read is a stage with none recorded
        return None


def advance(conn, run_id, index):
    """
    Queues stage `index` of the run, and each stage after it that has no units to wait for, or finishes the run
    when the last stage is one of them. Called by the unit that completes the stage before (the barrier), and by
    the starter for the first. A run that is no longer running (cancelling, cancelled, failed) queues nothing.
    """
    if runrecord._s(conn.hget(runrecord.run_key(run_id), "state")) != "running":
        return
    count = len(json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), "stages"))))
    while index < count:
        if _expand(conn, run_id, index):
            return
        index += 1
    runrecord.finish_run(conn, run_id, "done")


def _await_plan(conn, run_id, index, ensure, clock=time.monotonic, sleep=time.sleep, on_debug=None):
    """Stage `index`'s plan once its stage has been expanded (by the unit that completed the stage before it)."""
    last_ensure = float("-inf")
    while True:
        runrecord.touch(conn, run_id)
        if on_debug is not None:
            on_debug()
        plan = runrecord.read_plan(conn, run_id, index)
        if plan is not None:
            return plan
        state = runrecord._s(conn.hget(runrecord.run_key(run_id), "state"))
        if state not in runrecord.LIVE_STATES:
            error = runrecord._s(conn.hget(runrecord.run_key(run_id), "error"))
            raise units.RunStopped(error or f"the run {state} before stage {index} was queued")
        if ensure is not None and clock() - last_ensure >= units.ENSURE_SECONDS:
            last_ensure = clock()
            ensure()
        sleep(units.POLL_SECONDS)


def _await_conclusion(conn, run_id, index, ensure, clock=time.monotonic, sleep=time.sleep, on_debug=None):
    """What stage `index`'s conclusion found (`conclude_stage`), once the unit that completed the stage has made it."""
    last_ensure = float("-inf")
    while True:
        runrecord.touch(conn, run_id)
        if on_debug is not None:
            on_debug()
        found = runrecord.conclusion(conn, run_id, index)
        if found is not None:
            return found
        state = runrecord._s(conn.hget(runrecord.run_key(run_id), "state"))
        if state not in runrecord.LIVE_STATES:
            error = runrecord._s(conn.hget(runrecord.run_key(run_id), "error"))
            raise units.RunStopped(error or f"the run {state} before stage {index} was concluded")
        if ensure is not None and clock() - last_ensure >= units.ENSURE_SECONDS:
            last_ensure = clock()
            ensure()
        sleep(units.POLL_SECONDS)


def queue_of(kind):
    """The queue the units of a run of `kind` go on (PERF.150): the maintenance queue for a `MAINTENANCE` or `CHECK`
    run, the bulk queue for a `GENERATION` run."""
    return redisqueue.BULK_QUEUE if kind == runrecord.GENERATION else redisqueue.MAINTENANCE_QUEUE


def waiting_text(holder, position):
    """What a run waiting for the machine says (PERF.153; PERF.147's bar uses the same words)."""
    return f"Waiting at position {position} for run {holder}"


def report_wait(holder, position):
    """
    `runline.start_when_free`'s `on_wait`: shows that this run waits, where its bar is read (the Generate page's
    progress file) and where its job is (the API job's `meta`, which `api_jobs.status` reports as `waiting`), and
    says so in the log. `holder` `None` clears it.
    """
    if holder is not None:
        log.normal(waiting_text(holder, position))
        progress_file.report(0, None, description=waiting_text(holder, position), force=True)
    try:
        import rq
        job = rq.get_current_job()
        if job is not None:
            job.meta["waiting"] = None if holder is None else {"run": holder, "position": position,
                                                               "text": waiting_text(holder, position)}
            job.save_meta()
    except Exception:  # noqa: BLE001 -- showing the wait is a nicety, never a reason to fail the run
        pass


def stage_ensure(conn, url, kind, step_kind, workers_count=None):
    """
    The `ensure` for one stage of a run (PERF.220): `workers.ensure_workers` with the run's worker count held to what
    the stage's workers can take in memory (`workers.stage_worker_cap`), and, on the bulk queue, the live workers over
    that asked to finish their unit and exit (`workers.shed_bulk`), so the workers an earlier stage left do not all
    scatter at once. A stage with ordinary workers is not limited.
    """
    cap = workers.stage_worker_cap(step_kind)
    count = workers_count if cap is None else min(workers_count or cap, cap)

    def ensure():
        if cap is not None and queue_of(kind) == redisqueue.BULK_QUEUE:
            workers.shed_bulk(conn, cap)
        return workers.ensure_workers(queue_of(kind), conn, url, count=count)
    return ensure


def run(title, steps, on_plan, on_result, on_progress=None, workers_count=None, who=None, database=None, argv=None,
        seed=None, listing=None, url=None, kind=runrecord.GENERATION, rows=None, on_stage_done=None, galaxy_seed=None,
        build=None, finishes=None):
    """
    Runs the queued `steps` as one run: starts the run record (the machine is busy for it), queues the first stage,
    makes sure workers are running (`workers_count` of them at most: the run's `--workers`) and follows each stage to
    its end. The run ends in Redis as `done`, or `failed` or `cancelled` if this process stops on an error or an
    interrupt first.

    Args:
        steps (list[commands.Step]): The steps to queue, in order, parameters complete (`database` in them).
        on_plan (callable): `on_plan(index, plan)` when stage `index` has been queued, before its results arrive:
            `plan` is `{"kind", "units", "cost", "detail"}`, `detail` the command's JSON summary of the stage.
        on_result (callable): `on_result(index, unit, result)`, once for each unit as it finishes.
        on_progress (callable, optional): `on_progress(index, unit, done, estimate)` for the units still working.
        listing (list[commands.Step], optional): Every step of the run, inline ones too, as the wrapper's list
            (`steps_for`); stored in the run record. Default: `steps`.
        on_stage_done (callable, optional): `on_stage_done(index, conclusion)` when every unit of stage `index` is done
            and its command's `conclude` has run (`conclusion` is what it returned, `None` for a command with none).
        galaxy_seed (str, optional): The galaxy's seed as hex and `build` the code's version key, recorded with the run
            so a later `Finish run` can tell the galaxy and the code are still the ones it was started on (PERF.182);
            `finishes` is the id of the run this one carries on from.
        rows (dict, optional): `{stage key: row}`, the `generation_stage_runs` row the unit that completes each stage
            writes (`generation.stages.hand_over_open_stage`, `stages.row_for`); a stage with none stores no row.
        kind (str): `runrecord.GENERATION` for a run that changes the galaxy (it owns the machine and shows the
            visitors' banner) or `runrecord.MAINTENANCE` for one that only keeps the site in order, such as the
            cache warm-up (PERF.127).

    Returns:
        str: The run id.

    Raises:
        redisqueue.Unavailable: No Redis server answers.
        runline.LineFull: The waiting list holds `queue.max_waiting` entries already; the message names them.
        runrecord.Busy: Another run of that kind is running. A `GENERATION` run joins the waiting list (PERF.154) and
            waits for it up to `queue.wait_minutes`, then fails with the run named; a `MAINTENANCE` or `CHECK` run does
            not wait.
        units.RunStopped: A unit failed or the run was cancelled.
    """
    url = url or redisqueue.redis_url()
    conn = redisqueue.connect(url)
    try:
        conn.ping()
    except Exception as exc:  # noqa: BLE001 -- any connection failure
        raise redisqueue.Unavailable(f"no Redis server answers at {url} (config.json's redis.url); "
                                     f"start Redis to run this: {exc}") from exc
    by_key = {step.key: step for step in steps}
    recorded = [by_key.get(step.key, step) for step in (listing or steps)]
    times = stagetimes.lookup(conn, [step.key for step in steps])
    stages = [runrecord.Stage(step.key, step.label, workers=workers_count or 1,
                              expected_seconds=(expected_seconds(step, (rows or {}).get(step.key))
                                                or times.get(step.key, {}).get("seconds")),
                              worker_rate=times.get(step.key, {}).get("worker_rate"),
                              row=(rows or {}).get(step.key)) for step in steps]
    debug = (log.areas_spec() or "all") if log.console_debug() else ""     # a debug run: the workers relay their lines
    start = lambda: runrecord.start_run(       # noqa: E731
        conn, title, stages, who=who, database=database, argv=argv, seed=seed, steps=[step.as_dict() for step in recorded],
        kind=kind, disregard_missing_data=fetch.disregarding_missing_data(), debug=debug,
        job=os.environ.get(JOB_ENV) or None, galaxy_seed=galaxy_seed, build=build, finishes=finishes)
    if kind == runrecord.GENERATION:
        run_id = runline.start_when_free(conn, start, get_settings().queue.wait_minutes * 60, on_wait=report_wait,
                                         title=title, who=who, database=database,
                                         admitted=bool(os.environ.get(JOB_ENV)))
    else:
        run_id = start()
    ensures = [stage_ensure(conn, url, kind, step.kind, workers_count) for step in steps]
    relay = log.RelayReader(conn, runrecord.debug_key(run_id)) if debug else None
    on_debug = relay.drain if relay else None
    log.trace(f"run {run_id} ({title}) started: {len(steps)} stage(s) on the {queue_of(kind)} queue, "
              f"{workers_count or 'the default number of'} worker(s), kind {kind}, database {database}"
              + (f", debug areas {debug}" if debug else ""), area="queue")
    events = runrecord.EventHandler(conn, run_id)       # what this process logs is the run's live log (PERF.99)
    logging.getLogger("planetgen").addHandler(events)
    try:
        advance(conn, run_id, 0)
        for index in range(len(steps)):
            ensure = ensures[index]
            plan = _await_plan(conn, run_id, index, ensure, on_debug=on_debug)
            log.trace(f"run {run_id}: stage {index} ({steps[index].label}) queued as {plan['units']} unit(s) of total "
                      f"cost {plan['cost']}", area="queue")
            on_plan(index, plan)
            conclusion = None
            if plan["units"]:
                units.follow_stage(
                    conn, run_id, index, plan["units"], lambda unit, result, index=index: on_result(index, unit, result),
                    None if on_progress is None else
                    (lambda unit, done, estimate, index=index: on_progress(index, unit, done, estimate)),
                    ensure, on_debug=on_debug)
            if commands.spec(steps[index].kind).conclude is not None:
                conclusion = _await_conclusion(conn, run_id, index, ensure, on_debug=on_debug)
            if on_stage_done is not None:
                on_stage_done(index, conclusion)
        runrecord.finish_run(conn, run_id, "done")      # a stage of no units has no last unit to do it
        log.trace(f"run {run_id} finished: done", area="queue")
        if relay:
            relay.finish()
    except BaseException as exc:
        stopped = isinstance(exc, KeyboardInterrupt) or (       # Ctrl+C, or a SIGTERM (`cli.generate._exit_on_sigterm`)
            isinstance(exc, SystemExit) and isinstance(exc.code, int) and exc.code >= 128)
        full = isinstance(exc, redis.exceptions.OutOfMemoryError)
        message = (f"Redis is full (used memory is over its maxmemory, {exc}); the run stopped. Raise maxmemory "
                   f"in the Redis configuration or flush the page cache, then run again.") if full else \
            (str(exc) or type(exc).__name__)
        try:
            runrecord.finish_run(conn, run_id, "cancelled" if stopped else "failed", error=message)
        except redis.exceptions.RedisError:
            pass        # a full Redis refuses this write too; the original error is the one to show
        log.trace(f"run {run_id} finished: {'cancelled' if stopped else 'failed'} ({message})", area="queue")
        if galaxy_seed and kind == runrecord.GENERATION:
            log.error(f"Run {run_id} stopped before it was done. 'planetgen finish-run' carries it on from the first stage "
                      f"that did not finish.")
        if relay:
            relay.finish()
        if full:
            raise RuntimeError(message) from None
        raise
    finally:
        logging.getLogger("planetgen").removeHandler(events)
        from planetgen.queue import runcontrol      # it imports this file
        runcontrol.write_log(conn, run_id)
    return run_id
