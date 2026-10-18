# planetgen/queue/wrapper.py

"""
The wrapper of a `planetgen galaxy` or `planetgen plan` run (PERF.118, docs/design/redis-and-the-work-queue.md 5.12):
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
process that started the run still does it, between the queued stages, until PERF.106 to PERF.108 and PERF.121 to
PERF.128 give it units.
"""

import json
import time

from planetgen.queue import commands, redisqueue, runrecord, units, workers


def database_args(mysql_config):
    """The database a unit opens, as plain JSON: its name, host, port and user. Never the password; a worker reads
    that from its own settings."""
    return {"host": mysql_config.host, "port": mysql_config.port, "user": mysql_config.user,
            "database": mysql_config.database}


def steps_for(command, args):
    """The steps of a `galaxy` or `plan` run for the parsed `args` (`commands.build_steps`): the same list for the
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
    listing = steps_for(command, args) if command in ("galaxy", "plan") else []
    found = next((step for step in listing if step.key == key), None)
    step = commands.Step(key, commands.KIND_OF_STAGE[key], found.label if found else label,
                         {**(found.params if found else stages.stage_settings(key, args)), **params})
    if found is not None:
        return step, [step if item.key == key else item for item in listing]
    later = next((index for index, item in enumerate(listing) if item.key == "luminosity"), len(listing))
    return step, listing[:later] + [step] + listing[later:]


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
    expansion = commands.expand_function(step.kind)(step.params)
    cost = sum(unit_cost for _args, unit_cost in expansion.units)
    runrecord.set_stage_plan(conn, run_id, index, len(expansion.units), cost,
                             {"kind": step.kind, "units": len(expansion.units), "cost": cost,
                              "detail": expansion.plan})
    runrecord.begin_stage(conn, run_id, index)
    units.enqueue_units(conn, run_id, index, step.kind, expansion.units)
    return len(expansion.units)


def advance(conn, run_id, index):
    """
    Queues stage `index` of the run, and each stage after it that has no units to wait for, or finishes the run
    when the last stage is one of them. Called by the unit that completes the stage before (the barrier), and by
    the starter for the first.
    """
    count = len(json.loads(runrecord._s(conn.hget(runrecord.run_key(run_id), "stages"))))
    while index < count:
        if _expand(conn, run_id, index):
            return
        index += 1
    runrecord.finish_run(conn, run_id, "done")


def _await_plan(conn, run_id, index, ensure, clock=time.monotonic, sleep=time.sleep):
    """Stage `index`'s plan once its stage has been expanded (by the unit that completed the stage before it)."""
    last_ensure = float("-inf")
    while True:
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


def run(title, steps, on_plan, on_result, on_progress=None, workers_count=None, who=None, database=None, argv=None,
        seed=None, listing=None, url=None):
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

    Returns:
        str: The run id.

    Raises:
        redisqueue.Unavailable: No Redis server answers.
        runrecord.Busy: Another run owns the machine.
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
    run_id = runrecord.start_run(
        conn, title, [runrecord.Stage(step.key, step.label, workers=workers_count or 1) for step in steps],
        who=who, database=database, argv=argv, seed=seed, steps=[step.as_dict() for step in recorded])
    ensure = lambda: workers.ensure_workers(redisqueue.BULK_QUEUE, conn, url, count=workers_count)  # noqa: E731
    try:
        advance(conn, run_id, 0)
        for index in range(len(steps)):
            plan = _await_plan(conn, run_id, index, ensure)
            on_plan(index, plan)
            if plan["units"]:
                units.follow_stage(
                    conn, run_id, index, plan["units"], lambda unit, result, index=index: on_result(index, unit, result),
                    None if on_progress is None else
                    (lambda unit, done, estimate, index=index: on_progress(index, unit, done, estimate)),
                    ensure)
        runrecord.finish_run(conn, run_id, "done")      # a stage of no units has no last unit to do it
    except BaseException as exc:
        runrecord.finish_run(conn, run_id, "cancelled" if isinstance(exc, KeyboardInterrupt) else "failed",
                             error=str(exc) or type(exc).__name__)
        raise
    return run_id
