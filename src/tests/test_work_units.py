"""PERF.97: a stage of a run on the new work queue -- its units cut from the work alone, queued as JSON on the bulk queue,
run by workers, joined by the done set and followed by the process that started the run (`planetgen/queue/units.py`,
started by the wrapper, `planetgen/queue/wrapper.py`, PERF.118).
The autouse `_run_queue_inline` fixture gives each test a queue and key prefix of its own and runs the waiting units
in this process when `ensure_workers` is called."""

import pytest

from planetgen.queue import commands, redisqueue, runrecord, units, wrapper


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


def _handlers(monkeypatch, **named):
    """Registers a test command for each `kind=handler`: its units run `handler(args, report)`, and its stage is cut
    into `params["count"]` units of cost `1 + n`."""
    for kind in named:
        monkeypatch.setitem(commands.REGISTRY, kind, commands.CommandSpec(kind, "A test", unit=(__name__, "_square"),
                                                                          expand=(__name__, "_cut")))
    monkeypatch.setattr(commands, "unit_function", lambda kind: named[kind])
    monkeypatch.setattr(commands, "expand_function", lambda kind: _cut)


def _square(args, report):
    report(1, 2)
    return {"value": args["n"] ** 2}


def _cut(params):
    return commands.Expansion([({"n": n}, 1.0 + n) for n in range(params["count"])], {"count": params["count"]})


def _step(count, key="test", kind="test"):
    return commands.Step(key, kind, "A test stage", {"count": count})


def _run(monkeypatch, handler, count=5, on_progress=None, **kwargs):
    _handlers(monkeypatch, test=handler)
    results = {}
    run_id = wrapper.run("Test", [_step(count)], lambda index, plan: None,
                         lambda index, unit, result: results.__setitem__(unit, result),
                         None if on_progress is None else lambda index, unit, done, estimate: on_progress(unit, done, estimate),
                         **kwargs)
    return run_id, results


def test_the_cut_depends_on_the_work_not_on_the_workers():
    rows = [(shell, 0, 1.0 + (shell % 7)) for shell in range(1000)]
    chunks = units.static_chunks(rows, units=100)
    assert chunks == units.static_chunks(list(rows), units=100)
    assert [row for chunk in chunks for row in chunk] == rows                  # in order, nothing lost or repeated
    assert 90 <= len(chunks) <= 115
    weights = [sum(row[-1] for row in chunk) for chunk in chunks]
    assert max(weights) < 1.3 * (sum(weights) / len(weights)) + max(row[-1] for row in rows)


def test_a_unit_holds_at_most_max_rows_and_a_stage_at_most_max_units(monkeypatch):
    rows = [(shell, 0, 0.001) for shell in range(300)]
    assert max(len(chunk) for chunk in units.static_chunks(rows, units=1, max_rows=64)) == 64
    monkeypatch.setattr(units, "MAX_UNITS", 10)
    capped = units.static_chunks(rows, units=300, max_rows=1)
    assert len(capped) <= 10 and [row for chunk in capped for row in chunk] == rows
    assert units.static_chunks([]) == []


def test_every_unit_runs_once_and_the_run_ends_done(monkeypatch, conn):
    seen = []
    progress = []

    def handler(args, report):
        seen.append(args["n"])
        return _square(args, report)

    run_id, results = _run(monkeypatch, handler, on_progress=lambda unit, done, estimate: progress.append(unit))
    assert sorted(seen) == [0, 1, 2, 3, 4]
    assert results == {n: {"value": n * n} for n in range(5)}
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "done" and run["stages"][0]["done"] == 5 and run["stages"][0]["units"] == 5
    assert run["stages"][0]["cost_done"] == pytest.approx(sum(1.0 + n for n in range(5)))
    assert runrecord.active_run(conn) is None, "the machine is free again"
    assert not conn.hgetall(runrecord.progress_key(run_id, 0)), "a unit that finished reports nothing more"


def test_the_arguments_travel_as_json_and_carry_no_credentials(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=1)])
    units.enqueue_units(conn, run_id, 0, "test", [({"n": 3, "database": {"host": "h", "user": "u"}}, 1.0)])
    queue = redisqueue.queue(redisqueue.BULK_QUEUE, conn)
    job = queue.jobs[0]
    assert job.func_name == units.UNIT_FUNCTION and job.id == units.job_id(run_id, 0, 0)
    assert tuple(job.args[:4]) == (run_id, 0, 0, "test") and job.args[4]["n"] == 3
    assert job.serializer is redisqueue.JsonSerializer and "password" not in job.args[4].get("database", {})


def test_a_unit_found_done_has_its_result_and_a_repeated_unit_counts_once(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=2)])
    runrecord.begin_stage(conn, run_id, 0)
    units.enqueue_units(conn, run_id, 0, "test", [({"n": 1}, 1.0), ({"n": 2}, 1.0)])
    units.enqueue_units(conn, run_id, 0, "test", [({"n": 1}, 1.0)])         # unit 0 queued twice (a retry)
    from planetgen.queue import workers
    workers.ensure_workers(redisqueue.BULK_QUEUE, conn)
    got = {}
    units.follow_stage(conn, run_id, 0, 2, lambda unit, result: got.__setitem__(unit, result))
    assert got == {0: {"value": 1}, 1: {"value": 4}}
    assert runrecord.read_run(conn, run_id)["stages"][0]["done"] == 2


def test_a_failing_unit_fails_the_run_and_the_starter_says_so(monkeypatch, conn):
    def handler(args, report):
        if args["n"] == 2:
            raise ValueError("no stars here")
        return _square(args, report)

    _handlers(monkeypatch, test=handler)
    stopped = []
    with pytest.raises(units.RunStopped, match="unit 2 of stage 0 failed: ValueError: no stars here"):
        wrapper.run("T", [_step(4)], lambda index, plan: None, lambda index, unit, result: stopped.append(unit))
    run_id = runrecord.recent_runs(conn)[0]
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "failed" and "no stars here" in run["error"]
    assert runrecord.active_run(conn) is None


def test_the_units_of_a_run_that_is_no_longer_live_do_nothing(monkeypatch, conn):
    ran = []
    _handlers(monkeypatch, test=lambda args, report: ran.append(args["n"]) or {})
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=3)])
    units.enqueue_units(conn, run_id, 0, "test", [({"n": n}, 1.0) for n in range(3)])
    runrecord.request_cancel(conn, run_id)
    from planetgen.queue import workers
    workers.ensure_workers(redisqueue.BULK_QUEUE, conn)
    assert ran == [] and runrecord.read_run(conn, run_id)["stages"][0]["done"] == 0


def test_an_interrupt_in_the_starter_cancels_the_run_and_frees_the_machine(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)

    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(wrapper.workers, "ensure_workers", interrupted)      # Ctrl+C while the starter waits
    with pytest.raises(KeyboardInterrupt):
        wrapper.run("T", [_step(2)], lambda index, plan: None, lambda index, unit, result: None)
    run = runrecord.read_run(conn, runrecord.recent_runs(conn)[0])
    assert run["state"] == "cancelled" and runrecord.active_run(conn) is None


def test_a_second_run_is_refused_while_one_has_the_machine(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)
    first = runrecord.start_run(conn, "First", [runrecord.Stage("test", "T", units=1)])
    with pytest.raises(runrecord.Busy):
        wrapper.run("Second", [_step(1)], lambda index, plan: None, lambda index, unit, result: None)
    assert runrecord.active_run(conn) == first


def test_a_stage_of_no_units_ends_at_once(monkeypatch, conn):
    run_id, results = _run(monkeypatch, _square, count=0)
    assert results == {} and runrecord.read_run(conn, run_id)["state"] == "done"


def test_more_units_than_the_cap_are_refused(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)
    monkeypatch.setattr(units, "MAX_UNITS", 2)
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=3)])
    with pytest.raises(ValueError, match="at most 2 units"):
        units.enqueue_units(conn, run_id, 0, "test", [({}, 1.0)] * 3)
    with pytest.raises(ValueError, match="no unit kind"):
        units.enqueue_units(conn, run_id, 0, "nope", [({}, 1.0)])


def test_only_the_table_names_what_a_unit_runs():
    assert callable(commands.unit_function("scatter-mass")) and callable(commands.expand_function("scatter-mass"))
    for kind in ("os.system", "nope"):
        with pytest.raises(KeyError):
            commands.unit_function(kind)
    with pytest.raises(ValueError, match="inline"):
        commands.unit_function("scatter-luminosity")


# --- The massive-star scatter on the new queue (the stage PERF.97 moves, the wrapper queues, PERF.118) ----------

def _scatter_mass(mysql_config, workers, on_run_record):
    """Draws the toy galaxy's massive stars (8 solar masses and up) and returns the stars stored, sorted."""
    from planetgen.db import store
    from planetgen.galaxy import scatter_limits
    from planetgen.generation import run_plan
    from planetgen.generation import star_population as starPopulation
    from tests.test_bright_star_scatter import EXTENTS, _plan_args, _seed_galaxy
    _seed_galaxy(mysql_config)
    args = _plan_args(mysql_config, "--bright-stars-only", "--workers", str(workers))
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        conn.execute("DELETE FROM bright_stars")
        conn.commit()
    finally:
        conn.close()
    if on_run_record:
        run_plan._scatter_mass(args, mysql_config, 1000.0, 8.0, 77)
    else:
        limits = scatter_limits.ScatterLimits.of(skeleton, 1000.0, 8.0)
        tiered = scatter_limits.mass_pass(limits, starPopulation.MASS_PASS_MIN_LUMINOSITY_SOL)
        run_plan._scatter_shells(args, mysql_config, skeleton, EXTENTS, set(),
                                 starPopulation.MASS_PASS_MIN_LUMINOSITY_SOL, 77, "Massive stars",
                                 stats_kind="scatter-mass", task_kind="mass-stars", tiered=tiered)
    conn = store.get_connection(mysql_config)
    try:
        columns = ", ".join(store.BRIGHT_STAR_COLUMNS)
        return sorted(tuple(row.values()) for row in conn.execute(f"SELECT {columns} FROM bright_stars").fetchall())
    finally:
        conn.close()


@pytest.mark.fresh_database
def test_the_mass_scatter_gives_the_same_stars_at_any_worker_count_and_as_the_old_queue(mysql_config):
    one = _scatter_mass(mysql_config, 1, True)
    assert len(one) > 5, "the toy galaxy holds massive stars to compare"
    assert _scatter_mass(mysql_config, 4, True) == one
    assert _scatter_mass(mysql_config, 1, False) == one


@pytest.mark.fresh_database
def test_the_mass_scatter_records_its_run_with_a_unit_per_run_of_shells(mysql_config, monkeypatch):
    runs = []
    real = wrapper.run
    monkeypatch.setattr(wrapper, "run", lambda *args, **kwargs: runs.append(real(*args, **kwargs)) or runs[-1])
    _scatter_mass(mysql_config, 2, True)
    conn = redisqueue.connect()
    run = runrecord.read_run(conn, runs[0])
    assert run["state"] == "done" and run["database"] == mysql_config.database
    assert [stage["key"] for stage in run["stages"]] == ["mass"]
    assert run["stages"][0]["done"] == run["stages"][0]["units"] > 0
    plan = runrecord.read_plan(conn, runs[0], 0)
    assert plan["kind"] == "scatter-mass" and plan["units"] == run["stages"][0]["units"]
    assert [step["kind"] for step in runrecord.read_steps(conn, runs[0])] == ["scatter-mass"]


@pytest.mark.fresh_database
def test_a_unit_of_the_mass_scatter_runs_alone_from_its_json_arguments(mysql_config):
    import json

    from planetgen.db import store
    from planetgen.generation import run_plan
    from tests.test_bright_star_scatter import _seed_galaxy
    _seed_galaxy(mysql_config)
    params = {"database": wrapper.database_args(mysql_config), "min_luminosity_sol": 1000.0, "mass_limit_sol": 8.0,
              "seed": 77}
    expansion = run_plan.expand_scatter_mass(params)
    again = run_plan.expand_scatter_mass(json.loads(json.dumps(params)))
    assert expansion.units == again.units and expansion.plan == again.plan, "the same params give the same units"
    assert len(expansion.units) == len(expansion.plan["tasks"]) > 1
    for args, cost in expansion.units:
        assert json.loads(json.dumps(args)) == args and "password" not in args["database"]
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM bright_stars")
        conn.commit()
        stars = 0
        for args, _cost in expansion.units:
            result = commands.run_alone("scatter-mass", json.loads(json.dumps(args)))
            stars += sum(result["counts"].values())
        assert stars > 5
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == stars
    finally:
        conn.close()


def test_the_shells_are_cut_into_the_same_units_whatever_the_worker_count():
    from planetgen.generation import run_plan
    rows = [(shell, 3, 1.0 + shell % 5) for shell in range(300)]
    chunks = run_plan.static_shell_chunks(rows)
    assert chunks == run_plan.static_shell_chunks(list(rows))
    assert sorted(row for chunk in chunks for row in chunk) == rows
    assert max(len(chunk) for chunk in chunks) <= run_plan.MAX_SHELLS_PER_CHUNK
    weights = [sum(row[2] for row in chunk) for chunk in chunks]
    assert weights == sorted(weights, reverse=True), "the heaviest unit first"


def test_a_tiered_pass_travels_as_json_and_comes_back_equal():
    from planetgen.galaxy import scatter_limits
    from tests.test_gen_scatter_tiers import FLOOR, MASS, _skeleton
    skeleton = _skeleton()
    limits = scatter_limits.ScatterLimits.of(skeleton, FLOOR, MASS)
    for tiered in (scatter_limits.mass_pass(limits, 100.0), scatter_limits.luminosity_pass(limits),
                   scatter_limits.band_pass(limits, scatter_limits.ScatterLimits.of(skeleton, FLOOR * 10, MASS)),
                   scatter_limits.mass_pass(scatter_limits.ScatterLimits.of(skeleton, FLOOR, None), 100.0), None):
        import json
        sent = json.loads(json.dumps(scatter_limits.tiered_to_json(tiered)))
        assert scatter_limits.tiered_from_json(sent, skeleton) == tiered
