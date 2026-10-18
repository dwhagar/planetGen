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


def test_a_units_time_limit_is_a_multiple_of_its_stages_expected_time_and_never_none(monkeypatch, conn):
    assert units.unit_timeout(None) == units.unit_timeout(0) == units.UNIT_TIMEOUT_SECONDS
    assert units.unit_timeout(1000) == 10_000 and units.unit_timeout(1) == units.MIN_UNIT_TIMEOUT_SECONDS
    assert units.unit_timeout(10 ** 9) == units.UNIT_TIMEOUT_SECONDS
    _handlers(monkeypatch, test=_square)
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("a", "A", units=1), runrecord.Stage("b", "B", units=1)])
    units.enqueue_units(conn, run_id, 0, "test", [({"n": 1}, 1.0)], expected_seconds=2000)
    units.enqueue_units(conn, run_id, 1, "test", [({"n": 2}, 1.0)])
    jobs = {job.id: job for job in redisqueue.queue(redisqueue.BULK_QUEUE, conn).jobs}
    assert {job.timeout for job in jobs.values()} == {20_000, units.UNIT_TIMEOUT_SECONDS}


def test_a_stage_with_a_recorded_time_gets_it_as_its_expected_seconds_and_its_units_a_limit_from_it(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)
    monkeypatch.setattr(wrapper, "expected_seconds", lambda step, row: 3000.0 if row else None)
    row = {"command": "plan", "n": 1, "total": 1, "label": "A test stage", "settings": {"workers": 1}}
    seen, real = [], units.enqueue_units
    monkeypatch.setattr(units, "enqueue_units", lambda *args, **kwargs: seen.append(kwargs) or real(*args, **kwargs))
    run_id = wrapper.run("T", [_step(1)], lambda index, plan: None, lambda index, unit, result: None, rows={"test": row})
    stage = runrecord.read_run(conn, run_id)["stages"][0]
    assert stage["expected_seconds"] == 3000.0 and stage["row"] == row
    assert seen == [{"expected_seconds": 3000.0}]


def test_a_killed_worker_process_fails_its_run_at_once_naming_the_unit(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=3)])
    units.enqueue_units(conn, run_id, 0, "test", [({"n": n}, 1.0) for n in range(3)])
    job = redisqueue.queue(redisqueue.BULK_QUEUE, conn).jobs[1]
    units.horse_killed(job, 4242, -9, None)
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "failed" and "unit 1 of stage 0 was killed" in run["error"] and "4242" in run["error"]
    assert runrecord.active_run(conn) is None, "the machine is free at once"
    from planetgen.queue import workers
    workers.ensure_workers(redisqueue.BULK_QUEUE, conn)             # RQ's retry finds the run over and does nothing
    assert runrecord.read_run(conn, run_id)["stages"][0]["done"] == 0


def test_a_killed_job_that_is_not_a_unit_is_left_to_rq(conn):
    other = type("Job", (), {"func_name": "planetgen.somewhere.else", "args": ("x", 1, 2), "connection": conn})()
    units.horse_killed(other, 1, -9, None)               # no run to fail, nothing raised


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
    for kind in ("scatter-luminosity", "scatter-band", "scatter-phenomena"):
        assert callable(commands.unit_function(kind)) and callable(commands.expand_function(kind))
    with pytest.raises(ValueError, match="inline"):
        commands.unit_function("scatter-backfill")


# --- The massive-star scatter on the new queue (the stage PERF.97 moves, the wrapper queues, PERF.118) ----------

def _stored_stars(mysql_config):
    from planetgen.db import store
    conn = store.get_connection(mysql_config)
    try:
        columns = ", ".join(store.BRIGHT_STAR_COLUMNS)
        return sorted(tuple(row.values()) for row in conn.execute(f"SELECT {columns} FROM bright_stars").fetchall())
    finally:
        conn.close()


def _scatter_mass(mysql_config, workers):
    """Draws the toy galaxy's massive stars (8 solar masses and up) and returns the stars stored, sorted."""
    from planetgen.db import store
    from planetgen.generation import run_plan
    from tests.test_bright_star_scatter import _plan_args, _seed_galaxy
    _seed_galaxy(mysql_config)
    args = _plan_args(mysql_config, "--bright-stars-only", "--workers", str(workers))
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM bright_stars")
        conn.commit()
    finally:
        conn.close()
    run_plan._scatter_mass(args, mysql_config, 1000.0, 8.0, 77)
    return _stored_stars(mysql_config)


@pytest.mark.fresh_database
def test_the_mass_scatter_gives_the_same_stars_at_any_worker_count(mysql_config):
    one = _scatter_mass(mysql_config, 1)
    assert len(one) > 5, "the toy galaxy holds massive stars to compare"
    assert _scatter_mass(mysql_config, 4) == one


@pytest.mark.fresh_database
def test_the_mass_scatter_records_its_run_with_a_unit_per_run_of_shells(mysql_config, monkeypatch):
    runs = []
    real = wrapper.run
    monkeypatch.setattr(wrapper, "run", lambda *args, **kwargs: runs.append(real(*args, **kwargs)) or runs[-1])
    _scatter_mass(mysql_config, 2)
    conn = redisqueue.connect()
    run = runrecord.read_run(conn, runs[0])
    assert run["state"] == "done" and run["database"] == mysql_config.database
    assert [stage["key"] for stage in run["stages"]] == ["mass"]
    assert run["stages"][0]["done"] == run["stages"][0]["units"] > 0
    plan = runrecord.read_plan(conn, runs[0], 0)
    assert plan["kind"] == "scatter-mass" and plan["units"] == run["stages"][0]["units"]
    assert [step["kind"] for step in runrecord.read_steps(conn, runs[0])] == ["scatter-mass"]


def _scatter_then_band(mysql_config, workers):
    """Scatters the toy galaxy's bright stars (both passes), then adds the band down to 100 L_sun; returns the stars
    stored after each, sorted."""
    from planetgen.generation import run_plan
    from tests.test_bright_star_scatter import _plan_args, _seed_galaxy
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--workers", str(workers)))
    scattered = _stored_stars(mysql_config)
    run_plan.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", "100", "--workers", str(workers)))
    return scattered, _stored_stars(mysql_config)


@pytest.mark.fresh_database
def test_the_luminosity_pass_and_the_band_give_the_same_stars_at_any_worker_count_and_record_their_runs(
        mysql_config, monkeypatch):
    runs = []
    real = wrapper.run
    monkeypatch.setattr(wrapper, "run", lambda *args, **kwargs: runs.append(real(*args, **kwargs)) or runs[-1])
    scattered, banded = _scatter_then_band(mysql_config, 1)
    assert len(scattered) > 20 and len(banded) > len(scattered)
    conn = redisqueue.connect()
    assert [step["kind"] for run in runs for step in runrecord.read_steps(conn, run)] == [
        "scatter-mass", "scatter-luminosity", "scatter-band"]
    for run_id in runs:
        run = runrecord.read_run(conn, run_id)
        assert run["state"] == "done" and run["stages"][0]["done"] == run["stages"][0]["units"] > 0
    monkeypatch.setattr(wrapper, "run", real)
    assert _scatter_then_band(mysql_config, 4) == (scattered, banded)


@pytest.mark.fresh_database
def test_a_scatter_stage_is_stored_once_by_the_unit_that_completes_it_with_its_speeds(mysql_config, monkeypatch):
    from planetgen.db import store
    from planetgen.generation import run_plan, stages
    from planetgen.generation import stats as generation_stats
    from tests.test_bright_star_scatter import _plan_args, _seed_galaxy
    store.get_control_connection(mysql_config, ensure_schema=True).close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    _seed_galaxy(mysql_config)
    args = _plan_args(mysql_config, "--bright-stars-only", "--workers", "2")
    stages.begin(stages.for_command("plan", args), args, "plan")
    summary = run_plan.scatter_bright_stars(args)
    stages.finish()
    conn = store.get_control_connection(mysql_config)
    try:
        history = generation_stats.stage_history(conn)
        job = generation_stats.stage_history(conn, stages.JOB_KEY)
        kinds = {kind for kind, _workers, _bucket in generation_stats.GenerationStats(
            store.control_mysql_config(mysql_config)).buckets}
    finally:
        conn.close()
    by_key = {}
    for run in history:
        by_key.setdefault(run["stage_key"], []).append(run)
    assert {key: len(rows) for key, rows in by_key.items() if key in ("mass", "luminosity")} == {"mass": 1, "luminosity": 1}
    mass, light = by_key["mass"][0], by_key["luminosity"][0]
    assert mass["settings"]["workers"] == 2 and mass["seconds"] > 0 and not mass["skipped"]
    assert (mass["stage_n"], mass["stage_total"]) == (1, 2) and (light["stage_n"], light["stage_total"]) == (2, 2)
    assert mass["metrics"]["objects"] + light["metrics"]["objects"] == summary["total"] > 20
    assert mass["metrics"]["shells"] > 0 and 0 < mass["metrics"]["shells_modified"] <= mass["metrics"]["shells"]
    assert mass["metrics"]["seconds_per_shell"] == pytest.approx(mass["seconds"] / mass["metrics"]["shells"])
    assert {"scatter-mass", "scatter-bright"} <= kinds, "each unit's speed is recorded under its pass's kind"
    assert job and job[0]["metrics"]["shells"] > 0, "the whole job's totals are still the starter's"


@pytest.mark.fresh_database
def test_the_band_stage_is_stored_once_by_the_starter_because_it_tops_up_after_its_units(mysql_config, monkeypatch):
    from planetgen.db import store
    from planetgen.generation import run_plan, stages
    from planetgen.generation import stats as generation_stats
    from tests.test_bright_star_scatter import _plan_args, _seed_galaxy
    store.get_control_connection(mysql_config, ensure_schema=True).close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    run_plan.run_plan(_plan_args(mysql_config, "--bright-stars-down-to", "100", "--workers", "2"))
    conn = store.get_control_connection(mysql_config)
    try:
        bands = [run for run in generation_stats.stage_history(conn) if run["stage_key"] == "band"]
        kinds = {kind for kind, _workers, _bucket in generation_stats.GenerationStats(
            store.control_mysql_config(mysql_config)).buckets}
    finally:
        conn.close()
    assert len(bands) == 1 and bands[0]["seconds"] > 0 and bands[0]["settings"]["workers"] == 2
    assert "scatter-bright" in kinds, "the band's units still record their speeds"


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


def _scatter_phenomena(mysql_config, workers):
    """Scatters the toy galaxy's phenomena and returns the rows stored (no id), sorted."""
    import json

    from planetgen.db import store
    from planetgen.generation import run_plan
    from tests.test_bright_star_scatter import _plan_args, _seed_galaxy
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config, "--workers", str(workers)))
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT * FROM phenomenon_scatter").fetchall()
    finally:
        conn.close()
    return sorted(json.dumps({key: str(value) for key, value in row.items() if key != "id"}, sort_keys=True)
                  for row in rows)


@pytest.mark.fresh_database
def test_the_phenomena_scatter_gives_the_same_rows_at_any_worker_count_and_records_its_run(mysql_config, monkeypatch):
    runs = []
    real = wrapper.run
    monkeypatch.setattr(wrapper, "run", lambda *args, **kwargs: runs.append(real(*args, **kwargs)) or runs[-1])
    one = _scatter_phenomena(mysql_config, 1)
    assert len(one) > 50
    conn = redisqueue.connect()
    assert [step["kind"] for run in runs for step in runrecord.read_steps(conn, run)] == ["scatter-phenomena"]
    run = runrecord.read_run(conn, runs[0])
    assert run["state"] == "done" and run["stages"][0]["done"] == run["stages"][0]["units"] > 1
    monkeypatch.setattr(wrapper, "run", real)
    assert _scatter_phenomena(mysql_config, 4) == one


@pytest.mark.fresh_database
def test_a_unit_of_the_phenomena_scatter_runs_alone_from_its_json_arguments(mysql_config):
    import json

    from planetgen.db import store
    from planetgen.generation import run_plan
    from tests.test_bright_star_scatter import _seed_galaxy
    _seed_galaxy(mysql_config)
    params = {"database": wrapper.database_args(mysql_config), "min_mass_solar": 1.0, "seed": 77}
    expansion = run_plan.expand_scatter_phenomena(params)
    assert expansion.units == run_plan.expand_scatter_phenomena(json.loads(json.dumps(params))).units
    assert expansion.plan == json.loads(json.dumps(expansion.plan))
    shape = run_plan._PhenomenaShape.from_json(expansion.plan)
    assert shape.to_json() == expansion.plan and len(shape.tasks) == len(expansion.units) > 1
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM phenomenon_scatter")
        conn.commit()
        placed = 0
        for args, _cost in expansion.units:
            assert json.loads(json.dumps(args)) == args and "password" not in args["database"]
            result = commands.run_alone("scatter-phenomena", json.loads(json.dumps(args)))
            placed += sum(result["counts"].values())
        assert placed > 50
        assert conn.execute("SELECT COUNT(*) AS n FROM phenomenon_scatter").fetchone()["n"] == placed
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
