"""PERF.97: a stage of a run on the new work queue -- its units cut from the work alone, queued as JSON on the bulk queue,
run by workers, joined by the done set and followed by the process that started the run (`planetgen/queue/units.py`,
started by the wrapper, `planetgen/queue/wrapper.py`, PERF.118).
The autouse `_run_queue_inline` fixture gives each test a queue and key prefix of its own and runs the waiting units
in this process when `ensure_workers` is called."""

import threading

import pytest

from planetgen.queue import commands, redisqueue, runline, runrecord, units, wrapper
from planetgen.util.settings import get_settings


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


def test_a_fixed_piece_size_grows_by_whole_multiples_to_stay_within_the_cap(monkeypatch):
    """Boss, 2026-10-24: "unit 0 of stage 3 failed: a stage queues at most 4000 units, not 4051"."""
    assert units.coarsen(0, 128) == 128 and units.coarsen(1, 128) == 128
    assert units.coarsen(128 * units.MAX_UNITS, 128) == 128, "exactly the cap needs no growth"
    assert units.coarsen(128 * units.MAX_UNITS + 1, 128) == 256, "one piece over the cap doubles the piece"
    assert units.coarsen(128 * 4051, 128) == 256
    monkeypatch.setattr(units, "MAX_UNITS", 10)
    for count in (1, 99, 1280, 1281, 5000, 123457):
        size = units.coarsen(count, 128)
        assert size % 128 == 0 and -(-count // size) <= 10, (count, size)
        assert size == 128 or -(-count // (size - 128)) > 10, "no bigger than it has to be"


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
    with pytest.raises(units.RunStopped, match="unit 2 of stage A test stage failed: ValueError: no stars here"):
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


def test_a_second_run_is_refused_while_one_has_the_machine_when_it_may_not_wait(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)
    monkeypatch.setattr(get_settings().queue, "wait_minutes", 0)
    first = runrecord.start_run(conn, "First", [runrecord.Stage("test", "T", units=1)])
    with pytest.raises(runrecord.Busy):
        wrapper.run("Second", [_step(1)], lambda index, plan: None, lambda index, unit, result: None)
    assert runrecord.active_run(conn) == first


def test_a_second_run_waits_in_line_for_the_first_and_then_runs(monkeypatch, conn):
    """PERF.153: a generation run that finds the machine taken waits for it instead of failing Busy."""
    _handlers(monkeypatch, test=_square)
    first = runrecord.start_run(conn, "First", [runrecord.Stage("test", "T", units=1)])
    waits = []
    monkeypatch.setattr(wrapper, "report_wait", lambda holder, position: waits.append((holder, position)))
    timer = threading.Timer(1.5, lambda: runrecord.finish_run(conn, first, "done"))
    timer.start()
    try:
        results = {}
        run_id = wrapper.run("Second", [_step(2)], lambda index, plan: None,
                             lambda index, unit, result: results.__setitem__(unit, result))
    finally:
        timer.join()
    assert runrecord.read_run(conn, run_id)["state"] == "done" and len(results) == 2
    assert waits[0] == (first, 1) and waits[-1] == (None, None)


def test_a_runs_results_do_not_depend_on_the_place_it_started_from(monkeypatch, conn):
    """PERF.154: the waiting list only decides when a run starts. The same steps give the same units, the same results and
    the same stages whether the machine was free, or the run was next in line, or third (the generated sectors
    are a function of the seed and the cut of the work, never of the place in the line)."""
    _handlers(monkeypatch, test=_square)
    monkeypatch.setattr(wrapper, "report_wait", lambda holder, position: None)

    def run_from(ahead):
        """A run that finds the machine held and `ahead` other runs in front of it."""
        holder = runrecord.start_run(conn, "Holder", [runrecord.Stage("test", "T", units=1)]) if ahead is not None else None
        tickets = [runline.join(conn, runline.RUN, f"Ahead {number}") for number in range(ahead or 0)]
        timers = [threading.Timer(1.0 + 0.4 * number, runline.leave, (conn, ticket)) for number, ticket in enumerate(tickets)]
        if holder:
            timers.append(threading.Timer(1.2 + 0.4 * len(tickets), runrecord.finish_run, (conn, holder, "done")))
        for timer in timers:
            timer.start()
        try:
            results = {}
            run_id = wrapper.run("Same", [_step(4)], lambda index, plan: None,
                                 lambda index, unit, result: results.__setitem__(unit, result))
        finally:
            for timer in timers:
                timer.join()
        run = runrecord.read_run(conn, run_id)
        return results, [(stage["key"], stage["units"], stage["cost"]) for stage in run["stages"]]

    free = run_from(None)
    next_in_line = run_from(0)
    third_in_line = run_from(2)
    assert free[0] == {unit: {"value": unit ** 2} for unit in range(4)}
    assert next_in_line == free and third_in_line == free


def test_a_run_gives_up_after_the_wait_bound_and_names_the_run_it_waited_for(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)
    monkeypatch.setattr(get_settings().queue, "wait_minutes", 0.03)           # 1.8 seconds
    first = runrecord.start_run(conn, "First", [runrecord.Stage("test", "T", units=1)])
    with pytest.raises(runrecord.Busy) as refused:
        wrapper.run("Second", [_step(1)], lambda index, plan: None, lambda index, unit, result: None)
    assert first in str(refused.value) and "waited" in str(refused.value)
    assert runrecord.active_run(conn) == first


def test_a_maintenance_run_still_gets_busy_and_does_not_wait(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)
    first = runrecord.start_run(conn, "Warm", [runrecord.Stage("test", "T", units=1)], kind=runrecord.MAINTENANCE)
    with pytest.raises(runrecord.Busy):
        wrapper.run("Warm again", [_step(1)], lambda index, plan: None, lambda index, unit, result: None,
                    kind=runrecord.MAINTENANCE)
    assert runrecord.active_run(conn, runrecord.MAINTENANCE) == first


def test_a_stage_of_no_units_ends_at_once(monkeypatch, conn):
    run_id, results = _run(monkeypatch, _square, count=0)
    assert results == {} and runrecord.read_run(conn, run_id)["state"] == "done"


def test_a_units_time_limit_is_a_multiple_of_its_stages_expected_time_and_never_none(monkeypatch, conn):
    assert units.unit_timeout(None) == units.unit_timeout(0) == redisqueue.MAX_TIMEOUT_SECONDS
    assert units.unit_timeout(1000) == 10_000 and units.unit_timeout(1) == redisqueue.MIN_TIMEOUT_SECONDS
    assert units.unit_timeout(10 ** 9) == redisqueue.MAX_TIMEOUT_SECONDS
    _handlers(monkeypatch, test=_square)
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("a", "A", units=1), runrecord.Stage("b", "B", units=1)])
    units.enqueue_units(conn, run_id, 0, "test", [({"n": 1}, 1.0)], expected_seconds=2000)
    units.enqueue_units(conn, run_id, 1, "test", [({"n": 2}, 1.0)])
    jobs = {job.id: job for job in redisqueue.queue(redisqueue.BULK_QUEUE, conn).jobs}
    assert {job.timeout for job in jobs.values()} == {20_000, redisqueue.MAX_TIMEOUT_SECONDS}


def test_a_stage_with_a_recorded_time_gets_it_as_its_expected_seconds_and_its_units_a_limit_from_it(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)
    monkeypatch.setattr(wrapper, "expected_seconds", lambda step, row: 3000.0 if row else None)
    row = {"command": "plan", "n": 1, "total": 1, "label": "A test stage", "settings": {"workers": 1}}
    seen, real = [], units.enqueue_units
    monkeypatch.setattr(units, "enqueue_units", lambda *args, **kwargs: seen.append(kwargs) or real(*args, **kwargs))
    run_id = wrapper.run("T", [_step(1)], lambda index, plan: None, lambda index, unit, result: None, rows={"test": row})
    stage = runrecord.read_run(conn, run_id)["stages"][0]
    assert stage["expected_seconds"] == 3000.0 and stage["row"] == row
    assert seen == [{"queue_name": redisqueue.BULK_QUEUE, "expected_seconds": 3000.0}]


def test_a_killed_worker_process_fails_its_run_at_once_naming_the_unit(monkeypatch, conn):
    _handlers(monkeypatch, test=_square)
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=3)])
    units.enqueue_units(conn, run_id, 0, "test", [({"n": n}, 1.0) for n in range(3)])
    job = redisqueue.queue(redisqueue.BULK_QUEUE, conn).jobs[1]
    units.horse_killed(job, 4242, -9, None)
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "failed" and "unit 1 of stage T was killed" in run["error"] and "4242" in run["error"]
    assert runrecord.active_run(conn) is None, "the machine is free at once"
    from planetgen.queue import workers
    workers.ensure_workers(redisqueue.BULK_QUEUE, conn)             # RQ's retry finds the run over and does nothing
    assert runrecord.read_run(conn, run_id)["stages"][0]["done"] == 0


def test_a_unit_that_rq_fails_after_it_was_done_does_not_fail_the_run(monkeypatch, conn):
    """The Execution-not-found crash of 2026-10-23: RQ failed units whose sectors were saved and counted, when it cleaned
    up after them, and each failure ended the whole run."""
    _handlers(monkeypatch, test=_square)
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=3)])
    units.enqueue_units(conn, run_id, 0, "test", [({"n": n}, 1.0) for n in range(3)])
    from planetgen.queue import workers
    workers.ensure_workers(redisqueue.BULK_QUEUE, conn)
    assert runrecord.read_run(conn, run_id)["stages"][0]["done"] == 3
    job = type("Job", (), {"func_name": units.UNIT_FUNCTION, "args": (run_id, 0, 1, "test", {"n": 1}, 1.0)})()
    units.unit_failed(job, conn, ValueError, ValueError("Execution 9ed8 not found in Redis"), None)
    assert runrecord.read_run(conn, run_id)["state"] == "done", "the unit was done: the run is not failed"
    open_run = runrecord.start_run(conn, "T2", [runrecord.Stage("test", "T", units=3)])
    units.enqueue_units(conn, open_run, 0, "test", [({"n": n}, 1.0) for n in range(3)])
    job = type("Job", (), {"func_name": units.UNIT_FUNCTION, "args": (open_run, 0, 1, "test", {"n": 1}, 1.0)})()
    units.unit_failed(job, conn, ValueError, ValueError("no stars here"), None)
    assert runrecord.read_run(conn, open_run)["state"] == "failed", "a unit that was not done still fails the run"


def test_a_unit_found_done_is_not_run_again(monkeypatch, conn):
    ran = []
    _handlers(monkeypatch, test=lambda args, report: ran.append(args["n"]) or {})
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=2)])
    runrecord.begin_stage(conn, run_id, 0)
    units.enqueue_units(conn, run_id, 0, "test", [({"n": 1}, 1.0), ({"n": 2}, 1.0)])
    from planetgen.queue import workers
    workers.ensure_workers(redisqueue.BULK_QUEUE, conn)
    assert sorted(ran) == [1, 2]
    units.enqueue_units(conn, run_id, 0, "test", [({"n": 1}, 1.0)])         # RQ's retry of unit 0, after it was done
    workers.ensure_workers(redisqueue.BULK_QUEUE, conn)
    assert sorted(ran) == [1, 2]


def test_a_unit_survives_a_warm_that_outlasts_an_execution_record(monkeypatch, conn):
    """PERF.146's warm used to run after RQ created the unit's execution record, which expires 90 s after its last beat
    (the first beat is 30 s after the fork), so a warm of over a minute on a loaded host left a unit that finished
    before the first beat without its record, and RQ failed it with `ValueError: Execution ... not found in Redis`
    when it cleaned up. Here the warm drops every execution record of the job, as the expiry would."""
    from planetgen.cli import worker
    _handlers(monkeypatch, test=_square)
    seen = []

    def long_warm(job):
        for key in conn.scan_iter(match=f"rq:execution:{job.id}:*"):
            seen.append(key)
            conn.delete(key)

    monkeypatch.setattr(units, "warm_unit", long_warm)
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=1)])
    runrecord.begin_stage(conn, run_id, 0)
    units.enqueue_units(conn, run_id, 0, "test", [({"n": 3}, 1.0)])
    queue = redisqueue.queue(redisqueue.BULK_QUEUE, conn)
    worker.worker_class()([queue], connection=conn, serializer=redisqueue.JsonSerializer).work(
        burst=True, logging_level="WARNING")
    assert seen == [], "no execution record exists while the worker warms"
    assert queue.failed_job_registry.count == 0
    run = runrecord.read_run(conn, run_id)
    assert run["state"] != "failed" and run["stages"][0]["done"] == 1


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


def test_only_the_table_names_what_a_unit_runs(monkeypatch):
    assert callable(commands.unit_function("scatter-mass")) and callable(commands.expand_function("scatter-mass"))
    for kind in ("os.system", "nope"):
        with pytest.raises(KeyError):
            commands.unit_function(kind)
    for kind in ("scatter-luminosity", "scatter-band", "scatter-phenomena", "plan-skeleton"):
        assert callable(commands.unit_function(kind)) and callable(commands.expand_function(kind))
    monkeypatch.setitem(commands.REGISTRY, "an-inline-command", commands.CommandSpec("an-inline-command", "Runs here"))
    with pytest.raises(ValueError, match="inline"):
        commands.unit_function("an-inline-command")


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
    from tests.test_bright_star_scatter import _scatter_args, _seed_galaxy
    _seed_galaxy(mysql_config)
    args = _scatter_args(mysql_config, "--mass", "--phenomenon-min-mass", "8", "--bright-star-min-luminosity", "2500",
                         "--workers", str(workers))
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM bright_stars")
        conn.commit()
    finally:
        conn.close()
    run_plan.scatter_bright_stars(args, ("mass",))
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
    assert [stage["key"] for stage in run["stages"]] == ["mass", "block-lights"]
    assert run["stages"][0]["done"] == run["stages"][0]["units"] > 0
    plan = runrecord.read_plan(conn, runs[0], 0)
    assert plan["kind"] == "scatter-mass" and plan["units"] == run["stages"][0]["units"]
    assert [step["kind"] for step in runrecord.read_steps(conn, runs[0])] == ["scatter-mass", "refresh-block-lights"]


def _scatter_then_band(mysql_config, workers):
    """Scatters the toy galaxy's bright stars (both passes), then adds the band down to 100 L_sun; returns the stars
    stored after each, sorted."""
    from planetgen.generation import run_plan
    from tests.test_bright_star_scatter import _scatter_args, _seed_galaxy
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--workers", str(workers)))
    scattered = _stored_stars(mysql_config)
    run_plan.add_bright_star_band(_scatter_args(mysql_config, "--band", "--down-to", "100", "--workers", str(workers)))
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
        "scatter-mass", "scatter-luminosity", "refresh-block-lights", "scatter-band"]
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
    from tests.test_bright_star_scatter import _scatter_args, _seed_galaxy
    store.get_control_connection(mysql_config, ensure_schema=True).close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    _seed_galaxy(mysql_config)
    args = _scatter_args(mysql_config, "--mass", "--workers", "2")
    stages.begin([stages.SCATTER_STAGE["mass"], stages.SCATTER_STAGE["luminosity"], stages.BLOCK_LIGHTS_STAGE], args,
                 "scatter")
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
    assert (mass["stage_n"], mass["stage_total"]) == (1, 3) and (light["stage_n"], light["stage_total"]) == (2, 3)
    assert mass["metrics"]["objects"] + light["metrics"]["objects"] == summary["total"] > 20
    assert mass["metrics"]["shells"] > 0 and 0 < mass["metrics"]["shells_modified"] <= mass["metrics"]["shells"]
    assert mass["metrics"]["seconds_per_shell"] == pytest.approx(mass["seconds"] / mass["metrics"]["shells"])
    assert {"scatter-mass", "scatter-bright"} <= kinds, "each unit's speed is recorded under its pass's kind"
    assert job and job[0]["metrics"]["shells"] > 0, "the whole job's totals are still the starter's"


@pytest.mark.fresh_database
def test_the_band_stage_is_stored_once_by_the_starter_because_it_tops_up_after_its_units(mysql_config, monkeypatch):
    from planetgen.db import store
    from planetgen.generation import run_plan
    from planetgen.generation import stats as generation_stats
    from tests.test_bright_star_scatter import _scatter_args, _seed_galaxy
    store.get_control_connection(mysql_config, ensure_schema=True).close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    run_plan.run_scatter(_scatter_args(mysql_config, "--band", "--down-to", "100", "--workers", "2"))
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
    from tests.test_bright_star_scatter import _scatter_args, _seed_galaxy
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config, "--phenomena", "--workers", str(workers)))
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
    assert [step["kind"] for run in runs for step in runrecord.read_steps(conn, run)] == [
        "scatter-phenomena", "finish-phenomena", "draw-nebula-field", "refresh-block-lights"]
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
    assert shape.to_json() == {key: value for key, value in expansion.plan.items()
                                      if key not in ("filled", "outline_shells", "min_mass_solar")} \
        and len(shape.tasks) == len(expansion.units) > 1
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


def _skeleton_args(mysql_config, *extra):
    import argparse

    from planetgen.cli import generate as generate_cli
    from tests.test_galaxy_gen import _PLAN_SHAPE_ARGV, _mysql_argv
    parser = argparse.ArgumentParser(prefix_chars="-+")
    generate_cli.add_plan_arguments(parser)
    args = parser.parse_args(_PLAN_SHAPE_ARGV + list(extra) + _mysql_argv(mysql_config))
    generate_cli.validate_plan_args(args, parser)
    return args


@pytest.fixture
def control_database(mysql_config, monkeypatch):
    """The test database doubling as the control database (the plan draws the galaxy's naming key into it)."""
    from planetgen.db import store
    store.get_control_connection(mysql_config, ensure_schema=True).close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    return mysql_config


@pytest.mark.fresh_database
def test_the_skeleton_is_one_unit_whose_lines_the_starter_logs_and_whose_run_is_recorded(
        control_database, monkeypatch, capsys):
    from planetgen.generation import run_plan
    runs = []
    real = wrapper.run
    monkeypatch.setattr(wrapper, "run", lambda *args, **kwargs: runs.append(real(*args, **kwargs)) or runs[-1])
    summary = run_plan._skeleton_on_queue(_skeleton_args(control_database))
    out = capsys.readouterr().out
    assert summary["shell_count"] > 0 and "Building skeleton:" in out and "Galaxy seed " in out and "Naming key " in out
    skeleton = next(line for line in out.splitlines() if "Building skeleton:" in line)     # API.33: every number has its unit
    assert " solar masses " in skeleton and " pc, " in skeleton and " degrees," in skeleton and "spiral arms" in skeleton
    assert "k_norm" not in skeleton and "threshold_rho" not in skeleton and "_pc=" not in skeleton
    drew = next(line for line in out.splitlines() if line.startswith("Drew ") or " Drew " in line)
    assert drew.count(" systems") == 3
    conn = redisqueue.connect()
    run = runrecord.read_run(conn, runs[0])
    assert run["state"] == "done" and [stage["key"] for stage in run["stages"]] == ["skeleton"]
    assert run["stages"][0]["units"] == 1
    assert [step["kind"] for step in runrecord.read_steps(conn, runs[0])] == ["plan-skeleton"]


@pytest.mark.fresh_database
def test_the_skeleton_unit_runs_alone_from_its_json_arguments(control_database):
    import json

    from planetgen.db import store
    from planetgen.generation import run_plan
    args = _skeleton_args(control_database, "--seed", "00112233445566778899AABBCCDDEEFF")
    params = json.loads(json.dumps(run_plan.skeleton_params(args, control_database)))
    assert params["galaxy_seed"] == "00112233445566778899AABBCCDDEEFF" and "password" not in params["database"]
    result = commands.run_alone("plan-skeleton", params)
    assert result["summary"]["shell_count"] > 0 and any(line.startswith("Galaxy seed 00112233") for line in result["lines"])
    conn = store.get_connection(control_database)
    try:
        assert store.get_galaxy_seed(conn).hex().upper() == "00112233445566778899AABBCCDDEEFF"
    finally:
        conn.close()


@pytest.mark.fresh_database
def test_a_skeleton_the_unit_refuses_ends_the_plan_with_its_message(control_database, caplog):
    from planetgen.generation import run_plan
    args = _skeleton_args(control_database)
    args.reach_shells = -1
    with pytest.raises(SystemExit) as stopped, caplog.at_level("ERROR"):
        run_plan._skeleton_on_queue(args)
    assert stopped.value.code == 1 and "--reach-shells -1" in caplog.text


@pytest.mark.fresh_database
def test_the_skeleton_stage_is_stored_once_by_its_unit(control_database, monkeypatch):
    from planetgen.db import store
    from planetgen.generation import run_plan, stages
    from planetgen.generation import stats as generation_stats
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    args = _skeleton_args(control_database, "--workers", "2")
    stages.begin(stages.plan_stages(args), args, "plan")
    stages.enter("skeleton")
    run_plan._skeleton_on_queue(args)
    stages.finish()
    conn = store.get_control_connection(control_database)
    try:
        rows = [run for run in generation_stats.stage_history(conn) if run["stage_key"] == "skeleton"]
    finally:
        conn.close()
    assert len(rows) == 1 and rows[0]["seconds"] > 0 and rows[0]["settings"]["workers"] == 2


def test_the_batch_fills_are_units_of_id_ranges_recorded_as_their_commands(mysql_config, monkeypatch):
    from planetgen.generation import run_plan
    from tests.test_map_luminosity import _seed_galaxy
    from tests.test_habitability_fill import _fill_args
    from tests.test_phenomenon_scatter import _scatter_args
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    monkeypatch.setattr(run_plan, "MAP_LUMINOSITY_BATCH_IDS", 5000)
    monkeypatch.setattr(run_plan, "HABITABILITY_BATCH_SYSTEMS", 1)
    runs = []
    real = wrapper.run
    monkeypatch.setattr(wrapper, "run", lambda *args, **kwargs: runs.append(real(*args, **kwargs)) or runs[-1])
    for fill, command in ((run_plan.fill_map_luminosity, "fill-map-luminosity"),
                          (run_plan.fill_habitability, "rescore-habitability")):
        args = _fill_args(mysql_config)
        args.command = command      # each fill is its own command's run: its list is that command's one stage
        fill(args)
    conn = redisqueue.connect()
    kinds = [[step["kind"] for step in runrecord.read_steps(conn, run)] for run in runs]
    assert kinds == [["fill-map-luminosity", "refresh-block-lights"], ["rescore-habitability"]]     # a luminosity fill changes the pre-placed objects: the block lights follow (PERF.172)
    stage = runrecord.read_run(conn, runs[0])["stages"][0]
    top = run_plan.expand_fill_map_luminosity(
        {"database": wrapper.database_args(mysql_config)}).plan["top"]
    assert stage["units"] == -(-top // 5000) > 1


def test_a_batch_fill_unit_runs_alone_from_its_json_arguments(mysql_config):
    import json

    from planetgen.generation import run_plan
    from tests.test_map_luminosity import _seed_galaxy
    from tests.test_phenomenon_scatter import _scatter_args
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    params = json.loads(json.dumps({"database": wrapper.database_args(mysql_config)}))
    expansion = run_plan.expand_fill_map_luminosity(params)
    assert sum(cost for _args, cost in expansion.units) == expansion.plan["top"]
    result = commands.run_alone("fill-map-luminosity", json.loads(json.dumps(expansion.units[0][0])))
    assert result["ids"] == expansion.units[0][1] and result["filled"] == 0


def test_the_nebula_field_is_a_unit_per_slab_and_the_clean_up_one_unit(mysql_config, monkeypatch):
    import json

    from planetgen.generation import run_plan
    from tests.test_nebula_field_table import _seed_galaxy, _table, no_superbubbles, plenty_of_clouds  # noqa: F401
    from tests.test_phenomenon_scatter import _scatter_args
    _seed_galaxy(mysql_config)
    runs = []
    real = wrapper.run
    monkeypatch.setattr(wrapper, "run", lambda *args, **kwargs: runs.append(real(*args, **kwargs)) or runs[-1])
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    conn = redisqueue.connect()
    drawn = [run for run in runs if "draw-nebula-field" in [step["kind"] for step in runrecord.read_steps(conn, run)]]
    assert len(drawn) == 1, "the nebula field is a stage of the phenomena run"
    [field] = [stage for stage in runrecord.read_run(conn, drawn[0])["stages"] if stage["key"] == "nebula-draw"]
    slabs = field["units"]
    assert slabs > 1 and sum(len(_table(mysql_config, f"cell_i = {i}")) for i in range(-slabs, slabs)) == len(_table(mysql_config))
    params = {"database": wrapper.database_args(mysql_config)}
    expansion = run_plan.expand_draw_nebula_field(params)
    assert len(expansion.units) == slabs and expansion.plan["slabs"] == slabs
    unit = json.loads(json.dumps(expansion.units[slabs // 2][0]))
    result = commands.run_alone("draw-nebula-field", unit)
    assert result["slab"] == unit["i"] and result["count"] == 1 and set(result) >= {"written", "slab", "count", "drawn"}
    assert [drawn["slab"] for drawn in result["drawn"]] == [unit["i"]] and set(result["drawn"][0]) >= {"counts", "shells"}
    check = commands.run_alone("clean-nebula-field", {**params, "apply": False}, lambda done, estimate: None)
    assert check["deleted"] == 0 and check["lines"][0].startswith("Nebula field cleanup (dry run)")


def test_a_galaxy_of_more_than_the_unit_cap_of_slabs_makes_bigger_nebula_field_units(mysql_config, monkeypatch):
    """PERF.192: one unit per 50 pc slab made a stage `enqueue_units` refused once a galaxy had more than 4000 slabs; the cut
    now makes a unit a run of consecutive slabs, the same cut at any worker count, every slab in exactly one unit."""
    import json

    from planetgen import tuning
    from planetgen.galaxy import nebula_field
    from planetgen.generation import run_plan
    from planetgen.queue import units
    from tests.test_nebula_field_table import _seed_galaxy
    cut = run_plan.nebula_slab_cut(-1500, 1500)
    assert len(cut) == 3001 and all(count == 1 for _start, count in cut), "a galaxy under the cap keeps a unit per slab"
    cut = run_plan.nebula_slab_cut(-10_000, 10_000)
    assert len(cut) <= units.MAX_UNITS and cut == run_plan.nebula_slab_cut(-10_000, 10_000)
    assert [slab for start, count in cut for slab in range(start, start + count)] == list(range(-10_000, 10_001))
    assert max(count for _start, count in cut) > 1
    _seed_galaxy(mysql_config)
    edge = tuning.NEBULA_FIELD_CELL_PC
    monkeypatch.setattr(nebula_field, "field_extent_pc", lambda edge_pc, extents: (5000 * edge, 100.0))
    expansion = run_plan.expand_draw_nebula_field({"database": wrapper.database_args(mysql_config)})
    assert expansion.plan["slabs"] == 10_001 and len(expansion.units) <= units.MAX_UNITS
    assert sum(cost for _args, cost in expansion.units) == expansion.plan["slabs"]
    json.dumps(expansion.units)


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


def test_the_sectors_are_cut_into_the_same_units_whatever_the_worker_count(monkeypatch):
    import json

    from planetgen.generation import run_galaxy
    items = [value for sector in range(25) for value in (3, sector % 4 - 1, sector)]
    params = {"database": {"host": "h", "port": 1, "user": "u", "database": "d"}, "args": {"num_systems": 2},
              "edge_pc": 10.0, "items": items, "workers": 1}
    for workers in (1, 8):
        assert run_galaxy.expand_generate_sectors({**params, "workers": workers}).units \
            == run_galaxy.expand_generate_sectors(json.loads(json.dumps(params))).units
    monkeypatch.setattr(run_galaxy, "SECTOR_UNITS_MAX", 4)
    expansion = run_galaxy.expand_generate_sectors(params)
    assert len(expansion.units) == 4 and expansion.plan == {"sectors": 25}
    assert [value for unit_args, _cost in expansion.units for value in unit_args["items"]] == items, "in order, whole"
    assert sum(cost for _args, cost in expansion.units) == 25


def test_the_arguments_of_a_sector_unit_carry_no_credentials_and_round_trip():
    import argparse
    import json

    from planetgen.generation import run_galaxy
    args = argparse.Namespace(num_systems=2, density=None, who="x", mysql_password="secret", mysql_host="h",
                              mysql_port=3306, mysql_user="u", mysql_database="d", _drawn=object(), tags=["a", 1],
                              stamp=(1, 2))
    plain = run_galaxy._sector_args_json(args)
    assert plain == {"num_systems": 2, "density": None, "who": "x", "tags": ["a", 1], "stamp": [1, 2]}
    unit_args = json.loads(json.dumps({"args": plain, "database": {"host": "h", "port": 3306, "user": "u",
                                                                  "database": "d"}}))
    rebuilt = run_galaxy._sector_args_of(unit_args)
    assert rebuilt.num_systems == 2 and rebuilt.mysql_password is None and rebuilt.mysql_database == "d"
    assert "secret" not in json.dumps(unit_args)


@pytest.mark.fresh_database
def test_a_galaxy_run_stores_each_of_its_stages_once_and_follows_its_sector_units(mysql_config, monkeypatch, capsys):
    from planetgen.db import store
    from tests.bughunt_support import mysql_argv, run_cli
    from tests.test_galaxy_gen import _seed_skeleton
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1), (2, 1)])
    runs = []
    real = wrapper.run
    monkeypatch.setattr(wrapper, "run", lambda *args, **kwargs: runs.append(real(*args, **kwargs)) or runs[-1])
    run_cli("galaxy", ["--shell", "2", "--num-systems", "2", *mysql_argv(mysql_config)])
    out = capsys.readouterr().out
    conn = redisqueue.connect()
    records = [runrecord.read_run(conn, run_id) for run_id in runs]
    keys = [stage["key"] for record in records for stage in record["stages"]]
    assert keys.count("sectors") == 1 and "link" in keys and "settle" in keys
    assert all(record["state"] == "done" for record in records)
    sectors = next(stage for record in records for stage in record["stages"] if stage["key"] == "sectors")
    assert sectors["units"] >= 1
    db = store.get_connection(mysql_config)
    try:
        saved = db.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"]
        paths = db.execute("SELECT COUNT(*) AS n FROM sector_paths").fetchone()["n"]
    finally:
        db.close()
    assert saved > 0 and out.count("Saved sector '") == saved, "each sector logged once, as it was saved"
    assert paths > 0 and "Linked " in out


def _backfilled_rows(config, monkeypatch, unit_cells):
    """The bright stars a one-sector `galaxy` run's backfill places, cut into units of `unit_cells` sectors, and the
    units it was cut into."""
    from planetgen.db import store
    from planetgen.generation import run_galaxy
    from planetgen.galaxy.spherical import sector_position_pc
    from tests import first_fill
    from tests.test_bright_star_scatter import EDGE_PC, EXTENTS, SEED, SHAPE, _database_now
    first_fill.without_the_deepest_floor(monkeypatch)   # the rows below are the backfill's alone (GEN.633)
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8,      # a sparse toy galaxy: a quick draw
                            expected_system_count_at_density_1=100.0, config=config, cell_star_ceiling=1e12,
                            galaxy_seed=SEED)
    store.replace_galaxy_shells(EXTENTS, config=config)
    monkeypatch.setattr(run_galaxy, "BACKFILL_UNIT_CELLS", unit_cells)
    args = run_galaxy._default_generation_args(config=config)
    args.num_systems = 0
    address = (3, 0, 4)
    started = _database_now(config)
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    args.shell, args.band, args.wedge = address
    cut = run_galaxy.expand_backfill_bright_stars({"database": wrapper.database_args(config),
                                                   "since": started.isoformat(sep=" ")})
    summary = run_galaxy._backfill_on_queue(args, started)
    conn = store.get_connection(config)
    try:
        rows = [tuple(row.values()) for row in conn.execute(
            "SELECT shell_index, band_index, wedge_index, initial_mass_sol, luminosity_w, star_type, position_x_mpc,"
            " position_y_mpc, position_z_mpc, seed"
            " FROM bright_stars ORDER BY shell_index, band_index, wedge_index, initial_mass_sol, seed"
        ).fetchall()]
    finally:
        conn.close()
    return rows, cut, summary


@pytest.mark.fresh_database
def test_the_backfill_gives_the_same_stars_in_one_unit_or_many(mysql_config, monkeypatch):
    from tests.test_gen_resume import second_mysql_config  # noqa: F401  (the fixture, to be called below)
    from planetgen.db import store
    import tests.conftest as conftest
    import pymysql
    import uuid
    kwargs = conftest._test_server_kwargs()
    name = f"planetgen_test_{uuid.uuid4().hex[:16]}"
    server = pymysql.connect(**kwargs)
    try:
        with server.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{name}`")
        server.commit()
        other = store.MySQLConfig(database=name, **kwargs)
        whole, whole_cut, whole_summary = _backfilled_rows(mysql_config, monkeypatch, 10_000)
        size = max(1, whole_cut.plan["sectors"] // 4)
        many, many_cut, many_summary = _backfilled_rows(other, monkeypatch, size)
        store.close_pool(other)
    finally:
        with server.cursor() as cur:
            cur.execute(f"DROP DATABASE IF EXISTS `{name}`")
        server.commit()
        server.close()
    assert len(whole_cut.units) == 1 and whole_cut.plan["sectors"] >= 4 and whole_cut.plan == many_cut.plan
    assert 2 <= len(many_cut.units) <= 5 and max(len(unit[0]["items"]) for unit in many_cut.units) == size
    assert whole and whole == many and whole_summary == many_summary and whole_summary["stars"] == len(whole)
    flat = [item for unit, _cost in many_cut.units for item in unit["items"]]
    assert flat == [item for unit, _cost in whole_cut.units for item in unit["items"]], "the same sectors, in order"


@pytest.mark.fresh_database
def test_a_backfill_with_nothing_to_draw_is_one_empty_unit_and_stores_its_stage(mysql_config):
    from planetgen.generation import run_galaxy
    from tests.test_bright_star_scatter import _database_now, _seed_galaxy
    _seed_galaxy(mysql_config)
    cut = run_galaxy.expand_backfill_bright_stars({"database": wrapper.database_args(mysql_config),
                                                   "since": _database_now(mysql_config).isoformat(sep=" ")})
    assert [unit[0]["items"] for unit in cut.units] == [[]] and cut.plan == {
        "sectors": 0, "skip": "the run generated no sector"}
    result = commands.run_alone("scatter-backfill", cut.units[0][0])
    assert result == {"sectors": 0, "stars": 0, "visited": 0, "shells": {}}


def _says_whether_missing_data_is_disregarded(args, report):
    from planetgen.datasets import fetch
    return {"disregarding": fetch.disregarding_missing_data()}


@pytest.mark.parametrize("told", [True, False])
def test_a_unit_runs_with_the_switch_its_run_was_started_with(monkeypatch, told):
    """GEN.599: a worker is a process of its own, so the run record carries the option to every unit."""
    from planetgen.datasets import fetch
    fetch.set_disregard_missing_data(told)
    try:
        run_id, results = _run(monkeypatch, _says_whether_missing_data_is_disregarded, count=3)
    finally:
        fetch.set_disregard_missing_data(False)
    assert [result["disregarding"] for result in results.values()] == [told] * 3


def test_the_shells_of_a_scatter_are_weighed_in_processes_to_the_same_weights_as_a_loop(monkeypatch):
    from types import SimpleNamespace

    from planetgen.generation import bright_stars, run_plan

    def weigh(shape, shell_index, high_band, edge_pc, e_value, fractions, beyond=None):
        return shell_index * 1.5 + high_band + (0.0 if beyond is None else beyond), shell_index * 3.0

    monkeypatch.setattr(bright_stars, "shell_weight", weigh)
    skeleton = SimpleNamespace(shape="shape", edge_pc=100.0)
    shells = [(shell, shell % 7 - 1) for shell in range(run_plan.SHELL_WEIGHT_PARALLEL_MIN * 3 + 5)]
    for estimates in (None, {shell: shell / 4 for shell, _band in shells}):
        loop = [weigh(None, shell, band, None, None, None, beyond=None if estimates is None else estimates[shell])
                for shell, band in shells]
        monkeypatch.setattr(run_plan, "SHELL_WEIGHT_PARALLEL_MIN", 10 ** 9)
        assert run_plan._shell_weights(skeleton, shells, 1.0, {}, estimates) == loop
        monkeypatch.setattr(run_plan, "SHELL_WEIGHT_PARALLEL_MIN", 2)
        monkeypatch.setattr(run_plan.os, "cpu_count", lambda: 4)
        assert run_plan._shell_weights(skeleton, shells, 1.0, {}, estimates) == loop
    assert run_plan._weighing is None
