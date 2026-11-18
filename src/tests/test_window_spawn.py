"""OPS.240: a stage with more units than a window holds is cut into its first window of unit rows and a `spawn` row that holds the cursor; `planrows.spawn_window` writes the next window when the units written and
not finished fall under half of one. Real MySQL and Redis; the commands are the test's own."""

import json

import pytest

from planetgen.db import store
from planetgen.queue import commands, itemfile, planrows, redisqueue, runrecord, wrapper
from tests.test_plan_rows import _jobs

pytestmark = [pytest.mark.fresh_database]

WINDOW, TOTAL = 6, 20


def _later(start, count, **_):
    return [({"n": n}, 1.0) for n in range(start, min(start + count, TOTAL))]


RAN = []


def _unit(args, report):
    RAN.append(args["n"])
    return {"n": args["n"]}


def _cut(params):
    return commands.Expansion([({"n": n}, 1.0) for n in range(params["window"])], {"count": TOTAL}, total=TOTAL, cost=float(TOTAL), later={"made": True})


@pytest.fixture(autouse=True)
def _windowed(monkeypatch):
    monkeypatch.setitem(commands.REGISTRY, "windowed", commands.CommandSpec("windowed", "A test", unit=(__name__, "_unit"), expand=(__name__, "_cut"), window=(__name__, "_window")))
    monkeypatch.setattr(commands, "window_function", lambda kind: lambda later, start, count: _later(start, count))
    monkeypatch.setattr(planrows, "window_size", lambda workers: WINDOW)


def _started(mysql_config):
    conn = redisqueue.connect()
    step = commands.Step("a", "windowed", "A", {"database": wrapper.database_args(mysql_config)})
    run_id = runrecord.start_run(conn, "Windows", [runrecord.Stage("a", "A", units=TOTAL)], steps=[step.as_dict()])
    planrows.record_plan(conn, run_id, [step])
    first = _cut({"window": WINDOW})
    planrows.record_cut(conn, run_id, 0, "windowed", first.units, float(TOTAL), first.plan, total=TOTAL, later=first.later, window=WINDOW)
    planrows.release_units(conn, run_id, 0)
    return conn, run_id


def _rows(conn, run_id, mysql_config):
    ledger_run = int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))
    return _jobs(mysql_config, ledger_run)


def _units_of(rows):
    return [(int(job["item_key"]), job["state"]) for job in rows if job["command"] == "windowed" and job["item_key"]]


def _spawn(rows):
    return json.loads(next(job for job in rows if job["command"] == planrows.SPAWN)["arguments"])


def test_the_cut_writes_the_first_window_and_a_spawn_row_that_holds_the_cursor(mysql_config):
    conn, run_id = _started(mysql_config)
    rows = _rows(conn, run_id, mysql_config)
    assert _units_of(rows) == [(n, "ready") for n in range(WINDOW)]
    spawn = _spawn(rows)
    assert (spawn["next"], spawn["released"], spawn["total"], spawn["window"]) == (WINDOW, WINDOW, TOTAL, WINDOW)
    assert planrows.read_plan(conn, run_id, 0)["units"] == TOTAL, "the stage's units are all of them, written or not"
    assert next(job for job in rows if job["command"] == planrows.SPAWN)["state"] == "waiting"


def test_a_window_is_written_only_when_the_unfinished_rows_fall_under_half_of_one(mysql_config):
    conn, run_id = _started(mysql_config)
    done = runrecord.done_key(run_id, 0)
    assert planrows.spawn_window(conn, run_id, 0) == 0, "nothing is done: six rows wait, a window is not due"
    conn.sadd(done, 0, 1, 2)
    assert planrows.spawn_window(conn, run_id, 0) == 0, "three of six are not finished: that is not under half"
    conn.sadd(done, 3)
    assert planrows.spawn_window(conn, run_id, 0) == WINDOW
    rows = _rows(conn, run_id, mysql_config)
    assert _units_of(rows) == [(n, "ready") for n in range(2 * WINDOW)], "the next window is written and released in the order of the items"
    assert (_spawn(rows)["next"], _spawn(rows)["released"]) == (2 * WINDOW, 2 * WINDOW)
    assert planrows.spawn_window(conn, run_id, 0) == 0, "asked again at once it writes nothing: the rows written and not done are a whole window"


def test_a_unit_held_for_good_does_not_keep_a_window_from_being_written(mysql_config):
    conn, run_id = _started(mysql_config)
    conn.sadd(runrecord.done_key(run_id, 0), 0, 1)
    conn.hset(runrecord.held_key(run_id, 0), mapping={"2": "9", "3": "10"})
    assert planrows.spawn_window(conn, run_id, 0) == WINDOW


def test_the_last_window_is_the_rest_and_ends_the_spawn_row(mysql_config):
    conn, run_id = _started(mysql_config)
    done = runrecord.done_key(run_id, 0)
    written = 0
    for _ in range(5):
        conn.sadd(done, *range(0, planrows_next(conn, run_id, mysql_config)))
        written += planrows.spawn_window(conn, run_id, 0)
    rows = _rows(conn, run_id, mysql_config)
    assert [number for number, _state in _units_of(rows)] == list(range(TOTAL)) and written == TOTAL - WINDOW, "every item once, no gap"
    assert next(job for job in rows if job["command"] == planrows.SPAWN)["state"] == "done"


def planrows_next(conn, run_id, mysql_config):
    return _spawn(_rows(conn, run_id, mysql_config))["next"]


def test_a_window_written_twice_is_written_once(mysql_config):
    """A spawner stopped after the rows went in and before the cursor moved would write them again: the key of a row is unique, so the rows are the same ones."""
    conn, run_id = _started(mysql_config)
    conn.sadd(runrecord.done_key(run_id, 0), *range(WINDOW))
    control = planrows.control_connection(conn, run_id)
    try:
        run, stage = planrows.engine_run(conn, run_id, control), planrows._stage_key(conn, run_id, 0)
        from planetgen.queue import ledger
        ledger.add_jobs(control, run, [{"stage": stage, "command": "windowed", "item_key": str(WINDOW + n), "arguments": {"args": {"n": WINDOW + n}, "cost": 1.0}} for n in range(WINDOW)])
        control.commit()
    finally:
        control.close()
    planrows.spawn_window(conn, run_id, 0)
    assert [number for number, _state in _units_of(_rows(conn, run_id, mysql_config))] == list(range(2 * WINDOW))


def test_a_run_that_is_over_gets_no_window_and_its_spawn_row_is_canceled_with_its_item_file(mysql_config, tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    conn, run_id = _started(mysql_config)
    path = itemfile.path_of(run_id, "a")
    itemfile.write(path, [(1, 0, 0, 1, "")])
    conn.sadd(runrecord.done_key(run_id, 0), *range(WINDOW))
    conn.hset(runrecord.run_key(run_id), "state", "canceled")
    assert planrows.spawn_window(conn, run_id, 0) == 0
    rows = _rows(conn, run_id, mysql_config)
    assert len(_units_of(rows)) == WINDOW and next(job for job in rows if job["command"] == planrows.SPAWN)["state"] == "canceled"
    assert not tmp_path.joinpath("stage-items", f"{run_id}-a.items").exists()


def test_a_held_run_keeps_its_spawn_row_and_its_item_file(mysql_config, tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    conn, run_id = _started(mysql_config)
    path = itemfile.path_of(run_id, "a")
    itemfile.write(path, [(1, 0, 0, 1, "")])
    conn.hset(runrecord.run_key(run_id), "state", runrecord.HELD)
    assert planrows.spawn_window(conn, run_id, 0) == 0
    assert next(job for job in _rows(conn, run_id, mysql_config) if job["command"] == planrows.SPAWN)["state"] == "waiting"
    itemfile.Reader(path).close()


def test_the_item_file_goes_when_the_stage_is_released(mysql_config, tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    conn, run_id = _started(mysql_config)
    path = itemfile.path_of(run_id, "a")
    itemfile.write(path, [(1, 0, 0, 1, "")])
    planrows.mark_ended(conn, run_id, 0, planrows.CONCLUDE)
    itemfile.Reader(path).close()
    planrows.mark_ended(conn, run_id, 0, planrows.RELEASE)
    with pytest.raises(itemfile.ItemFileError, match="missing"):
        itemfile.Reader(path)


def test_a_run_goes_through_all_of_a_windowed_stage_one_window_after_another(mysql_config, monkeypatch):
    """The launcher's turn writes each next window as the units of the one before finish: every unit ran once, in the order of the items, and the stage's count is all of them."""
    RAN.clear()
    seen, results = [], {}
    run_id = wrapper.run("Windows", [commands.Step("a", "windowed", "A", {"database": wrapper.database_args(mysql_config)})], lambda index, plan: seen.append(plan["units"]),
                         lambda index, unit, result: results.__setitem__(unit, result))
    assert seen == [TOTAL] and sorted(results) == list(range(TOTAL)) and sorted(RAN) == list(range(TOTAL)) and len(RAN) == TOTAL
    conn = redisqueue.connect()
    rows = _rows(conn, run_id, mysql_config)
    assert {state for _number, state in _units_of(rows)} == {"done"} and len(_units_of(rows)) == TOTAL
    assert next(job for job in rows if job["command"] == planrows.SPAWN)["state"] == "done"
    assert runrecord.read_run(conn, run_id)["state"] == "done"


def test_a_stage_that_cannot_read_its_items_holds_the_run_and_a_retry_writes_the_window(mysql_config, monkeypatch):
    """The launcher's turn: a window that cannot be written (the item file is gone) is a lasting failure, not a stage with fewer units; the run is held and `Retry held` lets the spawner try again."""
    from planetgen.queue import held, launcher, memorygate
    conn, run_id = _started(mysql_config)
    conn.sadd(runrecord.done_key(run_id, 0), *range(WINDOW))
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        one = launcher.Launcher(control, None, read_host=lambda: memorygate.Host(10 ** 6, 10 ** 6, 64, 0.0), spawner=conn)
        monkeypatch.setattr(commands, "window_function", lambda kind: lambda later, start, count: (_ for _ in ()).throw(itemfile.ItemFileError("the item file x is missing")))
        assert one.spawn() == 0
        row = next(job for job in _rows(conn, run_id, mysql_config) if job["command"] == planrows.SPAWN)
        assert row["state"] == "failed" and "missing" in row["waited"]
        assert runrecord.run_state(conn, run_id) == runrecord.HELD
        assert len(_units_of(_rows(conn, run_id, mysql_config))) == WINDOW, "no unit was invented for the items it could not read"
        monkeypatch.setattr(commands, "window_function", lambda kind: lambda later, start, count: _later(start, count))
        held.retry_held(conn, run_id)
        assert next(job for job in _rows(conn, run_id, mysql_config) if job["command"] == planrows.SPAWN)["state"] == "waiting"
        one.spawned_at = None
        assert one.spawn() == WINDOW
    finally:
        control.close()
