"""OPS.216: the plan is rows. The cut of a stage writes its units and its own done mark to the job ledger in one transaction; everything that reads a plan reads the rows
(`queue/planrows.py`). Real MySQL and Redis."""

import pytest

from planetgen.db import store
from planetgen.queue import commands, planrows, redisqueue, runrecord, units, wrapper
from tests.test_queue_commands import _register, _units

pytestmark = [pytest.mark.fresh_database]


def _jobs(mysql_config, run):
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        return control.execute("SELECT * FROM engine_jobs WHERE run_id = ? ORDER BY id", (run,)).fetchall()
    finally:
        control.close()


def _run(monkeypatch, mysql_config, count=3):
    _register(monkeypatch, first=_units(count))
    return wrapper.run("Rows", [commands.Step("a", "first", "A", {"x": 1, "database": wrapper.database_args(mysql_config)})], lambda index, plan: None,
                       lambda index, unit, result: None)


def test_the_cut_of_a_stage_is_its_unit_rows_and_its_own_done_row(monkeypatch, mysql_config):
    run_id = _run(monkeypatch, mysql_config, 3)
    conn = redisqueue.connect()
    ledger_run = int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))
    jobs = _jobs(mysql_config, ledger_run)
    assert [(job["stage"], job["command"], job["item_key"]) for job in jobs] == [("a", "open", ""), ("a", "first", "0"), ("a", "first", "1"), ("a", "first", "2")]
    assert jobs[0]["state"] == "done"
    assert planrows.read_plan(conn, run_id, 0) == {"kind": "first", "units": 3, "cost": 3.0, "detail": {"count": 3}}


def test_a_cut_that_was_made_is_not_made_again_and_gives_the_same_row_ids(monkeypatch, mysql_config):
    run_id = _run(monkeypatch, mysql_config, 2)
    conn = redisqueue.connect()
    again = planrows.record_cut(conn, run_id, 0, "first", [({"n": 9}, 1.0)] * 5, 5.0, {"other": True})
    assert len(again) == 2, "the retry of a killed opening finds the cut it made: two units, whatever it was given"
    assert planrows.read_plan(conn, run_id, 0)["detail"] == {"count": 2}


def test_a_cut_that_fails_while_inserting_leaves_no_rows_and_no_done_mark(monkeypatch, mysql_config):
    from planetgen.queue import ledger
    _register(monkeypatch, first=_units(2))
    conn = redisqueue.connect()
    step = commands.Step("a", "first", "A", {"database": wrapper.database_args(mysql_config)})
    run_id = runrecord.start_run(conn, "Rows", [runrecord.Stage("a", "A", units=2)], steps=[step.as_dict()])
    real = ledger.add_jobs

    def broken(control, run, jobs):
        real(control, run, jobs)
        raise RuntimeError("the connection went away after the insert")

    monkeypatch.setattr(ledger, "add_jobs", broken)
    with pytest.raises(RuntimeError):
        planrows.record_cut(conn, run_id, 0, "first", [({}, 1.0), ({}, 1.0)], 2.0, {})
    monkeypatch.setattr(ledger, "add_jobs", real)
    ledger_run = int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))
    assert _jobs(mysql_config, ledger_run) == [], "the insert and the done mark are one transaction: nothing of the cut stands"
    assert planrows.read_plan(conn, run_id, 0) is None
    assert len(planrows.record_cut(conn, run_id, 0, "first", [({}, 1.0), ({}, 1.0)], 2.0, {})) == 2


def test_each_unit_is_handed_to_the_queue_with_its_row(monkeypatch, mysql_config):
    seen = []
    real = units.enqueue_units

    def spy(conn, run_id, stage, kind, unit_list, queue_name=None, expected_seconds=None, rows=None):
        seen.append(rows)
        return real(conn, run_id, stage, kind, unit_list, queue_name=queue_name, expected_seconds=expected_seconds, rows=rows)

    monkeypatch.setattr(units, "enqueue_units", spy)
    run_id = _run(monkeypatch, mysql_config, 3)
    conn = redisqueue.connect()
    ledger_run = int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))
    unit_ids = [job["id"] for job in _jobs(mysql_config, ledger_run) if job["command"] == "first"]
    assert seen == [unit_ids] and len(unit_ids) == 3


def test_a_stage_with_no_unit_still_has_its_plan_as_a_row(monkeypatch, mysql_config):
    _register(monkeypatch, first=_units(0))
    run_id = wrapper.run("None", [commands.Step("a", "first", "A", {"database": wrapper.database_args(mysql_config)})], lambda index, plan: None, lambda index, unit, result: None)
    assert planrows.read_plan(redisqueue.connect(), run_id, 0)["units"] == 0
