"""OPS.217: the fenced done mark. A unit takes an attempt in the ledger, every commit it fences makes is `WHERE attempt = mine`, its done mark is true in the galaxy before Redis hears of it, and a unit
another attempt took over rolls back and counts nothing (`planetgen/queue/attempts.py`)."""

import pytest

from planetgen.db import store
from planetgen.queue import attempts, commands, engine, ledger, planrows, redisqueue, runrecord, units, wrapper
from tests.test_work_units import _handlers, _square

pytestmark = pytest.mark.fresh_database


def _cut_with_database(database):
    return lambda params: commands.Expansion([({"n": n, "database": database}, 1.0) for n in range(params["count"])], {"count": params["count"]})


def _run(monkeypatch, mysql_config, count=2):
    _handlers(monkeypatch, test=_square)
    database = wrapper.database_args(mysql_config)
    monkeypatch.setattr(commands, "expand_function", lambda kind: _cut_with_database(database))
    return wrapper.run("Fence", [commands.Step("test", "test", "A test stage", {"count": count, "database": database})], lambda index, plan: None, lambda index, unit, result: None)


def _control(mysql_config):
    return store.get_control_connection(store.control_mysql_config(mysql_config))


def _ledger_run(run_id):
    return int(runrecord._s(redisqueue.connect().hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))


def test_each_unit_of_a_run_takes_an_attempt_and_is_marked_done_in_the_galaxy_and_the_ledger(monkeypatch, mysql_config):
    run_id = _run(monkeypatch, mysql_config, 2)
    control, galaxy = _control(mysql_config), store.get_connection(mysql_config)
    try:
        jobs = control.execute("SELECT id, state, attempt FROM engine_jobs WHERE run_id = ? AND command = 'test' ORDER BY id", (_ledger_run(run_id),)).fetchall()
        assert [(job["state"], job["attempt"]) for job in jobs] == [("done", 1), ("done", 1)]
        for job in jobs:
            item = ledger.item_state(galaxy, job["id"])
            assert item == {"attempt": 1, "cursor_key": None, "done": True}, "the done mark is in the galaxy's own database"
            assert control.execute("SELECT outcome FROM engine_attempts WHERE job_id = ?", (job["id"],)).fetchone()["outcome"] == "done"
    finally:
        control.close()
        galaxy.close()


def _two_attempts(mysql_config):
    """A job row and the galaxy connection, with attempt 1 taken and then attempt 2 (a retry started while the first was alive)."""
    control, galaxy = _control(mysql_config), store.get_connection(mysql_config)
    run = ledger.add_run(control, "generate", "g")
    ledger.add_jobs(control, run, [{"stage": "s", "command": "c", "item_key": "0"}])
    control.commit()
    row = control.execute("SELECT id FROM engine_jobs WHERE run_id = ?", (run,)).fetchone()["id"]
    config = mysql_config
    first = attempts.Attempt(control, row, ledger.begin_attempt(control, row, galaxy), config)
    second = attempts.Attempt(control, row, ledger.begin_attempt(control, row, galaxy), config)
    return control, galaxy, row, first, second


def test_a_replaced_attempts_commit_rolls_back_its_data_and_the_other_attempts_commit_lands(mysql_config):
    control, galaxy, row, first, second = _two_attempts(mysql_config)
    try:
        galaxy.execute("CREATE TABLE IF NOT EXISTS fence_probe (a INT)")
        galaxy.execute("INSERT INTO fence_probe VALUES (1)")
        with pytest.raises(attempts.Replaced):
            first.commit(galaxy)
        assert galaxy.execute("SELECT COUNT(*) AS n FROM fence_probe").fetchone()["n"] == 0, "the replaced attempt's data was rolled back"
        galaxy.execute("INSERT INTO fence_probe VALUES (2)")
        second.commit(galaxy)
        assert galaxy.execute("SELECT a FROM fence_probe").fetchall()[0]["a"] == 2
    finally:
        galaxy.execute("DROP TABLE IF EXISTS fence_probe")
        control.close()
        galaxy.close()


def test_only_the_attempt_that_still_owns_the_unit_sets_its_done_mark(mysql_config):
    control, galaxy, row, first, second = _two_attempts(mysql_config)
    try:
        with pytest.raises(attempts.Replaced):
            first.finish()
        assert ledger.item_state(galaxy, row)["done"] is False, "the replaced attempt marked nothing done"
        galaxy.commit()         # this connection read inside a transaction: start a new look at the row
        second.finish()
        galaxy.commit()
        assert ledger.item_state(galaxy, row)["done"] is True
    finally:
        control.close()
        galaxy.close()


def test_a_unit_whose_attempt_was_replaced_counts_nothing_and_ends_replaced(monkeypatch, mysql_config):
    run_id = runrecord.start_run(redisqueue.connect(), "Fence", [runrecord.Stage("test", "T", units=1)],
                                 steps=[commands.Step("test", "test", "T", {"database": wrapper.database_args(mysql_config)}).as_dict()])
    conn = redisqueue.connect()
    _handlers(monkeypatch, test=_square)
    database = wrapper.database_args(mysql_config)
    rows = planrows.record_cut(conn, run_id, 0, "test", [({"n": 1, "database": database}, 1.0)], 1.0, {})
    monkeypatch.setattr(attempts.Attempt, "finish", lambda self: (_ for _ in ()).throw(attempts.Replaced("another attempt has it")))
    assert units._do_unit(conn, run_id, 0, 0, "test", {"n": 1, "database": database}, 1.0, rows[0]) == {"replaced": True}
    assert not conn.sismember(runrecord.done_key(run_id, 0), 0), "the unit is not in the done set"
    control = _control(mysql_config)
    try:
        assert control.execute("SELECT outcome FROM engine_attempts WHERE job_id = ?", (rows[0],)).fetchone()["outcome"] == "replaced"
    finally:
        control.close()


def test_the_engines_commit_is_the_fence_when_the_unit_has_an_attempt(mysql_config):
    control, galaxy, row, first, second = _two_attempts(mysql_config)
    try:
        ctx = engine.Context("c", {"database": wrapper.database_args(mysql_config)}, lambda done, total: None, attempt=first)
        ctx._db = galaxy
        with pytest.raises(attempts.Replaced):
            ctx.commit()
        ctx.attempt = second
        ctx.commit()
        ctx._db = None
    finally:
        control.close()
        galaxy.close()


def test_a_unit_with_no_row_has_no_attempt_and_runs_as_before(monkeypatch):
    assert attempts.begin(redisqueue.connect(), "no-run", None, {"database": {}}) is None
