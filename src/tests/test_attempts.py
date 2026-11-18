"""OPS.217: the fenced done mark. A unit takes an attempt in the ledger, every commit it fences makes is `WHERE attempt = mine`, its done mark is true in the galaxy before Redis hears of it, and a unit
another attempt took over rolls back and counts nothing (`planetgen/queue/attempts.py`)."""

import pytest

from planetgen.db import store
from planetgen.queue import attemptlog, attempts, commands, engine, ledger, planrows, redisqueue, runrecord, units, wrapper
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


# --- OPS.223: the history of every attempt, and the debug log of a restarted one ---------------------------------------------------------------------


def _one_unit(monkeypatch, mysql_config, handler):
    """A run with one cut unit whose attempts the test takes by hand: `(redis, run id, unit args, the unit's ledger row)`."""
    conn = redisqueue.connect()
    database = wrapper.database_args(mysql_config)
    run_id = runrecord.start_run(conn, "History", [runrecord.Stage("test", "T", units=1)], steps=[commands.Step("test", "test", "T", {"database": database}).as_dict()])
    _handlers(monkeypatch, test=handler)
    args = {"n": 1, "database": database}
    rows = planrows.record_cut(conn, run_id, 0, "test", [(args, 1.0)], 1.0, {})
    return conn, run_id, args, rows[0]


def _boom(args, report):
    raise ValueError("the unit broke")


@pytest.fixture
def log_folder(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(tmp_path))
    return tmp_path


def _history(mysql_config, row):
    control = _control(mysql_config)
    try:
        return control.execute("SELECT attempt, outcome, peak_mb, cpu_seconds, error_kind, error_text, stack_file, ended_at FROM engine_attempts WHERE job_id = ? ORDER BY attempt",
                               (row,)).fetchall()
    finally:
        control.close()


def test_an_attempt_that_finished_keeps_what_it_cost(monkeypatch, mysql_config, log_folder):
    conn, run_id, args, row = _one_unit(monkeypatch, mysql_config, _square)
    units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row)
    (done,) = _history(mysql_config, row)
    assert done["outcome"] == "done" and done["ended_at"] is not None
    assert done["peak_mb"] is not None and done["peak_mb"] >= 0 and done["cpu_seconds"] is not None and done["cpu_seconds"] >= 0
    assert done["stack_file"] is None and done["error_kind"] is None


def test_an_attempt_that_failed_keeps_its_kind_its_traceback_and_a_stack_file(monkeypatch, mysql_config, log_folder):
    conn, run_id, args, row = _one_unit(monkeypatch, mysql_config, _boom)
    with pytest.raises(ValueError):
        units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row)
    (failed,) = _history(mysql_config, row)
    assert failed["outcome"] == "failed" and failed["error_kind"] == "ValueError"
    assert "the unit broke" in failed["error_text"] and "_boom" in failed["error_text"], "the text is the traceback"
    assert failed["stack_file"] == attemptlog.stack_path(row, 1) and "_boom" in open(failed["stack_file"]).read()
    assert failed["stack_file"].startswith(str(log_folder)), "the file is in the log folder's unit-debug folder"


def test_the_time_limit_is_a_timed_out_attempt():
    class JobTimeoutException(Exception):
        pass

    kind, text, outcome = attempts.failure(JobTimeoutException("limit"))
    assert (kind, outcome) == ("JobTimeoutException", "timed_out") and "limit" in text
    assert attempts.failure(KeyError("x"))[2] == "failed"


def test_the_next_attempt_closes_one_nobody_ended_as_lost(mysql_config):
    control, galaxy, row, first, second = _two_attempts(mysql_config)
    try:
        lost, running = _history(mysql_config, row)
        assert lost["outcome"] == "lost" and lost["error_kind"] == "lost" and lost["ended_at"] is not None, "attempt 1 had no end written when attempt 2 began"
        assert running["outcome"] == "running"
    finally:
        control.close()
        galaxy.close()


def test_a_killed_horse_closes_its_attempt_in_the_ledger_and_the_job_waits_for_its_retry(monkeypatch, mysql_config, log_folder):
    conn, run_id, args, row = _one_unit(monkeypatch, mysql_config, _square)
    control = _control(mysql_config)
    try:
        ledger.begin_attempt(control, row)
    finally:
        control.close()
    assert attempts.close_lost(conn, run_id, row, "killed_memory", error_kind="killed", error_text="status -9", stack_text="a stack") is True
    (killed,) = _history(mysql_config, row)
    assert killed["outcome"] == "killed_memory" and killed["error_text"] == "status -9" and open(killed["stack_file"]).read() == "a stack"
    control = _control(mysql_config)
    try:
        assert control.execute("SELECT state FROM engine_jobs WHERE id = ?", (row,)).fetchone()["state"] == "waiting"
    finally:
        control.close()
    assert attempts.close_lost(conn, run_id, row, "lost") is False, "nothing is left running"
    assert attempts.close_lost(conn, run_id, None, "lost") is False, "a unit with no row has no history"


def test_a_restarted_attempt_writes_a_debug_log_that_opens_with_how_the_one_before_it_ended(monkeypatch, mysql_config, log_folder):
    from planetgen.queue import runrecord as run_record

    conn, run_id, args, row = _one_unit(monkeypatch, mysql_config, _boom)
    monkeypatch.setattr(attemptlog, "enabled", lambda: True)
    with pytest.raises(ValueError):
        units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row)
    assert not (log_folder / attemptlog.FOLDER / f"job{row}-unit0-attempt1.log").exists(), "a first attempt writes no debug log"
    _handlers(monkeypatch, test=_square)
    units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row)
    path = log_folder / attemptlog.FOLDER / f"job{row}-unit0-attempt2.log"
    text = path.read_text()
    assert f"debug log of attempt 2 of unit 0 of stage 0 of run {run_id}" in text
    assert "attempt 1 ended failed" in text and "ValueError" in text and "the unit broke" in text, "the header says how the first attempt ended and its stack"
    assert "unit 0 beat" in text or "starts in worker" in text, "the attempt's own debug lines follow"
    assert any(f"debug log is {path}" in record.message for record in run_record.log_records(conn, run_id)), "the run log names the file"


def test_the_setting_off_writes_nothing_extra(monkeypatch, mysql_config, log_folder):
    conn, run_id, args, row = _one_unit(monkeypatch, mysql_config, _boom)
    monkeypatch.setattr(attemptlog, "enabled", lambda: False)
    with pytest.raises(ValueError):
        units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row)
    _handlers(monkeypatch, test=_square)
    units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row)
    names = [path.name for path in (log_folder / attemptlog.FOLDER).iterdir()]
    assert names == [f"job{row}-attempt1-stack.txt"], "only the failed attempt's stack, no debug log of the restart"


def test_the_setting_is_on_by_default_and_the_folder_is_bounded(monkeypatch, log_folder):
    from planetgen.util import settings

    assert settings.Settings.model_fields["queue"].annotation.model_fields["debug_restarts"].default is True
    monkeypatch.setattr(attemptlog, "KEEP_FILES", 3)
    for number in range(1, 7):
        assert attemptlog.keep_stack(1, number, f"stack {number}") is not None
    assert len(list((log_folder / attemptlog.FOLDER).iterdir())) <= 3, "the oldest files are removed as new ones are made"
