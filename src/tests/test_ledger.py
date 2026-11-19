"""OPS.215: the job engine's ledger (control schema v23, galaxy schema v261, planetgen/queue/ledger.py). Nothing else reads or writes the rows yet."""

import pytest

from planetgen.db import store
from planetgen.queue import ledger


@pytest.fixture
def control(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.execute("DELETE FROM engine_runs")
    conn.commit()
    try:
        yield conn
    finally:
        conn.execute("DELETE FROM engine_runs")
        conn.commit()
        conn.close()


@pytest.fixture
def galaxy(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        yield conn
    finally:
        conn.close()


def _plan(control, count=3, kind="generate", galaxy_name="g1"):
    run = ledger.add_run(control, kind, galaxy_name)
    ids = ledger.add_jobs(control, run, [{"stage": "fill", "command": "fill-unit", "item_key": f"part-{n}", "arguments": {"part": n}, "needs": ["scatter"],
                                          "makes": ["sectors"], "claim_mb": 900, "tokens": 2} for n in range(count)])
    control.commit()
    return run, ids


def _jobs(control, run):
    return control.execute("SELECT * FROM engine_jobs WHERE run_id = ? ORDER BY id", (run,)).fetchall()


def test_a_run_and_its_jobs_are_rows(control):
    run, added = _plan(control)
    row = control.execute("SELECT * FROM engine_runs WHERE id = ?", (run,)).fetchone()
    assert (row["kind"], row["galaxy"], row["state"], row["held_reason"]) == ("generate", "g1", "waiting", None)
    jobs = _jobs(control, run)
    assert added == 3 and [j["item_key"] for j in jobs] == ["part-0", "part-1", "part-2"]
    assert (jobs[0]["arguments"], jobs[0]["needs"], jobs[0]["makes"]) == ('{"part":0}', '["scatter"]', '["sectors"]')
    assert (jobs[0]["state"], jobs[0]["attempt"], jobs[0]["job_class"], jobs[0]["claim_mb"], jobs[0]["tokens"]) == ("waiting", 0, "bulk", 900, 2)


def test_adding_the_same_plan_again_adds_nothing(control):
    run, _ = _plan(control)
    again = ledger.add_jobs(control, run, [{"stage": "fill", "command": "fill-unit", "item_key": "part-1"}, {"stage": "fill", "command": "fill-unit", "item_key": "part-9"}])
    control.commit()
    assert again == 1 and len(_jobs(control, run)) == 4


def test_a_job_with_no_stage_or_an_unknown_field_is_refused(control):
    run = ledger.add_run(control, "generate")
    with pytest.raises(ValueError):
        ledger.add_jobs(control, run, [{"command": "x"}])
    with pytest.raises(ValueError):
        ledger.add_jobs(control, run, [{"stage": "s", "command": "x", "speed": 3}])


def test_a_run_is_stamped_when_it_starts_and_finishes_and_says_why_it_is_held(control):
    run, _ = _plan(control)
    ledger.set_run_state(control, run, "running")
    ledger.set_run_state(control, run, "held", "waiting for memory")
    row = control.execute("SELECT * FROM engine_runs WHERE id = ?", (run,)).fetchone()
    assert row["state"] == "held" and row["held_reason"] == "waiting for memory" and row["started_at"] is not None and row["finished_at"] is None
    ledger.set_run_state(control, run, "done")
    row = control.execute("SELECT * FROM engine_runs WHERE id = ?", (run,)).fetchone()
    assert row["held_reason"] is None and row["finished_at"] is not None
    with pytest.raises(ValueError):
        ledger.set_run_state(control, run, "paused")


def test_a_job_state_outside_the_set_is_refused_by_the_module_and_by_the_database(control):
    run, _ = _plan(control, 1)
    (job,) = _jobs(control, run)
    ledger.set_job_state(control, job["id"], "ready", waited="no memory", not_before=None)
    assert control.execute("SELECT state, waited FROM engine_jobs WHERE id = ?", (job["id"],)).fetchone() == {"state": "ready", "waited": "no memory"}
    with pytest.raises(ValueError):
        ledger.set_job_state(control, job["id"], "asleep")
    with pytest.raises(Exception):
        control.execute("UPDATE engine_jobs SET state = 'asleep' WHERE id = ?", (job["id"],))


def test_each_dispatch_raises_the_attempt_in_the_ledger_and_in_the_galaxy(control, galaxy):
    run, _ = _plan(control, 1)
    (job,) = _jobs(control, run)
    assert ledger.begin_attempt(control, job["id"], galaxy, pid=4242) == 1
    assert ledger.item_state(galaxy, job["id"]) == {"attempt": 1, "cursor_key": None, "done": False}
    ledger.end_attempt(control, job["id"], 1, "lost")
    control.commit()
    assert ledger.begin_attempt(control, job["id"], galaxy, pid=4243) == 2
    assert ledger.item_state(galaxy, job["id"])["attempt"] == 2
    rows = control.execute("SELECT attempt, pid, outcome, ended_at FROM engine_attempts WHERE job_id = ? ORDER BY attempt", (job["id"],)).fetchall()
    assert [(r["attempt"], r["pid"], r["outcome"], r["ended_at"] is not None) for r in rows] == [(1, 4242, "lost", True), (2, 4243, "running", False)]
    assert control.execute("SELECT state, attempt FROM engine_jobs WHERE id = ?", (job["id"],)).fetchone() == {"state": "running", "attempt": 2}


def test_a_run_is_running_from_the_start_of_its_first_job_and_a_held_or_ended_run_keeps_its_state(control):
    """OPS.281: a generation run stayed `waiting` with no `started_at` while its units ran."""
    run, _ = _plan(control, 2)
    first, second = _jobs(control, run)

    def row():
        return control.execute("SELECT state, started_at, finished_at FROM engine_runs WHERE id = ?", (run,)).fetchone()

    assert row() == {"state": "waiting", "started_at": None, "finished_at": None}
    ledger.begin_attempt(control, first["id"])
    began = row()
    assert began["state"] == "running" and began["started_at"] is not None and began["finished_at"] is None
    ledger.begin_attempt(control, second["id"])
    assert row()["started_at"] == began["started_at"], "the second job does not move the start"
    ledger.set_run_state(control, run, "held", "a unit failed")
    control.commit()
    ledger.end_attempt(control, second["id"], 1, "lost")
    ledger.begin_attempt(control, second["id"])
    assert row()["state"] == "held", "a held run is held until it is retried"


def test_an_attempt_keeps_what_it_cost_and_how_it_failed(control):
    run, _ = _plan(control, 1)
    (job,) = _jobs(control, run)
    attempt = ledger.begin_attempt(control, job["id"])
    ledger.end_attempt(control, job["id"], attempt, "failed", peak_mb=812.5, cpu_seconds=33.0, error_kind="ValueError", error_text="bad cell", stack_file="/tmp/s.txt")
    control.commit()
    row = control.execute("SELECT * FROM engine_attempts WHERE job_id = ?", (job["id"],)).fetchone()
    assert (row["outcome"], row["peak_mb"], row["cpu_seconds"], row["error_kind"], row["error_text"], row["stack_file"]) == ("failed", 812.5, 33.0, "ValueError", "bad cell", "/tmp/s.txt")
    with pytest.raises(ValueError):
        ledger.end_attempt(control, job["id"], attempt, "running")


def test_the_commit_of_the_attempt_the_row_holds_moves_the_cursor_and_sets_done(control, galaxy):
    run, _ = _plan(control, 1)
    (job,) = _jobs(control, run)
    attempt = ledger.begin_attempt(control, job["id"], galaxy)
    ledger.commit_item(galaxy, job["id"], attempt, cursor_key="cell:17")
    galaxy.commit()
    assert ledger.item_state(galaxy, job["id"]) == {"attempt": attempt, "cursor_key": "cell:17", "done": False}
    ledger.commit_item(galaxy, job["id"], attempt, done=True)        # a commit with no new cursor keeps the one it had
    galaxy.commit()
    assert ledger.item_state(galaxy, job["id"]) == {"attempt": attempt, "cursor_key": "cell:17", "done": True}


def test_a_replaced_attempt_touches_no_row_and_is_told_to_stop(control, galaxy):
    run, _ = _plan(control, 1)
    (job,) = _jobs(control, run)
    old = ledger.begin_attempt(control, job["id"], galaxy)
    new = ledger.begin_attempt(control, job["id"], galaxy)
    with pytest.raises(ledger.Replaced):
        ledger.commit_item(galaxy, job["id"], old, cursor_key="cell:99", done=True)
    galaxy.rollback()
    assert ledger.item_state(galaxy, job["id"]) == {"attempt": new, "cursor_key": None, "done": False}
    with pytest.raises(ledger.Replaced):
        ledger.commit_item(galaxy, job["id"] + 1000, 1)       # a job the supervisor never handed out has no row


def test_deleting_a_run_takes_its_jobs_and_attempts_with_it(control):
    run, _ = _plan(control, 2)
    for job in _jobs(control, run):
        ledger.begin_attempt(control, job["id"])
    control.execute("DELETE FROM engine_runs WHERE id = ?", (run,))
    control.commit()
    assert control.execute("SELECT COUNT(*) AS n FROM engine_jobs WHERE run_id = ?", (run,)).fetchone()["n"] == 0
    assert control.execute("SELECT COUNT(*) AS n FROM engine_attempts").fetchone()["n"] == 0


def test_the_writer_runs_queued_calls_in_order_and_survives_a_failed_one(mysql_config):
    opened = []

    def connect():
        conn = store.get_control_connection(mysql_config, ensure_schema=True)
        opened.append(conn)
        return conn

    made = []
    writer = ledger.Writer(connect)
    try:
        writer.submit(lambda conn: made.append(ledger.add_run(conn, "first")))
        writer.submit(lambda conn: conn.execute("INSERT INTO engine_runs (kind, state, created_at) VALUES ('bad', 'nonsense', NOW(6))"))
        writer.submit(lambda conn: made.append(ledger.add_run(conn, "third")))
        assert writer.flush(30.0)
    finally:
        writer.close()
    check = store.get_control_connection(mysql_config)
    try:
        kinds = [r["kind"] for r in check.execute("SELECT kind FROM engine_runs WHERE id IN (%s) ORDER BY id" % ",".join("?" * len(made)), tuple(made)).fetchall()]
        assert kinds == ["first", "third"] and len(opened) == 1
        assert check.execute("SELECT COUNT(*) AS n FROM engine_runs WHERE kind = 'bad'").fetchone()["n"] == 0
        check.execute("DELETE FROM engine_runs WHERE id IN (%s)" % ",".join("?" * len(made)), tuple(made))
        check.commit()
    finally:
        check.close()


OLD_SPELLING = "cancel" + "led"


def test_a_control_database_at_v23_gets_the_tables_again_with_the_state_spelled_canceled(mysql_config):
    """v24: the v23 tables, which held no rows, are dropped and made again, so the state `canceled` is allowed on a database that had the old check."""
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    try:
        for table in ("engine_attempts", "engine_jobs", "engine_runs"):
            conn.execute(f"DROP TABLE IF EXISTS {table}")
        conn.execute("CREATE TABLE engine_runs (id BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY, kind VARCHAR(64) NOT NULL, galaxy VARCHAR(64) NULL,"
                     " state VARCHAR(16) NOT NULL DEFAULT 'waiting', held_reason VARCHAR(255) NULL, created_at DATETIME(6) NOT NULL,"
                     " started_at DATETIME(6) NULL, finished_at DATETIME(6) NULL, CONSTRAINT chk_engine_runs_state CHECK (state IN ('waiting', 'running', 'held', 'done', '" + OLD_SPELLING + "')))")
        conn.execute("DELETE FROM control_schema_migrations")
        conn.execute("INSERT INTO control_schema_migrations (version) VALUES (23)")
        conn.commit()
        store._ensure_control_schema(conn)
        run = ledger.add_run(conn, "generate")
        ledger.set_run_state(conn, run, "canceled")
        conn.commit()
        assert conn.execute("SELECT state FROM engine_runs WHERE id = ?", (run,)).fetchone()["state"] == "canceled"
        assert conn.execute("SELECT MAX(version) AS v FROM control_schema_migrations").fetchone()["v"] == store.CONTROL_SCHEMA_VERSION
        conn.execute("DELETE FROM engine_runs")
        conn.commit()
    finally:
        conn.close()
