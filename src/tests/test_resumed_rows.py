"""OPS.222 part B: a run that carries on from a stopped one subtracts what the stopped run's rows say is done (`queue/resumed.py`). Real MySQL and Redis: the rows are the ledger's, the cursors
and the stored piece results are the galaxy's."""

import json

import pytest

from planetgen.db import store
from planetgen.queue import commands, ledger, planrows, redisqueue, resumed, runrecord, wrapper

pytestmark = pytest.mark.fresh_database

KIND = "resume-rows"
COUNT = 4


def listed(args, cursor):
    return list(range(args["count"]))[cursor:]


def work(args, piece, ctx):
    return None


@pytest.fixture(autouse=True)
def _command(monkeypatch):
    monkeypatch.setitem(commands.REGISTRY, KIND, commands.CommandSpec(KIND, "A test", expand=(__name__, "listed"), pieces=(__name__, "listed"), do=(__name__, "work"), resumes=True))
    monkeypatch.setitem(commands.REGISTRY, "restart-rows", commands.CommandSpec("restart-rows", "A test", expand=(__name__, "listed"), pieces=(__name__, "listed"), do=(__name__, "work")))


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture
def galaxy(mysql_config):
    connection = store.get_connection(mysql_config)
    yield connection
    connection.close()


def _units(mysql_config, count=COUNT):
    return [({"count": 6, "n": n, "database": wrapper.database_args(mysql_config)}, 1.0) for n in range(count)]


def _cut(conn, mysql_config, units, kind=KIND):
    step = commands.Step("s", kind, "S", {"database": wrapper.database_args(mysql_config)})
    run_id = runrecord.start_run(conn, "Rows", [runrecord.Stage("s", "S", units=len(units))], steps=[step.as_dict()])
    rows = planrows.record_cut(conn, run_id, 0, kind, units, float(len(units)), {})
    runrecord.finish_run(conn, run_id, "canceled")        # a run holds the machine while it lives: the next cut of the test is another run
    return run_id, rows


def _state(mysql_config, row):
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        return control.execute("SELECT state FROM engine_jobs WHERE id = ?", (row,)).fetchone()["state"]
    finally:
        control.close()


def _finished(conn, mysql_config, run_id, row, unit, result):
    """What a unit that was done leaves: its row done in the ledger and its result in the run's results."""
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        ledger.set_job_state(control, row, "done")
        control.commit()
    finally:
        control.close()
    conn.hset(runrecord.results_key(run_id, 0), unit, json.dumps(result))


def _checkpoint(galaxy, row, cursor, results):
    galaxy.execute("INSERT INTO engine_items (job_id, attempt, cursor_key, updated_at) VALUES (?, 1, ?, NOW(6))", (row, str(cursor)))
    for index, result in enumerate(results):
        if result is not None:
            galaxy.execute("INSERT INTO engine_piece_results (job_id, piece_index, result) VALUES (?, ?, ?)", (row, index, json.dumps(result)))
    galaxy.commit()


def _earlier_run(conn, mysql_config, galaxy):
    """A stopped run: unit 0 and 1 done, unit 2 stopped after two pieces, unit 3 never started."""
    units = _units(mysql_config)
    run_id, rows = _cut(conn, mysql_config, units)
    _finished(conn, mysql_config, run_id, rows[0], 0, {"n": 0})
    _finished(conn, mysql_config, run_id, rows[1], 1, {"n": 1})
    _checkpoint(galaxy, rows[2], 2, [{"piece": 0}, None])
    return run_id, rows, units


def test_the_units_the_earlier_run_finished_are_marked_done_and_their_results_are_found(conn, mysql_config, galaxy):
    prior, _rows, units = _earlier_run(conn, mysql_config, galaxy)
    new, rows = _cut(conn, mysql_config, units)
    found = resumed.carry(conn, new, 0, prior, KIND, rows, units)
    assert found.done == {0: ('{"n": 0}', 1.0), 1: ('{"n": 1}', 1.0)}
    assert [_state(mysql_config, row) for row in rows] == ["done", "done", "waiting", "waiting"]


def test_a_unit_that_stopped_partway_starts_after_the_pieces_it_finished(conn, mysql_config, galaxy):
    prior, _rows, units = _earlier_run(conn, mysql_config, galaxy)
    new, rows = _cut(conn, mysql_config, units)
    found = resumed.carry(conn, new, 0, prior, KIND, rows, units)
    assert found.resumed == 1
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        attempt = ledger.begin_attempt(control, rows[2], galaxy)
    finally:
        control.close()
    assert ledger.resume_point(galaxy, rows[2], attempt) == (2, [{"piece": 0}, None]), "the cursor and the stored results of the earlier unit"
    assert ledger.resume_point(galaxy, rows[3], 1) == (0, []), "a unit that never started has none to take"


def test_a_unit_whose_arguments_changed_is_not_taken_for_the_earlier_one(conn, mysql_config, galaxy):
    prior, _rows, units = _earlier_run(conn, mysql_config, galaxy)
    changed = [units[0], ({**units[1][0], "n": 99}, 1.0), ({**units[2][0], "count": 7}, 1.0), units[3]]
    new, rows = _cut(conn, mysql_config, changed)
    found = resumed.carry(conn, new, 0, prior, KIND, rows, changed)
    assert set(found.done) == {0} and found.resumed == 0, "only the unit that is the same unit is subtracted or resumed"
    assert [_state(mysql_config, row) for row in rows] == ["done", "waiting", "waiting", "waiting"]


def test_a_unit_the_rows_say_is_done_but_that_has_no_stored_result_runs_again(conn, mysql_config, galaxy):
    prior, _rows, units = _earlier_run(conn, mysql_config, galaxy)
    conn.hdel(runrecord.results_key(prior, 0), 1)
    new, rows = _cut(conn, mysql_config, units)
    assert set(resumed.carry(conn, new, 0, prior, KIND, rows, units).done) == {0}


def test_a_command_that_cannot_resume_takes_no_cursor(conn, mysql_config, galaxy):
    units = _units(mysql_config)
    prior, rows = _cut(conn, mysql_config, units, "restart-rows")
    _checkpoint(galaxy, rows[1], 3, [None, None, None])
    new, new_rows = _cut(conn, mysql_config, units, "restart-rows")
    found = resumed.carry(conn, new, 0, prior, "restart-rows", new_rows, units)
    assert found.resumed == 0 and ledger.resume_point(galaxy, new_rows[1], 1) == (0, [])


def test_a_run_with_no_rows_for_the_stage_gives_nothing_to_subtract(conn, mysql_config):
    units = _units(mysql_config)
    prior = runrecord.start_run(conn, "Old", [runrecord.Stage("s", "S", units=len(units))], steps=[commands.Step("s", KIND, "S", {"database": wrapper.database_args(mysql_config)}).as_dict()])
    runrecord.finish_run(conn, prior, "canceled")
    new, rows = _cut(conn, mysql_config, units)
    found = resumed.carry(conn, new, 0, prior, KIND, rows, units)
    assert found.done == {} and found.resumed == 0


def test_settling_counts_the_finished_units_in_the_stage(conn, mysql_config, galaxy):
    prior, _rows, units = _earlier_run(conn, mysql_config, galaxy)
    new, rows = _cut(conn, mysql_config, units)
    found = resumed.carry(conn, new, 0, prior, KIND, rows, units)
    runrecord.begin_stage(conn, new, 0)
    resumed.settle(conn, new, 0, found)
    assert conn.scard(runrecord.done_key(new, 0)) == 2
    assert sorted(int(runrecord._s(unit)) for unit in conn.hkeys(runrecord.results_key(new, 0))) == [0, 1]
