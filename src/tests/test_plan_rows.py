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
    assert [(job["stage"], job["command"], job["item_key"]) for job in jobs] == [("a", "open", ""), ("a", "first", "0"), ("a", "first", "1"), ("a", "first", "2"),
                                                                                  ("a", "conclude", ""), ("a", "finish", ""), ("a", "release", "")]       # OPS.218: the end of the stage is rows too
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


def test_each_unit_is_handed_to_the_launcher_with_its_row(monkeypatch, mysql_config):
    released = []
    real = planrows.release_units

    def spy(conn, run_id, index, skip=()):
        released.append(real(conn, run_id, index, skip=skip))
        return released[-1]

    monkeypatch.setattr(planrows, "release_units", spy)
    run_id = _run(monkeypatch, mysql_config, 3)
    conn = redisqueue.connect()
    ledger_run = int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))
    unit_rows = [job for job in _jobs(mysql_config, ledger_run) if job["command"] == "first"]
    assert released == [3] and len(unit_rows) == 3, "the rows of the three units were released to the supervisor's launcher (OPS.228), each unit its own row"


def test_a_stage_with_no_unit_still_has_its_plan_as_a_row(monkeypatch, mysql_config):
    _register(monkeypatch, first=_units(0))
    run_id = wrapper.run("None", [commands.Step("a", "first", "A", {"database": wrapper.database_args(mysql_config)})], lambda index, plan: None, lambda index, unit, result: None)
    assert planrows.read_plan(redisqueue.connect(), run_id, 0)["units"] == 0


def test_the_cut_writes_the_jobs_that_end_the_stage_with_an_edge_each_and_the_run_marks_them_done(monkeypatch, mysql_config):
    run_id = _run(monkeypatch, mysql_config, 2)
    conn = redisqueue.connect()
    ledger_run = int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))
    ends = [(job["command"], job["needs"], job["state"]) for job in _jobs(mysql_config, ledger_run) if job["item_key"] == "" and job["command"] != "open"]
    assert ends == [("conclude", '["units"]', "done"), ("finish", '["conclude"]', "done"), ("release", '["finish"]', "done")]
    assert planrows.ended(conn, run_id, 0) == {planrows.CONCLUDE, planrows.FINISH, planrows.RELEASE}
    assert len(planrows.record_cut(conn, run_id, 0, "first", [], 0.0, {})) == 2, "the end rows are not unit rows"


def test_a_kind_with_end_commands_has_them_as_the_commands_of_its_end_rows(monkeypatch, mysql_config):
    """OPS.260: the `conclude` and `finish` rows of a stage are named for the stage's own end commands, and `ended` still reports the roles."""
    _register(monkeypatch, first=_units(2))
    CALLS.clear()
    monkeypatch.setitem(commands.END_COMMANDS, ("conclude", "first"), commands.EndCommand("conclude", "first", (__name__, "_conclusion")))
    monkeypatch.setitem(commands.END_COMMANDS, ("finish", "first"), commands.EndCommand("finish", "first", (__name__, "_metrics")))
    run_id = wrapper.run("Ends", [commands.Step("a", "first", "A", {"x": 1, "database": wrapper.database_args(mysql_config)})], lambda index, plan: None, lambda index, unit, result: None)
    conn = redisqueue.connect()
    ledger_run = int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))
    ends = [(job["command"], job["needs"], job["state"]) for job in _jobs(mysql_config, ledger_run) if job["item_key"] == "" and job["command"] != "open"]
    assert ends == [("conclude-first", '["units"]', "done"), ("finish-first", '["conclude-first"]', "done"), ("release", '["finish-first"]', "done")]
    assert planrows.ended(conn, run_id, 0) == {planrows.CONCLUDE, planrows.FINISH, planrows.RELEASE}
    assert CALLS == ["conclude", "finish"], "each end command ran once, in order"
    assert commands.end_name("first", "conclude") == "conclude-first" and commands.end_name("nothing-like-it", "finish") == "finish"


CALLS = []


def _conclusion(params, plan, results, seconds):
    CALLS.append("conclude")
    return {"units": len(results)}


def _metrics(params, plan, results, seconds):
    CALLS.append("finish")
    return {}


def test_a_database_named_by_its_name_alone_opens_on_the_configured_server():
    """A warm-up's steps name the galaxy database by its name, a generation run's by a dict (`wrapper.database_args`): both open (TEST.771)."""
    from planetgen.db import store
    configured = store.MySQLConfig()
    by_name = store.MySQLConfig.from_args("a_galaxy")
    assert by_name.database == "a_galaxy" and (by_name.host, by_name.port, by_name.user) == (configured.host, configured.port, configured.user)
    by_dict = store.MySQLConfig.from_args({"host": "db.example", "port": 3307, "user": "gen", "database": "other"})
    assert (by_dict.host, by_dict.port, by_dict.user, by_dict.database) == ("db.example", 3307, "gen", "other")

def test_the_whole_run_is_rows_before_its_first_stage_opens(monkeypatch, mysql_config):
    """OPS.238: a run's steps are written as `open` rows that wait, each needing the stage before it and the stages that make what it needs, and the cut of a stage fills its row in."""
    conn = redisqueue.connect()
    database = wrapper.database_args(mysql_config)
    steps = [commands.Step("mass-clear", "clear-bright-stars", "Clear", {"database": database}), commands.Step("mass", "scatter-mass", "Mass", {"database": database}),
             commands.Step("luminosity", "scatter-luminosity", "Light", {"database": database})]
    run_id = runrecord.start_run(conn, "Plan", [runrecord.Stage(step.key, step.label) for step in steps], steps=[step.as_dict() for step in steps])
    planrows.record_plan(conn, run_id, steps)
    planrows.record_plan(conn, run_id, steps)                 # written once: a second call leaves the rows as they are
    ledger_run = int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))
    jobs = _jobs(mysql_config, ledger_run)
    assert [(job["stage"], job["command"], job["state"]) for job in jobs] == [("mass-clear", "open", "waiting"), ("mass", "open", "waiting"), ("luminosity", "open", "waiting")]
    import json
    assert [json.loads(job["needs"]) for job in jobs] == [[], ["after:mass-clear"], ["after:mass"]]
    assert [json.loads(job["makes"]) for job in jobs] == [["bright-stars-cleared"], ["bright-stars-mass"], ["bright-stars-luminosity"]]
    assert planrows.read_plan(conn, run_id, 1) is None, "a planned stage is not a cut stage"


def test_the_cut_of_a_planned_stage_fills_its_row_in_and_marks_it_done(monkeypatch, mysql_config):
    _register(monkeypatch, first=_units(2))
    conn = redisqueue.connect()
    step = commands.Step("a", "first", "A", {"database": wrapper.database_args(mysql_config)})
    run_id = runrecord.start_run(conn, "Rows", [runrecord.Stage("a", "A", units=2)], steps=[step.as_dict()])
    planrows.record_plan(conn, run_id, [step])
    assert len(planrows.record_cut(conn, run_id, 0, "first", [({}, 1.0), ({}, 1.0)], 2.0, {"count": 2})) == 2
    ledger_run = int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))
    opened = [job for job in _jobs(mysql_config, ledger_run) if job["command"] == "open"]
    assert len(opened) == 1 and opened[0]["state"] == "done"
    assert planrows.read_plan(conn, run_id, 0) == {"kind": "first", "units": 2, "cost": 2.0, "detail": {"count": 2}}


def _planned_run(conn, mysql_config, kinds):
    database = wrapper.database_args(mysql_config)
    steps = [commands.Step(kind, kind, kind, {"database": database}) for kind in kinds]
    run_id = runrecord.start_run(conn, "Plan", [runrecord.Stage(step.key, step.label) for step in steps], steps=[step.as_dict() for step in steps])
    return run_id, steps


def _ended(conn, run_id, steps):
    """Cuts and ends every stage of the planned run with no unit, as a stage that has nothing to wait for ends."""
    for index, step in enumerate(steps):
        planrows.record_cut(conn, run_id, index, step.kind, [], 0.0, {})
        for role in planrows.END_ROLES:
            planrows.mark_ended(conn, run_id, index, role)
    runrecord.finish_run(conn, run_id, "done")


def test_what_the_rows_say_exists_meets_the_need_of_a_later_run(monkeypatch, mysql_config):
    """OPS.238: a need that no step of the run makes is met when the rows of the galaxy say a stage that ended made it; otherwise the plan says it is unmet."""
    conn = redisqueue.connect()
    first, steps = _planned_run(conn, mysql_config, ["plan-skeleton", "clear-bright-stars"])
    assert planrows.record_plan(conn, first, steps) == {}
    _ended(conn, first, steps)
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        assert planrows.made_artifacts(control, mysql_config.database) == {"skeleton", "bright-stars-cleared"}
    finally:
        control.close()
    second, steps = _planned_run(conn, mysql_config, ["scatter-mass", "scatter-luminosity"])
    assert planrows.record_plan(conn, second, steps) == {}, "the skeleton and the cleared scatter were made by the run before"
    runrecord.finish_run(conn, second, "done")
    third, steps = _planned_run(conn, mysql_config, ["scatter-band"])
    assert planrows.record_plan(conn, third, steps) == {"scatter-band": ["bright-stars-luminosity"]}, "no stage that ended made the luminosity pass"
    ledger_run = int(runrecord._s(conn.hget(runrecord.run_key(third), planrows.ENGINE_RUN_FIELD)))
    import json
    assert json.loads(_jobs(mysql_config, ledger_run)[0]["arguments"])["unmet"] == ["bright-stars-luminosity"]


def test_a_new_skeleton_forgets_what_the_galaxy_had():
    from planetgen.queue import needs
    assert needs.exists_after([("plan-skeleton", ["skeleton"]), ("scatter-mass", ["bright-stars-mass"]), ("plan-skeleton", ["skeleton"])]) == {"skeleton"}
    assert needs.exists_after([("scatter-mass", ["bright-stars-mass"]), ("reset-galaxy", [])]) == set()

def test_a_run_that_names_no_database_keeps_its_plan_in_the_ledger_of_the_configured_control_database(monkeypatch, mysql_config):
    """OPS.271: `system --output` touches no galaxy database, and its plan, cut and end marks are rows of the ledger all the same (the control database the machine is configured with)."""
    _register(monkeypatch, first=_units(2))
    run_id = wrapper.run("Rows", [commands.Step("a", "first", "A", {"x": 1, "database": None})], lambda index, plan: None, lambda index, unit, result: None)
    conn = redisqueue.connect()
    assert conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD) is not None
    assert planrows.read_plan(conn, run_id, 0) == {"kind": "first", "units": 2, "cost": 2.0, "detail": {"count": 2}}
    assert planrows.ended(conn, run_id, 0) >= {planrows.CONCLUDE, planrows.FINISH, planrows.RELEASE}


def test_a_run_with_no_control_database_is_refused_and_not_run_without_a_ledger(monkeypatch):
    """OPS.271: there is no run without the ledger; a control database that does not answer is the run's error."""
    def refuse(*args):
        raise OSError("no control database")
    monkeypatch.setattr(planrows, "control_connection", refuse)
    _register(monkeypatch, first=_units(2))
    with pytest.raises(OSError, match="no control database"):
        wrapper.run("Rows", [commands.Step("a", "first", "A", {"x": 1, "database": None})], lambda index, plan: None, lambda index, unit, result: None)


def test_the_units_of_a_cut_stage_are_released_for_the_launcher_except_the_ones_skipped(monkeypatch, mysql_config):
    """OPS.228: `release_units` makes the stage's waiting unit rows ready (not the stage's open or end rows), leaves the numbers in `skip`, and the open row names the Redis run and the stage's index."""
    _register(monkeypatch, first=_units(3))
    conn = redisqueue.connect()
    step = commands.Step("a", "first", "A", {"database": wrapper.database_args(mysql_config)})
    run_id = runrecord.start_run(conn, "Rows", [runrecord.Stage("a", "A", units=3)], steps=[step.as_dict()])
    planrows.record_plan(conn, run_id, [step])
    planrows.record_cut(conn, run_id, 0, "first", [({"n": n}, 1.0) for n in range(3)], 3.0, {})
    assert planrows.release_units(conn, run_id, 0, skip=[1]) == 2
    assert planrows.release_units(conn, run_id, 0) == 1, "the unit that was skipped is released by a later call, and the others are not released twice"
    import json
    jobs = _jobs(mysql_config, int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD))))
    assert {job["item_key"]: job["state"] for job in jobs if job["item_key"]} == {"0": "ready", "1": "ready", "2": "ready"}
    opened = next(job for job in jobs if job["command"] == "open")
    assert json.loads(opened["arguments"])["run_id"] == run_id and json.loads(opened["arguments"])["index"] == 0
    assert all(job["state"] == "waiting" for job in jobs if job["command"] != "open" and not job["item_key"])
