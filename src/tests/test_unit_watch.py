"""OPS.138: the follower judges each running unit's beat, logs each change of state once, publishes the verdicts, and `GET /api/admin/work/<id>/units` shows them."""

import json

import pytest

from planetgen.queue import beatjudge, redisqueue, runrecord, unitbeat, unitwatch, units


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture
def run_id(conn):
    return runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=3)])


def _beat(conn, run_id, unit, *, age=1.0, still=0.0, done=5.0, estimate=100.0, waiting_on="", phase="p", running=2000.0):
    """Writes a beat record as a horse would, `age` seconds old, with the counter standing for `still` seconds."""
    now = unitbeat.server_time(conn)
    record = {"done": done, "estimate": estimate, "beat": now - age, "started": now - running, "phase": phase, "waiting_on": waiting_on, "restarts": 0,
              "moved": now - still}
    conn.hset(runrecord.progress_key(run_id, 0), unit, json.dumps(record))


def _lines(conn, run_id):
    return [fields["message"] for _, fields in runrecord.read_events(conn, run_id) if fields.get("kind") == "unit-state"]


def _watch(conn, run_id, errors):
    return unitwatch.UnitWatch(conn, run_id, 0, "test", record_error=lambda message, **fields: errors.append(message), stuck=False)


def test_a_unit_that_goes_quiet_is_logged_once_and_again_when_it_works(conn, run_id):
    errors = []
    watch = _watch(conn, run_id, errors)
    _beat(conn, run_id, 0)
    assert watch.check()[0].state == beatjudge.WORKING and _lines(conn, run_id) == [], "working from the start is not announced"
    for _ in range(3):
        _beat(conn, run_id, 0, still=2000.0)
        assert watch.check()[0].state == beatjudge.QUIET
    lines = _lines(conn, run_id)
    assert len(lines) == 1 and "unit 0 is quiet" in lines[0] and len(errors) == 1, "one line for three checks"
    _beat(conn, run_id, 0, still=1.0, done=9.0)
    assert watch.check()[0].state == beatjudge.WORKING
    assert len(_lines(conn, run_id)) == 2 and "works again" in _lines(conn, run_id)[1]


def test_a_dead_beat_is_dead_and_two_checks_in_a_row_agree(conn, run_id):
    watch = _watch(conn, run_id, [])
    _beat(conn, run_id, 1, age=300.0)
    _beat(conn, run_id, 2, still=2000.0)
    _beat(conn, run_id, 0)
    assert watch.check()[1].state == beatjudge.DEAD and watch.agreed() == [], "one check acts on nothing"
    _beat(conn, run_id, 1, age=330.0)
    _beat(conn, run_id, 2, still=2030.0)
    watch.check()
    assert watch.agreed() == [1, 2], "the working unit is never in it"


def test_a_unit_waiting_on_a_lock_is_waiting_not_quiet(conn, run_id):
    watch = _watch(conn, run_id, [])
    _beat(conn, run_id, 0, still=2000.0, waiting_on="sky lock")
    assert watch.check()[0].state == beatjudge.WAITING and watch.agreed() == []


def test_a_unit_that_finished_is_forgotten_and_a_new_attempt_is_announced_again(conn, run_id):
    watch = _watch(conn, run_id, [])
    _beat(conn, run_id, 0, still=2000.0)
    watch.check()
    conn.hdel(runrecord.progress_key(run_id, 0), 0)
    assert watch.check() == {} and unitwatch.read(conn, run_id, 0) == {}
    _beat(conn, run_id, 0, still=2000.0)
    watch.check()
    assert len(_lines(conn, run_id)) == 2, "a second attempt that goes quiet is a new change"


def test_the_verdicts_are_published_for_the_page_and_cleared_at_the_end_of_the_stage(conn, run_id):
    watch = _watch(conn, run_id, [])
    _beat(conn, run_id, 0, waiting_on="mysql", phase="drawing")
    watch.check()
    shown = unitwatch.read(conn, run_id, 0)[0]
    assert (shown["state"], shown["waiting_on"], shown["phase"], shown["done"], shown["estimate"]) == ("waiting", "mysql", "drawing", 5.0, 100.0)
    assert shown["beat_age_seconds"] >= 0.0 and shown["actionable"] is False
    watch.clear()
    assert unitwatch.read(conn, run_id, 0) == {}


def test_a_redis_that_cannot_answer_for_one_poll_changes_nothing(conn, run_id, monkeypatch):
    watch = _watch(conn, run_id, [])
    _beat(conn, run_id, 0, still=2000.0)
    first = watch.check()
    monkeypatch.setattr(unitbeat, "read", lambda *a: (_ for _ in ()).throw(ConnectionError("down")))
    assert watch.check() == first and len(_lines(conn, run_id)) == 1


def test_the_supervisor_judges_the_stage_in_flight_of_a_live_run_each_turn_and_drops_the_watch_when_it_ends(conn, run_id, monkeypatch):
    """OPS.208: the watch moved from the process that started the run to the supervisor (`supervisor.RunWatches`)."""
    from planetgen.queue import supervisor
    runrecord.begin_stage(conn, run_id, 0)
    runrecord.set_stage_plan(conn, run_id, 0, 1, 1.0, {"units": 1, "detail": {}})
    calls = []

    class _Watch:
        def check(self):
            calls.append("check")

        def clear(self):
            calls.append("clear")

    monkeypatch.setattr(runrecord, "read_steps", lambda c, r: [{"key": "test", "kind": "test"}])
    watches = supervisor.RunWatches(conn, make=lambda c, r, stage, kind: _Watch())
    assert supervisor.running_stage(conn, run_id) == 0
    assert watches.turn() == [(run_id, 0)] and watches.turn() == [(run_id, 0)] and calls == ["check", "check"]
    runrecord.unit_done(conn, run_id, 0, 0, seconds=1.0)
    assert supervisor.running_stage(conn, run_id) is None
    assert watches.turn() == [] and calls == ["check", "check", "clear"]


def test_a_failing_watch_does_not_stop_the_supervisors_turn(conn, run_id, monkeypatch):
    from planetgen.queue import supervisor
    runrecord.begin_stage(conn, run_id, 0)
    runrecord.set_stage_plan(conn, run_id, 0, 1, 1.0, {"units": 1, "detail": {}})
    monkeypatch.setattr(runrecord, "read_steps", lambda c, r: [{"key": "test", "kind": "test"}])

    class _Broken:
        def check(self):
            raise ConnectionError("down")

        def clear(self):
            pass

    watches = supervisor.RunWatches(conn, make=lambda c, r, stage, kind: _Broken())
    assert watches.turn() == [(run_id, 0)], "the failure is logged and the turn goes on"
