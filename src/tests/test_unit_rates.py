"""PERF.147: each command records how fast its units work (items a second), in the run and across runs, and every running unit's published view carries its own bar with an ETA."""

import json

import pytest

from planetgen.queue import redisqueue, runrecord, unitbeat, unitrates, unitwatch


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture
def run_id(conn):
    return runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=2)])


class _Clock:
    def __init__(self):
        self.now = 100.0

    def __call__(self):
        return self.now


def test_the_rate_of_a_command_is_the_items_over_the_seconds_its_finished_units_worked(conn, run_id):
    unitrates.note(conn, run_id, "scatter", 100, 20.0)
    unitrates.note(conn, run_id, "scatter", 300, 20.0)
    unitrates.note(conn, run_id, "link", 5, 10.0)
    assert unitrates.read(conn, run_id) == {"scatter": 10.0, "link": 0.5}


def test_a_unit_that_did_nothing_or_worked_under_a_second_is_not_noted(conn, run_id):
    unitrates.note(conn, run_id, "scatter", 0, 20.0)
    unitrates.note(conn, run_id, "scatter", 50, 0.2)
    unitrates.note(conn, run_id, "scatter", None, None)
    assert unitrates.read(conn, run_id) == {}


def test_a_unit_counts_the_seconds_it_worked_not_the_ones_it_waited(conn, run_id):
    clock = _Clock()
    beat = unitbeat.UnitBeat(conn, run_id, 0, 0, clock=clock)
    clock.now += 10
    beat(40, 100)
    with beat.waiting("gate"):
        clock.now += 30
    clock.now += 10
    beat(80, 100)
    assert beat.totals() == (80.0, 20.0)


def test_a_run_that_ends_done_teaches_the_server_and_a_failed_one_does_not(conn, run_id):
    unitrates.note(conn, run_id, "scatter", 100, 10.0)
    runrecord.finish_run(conn, run_id, "done")
    assert unitrates.persisted(conn) == {"scatter": 10.0}
    second = runrecord.start_run(conn, "T2", [runrecord.Stage("test", "T", units=1)])
    unitrates.note(conn, second, "scatter", 200, 10.0)
    runrecord.finish_run(conn, second, "done")
    assert unitrates.persisted(conn)["scatter"] == pytest.approx(15.0)          # half the stored rate stays
    third = runrecord.start_run(conn, "T3", [runrecord.Stage("test", "T", units=1)])
    unitrates.note(conn, third, "scatter", 1, 1.0)
    runrecord.finish_run(conn, third, "failed", error="x")
    assert unitrates.persisted(conn)["scatter"] == pytest.approx(15.0)


def test_the_rate_for_the_bar_is_this_runs_else_the_servers(conn, run_id):
    rate_of = unitrates.rate_of(conn, run_id)
    assert rate_of("scatter") is None
    conn.hset(unitrates.persisted_key(), "scatter", 4.0)
    assert rate_of("scatter") == 4.0
    unitrates.note(conn, run_id, "scatter", 100, 10.0)
    assert rate_of("scatter") == 10.0, "the current run counts for more than the stored figures"


def _beat(conn, run_id, unit, *, done, estimate, running):
    now = unitbeat.server_time(conn)
    conn.hset(runrecord.progress_key(run_id, 0), unit, json.dumps({"done": done, "estimate": estimate, "beat": now - 1.0, "started": now - running, "phase": "p",
                                                                   "waiting_on": "", "restarts": 0, "moved": now - 1.0}))


def test_each_published_unit_carries_its_bar_with_the_stored_rate_eta_at_the_start(conn, run_id):
    conn.hset(unitrates.persisted_key(), "test", 10.0)
    watch = unitwatch.UnitWatch(conn, run_id, 0, "test", rate_of=unitrates.rate_of(conn, run_id), stuck=False, record_error=lambda *a, **k: None)
    _beat(conn, run_id, 0, done=0.0, estimate=1000.0, running=2.0)
    _beat(conn, run_id, 1, done=0.0, estimate=100.0, running=2.0)
    watch.check()
    shown = unitwatch.read(conn, run_id, 0)
    assert shown[0]["progress"]["show"] is True and shown[0]["progress"]["source"] == "stored" and shown[0]["progress"]["eta_seconds"] == pytest.approx(100.0, abs=0.5)
    assert shown[1]["progress"]["show"] is False, "a unit expected to take 10 s shows no bar"
    _beat(conn, run_id, 1, done=40.0, estimate=100.0, running=20.0)
    watch.check()
    assert unitwatch.read(conn, run_id, 0)[1]["progress"]["show"] is True, "one that runs over 15 s gains its bar"


def test_without_a_stored_rate_the_published_unit_still_has_a_progress_record(conn, run_id):
    watch = unitwatch.UnitWatch(conn, run_id, 0, "test", stuck=False, record_error=lambda *a, **k: None)
    _beat(conn, run_id, 0, done=50.0, estimate=100.0, running=100.0)
    watch.check()
    progress = unitwatch.read(conn, run_id, 0)[0]["progress"]
    assert progress["show"] is True and progress["percent"] == 50.0 and progress["source"] == "live" and progress["eta_seconds"] == pytest.approx(100.0, abs=0.5)
