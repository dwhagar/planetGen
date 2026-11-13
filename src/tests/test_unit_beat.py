"""OPS.136: the beat record a horse writes for the unit it runs (`planetgen/queue/unitbeat.py`): a thread writes it on the Redis server clock while the work loop only sets numbers."""

import json
import threading
import time

import pytest

from planetgen.queue import beatjudge, redisqueue, runrecord, unitbeat, units


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture
def run_id(conn):
    return runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=1)])


class _Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


def _stored(conn, run_id, unit=0):
    return json.loads(conn.hget(runrecord.progress_key(run_id, 0), str(unit)))


def test_the_record_carries_the_counter_the_phase_the_wait_and_the_server_clock(conn, run_id):
    clock = _Clock()
    beat = unitbeat.UnitBeat(conn, run_id, 0, 0, clock=clock, resident=lambda: 123.456)
    clock.now += 5.0
    beat(10, 40)
    beat.phase("drawing tier 3")
    clock.now += 2.0
    with beat.waiting("sky lock"):
        waiting = beat.record()
    free = beat.record()
    assert (waiting["done"], waiting["estimate"], waiting["phase"], waiting["waiting_on"]) == (10.0, 40.0, "drawing tier 3", "sky lock")
    assert free["waiting_on"] == "" and free["resident_mb"] == 123.5
    server = unitbeat.server_time(conn)
    assert abs(free["beat"] - server) < 5.0, "the beat is on the Redis server clock"
    assert free["beat"] - free["started"] == pytest.approx(7.0, abs=0.01), "started is worked out from this process's clock"
    assert free["beat"] - free["moved"] == pytest.approx(2.0, abs=0.01), "moved is when the counter last changed"


def test_the_record_carries_the_time_spent_waiting_since_the_counter_last_changed(conn, run_id):
    """OPS.143: the judge subtracts it, so a unit that has just left a long wait is not quiet for the time it spent waiting."""
    clock = _Clock()
    beat = unitbeat.UnitBeat(conn, run_id, 0, 0, clock=clock)
    beat(1, 10)
    with beat.waiting("sky lock"):
        clock.now += 90.0
        assert beat.record()["waited"] == 0.0, "a wait still going on is reported as waiting_on, not yet as time waited"
    clock.now += 10.0
    with beat.waiting("gate"):
        clock.now += 30.0
    stored = beat.record()
    assert stored["waited"] == pytest.approx(120.0) and unitbeat.beat_of(stored).waited == pytest.approx(120.0)
    beat(2, 10)
    assert beat.record()["waited"] == 0.0, "a change of the counter starts the count again"


def test_moved_is_absent_until_the_counter_first_changes_and_follows_only_a_change(conn, run_id):
    clock = _Clock()
    beat = unitbeat.UnitBeat(conn, run_id, 0, 0, clock=clock)
    assert "moved" not in beat.record()
    beat(0, 10)
    assert "moved" not in beat.record(), "a call that leaves the counter as it was moves nothing"
    clock.now += 3.0
    beat(4, 10)
    clock.now += 6.0
    beat(4, 10)
    assert beat.record()["beat"] - beat.record()["moved"] == pytest.approx(6.0, abs=0.01)


def test_a_phase_that_begins_again_counts_as_a_restart(conn, run_id):
    beat = unitbeat.UnitBeat(conn, run_id, 0, 0)
    beat.phase("write batch")
    assert beat.record()["restarts"] == 0
    beat.phase("write batch")
    beat.phase("write batch")
    assert beat.record()["restarts"] == 2
    beat.phase("read")
    assert beat.record()["restarts"] == 0


def test_a_running_unit_has_a_fresh_beat_while_its_work_loop_is_frozen(conn, run_id):
    frozen = threading.Event()
    with unitbeat.UnitBeat(conn, run_id, 0, 0, beat_seconds=0.1, tick_seconds=0.02) as beat:
        beat(3, 9)
        first = _stored(conn, run_id) if conn.hexists(runrecord.progress_key(run_id, 0), "0") else None
        deadline = time.monotonic() + 5.0
        while time.monotonic() < deadline:                  # the work loop does nothing from here on
            stored = conn.hget(runrecord.progress_key(run_id, 0), "0")
            if stored and json.loads(stored)["beat"] > (first["beat"] if first else 0.0) + 0.25:
                break
            frozen.wait(0.05)
        second = _stored(conn, run_id)
    assert first is None or second["beat"] > first["beat"], "the beat moved while the counter did not"
    assert second["done"] == 3.0 and second["moved"] is not None


def test_the_counter_is_written_between_beats_at_a_slower_pace_than_it_moves(conn, run_id):
    writes = []
    real = conn.hset

    class _Spy:
        def __getattr__(self, name):
            return getattr(conn, name)

        def hset(self, key, field, value):
            writes.append(value)
            return real(key, field, value)

    with unitbeat.UnitBeat(_Spy(), run_id, 0, 0, beat_seconds=60.0, counter_write_seconds=0.2, tick_seconds=0.02) as beat:
        end = time.monotonic() + 1.0
        count = 0
        while time.monotonic() < end:
            count += 1
            beat(count, 10_000_000)             # the loop reports as fast as it can
            time.sleep(0.001)
    assert 2 <= len(writes) <= 12, f"{len(writes)} writes in a second of constant reporting"


def test_a_record_reads_back_as_a_beat_for_the_judge(conn, run_id):
    clock = _Clock()
    beat = unitbeat.UnitBeat(conn, run_id, 2, 7, clock=clock)
    beat.phase("p")
    beat(5, 20)
    conn.hset(runrecord.progress_key(run_id, 2), "7", json.dumps(beat.record()))
    conn.hset(runrecord.progress_key(run_id, 2), "8", "not json")
    found = unitbeat.read(conn, run_id, 2)
    assert set(found) == {7} and isinstance(found[7], beatjudge.Beat)
    assert (found[7].done, found[7].estimate, found[7].phase) == (5.0, 20.0, "p")


def test_the_follower_reads_the_counter_of_a_unit_that_has_one_and_skips_one_that_has_none(conn, run_id):
    conn.hset(runrecord.progress_key(run_id, 0), "0", json.dumps({"done": 2, "estimate": 8, "beat": 1.0}))
    conn.hset(runrecord.progress_key(run_id, 0), "1", json.dumps({"done": 0, "estimate": 0, "beat": 1.0}))
    assert units._progress(conn, run_id, 0) == {0: (2.0, 8.0)}


def test_code_deep_in_a_unit_tags_a_wait_through_the_module_and_does_nothing_outside_a_unit(conn, run_id):
    with unitbeat.waiting("sky lock"):          # no unit here: nothing happens
        unitbeat.phase("ignored")
    with unitbeat.UnitBeat(conn, run_id, 0, 0, beat_seconds=60.0) as beat:
        unitbeat.phase("charting")
        with unitbeat.waiting("sky lock"):
            assert beat.record()["waiting_on"] == "sky lock"
        assert beat.record()["waiting_on"] == "" and beat.record()["phase"] == "charting"
    with unitbeat.waiting("mysql"):
        pass
    assert unitbeat._current is None


def test_the_record_carries_the_evidence_of_work_besides_the_counter(conn, run_id):
    """OPS.193: processor seconds, bytes read and written, and the statement the unit is in, so a unit that does not move its counter can still be seen to work."""
    from planetgen.db import statementmark
    beat = unitbeat.UnitBeat(conn, run_id, 0, 0)
    quiet = beat.record()
    assert quiet["cpu_s"] >= 0.0 and quiet["io_read"] >= 0 and quiet["io_write"] >= 0
    assert "sql" not in quiet and "sql_since" not in quiet, "no statement out, no field"
    sum(n * n for n in range(300_000))
    with statementmark.labelled("clear the scatter"), statementmark.statement("DELETE FROM scatter_masses"):
        busy = beat.record()
    assert busy["sql"] == "clear the scatter"
    assert busy["beat"] - 5.0 < busy["sql_since"] <= busy["beat"], "the statement began on the server clock, a moment before the beat"
    assert busy["cpu_s"] >= quiet["cpu_s"]
    with statementmark.statement("SELECT   1\n  FROM dual"):
        assert beat.record()["sql"] == "SELECT 1 FROM dual"
    assert statementmark.running() is None


def test_the_statement_mark_names_the_longest_running_statement_and_clears_on_error():
    from planetgen.db import statementmark
    seen = []

    def other():
        with statementmark.statement("UPDATE b SET x = 1"):
            seen.append(statementmark.running())
            gate.wait(5)

    gate = threading.Event()
    with pytest.raises(RuntimeError):
        with statementmark.statement("SELECT 1"):
            thread = threading.Thread(target=other)
            thread.start()
            while not seen:
                time.sleep(0.01)
            assert statementmark.running()[0] == "SELECT 1", "the one that has waited longest"
            raise RuntimeError("boom")
    assert statementmark.running()[0] == "UPDATE b SET x = 1"
    gate.set()
    thread.join()
    assert statementmark.running() is None
    assert len(statementmark.shorten("x " * 100)) == statementmark.LABEL_CHARS
