"""OPS.135: each command records the longest time its counter stood still, in the run and across runs, and the judge's quiet limit reads it."""

import pytest

from planetgen.queue import beatjudge, commandgaps, redisqueue, runrecord, unitbeat, unitwatch


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture
def run_id(conn):
    return runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=1)])


class _Clock:
    def __init__(self):
        self.now = 100.0

    def __call__(self):
        return self.now


def _beat(conn, run_id, clock):
    return unitbeat.UnitBeat(conn, run_id, 0, 0, clock=clock)


def test_the_gap_is_the_longest_stand_between_two_changes_from_start_to_end(conn, run_id):
    clock = _Clock()
    beat = _beat(conn, run_id, clock)
    clock.now += 30
    beat(1, 10)             # 30 s from the start to the first change
    clock.now += 5
    beat(2, 10)
    clock.now += 90
    beat(3, 10)             # 90 s: the longest
    beat(3, 10)             # no change, no new gap
    clock.now += 20
    assert beat.longest_gap() == pytest.approx(90.0)
    clock.now += 100        # the unit stands 120 s at the end before it finishes
    assert beat.longest_gap() == pytest.approx(120.0)


def test_time_spent_waiting_on_a_lock_is_not_standing_still(conn, run_id):
    clock = _Clock()
    beat = _beat(conn, run_id, clock)
    clock.now += 10
    beat(1, 10)
    with beat.waiting("sky lock"):
        clock.now += 1800       # waits half an hour for another unit's draw
    clock.now += 40
    beat(2, 10)
    assert beat.longest_gap() == pytest.approx(40.0)


def test_the_run_keeps_the_largest_gap_of_each_command_and_a_done_run_teaches_the_server(conn, run_id):
    commandgaps.note(conn, run_id, "sky", 120.0)
    commandgaps.note(conn, run_id, "sky", 80.0)
    commandgaps.note(conn, run_id, "sector", 3.0)
    commandgaps.note(conn, run_id, "sector", 0.0)
    commandgaps.note(conn, run_id, "none", None)
    assert commandgaps.read(conn, run_id) == {"sky": 120.0, "sector": 3.0}
    assert commandgaps.persisted(conn) == {}, "nothing is learned before the run ends"
    runrecord.finish_run(conn, run_id, "done")
    assert commandgaps.persisted(conn) == {"sky": 120.0, "sector": 3.0}


def test_a_run_that_failed_teaches_nothing_and_a_shorter_gap_pulls_the_stored_one_down_by_half(conn):
    first = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=1)])
    commandgaps.note(conn, first, "sky", 400.0)
    runrecord.finish_run(conn, first, "failed", error="x")
    assert commandgaps.persisted(conn) == {}
    runrecord.finish_run(conn, first, "failed", error="x")
    second = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=1)])
    commandgaps.note(conn, second, "sky", 400.0)
    runrecord.finish_run(conn, second, "done")
    third = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=1)])
    commandgaps.note(conn, third, "sky", 100.0)
    runrecord.finish_run(conn, third, "done")
    assert commandgaps.persisted(conn)["sky"] == pytest.approx(250.0)
    fourth = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=1)])
    commandgaps.note(conn, fourth, "sky", 900.0)
    runrecord.finish_run(conn, fourth, "done")
    assert commandgaps.persisted(conn)["sky"] == pytest.approx(900.0), "a longer gap replaces it"


def test_the_judge_gets_the_larger_of_the_servers_gap_and_this_runs(conn, run_id):
    history = commandgaps.history_of(conn, run_id)
    assert history("sky").longest_gap == 0.0
    conn.hset(commandgaps.persisted_key(), "sky", 200.0)
    assert history("sky").longest_gap == 200.0
    commandgaps.note(conn, run_id, "sky", 500.0)
    assert history("sky").longest_gap == 500.0 and history("sector").longest_gap == 0.0
    assert beatjudge.quiet_limit_seconds(history("sky").longest_gap) == 2500.0, "5 times the longest gap, over the 5 minute floor"


def test_a_unit_stands_quiet_only_past_five_times_its_commands_longest_gap(conn, run_id):
    commandgaps.note(conn, run_id, "sky", 600.0)
    watch = unitwatch.UnitWatch(conn, run_id, 0, "sky", history_of=commandgaps.history_of(conn, run_id), stuck=False)
    now = unitbeat.server_time(conn)
    import json
    for unit, still in ((0, 2000.0), (1, 3500.0)):
        conn.hset(runrecord.progress_key(run_id, 0), unit, json.dumps({"done": 5.0, "estimate": 100.0, "beat": now, "started": now - 4000, "moved": now - still,
                                                                       "phase": "p", "waiting_on": "", "restarts": 0}))
    verdicts = watch.check()
    assert verdicts[0].state == beatjudge.WORKING and verdicts[1].state == beatjudge.QUIET, "limit is 3000 s for a command that once stood 600 s"
