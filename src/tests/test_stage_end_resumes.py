"""OPS.186 and OPS.187: the end of a stage can be killed or fail part way and is finished by the retry, once; a failure RQ reports after
an opening has cut its stage does not fail the run (`planetgen/queue/units.py`)."""

import types

import pytest
import rq

from planetgen.queue import jobprofiles, memorygate, planrows, redisqueue, runrecord, units, wrapper


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture
def ended(conn, monkeypatch):
    """A run whose stage 0 has just been completed by its only unit (so the end of the stage is owed), with the conclusion counted."""
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("a", "A", units=1), runrecord.Stage("b", "B", units=1)])
    runrecord.begin_stage(conn, run_id, 0)
    assert runrecord.unit_done(conn, run_id, 0, 0).advance
    concluded, queued = [], []
    def conclude(c, r, i, label=None):
        concluded.append(i)
        planrows.mark_ended(c, r, i, planrows.CONCLUDE)
        return 0.0

    monkeypatch.setattr(wrapper, "conclude_stage_job", conclude)
    monkeypatch.setattr(wrapper, "finish_stage_job", lambda c, r, i, took=0.0: planrows.mark_ended(c, r, i, planrows.FINISH))
    monkeypatch.setattr(units, "enqueue_open", lambda c, r, i: queued.append(i))
    return types.SimpleNamespace(run_id=run_id, concluded=concluded, queued=queued)


def _unit_job(run_id, retries_left):
    return types.SimpleNamespace(func_name=units.UNIT_FUNCTION, args=(run_id, 0, 0, "test", {}, 1.0), retries_left=retries_left)


def test_a_stage_end_that_cannot_queue_the_opening_tries_again_and_then_goes_on(conn, ended, monkeypatch):
    calls = []

    def hiccup(c, run_id, index):
        calls.append(index)
        if len(calls) < 3:
            raise ConnectionError("Redis went away")
        ended.queued.append(index)

    monkeypatch.setattr(units, "enqueue_open", hiccup)
    units._end_stage(conn, ended.run_id, 0, wrapper, sleep=lambda seconds: None)
    assert ended.queued == [1] and ended.concluded == [0]
    assert not units.stage_unended(conn, ended.run_id, 0)


def test_a_retry_after_a_kill_finishes_the_end_without_concluding_or_queuing_twice(conn, ended, monkeypatch):
    monkeypatch.setattr(units, "enqueue_open", lambda c, r, i: (_ for _ in ()).throw(ConnectionError("down")))
    with pytest.raises(ConnectionError):
        units._end_stage(conn, ended.run_id, 0, wrapper, sleep=lambda seconds: None)
    assert units.stage_unended(conn, ended.run_id, 0) and ended.concluded == [0]
    monkeypatch.setattr(units, "enqueue_open", lambda c, r, i: ended.queued.append(i))
    units._end_stage(conn, ended.run_id, 0, wrapper, sleep=lambda seconds: None)        # the retry: the unit is done, the end is not
    assert ended.concluded == [0] and ended.queued == [1]
    units._end_stage(conn, ended.run_id, 0, wrapper, sleep=lambda seconds: None)
    assert ended.queued == [1], "now it is over"


def test_a_second_attempt_does_not_end_a_stage_beside_the_first(conn, ended):
    conn.set(runrecord.ending_key(ended.run_id, 0), 1)
    units._end_stage(conn, ended.run_id, 0, wrapper)
    assert ended.concluded == [] and ended.queued == []


def test_the_repair_finds_every_finished_stage_that_did_not_advance_and_ends_it_once(conn, ended):
    assert units.unended_stages(conn, ended.run_id) == [0]
    assert units.repair_unended_stages(conn, ended.run_id, wrapper) == [0]
    assert ended.concluded == [0] and ended.queued == [1]
    assert units.repair_unended_stages(conn, ended.run_id, wrapper) == [] and ended.queued == [1]
