"""OPS.186 and OPS.187: the end of a stage can be killed or fail part way and is finished by the retry, once; a failure RQ reports after
an opening has cut its stage does not fail the run (`planetgen/queue/units.py`)."""

import types

import pytest
import rq

from planetgen.queue import jobprofiles, memorygate, redisqueue, runrecord, units, wrapper


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
    monkeypatch.setattr(wrapper, "conclude_and_finish_stage", lambda c, r, i, label=None: concluded.append(i))
    monkeypatch.setattr(units, "enqueue_open", lambda c, r, i, q: queued.append(i))
    return types.SimpleNamespace(run_id=run_id, concluded=concluded, queued=queued)


def _unit_job(run_id, retries_left):
    return types.SimpleNamespace(func_name=units.UNIT_FUNCTION, args=(run_id, 0, 0, "test", {}, 1.0), retries_left=retries_left)


def test_a_stage_end_that_cannot_queue_the_opening_tries_again_and_then_goes_on(conn, ended, monkeypatch):
    calls = []

    def hiccup(c, run_id, index, queue_name):
        calls.append(index)
        if len(calls) < 3:
            raise ConnectionError("Redis went away")
        ended.queued.append(index)

    monkeypatch.setattr(units, "enqueue_open", hiccup)
    units._end_stage(conn, ended.run_id, 0, wrapper, sleep=lambda seconds: None)
    assert ended.queued == [1] and ended.concluded == [0]
    assert not units.stage_unended(conn, ended.run_id, 0)


def test_a_retry_after_a_kill_finishes_the_end_without_concluding_or_queuing_twice(conn, ended, monkeypatch):
    monkeypatch.setattr(units, "enqueue_open", lambda c, r, i, q: (_ for _ in ()).throw(ConnectionError("down")))
    with pytest.raises(ConnectionError):
        units._end_stage(conn, ended.run_id, 0, wrapper, sleep=lambda seconds: None)
    assert units.stage_unended(conn, ended.run_id, 0) and ended.concluded == [0]
    monkeypatch.setattr(units, "enqueue_open", lambda c, r, i, q: ended.queued.append(i))
    monkeypatch.setattr(rq, "get_current_job", lambda: types.SimpleNamespace(connection=conn))
    result = units.run_unit(ended.run_id, 0, 0, "test", {}, 1.0)        # the retry: the unit is done, the end is not
    assert result == {"resumed": True} and ended.concluded == [0] and ended.queued == [1]
    assert units.run_unit(ended.run_id, 0, 0, "test", {}, 1.0) == {"skipped": True}, "now it is over"
    assert ended.queued == [1]


def test_a_second_attempt_does_not_end_a_stage_beside_the_first(conn, ended):
    conn.set(runrecord.ending_key(ended.run_id, 0), 1)
    units._end_stage(conn, ended.run_id, 0, wrapper)
    assert ended.concluded == [] and ended.queued == []


def test_a_unit_failure_with_the_end_unfinished_waits_for_its_retry_and_fails_the_run_when_none_is_left(conn, ended):
    units.unit_failed(_unit_job(ended.run_id, 1), conn, ConnectionError, ConnectionError("down"), None)
    assert runrecord.read_run(conn, ended.run_id)["state"] == "running", "the retry will finish the end"
    units.unit_failed(_unit_job(ended.run_id, 0), conn, ConnectionError, ConnectionError("down"), None)
    run = runrecord.read_run(conn, ended.run_id)
    assert run["state"] == "failed" and "opening of the next could not be queued" in run["error"]


def test_a_memory_kill_with_a_retry_left_frees_the_lease_of_the_end(conn, ended, monkeypatch):
    conn.set(runrecord.ending_key(ended.run_id, 0), 1)
    monkeypatch.setattr(memorygate, "kill_note", lambda c, worker: {"claimed_mb": 100.0, "peak_mb": 90.0, "kind": "test"})
    monkeypatch.setattr(memorygate, "grow", lambda *args: 200.0)
    monkeypatch.setattr(jobprofiles, "record", lambda *args, **kwargs: None)
    job = types.SimpleNamespace(connection=conn, id="j", func_name=units.UNIT_FUNCTION, args=(ended.run_id, 0, 0, "test", {}, 1.0), retries_left=1, worker_name="w")
    units.horse_killed(job, 4242, -9, None)
    assert not conn.exists(runrecord.ending_key(ended.run_id, 0)) and runrecord.read_run(conn, ended.run_id)["state"] == "running"


def _opening_job(run_id):
    return types.SimpleNamespace(func_name=units.OPEN_FUNCTION, args=(run_id, 1), retries_left=0)


def test_a_failure_rq_reports_after_an_opening_cut_its_stage_does_not_fail_the_run(conn, ended):
    units.open_failed(_opening_job(ended.run_id), conn, ValueError, ValueError("Execution not found in Redis"), None)
    assert runrecord.read_run(conn, ended.run_id)["state"] == "failed", "an opening that never cut its stage still fails the run"
    run_id = runrecord.start_run(conn, "T2", [runrecord.Stage("a", "A", units=1), runrecord.Stage("b", "B", units=1)])
    conn.set(runrecord.opened_key(run_id, 1), 1)
    units.open_failed(_opening_job(run_id), conn, ValueError, ValueError("Execution not found in Redis"), None)
    assert runrecord.read_run(conn, run_id)["state"] == "running"


def test_an_opening_whose_cut_is_done_is_not_cut_again_by_its_retry(conn, ended, monkeypatch):
    conn.set(runrecord.opened_key(ended.run_id, 1), 1)
    monkeypatch.setattr(rq, "get_current_job", lambda: types.SimpleNamespace(connection=conn))
    assert units.open_stage(ended.run_id, 1) == {"skipped": True}


def test_the_repair_finds_every_finished_stage_that_did_not_advance_and_ends_it_once(conn, ended):
    assert units.unended_stages(conn, ended.run_id) == [0]
    assert units.repair_unended_stages(conn, ended.run_id, wrapper) == [0]
    assert ended.concluded == [0] and ended.queued == [1]
    assert units.repair_unended_stages(conn, ended.run_id, wrapper) == [] and ended.queued == [1]
