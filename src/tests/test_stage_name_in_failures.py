"""PERF.221: the failed-unit and killed-unit messages name the stage the page shows ("Scatter the bright stars"), not the
internal stage number."""

import json
import types

import pytest

from planetgen.queue import redisqueue, runrecord, units, wrapper


@pytest.fixture
def run(redis_server, key_prefix, monkeypatch):
    """A running run of three stages in Redis under this test's own prefix."""
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    conn = redisqueue.connect(redis_server)
    stages = [{"label": "Plan the galaxy"}, {"label": "Scatter the massive stars"}, {"label": "Scatter the bright stars"}]
    conn.hset(runrecord.run_key("r1"), mapping={"state": "running", "kind": runrecord.GENERATION,
                                                "stages": json.dumps(stages)})
    return conn


def _job(conn, stage, func_name=units.UNIT_FUNCTION):
    args = ("r1", stage) if func_name == units.OPEN_FUNCTION else ("r1", stage, 101)
    return types.SimpleNamespace(args=args, connection=conn, func_name=func_name, id="job-1", worker_name="w1", retries_left=0)


def test_a_stage_is_named_by_its_label_and_by_its_number_only_when_the_run_cannot_say(run):
    assert wrapper.stage_name(run, "r1", 2) == "Scatter the bright stars"
    assert wrapper.stage_name(run, "r1", 7) == "7"                    # past the last stage
    assert wrapper.stage_name(run, "gone", 2) == "2"                  # a run Redis no longer holds


def test_a_failed_unit_names_the_stage_the_page_shows(run):
    units.unit_failed(_job(run, 2), run, ValueError, ValueError("no stars here"), None)
    error = run.hget(runrecord.run_key("r1"), "error").decode()
    assert error == "unit 101 of stage Scatter the bright stars failed: ValueError: no stars here"


def test_a_killed_unit_names_the_stage_the_page_shows(run):
    units.horse_killed(_job(run, 1), 4242, 9, None)
    error = run.hget(runrecord.run_key("r1"), "error").decode()
    assert error.startswith("unit 101 of stage Scatter the massive stars was killed: its worker process 4242 ended with status 9")


def test_a_kill_of_an_opening_names_the_stage_being_opened(run):
    """PERF.254: the cut of a stage is a job of its own, so a kill there is the opening's, not the unit's before it."""
    units.horse_killed(_job(run, 2, func_name=units.OPEN_FUNCTION), 4242, 9, None)
    error = run.hget(runrecord.run_key("r1"), "error").decode()
    assert error.startswith("the opening of stage Scatter the bright stars was killed: its worker process 4242 ended with status 9")


def test_the_opening_mark_is_cleared_when_the_opening_ends_or_fails(run, monkeypatch):
    assert runrecord.opening(run, "r1") is None
    seen = []

    def opened(conn, run_id, index):
        seen.append(runrecord.opening(conn, run_id))
        raise RuntimeError("cannot cut")
    monkeypatch.setattr(wrapper, "_open", opened)
    with pytest.raises(RuntimeError):
        wrapper._expand(run, "r1", 2)
    assert seen == [2] and runrecord.opening(run, "r1") is None


def test_an_opening_logs_each_step_with_the_memory_the_process_holds(monkeypatch):
    logged = []
    monkeypatch.setattr(wrapper.runrecord, "log", lambda conn, run_id, text, **fields: logged.append((text, fields)))
    readings = iter([1200.0, 1900.0, 1500.0])
    monkeypatch.setattr(wrapper.units, "resident_mb", lambda: next(readings))
    held = 0.0
    for text in ("is being cut into its units", "its preparation is done", "cut into 3 unit(s) after 4 s"):
        held = wrapper._say_opening(None, "r1", 9, text, held)
    assert held == 1900.0
    assert logged[0][0] == "Opening stage 10: is being cut into its units (this process holds 1,200 MB, at most 1,200 MB so far)"
    assert logged[2][0].endswith("(this process holds 1,500 MB, at most 1,900 MB so far)") and logged[2][1]["peak_mb"] == 1900
    assert all(fields["kind"] == "stage.opening" and fields["stage"] == 9 for _text, fields in logged)


def test_the_opening_asks_the_other_bulk_workers_to_go_and_the_run_waits_to_start_more(run, monkeypatch):
    """PERF.226: the stage before is complete, so the workers beside the one that opens the next stage only wait; they leave
    before the cut, and the run's `ensure` does not start them again until the cut is done."""
    shed = []
    monkeypatch.setattr(wrapper.workers, "shed_bulk", lambda conn, cap: shed.append(cap) or 0)
    seen = []
    monkeypatch.setattr(wrapper, "_open", lambda conn, run_id, index: seen.append(runrecord.opening(conn, run_id)) or 0)
    wrapper._expand(run, "r1", 2)
    assert shed == [1] and seen == [2]
    ensured = []
    clock = iter(range(0, 1000, 10))
    runrecord.begin_opening(run, "r1", 2)
    ticks = []

    def sleep(_seconds):
        ticks.append(1)
        if len(ticks) == 2:
            runrecord.end_opening(run, "r1")
        if len(ticks) > 5:
            raise RuntimeError("never ensured")
    monkeypatch.setattr(runrecord, "read_plan", lambda conn, run_id, index: {"units": 1} if ensured else None)
    result = wrapper._await_plan(run, "r1", 2, lambda: ensured.append(len(ticks)), clock=lambda: next(clock), sleep=sleep, touch=False)
    assert result == {"units": 1} and ensured and ensured[0] >= 2          # not while the opening lasted


def test_a_failed_unit_reaches_error_log_with_its_run_stage_unit_kind_and_traceback(run, monkeypatch):
    """OPS.162: RQ's own copy of the traceback has no run, stage, unit or kind; the unit's failure line carries them as fields."""
    from planetgen.util import log
    seen = []
    monkeypatch.setattr(log, "record_error", lambda component, message, exc_info=None, **fields: seen.append((component, message, exc_info, fields)))
    job = _job(run, 2)
    job.args = ("r1", 2, 101, "sectors", {}, None)
    try:
        raise ValueError("no stars here")
    except ValueError as exc:
        units.unit_failed(job, run, ValueError, exc, exc.__traceback__)
    (component, message, exc_info, fields), = seen
    assert component == "queue"
    assert message == "Unit 101 of stage Scatter the bright stars failed: ValueError: no stars here"
    assert exc_info[0] is ValueError and exc_info[2] is not None
    assert fields == {"run_id": "r1", "stage": 2, "stage_name": "Scatter the bright stars", "unit": 101, "kind": "sectors"}
