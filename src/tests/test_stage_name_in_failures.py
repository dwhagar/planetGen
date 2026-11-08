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
    return types.SimpleNamespace(args=("r1", stage, 101), connection=conn, func_name=func_name)


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


def test_a_kill_while_a_stage_is_being_opened_names_both_stages(run):
    """PERF.252: the last unit of a stage opens the next in its own work horse, so a kill there belongs to the stage being opened."""
    runrecord.begin_opening(run, "r1", 2)
    units.horse_killed(_job(run, 1), 4242, 9, None)
    error = run.hget(runrecord.run_key("r1"), "error").decode()
    assert error.startswith("unit 101 of stage Scatter the massive stars was killed while it was opening stage Scatter the bright stars: its worker process 4242")


def test_a_failure_while_a_stage_is_being_opened_names_both_stages(run):
    runrecord.begin_opening(run, "r1", 2)
    units.unit_failed(_job(run, 1), run, ValueError, ValueError("no galaxy"), None)
    error = run.hget(runrecord.run_key("r1"), "error").decode()
    assert error == "unit 101 of stage Scatter the massive stars failed while it was opening stage Scatter the bright stars: ValueError: no galaxy"


def test_the_opening_mark_is_cleared_when_the_opening_ends_or_fails(run, monkeypatch):
    assert runrecord.opening(run, "r1") is None and wrapper.opening_note(run, "r1", 1) == ""
    seen = []

    def opened(conn, run_id, index):
        seen.append(runrecord.opening(conn, run_id))
        raise RuntimeError("cannot cut")
    monkeypatch.setattr(wrapper, "_open", opened)
    with pytest.raises(RuntimeError):
        wrapper._expand(run, "r1", 2)
    assert seen == [2] and runrecord.opening(run, "r1") is None
    runrecord.begin_opening(run, "r1", 1)
    assert wrapper.opening_note(run, "r1", 1) == ""                 # a unit of the stage itself is not an opening
    runrecord.end_opening(run, "r1")


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
