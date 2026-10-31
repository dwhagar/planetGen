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
