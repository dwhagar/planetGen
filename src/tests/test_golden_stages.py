# tests/test_golden_stages.py

"""OPS.147: one digest per generation stage of the golden galaxy, and `scripts/train.py stages` checking what moved against what the train declared."""

import importlib.util
import json
import pathlib

import pytest

from planetgen.db import store

from tests import golden_seed, golden_stages
from tests.test_galaxy_gen import _plan_wide_galaxy
from tests.test_parallel_generation import make_database  # noqa: F401 -- a fixture

ROOT = pathlib.Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("train_for_stages", ROOT / "scripts" / "train.py")
train = importlib.util.module_from_spec(spec)
spec.loader.exec_module(train)


def _record(**rows):
    return {stage: (None if rows.get(stage) is None else {"full": rows[stage], "structural": "s", "rows": 3}) for stage in golden_stages.STAGES}


def test_the_stages_are_the_epoch_files_stages_in_order():
    assert golden_stages.STAGES == tuple(train.bump_version.GENERATION_STAGES)


def test_a_stage_moved_when_its_full_digest_differs_and_a_stage_nobody_measured_never_moves():
    before = _record(plan="a", scatter="b", fill="c", population=None)
    after = _record(plan="a", scatter="B", fill="C", population="p")
    assert golden_stages.moved(before, after) == ["scatter", "fill"]
    assert golden_stages.moved(before, before) == []


@pytest.mark.parametrize("moved, declared, failures, warnings", [
    (["fill"], ["fill"], [], []),
    (["fill", "population"], ["fill"], [], []),                       # fill covers the stage after it
    (["scatter", "fill", "population"], ["plan"], [], ["plan"]),      # plan covers all; it did not move: a warning
    (["plan", "fill"], ["fill"], ["plan"], []),                       # plan moved and nobody declared plan or anything before it
    (["scatter"], [], ["scatter"], []),                               # `generation-affects: none` with a moved stage
    ([], ["scatter", "fill"], [], ["scatter", "fill"]),
])
def test_a_moved_stage_must_be_covered_by_a_declared_one_and_a_declared_stage_that_stood_still_is_a_warning(moved, declared, failures, warnings):
    assert train.judge_stages(moved, declared) == (failures, warnings)


def test_the_stages_command_fails_a_moved_stage_nobody_declared_and_passes_a_covered_one(tmp_path, capsys):
    before, after = tmp_path / "before.json", tmp_path / "after.json"
    before.write_text(json.dumps(_record(plan="a", scatter="b", fill="c")))
    after.write_text(json.dumps(_record(plan="a", scatter="b", fill="C")))
    ok = train.main(["stages", "--before", str(before), "--after", str(after), "--affects", "fill"])
    assert ok == 0 and "fill        MOVED" in capsys.readouterr().out
    with pytest.raises(SystemExit) as stop:
        train.main(["stages", "--before", str(before), "--after", str(after), "--affects", "population"])
    assert stop.value.code == 1 and "FAIL fill moved" in capsys.readouterr().out


@pytest.mark.db
def test_the_stage_digests_of_a_stored_plan_are_stable_and_the_unmade_stages_are_empty(make_database, tmp_path):
    """A plan alone: the plan has rows; the scatter holds only the shape's columns; no sector means no fill; no facilities, so no population."""
    config = make_database()
    _plan_wide_galaxy(config, galaxy_seed=golden_seed.GALAXY_SEED)
    conn = store.get_connection(config)
    try:
        first = golden_seed.galaxy_stages(conn)
        second = golden_seed.galaxy_stages(conn)
        golden_stages.write(conn, tmp_path / "out.json")
    finally:
        conn.close()
    assert first == second == golden_stages.read(tmp_path / "out.json")
    assert set(first) == set(golden_stages.STAGES)
    assert first["plan"]["rows"] > 1 and first["fill"]["rows"] == 0 and first["population"] is None
