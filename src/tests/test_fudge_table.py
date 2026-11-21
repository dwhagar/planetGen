"""
The fudge table (GEN.1313, scripts/build_fudge_table.py): the committed file is what the script builds from the matrix, the class table and the registry; every quantity the class fit gives a
shoulder has a width; no width passes its cap; a width is the root sum of squares of the movements its rows make, each widened by the factor of its tier.
"""

import json
import math
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_fudge_table as builder  # noqa: E402


@pytest.fixture(scope="module")
def table():
    with open(builder.OUT, encoding="utf-8") as handle:
        return json.load(handle)


def test_the_committed_table_is_what_the_script_builds():
    assert open(builder.OUT, encoding="utf-8").read() == builder.make(), "run python scripts/build_fudge_table.py"


def test_the_build_is_reproducible():
    assert builder.make() == builder.make()


def test_every_shoulder_quantity_has_a_positive_width(table):
    from planetgen.physics import class_fit

    assert set(table["quantities"]) == set(class_fit.WIDTHS)
    for quantity, entry in table["quantities"].items():
        assert entry["width"] > 0.0 or not entry["edges"], quantity        # a quantity no class edge reads has no width to speak of
        assert entry["basis"] in ("measured", "fallback"), quantity
        assert (entry["tier"], entry["tweak"]) == (builder.TIER, builder.TWEAK)
        assert entry["unit"] == ("ln" if quantity in builder.LOG_QUANTITIES else "plain")


def test_a_measured_width_is_two_sigma_of_the_widened_movements(table):
    for quantity, entry in table["quantities"].items():
        if entry["basis"] != "measured":
            assert entry["rows"] == 0, quantity
            assert entry["sigma"] == 0.0, quantity
            continue
        assert entry["rows"] > 0, quantity
        assert entry["width"] == pytest.approx(builder.SIGMA_MULTIPLE * entry["sigma"], rel=1e-4, abs=1e-5), quantity
        assert sum(item["share"] for item in entry["contributors"]) <= 1.0 + 1e-3, quantity
        assert [item["share"] for item in entry["contributors"]] == sorted((item["share"] for item in entry["contributors"]), reverse=True), quantity


def test_no_width_in_use_passes_its_cap(table):
    for quantity, entry in table["quantities"].items():
        for edge in entry["edges"]:
            if edge["cap"] is None:
                assert edge["width"] == entry["width"]
                continue
            assert edge["width"] <= edge["cap"] + 1e-9, (quantity, edge)
            assert (edge["width"] == edge["cap"]) or not edge["binds"], (quantity, edge)        # a width that was cut is the cap


def test_tier_factor_widens_an_argued_range():
    assert builder.tier_factor("T1") == 1.0
    assert builder.tier_factor("T3") == 1.0
    assert builder.tier_factor("T4") == 1.5
    assert builder.tier_factor("T8") == 1.5
    assert builder.tier_factor("T9") == 2.0
    assert builder.tier_factor(None) == 2.0


def test_build_sums_independent_movements_in_quadrature():
    class Cls:
        def __init__(self, when):
            self.when = when

    classes = [Cls([{"quantity": "temperature_k", "min": 250.0}]), Cls([{"quantity": "temperature_k", "lt": 330.0}])]
    matrix = {"sample": {}, "rows": {"A1": {"quantities": {"temperature_k": {"rms": 3.0}}}, "B1": {"quantities": {"temperature_k": {"rms": 4.0}}}}}
    rows = [{"fill_id": "A1", "tier": "T2", "tweak": "L0"}, {"fill_id": "B1", "tier": "T9", "tweak": "L3"}]
    entry = builder.build(matrix, rows, classes)["quantities"]["temperature_k"]
    assert entry["sigma"] == pytest.approx(math.hypot(3.0, 8.0), rel=1e-5)
    assert entry["width"] == pytest.approx(2.0 * math.hypot(3.0, 8.0), rel=1e-5)
    assert entry["contributors"][0]["fill_id"] == "B1"
    assert entry["tuned_share"] == pytest.approx(64.0 / 73.0, abs=1e-3)
    assert all(edge["width"] <= 40.0 for edge in entry["edges"])        # half the 80 K gap between the two edges


def test_a_movement_past_the_half_gap_counts_as_the_half_gap():
    cell = {"rms": 50.0, "bodies": 10, "steps": [100.0, 2.0]}
    assert builder.clipped_rms(cell, None) == 50.0
    assert builder.clipped_rms(cell, 4.0) == pytest.approx(math.sqrt((16.0 + 4.0) / 10.0))
    assert builder.clipped_rms({"rms": 7.0}, 4.0) == 7.0


# --- The input stamp and the update.sh check (GEN.1317) ---------------------------------------------------------------------------

import copy  # noqa: E402
import subprocess  # noqa: E402


def _inputs():
    import stand_in_sensitivity

    with open(builder.MATRIX, encoding="utf-8") as handle:
        matrix = json.load(handle)
    return matrix, stand_in_sensitivity.read_rows()


def test_the_header_stamps_what_the_table_was_built_from(table):
    matrix, rows = _inputs()
    assert table["inputs"] == builder.input_hashes(matrix, rows)
    assert set(table["inputs"]) == {"registry", "class_table", "matrix", "matrix_epoch"}
    assert table["inputs"]["matrix_epoch"] == matrix["sample"]["epoch"], "the epoch is stamped in the matrix, not in the table: the table depends on files, and the epoch is numbered at merge"


def test_a_current_table_is_silent_about_its_files(table):
    ok, lines = builder.check_inputs(table=table, current_epoch=table["inputs"]["matrix_epoch"])
    assert ok and lines == []


def test_a_changed_input_is_named(table):
    matrix, rows = _inputs()
    other = copy.deepcopy(table)
    other["inputs"]["class_table"] = "0" * 64
    ok, lines = builder.check_inputs(table=other, current_epoch=table["inputs"]["matrix_epoch"])
    assert not ok and "class_table" in lines[0] and "build_fudge_table.py" in lines[0]
    assert "registry" not in lines[0] and "matrix," not in lines[0]
    # a registry row in scope that changes its tier moves the stamp; a row the matrix never measured does not
    measured = next(iter(matrix["rows"]))
    moved = [dict(row, tier="X9") if row["fill_id"] == measured else row for row in rows]
    assert builder.input_hashes(matrix, moved)["registry"] != table["inputs"]["registry"]
    unmeasured = [dict(row, tier="X9") if row["fill_id"] not in matrix["rows"] else row for row in rows]
    assert builder.input_hashes(matrix, unmeasured)["registry"] == table["inputs"]["registry"]


def test_an_old_matrix_is_reported_with_the_epochs_that_change_the_fill_stage(table):
    from planetgen.util import generator_epoch

    measured = table["inputs"]["matrix_epoch"]
    current = measured + 5
    ok, lines = builder.check_inputs(table=table, current_epoch=current)
    fill = len(generator_epoch.epochs_changing("fill", measured, current))
    assert ok, "an old matrix is reported, never refused"
    assert len(lines) == 1 and f"epoch {measured}" in lines[0] and f"5 epochs since, {fill} of them" in lines[0]


def test_a_missing_table_is_reported_not_raised(tmp_path, monkeypatch):
    monkeypatch.setattr(builder, "OUT", str(tmp_path / "none.json"))
    ok, lines = builder.check_inputs(current_epoch=1)
    assert not ok and "cannot be read" in lines[0]


def test_update_sh_runs_the_check_after_the_reference_pages_and_never_the_matrix():
    text = (ROOT / "update.sh").read_text(encoding="utf-8")
    step3 = text.split("log_stage 3 8")[1].split("log_stage 4 8")[0]
    assert step3.index("build_reference_pages") < step3.index("check_fudge_table")
    shell = (ROOT / "scripts" / "deploy-common.sh").read_text(encoding="utf-8")
    body = shell.split("check_fudge_table() {")[1].split("\n}\n")[0]
    assert "--check-inputs" in body and "stand_in_sensitivity" not in body and "--assemble" not in body


def _check_step(python, script):
    code = f'set -euo pipefail\nSCRIPT_DIR="{ROOT}"\nPYTHON="{python}"\nsource "$SCRIPT_DIR/scripts/deploy-common.sh"\n{script}\n'
    return subprocess.run(["bash", "-c", code], capture_output=True, text=True, timeout=120)


def test_update_sh_step_says_how_old_the_matrix_is_and_does_not_stop(tmp_path):
    result = _check_step(sys.executable, 'check_fudge_table; echo "failed=$FUDGE_TABLE_FAILED"')
    assert result.returncode == 0, result.stderr
    assert "failed=0" in result.stdout                       # an old matrix is a line, not a failure


def test_update_sh_step_warns_when_the_script_finds_other_files(tmp_path):
    fake = tmp_path / "python"
    fake.write_text("#!/bin/sh\necho 'the fudge table was built from other registry than the files now hold'\nexit 1\n")
    fake.chmod(0o755)
    result = _check_step(str(fake), 'check_fudge_table; echo "failed=$FUDGE_TABLE_FAILED"')
    assert result.returncode == 0, result.stderr
    assert "failed=1" in result.stdout and "built from other registry" in result.stderr
