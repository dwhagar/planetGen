# tests/test_fudge_table_update.py

"""GEN.1317: the fudge table in the data folder is checked against the generator epoch and the inputs at update, built again from the stored matrix when they differ, and
refused (like a missing dataset) when a bulk generation starts without a current one."""

import json
import pathlib
import sys

import pytest

from planetgen.generation import run_galaxy
from planetgen.physics import class_fit, fudge_table
from planetgen.util import generator_epoch

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_fudge_table as builder  # noqa: E402

pytestmark = pytest.mark.real_datasets


@pytest.fixture
def data_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_DATA_DIR", str(tmp_path))
    return tmp_path


def _quiet(text):
    return None


def test_a_missing_table_is_built_and_a_current_one_is_left_alone(data_dir):
    assert fudge_table.state()[0] == "missing"
    assert builder.update(say=_quiet) == 0
    assert fudge_table.state() == ("current", "")
    first = (data_dir / "fudge_table.json").read_bytes()
    stamp = (data_dir / "fudge_table.json").stat().st_mtime_ns
    said = []
    assert builder.update(say=said.append) == 0
    assert (data_dir / "fudge_table.json").stat().st_mtime_ns == stamp and (data_dir / "fudge_table.json").read_bytes() == first
    assert "nothing to do" in said[0]


def test_the_header_holds_the_epoch_and_the_hashes_of_the_inputs(data_dir):
    builder.update(say=_quiet)
    header = fudge_table.header_of(fudge_table.data_file())
    assert header["epoch"] == generator_epoch.generator_epoch()
    assert {key: header[key] for key in ("registry_hash", "matrix_hash", "class_hash")} == {key: value for key, value in fudge_table.current_inputs().items() if key != "epoch"}
    assert "matrix_epoch" in header


def test_a_stale_epoch_stamp_makes_update_build_the_table_again(data_dir, monkeypatch):
    builder.update(say=_quiet)
    older = generator_epoch.generator_epoch()
    monkeypatch.setattr(generator_epoch, "generator_epoch", lambda: older + 1)
    found, reason = fudge_table.state()
    assert found == "stale" and "epoch" in reason
    assert builder.update(say=_quiet) == 0
    assert fudge_table.state()[0] == "current" and fudge_table.header_of(fudge_table.data_file())["epoch"] == older + 1


def test_an_edited_registry_row_makes_the_table_stale(data_dir, tmp_path):
    plan = tmp_path / "plan.md"
    text = pathlib.Path(fudge_table.PLAN).read_text(encoding="utf-8")
    plan.write_text(text, encoding="utf-8")
    first = fudge_table.registry_hash(str(plan))
    assert first == fudge_table.registry_hash()
    row = next(line for line in text.splitlines() if line.startswith("| FT3 |"))
    plan.write_text(text.replace(row, row + " (edited)"), encoding="utf-8")
    assert fudge_table.registry_hash(str(plan)) != first
    plan.write_text(text.replace("\n# ", "\n## ", 1), encoding="utf-8")      # prose outside the registry rows changes nothing
    assert fudge_table.registry_hash(str(plan)) == first


def test_the_built_table_is_the_committed_one_with_a_header_and_is_reproducible(data_dir):
    builder.update(say=_quiet)
    built = json.loads((data_dir / "fudge_table.json").read_text(encoding="utf-8"))
    header = built.pop("header")
    assert header and built == json.loads(pathlib.Path(builder.OUT).read_text(encoding="utf-8"))
    builder.update(force=True, say=_quiet)
    assert json.loads((data_dir / "fudge_table.json").read_text(encoding="utf-8")) == {**built, "header": header}


def test_a_file_with_no_header_or_a_broken_one_is_stale(data_dir):
    (data_dir / "fudge_table.json").write_text(pathlib.Path(builder.OUT).read_text(encoding="utf-8"), encoding="utf-8")
    assert fudge_table.state()[0] == "stale"
    (data_dir / "fudge_table.json").write_text("not json", encoding="utf-8")
    assert fudge_table.state()[0] == "stale"


def test_a_bulk_generation_is_refused_without_a_current_table(data_dir):
    with pytest.raises(run_galaxy.DatasetsMissing) as refused:
        run_galaxy.require_datasets()
    assert "fudge table is missing" in str(refused.value) and "nothing was generated" in str(refused.value) and "update.sh" in str(refused.value)
    builder.update(say=_quiet)
    fudge_table.require_current()          # a current table passes the table's own gate


def test_the_class_fit_reads_the_fixture_until_a_current_table_is_built(data_dir):
    assert fudge_table.widths_file() == fudge_table.FIXTURE
    builder.update(say=_quiet)
    assert fudge_table.widths_file() == fudge_table.data_file()
    assert class_fit.FUDGE_FILE == fudge_table.FIXTURE
