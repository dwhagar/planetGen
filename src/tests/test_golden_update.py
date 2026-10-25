# tests/test_golden_update.py

"""OPS.56: `python -m tests.golden_update` bumps the generator epoch, and refuses where it must."""

import json
import os
import shutil
import subprocess

import pytest

from tests import golden_update

needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")


def _git(root, *args):
    subprocess.run(["git", "-c", "user.email=t@example.com", "-c", "user.name=t", *args], cwd=root, check=True,
                   capture_output=True)


def _repo(tmp_path):
    root = str(tmp_path)
    os.makedirs(os.path.join(root, golden_update.GOLDEN_DIR))
    with open(os.path.join(root, golden_update.EPOCHS_FILE), "w") as f:
        json.dump({"epochs": [{"epoch": 1, "first_release": None, "reason": "first"}]}, f)
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write("one")
    _git(root, "init", "-q", "-b", "main")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "base")
    return root


def _epochs(root):
    return json.load(open(os.path.join(root, golden_update.EPOCHS_FILE)))["epochs"]


@needs_git
def test_a_changed_golden_file_bumps_the_epoch_and_leaves_it_pending(tmp_path):
    root = _repo(tmp_path)
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write("two")
    assert golden_update.open_epoch("the radius draw moved", root, environ={}) == 2
    assert _epochs(root)[1] == {"epoch": 2, "first_release": None, "pending": True, "reason": "the radius draw moved"}


@needs_git
def test_a_new_golden_file_counts_as_a_change(tmp_path):
    root = _repo(tmp_path)
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "nebula.txt"), "w") as f:
        f.write("new")
    assert golden_update.open_epoch("a new golden case", root, environ={}) == 2


@needs_git
def test_it_refuses_with_no_golden_change_under_ci_or_without_a_reason(tmp_path):
    root = _repo(tmp_path)
    with pytest.raises(golden_update.Refused, match="nothing to bump"):
        golden_update.require_changed_golden(root)
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write("two")
    with pytest.raises(golden_update.Refused, match="CI"):
        golden_update.open_epoch("why", root, environ={"CI": "true"})
    with pytest.raises(golden_update.Refused, match="reason"):
        golden_update.open_epoch("  ", root, environ={})
    assert len(_epochs(root)) == 1


@needs_git
def test_a_second_run_in_the_same_change_does_not_bump_twice(tmp_path):
    root = _repo(tmp_path)
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write("two")
    assert golden_update.open_epoch("first", root, environ={}) == 2
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write("three")
    assert golden_update.open_epoch("again", root, environ={}) == 2
    assert [entry["epoch"] for entry in _epochs(root)] == [1, 2]


@needs_git
def test_with_no_golden_file_yet_the_epoch_bumps_on_a_reason_alone(tmp_path):
    root = _repo(tmp_path)
    _git(root, "rm", "-q", "-f", os.path.join(golden_update.GOLDEN_DIR, "sector.txt"))
    _git(root, "commit", "-q", "-m", "no golden files")
    assert golden_update.open_epoch("the maths moved", root, environ={}) == 2
    assert _epochs(root)[-1]["pending"] is True


def test_rewriting_the_golden_files_is_refused_under_ci(tmp_path):
    with pytest.raises(golden_update.Refused, match="under CI"):
        golden_update.write_golden(str(tmp_path), environ={"CI": "true"})


def test_the_main_command_stops_when_the_rewrite_is_refused(monkeypatch, capsys):
    monkeypatch.setenv("CI", "true")
    assert golden_update.main(["--reason", "why", "--affects", "fill"]) == 1
    assert "under CI" in capsys.readouterr().err


def _with_stages(root):
    path = os.path.join(root, golden_update.EPOCHS_FILE)
    data = json.load(open(path))
    data["stages"] = ["plan", "scatter", "fill", "population"]
    data["epochs"][0]["affects"] = []
    json.dump(data, open(path, "w"), indent=2)
    _git(root, "commit", "-q", "-am", "stages")


@needs_git
def test_a_new_epoch_records_the_stages_it_changed(tmp_path):
    root = _repo(tmp_path)
    _with_stages(root)
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write("two")
    assert golden_update.open_epoch("the radius draw moved", root, environ={}, affects=["scatter", "fill"]) == 2
    assert _epochs(root)[1]["affects"] == ["scatter", "fill"]
    text = open(os.path.join(root, golden_update.EPOCHS_FILE)).read()
    assert '"affects": ["scatter", "fill"]' in text


@needs_git
def test_a_new_epoch_without_stages_or_with_an_unknown_one_is_refused(tmp_path):
    root = _repo(tmp_path)
    _with_stages(root)
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write("two")
    with pytest.raises(golden_update.Refused, match="--affects"):
        golden_update.open_epoch("why", root, environ={})
    with pytest.raises(golden_update.Refused, match="one or more of"):
        golden_update.open_epoch("why", root, environ={}, affects=["sky"])
    with pytest.raises(golden_update.Refused, match="--affects"):
        golden_update.check_affects_given(None, root)
    golden_update.check_affects_given(["fill"], root)


@needs_git
def test_a_second_run_can_change_the_stages_of_the_pending_epoch(tmp_path):
    root = _repo(tmp_path)
    _with_stages(root)
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write("two")
    golden_update.open_epoch("why", root, environ={}, affects=["fill"])
    assert golden_update.open_epoch("why", root, environ={}, affects=["scatter", "fill"]) == 2
    assert len(_epochs(root)) == 2 and _epochs(root)[1]["affects"] == ["scatter", "fill"]
    golden_update.check_affects_given(None, root)    # the pending epoch already has them


def _stage_a_rewrite(monkeypatch, root, content="two"):
    """Stands in for the ten-minute rewrite: it writes a golden file that says which epoch it was made under, the
    way a sector's structure does."""
    def rewrite(where=None):
        epoch = _epochs(root)[-1]["epoch"]
        with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
            f.write(f"{content} epoch {epoch}")
    monkeypatch.setattr(golden_update, "write_golden", rewrite)


@needs_git
def test_one_run_writes_every_golden_file_under_the_new_epoch(tmp_path, monkeypatch):
    """TEST.404: the epoch went up after the files were written, so the sector section of the first run carried the old
    epoch and test_golden_seed failed until a second run."""
    root = _repo(tmp_path)
    _stage_a_rewrite(monkeypatch, root)
    assert golden_update.main(["--reason", "the radius draw moved"], root=root) == 0
    assert open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt")).read() == "two epoch 2"
    assert [entry["epoch"] for entry in _epochs(root)] == [1, 2]


@needs_git
def test_a_run_that_changes_no_golden_file_leaves_the_epochs_as_they_were(tmp_path, monkeypatch, capsys):
    root = _repo(tmp_path)
    before = open(os.path.join(root, golden_update.EPOCHS_FILE)).read()
    monkeypatch.setattr(golden_update, "write_golden", lambda where=None: None)
    assert golden_update.main(["--reason", "why"], root=root) == 1
    assert "nothing to bump" in capsys.readouterr().err
    assert open(os.path.join(root, golden_update.EPOCHS_FILE)).read() == before


@needs_git
def test_a_second_run_rewrites_the_files_under_the_epoch_the_change_already_holds(tmp_path, monkeypatch):
    root = _repo(tmp_path)
    _stage_a_rewrite(monkeypatch, root)
    assert golden_update.main(["--reason", "why"], root=root) == 0
    _stage_a_rewrite(monkeypatch, root, content="three")
    assert golden_update.main(["--reason", "why"], root=root) == 0
    assert open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt")).read() == "three epoch 2"
    assert [entry["epoch"] for entry in _epochs(root)] == [1, 2]
