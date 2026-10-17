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
    assert golden_update.bump("the radius draw moved", root, environ={}) == 2
    assert _epochs(root)[1] == {"epoch": 2, "first_release": None, "pending": True, "reason": "the radius draw moved"}


@needs_git
def test_a_new_golden_file_counts_as_a_change(tmp_path):
    root = _repo(tmp_path)
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "nebula.txt"), "w") as f:
        f.write("new")
    assert golden_update.bump("a new golden case", root, environ={}) == 2


@needs_git
def test_it_refuses_with_no_golden_change_under_ci_or_without_a_reason(tmp_path):
    root = _repo(tmp_path)
    with pytest.raises(golden_update.Refused, match="nothing to bump"):
        golden_update.bump("why", root, environ={})
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write("two")
    with pytest.raises(golden_update.Refused, match="CI"):
        golden_update.bump("why", root, environ={"CI": "true"})
    with pytest.raises(golden_update.Refused, match="reason"):
        golden_update.bump("  ", root, environ={})
    assert len(_epochs(root)) == 1


@needs_git
def test_a_second_run_in_the_same_change_does_not_bump_twice(tmp_path):
    root = _repo(tmp_path)
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write("two")
    assert golden_update.bump("first", root, environ={}) == 2
    with open(os.path.join(root, golden_update.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write("three")
    assert golden_update.bump("again", root, environ={}) == 2
    assert [entry["epoch"] for entry in _epochs(root)] == [1, 2]


@needs_git
def test_with_no_golden_file_yet_the_epoch_bumps_on_a_reason_alone(tmp_path):
    root = _repo(tmp_path)
    _git(root, "rm", "-q", "-f", os.path.join(golden_update.GOLDEN_DIR, "sector.txt"))
    _git(root, "commit", "-q", "-m", "no golden files")
    assert golden_update.bump("the maths moved", root, environ={}) == 2
    assert _epochs(root)[-1]["pending"] is True
