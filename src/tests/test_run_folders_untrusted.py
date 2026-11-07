# tests/test_run_folders_untrusted.py

"""
SEC.38: the `running/` and `stop-reasons/` folders of the state folder are open to every user (mode 1777) and the guard that reads them
(`update.sh` runs it as root) must not believe what a local user leaves there: a forged marker must not make the guard signal somebody
else's process, and a link planted under a pid's name must not make a write (as root) overwrite the file it points at. No database or
Redis is needed.
"""

import json
import os
import subprocess
import sys

import pytest

from planetgen.admin import stop_reason, updateguard


@pytest.fixture
def state(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    return tmp_path


@pytest.fixture
def sleeper():
    """A process of this user that lives for the test."""
    process = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    yield process
    process.kill()
    process.wait()


def _marker(kind, pid, folder=None):
    folder = folder or updateguard.running_dir()
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, f"{kind}-{pid}.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump({"kind": kind, "title": "t", "pid": pid, "since": "2026-11-07T00:00:00Z"}, handle)
    return path


def test_a_marker_the_process_of_its_pid_did_not_leave_is_not_a_run(state, sleeper, monkeypatch):
    path = _marker("orbits", sleeper.pid)
    source = updateguard.MarkerSource()
    assert [item["pid"] for item in source.items()] == [sleeper.pid]            # the owner of the marker is the owner of the process
    # the process runs as somebody else than the user who wrote the marker: the marker is forged, whoever it names
    monkeypatch.setattr(updateguard, "process_owner", lambda pid: os.getuid() + 1)
    assert source.items() == []
    assert os.path.exists(path) and sleeper.poll() is None                      # nothing was signalled, nothing removed


def test_a_cancel_never_reaches_a_process_a_forged_marker_names(state, sleeper, monkeypatch):
    _marker("orbits", sleeper.pid)
    monkeypatch.setattr(updateguard, "process_owner", lambda pid: os.getuid() + 1)
    sources = [updateguard.MarkerSource()]
    assert updateguard.list_running(sources) == []
    updateguard.cancel_all(sources, updateguard.list_running(sources), write=lambda *_: None)
    assert sleeper.poll() is None


def test_a_marker_that_does_not_name_its_own_pid_is_not_a_run(state, sleeper):
    folder = updateguard.running_dir()
    other = _marker("orbits", sleeper.pid, folder)
    os.rename(other, os.path.join(folder, "orbits-1.json"))                      # file name and pid inside disagree
    assert updateguard.MarkerSource().items() == []


def test_a_link_under_a_marker_name_is_not_a_run_and_is_not_followed(state, sleeper, tmp_path):
    real = _marker("orbits", sleeper.pid, str(tmp_path / "elsewhere"))
    folder = updateguard.running_dir()
    os.makedirs(folder, exist_ok=True)
    os.symlink(real, os.path.join(folder, f"orbits-{sleeper.pid}.json"))
    assert updateguard.MarkerSource().items() == []


def test_a_pipe_under_a_marker_name_does_not_block_the_guard(state):
    folder = updateguard.running_dir()
    os.makedirs(folder, exist_ok=True)
    os.mkfifo(os.path.join(folder, "warm-map-1.json"))
    assert updateguard.MarkerSource().items() == []


def test_a_run_registers_over_a_planted_link_without_writing_through_it(state, tmp_path):
    victim = tmp_path / "victim.txt"
    victim.write_text("keep me\n", encoding="utf-8")
    folder = updateguard.running_dir()
    os.makedirs(folder, exist_ok=True)
    marker = os.path.join(folder, f"orbits-{os.getpid()}.json")
    os.symlink(victim, marker)
    with updateguard.running("orbits", "Orbital update"):
        assert not os.path.islink(marker) and json.load(open(marker, encoding="utf-8"))["pid"] == os.getpid()
        assert [item["pid"] for item in updateguard.MarkerSource().items()] == [os.getpid()]
    assert victim.read_text(encoding="utf-8") == "keep me\n"
    assert not os.path.lexists(marker)


def test_a_reason_is_left_over_a_planted_link_without_writing_through_it(state, tmp_path):
    victim = tmp_path / "victim.txt"
    victim.write_text("keep me\n", encoding="utf-8")
    os.makedirs(stop_reason.folder())
    os.symlink(victim, os.path.join(stop_reason.folder(), "4242.txt"))
    assert stop_reason.leave(4242, "an update starts")
    assert victim.read_text(encoding="utf-8") == "keep me\n"
    path = os.path.join(stop_reason.folder(), "4242.txt")
    assert not os.path.islink(path) and open(path, encoding="utf-8").read() == "an update starts\n"
    assert os.stat(path).st_mode & 0o777 == 0o644           # the run, another user, reads it


def test_a_reason_is_left_over_a_file_somebody_else_made(state):
    os.makedirs(stop_reason.folder())
    path = os.path.join(stop_reason.folder(), "4243.txt")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("old\n")
    assert stop_reason.leave(4243, "new reason")
    assert stop_reason.take(4243) == "new reason"
    assert stop_reason.take(4243) is None


def test_a_run_reads_no_reason_through_a_link(state, tmp_path):
    secret = tmp_path / "secret.txt"
    secret.write_text("not a reason\n", encoding="utf-8")
    os.makedirs(stop_reason.folder())
    os.symlink(secret, os.path.join(stop_reason.folder(), "4244.txt"))
    assert stop_reason.take(4244) is None
    assert secret.read_text(encoding="utf-8") == "not a reason\n"


@pytest.mark.skipif(os.geteuid() != 0, reason="giving a file to another user needs root")
def test_a_run_reads_no_reason_a_third_user_wrote(state):
    os.makedirs(stop_reason.folder())
    path = os.path.join(stop_reason.folder(), "4245.txt")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("forged\n")
    os.chown(path, 12345, 12345)
    assert stop_reason.take(4245) is None
    assert not os.path.exists(path)
