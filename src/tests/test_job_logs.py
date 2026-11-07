# tests/test_job_logs.py

"""
OPS.101: the nightly whole-map build and the nightly orbit update append each run to a log of their own in the log folder, `tiles.log` and
`orbits.log`, the way `planetgen warm` appends to `warmup.log` (`planetgen/util/joblog.py`); `update.log` is the same for `update.sh`.
"""

import re

import pytest

from planetgen.cli import orbits, tiles
from planetgen.util import joblog, log, logpaths

LINE = re.compile(r"^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d [A-Za-z0-9+\-:]+ (INFO|WARNING|ERROR)\s+\w+: ")


@pytest.fixture
def folder(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(tmp_path))
    monkeypatch.delenv("INVOCATION_ID", raising=False)
    log.configure(log.NORMAL)
    return tmp_path


def _lines(path):
    return path.read_text(encoding="utf-8").splitlines()


def test_every_job_has_its_own_log_in_the_log_folder_and_the_rule_for_log_files_covers_it(folder):
    names = {"warmup.log", "tiles.log", "orbits.log", "update.log", "error.log"}
    assert set(logpaths.JOB_LOG_NAMES) == names
    for name in names:
        assert logpaths.job_log_path(name) == str(folder / name) and name.endswith(".log")
    with pytest.raises(ValueError):
        logpaths.job_log_path("smoke.log")


def test_a_run_is_bracketed_by_what_started_it_and_how_it_ended(folder):
    with joblog.JobLog("tiles.log", "tiles", "Whole-map build", started_by="systemd timer"):
        log.normal("level 3: 120 nodes")
    lines = _lines(folder / "tiles.log")
    assert len(lines) == 3 and all(LINE.match(line) for line in lines)
    assert lines[0].endswith("tiles: Whole-map build started (systemd timer)")
    assert lines[1].endswith("tiles: level 3: 120 nodes")
    assert re.search(r"tiles: Whole-map build finished after \d+:\d\d$", lines[2])


@pytest.mark.parametrize("raises, ending", [
    (SystemExit(3), "ended with exit status 3 after"),
    (KeyboardInterrupt(), "ended by Ctrl+C after"),
    (RuntimeError("boom"), "failed after 0:00: boom"),
])
def test_a_run_that_does_not_finish_says_how_it_ended(folder, raises, ending):
    with pytest.raises(type(raises)):
        with joblog.JobLog("orbits.log", "orbits", "Orbital update"):
            raise raises
    assert ending in _lines(folder / "orbits.log")[-1]


def test_a_returned_exit_status_is_said_in_the_closing_line(folder):
    with joblog.JobLog("tiles.log", "tiles", "Whole-map build") as record:
        record.status = 75
    assert _lines(folder / "tiles.log")[-1].endswith("ended with exit status 75 after 0:00")


def test_the_tiles_command_leaves_its_run_in_tiles_log(folder, monkeypatch):
    monkeypatch.setattr(tiles, "_build", lambda args, config: log.normal("Built 3 nodes") or 0)
    assert tiles.main(["--mysql-database", "x"]) == 0
    lines = _lines(folder / "tiles.log")
    assert lines[0].endswith("tiles: Whole-map build started (by hand)")
    assert lines[1].endswith("tiles: Built 3 nodes") and "finished" in lines[-1]


def test_a_tiles_run_that_failed_says_so_in_its_log(folder, monkeypatch):
    monkeypatch.setattr(tiles, "_build", lambda args, config: 1)
    assert tiles.main(["--mysql-database", "x"]) == 1
    assert "ended with exit status 1" in _lines(folder / "tiles.log")[-1]


def test_the_timer_is_named_when_systemd_started_the_run(folder, monkeypatch):
    monkeypatch.setenv("INVOCATION_ID", "abc")
    monkeypatch.setattr(tiles, "_build", lambda args, config: 0)
    tiles.main(["--mysql-database", "x"])
    assert _lines(folder / "tiles.log")[0].endswith("started (systemd timer)")


def test_the_orbits_command_leaves_its_run_in_orbits_log(folder, monkeypatch):
    monkeypatch.setattr(orbits, "_update", lambda args, config: log.normal("0.000100 years elapsed"))
    monkeypatch.setattr("sys.argv", ["planetgen.cli.orbits", "--mysql-database", "x"])
    orbits.main()
    lines = _lines(folder / "orbits.log")
    assert lines[0].endswith("orbits: Orbital update started (by hand)")
    assert lines[1].endswith("orbits: 0.000100 years elapsed") and "finished" in lines[-1]


def test_an_orbit_update_that_exits_non_zero_says_so(folder, monkeypatch):
    def refuse(args, config):
        raise SystemExit(75)

    monkeypatch.setattr(orbits, "_update", refuse)
    monkeypatch.setattr("sys.argv", ["planetgen.cli.orbits", "--mysql-database", "x"])
    with pytest.raises(SystemExit):
        orbits.main()
    assert "ended with exit status 75" in _lines(folder / "orbits.log")[-1]


def test_a_log_that_cannot_be_written_is_said_and_the_run_goes_on(tmp_path, monkeypatch, capsys):
    blocker = tmp_path / "file"
    blocker.write_text("")
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(blocker))
    log.configure(log.NORMAL)
    with joblog.JobLog("tiles.log", "tiles", "Whole-map build"):
        log.normal("still running")
    out = capsys.readouterr().out
    assert "Cannot write the whole-map build log" in out and "still running" in out
