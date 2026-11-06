# tests/test_warm_log.py

"""
`warmup.log` (OPS.96): `planetgen warm` appends each run to it itself, in local time with its zone, bracketed by what started it and
how it ended (`planetgen/cli/warm_log.py`). No database or Redis is needed.
"""

import re
from types import SimpleNamespace

import pytest

from planetgen.cli import warm as warm_cli
from planetgen.cli import warm_log, warm_progress
from planetgen.util import log, logpaths

LINE = re.compile(r"^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d [A-Za-z0-9+\-:]+ (INFO|WARNING|ERROR)\s+warm: ")


@pytest.fixture
def path(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(tmp_path))
    for name in (warm_log.STARTED_BY_ENV, warm_log.TIME_LIMIT_ENV):
        monkeypatch.delenv(name, raising=False)
    log.configure(log.NORMAL)
    return tmp_path / "warmup.log"


def _lines(path):
    return path.read_text(encoding="utf-8").splitlines()


def test_the_warm_up_log_is_warmup_log_in_the_log_folder(tmp_path, monkeypatch):
    monkeypatch.delenv("PLANETGEN_LOG_DIR", raising=False)
    assert logpaths.warmup_log_path({"log_dir": str(tmp_path)}) == str(tmp_path / "warmup.log")
    assert logpaths.warmup_log_path({}).endswith("/warmup.log")


def test_a_run_is_bracketed_and_every_line_has_the_local_time_and_zone(path):
    with warm_log.RunLog() as record:
        log.normal("galaxy-map: 1 warmed (0 already cached, 0 left out by the cap) in 4.0 s")
        record.say("a line only the bars' console printed")
        record.families = 5
    lines = _lines(path)
    assert len(lines) == 4 and all(LINE.match(line) for line in lines)
    assert lines[0].endswith("Cache warm-up started (by hand)")
    assert lines[1].endswith("warm: galaxy-map: 1 warmed (0 already cached, 0 left out by the cap) in 4.0 s")
    assert lines[2].endswith("a line only the bars' console printed")
    assert re.search(r"Cache warm-up finished after \d+:\d\d \(5 families\)$", lines[3])


def test_runs_are_appended_not_overwritten(path):
    for _ in range(2):
        with warm_log.RunLog():
            pass
    assert [line.count("started") for line in _lines(path)].count(1) == 2 and len(_lines(path)) == 4


def test_the_opening_line_says_what_started_the_run(path, monkeypatch):
    environ = {warm_log.STARTED_BY_ENV: "update.sh, part 2 of 2, the full warm-up", warm_log.TIME_LIMIT_ENV: "14400"}
    assert warm_log.started_by(environ=environ) == "update.sh, part 2 of 2, the full warm-up; time limit 4 hours"
    assert warm_log.started_by(environ={}) == "by hand"
    assert warm_log.started_by("orbits, systems", environ=environ) == "after a change: orbits, systems"
    assert warm_log.started_by(environ={warm_log.TIME_LIMIT_ENV: "5400"}) == "by hand; time limit 90 minutes"
    monkeypatch.setenv(warm_log.STARTED_BY_ENV, "update.sh, part 1 of 2")
    with warm_log.RunLog():
        pass
    assert _lines(path)[0].endswith("Cache warm-up started (update.sh, part 1 of 2)")


def test_a_failed_run_says_why_and_a_stopped_one_says_by_what(path):
    with pytest.raises(warm_cli.Failed):
        with warm_log.RunLog():
            raise warm_cli.Failed("no Redis server answers")
    with pytest.raises(warm_cli.Stopped):
        with warm_log.RunLog():
            raise warm_cli.Stopped("a SIGTERM that left no reason", 143)
    with pytest.raises(KeyboardInterrupt):
        with warm_log.RunLog():
            raise KeyboardInterrupt()
    ends = [line for line in _lines(path) if "started" not in line]
    assert re.search(r"ERROR\s+warm: Cache warm-up failed after \d+:\d\d: no Redis server answers$", ends[0])
    assert re.search(r"ended by a SIGTERM that left no reason after \d+:\d\d$", ends[1])
    assert re.search(r"ended by Ctrl\+C after \d+:\d\d$", ends[2])


def test_lines_printed_after_the_run_are_not_in_the_file(path):
    with warm_log.RunLog():
        pass
    log.normal("after the run")
    assert not any("after the run" in line for line in _lines(path))


def test_a_file_that_cannot_be_written_is_said_once_and_the_run_goes_on(tmp_path, capsys):
    blocked = tmp_path / "file"
    blocked.write_text("not a folder")
    with warm_log.RunLog(path=str(blocked / "warmup.log")) as record:
        record.say("still runs")
    out = capsys.readouterr()
    assert "Cannot write the warm-up log" in out.out + out.err


def test_the_bars_console_lines_reach_the_file_through_the_mirror():
    mirrored = []
    bars = warm_progress.WarmProgress(["lists"], lambda line: None, show=True, rich_bars=True, mirror=mirrored.append)
    with bars:
        bars.plan(0, {"units": 1, "cost": None})
        bars.result({"warmed": 3, "built": 3, "cached": 0, "left_out": 0, "noun": "pages", "seconds": 1.0,
                     "line": "lists: 3 warmed (0 already cached, 0 left out by the cap) in 1.0 s"})
    assert mirrored == ["lists: 3 warmed (0 already cached, 0 left out by the cap) in 1.0 s"]


def test_warm_once_keeps_a_record_of_the_run_and_the_reason_a_change_gave(path, monkeypatch):
    def run(args, reason, record):
        log.normal("one family warmed")
        record.families = 1

    monkeypatch.setattr(warm_cli, "_warm_run", run)
    warm_cli._warm_once(SimpleNamespace(), reason="orbits")
    lines = _lines(path)
    assert lines[0].endswith("Cache warm-up started (after a change: orbits)") and lines[1].endswith("one family warmed")
    assert lines[2].endswith("(1 family)") and "finished" in lines[2]
