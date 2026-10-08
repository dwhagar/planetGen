"""
ADM.24 and ADM.25: a failed run shows its error in full and, at a terminal,
waits for Enter; a worker's traceback reaches the run; a job runner's own
failure puts its traceback in the job's log.
"""

import io
import json
import os
import pickle
import sys

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.admin import activity_log
from planetgen.queue import work
from planetgen.web import job_runner


@pytest.fixture
def broken_handlers(monkeypatch):
    def explode(args):
        raise RuntimeError("the sector fill broke")

    for name in list(generate_cli._COMMAND_HANDLERS):
        monkeypatch.setitem(generate_cli._COMMAND_HANDLERS, name, explode)
    monkeypatch.setattr(activity_log, "event", lambda *a, **k: None)
    monkeypatch.setattr(sys, "argv", ["planetgen", "check-math"])


def test_an_unexpected_error_still_propagates_with_its_traceback(broken_handlers, monkeypatch):
    waited = []
    monkeypatch.setattr(generate_cli, "_wait_for_enter", lambda: waited.append(True))
    with pytest.raises(RuntimeError, match="the sector fill broke") as caught:
        generate_cli.main()
    assert caught.value.__traceback__ is not None  # Python prints it as the process exits
    assert waited == [True]


def test_a_database_error_prints_its_traceback_too(monkeypatch, capsys):
    import pymysql

    def refuse(args):
        raise pymysql.err.OperationalError(2003, "can't connect")

    for name in list(generate_cli._COMMAND_HANDLERS):
        monkeypatch.setitem(generate_cli._COMMAND_HANDLERS, name, refuse)
    monkeypatch.setattr(activity_log, "event", lambda *a, **k: None)
    monkeypatch.setattr(sys, "argv", ["planetgen", "check-math"])
    with pytest.raises(SystemExit):
        generate_cli.main()
    shown = capsys.readouterr()
    text = shown.out + shown.err
    assert "database error" in text and "Traceback (most recent call last)" in text


class _Terminal(io.StringIO):
    def isatty(self):
        return True


def test_a_failed_run_at_a_terminal_waits_for_enter(monkeypatch):
    asked = []
    monkeypatch.setattr(sys, "stdin", _Terminal())
    monkeypatch.setattr(sys, "stdout", _Terminal())
    monkeypatch.setattr("builtins.input", lambda prompt="": asked.append(prompt) or "")
    generate_cli._wait_for_enter()
    assert len(asked) == 1 and "Enter" in asked[0]


def test_a_failed_run_with_its_output_redirected_does_not_wait(monkeypatch):
    monkeypatch.setattr("builtins.input", lambda prompt="": pytest.fail("waited for input"))
    generate_cli._wait_for_enter()  # pytest's stdin/stdout are not terminals


def test_a_database_error_at_a_terminal_waits_after_showing_the_error(monkeypatch, capsys):
    import pymysql

    def refuse(args):
        raise pymysql.err.OperationalError(2003, "can't connect")

    for name in list(generate_cli._COMMAND_HANDLERS):
        monkeypatch.setitem(generate_cli._COMMAND_HANDLERS, name, refuse)
    monkeypatch.setattr(activity_log, "event", lambda *a, **k: None)
    monkeypatch.setattr(sys, "argv", ["planetgen", "check-math"])
    waited = []
    monkeypatch.setattr(generate_cli, "_wait_for_enter", lambda: waited.append(True))
    with pytest.raises(SystemExit):
        generate_cli.main()
    shown = capsys.readouterr()
    assert "Traceback (most recent call last)" in shown.out + shown.err
    assert waited == [True]


def test_a_workers_traceback_rides_along_with_its_exception():
    try:
        raise ValueError("bad payload")
    except ValueError as caught:
        exc = work._picklable_error(caught)
    again = pickle.loads(pickle.dumps(exc))
    assert isinstance(again, ValueError)
    assert any("Traceback in the worker" in note and "bad payload" in note for note in again.__notes__)


def test_a_worker_exception_that_cannot_be_pickled_keeps_its_traceback_text():
    class Unpicklable(Exception):
        def __init__(self):
            super().__init__("local")
            self.handle = lambda: None

    try:
        raise Unpicklable()
    except Unpicklable as caught:
        exc = work._picklable_error(caught)
    assert isinstance(exc, RuntimeError)
    assert "Unpicklable: local" in str(exc) and "Traceback in the worker" in str(exc)


def test_the_job_runner_logs_its_own_traceback(tmp_path):
    try:
        raise RuntimeError("runner broke")
    except RuntimeError:
        job_runner._log_traceback(str(tmp_path))
    text = (tmp_path / "output.log").read_text(encoding="utf-8")
    assert "Traceback (most recent call last)" in text and "RuntimeError: runner broke" in text
