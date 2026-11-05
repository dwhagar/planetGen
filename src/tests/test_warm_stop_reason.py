# tests/test_warm_stop_reason.py

"""
A warm-up told to stop says what ended it in one line (OPS.95): `planetgen/cli/warm.py` and the reasons the cancellers leave
(`planetgen/admin/stop_reason.py`). No database or Redis is needed.
"""

import os
import signal
import time
from types import SimpleNamespace

import pytest

from planetgen.admin import stop_reason, updateguard
from planetgen.cli import warm as warm_cli
from planetgen.queue import stopall


@pytest.fixture
def state(tmp_path, monkeypatch):
    """A state folder of its own; the signal handlers the run installs are put back after the test."""
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    monkeypatch.delenv(warm_cli.TIME_LIMIT_ENV, raising=False)
    handlers = signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGHUP)
    yield tmp_path
    signal.signal(signal.SIGTERM, handlers[0])
    signal.signal(signal.SIGHUP, handlers[1])


def _run_ended_by(monkeypatch, number, before=lambda: None):
    """`run_warm` with a warm-up that receives signal `number` from outside; returns the exit status."""
    monkeypatch.setattr(warm_cli, "_site_database", lambda: "galaxy")

    def warm_once(_args, reason=None):
        before()
        os.kill(os.getpid(), number)
        time.sleep(5)       # never reached: the handler raises

    monkeypatch.setattr(warm_cli, "_warm_once", warm_once)
    with pytest.raises(SystemExit) as ended:
        warm_cli.run_warm(SimpleNamespace(on_change=None))
    return ended.value.code


def _line(capsys):
    seen = capsys.readouterr()
    return (seen.out + seen.err).strip()


def test_a_sigterm_ends_the_warm_up_with_one_line_and_a_status_not_a_traceback(state, monkeypatch, capsys):
    status = _run_ended_by(monkeypatch, signal.SIGTERM)
    line = _line(capsys)
    assert status == 143
    assert "Traceback" not in line and "KeyboardInterrupt" not in line and len(line.splitlines()) == 1
    assert "The cache warm-up was ended by a SIGTERM that left no reason" in line
    assert "What it stored is kept: run `planetgen warm --progress` again and it continues" in line


def test_the_reason_the_canceller_left_is_what_the_line_says(state, monkeypatch, capsys):
    status = _run_ended_by(monkeypatch, signal.SIGTERM, before=lambda: stop_reason.leave(os.getpid(), "stop-all ended every command"))
    assert status == 143 and "ended by stop-all ended every command." in _line(capsys)
    assert stop_reason.take(os.getpid()) is None        # read once


def test_a_reason_left_for_a_pid_used_before_is_not_read(state, monkeypatch, capsys):
    stop_reason.leave(os.getpid(), "an old cancel")
    _run_ended_by(monkeypatch, signal.SIGTERM)
    assert "an old cancel" not in _line(capsys)


def test_update_sh_time_limit_is_named_when_the_limit_it_told_has_passed(state, monkeypatch, capsys):
    monkeypatch.setenv(warm_cli.TIME_LIMIT_ENV, "3")
    monkeypatch.setattr(warm_cli.time, "monotonic", lambda: _clock[0])
    _clock = [100.0]
    monkeypatch.setattr(warm_cli, "_started", 100.0)
    status = _run_ended_by(monkeypatch, signal.SIGTERM, before=lambda: _clock.__setitem__(0, 104.0))
    assert status == 143
    assert "ended by update.sh's time limit for this stage (3 s, WARM_STAGE_SECONDS)" in _line(capsys)


def test_a_sigterm_long_before_the_limit_is_not_blamed_on_it(state, monkeypatch, capsys):
    monkeypatch.setenv(warm_cli.TIME_LIMIT_ENV, "14400")
    _run_ended_by(monkeypatch, signal.SIGTERM)
    assert "time limit" not in _line(capsys)


def test_a_hangup_says_the_terminal_was_closed(state, monkeypatch, capsys):
    status = _run_ended_by(monkeypatch, signal.SIGHUP)
    assert status == 129 and "ended by the terminal it ran in was closed (hangup)" in _line(capsys)


def test_ctrl_c_is_named_and_exits_130(state, monkeypatch, capsys):
    monkeypatch.setattr(warm_cli, "_site_database", lambda: "galaxy")
    monkeypatch.setattr(warm_cli, "_warm_once", lambda _args, reason=None: (_ for _ in ()).throw(KeyboardInterrupt()))
    with pytest.raises(SystemExit) as ended:
        warm_cli.run_warm(SimpleNamespace(on_change=None))
    assert ended.value.code == 130 and "ended by Ctrl+C." in _line(capsys)


def test_the_update_guards_cancel_leaves_its_reason_before_it_signals(state, monkeypatch):
    sent = []
    monkeypatch.setattr(updateguard.os, "kill", lambda pid, number: sent.append((pid, number)))
    updateguard.MarkerSource().cancel({"pid": 424242})
    assert sent == [(424242, signal.SIGTERM)]
    assert "update guard cancelled it" in stop_reason.take(424242)


def test_stop_all_leaves_its_reason_for_each_process(state, monkeypatch):
    rounds = iter([[777001, 777002], []])
    monkeypatch.setattr(stopall, "process_table", lambda: {})
    monkeypatch.setattr(stopall, "find_targets", lambda table, protected: next(rounds))
    monkeypatch.setattr(stopall, "_signal", lambda pids, number: [])
    monkeypatch.setattr(stopall, "_alive", lambda pid: False)
    stopall.kill_processes(term_seconds=0.0)
    assert "stop-all ended every planetgen worker" in stop_reason.take(777001)
    assert "stop-all ended every planetgen worker" in stop_reason.take(777002)


def test_a_folder_that_cannot_be_written_only_means_a_run_says_less(tmp_path, monkeypatch):
    blocker = tmp_path / "file"
    blocker.write_text("x")
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(blocker / "below"))
    assert stop_reason.leave(1234, "anything") is False
    assert stop_reason.take(1234) is None
