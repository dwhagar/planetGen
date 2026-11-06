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
        warm_cli.run_warm(SimpleNamespace(on_change=None, follow=False))
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


def test_a_hangup_does_not_end_the_warm_up_it_only_lets_go_of_the_terminal(state, monkeypatch):
    """PERF.235: a lost terminal silences the process's output, says so in the run's log, and the warm-up goes on to its end."""
    from planetgen.util import log
    silenced, said, finished = [], [], []
    monkeypatch.setattr(warm_cli, "_silence", lambda fds=(1, 2): silenced.append(fds))
    monkeypatch.setattr(log, "normal", said.append)
    monkeypatch.setattr(warm_cli, "_site_database", lambda: "galaxy")

    def warm_once(_args, reason=None):
        os.kill(os.getpid(), signal.SIGHUP)
        os.kill(os.getpid(), signal.SIGHUP)     # a second hangup changes nothing
        finished.append(True)

    monkeypatch.setattr(warm_cli, "_warm_once", warm_once)
    warm_cli.run_warm(SimpleNamespace(on_change=None, follow=False))
    assert finished == [True], "the warm-up ran to its end"
    assert silenced == [(1, 2)], "standard output and error were pointed away from the closed terminal, once"
    assert len(said) == 1 and "terminal this warm-up ran in was closed" in said[0] and "planetgen warm --follow" in said[0]


def test_silence_points_the_descriptors_at_dev_null(tmp_path):
    target = os.open(str(tmp_path / "out"), os.O_WRONLY | os.O_CREAT)
    try:
        saved = os.dup(target)
        warm_cli._silence((target,))
        os.write(target, b"gone")
        os.dup2(saved, target)
        os.close(saved)
    finally:
        os.close(target)
    assert (tmp_path / "out").read_bytes() == b""


def test_ctrl_c_is_named_and_exits_130(state, monkeypatch, capsys):
    monkeypatch.setattr(warm_cli, "_site_database", lambda: "galaxy")
    monkeypatch.setattr(warm_cli, "_warm_once", lambda _args, reason=None: (_ for _ in ()).throw(KeyboardInterrupt()))
    with pytest.raises(SystemExit) as ended:
        warm_cli.run_warm(SimpleNamespace(on_change=None, follow=False))
    assert ended.value.code == 130 and "ended by Ctrl+C." in _line(capsys)


def test_the_update_guards_cancel_leaves_its_reason_before_it_signals(state, monkeypatch):
    sent = []
    monkeypatch.setattr(updateguard.os, "kill", lambda pid, number: sent.append((pid, number)))
    updateguard.MarkerSource().cancel({"pid": 424242})
    assert sent == [(424242, signal.SIGTERM)]
    assert "update guard canceled it" in stop_reason.take(424242)


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
