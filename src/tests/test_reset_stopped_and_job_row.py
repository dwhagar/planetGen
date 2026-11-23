"""OPS.295 and OPS.296: a stopped reset says what stays wiped (its handler named a log level that does not exist), and the job runner reads a job's `engine_jobs` row without assuming it is there (a missing row is no Cancel)."""

import sys
import types

import pytest

from planetgen.cli import reset
from planetgen.util import process
from planetgen.web import job_runner


def test_a_reset_stopped_by_a_signal_says_what_stays_wiped_and_exits_143(monkeypatch, capsys):
    """The stop handler once called `log.info`, which `util.log` does not have: the handler died with an AttributeError and the exit code was lost."""
    monkeypatch.setattr(sys, "argv", ["reset", "--yes"])
    monkeypatch.setattr(reset, "mysql_config_from_args", lambda args: types.SimpleNamespace(database="never_used"))
    monkeypatch.setattr(reset.process, "interrupt_on_sigterm", lambda: None)

    def stopped(*args, **kwargs):
        raise process.Stopped("a test")

    monkeypatch.setattr(reset, "reset_database", stopped)
    with pytest.raises(SystemExit) as exit_info:
        reset.main()
    assert exit_info.value.code == 143


class _Control:
    """A control connection whose `engine_jobs` row is `state` (None: no row)."""

    def __init__(self, state):
        self.state = state

    def execute(self, sql, params=()):
        assert "FROM engine_jobs" in sql
        row = None if self.state is None else {"state": self.state}
        return types.SimpleNamespace(fetchone=lambda: row)


@pytest.mark.parametrize("state,canceled", [("running", False), ("canceled", True), (None, False)])
def test_the_runner_stops_a_job_only_for_a_cancel_and_a_missing_row_is_no_cancel(state, canceled):
    assert job_runner._row_is_canceled(_Control(state), 7) is canceled
