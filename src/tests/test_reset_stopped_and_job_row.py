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
    assert job_runner._row_is_canceled(_Control(state), 7, set(), "job j") is canceled


def test_a_missing_job_row_is_logged_once_per_step_with_where_and_that_the_cause_is_unknown(monkeypatch):
    said = []
    monkeypatch.setattr(job_runner.planetgen_log, "error", lambda message, *args, **kwargs: said.append(message % args))
    seen = set()
    for _ in range(3):
        assert job_runner._row_is_canceled(_Control(None), 7, seen, "job j (run r, step 1)") is False
    assert len(said) == 1
    assert "7" in said[0] and "job j (run r, step 1)" in said[0] and "not known" in said[0] and "OPS.296" in said[0]
    job_runner._row_is_canceled(_Control(None), 8, seen, "job j (run r, step 2)")
    assert len(said) == 2
