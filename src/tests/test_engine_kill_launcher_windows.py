"""OPS.246 (the launcher killed between two windows, from the OPS.240 note of Lane 16): a launcher that dies in the middle of `planrows.spawn_window` (before or after each commit of the ledger) leaves a stage that the next launcher carries on, with every unit row written once and every unit run once. Every kill point of a window is tried, one at a time, over the whole stage."""

import pytest

from planetgen.queue import planrows, runrecord
from tests.test_window_spawn import RAN, TOTAL, WINDOW, _rows, _spawn, _started, _units_of  # noqa: F401 -- the autouse fixture that cuts the stage into windows

from tests.test_window_spawn import _windowed  # noqa: F401

pytestmark = [pytest.mark.fresh_database]


class Killed(BaseException):
    """A launcher process that the kernel killed: nothing may catch it by accident."""


class _Counted:
    """The control connection of a window, whose commits are the points a kill lands on."""

    def __init__(self, real, plan):
        self._real, self._plan = real, plan

    def commit(self):
        number, self._plan["seen"] = self._plan["seen"], self._plan["seen"] + 1
        mode = self._plan["kills"].get(number)
        if mode == "before":
            self._plan["fired"].append(number)
            raise Killed(f"commit {number} before")
        self._real.commit()
        if mode == "after":
            self._plan["fired"].append(number)
            raise Killed(f"commit {number} after")

    def __getattr__(self, name):
        return getattr(self._real, name)


def _all_windows(conn, run_id, plan, monkeypatch, config):
    """A launcher's spawning, turn after turn: a turn that is killed ends, and the next turn is the next launcher's. The units of every window written are finished as the stage's own workers would finish them."""
    real = planrows.control_connection
    monkeypatch.setattr(planrows, "control_connection", lambda conn_, run: _Counted(real(conn_, run), plan))
    written = 0
    for _turn in range(4 * TOTAL):
        try:
            written += planrows.spawn_window(conn, run_id, 0)
        except Killed:
            pass
        names = [key for key, _state in _units_of(_rows(conn, run_id, config))]
        if len(names) >= TOTAL and _spawn(_rows(conn, run_id, config))["next"] >= TOTAL and _spawn_state(conn, run_id, config) == "done":
            break
        if names:
            conn.sadd(runrecord.done_key(run_id, 0), *names)
    return written


def _spawn_state(conn, run_id, config):
    return next(job for job in _rows(conn, run_id, config) if job["command"] == planrows.SPAWN)["state"]


def _spawn_without_kills(mysql_config, monkeypatch):
    conn, run_id = _started(mysql_config)
    plan = {"seen": 0, "kills": {}, "fired": []}
    _all_windows(conn, run_id, plan, monkeypatch, mysql_config)
    runrecord.finish_run(conn, run_id, state="failed", error="the test is over with it")
    return plan["seen"]


def test_every_commit_of_a_window_has_a_kill_that_leaves_the_stage_written_whole(mysql_config, monkeypatch):
    with monkeypatch.context() as patch:
        commits = _spawn_without_kills(mysql_config, patch)
    assert commits > 0
    for number in range(commits):
        for mode in ("before", "after"):
            with monkeypatch.context() as patch:
                conn, run_id = _started(mysql_config)
                plan = {"seen": 0, "kills": {number: mode}, "fired": []}
                _all_windows(conn, run_id, plan, patch, mysql_config)
                rows = _rows(conn, run_id, mysql_config)
                runrecord.finish_run(conn, run_id, state="failed", error="the test is over with it")
                assert plan["fired"] == [number], f"commit {number} {mode} was not reached"
                assert [key for key, _state in _units_of(rows)] == list(range(TOTAL)), f"commit {number} {mode}: every unit row once, no gap"
                assert _spawn(rows)["next"] == TOTAL, f"commit {number} {mode}: the cursor is at the end"
