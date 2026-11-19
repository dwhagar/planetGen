"""OPS.190, checked against the new engine: a unit whose process is alive but silent (stopped, as a swapped-out or frozen one is) is not handed out a second time, so a retry never overlaps an old attempt; when the old process goes on it finishes the unit once. The late commit of an attempt that was replaced is refused by the fence (`test_attempts.py`: the replaced attempt's commit rolls back, the other lands, a replaced unit counts nothing and ends replaced). Real processes (`tests/test_engine_kill_processes.py`)."""

import os
import signal

import pytest

from tests.test_engine_kill_processes import _attempts, _each_piece_finished_once, _launcher, _state, _until, _units, engine  # noqa: F401 (the fixture)

pytestmark = [pytest.mark.fresh_database]


def test_a_job_whose_process_is_alive_but_silent_is_not_handed_out_again_and_finishes_once_when_it_goes_on(engine):  # noqa: F811
    _run, (job,), (marks,) = _units(engine, pieces=6, piece_seconds=0.3)
    one = _launcher(engine)
    one.turn()
    host = engine.forker.hosts["plain"]
    _until(lambda: os.path.exists(marks), 20, "the first piece")
    process = next(iter(host.jobs.values()))
    os.kill(process.pid, signal.SIGSTOP)
    try:
        for _ in range(12):
            one.recover()
            one.turn()
            assert _state(engine, job)["state"] == "running" and _state(engine, job)["attempt"] == 1, "no second attempt while the first process is alive"
    finally:
        os.kill(process.pid, signal.SIGCONT)
    _until(lambda: _state(engine, job)["state"] == "done", 60, "the job to finish once it goes on", step=one.turn)
    assert _attempts(engine, job) == [("done", None)]
    _each_piece_finished_once(marks, 6)
