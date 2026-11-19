"""OPS.246 (the transition missing from the table of `test_engine_kill_chain`): a run canceled while its job is killed, in either order, ends with the unit canceled and handed out no more. The pieces it had done stay done, none is done twice, and no later attempt is made (real processes, `tests/test_engine_kill_processes.py`)."""

import os
import signal

import pytest

from planetgen.queue import controls
from tests.test_engine_kill_processes import _attempts, _launcher, _pieces, _state, _until, _units, engine  # noqa: F401 (the fixture)

pytestmark = [pytest.mark.fresh_database]


@pytest.mark.parametrize("order", ["canceled first", "killed first"])
def test_a_run_canceled_while_its_job_is_killed_leaves_the_unit_canceled_and_never_handed_out_again(engine, order):  # noqa: F811
    run, (job,), (marks,) = _units(engine, pieces=10, piece_seconds=0.4)
    one = _launcher(engine)
    one.turn()
    host = engine.forker.hosts["plain"]
    _until(lambda: os.path.exists(marks), 20, "the first piece")
    process = next(iter(host.jobs.values()))

    def cancel():
        ended = controls.end_rows(engine.control, run)
        engine.control.commit()
        assert ended["rows"] >= 1 and ended["runs"] == 1

    if order == "canceled first":
        cancel()
        os.kill(process.pid, signal.SIGKILL)
    else:
        os.kill(process.pid, signal.SIGKILL)
        cancel()
    _until(lambda: job not in one.running, 30, "the launcher to settle the end of the killed job", step=one.turn)
    for _ in range(5):
        one.recover()
        one.turn()
    assert _state(engine, job)["state"] == "canceled", "the unit is canceled, not ready to be handed out"
    assert len(_attempts(engine, job)) == 1 and _attempts(engine, job)[0][0] is not None, f"one attempt, closed: no later attempt was made: {_attempts(engine, job)}"
    starts, dones = _pieces(marks)
    assert all(count == 1 for count in dones.values()) and len(dones) < 10, f"no piece was done twice and the unit did not run on: {dones}"
