"""OPS.228, slice A: the launcher (`planetgen/queue/launcher.py`) hands out the ready rows of the ledger by the dispatch rule, writes each attempt before it forks, and closes each attempt with how its process ended.
The templates are stood in for by a fake forker; the rows are real."""

import os
import signal

import pytest

from planetgen.db import store
from planetgen.queue import jobhost, launcher, ledger, memorygate

HOST = memorygate.Host(16000.0, 12000.0, 8, 0.2)


class FakeJob:
    def __init__(self, job_id, pid):
        self.id, self.pid, self.said, self.closed = job_id, pid, [], False

    def messages(self):
        said, self.said = self.said, []
        return said

    def close(self):
        self.closed = True


class FakeForker:
    def __init__(self):
        self.forked, self.exits, self.made = [], [], 0

    def fork(self, kind, runner, job=None):
        self.made += 1
        held = FakeJob(self.made, 1000 + self.made)
        self.forked.append((kind, runner, job, held))
        return "plain", held

    def poll(self, timeout=0.0):
        found, self.exits = self.exits, []
        return found


@pytest.fixture
def control(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.execute("DELETE FROM engine_runs")
    conn.commit()
    try:
        yield conn
    finally:
        conn.execute("DELETE FROM engine_runs")
        conn.commit()
        conn.close()


def _ready(control, count=3, mb=900):
    run = ledger.add_run(control, "generate", "g1")
    ledger.add_jobs(control, run, [{"stage": "fill", "command": "open", "arguments": {"index": 4, "plan": {}, "run_id": "20261115-x"}}]
                    + [{"stage": "fill", "command": "fill-unit", "item_key": str(n), "arguments": {"args": {"unit": n}, "cost": 1.5}, "claim_mb": mb} for n in range(count)])
    control.execute("UPDATE engine_jobs SET state = 'ready' WHERE run_id = ? AND item_key <> ''", (run,))
    control.commit()
    return run


def _state(control, run):
    return [(row["item_key"], row["state"], row["attempt"]) for row in control.execute("SELECT item_key, state, attempt FROM engine_jobs WHERE run_id = ? AND item_key <> '' ORDER BY id", (run,)).fetchall()]


def _launcher(control, forker, host=HOST):
    return launcher.Launcher(control, forker, read_host=lambda: host)


def test_every_ready_job_that_fits_is_handed_out_as_the_unit_job_forkedjob_takes(control):
    run = _ready(control)
    forker = FakeForker()
    started = _launcher(control, forker).turn()
    assert len(started) == 3
    assert [row[1] for row in _state(control, run)] == ["running"] * 3
    kind, runner, job, _held = forker.forked[0]
    assert (kind, runner) == ("fill-unit", launcher.RUNNER)
    assert job == {"run_id": "20261115-x", "stage": 4, "unit": 0, "kind": "fill-unit", "args": {"unit": 0}, "cost": 1.5, "row": started[0], "retries_left": 1}


def test_a_job_that_does_not_fit_waits_with_its_reason_on_the_row(control):
    run = _ready(control, count=3, mb=7000)
    forker = FakeForker()
    one = _launcher(control, forker)
    started = one.turn()
    assert 0 < len(started) < 3
    waiting = control.execute("SELECT waited FROM engine_jobs WHERE run_id = ? AND state = 'ready'", (run,)).fetchall()
    assert waiting and all(row["waited"] for row in waiting)


def test_a_unit_that_wrote_its_own_end_is_left_as_it_wrote_it(control):
    run = _ready(control, count=2)
    forker = FakeForker()
    one = _launcher(control, forker)
    first, second = one.turn()
    ledger.begin_attempt(control, first)
    ledger.end_attempt(control, first, 1, "done")
    ledger.set_job_state(control, first, "done")
    control.commit()
    forker.exits = [("plain", jobhost.Exit(1, 1001, code=0))]
    one.turn()
    assert [row[1] for row in _state(control, run)] == ["done", "running"]
    assert ledger.attempt_row(control, first, 1)["outcome"] == "done" and list(one.running) == [second]


def test_a_process_killed_before_it_ended_its_unit_has_the_unit_handed_out_once_more_then_failed(control):
    run = _ready(control, count=1)
    forker = FakeForker()
    one = _launcher(control, forker)
    (job,) = one.turn()
    ledger.begin_attempt(control, job)                                  # the unit took its attempt, then the process was killed
    forker.exits = [("plain", jobhost.Exit(1, 1001, signal=signal.SIGKILL))]
    started_again = one.turn()
    assert ledger.attempt_row(control, job, 1)["outcome"] == "lost" and "SIGKILL" in ledger.attempt_row(control, job, 1)["error_text"]
    assert started_again == [job] and forker.forked[1][2]["retries_left"] == 0, "the second hand-out has no retry left"
    ledger.begin_attempt(control, job)
    forker.exits = [("plain", jobhost.Exit(2, 1002, signal=signal.SIGKILL))]
    assert one.turn() == []
    assert _state(control, run) == [("0", "failed", 2)] and one.running == {}


def test_a_failure_the_job_said_is_the_reason_on_a_row_it_left_running(control):
    run = _ready(control, count=1)
    forker = FakeForker()
    one = _launcher(control, forker)
    (job,) = one.turn()
    ledger.begin_attempt(control, job)
    ledger.begin_attempt(control, job)                                  # a retry already used
    forker.forked[0][3].said.append({"event": "failed", "error": "ValueError: no star"})
    forker.exits = [("plain", jobhost.Exit(1, 1001, code=1))]
    one.turn()
    assert _state(control, run) == [("0", "failed", 2)]
    assert "ValueError: no star" in ledger.attempt_row(control, job, 2)["error_text"]


def test_a_job_that_is_running_is_not_handed_out_again(control):
    run = _ready(control, count=2)
    forker = FakeForker()
    one = _launcher(control, forker)
    one.turn()
    assert one.turn() == [] and len(forker.forked) == 2


def test_a_unit_that_found_its_run_gone_and_skipped_is_canceled_not_handed_out_again(control):
    run = _ready(control, count=1)
    forker = FakeForker()
    one = _launcher(control, forker)
    (job,) = one.turn()
    forker.forked[0][3].said = [{"event": "unit", "result": {"skipped": True}}]
    forker.exits = [("plain", jobhost.Exit(1, 1001, code=0))]
    assert one.turn() == []
    assert _state(control, run) == [("0", "canceled", 0)] and one.running == {}


def test_rows_left_running_by_a_supervisor_that_ended_are_handed_out_again_once_their_process_is_gone(control):
    run = _ready(control, count=2)
    control.execute("UPDATE engine_jobs SET state = 'running' WHERE run_id = ? AND item_key <> ''", (run,))
    control.commit()
    first, second = [row["id"] for row in control.execute("SELECT id FROM engine_jobs WHERE run_id = ? AND item_key <> '' ORDER BY id", (run,)).fetchall()]
    ledger.begin_attempt(control, second, pid=os.getpid())                  # alive: its own to end
    control.commit()
    forker = FakeForker()
    one = _launcher(control, forker)
    assert one.recover() == [first], "no attempt was ever taken: the supervisor ended between marking the row and forking it"
    assert [row[1] for row in _state(control, run)] == ["ready", "running"]
    assert one.recover() == [], "looked for once in RECOVER_SECONDS"


def test_a_job_whose_template_is_still_being_built_stays_ready_and_is_handed_out_when_it_is_up(control):
    run = _ready(control, count=1)

    class Warming(FakeForker):
        building = True

        def fork(self, kind, runner, job=None):
            if self.building:
                raise jobhost.NotReady("the template 'plain' is being built")
            return super().fork(kind, runner, job)

    forker = Warming()
    one = _launcher(control, forker)
    assert one.turn() == [] and _state(control, run) == [("0", "ready", 0)]
    forker.building = False
    assert len(one.turn()) == 1 and _state(control, run) == [("0", "running", 0)]


def _opening(control):
    run = ledger.add_run(control, "generate", "g1")
    ledger.add_jobs(control, run, [{"stage": "fill", "command": "opening", "arguments": {"run_id": "20261115-x", "index": 4}, "claim_mb": 1500}])
    control.execute("UPDATE engine_jobs SET state = 'ready' WHERE run_id = ?", (run,))
    control.commit()
    return run


def _opening_state(control, run):
    row = control.execute("SELECT state, attempt, waited FROM engine_jobs WHERE run_id = ?", (run,)).fetchone()
    return row["state"], row["attempt"], row["waited"]


def test_the_opening_of_a_stage_is_handed_out_as_an_opening_job_with_an_attempt_the_launcher_writes(control):
    run = _opening(control)
    forker = FakeForker()
    one = _launcher(control, forker)
    started = one.turn()
    assert len(started) == 1
    kind, _runner, job, held = forker.forked[0]
    assert kind == "opening" and job == {"opening": True, "run_id": "20261115-x", "stage": 4, "kind": "opening", "row": started[0], "retries_left": 1}
    assert _opening_state(control, run) == ("running", 1, None)
    attempt = control.execute("SELECT pid, outcome FROM engine_attempts WHERE job_id = ?", (started[0],)).fetchone()
    assert (attempt["pid"], attempt["outcome"]) == (held.pid, "running")
    held.said = [{"event": "unit", "result": {"opened": True}}]
    forker.exits = [("plain", jobhost.Exit(held.id, held.pid, 0, None))]
    one.turn()
    assert _opening_state(control, run)[0] == "done"
    assert control.execute("SELECT outcome FROM engine_attempts WHERE job_id = ?", (started[0],)).fetchone()["outcome"] == "done"


def test_an_opening_of_a_run_that_is_no_longer_live_is_canceled_and_its_attempt_closed(control):
    run = _opening(control)
    forker = FakeForker()
    one = _launcher(control, forker)
    started = one.turn()
    held = forker.forked[0][3]
    held.said = [{"event": "unit", "result": {"skipped": True}}]
    forker.exits = [("plain", jobhost.Exit(held.id, held.pid, 0, None))]
    one.turn()
    assert _opening_state(control, run)[0] == "canceled"
    assert control.execute("SELECT outcome FROM engine_attempts WHERE job_id = ?", (started[0],)).fetchone()["outcome"] == "done"


def test_a_stage_runs_no_more_units_at_once_than_its_run_was_asked_to_use(control):
    run = ledger.add_run(control, "generate", "g1")
    ledger.add_jobs(control, run, [{"stage": "fill", "command": "open", "arguments": {"index": 0, "plan": {}, "run_id": "20261115-x", "cap": 2}}]
                    + [{"stage": "fill", "command": "fill-unit", "item_key": str(n), "arguments": {"args": {}, "cost": 1.0}, "claim_mb": 500} for n in range(5)])
    control.execute("UPDATE engine_jobs SET state = 'done' WHERE run_id = ? AND item_key = ''", (run,))
    control.execute("UPDATE engine_jobs SET state = 'ready' WHERE run_id = ? AND item_key <> ''", (run,))
    control.commit()
    forker = FakeForker()
    one = _launcher(control, forker)
    assert len(one.turn()) == 2, "two of the five, the most the run was asked to use"
    assert len(one.turn()) == 0, "the two are still running"
    held = forker.forked[0][3]
    forker.exits = [("plain", jobhost.Exit(held.id, held.pid, 0, None))]
    held.said = [{"event": "unit", "result": {"seconds": 1.0}}]
    assert len(one.turn()) == 1, "one ended, one more starts"
