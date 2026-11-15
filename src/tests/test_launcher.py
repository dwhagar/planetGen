"""OPS.228, slice A: the launcher (`planetgen/queue/launcher.py`) hands out the ready rows of the ledger by the dispatch rule, writes each attempt before it forks, and closes each attempt with how its process ended.
The templates are stood in for by a fake forker; the rows are real."""

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
    ledger.add_jobs(control, run, [{"stage": "fill", "command": "fill-unit", "item_key": str(n), "arguments": {"unit": n}, "claim_mb": mb} for n in range(count)])
    control.execute("UPDATE engine_jobs SET state = 'ready' WHERE run_id = ?", (run,))
    control.commit()
    return run


def _state(control, run):
    return [(row["item_key"], row["state"], row["attempt"]) for row in control.execute("SELECT item_key, state, attempt FROM engine_jobs WHERE run_id = ? ORDER BY id", (run,)).fetchall()]


def _launcher(control, forker, host=HOST):
    return launcher.Launcher(control, forker, read_host=lambda: host)


def test_every_ready_job_that_fits_is_handed_out_with_its_attempt_written_first(control):
    run = _ready(control)
    forker = FakeForker()
    started = _launcher(control, forker).turn()
    assert len(started) == 3
    assert _state(control, run) == [("0", "running", 1), ("1", "running", 1), ("2", "running", 1)]
    kind, runner, job, _held = forker.forked[0]
    assert (kind, runner) == ("fill-unit", launcher.RUNNER) and job["unit"] == 0 and job["attempt"] == 1 and job["row"] == started[0]
    attempt = ledger.attempt_row(control, started[0], 1)
    assert attempt["outcome"] == "running" and attempt["pid"] == 1001


def test_a_job_that_does_not_fit_waits_with_its_reason_on_the_row(control):
    run = _ready(control, count=3, mb=7000)
    forker = FakeForker()
    one = _launcher(control, forker)
    started = one.turn()
    assert 0 < len(started) < 3
    waiting = control.execute("SELECT waited FROM engine_jobs WHERE run_id = ? AND state = 'ready'", (run,)).fetchall()
    assert waiting and all(row["waited"] for row in waiting)


def test_the_end_of_a_process_closes_its_attempt_and_moves_its_job(control):
    run = _ready(control, count=3)
    forker = FakeForker()
    one = _launcher(control, forker)
    ids = one.turn()
    done, failed, killed = ids
    forker.forked[1][3].said.append({"event": "failed", "error": "ValueError: no star"})
    forker.exits = [("plain", jobhost.Exit(1, 1001, code=0)), ("plain", jobhost.Exit(2, 1002, code=1)), ("plain", jobhost.Exit(3, 1003, signal=signal.SIGKILL))]
    one.turn()
    assert [row[1] for row in _state(control, run)] == ["done", "failed", "failed"]
    assert ledger.attempt_row(control, done, 1)["outcome"] == "done"
    assert (ledger.attempt_row(control, failed, 1)["outcome"], ledger.attempt_row(control, failed, 1)["error_text"]) == ("failed", "ValueError: no star")
    row = ledger.attempt_row(control, killed, 1)
    assert row["outcome"] == "lost" and "SIGKILL" in row["error_text"]
    assert one.running == {} and all(item[3].closed for item in forker.forked)


def test_a_job_that_is_running_is_not_handed_out_again(control):
    run = _ready(control, count=2)
    forker = FakeForker()
    one = _launcher(control, forker)
    one.turn()
    assert one.turn() == [] and len(forker.forked) == 2 and [row[2] for row in _state(control, run)] == [1, 1]
