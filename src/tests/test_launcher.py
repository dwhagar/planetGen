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
        self.forked, self.exits, self.made, self.templates = [], [], 0, []

    def fork(self, kind, runner, job=None, template=None):
        self.made += 1
        held = FakeJob(self.made, 1000 + self.made)
        self.forked.append((kind, runner, job, held))
        self.templates.append(template)
        return template or "plain", held

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

        def fork(self, kind, runner, job=None, template=None):
            if self.building:
                raise jobhost.NotReady("the template 'plain' is being built")
            return super().fork(kind, runner, job, template)

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


def test_a_launcher_with_redis_claims_what_the_gate_has_learned_not_the_minimum(control, redis_server, monkeypatch):
    from planetgen.queue import claims, redisqueue
    conn = redisqueue.connect(redis_server)
    monkeypatch.setattr(claims, "unit_table_mb", lambda kind: 3000.0)
    monkeypatch.setattr(claims, "unit_tokens", lambda kind: ("mysql-write",))
    run = _ready(control, count=1, mb=None)
    one = launcher.Launcher(control, FakeForker(), read_host=lambda: HOST, redis=conn)
    row = control.execute("SELECT id, run_id, stage, command, item_key, arguments, claim_mb, job_class, attempt FROM engine_jobs WHERE run_id = ? AND item_key <> ''", (run,)).fetchone()
    claim = one.claim(row)
    assert (claim.kind, claim.mb, claim.cpu, claim.tokens) == ("unit:fill-unit", 3000.0, claims.UNIT_CPU, ("mysql-write",))


def _redis_launcher(control, redis_server, monkeypatch, forker=None):
    from planetgen.queue import claims, redisqueue
    conn = redisqueue.connect(redis_server)
    monkeypatch.setattr(claims, "unit_table_mb", lambda kind: 2000.0)
    monkeypatch.setattr(claims, "unit_tokens", lambda kind: ())
    return conn, launcher.Launcher(control, forker or FakeForker(), read_host=lambda: HOST, redis=conn)


def _row(control, run):
    return control.execute("SELECT id, run_id, stage, command, item_key, arguments, claim_mb, job_class, attempt FROM engine_jobs WHERE run_id = ? AND item_key <> ''", (run,)).fetchone()


def test_a_unit_whose_process_was_killed_is_handed_out_again_alone_on_a_claim_grown_by_one_and_a_half(control, redis_server, monkeypatch):
    """OPS.236 (bug since OPS.235 deleted `units.horse_killed`): the retry of a job the kernel killed had the same claim and shared the host."""
    forker = FakeForker()
    _conn, one = _redis_launcher(control, redis_server, monkeypatch, forker)
    run = _ready(control, count=1, mb=None)
    (job,) = one.turn()
    first = one.claim(_row(control, run))
    assert first.mb == 2000.0 and not first.solo, "the first hand-out is an ordinary claim"
    ledger.begin_attempt(control, job)
    forker.exits = [("plain", jobhost.Exit(1, 1001, signal=signal.SIGKILL))]
    again = one.turn()
    assert again == [job], "the unit goes out once more"
    retry = one.claim(_row(control, run))
    assert retry.mb == pytest.approx(3000.0, rel=1e-6) and retry.solo, (retry.mb, retry.tokens)
    assert one.running[job].claim.mb == pytest.approx(3000.0, rel=1e-6) and one.running[job].claim.solo, "and it was admitted as the grown, solo claim"


def test_a_unit_that_failed_without_a_signal_or_was_lost_with_its_supervisor_is_not_grown_or_made_to_run_alone(control, redis_server, monkeypatch):
    forker = FakeForker()
    _conn, one = _redis_launcher(control, redis_server, monkeypatch, forker)
    run = _ready(control, count=2, mb=None)
    first, second = one.turn()
    ledger.begin_attempt(control, first)
    forker.exits = [("plain", jobhost.Exit(1, 1001, code=3))]
    one.turn()
    ledger.begin_attempt(control, second)
    forker.exits = [("plain", jobhost.Exit(2, 1002))]       # ended, and nothing says how: lost, not killed
    one.turn()
    for row in control.execute("SELECT id, run_id, stage, command, item_key, arguments, claim_mb, job_class, attempt FROM engine_jobs WHERE run_id = ? AND item_key <> ''", (run,)).fetchall():
        claim = one.claim(row)
        assert claim.mb == 2000.0 and not claim.solo, (row["item_key"], claim)


def test_a_killed_unit_that_is_killed_again_has_no_third_hand_out(control, redis_server, monkeypatch):
    forker = FakeForker()
    _conn, one = _redis_launcher(control, redis_server, monkeypatch, forker)
    run = _ready(control, count=1, mb=None)
    (job,) = one.turn()
    for number in (1, 2):
        ledger.begin_attempt(control, job)
        forker.exits = [("plain", jobhost.Exit(number, 1000 + number, signal=signal.SIGKILL))]
        started = one.turn()
    assert started == [] and _state(control, run) == [("0", "failed", 2)]


def _upkeep(control, count=4, name="warm-unit"):
    """OPS.229: a cache warm-up's rows (class `maintenance`, as `planrows.job_class` gives a run of its kind)."""
    run = ledger.add_run(control, "maintenance", "g1")
    ledger.add_jobs(control, run, [{"stage": "warm", "command": "open", "arguments": {"index": 0, "plan": {}, "run_id": "20261115-w"}, "job_class": "maintenance"}]
                    + [{"stage": "warm", "command": name, "item_key": str(n), "arguments": {"args": {"unit": n}, "cost": 1.0}, "claim_mb": 300, "job_class": "maintenance"} for n in range(count)])
    control.execute("UPDATE engine_jobs SET state = 'done' WHERE run_id = ? AND item_key = ''", (run,))
    control.execute("UPDATE engine_jobs SET state = 'ready' WHERE run_id = ? AND item_key <> ''", (run,))
    control.commit()
    return run


def test_the_units_of_a_warm_up_or_a_check_are_the_maintenance_lanes_and_fork_from_the_template_that_holds_the_app(control):
    upkeep, bulk = _upkeep(control, count=2), _ready(control, count=1)
    forker = FakeForker()
    one = _launcher(control, forker)
    assert len(one.turn()) == 3
    by_run = {job["row"]: template for (_kind, _runner, job, _held), template in zip(forker.forked, forker.templates)}
    classes = {row["id"]: (row["run_id"], row["job_class"]) for row in control.execute("SELECT id, run_id, job_class FROM engine_jobs WHERE item_key <> ''").fetchall()}
    assert {(classes[row][0] == upkeep, template) for row, template in by_run.items()} == {(True, "app"), (False, None)}
    assert {holder.claim.lane for holder in one.holders() if "warm" in holder.claim.kind} == {memorygate.LANE_MAINTENANCE}
    assert bulk and all(classes[row][1] == "bulk" for row, template in by_run.items() if template is None)


def test_the_maintenance_lane_runs_no_more_than_its_share_while_a_generation_run_is_live_and_what_the_gate_allows_otherwise(control, redis_server, monkeypatch):
    from planetgen.queue import redisqueue, runrecord
    conn = redisqueue.connect(redis_server)
    live = []
    monkeypatch.setattr(launcher.runrecord, "active_run", lambda connection, kind: "g-run" if live and kind == runrecord.GENERATION else None)
    import types
    monkeypatch.setattr(launcher, "get_settings", lambda: types.SimpleNamespace(queue=types.SimpleNamespace(maintenance_workers=2)))
    _upkeep(control, count=5)
    forker = FakeForker()
    one = launcher.Launcher(control, forker, read_host=lambda: HOST, redis=conn)
    live.append(True)
    assert len(one.turn()) == 2, "a generation run is live: the upkeep takes its two and no more"
    assert len(one.turn()) == 0
    live.clear()
    assert len(one.turn()) == 3, "no generation run: the three that are left go (the gate is the limit)"


def _page(control, steps=2):
    """A Generate page job's run (`pagerun.start`) with `steps` steps."""
    from planetgen.queue import pagerun
    job = {"id": "20261118-aaaa", "kind": "reset", "title": "A job", "admin": "boss", "database": None, "draft": None, "created_at": 1.0, "cwd": "/tmp", "root": "/tmp/jobs"}
    run = pagerun.start(control, job, [{"label": f"Step {n}", "argv": ["true"], "stages": None} for n in range(1, steps + 1)])
    return run, [row["id"] for row in control.execute("SELECT id FROM engine_jobs WHERE run_id = ? AND command = ? ORDER BY item_key", (run, pagerun.STEP)).fetchall()]


def test_a_step_of_a_generate_page_job_is_handed_out_as_a_page_step_with_an_attempt_the_launcher_writes(control):
    run, (first, second) = _page(control)
    forker = FakeForker()
    started = _launcher(control, forker).turn()
    assert started == [first], "only the step that is ready: the next one waits for it"
    _kind, _runner, job, held = forker.forked[0]
    assert job == {"page_step": True, "run": run, "row": first, "kind": "page-step", "retries_left": 0}
    attempt = control.execute("SELECT pid, outcome FROM engine_attempts WHERE job_id = ?", (first,)).fetchone()
    assert (attempt["pid"], attempt["outcome"]) == (held.pid, "running")


def test_a_page_step_whose_process_was_killed_is_failed_with_the_steps_after_it_canceled_and_never_handed_out_again(control):
    run, (first, second) = _page(control)
    forker = FakeForker()
    one = _launcher(control, forker)
    assert one.turn() == [first]
    forker.exits = [("plain", jobhost.Exit(1, 1001, signal=signal.SIGKILL))]
    assert one.turn() == []
    states = [row["state"] for row in control.execute("SELECT state FROM engine_jobs WHERE run_id = ? AND command = 'page-step' ORDER BY item_key", (run,)).fetchall()]
    assert states == ["failed", "canceled"] and len(forker.forked) == 1
    assert control.execute("SELECT state FROM engine_runs WHERE id = ?", (run,)).fetchone()["state"] == "held"
    assert ledger.attempt_row(control, first, 1)["outcome"] == "lost"


def test_a_page_step_left_running_by_a_supervisor_that_ended_is_failed_not_handed_out_again(control):
    run, (first, second) = _page(control)
    ledger.set_job_state(control, first, "running")
    ledger.begin_attempt(control, first, pid=2_000_000_000)
    control.commit()
    assert _launcher(control, FakeForker()).recover() == [first]
    states = [row["state"] for row in control.execute("SELECT state FROM engine_jobs WHERE run_id = ? AND command = 'page-step' ORDER BY item_key", (run,)).fetchall()]
    assert states == ["failed", "canceled"]
