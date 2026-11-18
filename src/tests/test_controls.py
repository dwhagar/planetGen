"""OPS.234: pause, Stop All, cancel and clear as requests to the supervisor (`planetgen/queue/controls.py`): the request list and its reply, the duty that takes them, the ledger's rows ended, a
job asked to stop through its pipe and killed when it will not, and the hand-outs held while the queues are paused. The last two tests use real job processes forked from a real template."""

import json
import signal
import threading
import time
import types

import pytest

from planetgen.db import store
from planetgen.queue import allowlist, controls, jobhost, launcher, ledger, memorygate, pause, planrows, redisqueue, runrecord, supervisor

R = "tests.template_runners"
HOST = memorygate.Host(16000.0, 12000.0, 8, 0.2)


@pytest.fixture
def conn(redis_server):
    connection = redisqueue.connect(redis_server)
    for key in connection.scan_iter(match=f"{runrecord.PREFIX}:*"):
        connection.delete(key)
    return connection


@pytest.fixture
def control(mysql_config):
    connection = store.get_control_connection(mysql_config, ensure_schema=True)
    connection.execute("DELETE FROM engine_runs")
    connection.commit()
    try:
        yield connection
    finally:
        connection.execute("DELETE FROM engine_runs")
        connection.commit()
        connection.close()


@pytest.fixture
def new_control(mysql_config):
    """What the duty opens when there is no launcher: a connection of its own."""
    return lambda: store.get_control_connection(mysql_config, ensure_schema=True)


class Clock:
    def __init__(self):
        self.t = 5000.0

    def __call__(self):
        return self.t


class FakeJob:
    """A job process as the supervisor holds it: what was sent it, whether it was killed."""

    def __init__(self, job_id, pid):
        self.id, self.pid, self.sent, self.killed = job_id, pid, [], False

    def send(self, message):
        self.sent.append(message)

    def kill(self, number=signal.SIGKILL):
        self.killed = True


class FakeLauncher:
    """The launcher as `Controls` sees it: `running` (job row id to `Started`-like) and the control connection; `ended` takes away what a test says has ended."""

    def __init__(self, control):
        self.control, self.running = control, {}

    def run(self, row_id, engine):
        job = FakeJob(row_id, 900 + row_id)
        self.running[row_id] = types.SimpleNamespace(row={"id": row_id, "run_id": engine}, job=job)
        return job


def _run(control, count=4):
    """A ledger run of one stage: an opening (done) and `count` units, the first two running with attempts, the rest ready."""
    engine = ledger.add_run(control, "generate", "g1")
    ledger.add_jobs(control, engine, [{"stage": "fill", "command": "open", "arguments": {"index": 4, "plan": {}, "run_id": "20261117-x"}}]
                    + [{"stage": "fill", "command": "fill-unit", "item_key": str(n), "arguments": {"args": {"unit": n}}, "claim_mb": 900} for n in range(count)])
    control.execute("UPDATE engine_jobs SET state = 'done' WHERE run_id = ? AND command = 'open'", (engine,))
    control.execute("UPDATE engine_jobs SET state = 'ready' WHERE run_id = ? AND item_key <> ''", (engine,))
    ids = [row["id"] for row in control.execute("SELECT id FROM engine_jobs WHERE run_id = ? AND item_key <> '' ORDER BY id", (engine,)).fetchall()]
    for job_id in ids[:2]:
        ledger.begin_attempt(control, job_id, pid=900 + job_id)
    ledger.set_run_state(control, engine, "running")
    control.commit()
    return engine, ids


def _states(control, engine):
    return [row["state"] for row in control.execute("SELECT state FROM engine_jobs WHERE run_id = ? AND item_key <> '' ORDER BY id", (engine,)).fetchall()]


def _alive(conn, now=1.0):
    supervisor.RedisStore(conn).publish(4242, now, now, 1)


def _duty(conn, running, new_control, clock=None, grace=10.0):
    launches = types.SimpleNamespace(launcher=running)
    return controls.Controls(conn, launches=launches, clock=clock or Clock(), grace=grace, control=new_control)


def _reply(conn, request_id):
    got = conn.lpop(controls.reply_key(request_id))
    return json.loads(runrecord._s(got)) if got is not None else None


# -- the request -------------------------------------------------------------------------------------------------------------

def test_no_request_is_left_when_no_supervisor_is_alive(conn):
    assert controls.ask(conn, "pause", "Dana") is None
    assert conn.llen(supervisor.requests_key()) == 0 and pause.state(conn)["paused"] is False


def test_a_request_is_left_when_one_is_alive_and_a_caller_that_waits_gets_the_reply(conn, new_control):
    _alive(conn)
    duty = _duty(conn, None, new_control)
    answers = []
    asking = threading.Thread(target=lambda: answers.append(controls.ask(conn, "pause", "Dana", wait=10)))
    asking.start()
    deadline = time.monotonic() + 10
    while asking.is_alive() and time.monotonic() < deadline:
        duty.turn()
        time.sleep(0.05)
    asking.join()
    assert answers == [{"ok": True}] and pause.state(conn)["by"] == "Dana"


def test_a_request_nobody_took_is_taken_back_so_the_callers_own_path_runs(conn):
    _alive(conn)
    assert controls.ask(conn, "pause", "Dana", wait=1) is None
    assert conn.llen(supervisor.requests_key()) == 0


def test_a_request_that_asks_for_nothing_the_supervisor_knows_is_refused_before_it_is_left(conn):
    with pytest.raises(ValueError):
        controls.ask(conn, "format-the-disk")


def test_a_bad_request_line_is_dropped_and_the_next_one_is_still_taken(conn, new_control):
    conn.rpush(supervisor.requests_key(), "not json", json.dumps({"id": "a", "op": "resume", "who": "Dana", "arguments": {}}))
    assert _duty(conn, None, new_control).turn() == 2
    assert _reply(conn, "a") == {"ok": False}


def test_pause_and_resume_write_and_clear_the_flag_the_engine_reads(conn, new_control):
    duty = _duty(conn, None, new_control)
    conn.rpush(supervisor.requests_key(), json.dumps({"id": "p", "op": "pause", "who": "Dana", "arguments": {}}))
    duty.turn()
    assert _reply(conn, "p") == {"ok": True} and pause.is_paused(conn)
    conn.rpush(supervisor.requests_key(), json.dumps({"id": "r", "op": "resume", "who": "Dana", "arguments": {}}))
    duty.turn()
    assert _reply(conn, "r") == {"ok": True} and not pause.is_paused(conn)


# -- the ledger's rows -------------------------------------------------------------------------------------------------------

def test_ending_a_runs_rows_cancels_what_is_not_over_closes_the_running_attempts_and_leaves_the_rest(control):
    engine, ids = _run(control)
    other, other_ids = _run(control, count=2)
    ended = controls.end_rows(control, engine)
    control.commit()
    assert ended == {"rows": 4, "attempts": 2, "runs": 1}
    assert _states(control, engine) == ["canceled"] * 4 and _states(control, other) == ["running", "running"], "another run is left as it was"
    assert control.execute("SELECT state FROM engine_runs WHERE id = ?", (engine,)).fetchone()["state"] == "canceled"
    assert control.execute("SELECT state FROM engine_runs WHERE id = ?", (other,)).fetchone()["state"] == "running"
    attempt = ledger.attempt_row(control, ids[0], 1)
    assert (attempt["outcome"], attempt["error_kind"]) == ("lost", "canceled")
    assert control.execute("SELECT state FROM engine_jobs WHERE run_id = ? AND command = 'open'", (engine,)).fetchone()["state"] == "done", "a job that is over stays as it ended"
    assert ledger.attempt_row(control, other_ids[0], 1)["outcome"] == "running"


def test_a_cancel_with_no_supervisor_alive_ends_the_rows_itself_and_a_run_with_none_is_nothing_to_stop(conn, control, new_control, monkeypatch):
    engine, _ids = _run(control)
    conn.hset(runrecord.run_key("20261117-x"), planrows.ENGINE_RUN_FIELD, engine)
    monkeypatch.setattr(planrows, "control_connection", lambda *_: new_control())
    assert controls.stop_run(conn, "20261117-x", "Dana") == {"rows": 4, "attempts": 2, "runs": 1}
    assert _states(control, engine) == ["canceled"] * 4
    assert controls.stop_run(conn, "20261117-none", "Dana") is None


def test_a_generate_page_jobs_run_is_stopped_by_its_ledger_run_with_no_supervisor_alive(conn, control, new_control, monkeypatch):
    """OPS.244: a page job has no Redis run record, so its cancel names the ledger run."""
    engine, _ids = _run(control)
    own = new_control()
    monkeypatch.setattr(store, "get_control_connection", lambda *_args, **_kwargs: own)
    assert controls.stop_engine_run(conn, engine, "Dana") == {"rows": 4, "attempts": 2, "runs": 1}
    assert _states(control, engine) == ["canceled"] * 4


def test_a_generate_page_jobs_run_is_stopped_through_the_supervisor_that_is_alive_and_not_waited_for(conn, control, new_control):
    engine, ids = _run(control)
    other, _other_ids = _run(control, count=2)
    held = FakeLauncher(control)
    mine = held.run(ids[0], engine)
    theirs = held.run(_other_ids[0], other)
    duty = _duty(conn, held, new_control)
    _alive(conn)
    assert controls.stop_engine_run(conn, engine, "Dana")["pending"] is True
    duty.turn()
    assert mine.sent == [{"op": "stop"}] and theirs.sent == []
    assert _states(control, engine) == ["canceled"] * 4 and "canceled" not in _states(control, other)


def test_a_stop_of_a_run_that_cannot_be_reached_never_raises(conn, monkeypatch):
    def refuse(*_args, **_kwargs):
        raise RuntimeError("no database")

    monkeypatch.setattr(store, "get_control_connection", refuse)
    assert controls.stop_engine_run(conn, 12345, "Dana") is None


# -- the jobs ----------------------------------------------------------------------------------------------------------------

def test_a_cancel_asks_each_job_of_the_run_to_stop_cancels_its_rows_and_leaves_other_runs_alone(conn, control, new_control):
    engine, ids = _run(control)
    other, _other_ids = _run(control, count=2)
    held = FakeLauncher(control)
    mine, theirs = held.run(ids[0], engine), held.run(ids[1], engine)
    third = held.run(_other_ids[0], other)
    duty = _duty(conn, held, new_control)
    conn.rpush(supervisor.requests_key(), json.dumps({"id": "c", "op": "stop_jobs", "who": "Dana", "arguments": {"run_id": "20261117-x"}}))
    conn.hset(runrecord.run_key("20261117-x"), planrows.ENGINE_RUN_FIELD, engine)
    duty.turn()
    assert mine.sent == theirs.sent == [{"op": "stop"}] and third.sent == []
    assert _states(control, engine) == ["canceled"] * 4 and "canceled" not in _states(control, other)
    assert _reply(conn, "c") is None, "a cancel does not wait for the jobs: nobody is waiting for the reply"
    held.running.pop(ids[0]), held.running.pop(ids[1])         # the launcher saw both end
    duty.turn()
    assert duty.pending == []


def test_stop_all_waits_for_every_job_and_answers_when_the_last_is_gone(conn, control, new_control):
    engine, ids = _run(control)
    other, other_ids = _run(control, count=2)
    held, clock = FakeLauncher(control), Clock()
    jobs = [held.run(ids[0], engine), held.run(other_ids[0], other)]
    duty = _duty(conn, held, new_control, clock)
    conn.rpush(supervisor.requests_key(), json.dumps({"id": "s", "op": "stop_all", "who": "Dana", "arguments": {}}))
    duty.turn()
    assert all(job.sent == [{"op": "stop"}] for job in jobs)
    assert _reply(conn, "s") is None
    held.running.pop(ids[0])
    clock.t += 1
    duty.turn()
    assert _reply(conn, "s") is None, "one job is still there"
    held.running.pop(other_ids[0])
    duty.turn()
    assert _reply(conn, "s") == {"ok": True, "jobs": 2, "rows": 6}
    assert _states(control, engine) == ["canceled"] * 4 and _states(control, other) == ["canceled"] * 2


def test_a_job_that_will_not_stop_is_killed_after_the_grace_and_is_given_one_more_for_the_launcher_to_see_it_end(conn, control, new_control):
    engine, ids = _run(control)
    held, clock = FakeLauncher(control), Clock()
    job = held.run(ids[0], engine)
    duty = _duty(conn, held, new_control, clock, grace=10.0)
    conn.rpush(supervisor.requests_key(), json.dumps({"id": "s", "op": "stop_all", "who": "Dana", "arguments": {}}))
    duty.turn()
    clock.t += 9
    duty.turn()
    assert not job.killed and _reply(conn, "s") is None
    clock.t += 2
    duty.turn()
    assert job.killed and _reply(conn, "s") is None
    clock.t += 11
    duty.turn()                                   # it never left the launcher's table: the duty gives it up and says how many it killed
    assert _reply(conn, "s")["killed"] == 1 and duty.pending == []


def test_a_supervisor_that_ends_answers_the_caller_that_waits(conn, control, new_control):
    engine, ids = _run(control)
    held = FakeLauncher(control)
    held.run(ids[0], engine)
    duty = _duty(conn, held, new_control)
    conn.rpush(supervisor.requests_key(), json.dumps({"id": "s", "op": "stop_all", "who": "Dana", "arguments": {}}))
    duty.turn()
    duty.clear()
    assert _reply(conn, "s")["jobs"] == 1


# -- the hand-outs while paused ----------------------------------------------------------------------------------------------

class Forker:
    def __init__(self):
        self.forked, self.exits, self.made = [], [], 0

    def fork(self, kind, runner, job=None, template=None):
        self.made += 1
        held = FakeJob(self.made, 1000 + self.made)
        held.messages = lambda: []
        held.close = lambda: None
        self.forked.append(held)
        return "plain", held

    def poll(self, timeout=0.0):
        found, self.exits = self.exits, []
        return found


def test_the_hand_outs_are_held_while_the_queues_are_paused_and_go_on_when_they_are_resumed(conn, control):
    engine, ids = _run(control)
    control.execute("UPDATE engine_jobs SET state = 'ready' WHERE run_id = ? AND item_key <> ''", (engine,))
    control.execute("UPDATE engine_attempts SET outcome = 'lost' WHERE outcome = 'running'")
    control.commit()
    forker = Forker()
    made = launcher.Launcher(control, forker, read_host=lambda: HOST)
    launches = supervisor.Launches(conn)
    launches.launcher = made          # the setting that forks jobs is not what this tests
    pause.pause("Dana", conn)
    assert launches.turn() == 0 and forker.forked == []
    assert set(_states(control, engine)) == {"ready"}
    pause.resume(conn)
    assert launches.turn() > 0 and forker.forked


# -- real processes ----------------------------------------------------------------------------------------------------------

@pytest.fixture
def host(monkeypatch):
    monkeypatch.setenv(allowlist.EXTRA_ENV, f"{R}:*")
    made = jobhost.JobHost(warm=[f"{R}.warm"])
    made.start(ready_seconds=60)
    yield made
    made.stop()


class RealLauncher(FakeLauncher):
    """What the launcher does with the ends the template reports: the job leaves `running`."""

    def __init__(self, control, host):
        super().__init__(control)
        self.host = host

    def ended(self):
        for exit_ in self.host.poll(0.05):
            for row_id, started in list(self.running.items()):
                if started.job.id == exit_.id:
                    del self.running[row_id]


def _until(check, seconds=30):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        if check():
            return True
        time.sleep(0.05)
    return False


def test_stop_all_stops_a_real_job_between_pieces_and_kills_one_that_will_not_stop(conn, control, new_control, host):
    """Two real job processes from a real template: one that obeys the stop request and one that sleeps through it. After Stop All's reply neither is alive."""
    engine, ids = _run(control)
    held = RealLauncher(control, host)
    polite, stubborn = host.fork(f"{R}.chats"), host.fork(f"{R}.lingers")
    for row_id, job in zip(ids[:2], (polite, stubborn)):
        held.running[row_id] = types.SimpleNamespace(row={"id": row_id, "run_id": engine}, job=job)
    duty = _duty(conn, held, new_control, time.time, grace=1.0)
    conn.rpush(supervisor.requests_key(), json.dumps({"id": "s", "op": "stop_all", "who": "Dana", "arguments": {}}))
    replies = []

    def settled():
        held.ended()
        duty.turn()
        replies.append(_reply(conn, "s"))
        return replies[-1] is not None

    assert _until(settled), "Stop All was never answered"
    reply = replies[-1]
    assert reply["ok"] and reply["jobs"] == 2 and reply["killed"] == 1
    assert _until(lambda: not jobhost._pid_alive(polite.pid) and not jobhost._pid_alive(stubborn.pid), 10)
    assert _states(control, engine) == ["canceled"] * 4
    assert control.execute("SELECT state FROM engine_runs WHERE id = ?", (engine,)).fetchone()["state"] == "canceled"


def test_the_workers_tile_counts_the_jobs_the_ledger_has_running_in_every_run(conn, control, mysql_config, monkeypatch):
    """OPS.235 slice 1d: `controls.running_units` is what the Queue page's workers tile shows now that no RQ worker is counted."""
    real = store.get_control_connection
    monkeypatch.setattr(store, "get_control_connection", lambda *a, **k: real(mysql_config, ensure_schema=True))
    base = controls.running_units(conn)
    _run(control)                       # two of its units are running
    _run(control, count=3)              # and two of this one's
    assert controls.running_units(conn) == base + 4
    monkeypatch.setattr(store, "get_control_connection", lambda *a, **k: 1 / 0)
    assert controls.running_units(conn) == 0                # a ledger that cannot be read counts none
