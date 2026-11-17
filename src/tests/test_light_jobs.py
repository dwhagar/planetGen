"""OPS.230: the site's jobs as light jobs (`queue/lightjobs.py`): a record in Redis the supervisor's duty hands out, the deadline of a job nobody waits for, the limits, the dead process, the
shared templates. The duty is driven with a fake clock and a fake forker; the last tests fork real jobs from a real template."""

import os
import signal
import time

import pytest

from planetgen.queue import jobhost, lightjobs, memorygate, redisqueue, runrecord

BIG_HOST = memorygate.Host(10 ** 6, 10 ** 6, 64, 0.0)


class Clock:
    def __init__(self, now=None):
        self.now = time.time() if now is None else now

    def __call__(self):
        return self.now


class FakeJob:
    def __init__(self, job_id, pid):
        self.id, self.pid, self.killed, self.closed = job_id, pid, False, False

    def kill(self, number=signal.SIGKILL):
        self.killed = number

    def close(self):
        self.closed = True

    def messages(self):
        return []


class FakeForker:
    """Jobs that run until the test ends them (`end`)."""

    def __init__(self):
        self.made, self.jobs, self.exits, self.forked = 0, {}, [], []

    def fork(self, kind, runner, job=None, owner=None):
        self.made += 1
        held = FakeJob(self.made, os.getpid())
        self.jobs[job["id"]] = held
        self.forked.append((kind, runner, dict(job), owner))
        return "fake", held

    def end(self, job_id, code=0, number=None):
        held = self.jobs[job_id]
        self.exits.append(("fake", jobhost.Exit(held.id, held.pid, code=None if number else code, signal=number)))

    def poll(self, timeout=0.0, owner=None):
        found, self.exits = self.exits, []
        return found

    def stop(self):
        pass


@pytest.fixture
def duty(redis_server, monkeypatch):
    monkeypatch.setattr(lightjobs, "ask_supervisor", lambda connection=None, url=None: None)      # the test drives the duty itself
    connection = redisqueue.connect()
    clock = Clock()
    forker = FakeForker()
    held = []
    light = lightjobs.LightJobs(connection, forker, read_host=lambda: BIG_HOST, clock=clock, held_mb=lambda pid: 0.0, also_held=lambda: held)
    light.test = (connection, clock, forker, held)
    return light


def _submit(connection, reference="math:pow", args=(2, 3), **options):
    return lightjobs.add_job(connection, reference, list(args), **options)


def test_a_pending_job_is_forked_with_what_it_needs_and_holds_its_claim_until_it_ends(duty):
    connection, clock, forker, _held = duty.test
    job_id = _submit(connection)
    assert duty.turn() == [job_id]
    kind, runner, job, owner = forker.forked[0]
    assert (kind, runner, owner) == (lightjobs.KIND, lightjobs.RUNNER, lightjobs.OWNER)
    assert job["id"] == job_id and job["function"] == "math:pow" and job["args"] == [2, 3] and job["nice"] == lightjobs.PROCESS_NICE[lightjobs.LANE_INTERACTIVE]
    assert lightjobs.status(connection, job_id)["state"] == "running" and duty.busy()
    assert [holder.claim.lane for holder in duty.holders()] == ["interactive"]
    lightjobs.finish(connection, job_id, "succeeded", result=8.0)          # the job process writes its own end
    forker.end(job_id, code=0)
    duty.turn()
    assert lightjobs.status(connection, job_id) == {"id": job_id, "state": "succeeded", "result": 8.0, "error": None, "error_status": None, "waiting": None}
    assert duty.holders() == [] and not duty.busy()


def test_a_job_whose_process_ended_without_writing_its_end_fails_naming_how_it_died(duty):
    connection, _clock, forker, _held = duty.test
    killed, quiet = _submit(connection), _submit(connection)
    duty.turn()
    forker.end(killed, number=signal.SIGKILL)
    forker.end(quiet, code=1)
    duty.turn()
    assert lightjobs.status(connection, killed)["error"] == "the job process ended on SIGKILL"
    assert lightjobs.status(connection, quiet)["error"] == "the job process exited with code 1 and said nothing"
    assert not connection.smembers(lightjobs.running_key())


def test_a_job_that_runs_past_its_time_limit_is_killed_and_fails_with_the_reason(duty):
    connection, clock, forker, _held = duty.test
    job_id = _submit(connection, expected_seconds=None)
    connection.hset(lightjobs.job_key(job_id), "limit", "100")
    duty.turn()
    clock.now += 99
    assert duty.limits() == []
    clock.now += 2
    assert duty.limits() == [job_id] and forker.jobs[job_id].killed == signal.SIGKILL
    forker.end(job_id, number=signal.SIGKILL)
    duty.turn()
    body = lightjobs.status(connection, job_id)
    assert body["state"] == "failed" and body["error"] == "the job ran past its time limit of 100 s and was ended"


def test_a_job_asked_to_stop_is_killed_and_a_pending_one_is_taken_off_the_line(duty):
    connection, _clock, forker, _held = duty.test
    running, waiting = _submit(connection), None
    duty.turn()
    waiting = _submit(connection)
    assert lightjobs.request_stop(connection, waiting) is True and lightjobs.status(connection, waiting)["error"] == lightjobs.CANCELED_TEXT
    assert lightjobs.request_stop(connection, running) is True
    assert duty.limits() == [running] and forker.jobs[running].killed
    assert lightjobs.request_stop(connection, "0" * 16) is False


def test_a_job_with_a_deadline_is_dropped_once_nobody_asks_after_it_and_its_keys_go(duty):
    connection, clock, _forker, _held = duty.test
    connection.set("the-claim", "x")
    job_id = _submit(connection, wanted_seconds=30, drop_keys=["the-claim"])
    clock.now += 20
    assert duty.drop() == []
    lightjobs.status(connection, job_id)                    # a visitor asks: the deadline moves
    wanted_until = float(connection.hget(lightjobs.job_key(job_id), "wanted_until"))
    assert wanted_until > time.time() + 25
    clock.now = wanted_until + 1
    assert duty.drop() == [job_id]
    body = lightjobs.status(connection, job_id, ask=False)
    assert body["state"] == "failed" and body["error"] == lightjobs.DROPPED_TEXT
    assert connection.get("the-claim") is None and not connection.zscore(lightjobs.pending_key(), job_id)


def test_a_job_that_began_or_has_no_deadline_is_never_dropped(duty):
    connection, clock, _forker, _held = duty.test
    started, never = _submit(connection, wanted_seconds=1), _submit(connection)
    duty.turn()                                             # both start at once: nothing to drop
    later = _submit(connection)
    clock.now += 10 ** 6
    assert duty.drop() == [] and lightjobs.status(connection, later)["state"] == "queued"
    assert lightjobs.status(connection, started, ask=False)["state"] == "running" and lightjobs.status(connection, never, ask=False)["state"] == "running"


def test_a_light_jobs_claim_is_its_api_functions_figure_on_its_lane(duty):
    from planetgen.queue import workers
    connection, _clock, _forker, _held = duty.test
    long = _submit(connection, "planetgen.queue.api_jobs:regenerate_sector", ({}, {}), expected_seconds=3600)
    claim = duty.claim_of(lightjobs._record(connection, long))
    assert claim == memorygate.Claim("api:planetgen.queue.api_jobs:regenerate_sector", float(workers.SECTOR_WORKER_MEMORY_MB), 1.0, "long")
    other = duty.claim_of({"function": "planetgen.nothing:here", "lane": "interactive"})
    assert other.kind == "fn:planetgen.nothing:here" and other.mb == 1000.0 and other.lane == "interactive"


def test_a_prefetch_chunk_is_a_maintenance_job_the_scheduler_weighs_as_the_sites_and_is_dropped_with_its_claims(duty):
    """OPS.230 slice B: the node prefetch's chunks are light jobs on the maintenance lane, claimed as the API function (`api:` kind, `scheduler.SITE_KINDS`), wanted for `CLAIM_SECONDS`, and a
    chunk nobody started in that time goes with the claims of its nodes."""
    from planetgen.queue import scheduler
    from planetgen.web.lib import tileprefetch
    connection, clock, _forker, _held = duty.test
    keys = [f"n{number}" for number in range(tileprefetch.CHUNK + 4)]
    for key in keys:
        connection.set(tileprefetch._claim_name("db", "gen", key), "1")
    assert tileprefetch._add_chunks(connection, "db", "gen", keys) == len(keys)
    assert lightjobs.pending_count(connection, tileprefetch.JOB) == 2 and lightjobs.pending_count(connection, "math:pow") == 0
    tickets, records = duty.tickets()
    assert [ticket.claim.lane for ticket in tickets] == ["maintenance", "maintenance"]
    assert {ticket.claim.kind for ticket in tickets} == {f"api:{tileprefetch.JOB}"} and set(scheduler.SITE_KINDS) == {f"api:{tileprefetch.JOB}"}
    assert scheduler.weight(tickets[0].claim) == scheduler.SITE_KINDS[f"api:{tileprefetch.JOB}"]
    job_id, record = next(iter(records.values()))
    assert float(record["wanted"]) == tileprefetch.CLAIM_SECONDS and int(record["limit"]) >= tileprefetch.MIN_JOB_SECONDS
    clock.now += tileprefetch.CLAIM_SECONDS + 1
    assert len(duty.drop()) == 2
    assert all(connection.get(tileprefetch._claim_name("db", "gen", key)) is None for key in keys)
    assert lightjobs.pending_count(connection, tileprefetch.JOB) == 0


def test_a_held_host_still_runs_the_sites_jobs_and_holds_the_upkeep_a_visit_asked_for(duty):
    """OPS.231: while the memory controller holds admission the interactive and long lanes go; the node prefetch (maintenance lane) waits, with the reason on its record."""
    connection, _clock, forker, _held = duty.test
    duty.hold = lambda: "admission is held: the host is short of memory"
    visitor = _submit(connection)
    prefetch = _submit(connection, "planetgen.web.lib.tileprefetch:build_nodes", ("db", []), lane="maintenance")
    assert duty.turn() == [visitor]
    record = lightjobs._record(connection, prefetch)
    assert record["state"] == "queued" and record["waited"].startswith("admission is held")
    duty.hold = lambda: None
    assert duty.turn() == [prefetch]


def test_the_lane_follows_the_expected_time_and_the_setting(redis_server, monkeypatch):
    monkeypatch.setattr(lightjobs, "ask_supervisor", lambda connection=None, url=None: None)
    from planetgen.util.settings import get_settings
    connection = redisqueue.connect()
    monkeypatch.setattr(get_settings().queue, "long_job_seconds", 5.0)
    lanes = [lightjobs.lane_of(expected) for expected in (None, float("nan"), 5, 6, 1e9)]
    assert lanes == ["interactive", "interactive", "interactive", "long", "long"]
    assert connection.hget(lightjobs.job_key(_submit(connection, expected_seconds=600)), "lane") == b"long"


def test_a_visitors_job_is_handed_out_at_once_while_the_bulk_lane_holds_the_host(duty):
    """The site's lanes are judged before the others and against the places kept for them (`dispatch.dispatch`)."""
    connection, _clock, _forker, held = duty.test
    duty.read_host = lambda: memorygate.Host(16_000, 3_000, 4, 0.5)          # the bulk lane has taken what it may
    held.append(memorygate.Holder("bulk-1", memorygate.Claim("unit:generate-system", 14_000, 1.0, memorygate.LANE_BULK), 1, "launcher", 14_000, 14_000))
    job_id = _submit(connection)
    assert duty.turn() == [job_id]


def test_a_job_that_cannot_start_yet_says_why_and_keeps_its_place(duty):
    connection, _clock, _forker, _held = duty.test
    job_id = _submit(connection)
    first = connection.zscore(lightjobs.pending_key(), job_id)

    def not_ready(*args, **options):
        raise jobhost.NotReady("the template 'plain' is being built")

    duty.forker.fork = not_ready
    assert duty.turn() == [] and lightjobs.waited(connection, job_id) == "the template 'plain' is being built"
    assert connection.zscore(lightjobs.pending_key(), job_id) == first and lightjobs.status(connection, job_id)["state"] == "queued"


def test_a_job_running_when_its_supervisor_ended_is_settled_by_the_next_one(duty):
    connection, clock, _forker, _held = duty.test
    gone, alive, nopid = _submit(connection), _submit(connection), _submit(connection)
    for job_id, pid in ((gone, 2 ** 22 + 99), (alive, os.getpid()), (nopid, None)):
        connection.sadd(lightjobs.running_key(), job_id)
        connection.zrem(lightjobs.pending_key(), job_id)
        connection.hset(lightjobs.job_key(job_id), mapping={"state": "running", "started": f"{clock.now:.3f}", **({} if pid is None else {"pid": pid})})
    assert duty.recover() == [gone]                         # the dead one is found at once, the one with no process is given a moment
    assert duty.recover() == []                             # and it is looked for every RECOVER_SECONDS, not at every turn
    clock.now += lightjobs.RECOVER_SECONDS + lightjobs.NO_PID_SECONDS
    assert duty.recover() == [nopid]
    assert lightjobs.status(connection, gone)["state"] == "failed" and lightjobs.status(connection, alive, ask=False)["state"] == "running"


def test_the_templates_tell_each_owner_only_of_its_own_jobs():
    class Host:
        def __init__(self, exits):
            self.exits = exits

        def poll(self, timeout=0.0):
            found, self.exits = self.exits, []
            return found

    templates = jobhost.TemplateSet()
    templates.hosts = {"plain": Host([jobhost.Exit(1, 11, code=0), jobhost.Exit(2, 12, code=0)])}
    templates._owners = {("plain", 1): None, ("plain", 2): lightjobs.OWNER}
    assert [(key, exit_.id) for key, exit_ in templates.poll(owner=lightjobs.OWNER)] == [("plain", 2)]
    assert [(key, exit_.id) for key, exit_ in templates.poll()] == [("plain", 1)]
    assert templates.poll() == [] and templates.poll(owner=lightjobs.OWNER) == []


def test_the_supervisor_has_something_to_watch_while_a_light_job_waits_or_runs(redis_server, monkeypatch):
    monkeypatch.setattr(lightjobs, "ask_supervisor", lambda connection=None, url=None: None)
    from planetgen.queue import supervisor
    connection = redisqueue.connect()
    store = supervisor.RedisStore(connection)
    assert not store.busy()
    job_id = _submit(connection)
    assert store.busy()
    assert lightjobs.cancel(connection, job_id) and not store.busy()


def test_an_update_and_stop_all_see_the_light_jobs(redis_server, monkeypatch):
    monkeypatch.setattr(lightjobs, "ask_supervisor", lambda connection=None, url=None: None)
    from planetgen.admin import updateguard
    from planetgen.queue import stopall
    connection = redisqueue.connect()
    job_id = _submit(connection)
    source = updateguard.QueuedJobSource.__new__(updateguard.QueuedJobSource)
    source._connection = lambda: connection
    listed = [item for item in source._light_jobs(connection)]
    assert [item["id"] for item in listed] == [job_id] and listed[0]["queued"] is True and listed[0]["holds"] == "the interactive lane"
    assert stopall._empty_queues(connection) >= 1
    assert lightjobs.status(connection, job_id, ask=False)["error"] == lightjobs.CANCELED_TEXT


def test_the_job_names_what_it_waits_for_to_its_status(redis_server, monkeypatch):
    monkeypatch.setattr(lightjobs, "ask_supervisor", lambda connection=None, url=None: None)
    connection = redisqueue.connect()
    job_id = _submit(connection)
    lightjobs.finish(connection, job_id, "failed", error="x")
    connection.hset(lightjobs.job_key(job_id), mapping={"state": "running", "waiting": '{"run": "r1", "position": 2, "text": "waiting"}'})
    assert lightjobs.status(connection, job_id, ask=False)["waiting"] == {"run": "r1", "position": 2, "text": "waiting"}
    lightjobs._current.update(id=job_id, url=None)
    try:
        lightjobs.report_waiting(None)
    finally:
        lightjobs._current.update(id=None, url=None)
    assert lightjobs.status(connection, job_id, ask=False)["waiting"] is None


# -- real processes ------------------------------------------------------------------------------------------------------------

@pytest.fixture
def forking(redis_server, monkeypatch):
    """A real template that forks the light jobs of this test, and the duty that hands them out."""
    from planetgen.queue import allowlist
    monkeypatch.setenv(allowlist.EXTRA_ENV, "tests.template_runners:*,tests.test_api_jobs:*,math:*,time:*")
    templates = jobhost.TemplateSet(ready_seconds=60)
    connection = redisqueue.connect()
    duty = lightjobs.LightJobs(connection, templates, read_host=lambda: BIG_HOST, runner="tests.template_runners.light_job_in_the_tests_queue", extra={"prefix": runrecord.PREFIX})
    yield duty
    duty.clear()


def _turns(duty, connection, job_id, states, seconds=60):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        duty.turn()
        if lightjobs.status(connection, job_id, ask=False)["state"] in states:
            return lightjobs.status(connection, job_id, ask=False)
        time.sleep(0.1)
    raise AssertionError(f"job {job_id} did not reach {states}: {lightjobs.status(connection, job_id, ask=False)}")


def test_a_light_job_runs_in_a_process_forked_from_the_template_and_writes_its_own_result(forking):
    connection = forking.connection
    job_id = _submit(connection, "math:pow", (2, 10))
    assert _turns(forking, connection, job_id, ("succeeded", "failed"))["result"] == 1024.0
    refused = _submit(connection, "tests.test_api_jobs:_refuse", (409,))
    body = _turns(forking, connection, refused, ("succeeded", "failed"))
    assert body["state"] == "failed" and body["error"] == "not this" and body["error_status"] == 409
    broke = _submit(connection, "tests.test_api_jobs:_raise_value_error", ())
    assert "boom" in _turns(forking, connection, broke, ("failed",))["error"]


def test_a_light_job_killed_from_outside_fails_naming_the_signal_and_one_over_its_limit_is_ended(forking):
    connection = forking.connection
    sleeping = _submit(connection, "time:sleep", (60,))
    _turns(forking, connection, sleeping, ("running",))
    os.kill(int(connection.hget(lightjobs.job_key(sleeping), "pid")), signal.SIGKILL)
    assert _turns(forking, connection, sleeping, ("failed",))["error"] == "the job process ended on SIGKILL"
    slow = _submit(connection, "time:sleep", (60,))
    connection.hset(lightjobs.job_key(slow), "limit", "1")
    body = _turns(forking, connection, slow, ("failed",))
    assert body["error"] == "the job ran past its time limit of 1 s and was ended"
