"""
OPS.143: every rule of the unit heartbeat, driven end to end by one fake clock.

`test_beat_judge.py` holds the judge's rules one at a time and `test_unit_watch.py`, `test_unit_stuck.py`, `test_job_beat.py` the follower, the stuck action and the job beats
one at a time. This file runs the whole chain (a horse that writes beats, `unitwatch.UnitWatch`, `unitstuck.Stuck`, `runrecord.is_stale`, the Queue page's job rows) over hours of
fake time with the Redis server clock replaced, the way a real galaxy would run, so the rules are checked in the order and the combinations a person would meet them:

* a unit that works, slows, waits, goes quiet and is stopped, with what the run log and `error.log` say at each change;
* the false kill: a counter that moves at ten times its estimate for six hours is flagged and never touched;
* a unit that dies, and a unit that is stuck twice (queued again once, then the run fails);
* the stale rule of OPS.141 for every mix of states, and the job kinds of OPS.142.
"""

import json

import pytest
import rq

from planetgen.queue import beatjudge, commands, jobbeat, pause, redisqueue, runrecord, unitbeat, units, unitwatch

POLL_SECONDS = 30.0
"""float: How often the follower polls in these runs (the real poll is about this)."""


class Clock:
    """The Redis server clock of the run, set by the test."""

    def __init__(self, now):
        self.now = now


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture(autouse=True)
def _no_hold_left(conn):
    pause.resume(conn)
    yield
    pause.resume(conn)


@pytest.fixture
def clock(conn, monkeypatch):
    fake = Clock(unitbeat.server_time(conn))
    monkeypatch.setattr(unitbeat, "server_time", lambda _conn: fake.now)
    return fake


@pytest.fixture
def run_id(conn, monkeypatch):
    monkeypatch.setitem(commands.REGISTRY, "heartbeatkind", commands.CommandSpec("heartbeatkind", "A test", unit=(__name__, "_unit"), expand=(__name__, "_unit")))
    run = runrecord.start_run(conn, "T", [runrecord.Stage("heartbeatkind", "Heartbeat", units=3)])
    runrecord.begin_stage(conn, run, 0)
    units.enqueue_units(conn, run, 0, "heartbeatkind", [({"n": n}, 1.0) for n in range(3)])
    yield run
    redisqueue.queue(redisqueue.BULK_QUEUE, conn).empty()


def _unit(args, report):
    return {}


class Horse:
    """What one unit's horse writes: a beat every step of the fake clock, and a counter that moves at `rate` items a second."""

    def __init__(self, conn, run_id, clock, unit=0, estimate=1000.0, phase="drawing"):
        self.conn, self.run_id, self.clock, self.unit = conn, run_id, clock, unit
        self.started, self.done, self.estimate, self.moved = clock.now, 0.0, estimate, None
        self.phase, self.waiting_on, self.restarts = phase, "", 0
        self.waited, self._wait_began = 0.0, None

    def write(self):
        record = {"done": self.done, "estimate": self.estimate, "beat": self.clock.now, "started": self.started, "phase": self.phase, "waiting_on": self.waiting_on,
                  "restarts": self.restarts, "waited": self.waited, "pid": 999_999_999, "host": "elsewhere"}
        if self.moved is not None:
            record["moved"] = self.moved
        self.conn.hset(runrecord.progress_key(self.run_id, 0), self.unit, json.dumps(record))

    def wait(self, tag):
        self.waiting_on, self._wait_began = tag, self.clock.now

    def end_wait(self):
        """The wait is over: the time it took is not time the counter stood still (the horse adds it to `waited`, which a counter change resets)."""
        self.waited += self.clock.now - self._wait_began
        self.waiting_on, self._wait_began = "", None

    def tick(self, seconds, rate, beating):
        self.clock.now += seconds
        if rate > 0:
            self.done += rate * seconds
            self.moved, self.waited = self.clock.now, 0.0
        if beating:
            self.write()


def run_for(seconds, horse, watch, rate=0.0, beating=True, step=10.0, others=()):
    """Advances the fake clock `seconds` in `step`s; every step the horse(s) write, every `POLL_SECONDS` the follower polls. Returns the verdicts of the last poll."""
    verdicts, last_poll = {}, horse.clock.now
    for _ in range(int(seconds / step)):
        horse.tick(step, rate, beating)
        for other in others:
            other.tick(step, 0.0, True)
        if horse.clock.now - last_poll >= POLL_SECONDS:
            verdicts, last_poll = watch.check(), horse.clock.now
    return verdicts


class Sink:
    """What `error.log` and the stop command were told."""

    def __init__(self):
        self.errors, self.stopped = [], []


def make_watch(conn, run_id, history=None, sink=None, solo=False):
    from planetgen.queue import unitstuck
    sink = sink or Sink()
    stuck = unitstuck.Stuck(conn, run_id, 0, dump=lambda raw: "STACK", stop=sink.stopped.append, record_error=lambda message, **fields: sink.errors.append(message),
                            solo=solo)
    watch = unitwatch.UnitWatch(conn, run_id, 0, "heartbeatkind", history_of=lambda kind: history or beatjudge.History(), record_error=lambda message, **fields: sink.errors.append(message),
                                stuck=stuck)
    stuck._forget = watch.forget
    return watch, sink


def lines(conn, run_id, kind="unit-state"):
    return [fields["message"] for _, fields in runrecord.read_events(conn, run_id) if fields.get("kind") == kind]


def state_of(watch, unit=0):
    return watch._current[unit]


# --- one unit through the rules, in order ----------------------------------------------------------------------

def test_a_unit_works_slows_waits_goes_quiet_and_is_stopped_each_change_logged_once(conn, clock, run_id):
    watch, sink = make_watch(conn, run_id, history=beatjudge.History(longest_gap=20.0, expected_seconds=1000.0))   # 1 item a second is the normal rate of 1,000
    horse = Horse(conn, run_id, clock)

    run_for(120, horse, watch, rate=1.0)
    assert state_of(watch).state == beatjudge.WORKING and lines(conn, run_id) == [], "a unit that works is not announced"

    run_for(480, horse, watch, rate=0.05)             # the count crawls: its average rate falls under a third of the expected one
    assert state_of(watch).state == beatjudge.SLOW and state_of(watch).severity == beatjudge.FLAG and not state_of(watch).actionable
    assert len(lines(conn, run_id)) == 1 and "is slow" in lines(conn, run_id)[0] and sink.errors == [], "slow is a flag in the run log, not an error"
    run_for(60, horse, watch, rate=0.05)
    assert len(lines(conn, run_id)) == 1, "a state is logged once, not at every poll"

    horse.wait("sky lock")                            # a wait on a lock is not stuck, whatever the counter does
    run_for(1500, horse, watch)
    assert state_of(watch).state == beatjudge.WAITING and not state_of(watch).over_limit and not state_of(watch).actionable
    run_for(600, horse, watch)
    assert state_of(watch).state == beatjudge.WAITING and state_of(watch).over_limit and state_of(watch).severity == beatjudge.FLAG
    assert watch.agreed() == [] and sink.stopped == [] and sink.errors == [], "a long wait is a flag, the holder is what gets judged"

    horse.end_wait()                                  # the wait ends after 35 minutes and the count does not move again
    run_for(240, horse, watch)
    assert state_of(watch).state in (beatjudge.WORKING, beatjudge.SLOW) and sink.stopped == [], "the 35 minutes of waiting are not 35 minutes of standing still"
    run_for(120, horse, watch)
    assert state_of(watch).state == beatjudge.QUIET and state_of(watch).actionable
    assert len(sink.stopped) == 1 and sink.stopped[0] == units.job_id(run_id, 0, 0), "two checks in a row said quiet, then the unit was stopped"
    assert any("is quiet" in line for line in lines(conn, run_id)) and any("is quiet" in message for message in sink.errors)
    assert lines(conn, run_id, "unit-stuck") and runrecord.read_run(conn, run_id)["state"] == "running", "the first stop queues the unit again; it does not fail the run"


# --- the false kill ----------------------------------------------------------------------------------------------

def test_a_counter_moving_at_ten_times_its_estimate_for_six_hours_is_flagged_and_never_touched(conn, clock, run_id):
    """The false kill of 11-08: a unit well over its estimate with the counter still moving is flagged for a person and never stopped, whatever the size."""
    watch, sink = make_watch(conn, run_id, history=beatjudge.History(longest_gap=20.0, expected_seconds=100.0))
    horse = Horse(conn, run_id, clock, estimate=100.0)
    seen = set()
    for _ in range(12):                                # half-hour blocks, 6 hours
        run_for(1800, horse, watch, rate=0.05, step=30.0)       # +90 items every half hour: 1,080 done of 100 expected at the end
        seen.add((state_of(watch).state, state_of(watch).severity))
        assert not state_of(watch).actionable and watch.agreed() == []
    assert horse.done > 10 * horse.estimate
    assert (beatjudge.LOOPING, beatjudge.FLAG) in seen and (beatjudge.LOOPING, beatjudge.SUSPECT) in seen
    assert sink.stopped == [] and watch.stuck.records() == {}
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "running" and not run.get("killed") and not runrecord.is_stale(run)
    assert lines(conn, run_id, "unit-stuck") == [] and lines(conn, run_id, "unit-requeued") == []
    assert len(sink.errors) == 1 and "looping" in sink.errors[0], "the suspect state reaches error.log once for a person to read"


def test_a_slow_unit_is_never_stopped_however_long_it_takes(conn, clock, run_id):
    watch, sink = make_watch(conn, run_id, history=beatjudge.History(longest_gap=60.0, expected_seconds=60.0))
    horse = Horse(conn, run_id, clock, estimate=100_000.0)
    run_for(4 * 3600, horse, watch, rate=0.01, step=30.0)       # 144 of 100,000 in four hours
    assert state_of(watch).state == beatjudge.SLOW and state_of(watch).projected_finish_seconds is not None
    assert sink.stopped == [] and runrecord.read_run(conn, run_id)["state"] == "running"


# --- dead --------------------------------------------------------------------------------------------------------

def test_a_unit_whose_beat_stops_is_dead_after_thirty_seconds_and_stopped_after_two_checks(conn, clock, run_id):
    watch, sink = make_watch(conn, run_id)
    horse = Horse(conn, run_id, clock)
    run_for(120, horse, watch, rate=1.0)
    assert state_of(watch).state == beatjudge.WORKING
    run_for(20, horse, watch, rate=1.0, beating=False, step=10.0)       # the process is gone: nothing is written any more
    assert state_of(watch).state == beatjudge.WORKING, "twenty seconds is two missed beats, not three"
    run_for(60, horse, watch, beating=False)
    assert state_of(watch).state == beatjudge.DEAD and state_of(watch).actionable
    assert len(sink.stopped) == 1, "a dead unit is acted on once two checks agree"
    assert any("is dead" in line for line in lines(conn, run_id)) and any("is dead" in message for message in sink.errors)


# --- stuck twice -----------------------------------------------------------------------------------------------

def _stop_and_requeue(conn, run_id, watch):
    """What the horse does once the stop is sent: the job leaves the queue, is marked stopped, and the next poll queues it again."""
    job = redisqueue.fetch_job(units.job_id(run_id, 0, 0), conn)
    redisqueue.queue(redisqueue.BULK_QUEUE, conn).remove(job.id)
    job.set_status(rq.job.JobStatus.STOPPED)
    watch.check()
    return job


def test_a_unit_stuck_twice_is_queued_again_once_and_then_the_run_fails(conn, clock, run_id):
    watch, sink = make_watch(conn, run_id)
    horse = Horse(conn, run_id, clock)
    run_for(60, horse, watch, rate=1.0)
    run_for(400, horse, watch)                                          # first stuck: two checks in a row after five minutes
    assert len(sink.stopped) == 1 and watch.stuck.records()[0]["attempts"] == 1
    _stop_and_requeue(conn, run_id, watch)
    assert lines(conn, run_id, "unit-requeued") and runrecord.read_run(conn, run_id)["state"] == "running"

    fresh = Horse(conn, run_id, clock)                                  # the new attempt starts afresh, works for a minute and then freezes the same way
    run_for(60, fresh, watch, rate=1.0)
    assert state_of(watch).state == beatjudge.WORKING, "the stopped attempt's history is not held against the new one"
    run_for(400, fresh, watch)
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "failed" and "stuck twice" in run["error"] and "STACK" in run["error"]
    assert len(sink.stopped) == 2 and any("stuck twice" in message for message in sink.errors)


# --- the stale rule ------------------------------------------------------------------------------------------------

def _states_to_run(conn, clock, run_id, conditions):
    """One unit per condition (`working`, `looping`, `waiting`, `quiet`, `dead`) for 800 s of fake time, polled every 30 s; returns the run record the Queue page would judge."""
    watch, _ = make_watch(conn, run_id)
    watch.stuck = None                                                  # judged and published only
    horses = [Horse(conn, run_id, clock, unit=index, estimate=100.0) for index in range(len(conditions))]
    for horse, condition in zip(horses, conditions):
        if condition == "looping":
            horse.done = 500.0
        if condition == "waiting":
            horse.wait("gate")
    for step in range(80):
        clock.now += 10.0
        for horse, condition in zip(horses, conditions):
            horse.tick(10.0, {"working": 0.1, "looping": 2.0}.get(condition, 0.0), beating=condition != "dead" or step < 4)
        if step % 3 == 2:
            watch.check()
    watch.check()
    return runrecord.read_run(conn, run_id)


@pytest.mark.parametrize("conditions,stale", [
    (("quiet",), True),
    (("dead",), True),
    (("quiet", "dead"), True),
    (("quiet", "working"), False),
    (("dead", "working"), False),
    (("waiting",), False),
    (("quiet", "waiting"), False),
    (("looping",), False),
    (("quiet", "looping"), False),
])
def test_a_run_is_stale_only_when_every_running_unit_is_dead_or_quiet(conn, clock, run_id, conditions, stale):
    run = _states_to_run(conn, clock, run_id, conditions)
    judged = sorted(run["stages"][0]["unit_states"].values())
    assert runrecord.is_stale(run, now=clock.now) is stale, f"{conditions}: judged {judged}"


def test_a_run_with_no_unit_judged_keeps_the_old_rule_and_a_long_job_with_a_moving_counter_is_not_stale(conn, clock, run_id):
    run = runrecord.read_run(conn, run_id)
    assert not runrecord.is_stale(run, now=run["started"] + 60.0) and runrecord.is_stale(run, now=run["started"] + 100_000.0), "no unit judged: nothing happened for hours"
    watch, _ = make_watch(conn, run_id)
    horse = Horse(conn, run_id, clock, estimate=100.0)
    run_for(3 * 3600, horse, watch, rate=0.02, step=30.0)
    assert not runrecord.is_stale(runrecord.read_run(conn, run_id), now=clock.now), "hours without a finished unit or an event, but the count moves"


# --- the other job kinds (OPS.142) -----------------------------------------------------------------------------------

def _job_record(conn, kind, job, now, *, age=1.0, still=0.0, done=5.0, estimate=100.0, running=900.0, waiting_on=""):
    conn.hset(jobbeat.hash_key(), f"{kind}:{job}", json.dumps({"kind": kind, "job": job, "label": job, "done": done, "estimate": estimate, "beat": now - age, "started": now - running,
                                                               "moved": now - still, "phase": "step", "waiting_on": waiting_on, "restarts": 0}))


@pytest.fixture
def jobs_clean(conn):
    conn.delete(jobbeat.hash_key())
    yield
    conn.delete(jobbeat.hash_key())


def test_api_web_and_tile_prefetch_jobs_are_judged_by_the_same_rules_on_the_same_clock(conn, jobs_clean):
    now = unitbeat.server_time(conn)
    _job_record(conn, "api", "working", now)
    _job_record(conn, "web", "quiet", now, still=600.0)
    _job_record(conn, "tile-prefetch", "dead", now, age=100.0)
    _job_record(conn, "api", "waiting", now, waiting_on="mysql", still=100.0)
    _job_record(conn, "web", "looping", now, done=1500.0, still=2.0)
    rows = {(row["kind"], row["job"]): row for row in jobbeat.states(conn, now=now)}
    assert rows[("api", "working")]["state"] == beatjudge.WORKING
    assert rows[("web", "quiet")]["state"] == beatjudge.QUIET and "600" in rows[("web", "quiet")]["reason"]
    assert rows[("tile-prefetch", "dead")]["state"] == beatjudge.DEAD
    assert rows[("api", "waiting")]["state"] == beatjudge.WAITING
    assert rows[("web", "looping")]["state"] == beatjudge.LOOPING


def test_a_job_that_is_dead_for_an_hour_leaves_the_list_and_a_working_one_never_does(conn, jobs_clean):
    now = unitbeat.server_time(conn)
    _job_record(conn, "api", "gone", now, age=jobbeat.FORGET_SECONDS + 5.0)
    _job_record(conn, "api", "long", now, running=10 * 3600.0, still=3.0, done=40.0)
    rows = {row["job"]: row for row in jobbeat.states(conn, now=now)}
    assert "gone" not in rows and rows["long"]["state"] == beatjudge.WORKING
    assert conn.hget(jobbeat.hash_key(), "api:gone") is None
