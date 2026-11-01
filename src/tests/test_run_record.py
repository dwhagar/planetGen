"""PERF.96: a run's record in Redis (`queue/runrecord.py`) and the time left in it (`queue/runeta.py`)."""

import json
import socket
import subprocess
import sys
import threading

import pytest

from planetgen.queue import progress_rate, redisqueue, runeta, runrecord
from planetgen.queue.runrecord import Stage
from planetgen.util import logrecord


@pytest.fixture
def conn(redis_server, key_prefix, monkeypatch):
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    return redisqueue.connect(redis_server)


def _stages():
    return [Stage("mass", "Mass scatter", units=4, cost=8.0), Stage("sectors", "Sectors", units=2)]


def test_a_second_run_gets_busy_until_the_first_ends(conn):
    first = runrecord.start_run(conn, "Fill", _stages(), who="boss", now=100)
    assert runrecord.active_run(conn) == first
    with pytest.raises(runrecord.Busy):
        runrecord.start_run(conn, "Another", _stages())
    assert runrecord.finish_run(conn, first, "done", now=200) is True
    assert runrecord.active_run(conn) is None
    assert runrecord.finish_run(conn, first, "failed") is False          # already over: left as it is
    assert runrecord.read_run(conn, first)["state"] == "done"
    assert runrecord.start_run(conn, "Another", _stages())


def test_the_record_reads_back_with_its_stages_and_counters(conn):
    run_id = runrecord.start_run(conn, "Fill shell 1", _stages(), who="boss", database="g1", argv=["galaxy"],
                                 seed=7, now=100)
    runrecord.begin_stage(conn, run_id, 0, now=101)
    run = runrecord.read_run(conn, run_id)
    assert run["title"] == "Fill shell 1" and run["state"] == "running" and run["current"] == 0
    assert [stage["key"] for stage in run["stages"]] == ["mass", "sectors"]
    assert run["stages"][0]["units"] == 4 and run["stages"][0]["cost_total"] == 8.0
    assert run["stages"][0]["started"] == 101 and run["stages"][1]["started"] is None
    assert run["stages"][1]["cost_total"] == 2.0                          # default: one per unit
    assert runrecord.read_run(conn, "20260101-000000-deadbeef") is None
    assert run_id in runrecord.recent_runs(conn)


def test_the_unit_that_completes_a_stage_is_the_only_one_to_advance(conn):
    run_id = runrecord.start_run(conn, "Fill", [Stage("mass", "Mass", units=4)], now=0)
    runrecord.begin_stage(conn, run_id, 0, now=0)
    outcomes = [runrecord.unit_done(conn, run_id, 0, index, seconds=1.0, now=float(index + 1)) for index in (0, 1, 2)]
    assert [o.advance for o in outcomes] == [False, False, False] and outcomes[-1].done == 3
    again = runrecord.unit_done(conn, run_id, 0, 1, now=9.0)               # a retried unit runs twice
    assert not again.counted and not again.advance and again.done == 3
    last = runrecord.unit_done(conn, run_id, 0, 3, now=5.0)
    assert last.counted and last.advance and last.done == 4
    assert not runrecord.unit_done(conn, run_id, 0, 3, now=6.0).advance    # and never twice
    stage = runrecord.read_run(conn, run_id)["stages"][0]
    assert stage["done"] == 4 and stage["cost_done"] == 4.0 and stage["ended"] == 5.0


def test_units_finishing_at_once_on_many_threads_advance_exactly_one(conn):
    run_id = runrecord.start_run(conn, "Fill", [Stage("mass", "Mass", units=40)], now=0)
    runrecord.begin_stage(conn, run_id, 0, now=0)
    advanced, lock = [], threading.Lock()

    def finish(first):
        for unit in range(first, 40, 8):
            if runrecord.unit_done(conn, run_id, 0, unit, now=1.0).advance:
                with lock:
                    advanced.append(unit)

    threads = [threading.Thread(target=finish, args=(n,)) for n in range(8)]
    [t.start() for t in threads]
    [t.join() for t in threads]
    assert len(advanced) == 1
    assert runrecord.read_run(conn, run_id)["stages"][0]["done"] == 40


def test_the_live_rate_in_redis_is_the_terminal_bars(conn):
    """The script's decayed sums are `progress_rate.DecayingRate`'s, so the page and the bar agree."""
    run_id = runrecord.start_run(conn, "Fill", [Stage("mass", "Mass", units=30)], now=1000.0)
    runrecord.begin_stage(conn, run_id, 0, now=1000.0)
    clock = [1000.0]
    reference = progress_rate.DecayingRate(clock=lambda: clock[0])
    for unit, at in enumerate([1002.0, 1002.0, 1005.5, 1006.0, 1020.0, 1020.5, 1021.0, 1060.0, 1075.0, 1075.0]):
        clock[0] = at
        runrecord.unit_done(conn, run_id, 0, unit, now=at)
        reference.add(1, now=at)
        stage = runrecord.read_run(conn, run_id)["stages"][0]
        assert stage["u"] / stage["t"] == pytest.approx(reference.live, rel=1e-9)


def test_the_event_stream_is_the_log_and_can_be_followed_and_kept(conn, tmp_path):
    run_id = runrecord.start_run(conn, "Fill", [Stage("mass", "Mass", units=1)], now=100.0)
    runrecord.log(conn, run_id, "scatter layer 3 saved", now=101.0, stage=0)
    runrecord.log(conn, run_id, "disk is slow", level="warning", now=102.0)
    events = runrecord.read_events(conn, run_id)
    assert [fields["message"] for _id, fields in events] == [
        "run started: Fill", "scatter layer 3 saved", "disk is slow"]
    assert runrecord.read_events(conn, run_id, after=events[1][0])[0][1]["message"] == "disk is slow"
    assert runrecord.read_events(conn, run_id, after=events[-1][0]) == []
    runrecord.finish_run(conn, run_id, "done", now=103.0)
    path = runrecord.write_log_file(conn, run_id, tmp_path / "jobs")
    text = open(path, encoding="utf-8").read()
    assert 'level=warning component=run msg="disk is slow"' in text and text.rstrip().endswith('msg="run done" event=run.done')
    assert text.splitlines()[0] == 'ts=1970-01-01T00:01:40.000Z level=info component=run msg="run started: Fill" event=run.start'
    assert [logrecord.parse_file_line(line).message for line in text.splitlines()] == [
        "run started: Fill", "scatter layer 3 saved", "disk is slow", "run done"]
    for n in range(25):
        (tmp_path / "jobs" / f"2020010{n % 10}-{n:06d}-aaaaaaaa.log").write_text("old")
    runrecord.write_log_file(conn, run_id, tmp_path / "jobs", keep=20)
    assert len(list((tmp_path / "jobs").glob("*.log"))) == 20


def test_cancel_marks_the_run_and_a_finished_run_expires_later(conn):
    run_id = runrecord.start_run(conn, "Fill", _stages(), now=0.0)
    assert not runrecord.is_cancelling(conn, run_id)
    assert runrecord.request_cancel(conn, run_id, now=1.0) and runrecord.is_cancelling(conn, run_id)
    assert runrecord.read_run(conn, run_id)["state"] == "cancelling"
    assert runrecord.finish_run(conn, run_id, "cancelled", now=2.0)
    assert not runrecord.request_cancel(conn, run_id)
    assert 0 < conn.ttl(runrecord.run_key(run_id)) <= runrecord.KEEP_SECONDS
    assert conn.ttl(runrecord.stage_key(run_id, 0)) > 0 and conn.ttl(runrecord.events_key(run_id)) > 0
    with pytest.raises(ValueError):
        runrecord.finish_run(conn, run_id, "weird")


def test_a_run_is_stale_when_nothing_has_happened_for_a_while(conn):
    run_id = runrecord.start_run(conn, "Fill", _stages(), now=1000.0)
    run = runrecord.read_run(conn, run_id)
    assert not runrecord.is_stale(run, now=1100.0) and runrecord.is_stale(run, now=1200.0)
    runrecord.begin_stage(conn, run_id, 0, now=1150.0)
    assert not runrecord.is_stale(runrecord.read_run(conn, run_id), now=1200.0)
    runrecord.finish_run(conn, run_id, "done", now=1300.0)
    assert not runrecord.is_stale(runrecord.read_run(conn, run_id), now=9999.0)


# --- the time left (pure) -----------------------------------------------------------------------------------

def _stage(**fields):
    base = {"key": "s", "label": "S", "units": 10, "cost": 10.0, "cost_total": 10.0, "done": 0, "cost_done": 0.0,
            "expected_seconds": None, "prior_rate": None, "worker_rate": None, "workers": 1, "started": None, "ended": None,
            "last": None, "last_end": None, "first_end": None, "u": None, "t": None, "longest_gap": None,
            "seconds": 0.0}
    base.update(fields)
    return base


def _run(*stages, state="running", started=0.0, current=0):
    return {"id": "r", "title": "Fill", "who": "", "database": "", "state": state, "error": None, "current": current,
            "created": started, "started": started, "updated": started, "finished": None, "stages": list(stages)}


def test_a_steady_run_estimates_its_time_left_from_the_live_rate():
    # 4 of 10 units in 40 s: 0.1 units a second, 6 left, so 60 s.
    stage = _stage(started=0.0, done=4, cost_done=4.0, u=4.0 * 0.6, t=40 * 0.6, last_end=40.0, longest_gap=10.0)
    eta = runeta.estimate(_run(stage), now=40.0)
    assert eta.seconds == pytest.approx(60.0, rel=0.02) and eta.source == "live" and not eta.from_earlier_runs
    assert eta.low == pytest.approx(eta.seconds * 0.88) and eta.high == pytest.approx(eta.seconds * 1.15)
    counted_down = runeta.estimate(_run(stage), now=45.0).seconds
    assert counted_down == pytest.approx(eta.seconds - 5.0, rel=0.02)


def test_before_any_unit_reports_the_recorded_rate_or_earlier_runs_carry_the_estimate():
    recorded = _stage(started=0.0, prior_rate=0.5)                         # 10 cost left at 0.5 a second
    eta = runeta.estimate(_run(recorded), now=3.0)
    assert eta.seconds == pytest.approx(20.0) and eta.source == "recorded" and eta.from_earlier_runs
    earlier = _stage(started=0.0, expected_seconds=300.0)                  # a phase that reports nothing
    eta = runeta.estimate(_run(earlier), now=100.0)
    assert eta.seconds == pytest.approx(200.0) and eta.source == "earlier runs" and eta.from_earlier_runs
    assert runeta.estimate(_run(earlier), now=1000.0).seconds == pytest.approx(30.0)     # over: a tenth is left


def test_a_phase_with_no_updates_and_no_history_uses_the_time_per_task_so_far_in_the_run():
    """Boss, 2026-10-15 05:30Z: the whole bar still shows a time remaining."""
    done = _stage(started=0.0, ended=100.0, done=10, cost_done=10.0, cost_total=10.0, u=1.0, t=1.0)   # 0.1/s so far
    silent = _stage(started=100.0, units=1, cost=5.0, cost_total=5.0)
    eta = runeta.estimate(_run(done, silent, current=1), now=120.0)
    assert eta.source == "this run so far" and not eta.partial
    assert eta.seconds == pytest.approx(5.0 / (10.0 / 120.0))
    nothing = runeta.estimate(_run(_stage(started=0.0), started=0.0), now=5.0)
    assert nothing.rough and nothing.source == "not timed yet" and not nothing.partial       # PERF.184: a guess, not a blank


def test_a_server_that_never_timed_a_stage_estimates_from_its_hardware_and_says_so(monkeypatch):
    """PERF.184 (Boss, 2026-10-24): no benchmark or history, the estimate must still be sensible."""
    monkeypatch.setattr(runeta.roughtime, "host_workers", lambda cores=None, memory_bytes=None: 8)
    started = _stage(started=0.0, units=16)                      # 16 units over 8 workers: 2 rounds of 20 s
    eta = runeta.estimate(_run(started), now=0.0)
    assert eta.rough and eta.source == "not timed yet" and eta.seconds == pytest.approx(40.0)
    assert (eta.low, eta.high) == pytest.approx((16.0, 120.0)), "a guess is wide"
    counted_down = runeta.estimate(_run(started), now=30.0)
    assert counted_down.seconds == pytest.approx(10.0)
    over = runeta.estimate(_run(started), now=500.0)
    assert over.seconds == pytest.approx(4.0), "past the guess, a tenth of it is left"
    later = runeta.estimate(_run(started, _stage(units=0, cost=0.0, cost_total=0.0)), now=0.0)
    assert later.seconds == pytest.approx(40.0 + 600.0 / 8), "a stage not expanded yet counts as a typical one"


def test_the_guess_gives_way_to_the_runs_own_pace_and_to_recorded_timings(monkeypatch):
    monkeypatch.setattr(runeta.roughtime, "host_workers", lambda cores=None, memory_bytes=None: 8)
    running = _stage(started=0.0, done=4, cost_done=4.0, u=4.0 * 0.6, t=40 * 0.6, last_end=40.0, longest_gap=10.0)
    eta = runeta.estimate(_run(running), now=40.0)
    assert eta.source == "live" and not eta.rough and eta.seconds == pytest.approx(60.0, rel=0.02)
    recorded = _stage(started=0.0, units=16, cost=16.0, cost_total=16.0, worker_rate=0.1)   # 0.1 cost per worker second
    eta = runeta.estimate(_run(recorded), now=0.0)
    assert eta.source == "recorded" and not eta.rough and eta.seconds == pytest.approx(16.0 / (0.1 * 8))
    earlier = runeta.estimate(_run(_stage(started=0.0, expected_seconds=300.0)), now=0.0)
    assert earlier.source == "earlier runs" and not earlier.rough


def test_later_stages_are_added_and_a_stage_with_no_history_uses_the_pace_of_the_run():
    now_stage = _stage(started=0.0, done=5, cost_done=5.0, u=5 * 0.6, t=50 * 0.6, last_end=50.0, longest_gap=10.0)
    later = _stage(expected_seconds=120.0)
    eta = runeta.estimate(_run(now_stage, later), now=50.0)
    assert eta.seconds == pytest.approx(50.0 + 120.0, rel=0.02) and not eta.partial
    unknown = runeta.estimate(_run(now_stage, _stage()), now=50.0)          # no history: this run's own pace (0.1/s)
    assert not unknown.partial and unknown.seconds == pytest.approx(50.0 + 100.0, rel=0.02)


def test_a_finished_run_has_none_left_and_a_stall_holds_the_estimate():
    done = _stage(started=0.0, ended=10.0, done=10, cost_done=10.0)
    assert runeta.estimate(_run(done, state="done"), now=500.0).seconds == 0.0
    assert runeta.estimate(_run(done), now=500.0).seconds == 0.0               # every stage over
    live = _stage(started=0.0, done=4, cost_done=4.0, u=4 * 0.6, t=40 * 0.6, last_end=40.0, longest_gap=10.0)
    held = runeta.estimate(_run(live), now=500.0)
    assert held.stalled and held.seconds == pytest.approx(60.0, rel=0.02)         # not counted down to zero


def test_the_public_view_has_what_a_banner_needs_and_nothing_private():
    stage = _stage(started=0.0, done=4, cost_done=4.0, u=4 * 0.6, t=40 * 0.6, last_end=40.0, longest_gap=10.0,
                   label="Mass scatter")
    run = _run(stage, _stage(label="Sectors", expected_seconds=60.0))
    run.update(who="boss", database="g1", updated=40.0)
    view = runeta.public_view(run, now=40.0)
    assert set(view) == {"active", "state", "title", "stage", "percent", "eta_seconds", "eta_low", "eta_high",
                         "eta_source", "from_earlier_runs", "partial", "stalled", "eta_rough", "updated"}
    assert view["active"] and view["stage"] == {"index": 0, "count": 2, "label": "Mass scatter", "done": 4, "units": 10}
    assert view["percent"] == 20.0 and view["eta_source"] == "live"
    assert "boss" not in json.dumps(view) and "g1" not in json.dumps(view)
    assert runeta.public_view(_run(stage, state="done"), now=1.0) == {"active": False}
    assert runeta.public_view(None, now=1.0) == {"active": False}


# --- a run whose process is gone does not keep the machine -----------------------------------------------------

def _dead_pid():
    child = subprocess.Popen([sys.executable, "-c", "pass"])
    child.wait()
    return child.pid


def test_a_run_whose_process_was_killed_no_longer_makes_the_next_run_busy(conn):
    first = runrecord.start_run(conn, "Fill", _stages(), now=100)
    conn.hset(runrecord.run_key(first), "owner", f"{socket.gethostname()}:{_dead_pid()}")     # kill -9, a restart
    second = runrecord.start_run(conn, "Another", _stages(), now=200)
    assert runrecord.active_run(conn) == second
    old = runrecord.read_run(conn, first)
    assert old["state"] == "failed" and "was gone" in old["error"]
    assert runrecord.is_cancelling(conn, first) is False and runrecord.finish_run(conn, first, "done") is False


def test_a_run_whose_process_lives_keeps_the_machine_even_when_its_sign_of_life_lapsed(conn):
    first = runrecord.start_run(conn, "Fill", _stages(), now=100)
    conn.delete(runrecord.alive_key(first))
    with pytest.raises(runrecord.Busy):
        runrecord.start_run(conn, "Another", _stages())
    assert runrecord.read_run(conn, first)["state"] == "running"


def test_a_run_of_another_host_is_dead_when_its_sign_of_life_has_run_out(conn):
    first = runrecord.start_run(conn, "Fill", _stages(), now=100)
    conn.hset(runrecord.run_key(first), "owner", "some-other-host:4242")
    with pytest.raises(runrecord.Busy):
        runrecord.start_run(conn, "Another", _stages())
    conn.delete(runrecord.alive_key(first))
    assert runrecord.start_run(conn, "Another", _stages())
    assert runrecord.read_run(conn, first)["state"] == "failed"


def test_a_run_recorded_before_runs_said_who_follows_them_is_dead_without_a_sign_of_life(conn):
    first = runrecord.start_run(conn, "Fill", _stages(), now=100)
    conn.hdel(runrecord.run_key(first), "owner")
    conn.delete(runrecord.alive_key(first))
    assert runrecord.start_run(conn, "Another", _stages())


def test_an_active_key_with_no_run_behind_it_is_freed(conn):
    conn.set(runrecord.active_key(), "20260101-000000-deadbeef")
    second = runrecord.start_run(conn, "Another", _stages())
    assert runrecord.active_run(conn) == second


def test_a_refused_run_leaves_no_sign_of_life_and_a_live_run_keeps_its_own(conn):
    first = runrecord.start_run(conn, "Fill", _stages())
    with pytest.raises(runrecord.Busy):
        runrecord.start_run(conn, "Another", _stages(), run_id="20260101-000000-feedface")
    assert not conn.exists(runrecord.alive_key("20260101-000000-feedface"))
    assert 0 < conn.ttl(runrecord.alive_key(first)) <= runrecord.ALIVE_SECONDS


def test_a_finished_stage_is_remembered_and_the_next_run_starts_from_it(conn):
    """PERF.184: from its second run on, a stage with no stored history has a time and a rate."""
    from planetgen.queue import stagetimes
    stage = Stage("warm-demo", "Warm demo", units=2, cost=2.0)
    first = runrecord.start_run(conn, "Run one", [stage], now=100.0)
    runrecord.begin_stage(conn, first, 0, now=100.0)
    runrecord.unit_done(conn, first, 0, 0, seconds=30.0, now=110.0)
    assert stagetimes.lookup(conn, ["warm-demo"]) == {}, "nothing is remembered before the stage ends"
    runrecord.unit_done(conn, first, 0, 1, seconds=30.0, now=130.0)
    found = stagetimes.lookup(conn, ["warm-demo", "never-run"])
    assert set(found) == {"warm-demo"} and found["warm-demo"]["seconds"] == pytest.approx(30.0)
    assert found["warm-demo"]["worker_rate"] == pytest.approx(2.0 / 60.0) and found["warm-demo"]["runs"] == 1
    runrecord.finish_run(conn, first, "done", now=131.0)
    second = runrecord.start_run(conn, "Run two", [Stage("warm-demo", "Warm demo", units=2, cost=2.0)], now=200.0)
    runrecord.begin_stage(conn, second, 0, now=200.0)
    runrecord.unit_done(conn, second, 0, 0, seconds=10.0, now=205.0)
    runrecord.unit_done(conn, second, 0, 1, seconds=10.0, now=210.0)
    again = stagetimes.lookup(conn, ["warm-demo"])["warm-demo"]
    assert again["runs"] == 2 and again["seconds"] == pytest.approx(0.5 * 30.0 + 0.5 * 10.0)


def test_a_stage_that_ended_at_once_is_not_remembered(conn):
    from planetgen.queue import stagetimes
    run_id = runrecord.start_run(conn, "Quick", [Stage("quick", "Quick", units=1)], now=100.0)
    runrecord.begin_stage(conn, run_id, 0, now=100.0)
    runrecord.unit_done(conn, run_id, 0, 0, seconds=0.01, now=100.1)
    assert stagetimes.lookup(conn, ["quick"]) == {}


def test_a_stage_keeps_the_largest_worker_memory_and_its_finish_line_says_it(conn):
    """PERF.226: each unit reports the most resident memory its process held; the stage keeps the largest, the finish line names it, a stage
    whose units reported none says nothing of it."""
    run_id = runrecord.start_run(conn, "Fill", _stages(), who="boss", now=100)
    runrecord.begin_stage(conn, run_id, 0, now=101)
    for unit, peak in enumerate((812.5, 1487.0, 640.0, 900.0)):
        runrecord.unit_done(conn, run_id, 0, unit, seconds=1.0, now=float(110 + unit), peak_mb=peak)
    run = runrecord.read_run(conn, run_id)
    assert run["stages"][0]["peak_mb"] == 1487.0
    finish = [message for _id, fields in runrecord.read_events(conn, run_id) if fields.get("kind") == logrecord.STAGE_FINISH
              for message in [fields["message"]]]
    assert finish and finish[0].endswith("(largest worker 1487 MB)")
    runrecord.begin_stage(conn, run_id, 1, now=200)
    runrecord.unit_done(conn, run_id, 1, 0, now=201.0)
    runrecord.unit_done(conn, run_id, 1, 1, now=202.0)
    last = [fields["message"] for _id, fields in runrecord.read_events(conn, run_id) if fields.get("kind") == logrecord.STAGE_FINISH][-1]
    assert "largest worker" not in last


def test_a_refusal_names_the_run_that_holds_the_machine(conn):
    first = runrecord.start_run(conn, "Check the data sets", _stages(), kind=runrecord.CHECK, now=100)
    with pytest.raises(runrecord.Busy) as refused:
        runrecord.start_run(conn, "Another", _stages(), kind=runrecord.CHECK)
    text = str(refused.value)
    assert f"run {first} is running" in text and "Check the data sets" in text and "check kind" in text
    assert "1970-01-01 00:01:40 UTC" in text and runrecord._owner() in text
