"""OPS.192: the second-level count of a unit (`ctx.step`, `ctx.advance`, a `Ticker` of its code) is in its beat, in the verdict the pages read, in the follower's call and on the run's second bar."""

import json

import pytest

from planetgen.generation import followed, run_common
from planetgen.queue import beatjudge, commands, engine, redisqueue, runrecord, unitbeat, units, unitwatch
from planetgen.util import ticker


class _Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture
def run_id(conn):
    return runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=2)])


def test_a_step_is_in_the_record_with_its_count_its_total_and_when_it_began():
    clock = _Clock()
    beat = unitbeat.UnitBeat(None, "r", 0, 0, clock=clock)
    assert "step" not in beat.record(now_server=2000.0)
    clock.now += 3.0
    beat.step("Removing the stored members", 4000)
    clock.now += 2.0
    beat.advance(500)
    beat.advance()
    made = beat.record(now_server=2000.0)
    assert (made["step"], made["step_done"], made["step_total"]) == ("Removing the stored members", 501.0, 4000.0)
    assert made["step_started"] == pytest.approx(2000.0 - 2.0)


def test_a_count_that_moves_starts_the_quiet_clock_again_and_a_new_step_starts_from_nothing():
    clock = _Clock()
    beat = unitbeat.UnitBeat(None, "r", 0, 0, clock=clock)
    beat.step("A", 10)
    clock.now += 100.0
    beat.advance(2)
    assert beat.record(now_server=5000.0)["moved"] == pytest.approx(5000.0)
    clock.now += 10.0
    beat.step("B")
    made = beat.record(now_server=5010.0)
    assert (made["step"], made["step_done"], made["step_total"]) == ("B", 0.0, 0.0)


def test_step_to_sets_the_count_and_the_total_where_a_loop_counts_for_itself():
    beat = unitbeat.UnitBeat(None, "r", 0, 0, clock=_Clock())
    beat.step("Rows", 100)
    beat.step_to(40)
    beat.step_to(60, 200)
    made = beat.record(now_server=1.0)
    assert (made["step_done"], made["step_total"]) == (60.0, 200.0)


def test_nothing_happens_to_a_step_outside_a_unit():
    unitbeat.step("x", 3)
    unitbeat.advance()
    unitbeat.step_to(1)
    assert ticker.sink is None


def test_a_ticker_of_the_units_code_is_its_second_level_count_and_the_sink_is_put_back(conn):
    beat = unitbeat.UnitBeat(conn, "r", 0, 0, clock=_Clock())
    assert ticker.sink is None
    with beat:
        pace = ticker.Ticker("Removing the stored members", 4000, say=lambda text: None, every=3600)
        pace.tick(1000)
        assert beat.record(now_server=1.0)["step"] == "Removing the stored members"
        pace.tick(2500)
        pace.finish(4000)
        made = beat.record(now_server=1.0)
        assert (made["step_done"], made["step_total"]) == (4000.0, 4000.0)
        ticker.Ticker("Rows of the next kind", say=lambda text: None).tick(5)
        assert beat.record(now_server=1.0)["step"] == "Rows of the next kind"
    assert ticker.sink is None


def test_a_ticker_that_changes_its_total_after_it_was_made_gives_the_new_total(conn):
    beat = unitbeat.UnitBeat(conn, "r", 0, 0, clock=_Clock())
    with beat:
        pace = ticker.Ticker("Columns", say=lambda text: None, every=3600)
        pace.tick(0)
        pace.total = 90
        pace.tick(3)
        assert beat.record(now_server=1.0)["step_total"] == 90.0


def test_ctx_step_and_advance_are_the_units_second_level_count(conn):
    beat = unitbeat.UnitBeat(conn, "r", 0, 0, clock=_Clock())
    ctx = engine.Context("k", {}, beat)
    with beat:
        ctx.step("Dropping the indexes", 3)
        ctx.advance()
        ctx.advance(1)
        made = beat.record(now_server=1.0)
    assert (made["step"], made["step_done"], made["step_total"]) == ("Dropping the indexes", 2.0, 3.0)


def test_the_judge_reads_a_step_from_a_stored_record():
    beat = unitbeat.beat_of({"beat": 10.0, "started": 1.0, "step": "Rows", "step_done": 5, "step_total": 20, "step_started": 4.0})
    assert (beat.step, beat.step_done, beat.step_total, beat.step_started) == ("Rows", 5.0, 20.0, 4.0)
    free = unitbeat.beat_of({"beat": 10.0, "started": 1.0})
    assert (free.step, free.step_total, free.step_started) == ("", 0.0, None)


def test_a_units_verdict_carries_its_step_with_its_own_bar_and_none_without_one():
    beat = beatjudge.Beat(beat=100.0, started=0.0, done=1.0, estimate=10.0, step="Rows", step_done=30.0, step_total=60.0, step_started=40.0)
    shown = unitwatch.step_view(beat, 100.0)
    assert (shown["name"], shown["done"], shown["total"]) == ("Rows", 30.0, 60.0)
    assert shown["progress"]["percent"] == 50.0 and shown["progress"]["rate"] == pytest.approx(0.5)
    assert unitwatch.step_view(beatjudge.Beat(beat=1.0, started=0.0), 1.0) is None


def _write(conn, run_id, unit, **record):
    conn.hset(runrecord.progress_key(run_id, 0), unit, json.dumps({"beat": 1.0, "done": 0, "estimate": 0, **record}))


def test_the_stages_count_is_the_first_working_units_step_summed_over_the_units_in_the_same_one(conn, run_id):
    _write(conn, run_id, 1, step="Rows", step_done=10, step_total=100)
    _write(conn, run_id, 0, step="Rows", step_done=5, step_total=50)
    _write(conn, run_id, 2, step="Indexes", step_done=1, step_total=4)
    assert units._step_now(conn, run_id, 0, set()) == {"name": "Rows", "done": 15.0, "total": 150.0}
    assert units._step_now(conn, run_id, 0, {0}) == {"name": "Rows", "done": 10.0, "total": 100.0}
    assert units._step_now(conn, run_id, 0, {0, 1}) == {"name": "Indexes", "done": 1.0, "total": 4.0}
    assert units._step_now(conn, run_id, 0, {0, 1, 2}) is None


def test_a_unit_with_no_step_gives_the_stage_none(conn, run_id):
    _write(conn, run_id, 0, done=3, estimate=9)
    assert units._step_now(conn, run_id, 0, set()) is None


def test_the_follower_hands_the_step_to_its_caller_each_poll_and_none_when_the_stage_ends(conn, run_id):
    _write(conn, run_id, 0, step="Rows", step_done=5, step_total=50)
    seen = []

    def on_result(unit, result):
        pass

    polls = []

    def sleep(seconds):
        polls.append(seconds)
        if len(polls) == 2:          # the unit finishes: its result and its done mark
            conn.hset(runrecord.results_key(run_id, 0), 0, json.dumps({}))
            conn.sadd(runrecord.done_key(run_id, 0), 0)
            conn.hdel(runrecord.progress_key(run_id, 0), 0)

    runrecord.begin_stage(conn, run_id, 0)
    units.follow_stage(conn, run_id, 0, 1, on_result, sleep=sleep, on_step=seen.append)
    assert seen[0] == {"name": "Rows", "done": 5.0, "total": 50.0} and seen[-1] is None


def test_the_runs_second_bar_follows_the_step_and_goes_when_none_is_left():
    with run_common._generation_progress(disable=True) as progress:
        main = progress.add_task("Stage", total=10)
        progress.main_task = main
        progress.show_step({"name": "Rows", "done": 5.0, "total": 50.0})
        detail = progress._detail()
        assert (detail["description"], detail["completed"], detail["total"]) == ("Rows", 5.0, 50.0)
        progress.show_step({"name": "Rows", "done": 20.0, "total": 50.0})
        assert progress._detail()["completed"] == 20.0
        progress.show_step({"name": "Indexes", "done": 1.0, "total": 0.0})
        assert progress._detail()["description"] == "Indexes" and progress._detail()["total"] is None
        progress.show_step(None)
        assert progress._detail() is None and progress.detail_task is None


def test_a_second_bar_another_step_holds_is_left_alone():
    with run_common._generation_progress(disable=True) as progress:
        held = progress.add_task("Shells", total=7)
        progress.detail_task = held
        progress.show_step({"name": "Rows", "done": 1.0, "total": 5.0})
        assert progress.detail_task == held and progress._detail()["description"] == "Shells"
        progress.show_step(None)
        assert progress.detail_task == held


def test_a_follower_draws_the_step_on_the_runs_display_and_one_without_a_display_does_nothing():
    follower = followed.Follower(None)
    follower.on_step({"name": "Rows", "done": 1.0, "total": 5.0})            # no display yet: nothing to draw on
    with run_common._generation_progress(disable=True) as progress:
        follower.progress = progress
        follower.on_step({"name": "Rows", "done": 1.0, "total": 5.0})
        assert progress._detail()["description"] == "Rows"
        follower.on_step(None)
        assert progress._detail() is None


def test_the_queue_page_row_says_the_count_the_rate_and_the_time_left():
    from planetgen.web import queue_page
    step = {"name": "Rows", "done": 30.0, "total": 60.0, "progress": {"show": True, "percent": 50.0, "rate": 0.5, "eta_seconds": 60.0, "source": "live"}}
    shown = queue_page._step_text(step)
    assert shown["name"] == "Rows" and shown["text"].startswith("Rows: 30 of 60, 2.0 s each, ") and shown["text"].endswith("left")
    assert shown["bar"]["percent"] == 50.0
    assert queue_page._step_text(None) is None
    assert queue_page._step_text({"name": "Rows", "done": 7.0, "total": 0.0, "progress": {}})["text"] == "Rows: 7"


def test_wrapper_run_and_follow_name_on_step():
    import inspect
    from planetgen.queue import wrapper
    assert "on_step" in inspect.signature(wrapper.run).parameters and "on_step" in inspect.signature(wrapper.follow).parameters
    assert "on_step" in inspect.signature(units.follow_stage).parameters and commands is not None
