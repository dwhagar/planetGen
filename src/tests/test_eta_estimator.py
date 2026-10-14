# tests/test_eta_estimator.py

"""PERF.33: the time left from the recorded rate blended with the live one."""

import pytest

from planetgen.generation import run_common
from planetgen.generation.stats import Bucket, Estimate, GenerationStats
from planetgen.queue import progress_rate, work
from planetgen.web import generate_page, queue_page


class Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


def test_without_a_recorded_rate_the_live_rate_is_used_as_before():
    assert progress_rate.blended_rate(2.0, 3, None) == 2.0
    assert progress_rate.blended_rate(None, 0, None) is None


def test_with_only_a_recorded_rate_it_is_the_estimate():
    assert progress_rate.blended_rate(None, 0, 4.0) == 4.0
    # The live rate is held back until enough units finished and time passed.
    assert progress_rate.blended_rate(1.0, 4, 4.0, elapsed=100.0) == 4.0
    assert progress_rate.blended_rate(1.0, 10, 4.0, elapsed=5.0) == 4.0


def test_the_live_rate_takes_over_as_units_finish():
    n = 15
    assert progress_rate.blended_rate(1.0, n, 3.0, elapsed=60.0) == pytest.approx(0.5 * 1.0 + 0.5 * 3.0)
    assert progress_rate.blended_rate(1.0, 1500, 3.0, elapsed=60.0) == pytest.approx(1.0, rel=0.02)


def test_the_time_constant_spans_many_tasks_however_slow_they_are():
    assert progress_rate.time_constant(None, 4) == 60.0
    assert progress_rate.time_constant(1.0, 4) == 60.0
    assert progress_rate.time_constant(120.0, 4) == pytest.approx(600.0)
    assert progress_rate.time_constant(120.0, 1) == pytest.approx(2400.0)


def test_the_range_is_the_middle_80_percent_around_the_estimate():
    low, high = progress_rate.eta_range(100.0)
    assert (low, high) == pytest.approx((88.0, 115.0))
    assert progress_rate.eta_range(None) is None


def test_a_bar_with_a_recorded_rate_has_an_eta_before_its_first_unit():
    clock = Clock()
    rate = progress_rate.DecayingRate(clock=clock, prior=2.0)
    assert rate.source == "recorded" and rate.rate == 2.0
    assert rate.eta(100) == pytest.approx(50.0)
    assert progress_rate.DecayingRate(clock=clock).eta(100) is None      # nothing recorded: unknown, as before


def test_a_bar_moves_from_the_recorded_rate_to_its_own_pace():
    clock = Clock()
    rate = progress_rate.DecayingRate(tau=1e9, clock=clock, prior=10.0)
    for _ in range(60):                       # the machine is five times slower than recorded: 2 units a second
        clock.now += 1.0
        rate.add(2)
    assert rate.live == pytest.approx(2.0, rel=0.05)
    weight = 60 / 75
    assert rate.rate == pytest.approx(weight * rate.live + (1 - weight) * 10.0, rel=0.01)
    assert 2.0 < rate.rate < 10.0


def test_a_stalled_bar_holds_its_estimate_rather_than_counting_it_down():
    clock = Clock()
    rate = progress_rate.DecayingRate(clock=clock, prior=1.0)
    clock.now += 30
    assert not rate.stalled()
    assert rate.eta(1000) == pytest.approx(999.0)       # held at the time the rest takes once the next unit is due
    clock.now += 600
    assert rate.stalled()
    assert rate.eta(1000) == pytest.approx(1000.0)


def test_the_pool_rate_and_task_seconds_follow_the_worker_count():
    stats = GenerationStats()
    for workers, seconds in ((1, 1.0), (5, 3.0)):
        bucket = Bucket("scatter", 3, workers=workers, samples=10, seconds_per_task=seconds * 10,
                        seconds_per_system=seconds)
        stats.buckets[("scatter", workers, 3)] = bucket
    assert stats.pool_rate("scatter", 1) == pytest.approx(1.0)
    assert stats.pool_rate("scatter", 5) == pytest.approx(5 / 3.0)
    assert stats.seconds_per_task("scatter", 5) == pytest.approx(30.0)
    assert stats.seconds_per_task("scatter", 3) == pytest.approx(20.0)      # blended between the neighbours
    assert stats.pool_rate("sector", 2) is None and stats.seconds_per_task("sector", 2) is None


def test_a_sector_bar_starts_from_the_recorded_rate_only_when_one_exists():
    measured = Estimate(workers=4, sectors=40, seconds=100.0, measured=True)
    rate, tau = run_common._sector_prior(measured)
    assert rate == pytest.approx(0.4) and tau == pytest.approx(60.0)
    assert run_common._sector_prior(Estimate(workers=4, sectors=40, seconds=100.0, measured=False)) is None
    assert run_common._sector_prior(Estimate(workers=4, sectors=0, seconds=0.0, measured=True)) is None


def test_the_generate_page_never_shows_dashes_for_the_time_left():
    assert generate_page.remaining_label(None, running=True) == "estimating the time left"
    assert generate_page.remaining_label(None, running=False) == ""
    assert generate_page.remaining_label(100.0) == "about 1 m 28 s to 1 m 54 s left"


def test_the_queue_page_shows_a_range_or_estimating():
    assert queue_page._eta_text(None, live=True) == "estimating"
    assert queue_page._eta_text(None, live=False) == ""
    assert queue_page._eta_text(100.0, live=True) == "1 m 28 s to 1 m 54 s"


def _node(done, done_seconds, planned, recorded, workers=2):
    return {"own_counts": {"done": (done, done_seconds)} if done else {}, "tasks_total": planned, "live": True,
            "workers": workers, "started_at": None, "finished_at": None, "children": [],
            "recorded_task_seconds": recorded}


def test_the_queue_blends_the_recorded_task_time_into_the_runs_own_pace():
    # 10 of 110 tasks done at 4 worker-seconds each; earlier runs took 2 s a task.
    node = _node(10, 40.0, 110, recorded=2.0)
    work._roll_up(node)
    n = 10
    rate = progress_rate.blended_rate(10 / 40.0, n, 1 / 2.0)
    assert node["totals"]["eta_seconds"] == pytest.approx(100 / rate / 2)
    own_only = _node(10, 40.0, 110, recorded=None)
    work._roll_up(own_only)
    assert own_only["totals"]["eta_seconds"] == pytest.approx(100 * 4.0 / 2)


def test_the_queue_has_a_time_left_before_the_first_task_with_a_recorded_rate():
    node = _node(0, 0.0, 100, recorded=3.0)
    work._roll_up(node)
    assert node["totals"]["eta_seconds"] == pytest.approx(100 * 3.0 / 2)
    none = _node(0, 0.0, 100, recorded=None)
    work._roll_up(none)
    assert none["totals"]["eta_seconds"] is None


def test_the_sector_bar_carries_the_recorded_rate_into_the_progress_bar():
    class Args:
        _sector_prior = (0.4, 90.0)

    with run_common._generation_progress(disable=True) as progress:
        with run_common.sector_bar(progress, Args(), "Sectors", 40) as bar:
            rate = progress.tasks[0].fields["rate"]
            assert bar.shown and rate.prior == 0.4 and rate.tau == 90.0
            assert rate.eta(40) == pytest.approx(100.0, rel=1e-3)

        class Unrecorded:
            pass

        with run_common.sector_bar(progress, Unrecorded(), "Sectors", 10):
            assert progress.tasks[1].fields["rate"].prior is None


def test_a_job_adds_the_time_of_each_step_not_started_yet():
    node = _node(0, 0.0, 0, recorded=None)
    node["unstarted_step_seconds"] = [60.0, 30.0]
    work._roll_up(node)
    assert node["totals"]["eta_seconds"] == pytest.approx(90.0)
    assert node["totals"]["eta_partial"] is False


def test_a_step_with_no_recorded_time_makes_the_job_estimate_partial():
    node = _node(0, 0.0, 0, recorded=None)
    node["unstarted_step_seconds"] = [60.0, None]
    work._roll_up(node)
    assert node["totals"]["eta_seconds"] == pytest.approx(60.0)
    assert node["totals"]["eta_partial"] is True
    nothing = _node(0, 0.0, 0, recorded=None)
    work._roll_up(nothing)
    assert nothing["totals"]["eta_seconds"] is None
    assert nothing["totals"]["eta_partial"] is False


def _job(**fields):
    base = {"steps": ["a", "b", "c"], "step": 1, "step_started_at": 1000.0, "step_estimates": [100.0, 60.0, 40.0],
            "elapsed_s": 50.0, "finished": False, "status": "running"}
    return {**base, **fields}


def test_a_one_step_job_has_no_overall_bar():
    assert generate_page.overall_view(_job(steps=["a"]), 10.0)["overall_shown"] is False


def test_the_overall_bar_adds_the_later_steps_to_the_running_ones_time_left():
    view = generate_page.overall_view(_job(), 30.0, now=1050.0)
    assert view["overall_shown"] and view["overall_text"] == "Whole job"
    # 30 s left on this step + 60 + 40 for the next two; 50 s already gone.
    assert view["overall_value"] == pytest.approx(50.0 / 180.0)
    assert view["overall_remaining_label"] == generate_page.remaining_label(130.0)


def test_the_overall_bar_falls_back_to_the_recorded_time_of_the_running_step():
    view = generate_page.overall_view(_job(), None, now=1030.0)
    # 100 s recorded less 30 s done, plus the later steps.
    assert view["overall_value"] == pytest.approx(50.0 / (50.0 + 170.0))


def test_an_unrecorded_later_step_is_counted_at_the_average_of_the_finished_ones():
    # Step 2 of 3 running, 20 s into it; step 1 took 50 - 20 = 30 s; step 3 has no record: 30 s stands in for it.
    job = _job(step=2, step_started_at=1030.0, elapsed_s=50.0, step_estimates=[30.0, 60.0, None])
    view = generate_page.overall_view(job, 25.0, now=1050.0)
    assert view["overall_value"] == pytest.approx(50.0 / (50.0 + 25.0 + 30.0))
    assert view["overall_remaining_label"] == generate_page.remaining_label(55.0)


def test_the_first_step_with_no_records_counts_each_later_step_as_long_as_its_own_total():
    job = _job(step_estimates=[None, None, None])
    view = generate_page.overall_view(job, 30.0, now=1050.0)
    # 50 s done + 30 s left = 80 s for this step, and the same for each of the other two.
    assert view["overall_value"] == pytest.approx(50.0 / (50.0 + 30.0 + 80.0 + 80.0))


def test_the_overall_time_left_only_equals_the_running_ones_on_the_last_step():
    job = _job(step=2, step_started_at=1030.0, elapsed_s=50.0, step_estimates=[30.0, None, None])
    middle = generate_page.overall_view(job, 25.0, now=1050.0)
    assert middle["overall_remaining_label"] != generate_page.remaining_label(25.0)
    last = generate_page.overall_view({**job, "step": 3}, 25.0, now=1050.0)
    assert last["overall_remaining_label"] == generate_page.remaining_label(25.0)


def test_the_overall_bar_counts_every_stage_of_the_job_not_every_step():
    stages = [{"n": n, "step": step, "label": f"s{n}", "skipped": reason}
              for n, (step, reason) in enumerate([(1, None), (2, None), (3, None), (3, None), (3, "not asked"), (3, None)], 1)]
    job = _job(steps=["a", "b", "c"], stages=stages, step=3, step_started_at=1000.0, elapsed_s=70.0,
               progress={"stage": {"index": 2, "started_at": 1060.0}}, stage_estimates=[10.0, 10.0, 20.0, 30.0, 0.0, 40.0])
    view = generate_page.overall_view(job, 5.0, now=1070.0)
    # Stage 4 of 6 runs (10 s in, 5 s left); stage 5 is skipped, stage 6 is expected at 40 s.
    assert view["overall_value"] == pytest.approx(70.0 / (70.0 + 5.0 + 40.0))


def test_the_overall_bar_says_estimating_without_any_time_for_the_running_step():
    view = generate_page.overall_view(_job(step_estimates=[]), None, now=1050.0)
    assert view["overall_value"] is None and view["overall_remaining_label"] == "estimating the time left"


def test_the_overall_bar_is_full_for_a_succeeded_job():
    view = generate_page.overall_view(_job(finished=True, status="succeeded"), None)
    assert view["overall_value"] == 1.0 and view["overall_remaining_label"] == ""
