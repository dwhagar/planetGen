# src/tests/test_task_progress.py

"""PERF.147: a task's own bar and ETA (`queue/taskprogress.py`), driven by fake clocks."""

import pytest

from planetgen.queue import taskprogress as tp


def test_a_task_that_ends_in_under_fifteen_seconds_shows_no_bar():
    for ran, done in ((0.0, 0), (5.0, 50), (15.0, 99)):
        assert not tp.progress(done, 100, 0.0, ran).show
    assert not tp.progress(100, 100, 0.0, 12.0).show, "a task that finished quickly never gained a bar"


def test_a_task_expected_to_take_longer_shows_its_bar_from_the_start_with_the_stored_rate_eta():
    task = tp.progress(0, 1000, 0.0, 0.0, stored_rate=10.0)       # 100 s expected
    assert task.show and task.source == tp.STORED and task.eta_seconds == 100.0 and task.expected_seconds == 100.0
    assert task.fraction == 0.0
    assert not tp.progress(0, 100, 0.0, 0.0, stored_rate=10.0).show, "10 s expected: no bar yet"


def test_a_task_that_overruns_fifteen_seconds_gains_the_bar_when_it_passes_it():
    assert not tp.progress(10, 100, 0.0, 15.0, stored_rate=10.0).show     # 10 s expected, 15 s in: not over yet
    assert tp.progress(10, 100, 0.0, 15.5, stored_rate=10.0).show
    assert tp.progress(10, 100, 0.0, 15.5).show, "with no recorded rate the bar comes when the task has run that long"


def test_without_a_recorded_rate_the_eta_is_the_running_rate_after_the_hold():
    early = tp.progress(5, 100, 0.0, 10.0)
    assert early.show is False and early.rate == 0.5 and early.source == tp.LIVE and early.eta_seconds == 190.0
    late = tp.progress(40, 100, 0.0, 40.0)
    assert late.source == tp.LIVE and late.rate == 1.0 and late.eta_seconds == 60.0


def test_the_eta_moves_from_the_stored_rate_to_the_running_rate_as_the_task_proceeds():
    stored = 10.0                                  # recorded: 10 items a second, 1000 items = 100 s
    # the machine is busy: the task really runs at 2 items a second
    etas, weights = [], []
    for done in (0, 20, 100, 200, 500):
        now = done / 2.0
        task = tp.progress(done, 1000, 0.0, now, stored_rate=stored)
        etas.append(task.eta_seconds)
        weights.append(task.weight)
    # held on the stored rate until HOLD_SECONDS (20 s) has passed: 100 items at 2 a second is 50 s in
    assert weights[0] == 0.0 and weights[1] == 0.0                     # 0 s and 10 s in
    assert weights[2] == pytest.approx(0.4)                            # 10% done of the 25% it takes to trust the running rate fully
    assert weights[3] == pytest.approx(0.8) and weights[4] == 1.0
    assert etas[1] == pytest.approx(98.0)                              # still the stored rate: 980 items at 10 a second
    assert etas[2] == pytest.approx(900 / (0.4 * 2.0 + 0.6 * 10.0))    # 900 items at the blended 6.8 a second
    assert etas[4] == pytest.approx(250.0)                             # the running rate alone: 500 left at 2 a second
    totals = [done / 2.0 + eta for done, eta in zip((0, 20, 100, 200, 500), etas)]
    assert totals == sorted(totals) and totals[0] == 100.0 and totals[-1] == pytest.approx(500.0)       # the guess of the whole time moves from the recorded 100 s to the real 500 s


def test_a_changed_rate_is_followed_in_both_directions():
    faster = tp.progress(500, 1000, 0.0, 25.0, stored_rate=2.0)        # recorded 2 a second, running 20 a second
    assert faster.source == tp.LIVE and faster.eta_seconds == pytest.approx(25.0)
    slower = tp.progress(100, 1000, 0.0, 100.0, stored_rate=20.0)      # recorded 20 a second, running 1 a second at 10% done
    assert slower.source == tp.BLENDED and slower.rate == pytest.approx(0.4 * 1.0 + 0.6 * 20.0)


def test_a_task_with_no_counter_has_no_percent_and_no_eta():
    task = tp.progress(0, 0, 0.0, 30.0)
    assert task.show and task.fraction is None and task.eta_seconds is None
    assert tp.progress(None, None, 0.0, 3.0) == tp.TaskProgress(show=False)


def test_a_finished_task_has_nothing_left_and_bad_numbers_are_ignored():
    assert tp.progress(100, 100, 0.0, 60.0).eta_seconds == 0.0
    assert tp.progress(250, 100, 0.0, 60.0).fraction == 1.0
    assert tp.progress(float("nan"), 100, 0.0, 60.0, stored_rate=float("inf")).eta_seconds is None
    assert tp.progress(-5, 100, 0.0, 60.0, stored_rate=0.0).fraction == 0.0


def test_the_dict_for_the_api_has_a_percent():
    shown = tp.as_dict(tp.progress(250, 1000, 0.0, 50.0, stored_rate=10.0))
    assert shown["percent"] == 25.0 and shown["show"] is True and set(shown) == {"show", "percent", "eta_seconds", "source", "rate", "expected_seconds"}


def test_the_bar_macro_shows_only_a_bar_that_is_to_be_shown(web_app=None):
    from jinja2 import Environment, FileSystemLoader
    import planetgen.web
    import os
    templates = os.path.join(os.path.dirname(planetgen.web.__file__), "templates")
    env = Environment(loader=FileSystemLoader(templates))
    macro = env.get_template("partials/task_bar.html").module.task_bar
    assert macro(None) == "" and macro({"show": False, "percent": 10.0, "eta_seconds": 5.0, "source": "live"}) == ""
    shown = macro(tp.as_dict(tp.progress(250, 1000, 0.0, 50.0, stored_rate=10.0)), "Unit 3")
    assert '<progress max="100" value="25.0" aria-label="Unit 3 progress">' in shown and "25%" in shown and "left" in shown
    assert "(from earlier runs)" in macro(tp.as_dict(tp.progress(0, 1000, 0.0, 0.0, stored_rate=10.0)))
