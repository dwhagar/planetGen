# tests/test_beat_judge.py

"""OPS.137: the heartbeat judge turns a unit's beat record, its command's history and a clock into one state; every rule, driven by a fake clock, and a moving counter is never judged for action."""

import pytest

from planetgen.queue import beatjudge
from planetgen.queue.beatjudge import Beat, History, judge

T0 = 1_000_000.0


def beat(now, **fields):
    """A fresh beat at `now` of a unit that started 600 s earlier and has moved its counter a minute ago."""
    values = {"beat": now, "started": now - 600.0, "done": 50.0, "estimate": 100.0, "moved": now - 60.0, "phase": "tier 3"}
    values.update(fields)
    return Beat(**values)


def test_a_unit_with_a_fresh_beat_and_a_moving_counter_is_working():
    verdict = judge(beat(T0), History(), T0)
    assert verdict.state == beatjudge.WORKING and not verdict.actionable and verdict.counter
    assert verdict.rate == pytest.approx(50.0 / 600.0)


def test_no_beat_for_three_intervals_is_dead_and_actionable():
    assert judge(beat(T0 - 29.0), History(), T0).state == beatjudge.WORKING
    verdict = judge(beat(T0 - 31.0), History(), T0)
    assert verdict.state == beatjudge.DEAD and verdict.actionable and "31" in verdict.reason


def test_a_beat_from_the_future_or_not_a_number_is_judged_without_crashing():
    assert judge(beat(T0 + 5.0), History(), T0).state == beatjudge.WORKING
    assert judge(beat(float("nan")), History(), T0).state == beatjudge.DEAD


def test_a_dead_unit_is_dead_even_when_it_was_waiting():
    assert judge(beat(T0 - 60.0, waiting_on="gate"), History(), T0).state == beatjudge.DEAD


def test_waiting_is_not_stuck_and_has_its_own_longer_limit():
    waiting = beat(T0, waiting_on="sky lock", moved=T0 - 1200.0)
    verdict = judge(waiting, History(), T0)
    assert verdict.state == beatjudge.WAITING and not verdict.actionable and not verdict.over_limit
    over = judge(beat(T0, waiting_on="sky lock", moved=T0 - 2400.0), History(), T0)
    assert over.state == beatjudge.WAITING and over.over_limit and not over.actionable and over.severity == beatjudge.FLAG


def test_a_fresh_beat_with_a_frozen_counter_is_quiet_after_five_minutes():
    assert judge(beat(T0, moved=T0 - 299.0), History(), T0).state == beatjudge.WORKING
    verdict = judge(beat(T0, moved=T0 - 301.0), History(), T0)
    assert verdict.state == beatjudge.QUIET and verdict.actionable and verdict.still_seconds == pytest.approx(301.0)


def test_the_quiet_limit_is_five_times_the_longest_gap_of_the_command_when_that_is_longer():
    history = History(longest_gap=240.0)                  # 5 x 240 = 1,200 s
    assert judge(beat(T0, moved=T0 - 1100.0), history, T0).state == beatjudge.WORKING
    assert judge(beat(T0, moved=T0 - 1300.0), history, T0).state == beatjudge.QUIET
    assert beatjudge.quiet_limit_seconds(None) == 300.0 and beatjudge.quiet_limit_seconds(0.0) == 300.0


def test_a_starved_host_doubles_the_quiet_limit():
    frozen = beat(T0, moved=T0 - 400.0)
    assert judge(frozen, History(), T0).state == beatjudge.QUIET
    assert judge(frozen, History(), T0, starved=True).state == beatjudge.WORKING
    assert judge(beat(T0, moved=T0 - 700.0), History(), T0, starved=True).state == beatjudge.QUIET


def test_a_unit_whose_counter_never_moved_is_timed_from_its_start():
    young = beat(T0, moved=None, done=0.0, started=T0 - 200.0)
    assert judge(young, History(), T0).state == beatjudge.WORKING
    assert judge(beat(T0, moved=None, done=0.0, started=T0 - 400.0), History(), T0).state == beatjudge.QUIET


def test_a_unit_with_no_counter_is_working_until_it_is_dead():
    blind = Beat(beat=T0, started=T0 - 3600.0)
    verdict = judge(blind, History(), T0)
    assert verdict.state == beatjudge.WORKING and not verdict.counter and not verdict.actionable
    assert judge(Beat(beat=T0 - 40.0, started=T0 - 3600.0), History(), T0).state == beatjudge.DEAD


def test_a_counter_over_three_times_its_estimate_is_a_looping_flag_and_never_actionable():
    verdict = judge(beat(T0, done=350.0), History(), T0)
    assert verdict.state == beatjudge.LOOPING and verdict.severity == beatjudge.FLAG and not verdict.actionable


def test_a_counter_over_ten_times_its_estimate_is_suspect_and_still_never_judged_for_killing():
    """The false kill of 11-08: a unit counting at ten times its estimate, with the counter moving, must not be acted on."""
    verdict = judge(beat(T0, done=1500.0, moved=T0 - 5.0), History(), T0)
    assert verdict.state == beatjudge.LOOPING and verdict.severity == beatjudge.SUSPECT
    assert not verdict.actionable and not beatjudge.agreed(verdict, verdict)


def test_the_count_going_backward_inside_one_phase_is_looping_but_a_new_phase_may_start_again():
    earlier = beat(T0 - 10.0, done=80.0)
    assert judge(beat(T0, done=20.0), History(), T0, before=earlier).state == beatjudge.LOOPING
    assert judge(beat(T0, done=20.0, phase="writing"), History(), T0, before=earlier).state == beatjudge.WORKING


def test_the_same_phase_beginning_again_too_often_is_looping():
    assert judge(beat(T0, restarts=3), History(), T0).state == beatjudge.WORKING
    assert judge(beat(T0, restarts=4), History(), T0).state == beatjudge.LOOPING


def test_a_rate_under_a_third_of_the_expected_one_is_slow_with_a_projected_finish():
    history = History(expected_seconds=100.0)             # 100 items in 100 s: 1 a second expected
    slow = beat(T0, done=10.0, estimate=100.0, started=T0 - 120.0, moved=T0 - 2.0)    # 0.083 a second
    verdict = judge(slow, history, T0)
    assert verdict.state == beatjudge.SLOW and verdict.severity == beatjudge.FLAG and not verdict.actionable
    assert verdict.projected_finish_seconds == pytest.approx(90.0 / (10.0 / 120.0))


def test_a_rate_near_the_expected_one_or_a_young_unit_is_not_slow():
    history = History(expected_seconds=100.0)
    assert judge(beat(T0, done=60.0, estimate=100.0, started=T0 - 120.0), history, T0).state == beatjudge.WORKING
    assert judge(beat(T0, done=1.0, estimate=100.0, started=T0 - 30.0, moved=T0 - 1.0), history, T0).state == beatjudge.WORKING
    assert judge(beat(T0, done=10.0, started=T0 - 120.0), History(), T0).state == beatjudge.WORKING     # no history, no expected rate


def test_only_dead_and_quiet_are_actionable():
    states = {judge(beat(T0 - 99.0), History(), T0).state: True, judge(beat(T0, moved=T0 - 999.0), History(), T0).state: True}
    assert set(states) == {beatjudge.DEAD, beatjudge.QUIET}
    assert beatjudge.ACTIONABLE == {beatjudge.DEAD, beatjudge.QUIET} and set(beatjudge.STATES) >= beatjudge.ACTIONABLE


def test_two_checks_in_a_row_must_agree_before_anything_is_done():
    quiet = judge(beat(T0, moved=T0 - 999.0), History(), T0)
    dead = judge(beat(T0 - 99.0), History(), T0)
    working = judge(beat(T0), History(), T0)
    assert beatjudge.agreed(quiet, quiet) and beatjudge.agreed(dead, dead)
    assert not beatjudge.agreed(None, quiet) and not beatjudge.agreed(working, quiet) and not beatjudge.agreed(quiet, dead) and not beatjudge.agreed(quiet, working)


def test_a_slow_unit_with_a_moving_counter_for_hours_is_never_acted_on():
    for hours in (1, 3, 6):
        now = T0 + hours * 3600.0
        long_running = Beat(beat=now, started=T0, done=10.0 * hours, estimate=1000.0, moved=now - 20.0, phase="tier 3")
        verdict = judge(long_running, History(expected_seconds=600.0), now)
        assert not verdict.actionable and verdict.state in (beatjudge.WORKING, beatjudge.SLOW)
