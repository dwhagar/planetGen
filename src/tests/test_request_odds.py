# tests/test_request_odds.py

"""GEN.1130: a request raises the odds of the wanted outcome and never forces it."""

import random

import pytest

from planetgen import tuning
from planetgen.util import draw, request_odds


def test_a_wanted_draw_is_always_kept_and_any_other_one_time_in_the_odds():
    rng = random.Random(4)
    assert all(request_odds.keeps(True, rng) for _ in range(200))
    kept = sum(request_odds.keeps(False, rng) for _ in range(20000))
    assert kept == pytest.approx(20000 / tuning.REQUEST_ODDS, rel=0.1)


def test_a_draw_takes_one_number_of_the_stream_wanted_or_not():
    with draw.bound("odds-a"):
        request_odds.keeps(True)
        after_wanted = draw.random()
    with draw.bound("odds-a"):
        request_odds.keeps(False)
        assert draw.random() == after_wanted


@pytest.mark.parametrize("chance", [0.001, 0.03, 0.2, 0.7])
def test_the_odds_of_a_chance_are_multiplied_and_the_chance_stays_below_one(chance):
    raised = request_odds.raised_chance(chance)
    assert raised / (1.0 - raised) == pytest.approx(tuning.REQUEST_ODDS * chance / (1.0 - chance))
    assert chance < raised < 1.0


def test_a_chance_of_nothing_or_of_everything_does_not_move():
    assert request_odds.raised_chance(0.0) == 0.0 and request_odds.raised_chance(1.0) == 1.0


def test_resampling_by_keeps_gives_the_raised_chance():
    rng = random.Random(9)
    chance, hits, runs = 0.05, 0, 6000
    for _ in range(runs):
        while True:
            met = rng.random() < chance
            if request_odds.keeps(met, rng):
                break
        hits += met
    assert hits / runs == pytest.approx(request_odds.raised_chance(chance), abs=0.03)
