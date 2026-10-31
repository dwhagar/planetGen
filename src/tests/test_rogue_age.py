# tests/test_rogue_age.py

"""
GEN.1112 RP5: a rogue planet's age, ejection delay and time since ejection (`planetgen.physics.rogue_age`): the parent's age from the star-formation `AgeLaw`, the delay's
power-law inverse cumulative, the identity age = delay + time since ejection, the window and the inputs refused.
"""

import math
import random

import pytest

from planetgen.galaxy.age_law import AgeLaw
from planetgen.physics import rogue_age as ra

OLD_HEAVY = AgeLaw((1.0, 10.0), (1.0, 3.0))
YOUNG_HEAVY = AgeLaw((1.0, 10.0), (3.0, 1.0))
FLAT = AgeLaw.uniform((1.0, 10.0))


def test_the_delay_is_log_uniform_between_the_minimum_and_the_parents_age():
    assert ra.ejection_delay_gy(0.0, 5.0) == pytest.approx(ra.EJECTION_DELAY_MIN_GY)
    assert ra.ejection_delay_gy(1.0, 5.0) == pytest.approx(5.0)
    assert ra.ejection_delay_gy(0.5, 5.0) == pytest.approx(math.sqrt(ra.EJECTION_DELAY_MIN_GY * 5.0))
    values = [ra.ejection_delay_gy(i / 100.0, 5.0) for i in range(101)]
    assert values == sorted(values)


def test_an_index_of_zero_is_uniform_in_time_and_two_is_steeper_than_one():
    assert ra.ejection_delay_gy(0.5, 5.0, index=0.0) == pytest.approx(0.5 * (ra.EJECTION_DELAY_MIN_GY + 5.0))
    assert ra.ejection_delay_gy(0.5, 5.0, index=2.0) < ra.ejection_delay_gy(0.5, 5.0, index=1.0) < ra.ejection_delay_gy(0.5, 5.0, index=0.0)


def test_a_parent_younger_than_the_minimum_delay_is_ejected_at_its_age():
    assert ra.ejection_delay_gy(0.3, 0.0005) == 0.0005


def test_an_age_is_delay_plus_time_since_ejection_and_a_delay_never_outlasts_the_age():
    rng = random.Random(3)
    for _ in range(500):
        got = ra.rogue_age(YOUNG_HEAVY, rng.random(), rng.random())
        assert got.age_gy == pytest.approx(got.delay_gy + got.since_ejection_gy)
        assert 1.0 <= got.age_gy <= 10.0 and ra.EJECTION_DELAY_MIN_GY <= got.delay_gy <= got.age_gy
        assert 0.0 <= got.lit_share <= 1.0


def test_the_ages_follow_the_star_formation_history():
    rng = random.Random(5)
    young = sorted(ra.draw_age(YOUNG_HEAVY, rng.random()) for _ in range(4000))
    old = sorted(ra.draw_age(OLD_HEAVY, rng.random()) for _ in range(4000))
    flat = sorted(ra.draw_age(FLAT, rng.random()) for _ in range(4000))
    assert young[2000] < flat[2000] < old[2000]
    assert flat[2000] == pytest.approx(5.5, abs=0.4)


def test_the_window_narrows_the_draw_and_the_minimum_delay_keeps_a_planet_from_being_older_than_nothing():
    assert 3.0 <= ra.draw_age(FLAT, 0.7, window=(3.0, 4.0)) <= 4.0
    law = AgeLaw.uniform((0.0, 1.0))
    assert ra.draw_age(law, 0.0) == pytest.approx(ra.EJECTION_DELAY_MIN_GY)


def test_a_record_gives_the_ages_directly():
    got = ra.from_record(4.0, 0.25)
    assert (got.age_gy, got.delay_gy, got.since_ejection_gy) == (4.0, 0.25, 3.75)
    assert got.lit_share == pytest.approx(0.0625)
    assert ra.from_record(4.0, 4.0).since_ejection_gy == 0.0


def test_most_rogues_were_ejected_early_in_their_parents_life():
    rng = random.Random(11)
    shares = sorted(ra.rogue_age(FLAT, rng.random(), rng.random()).lit_share for _ in range(4000))
    assert shares[2000] < 0.05                       # a log-uniform delay from 1 Myr: the median is a few percent of the age


def test_inputs_the_model_cannot_use_are_refused():
    with pytest.raises(ra.RogueAgeError):
        ra.ejection_delay_gy(1.5, 5.0)
    with pytest.raises(ra.RogueAgeError):
        ra.ejection_delay_gy(0.5, 0.0)
    with pytest.raises(ra.RogueAgeError):
        ra.draw_age(FLAT, -0.1)
    with pytest.raises(ra.RogueAgeError):
        ra.draw_age(FLAT, 0.5, window=(0.0, 0.0005))
    with pytest.raises(ra.RogueAgeError):
        ra.from_record(4.0, 5.0)
    with pytest.raises(ra.RogueAgeError):
        ra.from_record(0.0, 0.0)
