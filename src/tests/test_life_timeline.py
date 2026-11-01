# tests/test_life_timeline.py

"""GEN.750: `physics.life_timeline`, the five hard steps of life, their temperature clock and the oxygen the third one waits for."""

import math

import pytest

from planetgen import tuning
from planetgen.physics import life_timeline as timeline
from planetgen.util import draw, request_odds

UNIT = (1.0,) * timeline.STEPS


def test_the_reference_temperature_is_the_reference_pace():
    assert timeline.rate_factor(tuning.LIFE_REFERENCE_TEMPERATURE_K) == pytest.approx(1.0)
    assert timeline.rate_factor(None) == 1.0 and timeline.rate_factor(0.0) == 1.0


def test_a_warm_world_goes_faster_and_a_cold_one_slower_inside_the_range():
    low, high = tuning.LIFE_RATE_FACTOR_RANGE
    factors = [timeline.rate_factor(t) for t in (150.0, 250.0, 288.0, 300.0, 330.0, 900.0)]
    assert factors == sorted(factors)
    assert factors[0] == low and factors[-1] == high
    expected = math.exp(tuning.LIFE_ACTIVATION_ENERGY_EV * 1.602176634e-19 / 1.380649e-23 * (1 / 288.0 - 1 / 300.0))
    assert timeline.rate_factor(300.0) == pytest.approx(expected, rel=1e-3)


def test_the_pace_is_a_word_for_the_factor():
    assert timeline.pace_of(1.0) == "normal"
    assert timeline.pace_of(tuning.LIFE_PACE_FAST_ABOVE + 0.01) == "fast"
    assert timeline.pace_of(tuning.LIFE_PACE_SLOW_BELOW - 0.01) == "slow"


def test_the_unit_waits_are_exponential_with_mean_one_and_come_from_the_run():
    with draw.bound("timeline-waits"):
        waits = [timeline.unit_waits() for _ in range(4000)]
    assert all(w > 0.0 for row in waits for w in row)
    assert sum(w for row in waits for w in row) / (4000 * timeline.STEPS) == pytest.approx(1.0, abs=0.05)
    with draw.bound("timeline-waits"):
        assert timeline.unit_waits() == waits[0]


def test_each_step_waits_its_own_mean_scaled_and_sped_by_the_factor():
    ages = timeline.schedule(UNIT, 1.0, 0.0)
    gaps = [b - a for a, b in zip((tuning.LIFE_START_AGE_GYR, *ages), ages)]
    assert gaps == pytest.approx([mean * tuning.LIFE_STEP_SCALE for mean in tuning.LIFE_STEP_MEANS_GYR])
    fast = timeline.schedule(UNIT, 2.0, 0.0)
    assert (fast[-1] - tuning.LIFE_START_AGE_GYR) == pytest.approx((ages[-1] - tuning.LIFE_START_AGE_GYR) / 2.0)


def test_complex_cells_wait_for_the_oxygen_after_the_first_photosynthesis():
    now = timeline.schedule(UNIT, 1.0, 0.0)
    later = timeline.schedule(UNIT, 1.0, 3.0)
    assert later[:2] == now[:2]
    step = tuning.LIFE_STEP_MEANS_GYR[timeline.OXYGEN_STEP] * tuning.LIFE_STEP_SCALE
    assert later[timeline.OXYGEN_STEP] == pytest.approx(later[1] + 3.0 + step)
    assert later[-1] - now[-1] == pytest.approx(3.0)


def test_an_air_with_oxygen_from_the_start_does_not_wait():
    assert timeline.schedule(UNIT, 1.0, 5.0, oxygen_at_start=True) == timeline.schedule(UNIT, 1.0, 0.0)


def test_a_world_that_never_makes_oxygen_stops_at_photosynthesis():
    ages = timeline.schedule(UNIT, 1.0, None)
    assert ages[2:] == (None,) * (timeline.STEPS - 2) and None not in ages[:2]


def test_a_requested_civilization_moves_the_run_by_one_draw():
    with draw.bound("request-a"):
        timeline.raised_waits(4.0, 1.0, 1.0)
        after = draw.random()
    with draw.bound("request-a"):
        draw.random()
        assert draw.random() == after, "one draw of the run's stream however many tries it took"


def test_a_request_raises_the_odds_of_a_timeline_that_finishes_and_forces_none():
    """The odds of finishing inside the age are multiplied by `tuning.REQUEST_ODDS` (a chance of c becomes `raised_chance(c)`), and some timelines
    still do not finish: nothing is squeezed to fit."""
    age, factor, delay = 4.0, 1.0, 1.0
    runs = 1500
    with draw.bound("request-b"):
        plain = sum(timeline.schedule(timeline.unit_waits(), factor, delay)[-1] <= age for _ in range(runs)) / runs
        raised = [timeline.schedule(timeline.raised_waits(age, factor, delay), factor, delay)[-1] for _ in range(runs)]
    share = sum(last is not None and last <= age for last in raised) / runs
    expected = request_odds.raised_chance(plain)
    assert 4 * plain < share
    assert share == pytest.approx(expected, abs=0.04)
    assert share < 1.0


def test_a_request_on_a_world_that_cannot_fit_a_civilization_leaves_it_without_one():
    with draw.bound("request-c"):
        waits = timeline.raised_waits(1.0, 0.25, 0.5, tries=3)
        assert timeline.schedule(waits, 0.25, 0.5)[-1] > 1.0
        for _ in range(30):                        # it never makes oxygen: no timeline reaches a civilization, and the call still answers
            waits = timeline.raised_waits(5.0, 1.0, None)
        assert timeline.schedule(waits, 1.0, None)[-1] is None


def test_the_shares_of_the_stages_stay_where_the_tuning_says():
    """`tuning.LIFE_STEP_MEANS_GYR` (Boss's choice) give these shares (oxygenic life, 2 Gyr of delay, star ages of 3 to 13 Gyr): 41 in 100 reach animals and
    2 in 100 a civilization."""
    import random
    chance = random.Random(11)
    animals = civilizations = 0
    runs = 6000
    for _ in range(runs):
        age = chance.uniform(3.0, 13.0)
        ages = timeline.schedule(tuple(-math.log(1.0 - chance.random()) for _ in range(timeline.STEPS)), 1.0, 2.0)
        animals += ages[3] <= age
        civilizations += ages[4] <= age
    assert 0.34 < animals / runs < 0.48
    assert 0.010 < civilizations / runs < 0.035


@pytest.mark.parametrize("pace, factor", [("fast", tuning.LIFE_PACE_FAST_ABOVE + 0.01), ("normal", 1.0),
                                          ("slow", tuning.LIFE_PACE_SLOW_BELOW - 0.01)])
def test_every_pace_reaches_complex_cells_inside_the_age_of_the_universe(pace, factor):
    """GEN.873: on main the 'slow' pace could never pass photosynthesis. A median timeline, with the oxygen coming a Gyr after the first
    photosynthesis, reaches complex cells and multicellularity inside 13.8 Gyr at every pace."""
    assert timeline.pace_of(factor) == pace
    median = (math.log(2.0),) * timeline.STEPS
    ages = timeline.schedule(median, factor, 1.0)
    assert ages[timeline.OXYGEN_STEP] < tuning.UNIVERSE_AGE_GY
    assert ages[timeline.OXYGEN_STEP + 1] < tuning.UNIVERSE_AGE_GY
