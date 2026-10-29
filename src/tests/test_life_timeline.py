# tests/test_life_timeline.py

"""GEN.750: `physics.life_timeline`, the five hard steps of life, their temperature clock and the oxygen the third one waits for."""

import math

import pytest

from planetgen import tuning
from planetgen.physics import life_timeline as timeline
from planetgen.util import draw

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


def test_a_forced_civilization_finishes_inside_the_age_and_moves_the_run_by_one_draw():
    with draw.bound("forced-a"):
        waits = timeline.forced_waits(4.0, 1.0, 1.0)
        after = draw.random()
    ages = timeline.schedule(waits, 1.0, 1.0)
    assert ages[-1] <= 4.0 + 1e-9
    with draw.bound("forced-a"):
        draw.random()
        assert draw.random() == after, "one draw of the run's stream however many tries it took"


def test_a_forced_timeline_that_cannot_fit_is_squeezed_to_fit():
    with draw.bound("forced-b"):
        waits = timeline.forced_waits(1.0, 0.25, 0.5, tries=3)
    assert timeline.schedule(waits, 0.25, 0.5)[-1] <= 1.0 + 1e-6


def test_a_forced_timeline_has_none_when_the_oxygen_never_comes_or_comes_too_late():
    with draw.bound("forced-c"):
        assert timeline.forced_waits(5.0, 1.0, None) is None
        assert timeline.forced_waits(5.0, 1.0, 5.0) is None
        assert timeline.forced_waits(5.0, 1.0, None, oxygen_at_start=True) is not None


def test_the_shares_of_the_stages_stay_where_the_tuning_says():
    """`tuning.LIFE_STEP_SCALE` is chosen for these shares (oxygenic life, 2 Gyr of delay, star ages of 3 to 13 Gyr): 24 in 100 reach animals and
    4 in 100 a civilization."""
    import random
    chance = random.Random(11)
    animals = civilizations = 0
    runs = 6000
    for _ in range(runs):
        age = chance.uniform(3.0, 13.0)
        ages = timeline.schedule(tuple(-math.log(1.0 - chance.random()) for _ in range(timeline.STEPS)), 1.0, 2.0)
        animals += ages[3] <= age
        civilizations += ages[4] <= age
    assert 0.18 < animals / runs < 0.30
    assert 0.02 < civilizations / runs < 0.07
