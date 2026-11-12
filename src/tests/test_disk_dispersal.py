# tests/test_disk_dispersal.py

import math

from planetgen.physics import disk_dispersal as dd

SAMPLES = 4000


def _times(companion_au=None, mass=1.0):
    return [dd.disk_dispersal_time_yr(mass, ("star", i), companion_au) for i in range(SAMPLES)]


def test_curve_inverse_round_trips():
    for age in (0.0, 0.4e6, 1.0e6, 2.0e6, 3.0e6, 4.5e6, 5.5e6, 8.0e6):
        assert math.isclose(dd.survival_to_time_yr(dd.disk_share(age)), age, rel_tol=1e-9, abs_tol=1.0)


def test_curve_passes_the_survey_points():
    assert dd.disk_share(0.0) == 1.0
    assert math.isclose(dd.disk_share(1.0e6), 0.80)
    assert math.isclose(dd.disk_share(3.0e6), 0.50)
    assert math.isclose(dd.disk_share(5.5e6), 0.10)
    assert dd.disk_share(10.0e6) == 0.0 and dd.disk_share(50.0e6) == 0.0


def test_drawn_times_follow_the_curve_and_are_gone_by_ten_myr():
    times = _times()
    assert max(times) <= 10.0e6
    for age in (1.0e6, 3.0e6, 5.5e6):
        held = sum(1 for t in times if t > age) / SAMPLES
        assert abs(held - dd.disk_share(age)) < 0.03


def test_same_star_same_time_and_age_never_enters_the_draw():
    assert dd.disk_dispersal_time_yr(0.5, "a") == dd.disk_dispersal_time_yr(0.5, "a")
    assert dd.disk_dispersal_time_yr(0.5, "a") != dd.disk_dispersal_time_yr(0.5, "b")
    assert dd.disk_present(0.1e6, 0.5, "a") and not dd.disk_present(11.0e6, 0.5, "a")


def test_a_close_companion_cuts_the_gas_early_and_a_far_one_does_not():
    alone, close, far = _times(), _times(10.0), _times(500.0)
    assert far == alone
    assert all(c <= a for c, a in zip(close, alone))
    gone_in_a_million = sum(1 for t in close if t < dd.COMPANION_EARLY_TIME_YR) / SAMPLES
    assert gone_in_a_million > 0.6
    assert sum(1 for t in alone if t < dd.COMPANION_EARLY_TIME_YR) / SAMPLES < 0.25
