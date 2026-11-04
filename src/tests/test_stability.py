# tests/test_stability.py

"""
GEN.261: the stability criteria, the Kozai-Lidov screen and the Peters merger
time (`planetgen/physics/stability.py`). The expected numbers are the unit-test
table of docs/design/systems/exotic-environments-planets-and-compact-binaries.md 2.3
and the worked examples of docs/design/systems/multi-star-systems.md 3.3 and 9.
"""

import math

import pytest

from planetgen.physics import stability


@pytest.mark.parametrize("q_out, e_out, e_in, i_deg, ma01, eq4", [
    (1e-9, 0.0, 0.0, 0.0, 2.80, 2.40),
    (1.0, 0.5, 0.0, 0.0, 9.98, 7.28),
    (0.5, 0.0, 0.0, 0.0, 3.29, 2.82),
    (0.5, 0.0, 0.0, 90.0, 2.80, 3.19),
    (0.5, 0.0, 0.0, 180.0, 2.31, 2.12),
    (0.5, 0.0, 0.6, 0.0, 3.29, 3.74),
])
def test_the_criteria_reproduce_the_unit_test_table(q_out, e_out, e_in, i_deg, ma01, eq4):
    i = math.radians(i_deg)
    assert stability.mardling_aarseth_axis_ratio(q_out, e_out, i) == pytest.approx(ma01, abs=0.01)
    assert stability.eq4_axis_ratio(q_out, e_in, e_out, i) == pytest.approx(eq4, abs=0.01)


def test_the_effective_inner_eccentricity_is_the_kozai_average_inside_the_window_only():
    assert stability.effective_inner_eccentricity(0.1, math.radians(90)) == pytest.approx(0.5)
    assert stability.effective_inner_eccentricity(0.1, math.radians(20)) == pytest.approx(0.1)
    assert stability.effective_inner_eccentricity(0.7, math.radians(90)) == pytest.approx(0.7)


def test_the_switch_picks_the_criterion_and_refuses_an_unknown_one():
    args = (1.0, 1.0, 0.5, 0.0, 0.0, 0.0)
    ma01 = stability.required_axis_ratio(*args, method="ma01")
    eq4 = stability.required_axis_ratio(*args, method="eq4")
    ek95 = stability.required_axis_ratio(*args, method="ek95")
    assert len({round(ma01, 3), round(eq4, 3), round(ek95, 3)}) == 3
    assert stability.required_axis_ratio(*args) == eq4
    with pytest.raises(ValueError):
        stability.required_axis_ratio(*args, method="nope")


def test_auto_uses_mardling_aarseth_outside_the_fitted_range_where_its_own_range_holds():
    small_pair = (1.0, 0.001, 3.0, 0.0, 0.0, 0.0)     # q_in 0.001 is outside Equation 4; q_out 3 is inside Mardling and Aarseth's range
    assert stability.criterion(1.0, 0.001, 3.0) == "ma01"
    assert stability.required_axis_ratio(*small_pair) == pytest.approx(stability.required_axis_ratio(*small_pair, method="ma01"))
    assert not stability.eq4_in_range(0.001, 0.5, 0.1)
    assert stability.eq4_in_range(0.5, 0.5, 0.1)


def test_mardling_aarseth_is_not_used_past_the_outer_mass_ratio_it_is_given_for():
    """GEN.1028: past q_out = 5 and Equation 4's range, auto takes Eggleton and Kiseleva and a named ma01 is refused."""
    extreme = (1.0, 0.001, 100.0, 0.0, 0.0, 0.0)
    assert stability.criterion(1.0, 0.001, 100.0) == "ek95"
    assert stability.required_axis_ratio(*extreme) == pytest.approx(stability.required_axis_ratio(*extreme, method="ek95"))
    with pytest.raises(ValueError):
        stability.required_axis_ratio(*extreme, method="ma01")
    assert stability.ma01_in_range(5.0) and not stability.ma01_in_range(5.01)


def test_the_moon_of_a_planet_about_a_star_is_not_failed_by_a_criterion_outside_its_range():
    """The Earth and the Moon (q_in 0.0123, q_out 3e5): Equation 4 stops at q_out 100, and Mardling and Aarseth at 5 would have asked
    for an axis ratio of 430, more than the real 390; Eggleton and Kiseleva's 7 holds."""
    earth, moon, sun = 3.0e-6, 3.7e-8, 1.0
    assert stability.criterion(earth, moon, sun) == "ek95"
    assert stability.is_triple_stable(earth, moon, sun, 0.00257, 1.0, 0.055, 0.017, 0.09)


def test_a_triple_is_stable_just_outside_its_criterion_and_not_just_inside():
    ratio = stability.required_axis_ratio(1.0, 1.0, 1.0, 0.1, 0.2, 0.3, "ma01")
    assert stability.is_triple_stable(1, 1, 1, 1.0, ratio * 1.01, 0.1, 0.2, 0.3, "ma01")
    assert not stability.is_triple_stable(1, 1, 1, 1.0, ratio * 0.99, 0.1, 0.2, 0.3, "ma01")


def test_the_largest_inner_axis_is_just_stable():
    a_in = stability.max_inner_axis(1.0, 0.8, 0.6, 100.0, 0.0, 0.3, 0.4)
    assert stability.is_triple_stable(1.0, 0.8, 0.6, a_in * 0.99, 100.0, 0.0, 0.3, 0.4)
    assert not stability.is_triple_stable(1.0, 0.8, 0.6, a_in * 1.05, 100.0, 0.0, 0.3, 0.4)


def test_kepler_helpers_round_trip():
    assert stability.orbital_period_yr(1.0, 1.0) == pytest.approx(1.0)
    assert stability.axis_au(stability.orbital_period_yr(37.0, 2.5), 2.5) == pytest.approx(37.0)
    assert stability.mutual_inclination(0.3, 1.0, 0.3, 1.0) == pytest.approx(0.0)
    assert stability.mutual_inclination(0.0, 0.0, math.pi / 2, 0.0) == pytest.approx(math.pi / 2)


def test_v404_cygni_is_stable_by_a_wide_margin_and_its_kozai_cycle_is_quenched():
    period_in = 6.5 / 365.25
    a_in = stability.axis_au(period_in, 9.7)
    a_out = stability.axis_au(7.0e4, 10.6)
    assert a_in == pytest.approx(0.145, abs=0.002)
    assert a_out == pytest.approx(3731, rel=0.01)
    assert a_out / a_in == pytest.approx(25_700, rel=0.01)
    for e_out in (0.0, 0.5, 0.9):
        assert stability.is_triple_stable(9.0, 0.7, 0.9, a_in, a_out, 0.0, e_out, 0.0)
    assert stability.gr_precession_period_yr(a_in, 0.0, 9.7, period_in) == pytest.approx(9000, rel=0.05)
    assert stability.kozai_timescale_yr(period_in, 7.0e4, 10.6, 0.9, 0.0) == pytest.approx(6.9e11, rel=0.05)
    assert stability.kozai_verdict(9.0, 0.7, 0.9, a_in, a_out, 0.0, 0.0, math.radians(60), 1.0e10, 0.01) == "safe"


def test_the_kozai_screen_flags_a_cycle_that_reaches_the_stars():
    kwargs = dict(m1=1.0, m2=1.0, m3=1.0, a_in=1.0, a_out=30.0, e_in=0.0, e_out=0.3, i=math.radians(70), age_yr=1e10)
    assert stability.kozai_verdict(**kwargs, extent_sum_au=0.05) == "unsafe"
    assert stability.kozai_verdict(**kwargs, extent_sum_au=0.0001) == "kl_active"
    assert stability.kozai_verdict(**{**kwargs, "i": math.radians(5)}, extent_sum_au=0.05) == "safe"
    assert stability.kozai_verdict(**{**kwargs, "age_yr": 1e3}, extent_sum_au=0.05) == "safe"


def test_the_quadrupole_maximum_eccentricity_is_the_textbook_value():
    assert stability.kozai_max_eccentricity(0.0, math.radians(90)) == pytest.approx(math.sqrt(1 - 5 / 3 * 0.0))
    assert stability.kozai_max_eccentricity(0.0, math.radians(60)) == pytest.approx(math.sqrt(1 - 5 / 3 * 0.25))
    assert stability.kozai_max_eccentricity(0.2, math.radians(10)) == pytest.approx(0.2)


def test_the_double_pulsar_merges_in_85_million_years_and_fails_an_old_system():
    period_yr = 2.45 / 24.0 / 365.25
    a = stability.axis_au(period_yr, 1.34 + 1.25)
    assert a == pytest.approx(0.0059, abs=0.0003)
    assert stability.peters_merger_time_yr(1.34, 1.25, a, 0.088) == pytest.approx(85.0e6, rel=0.05)
    assert not stability.survives_gravitational_waves(1.34, 1.25, a, 0.088, 10.0e9)
    assert stability.survives_gravitational_waves(1.34, 1.25, a, 0.088, 10.0e6)


@pytest.mark.parametrize("m1, m2, a_au, e, expected_gyr", [
    (1.4, 1.4, 0.022, 0.0, 13.8),
    (1.4, 8.0, 0.046, 0.0, 13.8),
    (10.0, 10.0, 0.096, 0.0, 13.8),
    (1.4, 1.4, 0.094, 0.9, 13.8),
])
def test_the_hubble_time_axes_of_the_note_merge_in_a_hubble_time(m1, m2, a_au, e, expected_gyr):
    assert stability.peters_merger_time_yr(m1, m2, a_au, e) / 1.0e9 == pytest.approx(expected_gyr, rel=0.05)
