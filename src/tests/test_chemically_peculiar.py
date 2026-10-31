# tests/test_chemically_peculiar.py

"""GEN.1027: the Ap and Am stars (`planetgen/physics/chemically_peculiar.py`)."""

import pytest

from planetgen.physics import chemically_peculiar as cp


def test_the_shares_of_each_range_are_the_stand_in_shares():
    us = [(k + 0.5) / 10000.0 for k in range(10000)]
    in_both = [cp.kind_of(8500.0, u) for u in us]
    assert in_both.count(cp.AP) / len(us) == pytest.approx(cp.AP_SHARE, abs=1.0e-3)
    assert in_both.count(cp.AM) / len(us) == pytest.approx(cp.AM_SHARE, abs=1.0e-3)
    assert cp.AP not in [cp.kind_of(7200.0, u) for u in us]    # cooler than the Ap range
    assert cp.AM not in [cp.kind_of(11000.0, u) for u in us]    # hotter than the Am range
    assert {cp.kind_of(5000.0, u) for u in us} == {None} and {cp.kind_of(30000.0, u) for u in us} == {None}


def test_a_star_is_never_both():
    for temperature in (7000.0, 7500.0, 8500.0, 9500.0, 13000.0):
        kinds = [cp.kind_of(temperature, k / 1000.0) for k in range(1000)]
        assert set(kinds) <= {None, cp.AP, cp.AM}


def test_an_ap_star_rotates_in_days_inside_its_range_and_an_am_star_slower_than_100_km_s():
    periods = [cp.rotation_period_hours(cp.AP, 2.0, (k + 0.5) / 1000.0) for k in range(1000)]
    assert min(periods) >= 0.5 * 24.0 - 1.0e-9 and max(periods) <= 30.0 * 24.0 + 1.0e-9
    ordered = sorted(periods)
    assert ordered[500] / 24.0 == pytest.approx(cp.AP_PERIOD_MEDIAN_DAYS, rel=0.2)
    speed_low = 2.0 * 3.141592653589793 * 1.8 * cp.SOLAR_RADIUS_KM / (cp.rotation_period_hours(cp.AM, 1.8, 0.0) * 3600.0)
    speed_high = 2.0 * 3.141592653589793 * 1.8 * cp.SOLAR_RADIUS_KM / (cp.rotation_period_hours(cp.AM, 1.8, 1.0) * 3600.0)
    assert speed_low == pytest.approx(cp.AM_SPEED_RANGE_KMS[0]) and speed_high == pytest.approx(cp.AM_SPEED_RANGE_KMS[1])


def test_the_ap_field_is_kilogauss_scale_and_an_am_star_has_none():
    low, high = cp.AP_POLAR_FIELD_GAUSS
    assert cp.polar_field_gauss(cp.AP, 0.0) == pytest.approx(low) and cp.polar_field_gauss(cp.AP, 1.0) == pytest.approx(high)
    assert low < cp.polar_field_gauss(cp.AP, 0.5) < high
    assert cp.polar_field_gauss(cp.AM, 0.5) is None
    with pytest.raises(ValueError):
        cp.polar_field_gauss("Fm", 0.5)
    with pytest.raises(ValueError):
        cp.rotation_period_hours("Fm", 1.0, 0.5)
