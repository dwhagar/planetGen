# tests/test_saturation_time.py

"""GEN.493: the X-ray saturation time by stellar mass and rotation percentile (`physics/saturation_time.py`), read off [R2183: Fig. 10]."""

import pytest

from planetgen.physics import saturation_time as sat


def test_at_one_solar_mass_it_agrees_with_tu_2015():
    """Tu et al. 2015: 5.7, 23 and 226 Myr at the 10th, 50th and 90th percentiles; the figure's 5th, 50th and 95th are 5.2, 22.8 and 297."""
    assert sat.saturation_time_myr(1.0, 50.0) == pytest.approx(23.0, rel=0.05)
    assert sat.saturation_time_myr(1.0, 10.0) == pytest.approx(5.7, rel=0.3)
    assert sat.saturation_time_myr(1.0, 90.0) == pytest.approx(226.0, rel=0.4)


def test_a_lighter_star_stays_saturated_longer():
    for percentile in (5.0, 50.0, 95.0):
        values = [sat.saturation_time_myr(m, percentile) for m in (0.15, 0.3, 0.5, 0.8, 1.1)]
        assert values == sorted(values, reverse=True)


def test_a_faster_rotator_stays_saturated_longer_and_the_percentile_is_monotone():
    for mass in (0.15, 0.4, 0.9, 1.0):
        values = [sat.saturation_time_myr(mass, p) for p in (2.0, 5.0, 25.0, 50.0, 75.0, 95.0, 98.0)]
        assert values == sorted(values)
    assert sat.saturation_time_myr(0.5, 50.0) == pytest.approx(491.0)


def test_the_published_curves_are_reproduced_at_their_anchors():
    for mass, row in sat.ANCHORS_MYR.items():
        for percentile, value in zip(sat.PERCENTILES, row):
            assert sat.saturation_time_myr(mass, percentile) == pytest.approx(value, rel=1e-9)


def test_the_mid_to_late_m_dwarfs_stay_saturated_for_gigayears_as_medina_2022_finds():
    """[R3017]: the saturated to unsaturated flaring transition of 0.1 to 0.3 solar mass stars is at 2.4 +- 0.3 Gyr on average (kinematic ages)."""
    masses = [0.1 + 0.01 * i for i in range(21)]
    weights = [m ** -1.3 for m in masses]       # a Kroupa-like mass function below half a solar mass
    mean = sum(w * sat.saturation_time_myr(m, 50.0) for m, w in zip(masses, weights)) / sum(weights)
    assert 1500 < mean < 3300


def test_outside_the_models_range_it_takes_the_nearest_anchor_and_a_bad_percentile_raises():
    assert sat.saturation_time_myr(0.05, 50.0) == sat.saturation_time_myr(0.1, 50.0)
    assert sat.saturation_time_myr(2.0, 50.0) == sat.saturation_time_myr(1.2, 50.0)
    with pytest.raises(ValueError):
        sat.saturation_time_myr(1.0, 120.0)
