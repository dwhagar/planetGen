# tests/test_envelope_escape.py

"""GEN.822: the share of a starting hydrogen-helium envelope a close-in planet has lost to its star's XUV (`planetgen.physics.envelope_escape`)."""

import pytest

from planetgen.physics import envelope_escape as ee

SUN = dict(host_mass_solar=1.0, luminosity_solar=1.0, age_gy=5.0)


def _lost(core_mass, f_env0, insolation, **changes):
    return ee.lost_fraction(core_mass, f_env0, insolation, **{**SUN, **changes})


def test_a_planet_at_the_earths_distance_keeps_almost_all_of_its_envelope():
    assert _lost(5.0, 0.05, 1.0) < 0.05


def test_the_share_lost_grows_with_the_insolation_and_falls_with_the_mass_of_the_core():
    by_insolation = [_lost(6.0, 0.05, s) for s in (1.0, 10.0, 30.0, 100.0)]
    assert by_insolation == sorted(by_insolation) and by_insolation[0] < 0.1 and by_insolation[-1] > 0.9
    by_mass = [_lost(m, 0.05, 30.0) for m in (2.0, 4.0, 8.0, 16.0, 30.0)]
    assert by_mass == sorted(by_mass, reverse=True) and by_mass[0] == 1.0 and by_mass[-1] < 0.1


def test_a_planet_under_strong_irradiation_is_stripped_to_its_bare_core():
    assert _lost(3.0, 0.05, 200.0) == 1.0
    assert ee.envelope_left_earth(3.0, 0.05, 200.0, **SUN) == 0.0


def test_a_planet_younger_than_the_onset_has_lost_nothing_and_an_older_one_has_lost_at_least_as_much():
    assert _lost(4.0, 0.05, 30.0, age_gy=ee.ONSET_GY) == 0.0
    shares = [_lost(5.0, 0.05, 20.0, age_gy=age) for age in (1.5, 3.0, 5.0, 10.0)]
    assert shares == sorted(shares)


def test_a_thicker_starting_envelope_swells_the_planet_and_loses_more_mass():
    thin = 0.05 * 8.0 / 0.95 - ee.envelope_left_earth(8.0, 0.05, 30.0, **SUN)
    thick = 0.2 * 8.0 / 0.8 - ee.envelope_left_earth(8.0, 0.2, 30.0, **SUN)
    assert thick > thin > 0.0


def test_a_star_of_another_mass_gives_a_share_between_nothing_and_everything():
    for mass in (0.4, 0.8, 1.3):
        assert 0.0 <= _lost(6.0, 0.05, 20.0, host_mass_solar=mass, luminosity_solar=mass ** 3.5) <= 1.0


def test_the_starting_envelope_is_held_to_the_range_of_the_radius_law():
    assert _lost(5.0, 0.9, 20.0) == _lost(5.0, ee.MAX_ENVELOPE, 20.0)
    assert ee.envelope_left_earth(5.0, 0.0, 20.0, **SUN) == 0.0 and _lost(5.0, 0.0, 20.0) == 0.0


def test_the_equilibrium_temperature_is_255_kelvin_at_one_earth_insolation():
    assert ee.equilibrium_temperature_k(1.0) == pytest.approx(255.0)
    assert ee.equilibrium_temperature_k(16.0) == pytest.approx(510.0)


def test_a_given_equilibrium_temperature_replaces_the_one_from_the_insolation():
    assert _lost(6.0, 0.05, 30.0, teq_k=2000.0) != _lost(6.0, 0.05, 30.0)


def test_the_envelope_left_and_the_share_lost_agree():
    start = 0.05 * 6.0 / (1.0 - 0.05)
    left = ee.envelope_left_earth(6.0, 0.05, 15.0, **SUN)
    assert 0.0 < left < start
    assert _lost(6.0, 0.05, 15.0) == pytest.approx(1.0 - left / start)
