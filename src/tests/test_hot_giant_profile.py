# tests/test_hot_giant_profile.py

"""GEN.1577: the Guillot 2010 semi-grey profile of a hot giant's atmosphere (`physics/hot_giant_profile.py`), held to the worked example of the research (part 03 section 5.1)."""

import math

import pytest

from planetgen.physics import constants, hot_giant_profile as hgp

SIGMA = constants.STEFAN_BOLTZMANN_CONSTANT
HD209458 = dict(equilibrium_k=1469.0, internal_flux_w_m2=SIGMA * 100.0 ** 4, gravity_m_s2=10.0, gamma=0.4)
"""dict: Guillot's HD 209458 b example: kappa_th 1e-2 cm2/g (the module's value), kappa_v 4e-3 (gamma 0.4), T_int 100 K, T_eq 1469 K, g 1000 cm/s2, f = 1/4."""


@pytest.mark.parametrize("pressure_bar, expected_k", [(0.0667, 1473.0), (1.0, 1647.0), (10.0, 1648.0)])
def test_the_worked_example_of_the_research(pressure_bar, expected_k):
    assert hgp.temperature_k(pressure_bar, **HD209458) == pytest.approx(expected_k, abs=1.0)


def test_the_top_is_the_skin_and_the_deep_end_the_isotherm():
    assert hgp.temperature_k(1.0e-9, **HD209458) == pytest.approx(hgp.skin_k(1469.0, HD209458["internal_flux_w_m2"], gamma=0.4), rel=1e-6)
    deep = dict(HD209458, internal_flux_w_m2=0.0)
    isotherm = (0.75 * 1469.0 ** 4 * 4 * 0.25 * (2.0 / 3.0 + 1.0 / (0.4 * math.sqrt(3.0)))) ** 0.25
    assert hgp.temperature_k(1.0e3, **deep) == pytest.approx(isotherm, rel=1e-6)


def test_the_profile_rises_with_depth_for_a_modest_gamma_and_inverts_for_a_large_one():
    pressures = [10.0 ** (exponent / 2.0) for exponent in range(-12, 7)]            # 1e-6 to 1e3 bar
    modest = [hgp.temperature_k(p, **dict(HD209458, gamma=0.4)) for p in pressures]
    assert modest == sorted(modest)
    inverted = [hgp.temperature_k(p, **dict(HD209458, gamma=8.0, internal_flux_w_m2=0.0)) for p in pressures]
    assert inverted[0] > inverted[-1]


def test_the_internal_heat_warms_only_the_deep_end():
    cool = dict(HD209458, internal_flux_w_m2=0.0)
    assert hgp.temperature_k(1.0e-5, **HD209458) == pytest.approx(hgp.temperature_k(1.0e-5, **cool), rel=0.01)
    assert hgp.temperature_k(1.0e3, **HD209458) > 1.01 * hgp.temperature_k(1.0e3, **cool)


def test_a_hotter_planet_has_a_hotter_profile_and_the_dayside_is_hotter_than_the_globe():
    assert hgp.temperature_k(1.0, **dict(HD209458, equilibrium_k=2000.0)) > hgp.temperature_k(1.0, **HD209458)
    dayside = hgp.temperature_k(1.0, **HD209458, redistribution=0.5)
    assert dayside == pytest.approx(hgp.temperature_k(1.0, **HD209458) * 2.0 ** 0.25, rel=0.02)


def test_gamma_follows_the_irradiation_when_not_given_and_the_function_matches_temperature_k():
    t_irr = hgp.irradiation_temperature_k(1469.0)
    assert t_irr == pytest.approx(1469.0 * math.sqrt(2.0))
    assert hgp.opacity_ratio(2000.0) == pytest.approx(0.6)
    no_gamma = dict(HD209458)
    del no_gamma["gamma"]
    profile = hgp.profile_function(**no_gamma)
    assert profile(1.0) == hgp.temperature_k(1.0, **no_gamma)


def test_bad_inputs_are_refused():
    for bad in (dict(pressure_bar=0.0), dict(pressure_bar=-1.0)):
        with pytest.raises(ValueError):
            hgp.temperature_k(bad["pressure_bar"], **HD209458)
    with pytest.raises(ValueError):
        hgp.temperature_k(1.0, **dict(HD209458, gravity_m_s2=0.0))
    with pytest.raises(ValueError):
        hgp.temperature_k(1.0, **dict(HD209458, gamma=-1.0))
    with pytest.raises(ValueError):
        hgp.internal_temperature_k(-1.0)
