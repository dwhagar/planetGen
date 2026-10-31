"""
The branches of a climate with an ice-albedo feedback (GEN.1046, `physics/climate_branches.py`): the albedo law of the ice, the fixed points of a
response, their stability, and the window of flux in which two stable climates exist. The response here is a grey toy (a held greenhouse factor
over the effective temperature), a test target and not an input of the code.
"""

import numpy as np
import pytest

from planetgen.physics import climate_branches as cb, constants

SIGMA = constants.STEFAN_BOLTZMANN_CONSTANT
EARTH_FLUX = 1361.0


def _response(flux, greenhouse=0.63):
    return lambda albedo: (flux * (1.0 - albedo) / 4.0 / SIGMA) ** 0.25 * (1.0 + greenhouse) ** 0.25


def test_the_law_is_ice_free_above_covered_below_and_a_straight_line_between():
    law = cb.IceAlbedoLaw()
    assert law.albedo(300.0) == law.ice_free and law.albedo(200.0) == law.ice_covered
    middle = 0.5 * (law.free_above_k + law.covered_below_k)
    assert law.ice_fraction(middle) == pytest.approx(0.5)
    assert law.albedo(middle) == pytest.approx(0.5 * (law.ice_free + law.ice_covered))
    assert [law.albedo(t) for t in (270.0, 272.0, 274.0)] == sorted((law.albedo(t) for t in (270.0, 272.0, 274.0)), reverse=True)


@pytest.mark.parametrize("kwargs", [{"ice_free": 0.7, "ice_covered": 0.5}, {"ice_covered": 1.0}, {"free_above_k": 250.0, "covered_below_k": 260.0}])
def test_a_law_that_is_not_a_ramp_from_dark_to_bright_and_warm_to_cold_is_refused(kwargs):
    with pytest.raises(ValueError):
        cb.IceAlbedoLaw(**kwargs)


def test_an_earth_twin_holds_a_snowball_a_warm_world_and_an_unstable_one_between():
    branches = cb.fixed_points(_response(EARTH_FLUX))
    assert [branch.stable for branch in branches] == [True, False, True]
    snowball, between, warm = branches
    assert snowball.temperature_k < cb.ICE_COVERED_BELOW_K and snowball.ice_fraction == 1.0 and snowball.albedo == cb.ICE_COVERED_ALBEDO
    assert warm.temperature_k > cb.ICE_FREE_ABOVE_K and warm.ice_fraction == 0.0 and warm.albedo == cb.ICE_FREE_ALBEDO
    assert snowball.temperature_k < between.temperature_k < warm.temperature_k
    assert 0.0 < between.ice_fraction < 1.0


def test_every_branch_is_a_fixed_point_of_its_own_loop():
    law = cb.IceAlbedoLaw()
    response = _response(EARTH_FLUX)
    for branch in cb.fixed_points(response, law):
        assert response(law.albedo(branch.temperature_k)) == pytest.approx(branch.temperature_k, abs=1.0e-6)
        assert branch.albedo == pytest.approx(law.albedo(branch.temperature_k))


def test_a_cold_world_has_only_the_snowball_and_a_hot_one_only_the_warm_branch():
    cold = cb.fixed_points(_response(900.0))
    hot = cb.fixed_points(_response(2000.0))
    assert len(cold) == 1 and cold[0].ice_fraction == 1.0 and cold[0].stable
    assert len(hot) == 1 and hot[0].ice_fraction == 0.0 and hot[0].stable


def test_a_world_whose_temperature_does_not_depend_on_the_albedo_has_one_branch():
    branches = cb.fixed_points(lambda albedo: 280.0)
    assert [(round(b.temperature_k, 6), b.stable) for b in branches] == [(280.0, True)]


def test_the_stable_branches_are_the_coldest_and_the_warmest():
    stable = cb.stable_branches(_response(EARTH_FLUX))
    assert len(stable) == 2 and stable[0].temperature_k < stable[1].temperature_k
    assert cb.is_bistable(_response(EARTH_FLUX)) and not cb.is_bistable(_response(900.0))


def test_the_bistable_window_brackets_the_flux_of_an_earth_twin_and_is_none_outside_every_flux_tried():
    fluxes = np.arange(900.0, 2200.0, 20.0)
    low, high = cb.bistable_window(_response_at, fluxes)
    assert low < EARTH_FLUX < high
    assert not cb.is_bistable(_response_at(low - 40.0)) and not cb.is_bistable(_response_at(high + 40.0))
    assert cb.bistable_window(_response_at, np.array([800.0, 850.0])) is None


def _response_at(flux):
    return _response(flux)


def test_a_snowball_response_that_falls_with_the_albedo_more_gently_has_a_narrower_window():
    steep = cb.bistable_window(_response_at, np.arange(900.0, 2200.0, 10.0))
    gentle_law = cb.IceAlbedoLaw(ice_covered=0.45)
    gentle = cb.bistable_window(_response_at, np.arange(900.0, 2200.0, 10.0), gentle_law)
    assert gentle is None or gentle[1] - gentle[0] < steep[1] - steep[0]
