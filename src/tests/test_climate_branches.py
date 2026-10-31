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


def test_the_law_is_ice_free_above_covered_below_and_a_straight_line_between_when_the_exponent_is_one():
    law = cb.IceAlbedoLaw(exponent=1.0)
    assert law.albedo(law.free_above_k + 5.0) == law.ice_free and law.albedo(law.covered_below_k - 5.0) == law.ice_covered
    middle = 0.5 * (law.free_above_k + law.covered_below_k)
    assert law.ice_fraction(middle) == pytest.approx(0.5)
    assert law.albedo(middle) == pytest.approx(0.5 * (law.ice_free + law.ice_covered))
    assert [law.albedo(t) for t in (270.0, 272.0, 274.0)] == sorted((law.albedo(t) for t in (270.0, 272.0, 274.0)), reverse=True)


def test_a_cubic_ice_cover_grows_slowly_at_first_and_fast_when_the_ice_line_reaches_the_tropics():
    law = cb.IceAlbedoLaw()
    quarter, half, three_quarters = (law.free_above_k - q * (law.free_above_k - law.covered_below_k) for q in (0.25, 0.5, 0.75))
    assert law.ice_fraction(quarter) == pytest.approx(0.25 ** 3) and law.ice_fraction(half) == pytest.approx(0.125) and law.ice_fraction(three_quarters) == pytest.approx(0.75 ** 3)


def test_the_surface_albedo_feedback_at_the_earths_temperature_is_inside_the_range_of_the_assessment():
    """CC16: IPCC AR6 (2021) table 7.10 gives 0.35 W/m2 per K, very likely 0.10 to 0.60; the Earth's incoming flux 1361 W/m2 over four."""
    law = cb.IceAlbedoLaw()
    feedback = -(law.albedo(288.5) - law.albedo(287.5)) * EARTH_FLUX / 4.0
    assert 0.10 < feedback < 0.60 and feedback == pytest.approx(0.35, abs=0.1)


@pytest.mark.parametrize("kwargs", [{"ice_free": 0.7, "ice_covered": 0.5}, {"ice_covered": 1.0}, {"free_above_k": 250.0, "covered_below_k": 260.0}, {"exponent": 0.0}])
def test_a_law_that_is_not_a_ramp_from_dark_to_bright_and_warm_to_cold_is_refused(kwargs):
    with pytest.raises(ValueError):
        cb.IceAlbedoLaw(**kwargs)


FAINT_SUN_FLUX = 1150.0
"""float: A flux 15 percent under the Earth's, inside the bistable window of the toy (a faint young Sun, as in the snowball Earth of Hoffman 1998)."""


def test_an_earth_twin_under_a_fainter_sun_holds_a_snowball_a_warm_world_and_an_unstable_one_between():
    branches = cb.fixed_points(_response(FAINT_SUN_FLUX))
    assert [branch.stable for branch in branches] == [True, False, True]
    snowball, between, warm = branches
    assert snowball.temperature_k < cb.ICE_COVERED_BELOW_K and snowball.ice_fraction == 1.0 and snowball.albedo == cb.ICE_COVERED_ALBEDO
    assert warm.temperature_k > snowball.temperature_k + 20.0 and warm.ice_fraction < 0.25 and warm.albedo < cb.ICE_FREE_ALBEDO + 0.08
    assert snowball.temperature_k < between.temperature_k < warm.temperature_k
    assert 0.0 < between.ice_fraction < 1.0


def test_every_branch_is_a_fixed_point_of_its_own_loop():
    law = cb.IceAlbedoLaw()
    response = _response(FAINT_SUN_FLUX)
    for branch in cb.fixed_points(response, law):
        assert response(law.albedo(branch.temperature_k)) == pytest.approx(branch.temperature_k, abs=1.0e-6)
        assert branch.albedo == pytest.approx(law.albedo(branch.temperature_k))


def test_a_cold_world_has_only_the_snowball_and_a_hot_one_only_the_warm_branch():
    cold = cb.fixed_points(_response(900.0))
    hot = cb.fixed_points(_response(2000.0))
    assert len(cold) == 1 and cold[0].ice_fraction == 1.0 and cold[0].stable
    assert len(hot) == 1 and hot[0].ice_fraction == 0.0 and hot[0].stable


def test_the_earth_at_its_own_flux_has_only_the_warm_branch():
    branches = cb.fixed_points(_response(EARTH_FLUX))
    assert len(branches) == 1 and branches[0].stable and branches[0].temperature_k > 280.0


def test_a_world_whose_temperature_does_not_depend_on_the_albedo_has_one_branch():
    branches = cb.fixed_points(lambda albedo: 280.0)
    assert [(round(b.temperature_k, 6), b.stable) for b in branches] == [(280.0, True)]


def test_the_stable_branches_are_the_coldest_and_the_warmest():
    stable = cb.stable_branches(_response(FAINT_SUN_FLUX))
    assert len(stable) == 2 and stable[0].temperature_k < stable[1].temperature_k
    assert cb.is_bistable(_response(FAINT_SUN_FLUX)) and not cb.is_bistable(_response(900.0)) and not cb.is_bistable(_response(EARTH_FLUX))


def test_the_bistable_window_of_an_earth_twin_lies_under_the_present_flux_and_is_a_tenth_wide():
    """W20 asks for a window 10 to 20 percent of the flux wide for an Earth twin; the toy holds the greenhouse of the Earth fixed, which the vapour and cloud of
    the full response will move."""
    fluxes = np.arange(900.0, 2200.0, 10.0)
    low, high = cb.bistable_window(_response_at, fluxes)
    assert 0.75 * EARTH_FLUX < low < high < EARTH_FLUX and low < FAINT_SUN_FLUX < high
    assert 0.07 < (high - low) / EARTH_FLUX < 0.20
    assert not cb.is_bistable(_response_at(low - 40.0)) and not cb.is_bistable(_response_at(high + 40.0))
    assert cb.bistable_window(_response_at, np.array([800.0, 850.0])) is None


def _response_at(flux):
    return _response(flux)


def test_a_snowball_that_reflects_less_has_a_narrower_window():
    fluxes = np.arange(900.0, 2200.0, 10.0)
    bright = cb.bistable_window(_response_at, fluxes)
    dull = cb.bistable_window(_response_at, fluxes, cb.IceAlbedoLaw(ice_covered=0.50))
    assert dull is None or dull[1] - dull[0] < bright[1] - bright[0]
