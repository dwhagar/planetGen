"""
The greenhouse variable and the flux behind it (GEN.443, `planetgen.physics.greenhouse`). The surface temperature itself is `greenhouse_table`'s (GEN.850,
`test_greenhouse_table.py`); Earth and the Sun are test targets, never inputs of the code.
"""

import pytest

from planetgen.physics import greenhouse as gh


def test_absorbed_flux_of_the_sun_at_one_au():
    assert gh.absorbed_flux_w_m2(1361.17, 0.30) == pytest.approx(238.2, rel=0.001)


def test_a_tidally_locked_world_with_a_hot_dayside_absorbs_the_larger_share():
    assert gh.absorbed_flux_w_m2(1361.17, 0.30, redistribution=0.5) == pytest.approx(2.0 * gh.absorbed_flux_w_m2(1361.17, 0.30))


def test_the_effective_temperature_is_the_one_that_radiates_the_flux():
    assert gh.effective_temperature_k(238.2) == pytest.approx(254.6, abs=0.2)


def test_the_factor_is_the_fourth_power_of_the_temperature_ratio_less_one():
    assert gh.greenhouse_factor(288.0, 255.0) == pytest.approx((288.0 / 255.0) ** 4 - 1.0)
    assert gh.greenhouse_factor(255.0, 255.0) == 0.0
