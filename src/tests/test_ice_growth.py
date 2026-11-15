# tests/test_ice_growth.py

"""GEN.1728: how thick a world's ice and the snow on it are (`physics/ice_growth.py`): the degree days of a sine year, Stefan's law with snow, the perennial balance, the snow and the whole lid."""

import math

import pytest

from planetgen.physics import continental_weathering, ice_growth as ig, ice_light


def test_the_degree_days_of_a_sine_year_match_a_direct_sum():
    for mean, amplitude in ((258.0, 18.0), (270.0, 10.0), (273.15, 5.0), (250.0, 40.0)):
        steps = 20000
        direct = sum(max(ig.FREEZING_K - (mean + amplitude * math.cos(2.0 * math.pi * i / steps)), 0.0) for i in range(steps)) * ig.DAYS_PER_YEAR / steps
        assert ig.freezing_degree_days(mean, amplitude) == pytest.approx(direct, rel=1e-4)


def test_a_year_that_never_freezes_has_no_degree_days_and_one_that_never_thaws_has_them_all():
    assert ig.freezing_degree_days(300.0, 10.0) == 0.0 and ig.frozen_share(300.0, 10.0) == 0.0
    assert ig.freezing_degree_days(240.0, 10.0) == pytest.approx(33.15 * ig.DAYS_PER_YEAR) and ig.frozen_share(240.0, 10.0) == 1.0
    assert ig.frozen_share(ig.FREEZING_K, 10.0) == pytest.approx(0.5)
    assert ig.freezing_degree_days(250.0, 0.0) == pytest.approx(23.15 * ig.DAYS_PER_YEAR)


def test_more_cold_grows_more_ice_and_snow_slows_it():
    values = [ig.winter_growth_m(u) for u in (0.0, 100.0, 500.0, 2000.0, 8000.0)]
    assert values[0] == pytest.approx(0.0) and values == sorted(values)
    assert ig.winter_growth_m(2000.0, 0.3) < ig.winter_growth_m(2000.0, 0.1) < ig.winter_growth_m(2000.0)


def test_two_thousand_degree_days_make_a_first_year_arctic_lid():
    # first-year Arctic sea ice is 1.3 to 2 m thick after 2000 to 3000 degree days; the law with the air's resistance gives 1.4 m, and 15 cm of snow takes off about 0.8 m
    assert 1.2 < ig.winter_growth_m(2000.0) < 1.8
    assert 0.5 < ig.winter_growth_m(2000.0, 0.15) < 1.0


def test_stefans_law_without_resistance_is_a_square_root():
    k = ig.ice_conductivity_w_m_k(ig.FREEZING_K - 10.0)
    h = ig.winter_growth_m(1000.0, 0.0, air_coupling_w_m2_k=1e12)
    assert h == pytest.approx(math.sqrt(2.0 * k * 1000.0 * 86400.0 / (ig.ICE_DENSITY_KG_M3 * ig.LATENT_HEAT_J_KG)), rel=1e-3)
    assert ig.winter_growth_m(4000.0, 0.0, air_coupling_w_m2_k=1e12) == pytest.approx(2.0 * h, rel=1e-3)


def test_the_snow_conductivity_is_the_measured_one():
    assert ig.snow_conductivity_w_m_k(300.0) == pytest.approx(0.23, abs=0.01)
    assert ig.snow_conductivity_w_m_k(250.0) < ig.snow_conductivity_w_m_k(320.0)


def test_the_standard_arctic_case_gives_the_multiyear_thickness_of_the_model_that_made_it():
    # [R5271]: 2.88 m with 2 W per m2 from the water, 40 cm of ablation and 34 cm of snow ([R5270]) under a mean surface of about 256 K
    assert ig.perennial_thickness_m(256.0, 2.0, 0.34, 0.4) == pytest.approx(2.88, abs=0.5)


def test_the_perennial_lid_thickens_with_cold_and_thins_with_snow_water_heat_and_wear():
    base = ig.perennial_thickness_m(250.0)
    assert ig.perennial_thickness_m(240.0) > base
    assert ig.perennial_thickness_m(250.0, snow_m=0.3) < base
    assert ig.perennial_thickness_m(250.0, ocean_flux_w_m2=4.0) < base
    assert ig.perennial_thickness_m(250.0, ablation_m_yr=0.0) > base
    assert ig.perennial_thickness_m(280.0) == 0.0


def test_the_conduction_limit_is_the_hydrosphere_shell():
    from planetgen.physics import hydrosphere
    assert ig.conduction_limit_m(240.0, 0.09) == pytest.approx(hydrosphere.ice_shell_thickness_km(0.09, 240.0) * 1000.0)


def test_the_precipitation_share_is_the_weathering_laws_earth_identity():
    assert ig.PRECIPITATION_ENERGY_SHARE == pytest.approx(0.32, abs=0.02)
    earth = ig.snowfall_m_yr(continental_weathering.EARTH_ABSORBED_FLUX_W_M2, continental_weathering.REFERENCE_K)
    assert earth == pytest.approx(continental_weathering.PRECIPITATION_REF_M_YR, rel=1e-6)


def test_the_snow_follows_the_starlight_and_the_length_of_the_freezing_season():
    assert ig.snowfall_m_yr(0.0, 255.0) == 0.0
    assert ig.snowfall_m_yr(200.0, 255.0) > ig.snowfall_m_yr(100.0, 255.0)
    assert ig.snow_depth_m(100.0, 300.0, 10.0) == 0.0
    assert ig.snow_depth_m(100.0, 265.0, 20.0) < ig.snow_depth_m(100.0, 250.0, 20.0)
    assert ig.snowfall_m_yr(100.0, 20.0) == 0.0


def test_the_arctic_snow_is_within_a_factor_of_two_of_the_measured_depth():
    # [R5270]: 34 cm on multiyear Arctic ice; about 75 W per m2 of starlight absorbed at 255 K. The law has no dry cold air, so it errs high (IG6 says so)
    depth = ig.snow_depth_m(75.0, 255.0, 17.0)
    assert 0.34 / 2.0 < depth < 0.34 * 2.0


def test_a_year_that_never_freezes_has_no_lid():
    assert ig.lid_of(290.0, 10.0, 200.0) == ig.Lid(0.0, 0.0, False)


def test_a_cool_world_wears_a_seasonal_lid_and_a_cold_one_a_perennial_one():
    seasonal = ig.lid_of(268.0, 20.0, 200.0)
    assert not seasonal.perennial and 0.1 < seasonal.ice_m < 3.0 and seasonal.snow_m > 0.0
    perennial = ig.lid_of(240.0, 8.0, 100.0)
    assert perennial.perennial and perennial.ice_m > seasonal.ice_m


def test_the_lid_feeds_the_light_law():
    lid = ig.lid_of(262.0, 18.0, 150.0)
    assert 0.0 < ice_light.transmission(lid.ice_m, lid.snow_m) < 0.05


def test_a_colder_world_has_a_darker_lid():
    light = [ice_light.transmission(*(lambda lid: (lid.ice_m, lid.snow_m))(ig.lid_of(t, 10.0, 100.0))) for t in (262.0, 250.0, 235.0)]
    assert light == sorted(light, reverse=True)


@pytest.mark.parametrize("call", [
    lambda: ig.freezing_degree_days(250.0, -1.0),
    lambda: ig.winter_growth_m(-5.0),
    lambda: ig.winter_growth_m(10.0, -0.1),
    lambda: ig.perennial_thickness_m(250.0, 0.0),
    lambda: ig.conduction_limit_m(250.0, -1.0),
    lambda: ig.snowfall_m_yr(-1.0, 250.0),
    lambda: ig.snow_conductivity_w_m_k(0.0),
])
def test_bad_inputs_are_refused(call):
    with pytest.raises(ig.IceGrowthError):
        call()
