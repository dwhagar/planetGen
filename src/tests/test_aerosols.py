"""
Haze, dust and cloud (GEN.414, `planetgen.physics.aerosols`): the condensation level and the cloud limit, the methane haze gate and
depth, the dust draw, the layers handed to the radiation kernel, and the fixed number of draws.
"""

import math
import random

import numpy as np
import pytest

from planetgen.datasets import fetch
from planetgen.physics import aerosols, air_column, atmospheric_radiation as ar, constants, stellar_spectrum as ss


def _air(partials, gravity):
    return air_column.from_partials(partials, gravity)


EARTH = _air({"N2": 78084.0, "O2": 20946.0, "Ar": 934.0, "CO2": 40.0, "H2O": 1000.0}, 9.80665)
MARS = _air({"CO2": 570.0, "N2": 15.0, "Ar": 18.0, "H2O": 0.2}, 3.71)
VENUS = _air({"CO2": 8.9e6, "N2": 3.4e5, "H2O": 180.0, "SO2": 1e3}, 8.87)
TITAN = _air({"N2": 133000.0, "CH4": 7000.0}, 1.35)
DRY_EARTH = _air({"N2": 78084.0, "O2": 20946.0, "Ar": 934.0, "CO2": 40.0}, 9.80665)
HOT_JUPITER_LIKE = _air({"H2": 86000.0, "He": 14000.0}, 24.8)


class Counting:
    def __init__(self, seed):
        self.inner = random.Random(seed)
        self.calls = 0

    def random(self):
        self.calls += 1
        return self.inner.random()


def test_the_profile_is_a_dry_adiabat_under_an_isothermal_skin():
    pressure, temperature = aerosols.profile(EARTH, 288.0)
    assert temperature[0] == 288.0 and pressure[0] == EARTH.pressure_pa
    assert temperature[-1] == pytest.approx(288.0 * aerosols.SKIN_FACTOR)
    assert np.all(np.diff(temperature) <= 1e-9) and np.all(np.diff(pressure) < 0.0)
    height_km = 8.314 * 288.0 / (EARTH.mean_molar_mass_kg_mol * 9.80665) / 1e3
    lapse = (temperature[0] - temperature[1]) / (EARTH.pressure_pa - pressure[1]) * EARTH.pressure_pa / 288.0
    assert height_km == pytest.approx(8.4, rel=0.05) and 0.2 < lapse < 0.35        # kappa R / cp, about 0.286 for air


def test_the_saturation_curve_of_water_and_the_ceiling_of_a_supercritical_gas():
    assert aerosols.saturation_pa("H2O", 373.0) == pytest.approx(1e5, rel=0.02)
    assert aerosols.saturation_pa("H2O", 273.16) == pytest.approx(611.657, rel=0.001)     # the water triple point, from the checked curve (GEN.452)
    assert aerosols.saturation_pa("NH3", 195.49) == pytest.approx(6053.0, rel=0.03)       # the gas table's triple point
    assert math.isinf(aerosols.saturation_pa("CO2", 400.0))
    assert aerosols.saturation_pa("H2O", 300.0) < aerosols.saturation_pa("H2O", 320.0)


def test_earth_air_makes_a_water_cloud_with_the_vapour_of_a_moist_column():
    condensate = aerosols.condensate_columns_kg_m2(EARTH, 288.0)
    assert set(condensate) == {"H2O"}
    assert 20.0 < condensate["H2O"] < 60.0                    # the column of 1 percent vapour at 288 K above its condensation level
    assert aerosols.condensate_columns_kg_m2(DRY_EARTH, 288.0) == {}


def test_thin_cold_air_makes_a_thin_ice_cloud_and_hot_thick_air_none():
    thin = aerosols.condensate_columns_kg_m2(MARS, 215.0)
    assert set(thin) == {"H2O"} and thin["H2O"] < 0.1
    assert aerosols.condensate_columns_kg_m2(VENUS, 737.0).get("CO2") is None             # supercritical: nothing condenses
    assert aerosols.condensate_columns_kg_m2(TITAN, 94.0)["CH4"] > 10.0


def test_the_cloud_depth_limit_is_three_q_m_over_four_rho_r():
    mass = aerosols.condensate_columns_kg_m2(EARTH, 288.0)["H2O"]
    limit = aerosols.cloud_depth_limit(EARTH, 288.0, 10.0)
    assert limit == pytest.approx(3.0 * 2.0 * mass / (4.0 * 997.0 * 10e-6), rel=1e-9)
    assert aerosols.cloud_depth_limit(EARTH, 288.0, 20.0) == pytest.approx(limit / 2.0)


def test_the_cloud_draw_is_inside_its_bracket_and_zero_without_a_condensate():
    for seed in range(60):
        tau, cover = aerosols.draw_cloud(EARTH, 288.0, random.Random(seed))
        assert aerosols.CLOUD_TAU_FLOOR <= tau <= aerosols.cloud_depth_limit(EARTH, 288.0, 1.0)
        assert 0.1 <= cover <= 1.0
    assert aerosols.draw_cloud(DRY_EARTH, 288.0, random.Random(3)) == (0.0, 0.0)
    depths = sorted(aerosols.draw_cloud(EARTH, 288.0, random.Random(s))[0] for s in range(400))
    assert 5.0 < depths[200] < 200.0                           # a median in the range of Earth's bins (5, 10, 20) and Venus (20 to 30)


def test_the_methane_haze_gate_follows_the_ratio_and_the_hydrogen():
    assert aerosols.haze_gate(EARTH) == 0.0
    assert aerosols.haze_gate(TITAN) == 1.0
    onset = _air({"CO2": 90000.0, "CH4": 9000.0, "N2": 1000.0}, 9.8)
    thick = _air({"CO2": 80000.0, "CH4": 16000.0, "N2": 4000.0}, 9.8)
    middle = _air({"CO2": 85000.0, "CH4": 12750.0, "N2": 2250.0}, 9.8)
    assert aerosols.haze_gate(onset) == pytest.approx(0.0, abs=1e-9)
    assert aerosols.haze_gate(thick) == pytest.approx(1.0)
    assert aerosols.haze_gate(middle) == pytest.approx(0.5, abs=0.01)
    hydrogen = _air({"H2": 50000.0, "CH4": 5000.0, "N2": 45000.0}, 9.8)
    assert aerosols.haze_gate(hydrogen) == pytest.approx(aerosols.HAZE_H2_CUT)


def test_the_titan_like_haze_is_the_order_of_titans_and_never_a_preset():
    depths = np.array([_haze(TITAN, 94.0, s).tau for s in range(300)])
    assert 0.003 < depths.min() < 0.1 and 1.0 < depths.max() < 30.0                  # Titan's 9.6 at 500 nm is a target in the span's order
    assert len(set(np.round(depths, 6))) > 250
    assert _haze(EARTH, 288.0, 1) is aerosols.NO_HAZE
    assert _haze(HOT_JUPITER_LIKE, 300.0, 1) is aerosols.NO_HAZE                    # no methane, no haze


def _haze(column, temperature, seed, star_temperature_k=5772.0):
    return aerosols.haze_in(column, temperature, aerosols.draw_haze_draws(random.Random(seed)), star_temperature_k)


def test_a_haze_takes_a_share_of_the_sunlight_that_grows_with_its_depth_and_sits_aloft():
    hazes = [_haze(TITAN, 94.0, s) for s in range(60)]
    ordered = sorted(hazes, key=lambda haze: haze.tau)
    shares = [haze.absorbed_share for haze in ordered]
    assert shares == sorted(shares) and 0.0 < shares[0] and shares[-1] <= aerosols.HAZE_SHARE_CAP
    assert all(0.0 < haze.pressure_pa < TITAN.pressure_pa for haze in hazes)
    assert _haze(TITAN, 94.0, 2, star_temperature_k=None).absorbed_share == 0.0       # no star, no sunlight to take
    assert aerosols.haze_absorbed_share(0.0, 5772.0) == 0.0
    assert aerosols.haze_absorbed_share(1.0e4, 5772.0) == aerosols.HAZE_SHARE_CAP


def test_a_cooler_star_has_more_of_its_light_where_the_haze_absorbs_less():
    """The tholin absorbs in the blue (Khare): a red star's light is less absorbed than a hot star's."""
    assert aerosols.haze_absorbed_share(1.0, 3000.0) < aerosols.haze_absorbed_share(1.0, 5772.0) < aerosols.haze_absorbed_share(1.0, 9000.0)


def test_a_thick_haze_takes_its_sunlight_higher_up_than_a_thin_one_of_the_same_mass_profile():
    draws = aerosols.draw_haze_draws(random.Random(2))
    thin = aerosols.haze_in(TITAN, 94.0, aerosols.HazeDraws(draws.kzz_m2_s, 0.01, draws.radius_m), 5772.0)
    thick = aerosols.haze_in(TITAN, 94.0, aerosols.HazeDraws(draws.kzz_m2_s, 1.0, draws.radius_m), 5772.0)
    assert thick.tau > thin.tau and thick.pressure_pa < thin.pressure_pa


def test_the_haze_draws_are_a_stream_of_the_bodys_own_and_the_air_makes_the_same_haze_each_time():
    class Body:
        is_moon, mass, radius, distance = False, 5.97e24, 6371.0, 1.0

    first = aerosols.draw_haze_draws(aerosols.haze_stream_of(Body))
    assert first == aerosols.draw_haze_draws(aerosols.haze_stream_of(Body))
    Body.distance = 2.0
    assert first != aerosols.draw_haze_draws(aerosols.haze_stream_of(Body))
    assert _haze(TITAN, 94.0, 7) == _haze(TITAN, 94.0, 7)


def test_the_haze_supply_grows_with_the_mixing_and_not_with_the_methane_share():
    low = aerosols.haze_supply_kg_m2_s(TITAN, 94.0, 1e3)
    assert aerosols.haze_supply_kg_m2_s(TITAN, 94.0, 1e4) == pytest.approx(10.0 * low)
    richer = _air({"N2": 126000.0, "CH4": 14000.0}, 1.35)
    assert aerosols.haze_supply_kg_m2_s(richer, 94.0, 1e3) == pytest.approx(low, rel=0.1)      # the methane splits where its column is one optical depth
    assert aerosols.haze_supply_kg_m2_s(EARTH, 288.0, 1e3) == 0.0


def test_dust_needs_a_dry_ground_and_a_wind_over_the_threshold():
    drawn = [aerosols.draw_dust(True, random.Random(s)) for s in range(600)]
    share = sum(1 for d in drawn if d > 0.0) / len(drawn)
    assert share == pytest.approx(0.5, abs=0.07)                                      # u*/u*_t log-uniform on 0.3 to 3 is above 1 half the time
    assert all(0.5 <= d <= 1.5 for d in drawn if d > 0.0)
    assert all(aerosols.draw_dust(False, random.Random(s)) == 0.0 for s in range(50))


@pytest.mark.parametrize("column, temperature", [(EARTH, 288.0), (MARS, 215.0), (TITAN, 94.0), (None, 0.0)])
def test_every_body_takes_the_same_number_of_draws(column, temperature):
    rng = Counting(11)
    aerosols.draw_aerosols(column, temperature, rng, aerosols.NO_HAZE)
    assert rng.calls == aerosols.DRAWS


def test_the_same_stream_draws_the_same_aerosols():
    first = aerosols.draw_aerosols(EARTH, 288.0, random.Random(5), aerosols.NO_HAZE, dry_surface=False)
    assert first == aerosols.draw_aerosols(EARTH, 288.0, random.Random(5), aerosols.NO_HAZE, dry_surface=False)
    assert first.tau_dust == 0.0 and first.tau_haze == 0.0 and first.tau_cloud > 0.0
    assert aerosols.draw_aerosols(None, 0.0, random.Random(5), aerosols.NO_HAZE) == aerosols.NONE


def test_the_layers_carry_the_materials_optics():
    centres = np.array([300.0, 450.0, 550.0, 650.0, 800.0])
    aerosol = aerosols.Aerosols(tau_haze=2.0, tau_dust=0.0, tau_cloud=10.0, cloud_cover=0.6, aerosol_draw=0.5)
    haze, cloud = aerosols.layers(aerosol, centres)
    assert cloud.tau == pytest.approx(np.full(5, 10.0)) and cloud.g == pytest.approx(np.full(5, 0.85)) and cloud.omega == pytest.approx(np.ones(5))
    assert haze.tau[2] == pytest.approx(2.0)
    assert np.all(np.diff(haze.tau) < 0.0)                                                    # extinction falls with wavelength
    exponent = -math.log(haze.tau[3] / haze.tau[1]) / math.log(650.0 / 450.0)
    assert exponent == pytest.approx(aerosols._haze_exponent(0.5), rel=1e-6) and 0.97 <= exponent <= 2.34
    assert haze.omega[2] == pytest.approx(aerosols.HAZE_OMEGA_550, abs=1e-9)
    assert haze.omega[1] < haze.omega[2] < haze.omega[3]                                       # tholin absorbs blue: darker in the blue
    dust = aerosols.layers(aerosols.Aerosols(0.0, 1.0, 0.0, 0.0, 0.0), centres)[0]
    assert dust.tau == pytest.approx(np.ones(5))
    assert dust.omega[1] == pytest.approx(0.7725, abs=1e-6) and dust.omega[3] == pytest.approx(0.9545, abs=1e-6)       # Wolff 2009, 1.2 micron: 450 and 650 nm
    assert dust.omega[0] < dust.omega[1] < dust.omega[3]
    assert aerosols.layers(aerosols.NONE, centres) == (None, None)


def test_dust_and_haze_are_mixed_by_depth():
    centres = np.array([550.0])
    both, _ = aerosols.layers(aerosols.Aerosols(1.0, 1.0, 0.0, 0.0, 0.0), centres)
    haze, _ = aerosols.layers(aerosols.Aerosols(1.0, 0.0, 0.0, 0.0, 0.0), centres)
    dust, _ = aerosols.layers(aerosols.Aerosols(0.0, 1.0, 0.0, 0.0, 0.0), centres)
    assert both.tau[0] == pytest.approx(2.0)
    scattering = haze.tau * haze.omega + dust.tau * dust.omega
    assert both.omega[0] == pytest.approx(scattering[0] / 2.0)
    assert min(haze.g[0], dust.g[0]) <= both.g[0] <= max(haze.g[0], dust.g[0])


def test_a_cloud_dims_the_ground_and_brightens_the_planet():
    flux = ss.bolometric_flux_w_m2(constants.SOLAR_LUMINOSITY, 1.0)
    options = dict(temperature_k=5772.0, flux_w_m2=flux, uv_c_share=0.256, ground_albedo=0.15, mu0=0.8)
    with fetch.disregard_missing_data():
        clear = ar.surface_radiation(DRY_EARTH, **options)
        centres = clear.centres_nm
        _, cloud = aerosols.layers(aerosols.Aerosols(0.0, 0.0, 10.0, 1.0, 0.5), centres)
        overcast = ar.surface_radiation(DRY_EARTH, low_aerosol=cloud, cloud_cover=1.0, **options)
        half = ar.surface_radiation(DRY_EARTH, low_aerosol=cloud, cloud_cover=0.5, **options)
    visible = lambda s: s.band_ground_w_m2("visible")
    assert visible(overcast) < visible(half) < visible(clear)
    assert overcast.direct_w_m2.sum() < 0.1 * clear.direct_w_m2.sum()
    assert overcast.outgoing_w_m2.sum() > clear.outgoing_w_m2.sum()


def _body(**changes):
    from types import SimpleNamespace
    values = dict(is_moon=False, mass=5.97e24, radius=6371.0, distance=1.0, surface_temperature=288.0, atmospheric_pressure=101325.0,
                  planet_class="M")
    return SimpleNamespace(**{**values, **changes})


def test_the_aerosol_stream_is_keyed_on_the_bodys_identity_and_conditions_never_its_class_letter():
    """GEN.917, GEN.915: two bodies that differ only in their class letter draw the same aerosols; one that differs in a condition
    the draws read (the air pressure here, the surface temperature there) draws from another stream."""
    def first(body):
        return aerosols.stream_of(body).random()

    assert first(_body(planet_class="M")) == first(_body(planet_class="D")) == first(_body(planet_class=None))
    assert first(_body()) != first(_body(atmospheric_pressure=5e5))
    assert first(_body()) != first(_body(surface_temperature=250.0))
