"""
The grey infrared optical depth of an air from its mixture (GEN.703, `planetgen.physics.infrared_opacity`). Earth, Venus, Mars,
Titan and the Wordsworth and Pierrehumbert scenarios are test targets from docs/design/composition/atmosphere-composition-parts/
03-infrared-greenhouse-thermal.md, never inputs of the code.
"""

import math

import numpy as np
import pytest

from planetgen.physics import gases, mixture
from planetgen.physics import infrared_opacity as ir

EARTH = {"N2": 78084.0, "O2": 20946.0, "Ar": 934.0, "CO2": 41.0, "H2O": 1300.0}
"""dict: Earth's air with the vapour pressure of its surface (0.77 of saturation at 288 K, 1.3 kPa)."""
VENUS = {"CO2": 0.965 * 92.0e5, "N2": 0.035 * 92.0e5}
MARS = {"CO2": 600.0 * 0.953, "N2": 600.0 * 0.027, "Ar": 600.0 * 0.016}
TITAN = {"N2": 1.47e5 * 0.95, "CH4": 1.47e5 * 0.05}


def _part(mix, temperature, gravity, part):
    return ir.optical_depth(mix, temperature, gravity).by_part[part]


def test_planck_fraction_matches_a_quadrature():
    xs = np.linspace(1e-9, 60.0, 600001)
    density = xs ** 3 / np.expm1(xs)
    cumulative = np.cumsum((density[1:] + density[:-1]) / 2.0 * np.diff(xs))
    cumulative = np.concatenate([[0.0], cumulative]) * 15.0 / math.pi ** 4
    for x in (0.1, 0.5, 1.0, 1.99, 2.01, 3.0, 5.0, 10.0):
        index = int(np.searchsorted(xs, x))
        assert float(ir.planck_fraction_below(x)) == pytest.approx(float(cumulative[index]), abs=2e-5)
    assert float(ir.planck_fraction_below(0.0)) == 0.0
    assert float(ir.planck_fraction_below(80.0)) == pytest.approx(1.0, abs=1e-12)


def test_the_weights_are_a_distribution_and_peak_where_wien_says():
    for temperature in (50.0, 288.0, 735.0, 3000.0):
        weights = ir.planck_weights(temperature)
        assert weights.sum() == pytest.approx(1.0, abs=1e-12)
        assert (weights >= 0.0).all()
    weights = ir.planck_weights(288.0)
    peak = ir.CENTRES_CM1[int(np.argmax(weights))]
    assert peak == pytest.approx(2.821 * 288.0 / ir.SECOND_RADIATION_CM_K, abs=ir.BIN_CM1)


def test_a_saturated_band_widens_with_the_logarithm_of_its_column():
    """The logarithm of the carbon dioxide forcing law: each doubling of a saturated band's gas adds the same."""
    steps = []
    for pressure in (10.0, 20.0, 40.0, 80.0):
        lower = _part({"N2": 1.0e5, "CO2": pressure}, 288.0, 9.81, "bands")
        upper = _part({"N2": 1.0e5, "CO2": 2.0 * pressure}, 288.0, 9.81, "bands")
        steps.append(upper - lower)
    assert min(steps) > 0.0
    totals = [_part({"N2": 1.0e5, "CO2": pressure}, 288.0, 9.81, "bands") for pressure in (10.0, 20.0, 40.0, 80.0)]
    assert all(step < total for step, total in zip(steps, totals))          # twice the gas adds less than the gas had: it saturates
    assert max(steps) < 2.0 * min(steps)            # and the add stays of one size over a factor 8 of column (a logarithm, not a power)


def test_a_doubling_of_carbon_dioxide_gives_the_measured_forcing_of_a_doubling():
    """F = 5.32 ln 2 = 3.7 W/m2 (Byrne and Goldblatt 2014a, READ): the wing e-fold is set from the logarithm of this law, and the
    change in outgoing flux from the extra width of the bands, between the surface and the stratosphere (216.65 K), is the same size.
    The doubling is taken from 280 ppm (the law is the one of that state, 200 to 10000 ppm), where the forcing is a few percent apart."""
    surface, stratosphere = 288.0, 216.65
    flux_by_bin = 5.670374e-8 * (surface ** 4 * ir.planck_weights(surface) - stratosphere ** 4 * ir.planck_weights(stratosphere))
    for fraction in (140e-6, 280e-6, 560e-6, 1120e-6):
        column = fraction * 1.0e5 * gases.gas("CO2").molar_mass / 0.0289 / 9.81
        change = np.exp(-ir._band_depth("CO2", column, 1.0)) - np.exp(-ir._band_depth("CO2", 2.0 * column, 1.0))
        assert float((flux_by_bin * change).sum()) == pytest.approx(3.7, rel=0.25)


def test_bands_above_the_planck_peak_do_nothing_on_a_cold_world_and_much_on_a_hot_one():
    cold = _part(TITAN, 94.0, 1.35, "bands")
    hot = _part({"N2": 1.0e5, "CH4": 5.0e3}, 700.0, 9.81, "bands")
    assert cold < 1e-3
    assert hot > 0.2


def test_a_pure_carbon_dioxide_air_is_of_order_unity_at_one_bar_and_carbon_monoxide_is_thin():
    """Koll 2019's statements, which set the far-wing coefficients: CO2 about 3 at 1 bar, CO thin even at 1 bar."""
    co2 = _part({"CO2": 1.0e5}, 288.0, 9.81, "far_wing")
    co = ir.optical_depth({"CO": 1.0e5}, 288.0, 9.81).tau0
    assert 2.0 < co2 < 4.0
    assert co < 1.0


def test_a_dense_carbon_dioxide_air_is_a_runaway_like_greenhouse_with_no_ceiling():
    taus = [ir.optical_depth({"CO2": pressure}, 700.0, 8.87).tau0 for pressure in (1.0e5, 1.0e6, 1.0e7)]
    assert taus == sorted(taus)
    assert taus[-1] > 1.0e3
    assert all(math.isfinite(t) for t in taus)


def test_the_water_continuum_closes_the_window_of_a_humid_air_and_grows_faster_than_the_vapour():
    window = int(np.argmin(np.abs(ir.CENTRES_CM1 - 1000.0)))

    def depth(fraction, temperature):
        partials = {"H2O": fraction * 1.0e5, "N2": (1.0 - fraction) * 1.0e5}
        molar_mass = mixture.mean_molar_mass_kg_mol(partials)
        vapour = ir.vapour_profile(partials, temperature)
        return float(ir._water_continuum_depth(vapour, molar_mass, 9.81)[window])

    assert depth(2.0e-3, 288.0) / depth(1.0e-3, 288.0) > 2.5      # the self continuum goes as the vapour squared; the foreign continuum as the vapour, which dilutes the rise
    assert depth(4.0e-3, 270.0) > depth(4.0e-3, 310.0)       # the same vapour absorbs more when colder (the self temperature exponent)
    assert 0.02 < depth(1.0e-2, 288.0) < 2.0                  # a humid window at 10 um: a tenth of a unit, as measured
    assert _part({"N2": 1.0e5}, 288.0, 9.81, "water_continuum") == 0.0


@pytest.mark.parametrize("label,mix,temperature,gravity,expected,tolerance", [
    ("3.3 bar 90 percent N2 + 10 percent H2 (Wordsworth and Pierrehumbert 2013)", {"N2": 0.9 * 3.3e5, "H2": 0.1 * 3.3e5}, 273.0, 9.81, 1.26, 0.25),
    ("1.1 bar of the same mix", {"N2": 0.9 * 1.1e5, "H2": 0.1 * 1.1e5}, 273.0, 9.81, 0.25, 0.3),
    ("1 bar N2, 288 K", {"N2": 1.0e5}, 288.0, 9.81, 0.040, 0.5),
    ("3 bar CO2, 90 percent with 10 percent H2, 250 K (Ramirez 2014)", {"CO2": 0.9 * 3.0e5, "H2": 0.1 * 3.0e5}, 250.0, 3.71, 2.4, 0.5),
    ("10 bar of 85 percent H2 in He, 300 K (hycean-like)", {"H2": 0.85e6, "He": 0.15e6}, 300.0, 12.0, 6.5 + 3.6, 0.5),
])
def test_collision_induced_depths_match_the_table_near_the_planck_peak(label, mix, temperature, gravity, expected, tolerance):
    assert _part(mix, temperature, gravity, "collision_induced") == pytest.approx(expected, rel=tolerance), label


def test_collision_induced_depths_are_within_a_factor_of_three_for_the_cold_cases():
    """Titan N2 at 94 K: the table gives 1.24 (Koll's definition); the saturating form gives 3.3. A thin Mars-like CIA is exact."""
    titan = _part({"N2": 1.47e5}, 94.0, 1.35, "collision_induced")
    assert 1.24 / 3.0 < titan < 1.24 * 3.0
    assert _part({"CO2": 6.4e2}, 220.0, 3.71, "collision_induced") == pytest.approx(8.5e-5, rel=0.5)


def test_pair_coefficients_come_from_the_table_and_a_missing_pair_is_none():
    assert ir.pair_coefficient(("N2", "N2"), 300.0) == pytest.approx(1.9e-46)
    assert ir.pair_coefficient(("N2", "H2"), 300.0) == ir.pair_coefficient(("H2", "N2"), 300.0)
    assert ir.pair_coefficient(("N2", "O2"), 300.0) is None
    assert 3.9e-46 > ir.pair_coefficient(("N2", "N2"), 225.0) > 2.6e-46


def test_earth_venus_mars_and_titan_land_on_their_equivalent_depths():
    """The targets of part 03 section 5.4: Earth about 1.8, Venus about 2e5 (n 2) or 400 (n 1), Mars 0.1 to 0.25, Titan 1.6 to 5.6."""
    assert ir.optical_depth(EARTH, 288.0, 9.81).tau0 == pytest.approx(1.82, rel=0.35)
    assert 0.05 < ir.optical_depth(MARS, 215.0, 3.71).tau0 < 0.3
    venus = ir.optical_depth(VENUS, 735.0, 8.87).tau0
    assert 1.0e4 < venus < 1.0e6
    assert 1.0 < ir.optical_depth(TITAN, 94.0, 1.35).tau0 < 12.0


def test_the_parts_and_the_gases_explain_the_whole():
    opacity = ir.optical_depth(EARTH, 288.0, 9.81)
    assert set(opacity.by_part) == {"bands", "far_wing", "water_continuum", "collision_induced", "clouds"}
    assert opacity.tau0 <= sum(opacity.by_part.values()) + 1e-9           # overlap only ever lowers the sum
    assert opacity.tau0 >= max(opacity.by_part.values()) - 1e-12
    assert set(opacity.by_gas) == {"H2O", "CO2"}
    assert 0.0 < opacity.condensable_share < 1.0
    assert opacity.condensable_share == pytest.approx(opacity.by_gas["H2O"] / sum(opacity.by_gas.values()))


def test_lapse_scaling_runs_from_condensable_to_dry_and_hydrogen():
    wet = ir.lapse_scaling(ir.optical_depth({"H2O": 2.0e4, "N2": 8.0e4}, 288.0, 9.81))
    dry = ir.lapse_scaling(ir.optical_depth({"CO2": 1.0e5}, 288.0, 9.81))
    hydrogen = ir.lapse_scaling(ir.optical_depth({"H2": 8.6e4, "He": 1.4e4, "CH4": 20.0}, 120.0, 24.8))
    assert 0.6 <= wet < 0.75
    assert dry == ir.LAPSE_SCALING_DRY
    assert hydrogen == ir.LAPSE_SCALING_HYDROGEN


def test_a_gas_with_no_listed_band_still_absorbs_in_the_fingerprint_region():
    listed = {row.formula for row in gases.GASES if row.ir == "ir" and not (row.infrared and row.infrared.bands)}
    assert listed                     # the table has such gases
    formula = min(listed)
    assert ir.optical_depth({"N2": 1.0e5, formula: 1.0e3}, 300.0, 9.81).by_part["bands"] > 0.0


def test_no_gas_of_the_table_gives_a_nan_or_a_negative_depth():
    for row in gases.GASES:
        for temperature in (60.0, 300.0, 1500.0):
            opacity = ir.optical_depth({"N2": 1.0e5, row.formula: 1.0e3}, temperature, 9.81)
            assert math.isfinite(opacity.tau0) and opacity.tau0 >= 0.0, (row.formula, temperature)
            assert all(math.isfinite(v) and v >= 0.0 for v in opacity.by_part.values()), (row.formula, temperature)


def test_a_very_thick_air_does_not_underflow():
    opacity = ir.optical_depth({"CO2": 9.2e9}, 735.0, 8.87)
    assert math.isfinite(opacity.tau0) and opacity.tau0 > 1.0e5


def test_bad_input_is_refused():
    with pytest.raises(ValueError):
        ir.optical_depth({}, 288.0, 9.81)
    with pytest.raises(ValueError):
        ir.optical_depth(EARTH, 0.0, 9.81)
    with pytest.raises(ValueError):
        ir.optical_depth(EARTH, 288.0, -1.0)


def test_the_vapour_thins_with_height_and_holds_about_the_measured_column_for_the_earth():
    """Manabe and Wetherald's humidity falls linearly with pressure; the Earth's column of vapour is about 25 kg/m2 (precipitable water)."""
    vapour = ir.vapour_profile(EARTH, 288.0)
    assert vapour.surface_fraction == pytest.approx(1300.0 / sum(EARTH.values()))
    assert (np.diff(vapour.fraction) <= 1e-15).all() and vapour.fraction[-1] < 1e-3 * vapour.surface_fraction
    molar_mass = mixture.mean_molar_mass_kg_mol(EARTH)
    column = float(np.sum(vapour.fraction * vapour.pressure_pa * vapour.log_step)) * gases.gas("H2O").molar_mass / molar_mass / 9.81
    assert 15.0 < column < 35.0
    well_mixed = vapour.surface_fraction * sum(EARTH.values()) * gases.gas("H2O").molar_mass / molar_mass / 9.81
    assert well_mixed > 2.0 * column                 # a well-mixed column of the surface's vapour would hold three times the water


def test_a_vapour_over_a_surface_above_the_critical_temperature_stands_as_it_is_until_the_cold_condenses_it():
    mix = {"CO2": 9.0e6, "H2O": 270.0}
    vapour = ir.vapour_profile(mix, 735.0)
    assert vapour.cloud_cover == 0.0
    assert vapour.fraction[0] == pytest.approx(vapour.surface_fraction)
    assert vapour.fraction[0] >= vapour.fraction[-1]


def test_an_air_with_no_vapour_has_no_profile_no_continuum_and_no_cloud():
    opacity = ir.optical_depth({"N2": 1.0e5, "CO2": 40.0}, 288.0, 9.81)
    assert ir.vapour_profile({"N2": 1.0e5}, 288.0) is None
    assert opacity.by_part["water_continuum"] == 0.0 and opacity.by_part["clouds"] == 0.0


def test_the_cloud_scales_with_the_surface_humidity_and_adds_its_depth_to_the_rest():
    saturation_pa = ir.saturation.pressure_pa("H2O", 288.0)
    earth = ir.optical_depth(dict(EARTH, H2O=0.77 * saturation_pa), 288.0, 9.81)
    half = ir.optical_depth(dict(EARTH, H2O=0.385 * saturation_pa), 288.0, 9.81)
    assert earth.by_part["clouds"] == pytest.approx(ir.CLOUD_LONGWAVE_TAU0)
    assert half.by_part["clouds"] == pytest.approx(ir.CLOUD_LONGWAVE_TAU0 / 2.0)
    assert earth.tau0 > earth.by_part["clouds"] > 0.0


def test_the_layered_continuum_converges_with_the_number_of_layers(monkeypatch):
    partials = dict(EARTH)
    window = int(np.argmin(np.abs(ir.CENTRES_CM1 - 1000.0)))
    molar_mass = mixture.mean_molar_mass_kg_mol(partials)
    coarse = float(ir._water_continuum_depth(ir.vapour_profile(partials, 288.0), molar_mass, 9.81)[window])
    monkeypatch.setattr(ir, "VAPOUR_LAYERS", 2 * ir.VAPOUR_LAYERS)
    fine = float(ir._water_continuum_depth(ir.vapour_profile(partials, 288.0), molar_mass, 9.81)[window])
    assert coarse == pytest.approx(fine, rel=5e-3)


def test_a_vapour_over_a_surface_so_cold_that_its_saturation_pressure_underflows_is_all_condensate_and_gives_a_profile():
    vapour = ir.vapour_profile({"N2": 1.0e5, "H2O": 1.0e-3}, 5.0)
    assert ir.saturation.pressure_pa("H2O", 5.0) == 0.0
    assert np.isfinite(vapour.fraction).all() and vapour.fraction.max() <= vapour.surface_fraction
    assert math.isfinite(ir.optical_depth({"N2": 1.0e5, "H2O": 1.0e-3}, 5.0, 9.81).tau0)
