"""
The one-dimensional radiative-convective model (GEN.715, GEN.717, GEN.718, `planetgen.physics.radiative_convective`) on the offline fixture
tables of `tests/ir_datasets_support.py`: the water curve, the two-stream flux sweep, the adiabat (dry and moist), the humidity profile and
the balance the solver finds. The numbers of Earth, Venus, Mars and Titan against the real tables are `test_radiative_convective_real.py`.
"""

import math

import numpy as np
import pytest

from planetgen.datasets import fetch
from planetgen.physics import constants
from planetgen.physics import line_opacity as ir
from planetgen.physics import radiative_convective as rc
from tests import ir_datasets_support as support

SIGMA = constants.STEFAN_BOLTZMANN_CONSTANT
CARBON = {"N2": 9.0e4, "CO2": 1.0e4}
EARTHLIKE = {"N2": 7.8e4, "O2": 2.1e4, "CO2": 41.0, "H2O": 1.0e3}


@pytest.fixture(autouse=True)
def tables(tmp_path, monkeypatch):
    support.install(monkeypatch, tmp_path / "cache")
    yield
    fetch._STAND_INS.clear()


# ---- water ----

def test_water_boils_at_a_bar_and_freezes_near_its_triple_point():
    assert float(rc.saturation_pressure_pa(373.15)) == pytest.approx(101325.0, rel=0.02)
    assert float(rc.saturation_pressure_pa(273.16)) == pytest.approx(611.7, rel=0.02)
    assert float(rc.saturation_pressure_pa(300.0)) == pytest.approx(3536.0, rel=0.03)


def test_the_saturation_curve_rises_with_temperature_takes_arrays_and_is_infinite_past_the_critical_point():
    values = rc.saturation_pressure_pa(np.array([250.0, 300.0, 350.0, 400.0]))
    assert values.shape == (4,) and np.all(np.diff(values) > 0.0)
    assert math.isinf(float(rc.saturation_pressure_pa(rc.WATER_TABLE_TOP_K + 10.0)))
    assert np.isfinite(float(rc.saturation_pressure_pa(5.0)))


def test_the_latent_heat_of_water_follows_clausius_and_clapeyron():
    assert float(rc.latent_heat_j_mol(373.15)) == pytest.approx(40.66e3, rel=0.05)
    assert float(rc.latent_heat_j_mol(300.0)) == pytest.approx(43.99e3, rel=0.05)
    assert float(rc.latent_heat_j_mol(300.0)) > float(rc.latent_heat_j_mol(400.0))


# ---- the flux ----

def _isothermal(temperature, layers=6):
    return np.full((1, layers + 1), float(temperature))


def test_a_clear_column_sends_the_ground_flux_through_and_gets_nothing_back():
    emitted = ir.band_emission(np.array([300.0]))[0].sum()
    up, down = rc.sweep(_isothermal(300.0), np.zeros((6, ir.BANDS, ir.GRID_POINTS)), np.array([300.0]))
    assert up == pytest.approx(emitted) and not down.any()
    assert emitted == pytest.approx(SIGMA * 300.0 ** 4, rel=0.01)


def test_an_isothermal_column_in_equilibrium_with_the_ground_keeps_the_same_upward_flux_at_every_depth():
    rng = np.random.default_rng(5)
    depth = rng.lognormal(0.0, 1.5, (6, ir.BANDS, ir.GRID_POINTS))
    emitted = ir.band_emission(np.array([260.0]))[0].sum()
    up, down = rc.sweep(_isothermal(260.0), depth, np.array([260.0]))
    assert up[0] == pytest.approx(np.full(7, emitted), rel=1e-9)
    assert down[0, 0] == 0.0 and np.all(np.diff(down[0]) >= 0.0) and down[0, -1] <= emitted * (1.0 + 1e-9)


def test_a_thick_column_radiates_from_its_top_and_returns_what_it_emits_at_the_bottom():
    thick = np.full((6, ir.BANDS, ir.GRID_POINTS), 200.0)
    levels = np.linspace(200.0, 300.0, 7)[None, :]
    up, down = rc.sweep(levels, thick, np.array([400.0]))
    assert up[0, 0] == pytest.approx(ir.band_emission(np.array([levels[0, 0]]))[0].sum(), rel=0.12)          # the top layer's own temperature
    assert down[0, -1] == pytest.approx(ir.band_emission(np.array([levels[0, -1]]))[0].sum(), rel=0.08)       # the bottom layer's


def test_a_warmer_ground_sends_more_up_and_more_depth_lets_less_out():
    depth = np.full((6, ir.BANDS, ir.GRID_POINTS), 0.5)
    cold, _ = rc.sweep(_isothermal(250.0), depth, np.array([250.0]))
    warm, _ = rc.sweep(_isothermal(250.0), depth, np.array([300.0]))
    thick, _ = rc.sweep(_isothermal(250.0), 4.0 * depth, np.array([300.0]))
    assert warm[0, 0] > cold[0, 0] and warm[0, 0] > thick[0, 0]


# ---- the adiabat ----

def test_a_dry_adiabat_follows_the_heat_capacity_of_its_gas():
    column = rc.Column({"N2": 1.0e5}, 9.81)
    adiabat = column.adiabat(300.0)
    assert adiabat.shape == (column.layers,) and np.all(np.diff(adiabat) > 0.0) and adiabat[-1] < 300.0
    kappa = math.log(300.0 / adiabat[-1]) / math.log(column.surface_pressure_pa / column.level_pressure[-2])
    assert kappa == pytest.approx(2.0 / 7.0, rel=0.02)          # R / cp of a diatomic gas at 300 K


def test_the_more_atoms_a_molecule_has_the_less_the_air_cools_with_height():
    top = {gas: rc.Column({gas: 1.0e5}, 9.81).adiabat(300.0)[0] for gas in ("He", "N2", "CO2")}
    assert top["He"] < top["N2"] < top["CO2"]          # R / cp is 0.4, 0.29 and 0.22


def test_a_saturated_air_cools_more_slowly_with_height_than_a_dry_one():
    moist = rc.Column({"N2": 8.0e4, "H2O": 2.0e4}, 9.81)
    dry = rc.Column({"N2": 1.0e5}, 9.81)
    assert moist._saturated(300.0, 1.0e5) and not moist._saturated(340.0, 1.0e5)
    wet_profile, dry_profile = moist.adiabat(300.0), dry.adiabat(300.0)
    assert np.all(wet_profile[10:] > dry_profile[10:] - 1e-6)          # latent heat releases as the vapour condenses
    assert wet_profile[-1] > dry_profile[-1]
    assert np.all(np.isfinite(wet_profile)) and np.all(np.diff(wet_profile) > 0.0)


def test_methane_condenses_in_a_cold_hydrogen_air_and_slows_its_cooling_with_height():
    """GEN.1030: with no water in the air the condensable is methane, whose condensation near 80 K (the clouds of Uranus and Neptune) shallows the adiabat."""
    methane = rc.Column({"H2": 0.80e5, "He": 0.175e5, "CH4": 0.025e5}, 11.15)
    plain = rc.Column({"H2": 0.8205e5, "He": 0.1795e5}, 11.15)
    assert methane.condensable_surface == pytest.approx(0.025) and plain.condensable_surface == 0.0
    assert methane._saturated(80.0, 1.0e5) and not methane._saturated(100.0, 1.0e5)
    wet_profile, dry_profile = methane.adiabat(80.0), plain.adiabat(80.0)
    assert wet_profile[-1] > dry_profile[-1] + 1.0e-6 and np.all(np.isfinite(wet_profile)) and np.all(np.diff(wet_profile) > 0.0)


def test_water_takes_the_adiabat_before_methane_in_an_air_with_both():
    both = rc.Column({"N2": 0.7e5, "H2O": 0.2e5, "CH4": 0.1e5}, 9.81)
    assert both.condensable_surface == pytest.approx(0.2) and both._saturation is rc.saturation_pressure_pa


def test_the_methane_curve_is_the_saturation_module_and_ends_below_the_critical_point():
    from planetgen.physics import saturation
    assert float(rc.methane_saturation_pressure_pa(100.0)) == pytest.approx(saturation.pressure_pa("CH4", 100.0), rel=1.0e-3)
    assert np.isinf(rc.methane_saturation_pressure_pa(rc.METHANE_TABLE_TOP_K + 1.0))
    assert float(rc.methane_latent_heat_j_mol(110.0)) == pytest.approx(8.2e3, rel=0.05)          # about 8.2 kJ/mol at its boiling point


def test_cold_hydrogen_in_ortho_and_para_equilibrium_holds_more_heat_than_hot_hydrogen():
    cp = rc.equilibrium_hydrogen_cp_over_r
    assert float(cp(5.0)) == pytest.approx(2.5, abs=1.0e-6)          # all para, no rotation: translation alone
    assert 4.3 < float(cp(60.0)) < 4.5                                # the rotational bump of the ortho-para conversion
    assert float(cp(60.0)) > float(cp(100.0)) > float(cp(200.0))
    assert float(cp(2000.0)) == pytest.approx(3.5, abs=0.02)          # rotation fully excited
    assert float(cp(200.0)) == pytest.approx(rc.species_thermo.gas("H2").cp_over_r(200.0), rel=0.01)          # joins the NASA fit


def test_the_heat_capacity_of_hydrogen_on_the_grid_is_the_equilibrium_one_below_200_k_and_the_fit_above():
    values = rc._cp_over_r("H2")
    cold = rc._CP_GRID < rc.EQUILIBRIUM_HYDROGEN_BELOW_K
    assert values[cold].max() > 4.0 and np.all(np.diff(values[~cold][:20]) > 0.0)
    assert values[~cold][0] == pytest.approx(rc.species_thermo.gas("H2").cp_over_r(200.0))
    assert np.all(rc._cp_over_r("N2") == rc._cp_over_r.__wrapped__("N2"))          # no other gas is touched


def test_the_adiabat_is_continuous_and_rises_with_the_surface_temperature_across_the_saturation_corner():
    """A step that started on the dry branch and crossed into saturation used to take the saturated slope for its whole length (the top of
    the column came out 10 K too warm just below the corner and then fell). The slope of the top temperature has a kink at the corner, nothing more."""
    column = rc.Column({"N2": 8.0e4, "H2O": 2.0e4}, 9.81)
    temperatures = np.linspace(318.0, 340.0, 89)          # the vapour saturates at the surface below about 333 K
    top = np.array([column.adiabat(t)[0] for t in temperatures])
    assert np.all(np.diff(top) > 0.0) and np.max(np.diff(top)) < 3.0


# ---- the humidity ----

def test_the_vapour_follows_manabe_and_wetherald_and_is_cold_trapped():
    column = rc.Column(EARTHLIKE, 9.81)
    levels = np.linspace(215.0, 288.0, column.layers + 1)
    layer = 0.5 * (levels[1:] + levels[:-1])
    water = column.water_profile(layer, 288.0)
    assert water.shape == (column.layers,) and np.all(water >= 0.0) and np.all(np.diff(water) >= -1e-18)
    assert water[-1] <= column.water_surface * 1.0001           # the surface holds what the mixture gives when it is unsaturated
    assert water[0] < 1.0e-3 * water[-1]                         # the cold trap and the humidity going to nothing aloft
    saturated = rc.saturation_pressure_pa(layer)
    assert np.all(water * column.layer_pressure <= saturated * (1.0 + 1e-9))


def test_a_dry_air_holds_no_vapour_and_a_pure_steam_air_is_all_vapour():
    layer = np.linspace(200.0, 300.0, 40)
    assert not rc.Column({"N2": 1.0e5}, 9.81).water_profile(layer, 300.0).any()
    assert np.all(rc.Column({"H2O": 1.0e6}, 9.81).water_profile(layer, 600.0) == 1.0)


def test_the_vapour_can_be_held_fixed_as_a_whole_profile():
    held = np.linspace(1.0e-6, 1.0e-2, rc.LAYERS)
    column = rc.Column(EARTHLIKE, 9.81, water_fraction=held)
    assert column.water_profile(np.full(rc.LAYERS, 250.0), 288.0) is column.fixed_water
    assert column.water_profile(np.full(rc.LAYERS, 250.0), 288.0) == pytest.approx(held)


# ---- the cloud ----

def test_the_cloud_keeps_in_a_share_of_the_flux_in_proportion_to_the_surface_humidity_and_is_nil_without_vapour():
    from planetgen.physics import saturation

    assert rc.Column(CARBON, 9.81).cloud_share(288.0) == 0.0
    at_earths = rc.Column({"N2": 1.0e5 - 1.0e3, "H2O": 1.0e3}, 9.81)
    humidity = 1.0e3 / saturation.pressure_pa("H2O", 288.0)
    assert at_earths.cloud_share(288.0) == pytest.approx(rc.CLOUD_LONGWAVE_SHARE * humidity / rc.EARTH_HUMIDITY)
    assert at_earths.cloud_share(300.0) < at_earths.cloud_share(288.0)                 # a drier surface, less cloud
    assert at_earths.cloud_share(200.0) == pytest.approx(rc.CLOUD_LONGWAVE_SHARE / rc.EARTH_HUMIDITY)     # a saturated surface, the most


def test_the_cloud_warms_the_ground_because_it_keeps_in_a_share_of_the_flux(monkeypatch):
    clear = rc.solve(EARTHLIKE, 9.81, 240.0)
    monkeypatch.setattr(rc, "CLOUD_LONGWAVE_SHARE", 0.0)
    bare = rc.solve(EARTHLIKE, 9.81, 240.0)
    assert clear.converged and bare.converged
    assert clear.surface_temperature_k > bare.surface_temperature_k + 3.0


# ---- the sunlight aloft ----

def test_the_sunlight_aloft_is_a_share_of_the_flux_centred_on_the_haze():
    column = rc.Column(CARBON, 9.81)
    heating = rc.sunlight_aloft(column, 100.0, 0.4, 2.0e3)
    assert heating.sum() == pytest.approx(40.0)
    assert column.layer_pressure[int(np.argmax(heating))] == pytest.approx(2.0e3, rel=0.15)
    assert not rc.sunlight_aloft(column, 100.0, 0.0, 2.0e3).any()
    for share in (-0.1, 1.0, 1.5):
        with pytest.raises(ValueError):
            rc.sunlight_aloft(column, 100.0, share, 2.0e3)


# ---- the balance ----

def test_an_air_with_nothing_that_absorbs_in_the_infrared_has_the_effective_temperature():
    solution = rc.solve({"Ar": 1.0e5}, 9.81, 240.0)
    effective = (240.0 / SIGMA) ** 0.25
    assert solution.converged and solution.surface_temperature_k == pytest.approx(effective)
    assert np.all(solution.level_temperature_k == effective) and solution.outgoing_w_m2 == 240.0


def test_the_solution_sends_out_what_it_absorbs_and_is_warmer_than_without_air():
    solution = rc.solve(CARBON, 9.81, 240.0)
    effective = (240.0 / SIGMA) ** 0.25
    assert solution.converged
    assert solution.outgoing_w_m2 == pytest.approx(240.0, rel=0.02)
    assert solution.surface_temperature_k > effective
    assert solution.level_pressure_pa.shape == solution.level_temperature_k.shape == (rc.LAYERS + 1,)
    assert np.all(np.isfinite(solution.level_temperature_k)) and solution.level_temperature_k.min() > 0.0
    assert solution.level_pressure_pa[0] < solution.tropopause_pressure_pa < solution.level_pressure_pa[-1]
    assert solution.level_temperature_k[-1] == solution.surface_temperature_k


def test_below_the_tropopause_the_profile_is_the_adiabat_and_above_it_the_air_is_colder_than_the_tropopause():
    solution = rc.solve(CARBON, 9.81, 240.0)
    column = rc.Column(CARBON, 9.81)
    adiabat = column.adiabat(solution.surface_temperature_k)
    below = solution.level_pressure_pa[:-1] >= solution.tropopause_pressure_pa
    assert solution.level_temperature_k[:-1][below] == pytest.approx(adiabat[below], rel=1e-6)
    assert solution.level_temperature_k[0] < solution.tropopause_temperature_k


def test_more_gas_and_more_sunlight_each_warm_the_surface():
    base = rc.solve({"N2": 9.0e4, "CO2": 1.0e3}, 9.81, 240.0).surface_temperature_k
    more_gas = rc.solve({"N2": 9.0e4, "CO2": 1.0e4}, 9.81, 240.0).surface_temperature_k
    more_sun = rc.solve({"N2": 9.0e4, "CO2": 1.0e3}, 9.81, 300.0).surface_temperature_k
    assert more_gas > base and more_sun > base


def test_a_share_of_the_sunlight_taken_high_up_cools_the_ground():
    clear = rc.solve(CARBON, 9.81, 240.0)
    hazy = rc.solve(CARBON, 9.81, 240.0, aloft_share=0.5, aloft_pressure_pa=3.0e3)
    assert hazy.converged and hazy.outgoing_w_m2 == pytest.approx(240.0, rel=0.02)
    assert hazy.surface_temperature_k < clear.surface_temperature_k          # the anti-greenhouse of a haze


def test_a_solution_run_with_its_own_vapour_held_fixed_gives_the_same_surface():
    first = rc.solve(EARTHLIKE, 9.81, 240.0)
    again = rc.solve(EARTHLIKE, 9.81, 240.0, water_fraction=first.water_fraction)
    assert again.surface_temperature_k == pytest.approx(first.surface_temperature_k, abs=0.3)


def test_a_different_number_of_layers_gives_the_same_physics():
    coarse = rc.solve(CARBON, 9.81, 240.0, layers=24)
    fine = rc.solve(CARBON, 9.81, 240.0, layers=48)
    assert coarse.level_temperature_k.shape == (25,) and fine.level_temperature_k.shape == (49,)
    assert coarse.surface_temperature_k == pytest.approx(fine.surface_temperature_k, abs=3.0)


def test_bad_input_is_refused():
    for flux in (0.0, -5.0, math.nan, math.inf):
        with pytest.raises(ValueError):
            rc.solve(CARBON, 9.81, flux)
    with pytest.raises(ValueError):
        rc.solve(CARBON, 0.0, 240.0)
    with pytest.raises(ValueError):
        rc.solve({}, 9.81, 240.0)
    with pytest.raises(ValueError, match="aloft_pressure_pa"):
        rc.solve(CARBON, 9.81, 240.0, aloft_share=0.3)


def test_a_missing_dataset_refuses_the_solve_and_a_disregarding_run_takes_the_stand_ins(tmp_path, monkeypatch):
    support.forget(monkeypatch, tmp_path / "empty")
    with pytest.raises(fetch.DataMissing):
        rc.solve(CARBON, 9.81, 240.0)
    with fetch.disregard_missing_data():
        solution = rc.solve(CARBON, 9.81, 240.0)
        assert solution.surface_temperature_k > (240.0 / SIGMA) ** 0.25
        assert fetch.stand_ins_used()


def test_the_warming_is_the_model_in_the_repos_variable():
    warm = rc.warming(CARBON, 9.81, 240.0)
    solution = rc.solve(CARBON, 9.81, 240.0)
    assert warm.converged and warm.solution.surface_temperature_k == pytest.approx(solution.surface_temperature_k)
    assert warm.surface_temperature_k == warm.solution.surface_temperature_k
    assert warm.effective_temperature_k == pytest.approx((240.0 / SIGMA) ** 0.25)
    assert warm.factor == pytest.approx((warm.surface_temperature_k / warm.effective_temperature_k) ** 4 - 1.0)
    assert warm.factor > 0.0 and warm.tropopause_pressure_pa == warm.solution.tropopause_pressure_pa
    assert warm.tropopause_temperature_k < warm.surface_temperature_k


def test_a_haze_that_takes_sunlight_high_up_lowers_the_warming_factor():
    clear = rc.warming(CARBON, 9.81, 240.0)
    hazy = rc.warming(CARBON, 9.81, 240.0, aloft_share=0.5, aloft_pressure_pa=2.0e3)
    assert hazy.factor < clear.factor and hazy.converged


def test_the_factor_never_goes_negative_even_when_the_model_does_not_converge():
    warm = rc.warming({"N2": 5.0e4, "CO2": 1.0e4}, 9.81, 240.0, aloft_share=0.999, aloft_pressure_pa=2.0e3)
    assert warm.factor >= 0.0 and math.isfinite(warm.surface_temperature_k)


# ---- the heat capacity of every gas of the table (GEN.852, row IR9) ----

def test_a_gas_the_database_writes_under_another_name_takes_its_nasa_fit():
    from planetgen.physics import species_thermo
    ocs = species_thermo.species("COS")
    low, high = ocs.range_k
    at = int(np.searchsorted(rc._CP_GRID, 300.0))
    assert rc._cp_over_r("OCS")[at] == pytest.approx(ocs.cp_over_r(min(max(rc._CP_GRID[at], low), high)))
    assert rc._cp_over_r("C2H2")[at] == pytest.approx(species_thermo.species("C2H2,acetylene").cp_over_r(rc._CP_GRID[at]))


def test_a_gas_without_a_fit_takes_the_table_heat_capacities_held_at_the_ends():
    cp = rc._cp_over_r("HCl")
    at = int(np.searchsorted(rc._CP_GRID, 300.0))
    assert rc._CP_GRID[at] == pytest.approx(300.0)
    assert cp[at] * rc.R == pytest.approx(29.14, rel=1e-3)
    assert cp[-1] * rc.R == pytest.approx(31.63, rel=1e-3)
    assert cp[0] * rc.R == pytest.approx(29.12, rel=1e-3)


def test_a_gas_with_no_heat_capacity_at_all_counts_its_atoms():
    assert rc._cp_over_r("HC3N")[0] == pytest.approx(4.0 + 0.5 * (3 * 5 - 6))
    assert rc._cp_over_r("NH4Cl")[0] == pytest.approx(4.0 + 0.5 * (3 * 6 - 6))


def test_every_gas_of_the_table_has_a_heat_capacity_on_the_grid():
    import json
    from planetgen.physics import gases
    with open(gases.DATA_FILE, encoding="utf-8") as handle:
        rows = json.load(handle)
    rows = rows["gases"] if isinstance(rows, dict) else rows
    assert len(rows) > 100
    for row in rows:
        cp = rc._cp_over_r(row["formula"])
        assert cp.shape == rc._CP_GRID.shape and (cp > 0).all(), row["formula"]


def test_an_air_with_carbonyl_sulphide_makes_a_column():
    column = rc.Column({"N2": 0.8e5, "CO2": 0.1e5, "OCS": 50.0}, 9.8)
    assert column.cp_dry[0] > 0

