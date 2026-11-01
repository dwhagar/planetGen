"""
The air column and the Rayleigh optical depth of an atmosphere (GEN.417, `planetgen.physics.air_column` and
`planetgen.physics.rayleigh`): the numbers the research checks them against, and the way a body's stored mixture becomes a
column.
"""

import types

import pytest

from planetgen.physics import air_column, atmosphere, gases, rayleigh

from tests.air_support import fold_air

EARTH_G = 9.80665
EARTH_DRY_AIR_PA = {"N2": 0.78084 * 101325.0, "O2": 0.20946 * 101325.0, "Ar": 0.00934 * 101325.0, "CO2": 4e-4 * 101325.0}


def _earth_air():
    return air_column.from_partials(EARTH_DRY_AIR_PA, EARTH_G)


def test_an_earth_air_column_is_1033_grams_of_air_with_the_molar_mass_of_dry_air():
    column = _earth_air()
    assert column.column_g_cm2 == pytest.approx(1033.27, rel=1e-4)
    assert column.mean_molar_mass_kg_mol == pytest.approx(0.028964, rel=1e-3)
    assert column.total_number_column_cm2() == pytest.approx(2.148e25, rel=2e-3)
    assert column.fraction("N2") == pytest.approx(0.78084, abs=1e-4)
    assert sum(column.number_column_cm2(formula) for formula in column.partials_pa) == pytest.approx(column.total_number_column_cm2())


def test_a_column_needs_pressure_and_gravity():
    with pytest.raises(air_column.EmptyColumn):
        air_column.from_partials({"N2": 0.0}, EARTH_G)
    with pytest.raises(air_column.EmptyColumn):
        air_column.from_partials({"N2": 1e5}, 0.0)
    with pytest.raises(gases.UnknownGas):
        air_column.from_partials({"XX9": 1e5}, EARTH_G)


def test_earth_air_scatters_with_a_depth_of_0_097_at_550_nm():
    # [R0575: Table 4] gives 0.09721 for the standard atmosphere; the cross-section table is 0.3 percent below it
    assert rayleigh.optical_depth(_earth_air(), 550.0) == pytest.approx(0.09721, rel=0.01)


def test_the_blue_scatters_twelve_times_more_than_the_green():
    column = _earth_air()
    assert rayleigh.optical_depth(column, 300.0) / rayleigh.optical_depth(column, 550.0) == pytest.approx(12.6, rel=0.04)
    assert rayleigh.optical_depth(column, 450.0) / rayleigh.optical_depth(column, 650.0) == pytest.approx(4.4, rel=0.05)


def test_thick_carbon_dioxide_air_has_the_depth_of_the_research_worked_example():
    venus_like = air_column.from_partials({"CO2": 0.965 * 9.2e6, "N2": 0.035 * 9.2e6}, 8.87)
    assert rayleigh.optical_depth(venus_like, 550.0) == pytest.approx(16.4, rel=0.02)


def test_a_giant_reference_level_of_hydrogen_and_helium():
    # part 03 7.1: 0.089 at 550 nm for 1 bar of 86/14 H2/He at 24.8 m/s2 (a different He share by mass or number moves it 4 percent)
    giant = air_column.from_partials({"H2": 0.86e5, "He": 0.14e5}, 24.8)
    assert rayleigh.optical_depth(giant, 550.0) == pytest.approx(0.089, rel=0.06)


def test_the_fit_goes_through_the_three_tabulated_cross_sections_of_every_gas():
    for gas in gases.GASES:
        if gas.rayleigh is None:
            continue
        for wavelength, sigma in gas.rayleigh.sigma_cm2.items():
            assert rayleigh.cross_section_cm2(gas.formula, wavelength) == pytest.approx(sigma, rel=1e-6), gas.formula


def test_the_dispersion_terms_follow_the_research_ordering():
    assert 0.002 < rayleigh.dispersion_um2("He") < 0.003
    assert 0.005 < rayleigh.dispersion_um2("N2") < 0.007
    assert 0.010 < rayleigh.dispersion_um2("NH3") < 0.013
    assert rayleigh.dispersion_um2("He") < rayleigh.dispersion_um2("N2") < rayleigh.dispersion_um2("NH3")


def test_cross_sections_fall_with_the_fourth_power_of_the_wavelength_and_a_gas_without_a_row_is_refused():
    ratio = rayleigh.cross_section_cm2("Ar", 400.0) / rayleigh.cross_section_cm2("Ar", 800.0)
    assert 16.0 < ratio < 18.5
    missing = next(gas.formula for gas in gases.GASES if gas.rayleigh is None)
    with pytest.raises(rayleigh.NoRayleigh):
        rayleigh.cross_section_cm2(missing, 550.0)


def _body(**fields):
    base = {"planet_class": "M", "gravity": 1.0, "atmospheric_pressure": 101325.0}
    base.update(fields)
    return types.SimpleNamespace(**fold_air(base))


def test_a_bodys_stored_partials_become_its_column():
    body = _body(p_n2_kpa=79.12, p_o2_kpa=21.22, p_ar_kpa=0.95, p_co2_kpa=0.04)
    column = air_column.column_of(body)
    assert column.pressure_pa == pytest.approx(101330.0)
    assert column.gravity_m_s2 == pytest.approx(EARTH_G)
    assert rayleigh.optical_depth(column, 550.0) == pytest.approx(0.0972, rel=0.01)


def test_the_stored_balance_gas_of_a_giant_template_is_helium():
    giant = _body(planet_class="J", gravity=2.53, atmospheric_pressure=1e5, gas_partials_pa={"H2": 86000.0, "CH4": 300.0, "He": 13700.0})
    mixture = air_column.partials_pa(giant)
    assert mixture["H2"] == pytest.approx(86000.0)
    assert mixture["He"] == pytest.approx(1e5 - 86300.0)
    assert air_column.column_of(giant).pressure_pa == pytest.approx(1e5)


def test_sodium_is_left_out_of_the_column_for_want_of_a_cross_section():
    hot = _body(planet_class="B", gas_partials_pa={"N2": 80000.0, "Na": 20000.0})
    assert air_column.partials_pa(hot) == {"N2": 80000.0}


def test_a_body_without_air_or_gravity_has_no_column():
    assert air_column.column_of(_body(atmospheric_pressure=0.0)) is None
    assert air_column.column_of(_body(gravity=0.0, p_n2_kpa=80.0)) is None


def test_a_mixture_with_a_trace_of_an_uncovered_gas_has_no_depth():
    uncovered = next(gas.formula for gas in gases.GASES if gas.rayleigh is None and gas.molar_mass < 0.1)
    air = air_column.from_partials({**EARTH_DRY_AIR_PA, uncovered: 0.01 * 101325.0}, EARTH_G)
    with pytest.raises(rayleigh.NoRayleigh):
        rayleigh.optical_depth(air, 550.0)
