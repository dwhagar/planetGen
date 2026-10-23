# tests/test_saturation.py

"""GEN.470: the saturation-pressure function, against the gas table's triple points and the checked inversions of
docs/design/atmosphere-composition-parts/06-condensates-clouds-hazes-dust-minerals.md section 2.4."""

import json
import math

import pytest

from planetgen.physics import gases, saturation

# The species whose triple-point pressure on the highest-temperature fit is further than 3 percent from the gas table's,
# with the ratio (the module docstring says the same).
OFF_BY_OVER_THREE_PERCENT = {"SO2": 1.038, "C6H6": 1.071, "C2H2": 0.964}
NO_TRIPLE_PRESSURE = {"HCN"}


def _triple_ratio(formula):
    row = gases.gas(formula)
    return saturation.pressure_pa(formula, row.triple_point_k) / row.triple_point_pa


def _invert(formula, pressure_pa, low=6.0, high=400.0):
    """The temperature at which `formula` saturates at `pressure_pa`, by bisection on the monotone part."""
    for _ in range(200):
        middle = 0.5 * (low + high)
        if saturation.pressure_pa(formula, middle) < pressure_pa:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def test_every_species_is_in_the_gas_table_with_a_triple_point():
    assert saturation.formulas()
    for formula in saturation.formulas():
        assert gases.has(formula)
        assert saturation.triple_point_k(formula) is not None


@pytest.mark.parametrize("formula", sorted(set(saturation.formulas()) - {"H2O"} - NO_TRIPLE_PRESSURE))
def test_the_triple_point_pressure_is_the_gas_tables(formula):
    ratio = _triple_ratio(formula)
    if formula in OFF_BY_OVER_THREE_PERCENT:
        assert ratio == pytest.approx(OFF_BY_OVER_THREE_PERCENT[formula], abs=0.002)
    else:
        assert abs(ratio - 1.0) < 0.03


def test_carbon_dioxide_is_within_a_percent_on_the_upper_fit_and_four_under_on_the_lower():
    assert _triple_ratio("CO2") == pytest.approx(1.005, abs=0.002)
    lower = saturation.species("CO2").curves[0]
    row = gases.gas("CO2")
    assert saturation._fray_schmitt_pa(lower.coefficients, row.triple_point_k) / row.triple_point_pa == pytest.approx(0.961, abs=0.002)


def test_hydrogen_cyanide_has_no_triple_point_pressure_to_check_against():
    assert gases.gas("HCN").triple_point_pa is None
    assert saturation.pressure_pa("HCN", 259.86) == pytest.approx(1.934e4, rel=0.001)


def test_water_ice_is_the_iapws_form():
    assert saturation.pressure_pa("H2O", 273.16) == pytest.approx(611.657, rel=1e-9)
    assert saturation.pressure_pa("H2O", 250.0) == pytest.approx(76.03, rel=1e-3)
    assert saturation.pressure_pa("H2O", 200.0) == pytest.approx(0.1626, rel=1e-3)
    assert saturation.saturation("H2O", 250.0).curve_id == "H2O-FW"


# The temperatures at which each ice saturates at 1e-9, 1e-5, 1e-3 and 0.1 bar, from the table of part 06 section 2.4
# (the inversions of the same polynomials, so this tests the code and the transcription of the coefficients together;
# the table used the beta fit of N2 and CO throughout, here the alpha fit serves below the switch, so the coldest N2 and CO
# entries are the alpha fit's: 26.8 K against 26.3 K, 28.8 K against 29.1 K).
CHECKED_INVERSIONS = {
    "Ne": (8.4, 12.2, 15.6, 21.7), "N2": (26.8, 36.7, 46.0, 62.2), "CO": (28.8, 40.0, 49.6, 65.8),
    "Ar": (30.4, 42.6, 53.4, 71.7), "O2": (31.6, 42.9, 53.24, None), "CH4": (37.4, 52.7, 66.3, 89.5),
    "Kr": (42.3, 59.1, 73.8, 98.7), "Xe": (58.9, 82.1, 102.4, 136.6), "C2H6": (67.8, 89.9, None, None),
    "C2H4": (64.6, 85.9, 103.1, None), "H2S": (80.3, 110.5, 136.1, 177.1), "CO2": (86.1, 114.0, 136.3, 170.1),
    "NH3": (101.5, 134.6, 161.0, None), "SO2": (111.9, 146.7, 174.3, None), "HCN": (125.2, 165.8, 198.9, 250.1),
    "C6H6": (148.7, 196.2, 233.5, None), "H2O": (161.1, 212.6, 252.8, None),
}


@pytest.mark.parametrize("formula", sorted(CHECKED_INVERSIONS))
def test_the_checked_inversions_are_reproduced(formula):
    for pressure_bar, expected in zip((1e-9, 1e-5, 1e-3, 0.1), CHECKED_INVERSIONS[formula]):
        if expected is not None:
            assert _invert(formula, pressure_bar * 1e5) == pytest.approx(expected, abs=0.15), (formula, pressure_bar)


def test_two_fit_species_switch_where_the_fits_cross_or_come_closest():
    for formula in ("N2", "CO", "O2", "CH3OH", "Ne", "CO2"):
        row = saturation.species(formula)
        assert len(row.curves) == (3 if formula == "O2" else 2)
        for low, high in zip(row.curves, row.curves[1:]):
            switch = high.from_k
            jump = abs(math.log(saturation._fray_schmitt_pa(low.coefficients, switch)
                                / saturation._fray_schmitt_pa(high.coefficients, switch)))
            assert jump < (0.02 if row.selection == "closest-approach" else 0.005), formula
            assert saturation.saturation(formula, switch - 0.01).curve_id == low.curve_id
            assert saturation.saturation(formula, switch).curve_id == high.curve_id


def test_the_switch_points_are_near_the_known_solid_transitions():
    transitions = {"N2": (35.61,), "CO": (61.55,), "O2": (23.88, 43.8), "CH3OH": (157.4,)}
    for formula, expected in transitions.items():
        found = [curve.from_k for curve in saturation.species(formula).curves[1:]]
        assert len(found) == len(expected)
        for switch, transition in zip(found, expected):
            assert abs(switch - transition) < 1.0, formula


def test_each_curve_rises_with_temperature_from_one_kelvin_to_the_triple_point():
    for formula in saturation.formulas():
        top = saturation.triple_point_k(formula)
        previous = 0.0
        t = 1.0
        while t <= top:
            value = saturation.pressure_pa(formula, t)
            assert math.isfinite(value) and value >= previous, (formula, t)
            previous = value
            t += 0.5


def test_each_floor_is_where_the_polynomial_turns_over():
    for formula in saturation.formulas():
        for curve in saturation.species(formula).curves:
            if curve.floor_k is None or curve.from_k > 0 or curve.floor_k < 2.0:
                continue
            steps = [curve.floor_k - 1.0 + 0.05 * n for n in range(21)]
            values = [saturation._fray_schmitt_pa(curve.coefficients, t) for t in steps]
            fell = [steps[n] for n in range(1, 21) if values[n] < values[n - 1]]
            assert fell, formula
            assert max(fell) < curve.floor_k + 0.01, formula
            assert max(fell) > curve.floor_k - 0.2, formula


def test_status_follows_the_temperature():
    assert saturation.saturation("N2", 40.0).status == saturation.COLD_UNVERIFIED
    assert saturation.saturation("N2", 55.0).status == saturation.OK
    assert saturation.saturation("N2", 70.0).status == saturation.OK
    assert saturation.saturation("N2", 70.0).phase == saturation.LIQUID
    assert saturation.saturation("N2", 5.0).status == saturation.BELOW_FIT
    assert saturation.saturation("CO2", 150.0).status == saturation.OK
    assert saturation.saturation("CO2", 230.0).phase == saturation.LIQUID
    assert saturation.saturation("C2H2", 200.0).status == saturation.ABOVE_TRIPLE
    assert saturation.saturation("C2H2", 200.0).phase == saturation.SOLID
    assert saturation.saturation("H2O", 273.16).status == saturation.OK
    assert saturation.saturation("H2O", 300.0).status == saturation.OK
    assert saturation.saturation("H2O", 40.0).status == saturation.COLD_UNVERIFIED
    assert saturation.saturation("N2", 40.0).status in saturation.STATUSES


def test_the_cold_limit_is_fifty_kelvin():
    assert saturation.COLD_UNVERIFIED_BELOW_K == 50.0
    assert saturation.saturation("CH4", 49.9).status == saturation.COLD_UNVERIFIED
    assert saturation.saturation("CH4", 50.0).status == saturation.OK


def test_an_unknown_species_and_a_bad_temperature_are_refused():
    with pytest.raises(saturation.UnknownSpecies):
        saturation.saturation("He", 3.0)
    assert not saturation.has("He")
    for bad in (0.0, -3.0):
        with pytest.raises(ValueError):
            saturation.saturation("N2", bad)


def test_the_data_file_labels_its_stand_ins_and_the_unread_paper():
    with open(saturation.DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    assert "STAND-IN" in data["basis"]
    assert "2ND-HAND" in data["basis"]
    assert data["cold_unverified_below_k"] == 50.0
    assert "H2S-2" in data["species"]["H2S"]["note"]


def test_pressure_is_in_pascals_not_bar():
    assert saturation.pressure_pa("CO2", 194.7) == pytest.approx(1.01e5, rel=0.03)


LIQUID_FORMULAS = sorted(f for f in saturation.formulas() if saturation.species(f).liquid is not None)

# solid pressure over liquid pressure at the triple point, the "Err at Tt" of part 01 section 7.1 (the pieces do not meet)
JOINS = {"N2": 0.9726, "CH4": 1.0165, "H2S": 1.0165, "NH3": 1.0190, "SO2": 1.0381, "HCN": 1.0333, "C6H6": 1.0707}


def test_every_species_but_acetylene_has_a_liquid_branch():
    assert set(saturation.formulas()) - set(LIQUID_FORMULAS) == {"C2H2"}
    assert saturation.species("HCN").liquid.form == "dippr101"


@pytest.mark.parametrize("formula", LIQUID_FORMULAS)
def test_the_liquid_boils_at_the_gas_tables_boiling_point(formula):
    row = gases.gas(formula)
    if row.boiling_point_is_sublimation:
        pytest.skip("the boiling point of the gas table is a sublimation point")
    assert saturation.pressure_pa(formula, row.boiling_point_k) == pytest.approx(101325.0, rel=0.01)
    assert saturation.saturation(formula, row.boiling_point_k).phase == saturation.LIQUID


@pytest.mark.parametrize("formula", [f for f in LIQUID_FORMULAS if f != "HCN"])
def test_the_critical_point_is_the_gas_tables(formula):
    liquid = saturation.species(formula).liquid
    row = gases.gas(formula)
    assert liquid.t_max_k == pytest.approx(row.critical_temperature_k, abs=0.02)
    assert saturation._liquid_pa(liquid, liquid.t_max_k * (1 - 1e-12)) == pytest.approx(liquid.pc_pa, rel=1e-6)
    assert liquid.pc_pa == pytest.approx(row.critical_pressure_pa, rel=0.005)


def test_water_liquid_is_the_steam_table():
    assert saturation.pressure_pa("H2O", 300.0) == pytest.approx(3536.8, rel=0.003)
    assert saturation.pressure_pa("H2O", 373.124) == pytest.approx(101325.0, rel=0.003)


@pytest.mark.parametrize("formula", LIQUID_FORMULAS)
def test_the_two_branches_meet_within_the_published_spread(formula):
    tt = saturation.triple_point_k(formula)
    below = saturation.saturation(formula, tt)
    above = saturation.saturation(formula, tt + 1e-6)
    assert (below.phase, above.phase) == (saturation.SOLID, saturation.LIQUID)
    ratio = below.pressure_pa / above.pressure_pa
    assert ratio == pytest.approx(JOINS.get(formula, 1.0), abs=0.005 if formula in JOINS else 0.01)


@pytest.mark.parametrize("formula", LIQUID_FORMULAS)
def test_the_liquid_rises_to_the_critical_point_and_then_there_is_no_saturation(formula):
    liquid = saturation.species(formula).liquid
    tt = saturation.triple_point_k(formula)
    previous = 0.0
    t = tt + 0.5
    while t < liquid.t_max_k:
        value = saturation.pressure_pa(formula, t)
        assert value > previous, (formula, t)
        previous = value
        t += (liquid.t_max_k - tt) / 100.0
    hot = saturation.saturation(formula, liquid.t_max_k + 1.0)
    assert hot.status == saturation.ABOVE_CRITICAL
    assert math.isinf(hot.pressure_pa)


def test_hydrogen_is_a_clausius_clapeyron_solid_below_its_triple_point():
    row = gases.gas("H2")
    assert saturation.pressure_pa("H2", row.triple_point_k) == pytest.approx(row.triple_point_pa, rel=1e-9)
    assert saturation.saturation("H2", 10.0).status == saturation.COLD_UNVERIFIED
    assert saturation.pressure_pa("H2", 10.0) < 0.2 * row.triple_point_pa
    assert saturation.saturation("H2", 10.0).curve_id == "H2-CC"


def test_solid_carbon_dioxide_is_as_cold_as_its_latent_heat_allows():
    # 1e-9 bar from 194.7 K at 1 bar with the sublimation heat of 25.2 to 27.2 kJ/mol (gas table): 83.5 to 87.2 K
    low, high = gases.gas("CO2").latent_heat_sublimation_kj_mol
    for heat, expected in ((low, 83.5), (high, 87.2)):
        slope = heat * 1000.0 / saturation.GAS_CONSTANT_J_MOL_K
        t = 1.0 / (1.0 / 194.7 + math.log(1e9) / slope)
        assert t == pytest.approx(expected, abs=0.1)
    assert 83.0 < _invert("CO2", 1e-9 * 1e5) < 88.0


# Part 01 section 7.4 (T8): the condensation temperature, K, at 100 kPa, 10 kPa, 1 kPa, 100 Pa, 10 Pa, 1 Pa and 1 mPa
# (COMPUTED by the Atmosphere research lane from the same curves; carbon dioxide is not here, see the next test).
CHECKED_CONDENSATION_K = {
    "H2": (20.3, 14.5, 11.4, 9.4, 8.0, 7.0, 5.0), "N2": (77.2, 62.2, 52.8, 46.0, 40.8, 36.7, 28.7),
    "CO": (81.5, 65.9, 56.6, 49.6, 44.3, 40.0, 31.0), "Ar": (87.2, 71.7, 61.2, 53.4, 47.4, 42.6, 32.7),
    "O2": (90.1, 72.7, 61.3, 53.2, 47.4, 42.9, 33.8), "CH4": (111.5, 89.5, 76.1, 66.3, 58.7, 52.7, 40.3),
    "SO2": (262.8, 221.0, 192.7, 174.3, 159.3, 146.7, 118.9), "H2O": (372.8, 318.9, 280.1, 252.8, 231.0, 212.6, 171.5),
    "H2S": (212.6, 177.1, 153.9, 136.1, 122.0, 110.5, 86.2), "NH3": (239.6, 202.0, 178.6, 161.0, 146.6, 134.6, 108.1),
    "C2H6": (184.3, 150.4, 127.9, 111.8, 99.5, 89.9, 72.2), "HCN": (298.5, 250.1, 221.4, 198.9, 180.8, 165.8, 133.2),
}
CHECKED_PRESSURES_PA = (1e5, 1e4, 1e3, 1e2, 10.0, 1.0, 1e-3)


@pytest.mark.parametrize("formula", sorted(CHECKED_CONDENSATION_K))
def test_the_condensation_temperatures_are_the_checked_table(formula):
    for pressure, expected in zip(CHECKED_PRESSURES_PA, CHECKED_CONDENSATION_K[formula]):
        assert saturation.condensation_temperature_k(formula, pressure) == pytest.approx(expected, abs=0.12), (formula, pressure)


def test_carbon_dioxide_condenses_where_its_latent_heat_puts_it():
    # part 01's table used CO2-2 (3.9 percent under to 0.5 percent over at the triple point decided it), which is 4 K cold
    # at 100 Pa and 17 K at 1 mPa; the sublimation heat of the gas table (25.2 to 27.2 kJ/mol) from 194.7 K at 1 bar
    # brackets the answer instead
    low, high = gases.gas("CO2").latent_heat_sublimation_kj_mol
    for pressure in CHECKED_PRESSURES_PA[1:]:
        ln_ratio = math.log(1e5 / pressure)
        cold = 1.0 / (1.0 / 194.7 + ln_ratio * saturation.GAS_CONSTANT_J_MOL_K / (low * 1000.0))
        warm = 1.0 / (1.0 / 194.7 + ln_ratio * saturation.GAS_CONSTANT_J_MOL_K / (high * 1000.0))
        assert cold - 1.5 < saturation.condensation_temperature_k("CO2", pressure) < warm + 1.5, pressure
    assert saturation.condensation_temperature_k("CO2", 1e5) == pytest.approx(194.5, abs=0.5)


def test_the_condensation_temperature_is_the_inverse_and_has_limits():
    for formula in ("N2", "H2O", "CH4", "C2H2"):
        t = saturation.condensation_temperature_k(formula, 500.0)
        assert saturation.pressure_pa(formula, t) == pytest.approx(500.0, rel=1e-9)
    assert saturation.condensation_temperature_k("N2", 1e8) is None
    with pytest.raises(ValueError):
        saturation.condensation_temperature_k("N2", 0.0)
