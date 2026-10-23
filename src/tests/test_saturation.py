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
    assert saturation.saturation("N2", 70.0).status == saturation.ABOVE_TRIPLE
    assert saturation.saturation("N2", 5.0).status == saturation.BELOW_FIT
    assert saturation.saturation("CO2", 150.0).status == saturation.OK
    assert saturation.saturation("CO2", 230.0).status == saturation.ABOVE_TRIPLE
    assert saturation.saturation("H2O", 273.16).status == saturation.OK
    assert saturation.saturation("H2O", 300.0).status == saturation.ABOVE_TRIPLE
    assert saturation.saturation("H2O", 40.0).status == saturation.COLD_UNVERIFIED
    assert saturation.saturation("N2", 40.0).status in saturation.STATUSES


def test_the_cold_limit_is_fifty_kelvin():
    assert saturation.COLD_UNVERIFIED_BELOW_K == 50.0
    assert saturation.saturation("CH4", 49.9).status == saturation.COLD_UNVERIFIED
    assert saturation.saturation("CH4", 50.0).status == saturation.OK


def test_an_unknown_species_and_a_bad_temperature_are_refused():
    with pytest.raises(saturation.UnknownSpecies):
        saturation.saturation("H2", 20.0)
    assert not saturation.has("H2")
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
