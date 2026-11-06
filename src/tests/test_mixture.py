"""
The mixture accessor (GEN.411's first slice, `planetgen.physics.mixture`).
"""

from types import SimpleNamespace

import numpy as np
import pytest

from planetgen.physics import atmosphere, constants, elements, gases, mixture, species_thermo

from tests.air_support import air


def _air(**kpa):
    return SimpleNamespace(**air(**{f"p_{gas}_kpa": value for gas, value in kpa.items()}))


def test_every_stored_species_is_a_row_of_the_gas_table():
    for key in atmosphere.SPECIES:
        assert gases.has(gases.stored(key).formula), key
        assert gases.stored(key).formula.lower() == key


def test_the_partials_come_back_in_pascal_under_the_table_formulas():
    mix = mixture.partial_pressures_pa(_air(n2=78.084, o2=20.946, ar=0.934, co2=0.042))
    assert mix == {"N2": pytest.approx(78084.0), "O2": pytest.approx(20946.0), "Ar": pytest.approx(934.0),
                   "CO2": pytest.approx(42.0)}
    assert mixture.total_pressure_pa(mix) == pytest.approx(100006.0)


def test_an_airless_body_and_an_old_body_have_an_empty_mixture():
    assert mixture.partial_pressures_pa(_air()) == {}
    assert mixture.partial_pressures_pa(SimpleNamespace()) == {}
    with pytest.raises(ValueError, match="empty"):
        mixture.mole_fractions({})


def test_the_mean_molar_mass_is_the_mole_weighted_sum():
    mix = {"N2": 78084.0, "O2": 20946.0, "Ar": 934.0}
    fractions = mixture.mole_fractions(mix)
    assert sum(fractions.values()) == pytest.approx(1.0)
    expected = sum(fraction * elements.molar_mass(formula) for formula, fraction in fractions.items())
    assert mixture.mean_molar_mass_kg_mol(mix) == pytest.approx(expected)
    assert mixture.mean_molar_mass_kg_mol(mix) == pytest.approx(0.0289, abs=1e-4)
    assert mixture.mean_molar_mass_kg_mol({"CO2": 1.0}) == pytest.approx(0.044009, abs=1e-5)


def test_any_gas_of_the_table_may_be_in_a_mixture():
    # helium and ammonia are balance gases today; the accessor already takes them
    assert mixture.mean_molar_mass_kg_mol({"H2": 86.0, "He": 14.0}) == pytest.approx(
        0.86 * gases.gas("H2").molar_mass + 0.14 * gases.gas("He").molar_mass)
    with pytest.raises(gases.UnknownGas):
        mixture.mean_molar_mass_kg_mol({"unobtainium": 1.0})


def test_a_trace_gas_with_no_nasa_fit_has_a_heat_capacity_ratio_from_the_gas_table():
    """TEST.604: CH3Cl has a gas-table row and no NASA Glenn fit; the ratio of an air that holds it used to raise KeyError."""
    assert "CH3Cl" not in species_thermo.BY_FORMULA
    air = {"N2": 0.78e5, "O2": 0.21e5, "CH3Cl": 1.0e3}
    assert 1.2 < mixture.heat_capacity_ratio(air, 288.0) < 1.5
    pure = mixture.heat_capacity_ratio({"CH3Cl": 1.0e5}, 288.0)
    row = gases.gas("CH3Cl")
    cp = np.interp(288.0, sorted(row.cp_j_mol_k), [row.cp_j_mol_k[t] for t in sorted(row.cp_j_mol_k)]) / constants.R
    assert pure == pytest.approx(cp / (cp - 1.0))
    assert mixture.heat_capacity_ratio({"N2": 1.0e5}, 288.0) == pytest.approx(mixture.heat_capacity_ratio({"N2": 1.0e5, "CH3Cl": 0.0}, 288.0))


def test_a_gas_with_neither_a_fit_nor_table_heat_capacities_takes_the_rigid_rotor_count(monkeypatch):
    class Row:
        cp_j_mol_k = {100: None, 300: None, 1000: None}
    monkeypatch.setattr(mixture.gases, "gas", lambda formula: Row())
    assert mixture.cp_over_r_without_fit("CH3Cl", 300.0) == pytest.approx(4.0 + 0.5 * (3 * 5 - 6))
    assert mixture.cp_over_r_without_fit("HCl", 300.0) == pytest.approx(3.5)
    assert list(mixture.cp_over_r_without_fit("He", np.array([100.0, 300.0]))) == [2.5, 2.5]
