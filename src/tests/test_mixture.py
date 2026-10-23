"""
The mixture accessor (GEN.411's first slice, `planetgen.physics.mixture`).
"""

from types import SimpleNamespace

import pytest

from planetgen.physics import atmosphere, elements, gases, mixture


def _air(**kpa):
    values = {f"p_{gas}_kpa": 0.0 for gas in atmosphere.SPECIES}
    values.update({f"p_{gas}_kpa": value for gas, value in kpa.items()})
    return SimpleNamespace(**values)


def test_every_stored_species_is_a_row_of_the_gas_table():
    for key, formula in mixture.STORED_FORMULAS.items():
        assert gases.has(formula), key
        assert gases.gas(formula).formula.lower() == key


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
