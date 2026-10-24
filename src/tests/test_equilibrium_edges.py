# tests/test_equilibrium_edges.py

"""
The regime edges of a hot atmosphere (GEN.518, GEN.672): the temperatures the one equilibrium solver (`gas_equilibrium`) gives for the
methane-carbon monoxide line, the water and the hydrogen dissociation are those of the independent computation recorded in
tunable-gaps-sources/group4-giants-bd-clouds.md section 2.4, and the recipe rows h1 and h3 are half open at them.
"""

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import equilibrium_edges as edges
from planetgen.physics import gas_equilibrium, species_thermo
from planetgen.physics.atmosphere_recipes import RECIPES

SOLAR = edges.reference_budget(1.0, 0.55)
TABLE_SLACK = 0.07
# The points of the recipe tables (rows h1 to h3) were computed with the first GEN.518 module; the one solver gives the independent
# record's values, up to 16 K (0.6 percent) lower for the water edge, which moves the probability at the edge by up to 0.06 of the 100 K window.


def test_water_forms_from_its_elements_with_the_known_formation_free_energy():
    # H2O -> H2 + 1/2 O2 at 1500 K: the standard free energy of formation of water is -13.2 RT there
    delta = (species_thermo.g_over_rt("H2O", 1500.0) - species_thermo.g_over_rt("H2", 1500.0)
             - 0.5 * species_thermo.g_over_rt("O2", 1500.0))
    assert delta == pytest.approx(-13.2, abs=0.3)


def test_the_reference_budget_scales_the_metals_and_sets_the_carbon_to_oxygen():
    solar = edges.reference_budget(1.0)
    tenfold = edges.reference_budget(10.0, 0.55)
    assert tenfold["He"] == solar["He"] and tenfold["H"] == 1.0
    assert tenfold["O"] == pytest.approx(10.0 * solar["O"])
    assert tenfold["C"] == pytest.approx(0.55 * tenfold["O"])
    assert set(solar) == {"H", "He", "C", "N", "O", "S", "Ar"}


def test_methane_holds_the_carbon_when_cool_and_carbon_monoxide_when_hot():
    cool, hot = gas_equilibrium.solve(SOLAR, 800.0, 1.0), gas_equilibrium.solve(SOLAR, 1800.0, 1.0)
    assert cool.mole_fractions["CH4"] > 100 * cool.mole_fractions["CO"]
    assert hot.mole_fractions["CO"] > 1000 * hot.mole_fractions["CH4"]
    assert cool.mole_fractions["NH3"] > cool.mole_fractions["N2"] * 0.2
    assert hot.mole_fractions["N2"] > 100 * hot.mole_fractions["NH3"]


def test_the_pressure_favours_the_molecules():
    low, high = gas_equilibrium.solve(SOLAR, 2500.0, 0.001), gas_equilibrium.solve(SOLAR, 2500.0, 1.0)
    assert low.mole_fractions["H"] > high.mole_fractions["H"] and high.mole_fractions["H2O"] > low.mole_fractions["H2O"]


@pytest.mark.parametrize("pressure, metallicity, expected", [
    (0.1, 1.0, 966), (1.0, 1.0, 1156), (10.0, 1.0, 1438), (0.1, 10.0, 892), (1.0, 10.0, 1052)])
def test_the_methane_carbon_monoxide_line_matches_the_independent_computation(pressure, metallicity, expected):
    assert edges.ch4_co_temperature(pressure, edges.reference_budget(metallicity, 0.55)) == pytest.approx(expected, abs=6)


@pytest.mark.parametrize("pressure, ninety, half, atomic", [
    (0.001, 2331, 2567, 2574), (0.01, 2542, 2843, 2884), (0.1, 2792, 3182, 3278), (1.0, 3095, 3609, 3794)])
def test_the_water_and_hydrogen_dissociation_temperatures_match_the_independent_computation(pressure, ninety, half, atomic):
    assert edges.h2o_survival_temperature(pressure, SOLAR, 0.9) == pytest.approx(ninety, rel=0.01)
    assert edges.h2o_survival_temperature(pressure, SOLAR, 0.5) == pytest.approx(half, rel=0.01)
    assert edges.h2_dissociation_temperature(pressure, SOLAR) == pytest.approx(atomic, rel=0.01)


def test_the_hot_jupiter_regime_edges_of_the_recipes_are_the_computed_ones():
    # the soft window is half at its edge: row h1 has 0.5 at the methane line, row h3 has 0.5 at the 10 percent water line
    h1 = RECIPES.entry("h1")
    for bar in (0.01, 0.1, 1.0):
        edge = edges.ch4_co_temperature(bar, SOLAR)
        probability = dict(RECIPES.probabilities("h", Subject("J", {"zone": "h", "tday_k": edge, "photosphere_bar": bar})))
        assert probability["h1"] == pytest.approx(0.5, abs=TABLE_SLACK), (bar, edge)
        edge = edges.h2o_survival_temperature(bar, SOLAR, 0.9)
        probability = dict(RECIPES.probabilities("h", Subject("J", {"zone": "h", "tday_k": edge, "photosphere_bar": bar})))
        assert probability["h3"] == pytest.approx(0.5, abs=TABLE_SLACK), (bar, edge)
    assert h1.state == "partly verified"


def test_a_carbon_rich_water_edge_or_a_bad_budget_is_refused():
    with pytest.raises(gas_equilibrium.EquilibriumError, match="no oxygen"):
        edges.h2o_survival_temperature(1.0, edges.reference_budget(1.0, 1.2))
    with pytest.raises(gas_equilibrium.EquilibriumError, match="not positive"):
        edges.reference_budget(0.0)
    with pytest.raises(gas_equilibrium.EquilibriumError, match="not positive"):
        edges.reference_budget(1.0, -1.0)


def test_the_solution_does_not_depend_on_the_order_of_the_species():
    names = ["H2", "H2O", "CO", "CH4", "CO2", "N2", "NH3", "H", "O2", "He", "Ar", "OH", "HCN", "H2S", "S2", "SO2", "O"]
    forward = gas_equilibrium.solve(SOLAR, 1500.0, 0.1, species=names)
    backward = gas_equilibrium.solve(SOLAR, 1500.0, 0.1, species=names[::-1])
    assert forward.mole_fractions["CO"] == pytest.approx(backward.mole_fractions["CO"], rel=1e-8)
