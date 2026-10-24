# tests/test_equilibrium.py

"""
The gas-phase equilibrium (GEN.518): the free energy of a mixture is minimised under the conservation of each element, and the
edges it gives are those of the independent computation recorded in tunable-gaps-sources/group4-giants-bd-clouds.md section 2.4.
"""

import pytest

from planetgen.physics import equilibrium as eq
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.generation.possibility_table import Subject

SOLAR = eq.element_budget(1.0, 0.55)


def test_the_table_is_the_nasa_polynomials_with_a_known_formation_free_energy():
    assert len(eq.SPECIES) >= 50 and "4513" in open(eq.DATA_FILE, encoding="utf-8").read()
    # H2O -> H2 + 1/2 O2 at 1500 K: the standard free energy of formation of water is -13.2 RT there (the table's value)
    delta = eq.SPECIES["H2O"].gibbs_rt(1500.0) - eq.SPECIES["H2"].gibbs_rt(1500.0) - 0.5 * eq.SPECIES["O2"].gibbs_rt(1500.0)
    assert delta == pytest.approx(-13.2, abs=0.3)


def test_elements_are_conserved_and_the_fractions_add_to_one():
    solution = eq.solve(1200.0, 1.0, SOLAR)
    assert sum(solution.moles.values()) == pytest.approx(solution.total)
    for symbol, atoms in SOLAR.items():
        held = sum(n * eq.SPECIES[s].composition.get(symbol, 0) for s, n in solution.moles.items())
        assert held == pytest.approx(atoms, rel=1e-8), symbol
    assert sum(solution.fraction(s) for s in solution.moles) == pytest.approx(1.0)


def test_methane_holds_the_carbon_when_cool_and_carbon_monoxide_when_hot():
    cool, hot = eq.solve(800.0, 1.0, SOLAR), eq.solve(1800.0, 1.0, SOLAR)
    assert cool.fraction("CH4") > 100 * cool.fraction("CO") and hot.fraction("CO") > 1000 * hot.fraction("CH4")
    assert cool.fraction("NH3") > cool.fraction("N2") * 0.2 and hot.fraction("N2") > 100 * hot.fraction("NH3")


def test_the_pressure_favours_the_molecules():
    low, high = eq.solve(2500.0, 0.001, SOLAR), eq.solve(2500.0, 1.0, SOLAR)
    assert low.fraction("H") > high.fraction("H") and high.fraction("H2O") > low.fraction("H2O")


@pytest.mark.parametrize("pressure, metallicity, expected", [
    (0.1, 1.0, 966), (1.0, 1.0, 1156), (10.0, 1.0, 1438), (0.1, 10.0, 892), (1.0, 10.0, 1052)])
def test_the_methane_carbon_monoxide_line_matches_the_independent_computation(pressure, metallicity, expected):
    assert eq.ch4_co_temperature(pressure, eq.element_budget(metallicity, 0.55)) == pytest.approx(expected, abs=6)


@pytest.mark.parametrize("pressure, ninety, half, atomic", [
    (0.001, 2331, 2567, 2574), (0.01, 2542, 2843, 2884), (0.1, 2792, 3182, 3278), (1.0, 3095, 3609, 3794)])
def test_the_water_and_hydrogen_dissociation_temperatures_match_the_independent_computation(pressure, ninety, half, atomic):
    assert eq.h2o_survival_temperature(pressure, SOLAR, 0.9) == pytest.approx(ninety, rel=0.01)
    assert eq.h2o_survival_temperature(pressure, SOLAR, 0.5) == pytest.approx(half, rel=0.01)
    assert eq.h2_dissociation_temperature(pressure, SOLAR) == pytest.approx(atomic, rel=0.01)


def test_the_hot_jupiter_regime_edges_of_the_recipes_are_the_computed_ones():
    s = Subject("J", {"zone": "h", "tday_k": 0.0, "photosphere_bar": 1.0})
    # the soft window is half at its edge: row h1 has 0.5 at the methane line, row h3 has 0.5 at the 10 percent water line
    h1 = RECIPES.entry("h1")
    for bar in (0.01, 0.1, 1.0):
        edge = eq.ch4_co_temperature(bar, SOLAR)
        probability = dict(RECIPES.probabilities("h", Subject("J", {"zone": "h", "tday_k": edge, "photosphere_bar": bar})))
        assert probability["h1"] == pytest.approx(0.5, abs=0.03), (bar, edge)
        edge = eq.h2o_survival_temperature(bar, SOLAR, 0.9)
        probability = dict(RECIPES.probabilities("h", Subject("J", {"zone": "h", "tday_k": edge, "photosphere_bar": bar})))
        assert probability["h3"] == pytest.approx(0.5, abs=0.03), (bar, edge)
    assert h1.state == "partly verified" and s.quantities["photosphere_bar"] == 1.0


def test_a_budget_without_a_species_for_an_element_or_a_carbon_rich_water_edge_is_refused():
    with pytest.raises(eq.EquilibriumError, match="no species carries"):
        eq.solve(1000.0, 1.0, {"H": 1.0, "Ba": 1e-6})
    with pytest.raises(eq.EquilibriumError, match="no oxygen"):
        eq.h2o_survival_temperature(1.0, eq.element_budget(1.0, 1.2))
    with pytest.raises(eq.EquilibriumError, match="unknown species"):
        eq.solve(1000.0, 1.0, SOLAR, species=["H2", "XeF9"])


def test_the_solution_does_not_depend_on_the_order_of_the_species():
    names = list(eq.DEFAULT_SPECIES)
    forward = eq.solve(1500.0, 0.1, SOLAR, species=names)
    backward = eq.solve(1500.0, 0.1, SOLAR, species=names[::-1])
    assert forward.fraction("CO") == pytest.approx(backward.fraction("CO"), rel=1e-8)
