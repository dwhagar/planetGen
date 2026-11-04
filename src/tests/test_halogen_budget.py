# tests/test_halogen_budget.py

"""GEN.524 part 2: the chlorine and fluorine of the outgassing (`mixture_rocky`), which replace the presence draws of rows m5 and m6."""

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import mixture_draw as md, mixture_rocky as mr, mixture_trace, saturation, surface_pressure
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw

from tests.test_atmosphere_mixture_evolved import BASE

CL_AMU, F_AMU = 35.45, 18.998


def _inventory(**changes):
    subject = Subject("terrestrial", {**BASE, **changes})
    constants = {}
    with draw.bound(3):
        recipe = RECIPES.draw_all(subject, categories=mr.ROCKY_CATEGORIES)
        return mr.rocky_inventory(subject, recipe, constants), constants


def test_the_rows_for_hydrogen_chloride_and_fluoride_are_retired():
    assert not {"HCl", "HF"} & {gas for row in mixture_trace.ABUNDANCES.values() for gas in row}


def test_the_inventory_holds_chlorine_and_fluorine_inside_their_brackets_and_follows_the_hosts_iron():
    solar, constants = _inventory(metallicity_fe_h=0.0)
    metal_rich, _ = _inventory(metallicity_fe_h=0.3)
    low, high, _ = mr.INVENTORY["chlorine_wt_percent"]
    assert low <= constants["chlorine_wt_percent"] <= high
    assert solar["Cl"] * CL_AMU == pytest.approx(constants["chlorine_wt_percent"] * 10.0 ** mr.STELLAR_CL_OVER_FE_DEX, rel=1e-3)
    assert solar["F"] * F_AMU == pytest.approx(constants["fluorine_wt_percent"] * 10.0 ** mr.STELLAR_F_OVER_FE_DEX, rel=1e-3)
    assert metal_rich["Cl"] / solar["Cl"] == pytest.approx(10.0 ** 0.3, rel=1e-6)
    assert metal_rich["F"] / solar["F"] == pytest.approx(10.0 ** 0.3, rel=1e-6)


def test_a_body_with_no_star_keeps_the_reference_halogens():
    bare, constants = _inventory()
    assert bare["Cl"] * CL_AMU == pytest.approx(constants["chlorine_wt_percent"] * 10.0 ** mr.STELLAR_CL_OVER_FE_DEX, rel=1e-3)


def test_the_halogens_are_cut_to_the_hydrogen_they_can_pair_with_in_proportion():
    atoms = {"H": 1.0, "O": 0.5, "Cl": 0.3, "F": 0.1}
    constants = {"halogen_per_hydrogen_max": 0.1}
    capped = mr.cap_halogens(atoms, constants)
    assert capped["Cl"] + capped["F"] == pytest.approx(0.1)
    assert capped["Cl"] / capped["F"] == pytest.approx(3.0)
    assert capped["H"] == 1.0 and capped["O"] == 0.5
    assert mr.cap_halogens({"H": 1.0, "Cl": 0.01}, constants) == {"H": 1.0, "Cl": 0.01}      # under the cap: untouched
    assert "Cl" not in mr.cap_halogens({"C": 1.0, "Cl": 0.2}, constants)                       # no hydrogen: the halogens stay in the crust
    assert mr.cap_halogens({"H": 1.0}, constants) == {"H": 1.0}


def test_the_cap_is_drawn_once_inside_its_bracket():
    atoms = {"H": 1.0, "Cl": 5.0}
    constants = {}
    with draw.bound(5):
        mr.cap_halogens(atoms, constants)
        first = constants["halogen_per_hydrogen_max"]
        mr.cap_halogens(atoms, constants)
    low, high, _ = mr.INVENTORY["halogen_per_hydrogen_max"]
    assert low <= first <= high and constants["halogen_per_hydrogen_max"] == first


def test_the_degassing_efficiencies_carry_chlorine_and_fluorine_inside_the_bracket():
    constants = {}
    with draw.bound(1):
        efficiencies = surface_pressure.draw_efficiencies(constants)
    low, high, _ = surface_pressure.DEGAS_EFFICIENCY_HALOGEN
    assert low <= efficiencies["Cl"] <= high and low <= efficiencies["F"] <= high


def test_a_wet_outgassed_air_holds_hydrogen_chloride_and_fluoride_as_the_gases_and_no_atomic_halogen():
    seen = set()
    for seed in range(40):
        with draw.bound(seed):
            result = md.draw_mixture(Subject("terrestrial", {**BASE, "ocean_present": False, "metallicity_fe_h": 0.0}))
        x = result.mole_fractions
        seen |= {gas for gas in ("HCl", "HF") if gas in x}
        assert sum(x.get(gas, 0.0) for gas in ("Cl", "F", "Cl2", "F2")) < 1e-3          # the hydrogen the cap leaves holds the halogens as HCl and HF
    assert seen == {"HCl", "HF"}


def test_the_halides_condense_in_the_cold_and_stay_gas_in_the_warm():
    for formula in ("HCl", "HF"):
        assert saturation.has(formula)
        assert saturation.pressure_pa(formula, 68.0) < saturation.pressure_pa(formula, 200.0) < saturation.pressure_pa(formula, 250.0)
