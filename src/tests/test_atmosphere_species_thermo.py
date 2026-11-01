# tests/test_atmosphere_species_thermo.py

"""The NASA Glenn species thermodynamics (GEN.451, GEN.454, GEN.499): the data file, the loader and the free energies."""

import json
import math

import pytest

from planetgen.physics import elements, gases, species_thermo as st
from planetgen.physics.elements import DataError


def test_the_file_holds_the_gases_and_the_condensed_phases_the_chemistry_needs():
    assert len(st.SPECIES) == 104
    assert len(st.GASES) == 84 and len(st.CONDENSED) == 20
    for name in ("H2", "H2O", "CO", "CO2", "CH4", "NH3", "N2", "H2S", "SO2", "Na", "K", "SiO", "TiO", "VO", "Fe", "Mg"):
        assert name in st.GASES
    for name in ("H2O(cr)", "H2O(L)", "C(gr)", "Fe(L)", "SiO2(L)"):
        assert name in st.CONDENSED


LAVA_VAPOUR_NAMES = {"AL": "Al", "ALO": "AlO", "AL2O": "Al2O", "ALO2": "AlO2", "AL2O2": "Al2O2", "K2": "K2", "KO": "KO", "K2O": "K2O", "Na2O": "Na2O"}
"""dict: The lava vapour's aluminium and alkali species (GEN.751, GEN.818) by their name in the NASA file (which spells aluminium AL) and in the gas table."""


def test_the_aluminium_and_alkali_species_of_the_lava_vapour_are_in_the_table_under_the_gas_tables_formulas():
    for name, formula in LAVA_VAPOUR_NAMES.items():
        assert name in st.GASES and st.SPECIES[name].formula == formula
        assert gases.has(formula), formula
        row = gases.gas(formula)
        assert row.molar_mass == pytest.approx(st.SPECIES[name].molar_mass, abs=1e-5)
        # the Cp of the table is the database's, at 300 and 1000 K
        for temperature in (300, 1000):
            assert row.cp_j_mol_k[temperature] == pytest.approx(st.SPECIES[name].cp_over_r(float(temperature)) * 8.314462618, abs=0.01)


def test_every_species_has_elements_a_molar_mass_that_adds_up_and_joined_fits():
    for row in st.SPECIES.values():
        assert row.elements
        assert row.molar_mass == pytest.approx(sum(n * elements.element(s).atomic_weight for s, n in row.elements.items()) * 1e-3)
        for earlier, later in zip(row.fits, row.fits[1:]):
            assert earlier.t_max_k == later.t_min_k
        assert len(row.fits[0].a) == 7 and len(row.fits[0].b) == 2


def test_the_gas_table_names_agree_with_the_database_for_the_gases_both_hold():
    table = {row.formula: row for row in gases.GASES}
    shared = [st.SPECIES[name] for name in st.GASES if st.SPECIES[name].formula in table]
    assert len(shared) >= 30
    for row in shared:
        assert row.molar_mass == pytest.approx(table[row.formula].molar_mass, abs=1.5e-5), row.name
    assert st.species("COS").formula == "OCS" and st.species("C2H2,acetylene").formula == "C2H2"


def test_the_free_energy_matches_the_formation_data_at_298_k():
    # H2O: the enthalpy of formation of the gas is -241.826 kJ/mol; H/RT at 298.15 K is that over R T.
    row = st.species("H2O")
    assert row.h_over_rt(298.15) * st.GAS_CONSTANT * 298.15 == pytest.approx(-241826.0, rel=2e-4)
    assert st.species("O2").h_over_rt(298.15) == pytest.approx(0.0, abs=1e-3)
    assert st.species("O2").s_over_r(298.15) * st.GAS_CONSTANT == pytest.approx(205.15, abs=0.05)


@pytest.mark.parametrize("temperature, expected", [(1500.0, 5.311), (2173.0, 2.303)])
def test_the_co_oxidation_constant_is_the_value_of_the_research(temperature, expected):
    assert st.log10_k({"CO2": 1.0, "CO": -1.0, "O2": -0.5}, temperature) == pytest.approx(expected, abs=0.002)


def test_the_methane_constant_is_zahnle_and_marleys():
    # [R4863]: K = pCH4 pH2O / (pCO pH2^3) = 5.24e-14 exp(27285/T), in atm^-2; our log10 K is in bar^-2.
    for temperature in (700.0, 1000.0, 1300.0):
        ours = st.log10_k({"CH4": 1.0, "H2O": 1.0, "CO": -1.0, "H2": -3.0}, temperature)
        theirs = math.log10(5.24e-14 * math.exp(27285.0 / temperature) * 1.01325 ** 2)   # atm^-2 to bar^-2
        assert ours == pytest.approx(theirs, abs=0.12), temperature


def test_a_temperature_outside_the_fits_is_refused():
    with pytest.raises(ValueError, match="no thermodynamic fit"):
        st.g_over_rt("H2O", 100.0)
    with pytest.raises(ValueError, match="no thermodynamic fit"):
        st.g_over_rt("H2O(L)", 800.0)


def test_the_heat_capacity_of_an_ideal_diatomic_gas_is_near_7_2_r_at_300_k():
    assert st.species("N2").cp_over_r(300.0) == pytest.approx(3.5, abs=0.02)


def test_the_loader_refuses_a_row_that_does_not_add_up(tmp_path, monkeypatch):
    document = json.load(open(st.DATA_FILE, encoding="utf-8"))
    document["species"]["H2"]["molar_mass_g_mol"] = 3.0
    path = tmp_path / "species.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    monkeypatch.setattr(st, "DATA_FILE", str(path))
    with pytest.raises(DataError, match="molar mass"):
        st._load()
    document["species"]["H2"]["molar_mass_g_mol"] = 2.016
    document["species"]["H2"]["intervals"][1]["t_min_k"] = 1200.0
    path.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(DataError, match="do not join"):
        st._load()


def test_the_build_script_reproduces_the_committed_file(tmp_path):
    # The raw database is not in the repository: build a small thermo.inp of three species in its format and check the reader.
    import importlib.util
    import os
    spec = importlib.util.spec_from_file_location("build_species_thermo", os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "build_species_thermo.py"))
    build = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(build)
    h2 = st.species("H2")
    fit = h2.fits[0]
    text = (
        "!comment\n"
        "H2                Ref-Elm. Gurvich,1978 pt1 p103 pt2 p31.\n"
        " 1 tpis78 H   2.00    0.00    0.00    0.00    0.00 0    2.0158800          0.000\n"
        f"{fit.t_min_k:11.3f}{fit.t_max_k:11.3f}7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8467.340\n"
        + "".join(f"{value:16.9E}".replace("E", "D") for value in fit.a[:5]) + "\n"
        + "".join(f"{value:16.9E}".replace("E", "D") for value in fit.a[5:]) + " " * 16
        + "".join(f"{value:16.9E}".replace("E", "D") for value in fit.b) + "\n"
    )
    path = tmp_path / "thermo.inp"
    path.write_text(text, encoding="utf-8")
    read = build.read_species(str(path))
    assert read["H2"]["elements"] == {"H": 2.0}
    assert read["H2"]["intervals"][0]["a"] == pytest.approx(list(fit.a), rel=1e-8)
    assert read["H2"]["intervals"][0]["b"] == pytest.approx(list(fit.b), rel=1e-8)


HALOGEN_NAMES = {"F": "F", "F2": "F2", "CL": "Cl", "CL2": "Cl2", "HCL": "HCl", "HF": "HF"}
"""dict: The halogen gases of the outgassed air (GEN.524 part 2) by their name in the NASA file (which spells chlorine CL) and in the gas table."""


def test_the_halogen_species_are_in_the_table_and_hydrogen_chloride_leads_the_chlorine_of_a_wet_air():
    from planetgen.physics import gas_equilibrium
    for name, formula in HALOGEN_NAMES.items():
        assert name in st.GASES and st.SPECIES[name].formula == formula
    for formula in ("HCl", "HF"):
        assert gases.has(formula)
    for temperature in (500.0, 1000.0, 1500.0):
        result = gas_equilibrium.solve({"H": 2.0, "O": 1.0, "C": 0.3, "Cl": 0.002, "F": 0.001}, temperature, 100.0)
        assert result.moles["HCL"] > 0.99 * 0.002          # the dominant chlorine gas below 1500 K [R3945]
        assert result.moles["HF"] > 0.99 * 0.001
