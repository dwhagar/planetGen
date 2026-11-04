"""GEN.438: the table of the elements and their isotopes (`planetgen/physics/elements.py`, its data files in
`physics/data/`, built by `scripts/build_element_tables.py`). Atomic physics only: the molar masses the gas, materials
and stellar tables use, and the isotope constants the radiogenic heat and the ground's dose read."""

import json
import math
import os
import shutil

import pytest

from planetgen.physics import elements, habitability_world

KEV_J = 1.602176634e-16
"""float: One keV in joules."""

STORED_SPECIES_FORMULAS = {"o2": "O2", "co2": "CO2", "co": "CO", "n2": "N2", "ar": "Ar", "h2": "H2", "h2o": "H2O",
                           "ch4": "CH4", "h2s": "H2S", "so2": "SO2"}
"""dict: The formula of each species the atmosphere stores (`atmosphere.SPECIES`)."""


# --- The data files state where they come from ----------------------------------------------------------------

@pytest.mark.parametrize("name", ["elements.json", "isotopes.json", "decay_series.json"])
def test_every_data_file_states_its_sources_and_when_it_was_read(name):
    with open(os.path.join(elements.DATA_DIRECTORY, name), encoding="utf-8") as handle:
        data = json.load(handle)
    assert data["title"] and data["retrieved"] and "licence" not in data
    assert data["sources"] and all(len(source) > 30 for source in data["sources"])
    assert any("doi" in source or "[R" in source or "ciaaw" in source.lower() or "nist" in source.lower() for source in data["sources"])


@pytest.mark.parametrize("field", ["sources", "retrieved"])
def test_the_loader_refuses_a_file_that_states_no_source(tmp_path, monkeypatch, field):
    shutil.copytree(elements.DATA_DIRECTORY, tmp_path / "data")
    path = tmp_path / "data" / "elements.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    del data[field]
    path.write_text(json.dumps(data), encoding="utf-8")
    monkeypatch.setattr(elements, "DATA_DIRECTORY", str(tmp_path / "data"))
    with pytest.raises(elements.DataError, match=field):
        elements._load_elements()


def test_the_data_files_are_what_the_build_script_makes_from_the_published_tables():
    """The script's parsers read the same layouts the files came from (the source tables themselves are not committed)."""
    import importlib.util
    path = os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "build_element_tables.py")
    spec = importlib.util.spec_from_file_location("build_element_tables", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module._number("1.00782503223(9)") == pytest.approx((1.00782503223, 9e-11))
    assert module._number("0.999885(70)") == pytest.approx((0.999885, 7e-5))
    assert module._number("") == (None, None)
    assert module._modes("IS=0.0117 1;B-=89.28 13;B+=10.72 13") == ([["B-", 89.28], ["B+", 10.72]], 0.0117)
    assert module._modes("A~100;B-=0.276 4") == ([["A", 99.724, "~"], ["B-", 0.276]], None)
    row = ("040 0190   40K    -33535.50       0.06                                  1.248 Gy 0.003  4-*           17      "
           "    1935 IS=0.0117 1;B-=89.28 13;B+=10.72 13")
    state = module.read_nubase(row)[(19, 40)]
    assert state["half_life_s"] == pytest.approx(1.248 * 31557600e9) and state["mass_excess_kev"] == -33535.50
    assert module.read_nubase("206 0820   206Pb  -23785.5        1.1                                 stbl       >2.5Zy 0+"
                              "            08          1927 IS=24.1 30;A ?")[(82, 206)]["stable"]


# --- The elements ----------------------------------------------------------------------------------------------

def test_the_table_holds_the_118_elements_in_order_with_their_symbols():
    assert [item.z for item in elements.ELEMENTS] == list(range(1, 119))
    assert len({item.symbol for item in elements.ELEMENTS}) == 118 and len({item.name for item in elements.ELEMENTS}) == 118
    assert [elements.element(z).symbol for z in (1, 2, 6, 8, 26, 92, 118)] == ["H", "He", "C", "O", "Fe", "U", "Og"]
    assert elements.element("iron") is elements.element("Fe") is elements.element(26)
    for bad in ("Xx", 0, 119, "unobtainium"):
        with pytest.raises(elements.UnknownElement):
            elements.element(bad)


def test_the_atomic_weights_are_ciaaws_to_five_digits_and_an_interval_element_has_one_value_and_its_interval():
    known = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "Mg": 24.305, "Si": 28.085, "S": 32.06, "Cl": 35.45,
             "Fe": 55.845, "Pb": 207.2, "U": 238.03}
    for symbol, weight in known.items():
        assert elements.element(symbol).atomic_weight == weight
    for symbol in ("H", "Li", "B", "C", "N", "O", "Mg", "Si", "S", "Cl", "Ar", "Br", "Tl"):
        item = elements.element(symbol)
        low, high = item.standard_interval
        assert low - item.atomic_weight_uncertainty <= item.atomic_weight <= high + item.atomic_weight_uncertainty, symbol
        assert item.standard_weight is None
    assert elements.element("Fe").standard_interval is None and elements.element("Fe").standard_weight == 55.845


def test_an_element_without_a_standard_weight_names_its_longest_lived_isotope_and_has_no_molar_mass():
    technetium = elements.element("Tc")
    assert technetium.atomic_weight is None and technetium.longest_lived_mass_number in (97, 98)
    with pytest.raises(elements.DataError, match="no standard atomic weight"):
        technetium.molar_mass
    assert elements.element("Pu").longest_lived_mass_number == 244
    assert elements.element("U").longest_lived_mass_number is None and elements.element("U").molar_mass == pytest.approx(0.23803)


@pytest.mark.parametrize("symbol", [item.symbol for item in elements.ELEMENTS if item.has_natural_isotopes])
def test_the_natural_abundances_add_up_to_one_and_give_the_atomic_weight(symbol):
    """The isotope table (NIST, from IUPAC) and the weights table (CIAAW 2021) are two sources; they agree to 3e-4 of the
    weight (selenium and molybdenum differ most, by the 2013 and 2019 revisions of their weights)."""
    item = elements.element(symbol)
    natural = [isotope for isotope in elements.isotopes_of(symbol) if isotope.abundance is not None]
    if not natural:
        pytest.skip("no isotope of this element has a standard abundance")
    assert math.fsum(isotope.abundance for isotope in natural) == pytest.approx(1.0, abs=1e-4)
    mean = math.fsum(isotope.abundance * isotope.mass_u for isotope in natural)
    low, high = item.standard_interval or (item.atomic_weight, item.atomic_weight)
    assert low * (1 - 3e-4) <= mean <= high * (1 + 3e-4)


def test_an_isotopes_molar_mass_is_its_atomic_mass_and_a_radioactive_one_has_a_half_life():
    fe56 = elements.isotope("Fe", 56)
    assert fe56 is elements.isotope("Fe-56") and fe56.molar_mass == pytest.approx(0.055934936, rel=1e-8)
    assert fe56.stable and fe56.abundance == pytest.approx(0.91754, abs=1e-4)
    for isotope in elements.ISOTOPES:
        assert isotope.mass_u == pytest.approx(isotope.a, rel=2e-2) and isotope.name == f"{elements.element(isotope.z).symbol}-{isotope.a}"
        if not isotope.stable and isotope.abundance is None:
            continue
        assert isotope.stable == (isotope.half_life_s is None) or isotope.half_life_s > 3e14
    with pytest.raises(elements.DataError):
        elements.isotope("Fe-56").decay_constant
    with pytest.raises(KeyError):
        elements.isotope("Fe-1")


# --- Molar masses ----------------------------------------------------------------------------------------------

def test_a_formula_is_read_with_groups_and_hydrates():
    assert elements.parse_formula("H2O") == {"H": 2, "O": 1}
    assert elements.parse_formula("Ca3(PO4)2") == {"Ca": 3, "P": 2, "O": 8}
    assert elements.parse_formula("Fe[OH]3") == {"Fe": 1, "O": 3, "H": 3}
    assert elements.parse_formula("CaSO4.2H2O") == elements.parse_formula("CaSO4·2H2O") == \
        {"Ca": 1, "S": 1, "O": 6, "H": 4}
    assert elements.parse_formula("(NH4)2SO4") == {"N": 2, "H": 8, "S": 1, "O": 4}
    assert elements.parse_formula("Ar") == {"Ar": 1}
    for bad in ("", "Xx2", "CO2)", "(CO2", "co2", "C0O"):
        with pytest.raises(ValueError):
            elements.parse_formula(bad)


def test_a_molar_mass_is_the_sum_of_the_atoms_in_kg_per_mol():
    assert elements.molar_mass("H2O") == pytest.approx(2 * 0.001008 + 0.015999)
    assert elements.molar_mass("CaSO4.2H2O") == pytest.approx(0.172164, abs=1e-6)
    assert elements.molar_mass("O2") == pytest.approx(0.031998)
    with pytest.raises(elements.DataError):
        elements.molar_mass("TcO4")


def test_the_molar_mass_of_every_compound_in_the_stored_species_is_the_sum_of_its_elements():
    """The one table of compounds in the generator today is `habitability_world.MOLAR_KG` (the ten species the atmosphere
    stores). Stored to four or five digits from older weights, they agree with the elements to 2e-4 (sulfur dioxide is the
    furthest, 0.06407 against 0.064058). The gas and materials tables of GEN.410 and GEN.427 take their masses from
    `molar_mass`, and add their compounds to this check."""
    from planetgen.physics import atmosphere
    assert set(STORED_SPECIES_FORMULAS) == set(atmosphere.SPECIES) == set(habitability_world.MOLAR_KG)
    for species, formula in STORED_SPECIES_FORMULAS.items():
        assert habitability_world.MOLAR_KG[species] == pytest.approx(elements.molar_mass(formula), rel=2.5e-4), species
    assert habitability_world.CARBON_KG_MOL == pytest.approx(elements.element("C").molar_mass, rel=1e-9)
    assert habitability_world.NITROGEN_KG_MOL == pytest.approx(elements.element("N").molar_mass, rel=1e-9)


# --- The radioactive isotopes ----------------------------------------------------------------------------------

LONG_LIVED = ("U-238", "U-235", "Th-232", "K-40")


@pytest.mark.parametrize("name", LONG_LIVED)
def test_the_half_lives_agree_with_dyes_decay_constants_to_half_a_percent(name):
    isotope = elements.isotope(name)
    assert isotope.decay_constant == pytest.approx(isotope.dye_2012["decay_constant_per_s"], rel=5e-3)
    assert isotope.half_life_s == pytest.approx(
        {"U-238": 4.463e9, "U-235": 7.04e8, "Th-232": 1.40e10, "K-40": 1.248e9}[name] * 365.25 * 86400, rel=1e-3)


@pytest.mark.parametrize("name,tolerance", [("U-238", 2e-3), ("U-235", 3e-3), ("Th-232", 2e-3), ("K-40", 3e-2)])
def test_the_heat_per_kilogram_reproduces_dyes_table(name, tolerance):
    """Avogadro's constant times the decay constant times the energy per decay over the molar mass is Dye's equation 6; his
    heat constants (Table 2) are reproduced within 0.3 percent for uranium and thorium and 2.2 percent for potassium-40,
    whose energy per decay (Lugaro et al., the newer source, used here) is 4.5 percent below his."""
    isotope = elements.isotope(name)
    assert elements.heat_per_kg_of_isotope(name) == pytest.approx(isotope.dye_2012["heat_uw_per_kg"] * 1e-6, rel=tolerance)
    assert isotope.decay_energy_j == pytest.approx(isotope.dye_2012["q_heat_pj"] * 1e-12, rel=5e-2)


def test_dyes_energies_add_up_and_the_heat_is_the_energy_less_the_neutrinos():
    for name in LONG_LIVED:
        dye = elements.isotope(name).dye_2012
        assert dye["q_pj"] - dye["q_neutrino_pj"] == pytest.approx(dye["q_heat_pj"], abs=2e-3)


def test_short_lived_radionuclides_have_a_half_life_and_the_energy_of_a_decay():
    expected_half_life_years = {"Al-26": 7.17e5, "Fe-60": 2.62e6, "Cl-36": 3.01e5, "Be-10": 1.387e6, "Mn-53": 3.7e6,
                                "Pd-107": 6.5e6, "Hf-182": 8.9e6, "I-129": 1.614e7, "Cm-247": 1.56e7, "Sm-146": 6.8e7,
                                "Pu-244": 8.0e7, "Ca-41": 9.94e4, "Pb-205": 1.7e7, "Nb-92": 3.47e7}
    for name, years in expected_half_life_years.items():
        isotope = elements.isotope(name)
        assert isotope.half_life_s / (365.25 * 86400) == pytest.approx(years, rel=0.02), name
        assert isotope.decay_energy_j > 0
    al26 = elements.isotope("Al-26")
    assert al26.decay_energy_j / KEV_J / 1000 == pytest.approx(3.12, rel=2e-3), "Lichtenberg et al.'s 3.12 MeV per decay"
    assert 1 / al26.decay_constant / (1e6 * 365.25 * 86400) == pytest.approx(1.03, rel=5e-3), "the mean life, 1.03 Myr"


def test_the_activity_of_a_ppm_of_each_element_is_the_known_figure():
    """Specific activity per ppm by mass of the element (Bq/kg per ppm): uranium 12.35 (uranium-238 alone), thorium 4.06, and
    per weight percent of potassium 316 (the figures part E of docs/design/composition/body-composition-research.md computed)."""
    uranium = elements.specific_activity("U-238") * elements.isotopic_mass_fraction("U-238") * 1e-6
    thorium = elements.specific_activity("Th-232") * elements.isotopic_mass_fraction("Th-232") * 1e-6
    potassium = elements.specific_activity("K-40") * elements.isotopic_mass_fraction("K-40") * 1e-2
    assert uranium == pytest.approx(12.35, rel=5e-3)
    assert thorium == pytest.approx(4.06, rel=5e-3)
    assert potassium == pytest.approx(316.0, rel=1e-2)


# --- The decay chains ------------------------------------------------------------------------------------------

@pytest.mark.parametrize("parent,end,alphas", [("U-238", "Pb-206", 8), ("U-235", "Pb-207", 7), ("Th-232", "Pb-208", 6)])
def test_each_chain_runs_from_its_parent_to_a_stable_lead_isotope_and_every_path_releases_the_same(parent, end, alphas):
    chain = elements.series(parent)
    assert chain.parent == parent and chain.end == end and chain.members[0] == parent and chain.members[-1] == end
    assert elements.isotope(end).stable
    for name in chain.members[:-1]:
        assert elements.isotope(name).half_life_s, name
        shares = {}
        for step in chain.daughters(name):
            shares[step.daughter] = shares.get(step.daughter, 0.0) + step.percent
        assert math.fsum(shares.values()) == pytest.approx(100.0, abs=0.02), f"the branches of {name} account for all its decays"
    assert not chain.daughters(end)

    def paths(name, taken):
        steps = chain.daughters(name)
        if not steps:
            yield taken
        for step in steps:
            yield from paths(step.daughter, taken + [step])

    totals = {round(math.fsum(step.q_kev for step in path), 3) for path in paths(parent, [])}
    assert len(totals) == 1, "the Q values telescope: every branch of a chain releases the same in all"
    assert {sum(1 for step in path if step.mode == "alpha") for path in paths(parent, [])} == {alphas}
    total_j = totals.pop() * KEV_J
    dye = elements.isotope(parent).dye_2012
    assert total_j == pytest.approx(dye["q_pj"] * 1e-12, rel=2e-3), "the chain's Q against Dye's"
    lugaro = elements.isotope(parent).decay_energy_j
    assert lugaro == pytest.approx(total_j - dye["q_neutrino_pj"] * 1e-12, rel=1.5e-2), "the energy less the neutrinos"


def test_the_radon_isotopes_the_ground_dose_reads_are_in_the_chains():
    assert "Rn-222" in elements.series("U-238").members and "Rn-220" in elements.series("Th-232").members
    assert elements.isotope("Rn-222").half_life_s == pytest.approx(3.8215 * 86400, rel=1e-4)
    assert elements.isotope("Rn-220").half_life_s == pytest.approx(55.6, rel=1e-2)
    assert elements.isotope("Ra-226").half_life_s / (365.25 * 86400) == pytest.approx(1600, rel=1e-3)
