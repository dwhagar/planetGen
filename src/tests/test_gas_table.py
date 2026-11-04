"""GEN.410: the table of the gases and compounds that can be in an atmosphere (`planetgen/physics/gases.py`, its data
file `physics/data/gases.json`, built by `scripts/build_gas_table.py` from the research parts 01, 02 and 03). The audit at
the end fails when a gas the code or the docs name is not a row of the table."""

import json
import os
import re

import pytest

from planetgen import tuning
from planetgen.physics import atmosphere, elements, gases, habitability_world
from planetgen.physics.elements import DataError

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
TOXICITY_NOTE = os.path.join(REPO, "docs", "design", "research", "imported", "atmospheric-toxicity.md")


# --- The file states where each row comes from -----------------------------------------------------------------

def test_every_row_has_a_source_and_a_state():
    assert len(gases.GASES) >= 72
    for row in gases.GASES:
        assert row.source and row.state in ("verified", "partly verified", "open", "stand-in"), row.formula
        if row.state != "verified":
            assert row.limit, f"{row.formula} is {row.state} and names no limit"
        if row.rayleigh and row.rayleigh.state != "verified":
            assert row.rayleigh.limit, row.formula


def test_the_formulas_are_unique_and_every_name_is_set():
    formulas = [row.formula for row in gases.GASES]
    assert len(formulas) == len(set(formulas))
    names = [row.name for row in gases.GASES]
    assert all(names) and len(names) == len(set(names))


def test_the_molar_masses_come_from_the_element_table():
    for row in gases.GASES:
        assert abs(row.molar_mass * 1000.0 - row.research_molar_mass) <= gases.MOLAR_MASS_TOLERANCE_G_MOL, row.formula
    assert gases.gas("H2O").molar_mass == pytest.approx(elements.molar_mass("H2O"), rel=0, abs=0)
    assert gases.gas("(CH3)2S").molar_mass == pytest.approx(0.06213, abs=1e-6)
    # An isotope species is a sum of isotope masses; radon has no standard weight.
    assert gases.gas("3He").molar_mass == pytest.approx(3.0160293e-3, abs=1e-8)
    assert gases.gas("Rn").molar_mass == pytest.approx(0.2220176, abs=1e-6)
    assert gases.gas("HD").isotopes == {"H-1": 1, "H-2": 1}


def test_the_ten_stored_partials_are_rows_and_agree_with_the_stored_molar_masses():
    for key in gases.STORED_PARTIALS:
        row = gases.stored(key)
        assert row.molar_mass == pytest.approx(habitability_world.MOLAR_KG[key], rel=2e-3), key
    with pytest.raises(gases.UnknownGas):
        gases.stored("he")


def test_an_unknown_gas_is_refused_with_the_way_to_add_it():
    with pytest.raises(gases.UnknownGas, match="gases.json"):
        gases.gas("XeF2")
    assert gases.has("CO2") and not gases.has("co2x")


def test_the_flags_follow_the_rules_of_part_01():
    assert gases.gas("N2").ir == "cia" and gases.gas("O2").ir == "cia" and gases.gas("H2").ir == "cia"
    assert gases.gas("He").ir == "none" and gases.gas("Ar").ir == "none"
    assert gases.gas("CO2").is_infrared_active and gases.gas("CH4").is_infrared_active
    assert not gases.gas("N2").is_infrared_active
    assert gases.gas("CO2").boiling_point_is_sublimation and not gases.gas("CH4").boiling_point_is_sublimation
    assert gases.gas("O3").absorbs_near_uv and not gases.gas("N2").absorbs_near_uv
    # Volatility classes follow the boiling point (UV below 30 K, SV 30 to 150, V 150 to 350, SmV 350 to 1000, Ref above).
    for row in gases.GASES:
        if row.volatility and row.boiling_point_k and row.volatility != "Ref":
            low, high = {"UV": (0, 30), "SV": (30, 150), "V": (150, 350), "SmV": (350, 1000)}[row.volatility]
            assert low <= row.boiling_point_k < high or row.boiling_point_is_sublimation, row.formula


def test_the_phase_points_are_ordered():
    for row in gases.GASES:
        if row.triple_point_k and row.critical_temperature_k:
            assert row.triple_point_k < row.critical_temperature_k, row.formula
        if row.boiling_point_k and row.critical_temperature_k:
            assert row.boiling_point_k < row.critical_temperature_k, row.formula
        if row.gamma_300 is not None:
            assert 1.0 < row.gamma_300 <= 5 / 3 + 1e-3, row.formula


# --- Rayleigh: the reference temperature of each refractivity fit (GEN.408) ------------------------------------------

def test_each_refractivity_fit_stores_the_temperature_it_was_made_at():
    for formula in ("N2", "Ar", "CO2", "CH4", "N2O"):
        assert gases.gas(formula).rayleigh.reference_temperature_k == 288.15, formula
    for formula in ("O2", "CO", "H2", "NH3"):
        assert gases.gas(formula).rayleigh.reference_temperature_k == 273.15, formula
    # The polarizability route has no reference temperature, and neither has water vapor's (unknown).
    for formula in ("H2S", "SO2", "C2H6", "HCN", "H2O"):
        assert gases.gas(formula).rayleigh.reference_temperature_k is None, formula


def test_the_rayleigh_cross_sections_follow_the_research_values():
    n2 = gases.gas("N2").rayleigh
    assert n2.sigma_cm2[550] == pytest.approx(4.635e-27, rel=1e-3)
    assert n2.sigma_cm2[450] > n2.sigma_cm2[550] > n2.sigma_cm2[650]
    # CO takes the Allen/Cox form and CH4 the He 2021 dispersion, not the coded older form that is 30 percent high.
    assert "Allen/Cox" in gases.gas("CO").rayleigh.research_row
    assert "He 2021" in gases.gas("CH4").rayleigh.research_row
    assert gases.gas("CO").rayleigh.sigma_cm2[550] == pytest.approx(5.734e-27, rel=1e-3)
    assert gases.gas("CH4").rayleigh.sigma_cm2[550] == pytest.approx(9.905e-27, rel=1e-3)
    # H2 and He against N2 (ratios of part 02: 0.217 by the refractivity, 0.0132).
    assert gases.gas("H2").rayleigh.sigma_cm2[550] / n2.sigma_cm2[550] == pytest.approx(0.217, abs=0.002)
    assert gases.gas("He").rayleigh.sigma_cm2[550] / n2.sigma_cm2[550] == pytest.approx(0.0132, abs=0.0003)
    for row in gases.GASES:
        if row.rayleigh:
            assert row.rayleigh.king_factor >= 1.0, row.formula


def test_water_vapor_uses_the_vapor_depolarization_not_the_liquids():
    water = gases.gas("H2O").rayleigh
    assert water.king_factor == pytest.approx(1.0010, abs=1e-4)
    # No reference temperature is known for its refractivity, so the generic polarizability value is used (9 to 21 percent
    # above the n-based 3.11e-27 and 2.80e-27).
    assert water.sigma_cm2[550] == pytest.approx(3.38e-27, rel=1e-3)
    assert water.state == "partly verified" and "temperature" in water.limit


def test_the_gases_whose_numbers_are_not_found_stay_marked():
    assert gases.gas("NH3").rayleigh.king_factor_assumed and gases.gas("NH3").rayleigh.state == "partly verified"
    for formula in ("H2S", "SO2"):
        assert gases.gas(formula).rayleigh.state == "partly verified"
        assert "polarizability" in gases.gas(formula).rayleigh.limit
    assert gases.gas("CO").rayleigh.state == "partly verified"
    assert gases.gas("N2").rayleigh.state == "verified"


def test_the_henry_law_returns_the_tabulated_value_at_298_and_follows_d():
    for row in gases.GASES:
        if row.henry is None:
            continue
        assert row.henry.kh(298.15) == pytest.approx(row.henry.kh_298, rel=1e-12)
        assert row.henry.kh_range[0] <= row.henry.kh_298 <= row.henry.kh_range[1] or row.formula == "O3"
    co2 = gases.gas("CO2").henry
    # kH(T) = kH(298.15) exp(d (1/T - 1/298.15)): CO2 is about 2.4 times more soluble at 273.15 K (d = 2400 K).
    assert co2.kh(273.15) == pytest.approx(0.034 * 2.718281828 ** (2400 * (1 / 273.15 - 1 / 298.15)), rel=1e-6)
    assert co2.kh(273.15) > co2.kh(298.15) > co2.kh(323.15)
    assert co2.dissociates and not gases.gas("N2").henry.dissociates


def test_the_henry_law_is_refused_outside_the_range_the_research_gives():
    with pytest.raises(gases.NoData, match="273.15"):
        gases.gas("CO2").henry.kh(400.0)
    with pytest.raises(gases.NoData):
        gases.gas("CO2").henry.kh(200.0)


def test_the_henry_entries_name_where_the_sources_disagree():
    assert gases.gas("N2").henry.state == "verified"
    assert gases.gas("H2S").henry.state == "partly verified" and "factor" in gases.gas("H2S").henry.limit
    assert gases.gas("PH3").henry.state == "partly verified" and "one compiled value" in gases.gas("PH3").henry.limit
    assert gases.gas("HCl").henry.state == "partly verified"


def test_every_henry_species_is_a_row_with_the_research_count():
    assert sum(1 for row in gases.GASES if row.henry) == 35
    assert sum(1 for row in gases.GASES if row.viscosity) == 38


def test_chapman_enskog_reproduces_the_research_viscosities_at_300_k():
    checked = 0
    for row in gases.GASES:
        if row.viscosity is None or row.viscosity.polar:
            continue
        assert row.viscosity_upa_s(300.0) == pytest.approx(row.viscosity.research_mu_300_upa_s, abs=0.011), row.formula
        checked += 1
    assert checked == 37


def test_the_viscosity_agrees_with_nist_where_the_research_compared_it():
    for formula, nist in (("N2", 17.88), ("Ar", 22.72), ("CO2", 14.99), ("CH4", 11.12)):
        assert gases.gas(formula).viscosity_upa_s(300.0) == pytest.approx(nist, rel=0.02)
    assert gases.gas("N2").viscosity.state == "verified"
    assert gases.gas("NO").viscosity.state == "partly verified" and "NIST" in gases.gas("NO").viscosity.limit


def test_the_viscosity_rises_with_temperature_for_a_gas_and_refuses_water_vapor():
    assert gases.gas("N2").viscosity_upa_s(600.0) > gases.gas("N2").viscosity_upa_s(300.0)
    with pytest.raises(gases.NoData, match="polar"):
        gases.gas("H2O").viscosity_upa_s(300.0)
    assert gases.gas("H2O").viscosity.nist_300_upa_s == 9.77
    assert gases.gas("H2O").viscosity.state == "partly verified"


def test_the_viscosity_is_refused_outside_the_neufeld_range_and_for_a_gas_without_constants():
    with pytest.raises(gases.NoData, match="Neufeld"):
        gases.gas("N2").viscosity_upa_s(1.0)
    with pytest.raises(gases.NoData, match="no viscosity constants"):
        gases.gas("Na").viscosity_upa_s(1000.0)


def test_the_loader_refuses_a_henry_or_viscosity_entry_that_is_not_physical(broken):
    def zero_henry(data):
        next(row for row in data["gases"] if row["formula"] == "CO2")["henry"]["kh_298_mol_kg_bar"] = 0.0

    def zero_sigma(data):
        next(row for row in data["gases"] if row["formula"] == "N2")["viscosity"]["sigma_angstrom"] = 0.0

    def partly_without_limit(data):
        entry = next(row for row in data["gases"] if row["formula"] == "N2")["viscosity"]
        entry["state"], entry["limit"] = "partly verified", ""

    for edit in (zero_henry, zero_sigma, partly_without_limit):
        with pytest.raises(DataError):
            broken(edit)()


# --- The infrared bands and the ultraviolet bins (the third slice) -------------------------------------------------

def test_every_infrared_band_is_a_number_the_research_cell_states():
    listed = [row for row in gases.GASES if row.infrared]
    assert len(listed) == 37
    for row in listed:
        entry = row.infrared
        text = entry.fundamentals_text.replace(",", "")
        for band in entry.bands:
            for value in (band.wavenumber_cm1, *(band.range_cm1 or ())):
                assert value in (None, 0.0) or f"{value:g}" in text, f"{row.formula}: {value} is not in {entry.fundamentals_text!r}"


def test_the_band_centres_agree_with_the_wavelengths_the_research_prints_beside_them():
    water = {band.label: band for band in gases.gas("H2O").infrared.bands}
    assert water["nu2 bend"].wavelength_um == pytest.approx(6.27, abs=0.005)
    assert water["nu1"].wavelength_um == pytest.approx(2.73, abs=0.005)
    methane = {band.label: band for band in gases.gas("CH4").infrared.bands}
    assert methane["nu4"].wavenumber_cm1 == 1306.2 and methane["nu4"].wavelength_um == pytest.approx(7.66, abs=0.005)
    carbon_dioxide = gases.gas("CO2").infrared
    assert [band.wavenumber_cm1 for band in carbon_dioxide.bands if band.label == "hot band"] == [960, 1063]
    assert all(band.approximate for band in carbon_dioxide.bands), "the CO2 bands are recalled, not read"
    assert not any(band.approximate for band in gases.gas("CH4").infrared.bands)
    assert gases.gas("CH3OH").infrared.bands[2].wavenumber_cm1 is None and gases.gas("CH3OH").infrared.bands[2].range_cm1 == (1455.0, 1477.0)


def test_the_saturated_absorption_rank_is_the_one_the_research_states():
    ranks = {row.formula: row.infrared.saturated_rank for row in gases.GASES if row.infrared and row.infrared.saturated_rank}
    assert ranks == {"H2O": 11, "CO2": 16, "CH4": 27, "C2H2": 6, "SO2": 8, "CH3OH": 1}


def test_a_gas_without_fundamentals_has_no_band_and_says_why():
    for formula in ("N2", "O2", "H2", "He", "Ne", "Ar", "Kr", "Xe", "TiO", "Na", "HNO3"):
        entry = gases.gas(formula).infrared
        assert entry.bands == () and entry.strength_text, formula
    assert "collision-induced" in gases.gas("N2").infrared.fundamentals_text
    assert gases.gas("CH4").infrared.state == "verified" and gases.gas("CH4").infrared.limit == ""
    assert gases.gas("H2O").infrared.state == "partly verified" and "rotation edge RECALLED" in gases.gas("H2O").infrared.limit


def test_recalled_band_positions_are_partly_verified_with_their_limit():
    for formula in ("CO2", "O3", "N2O", "SO2", "CO", "OCS", "SiO"):
        entry = gases.gas(formula).infrared
        assert entry.state == "partly verified" and "RECALLED" in entry.limit, formula


def test_the_ultraviolet_bins_run_from_10_to_1100_nm_one_after_another():
    bins = gases.UV_BINS_NM
    assert len(bins) == 35 and bins[0] == (10.0, 30.0) and bins[-1] == (800.0, 1100.0)
    assert all(a[1] == b[0] for a, b in zip(bins, bins[1:]))


def test_the_ultraviolet_values_are_the_ones_the_research_table_prints():
    assert gases.gas("O3").uv.at(255.0) == 1.1e-17                 # the Hartley band, 240 to 260 nm
    assert gases.gas("O3").uv.at(620.0) == 3.8e-21                 # the Chappuis band
    assert gases.gas("O2").uv.at(145.0) == 1.4e-17                 # the Schumann-Runge continuum
    assert gases.gas("NO2").uv.at(420.0) == 5.6e-19                # the only broad visible band
    assert gases.gas("CH4").uv.at(900.0) == 1.6e-24
    assert gases.gas("CH2O").uv.at(310.0) == 2.8e-20               # the research's H2CO is the table's CH2O
    assert gases.gas("SO2").uv.at(30.0) == 3.5e-17                 # the edge of a bin belongs to the bin above it
    assert gases.gas("SO2").uv.at(1100.0) == 0.0


def test_a_bin_outside_the_data_is_none_and_counts_as_zero():
    h2 = gases.gas("H2").uv
    assert h2.coverage_nm == (0.1, 119.9)
    assert h2.cross_section_cm2[gases.UV_BINS_NM.index((120.0, 123.0))] is None and h2.at(121.6) == 0.0
    assert gases.gas("CO2").uv.partial_bin[gases.UV_BINS_NM.index((210.0, 220.0))] is True


def test_every_cell_is_missing_exactly_where_the_bin_lies_outside_the_coverage():
    listed = [row for row in gases.GASES if row.uv]
    assert len(listed) == 22
    for row in listed:
        low, high = row.uv.coverage_nm
        for (bin_low, bin_high), value in zip(gases.UV_BINS_NM, row.uv.cross_section_cm2):
            assert (value is None) == (bin_high <= low or bin_low >= high), (row.formula, bin_low, bin_high)


def test_the_ultraviolet_search_for_the_longest_absorbing_wavelength():
    assert gases.gas("SO2").uv.longest_wavelength_above(1e-20) == 320.0       # the 300 to 320 nm bin is the last at 1e-20
    assert gases.gas("NO2").uv.longest_wavelength_above(1e-20) == 650.0
    assert gases.gas("N2").uv.longest_wavelength_above(1e-10) is None
    assert gases.gas("H2O").uv.longest_wavelength_above(1e-18) < 200.0, "water vapour does not shield the ground ultraviolet"


def test_the_limits_the_research_names_for_the_ultraviolet_bins_are_kept():
    assert "k-distribution" in gases.gas("O2").uv.limit and gases.gas("O2").uv.state == "partly verified"
    assert "[R2258]" in gases.gas("CH4").uv.limit
    assert gases.gas("O3").uv.state == "verified" and gases.gas("O3").uv.limit == ""
    assert "[R2040] below 52.7 nm" in gases.gas("O3").uv.origin
    assert gases.gas("O2").uv.origin == "[R1873] or [R2040]"


def test_the_loader_refuses_an_infrared_or_ultraviolet_entry_that_is_not_physical(broken):
    def row(data, formula):
        return next(item for item in data["gases"] if item["formula"] == formula)

    def no_centre(data):
        row(data, "CH4")["infrared"]["bands"][0]["wavenumber_cm1"] = None

    def negative_centre(data):
        row(data, "CH4")["infrared"]["bands"][0]["wavenumber_cm1"] = -1.0

    def wrong_band_count(data):
        row(data, "N2")["infrared"]["bands"] = row(data, "CH4")["infrared"]["bands"]

    def short_values(data):
        row(data, "O3")["uv"]["cross_section_cm2"].pop()

    def negative_value(data):
        row(data, "O3")["uv"]["cross_section_cm2"][3] = -1e-17

    def bins_apart(data):
        data["uv_bins_nm"][4][1] += 1.0

    def no_bins(data):
        data.pop("uv_bins_nm")

    for edit in (no_centre, negative_centre, wrong_band_count, short_values, negative_value, bins_apart, no_bins):
        with pytest.raises(DataError):
            broken(edit)()


# --- The loader refuses a file or a row that does not say where it comes from -------------------------------------

@pytest.fixture
def broken(tmp_path, monkeypatch):
    """Returns `edit(function)`: copies the data file with `function(data)` applied and makes the loader read it."""
    original = gases.DATA_FILE

    def edit(function):
        data = json.loads(open(original, encoding="utf-8").read())
        function(data)
        path = tmp_path / "gases.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        monkeypatch.setattr(gases, "DATA_FILE", str(path))
        return gases._load
    return edit


@pytest.mark.parametrize("field", ["sources", "retrieved", "title"])
def test_the_loader_refuses_a_file_that_states_no_source(broken, field):
    with pytest.raises(DataError, match=field):
        broken(lambda data: data.pop(field))()


@pytest.mark.parametrize("field", ["source", "state"])
def test_the_loader_refuses_a_row_without_a_source_or_state(broken, field):
    with pytest.raises(DataError, match=field):
        broken(lambda data: data["gases"][3].pop(field))()


def test_the_loader_refuses_a_partly_verified_row_without_its_limit(broken):
    def edit(data):
        data["gases"][0]["state"] = "partly verified"
        data["gases"][0]["limit"] = ""
    with pytest.raises(DataError, match="names no limit"):
        broken(edit)()


def test_the_loader_refuses_a_duplicate_formula_and_an_unreadable_one(broken):
    with pytest.raises(DataError, match="twice"):
        broken(lambda data: data["gases"][1].update(formula=data["gases"][0]["formula"]))()
    with pytest.raises(DataError, match="molar mass"):
        broken(lambda data: data["gases"][2].update(formula="Qq9"))()


def test_the_loader_refuses_a_molar_mass_that_the_element_table_does_not_give(broken):
    with pytest.raises(DataError, match="more than"):
        broken(lambda data: data["gases"][2].update(molar_mass_research_g_mol=data["gases"][2]["molar_mass_research_g_mol"] + 1))()


def test_the_data_file_is_what_the_builder_makes():
    """The checked-in file is the research tables' values and nothing else (`scripts/build_gas_table.py --check`)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("build_gas_table", os.path.join(REPO, "scripts", "build_gas_table.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not os.path.isdir(module.PARTS):
        pytest.skip("the research parts are not in this checkout")
    expected = json.dumps(module.build(), indent=1, ensure_ascii=False) + "\n"
    with open(gases.DATA_FILE, encoding="utf-8") as handle:
        assert handle.read() == expected


# --- The audit: every gas the code or the docs name is a row (GEN.410) -----------------------------------------------

NAMES = {row.name: row.formula for row in gases.GASES}

NOT_FORMULAS = {"OH", "NUV", "UV", "ROS", "IW", "QFM", "FTT", "XUV", "FUV", "CIA", "MCP", "NOx", "CME", "CNS", "S0094"}
"""set: Tokens of the toxicity note that look like formulas and are not gases (the hydroxyl radical is an intermediate; the
rest are acronyms, a journal article number, and the NOx family, which is NO and NO2)."""


def _names_in(text):
    """The gas names an atmosphere text lists: its "mix of" list and its "traces of" list, without the adjectives."""
    found = []
    mix, _, traces = re.sub(r"^a .*?mix of (mostly )?", "", text).partition(", with traces of ")
    for part in (mix, traces):
        found.extend(name for name in re.split(r", and |, | and ", part) if name)
    return found


def test_no_planet_class_carries_an_atmosphere_text():
    """GEN.415: a class names no gas and does not say whether it has air; the air is drawn from the body and the text is written from
    the stored mixture."""
    for letter, class_data in tuning.PLANET_CLASSES.items():
        assert "atmosphere" not in class_data, letter
        assert "air" not in class_data, letter


def test_every_gas_the_atmosphere_text_names_is_a_row_of_the_table():
    """GEN.415: the text is written from the partial pressures, so a gas it names is a gas of the stored mixture, named as the table
    names it (one spelling: "water vapor", never "water vapour"), whichever gases the air holds."""
    for start in range(0, len(gases.GASES) - 5):
        rows = gases.GASES[start:start + 6]
        partials = {row.formula: 10.0 ** (3 - 0.8 * i) for i, row in enumerate(rows)}
        text = atmosphere.describe(partials)
        named = _names_in(text)
        assert named, text
        assert set(named) <= set(NAMES), (text, sorted(set(named) - set(NAMES)))
        assert named[0] == rows[0].name, text
    assert atmosphere.describe({"H2O": 1000.0}) == "a humid mix of mostly water vapor"


def _formulas_in(text):
    """The chemical formulas written in the toxicity note (LaTeX `\\text{SO}_2` and plain `CO2`)."""
    text = re.sub(r"\\text\{([^{}]*)\}", r"\1", text)
    text = re.sub(r"_\{?(\d+)\}?", r"\1", text)
    text = re.sub(r"\^\{?[+-]\}?", "", text)
    found = set()
    for token in re.findall(r"(?<![A-Za-z0-9])((?:[A-Z][a-z]?\d*){1,6}(?:\([A-Za-z0-9]+\)\d*)?)(?![A-Za-z0-9])", text):
        if len(token) < 2 or token in NOT_FORMULAS:
            continue
        try:
            elements.parse_formula(token)
        except ValueError:
            continue
        if re.search(r"\d", token) or re.search(r"[A-Z][a-z]?[A-Z]", token):
            found.add(token)
    return found


def test_every_gas_the_toxicity_note_names_is_in_the_table():
    with open(TOXICITY_NOTE, encoding="utf-8") as handle:
        formulas = _formulas_in(handle.read())
    assert {"CO2", "CO", "SO2", "H2S", "CH4"} <= formulas
    missing = sorted(formula for formula in formulas if not gases.has(formula))
    assert not missing, f"gases named in docs/design/research/imported/atmospheric-toxicity.md and missing from the table: {missing}"
