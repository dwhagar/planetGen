#!/usr/bin/env python3
"""
Builds the gas and compound constants file (GEN.410) from the atmosphere research tables it is typed in from.

    python scripts/build_gas_table.py            # write src/planetgen/physics/data/gases.json
    python scripts/build_gas_table.py --check    # build in memory and fail if the committed file differs

The file is `src/planetgen/physics/data/gases.json`. The numbers are read from the markdown tables of the research parts
(`docs/design/composition/atmosphere-composition-parts/`): T1 identity, T2 phase points and T3 heat capacity of part 01, and the
Rayleigh tables 2.1a and 2.1b of part 02; part 01 sections 8 (Henry constants) and 9 (viscosity) add two more blocks; part 03
table 2.1 gives the infrared bands and part 02 section 2.3 the ultraviolet and visible absorption bins (GEN.410 third slice). The names, the choice of which Rayleigh row is adopted for a gas, and the
reference temperature of each refractivity fit are typed in below, each with the reason. Nothing is fitted or invented
here: a number that the tables do not give is `null` and its block says so in its state.
"""

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARTS = os.path.join(HERE, "..", "docs", "design", "composition", "atmosphere-composition-parts")
OUT = os.path.join(HERE, "..", "src", "planetgen", "physics", "data", "gases.json")

PART_01 = "01-gas-table-thermodynamics.md"
PART_02 = "02-gas-optics-uv-visible.md"
PART_03 = "03-infrared-greenhouse-thermal.md"

RETRIEVED = "2026-10-17"
"""str: The day the research parts were merged (the tables' own date)."""

# A tag letter after a number in part 01 says where the value came from (section 0 of part 01).
PHASE_TAGS = {"h": "heos_constants", "e": "NIST EOS table", "n": "NIST WebBook median", "y": "Yaws", "c": "CRC"}
CP_TAGS = {"J": "JANAF 1998", "N": "NASA-7", "S": "Shomate", "f": "NIST fluid at 1e-4 MPa",
           "x": "NASA-7 extrapolated below 200 K (unreliable)", "C": "computed 5R/2", "P": "Poling polynomial"}

# Names of the species (typed in; part 01 gives formulas only), spelled as the code spells them (the gas table's names are the one source of every gas name in the text, GEN.415: American, "water vapor"). The key is the formula as part 01 writes it.
NAMES = {
    "H2": "hydrogen", "HD": "hydrogen deuteride", "He": "helium", "3He": "helium-3", "Ne": "neon", "Ar": "argon",
    "Kr": "krypton", "Xe": "xenon", "Rn": "radon", "N2": "nitrogen", "O2": "oxygen", "O3": "ozone", "H2O": "water vapor",
    "HDO": "semi-heavy water", "H2O2": "hydrogen peroxide", "CO": "carbon monoxide", "CO2": "carbon dioxide",
    "CH4": "methane", "C2H2": "acetylene", "C2H4": "ethylene", "C2H6": "ethane", "C3H8": "propane", "C4H10": "butane",
    "C6H6": "benzene", "HCN": "hydrogen cyanide", "CH2O": "formaldehyde", "CH3OH": "methanol", "NH3": "ammonia",
    "N2H4": "hydrazine", "N2O": "nitrous oxide", "NO": "nitric oxide", "NO2": "nitrogen dioxide",
    "N2O4": "dinitrogen tetroxide", "HNO3": "nitric acid", "HNCO": "isocyanic acid", "C3O2": "carbon suboxide",
    "CN": "cyanogen radical", "C2N2": "cyanogen", "CH3CN": "acetonitrile", "HC3N": "cyanoacetylene",
    "C4H2": "diacetylene", "PH3": "phosphine", "H2S": "hydrogen sulfide", "SO2": "sulfur dioxide",
    "SO3": "sulfur trioxide", "H2SO4": "sulfuric acid", "OCS": "carbonyl sulfide", "CS2": "carbon disulfide",
    "S2": "disulfur", "S8": "octasulfur", "SO": "sulfur monoxide", "CH3SH": "methanethiol",
    "(CH3)2S": "dimethyl sulfide", "HCl": "hydrogen chloride", "HF": "hydrogen fluoride", "Cl2": "chlorine",
    "CH3Cl": "chloromethane", "HBr": "hydrogen bromide", "HClO4": "perchloric acid", "CF4": "tetrafluoromethane",
    "SF6": "sulfur hexafluoride", "SiH4": "silane", "SiO": "silicon monoxide", "Na": "sodium",
    "K": "potassium", "Ca": "calcium", "Mg": "magnesium", "Fe": "iron",
    "TiO": "titanium monoxide", "NaCl": "sodium chloride vapor", "KCl": "potassium chloride vapor",
    "NH4Cl": "ammonium chloride", "FeH": "iron hydride", "CrH": "chromium hydride", "AlH": "aluminum monohydride",
    "VO": "vanadium monoxide", "CaOH": "calcium hydroxide", "NH4SH": "ammonium hydrosulfide",
    "H": "atomic hydrogen", "O": "atomic oxygen", "OH": "hydroxyl radical", "C": "atomic carbon", "N": "atomic nitrogen",
    "S": "atomic sulfur", "Si": "atomic silicon", "Ti": "atomic titanium", "V": "atomic vanadium", "Cr": "atomic chromium",
    "P2": "diphosphorus", "FeO": "iron monoxide", "MgO": "magnesium monoxide", "SiO2": "silicon dioxide vapor",
    "TiO2": "titanium dioxide vapor", "NaOH": "sodium hydroxide vapor", "KOH": "potassium hydroxide vapor",
    "Na2": "disodium", "NaO": "sodium monoxide", "Al": "atomic aluminum", "AlO": "aluminum monoxide", "Al2O": "dialuminum monoxide", "AlO2": "aluminum dioxide", "Al2O2": "dialuminum dioxide", "K2": "dipotassium", "KO": "potassium monoxide", "K2O": "potassium oxide vapor", "Na2O": "sodium oxide vapor", "MgOH": "magnesium hydroxide radical", "Mg(OH)2": "magnesium hydroxide vapor",
    "Fe(OH)2": "iron hydroxide vapor", "CaO": "calcium monoxide", "SiS": "silicon monosulfide", "HCOOH": "formic acid",
}

# Species whose molar mass is a sum of isotope masses, not of standard atomic weights (part 01 T1: "isotope masses for
# 3He, HD, HDO"; radon has no standard weight). The values are isotope names of `physics.elements`.
ISOTOPES = {
    "3He": {"He-3": 1},
    "HD": {"H-1": 1, "H-2": 1},
    "HDO": {"H-1": 1, "H-2": 1, "O-16": 1},
    "Rn": {"Rn-222": 1},
}

# Which Rayleigh row is adopted for each gas, and the reference temperature its refractivity fit was made at (GEN.408 section 2:
# the wrong one is 1.1 percent off for O2 and about 8 percent for codes that mix conventions). The label is the first cell of the
# row of part 02 table 2.1a or 2.1b. `ref_k` is None where the part names no reference temperature.
RAYLEIGH_ROWS = {
    "N2": ("a", "N2", 288.15, "n-based, report 2.2 re-evaluated"),
    "O2": ("a", "O2", 273.15, "n-based; measured 4.55e-27 at 532 nm against n-based 4.642e-27 (He 2021 Table 2)"),
    "Ar": ("a", "Ar", 288.15, "n-based"),
    "CO2": ("a", "CO2", 288.15, "n-based; He 2021 measured 13.3e-27 at 532 nm"),
    "CH4": ("a", "CH4 (He 2021 dispersion)", 288.15,
            "He 2021 dispersion, replacing the older coded form that is 30 percent above the measured value"),
    "CO": ("a", "CO (Allen/Cox form)", 273.15, "Allen/Cox form n - 1 = 32.7e-5 (1 + 8.1e-3/lambda^2) at 273.15 K; no measurement reached"),
    "H2": ("a", "H2, route A: Peck and Huang n-based", 273.15,
           "Peck and Huang 1977 at 273.15 K; Dalgarno and Williams route gives 0.202 of N2 against 0.217: draw 0.20 to 0.22"),
    "He": ("a", "He", None, "Mansfield and Peck n-based; the report formula gives 6.04e-29"),
    "Ne": ("a", "Ne", None, "Bideau-Mehu fit (650 nm extrapolated); Cuthbertson gives 2.263e-28, spread 3.5 percent"),
    "Kr": ("a", "Kr", None, "Bideau-Mehu fit; spread 0.8 percent"),
    "Xe": ("a", "Xe", None, "Bideau-Mehu fit; spread 4.0 percent"),
    "N2O": ("a", "N2O", 288.15, "He 2021 Table 1; 18.3e-27 measured at 532 nm"),
    "NH3": ("a", "NH3, Cuthbertson", 273.15, "Cuthbertson fit (0.48 to 0.67 um); the generic formula gives 6.65e-27"),
    "H2O": ("a", "H2O vapour (generic check only)", None, "see the note on the water vapour row"),
    "H2S": ("b", "H2S", None, "polarizability route, each value about +-8 percent; Bogaard 1978 not reached"),
    "SO2": ("b", "SO2", None, "polarizability route, each value about +-8 percent; Bogaard 1978 not reached"),
    "C2H6": ("b", "C2H6", None, "polarizability route"),
    "C2H4": ("b", "C2H4", None, "polarizability route"),
    "HCN": ("b", "HCN", None, "polarizability route"),
    "C2H2": ("b", "C2H2", None, "polarizability route"),
    "OCS": ("b", "OCS", None, "polarizability route"),
    "HCl": ("b", "HCl", None, "polarizability route"),
    "CH2O": ("b", "H2CO", None, "polarizability route"),
    "O3": ("b", "O3, CCCBDB alpha with DFT rho", None,
           "polarizability route; the Ranjan and Sasselov values give 1.64e-26 at 550 nm: use 1.5e-26 to 1.75e-26"),
}

WATER_SIGMA_550_REPORT = 3.11e-27
"""float: n-based sigma of water vapor at 550 nm from n - 1 = 2.362e-4, if the reference temperature is 288.15 K (part 02)."""
WATER_SIGMA_550_REPORT_273 = 2.80e-27


def read_table(filename, header_start):
    """The rows (lists of stripped cells) of the first markdown table whose header line starts with `header_start`."""
    with open(os.path.join(PARTS, filename), encoding="utf-8") as handle:
        lines = handle.read().split("\n")
    for index, line in enumerate(lines):
        if line.startswith(header_start):
            break
    else:
        raise SystemExit(f"{filename}: no table starting {header_start!r}")
    rows = []
    for line in lines[index + 2:]:
        if not line.startswith("|"):
            break
        rows.append([cell.strip() for cell in line.strip().strip("|").split("|")])
    return rows


NUMBER = r"-?\d+\.?\d*(?:e[+-]?\d+)?"
"""str: A number in the tables, with an optional exponent."""


def tagged(cell, tags):
    """`(value, tag letter or None, sublimation flag)` of a cell like `13.957h`, `194.67y sub`; `(None, None, False)` for a dash."""
    if cell in ("-", "--", ""):
        return None, None, False
    match = re.fullmatch(rf"({NUMBER})([a-zA-Z])?( sub)?", cell)
    if not match:
        raise SystemExit(f"cannot read the tagged number {cell!r}")
    tag = match.group(2)
    if tag is not None and tag not in tags:
        raise SystemExit(f"unknown tag {tag!r} in {cell!r}")
    return float(match.group(1)), tag, bool(match.group(3))


def with_range(cell):
    """`(value, low, high, count text)` of `33.15 [33.15-33.20/5]`; all None for `- [-]`."""
    if cell.startswith("-"):
        return None, None, None, None
    match = re.fullmatch(rf"({NUMBER}) \[({NUMBER})-({NUMBER})/(\S+)\]", cell)
    if not match:
        raise SystemExit(f"cannot read the ranged number {cell!r}")
    return float(match.group(1)), float(match.group(2)), float(match.group(3)), match.group(4)


def span(cell):
    """A number, or a `lo-hi` pair as `[lo, hi]`, or None for a dash (kJ/mol columns)."""
    if cell in ("-", "--"):
        return None
    match = re.fullmatch(rf"({NUMBER})(?:[a-zA-Z])?", cell)
    if match:
        return float(match.group(1))
    match = re.fullmatch(rf"({NUMBER})-({NUMBER})", cell)
    if match:
        return [float(match.group(1)), float(match.group(2))]
    raise SystemExit(f"cannot read the number or span {cell!r}")


def dh(cell):
    """`(value or [lo, hi], tag)` of an enthalpy cell such as `40.65e` or `7.6-8.1`."""
    if cell in ("-", "--"):
        return None, None
    match = re.fullmatch(rf"({NUMBER})([a-zA-Z])", cell)
    if match:
        return float(match.group(1)), match.group(2)
    return span(cell), None


def build_identity(rows):
    out = {}
    for species, cas, mass, ir, vol, threshold, nuv, status in rows:
        match = re.fullmatch(r"(\d+) \((.+)\)", threshold)
        out[species] = {
            "cas": cas,
            "molar_mass_research_g_mol": float(mass),
            "ir": {"CIA": "cia", "IR": "ir", "IR(weak)": "weak", "-": "none"}[ir],
            "volatility": {"?": None}.get(vol, vol),
            "photolysis_threshold_nm": float(match.group(1)) if match else None,
            "photolysis_channel": match.group(2) if match else None,
            "absorbs_near_uv": nuv == "Y",
            "identity_status": status,
        }
    return out


def build_phase(rows):
    out = {}
    for species, tt, pt, tb, tc, pc, omega, status in rows:
        tt_v, tt_tag, _ = tagged(tt, PHASE_TAGS)
        pt_v, pt_tag, _ = tagged(pt, PHASE_TAGS)
        tb_v, tb_tag, tb_sub = tagged(tb, PHASE_TAGS)
        tc_v, tc_lo, tc_hi, tc_n = with_range(tc)
        pc_v, pc_lo, pc_hi, pc_n = with_range(pc)
        out[species] = {
            "triple_point_k": tt_v, "triple_point_k_source": PHASE_TAGS.get(tt_tag),
            "triple_point_pa": pt_v, "triple_point_pa_source": PHASE_TAGS.get(pt_tag),
            "boiling_point_k": tb_v, "boiling_point_k_source": PHASE_TAGS.get(tb_tag),
            "boiling_point_is_sublimation": tb_sub,
            "critical_temperature_k": tc_v, "critical_temperature_range_k": [tc_lo, tc_hi] if tc_v is not None else None,
            "critical_temperature_n": tc_n,
            "critical_pressure_pa": pc_v * 1e6 if pc_v is not None else None,
            "critical_pressure_range_pa": [pc_lo * 1e6, pc_hi * 1e6] if pc_v is not None else None,
            "critical_pressure_n": pc_n,
            "acentric_factor": None if omega == "-" else float(omega),
            "phase_status": status,
        }
    return out


def build_thermal(rows):
    out = {}
    for species, cp100, cp300, cp1000, gamma, spread, dhvap, dhsub, status in rows:
        cps = [tagged(cell, CP_TAGS) for cell in (cp100, cp300, cp1000)]
        vap, vap_tag = dh(dhvap)
        sub, sub_tag = dh(dhsub)
        out[species] = {
            "cp_j_mol_k": {"100": cps[0][0], "300": cps[1][0], "1000": cps[2][0]},
            "cp_source": {"100": CP_TAGS.get(cps[0][1]), "300": CP_TAGS.get(cps[1][1]), "1000": CP_TAGS.get(cps[2][1])},
            "gamma_300": None if gamma == "-" else float(gamma),
            "latent_heat_vaporisation_kj_mol": vap, "latent_heat_vaporisation_source": PHASE_TAGS.get(vap_tag),
            "latent_heat_sublimation_kj_mol": sub,
            "thermal_status": status,
        }
    return out


def build_extras(rows):
    out = {}
    for species, mass, cp100, cp300, cp1000, other, status in rows:
        out[species] = {
            "molar_mass_research_g_mol": float(mass),
            "cp_j_mol_k": {"100": None if cp100 == "-" else float(cp100), "300": None if cp300 == "-" else float(cp300),
                           "1000": None if cp1000 == "-" else float(cp1000)},
            "note": other, "thermal_status": status,
        }
    return out


def sigma_row(cells, kind):
    """The adopted Rayleigh numbers of a row of table 2.1a (`kind` a) or 2.1b (`kind` b)."""
    if kind == "a":
        label, s450, s550, s650, ratio, fk, source, status = cells
    else:
        label, alpha, rho, fk, s450, s550, s650, ratio, brange, status = cells
        source = f"alpha {alpha} A^3, rho_n {rho}; sigma(550) at B = 0.004 to 0.0075 um^2: {brange}"
    return {
        "sigma_cm2": {"450": float(s450), "550": float(s550), "650": float(s650)},
        "ratio_to_n2_550": float(ratio),
        "king_factor_550": float(re.match(NUMBER, fk).group(0)),
        "king_factor_assumed": "assumed" in fk,
        "research_status": status,
        "research_row": label,
        "research_note": source,
    }


def build_rayleigh(rows_a, rows_b):
    by_label = {"a": {r[0]: r for r in rows_a}, "b": {r[0]: r for r in rows_b}}
    out = {}
    for species, (kind, label, ref_k, note) in RAYLEIGH_ROWS.items():
        entry = sigma_row(by_label[kind][label], kind)
        entry["reference_temperature_k"] = ref_k
        entry["note"] = note
        if species == "H2O":
            # The refractivity n - 1 = 2.362e-4 at 550 nm has no known reference temperature, so the generic polarizability
            # route is used (GEN.410 merged research; part 02 also prints the n-based values, kept in the entry).
            entry["note"] = ("generic polarizability route (alpha 1.501 A^3), because the refractivity n - 1 = 2.362e-4 at 550 nm has "
                             "no known reference temperature: n-based sigma(550) is 3.11e-27 if 288.15 K or 2.80e-27 if 273.15 K, "
                             "so the generic value is 9 to 21 percent above; F_K 1.0010 from the vapor depolarization 2.99e-4 "
                             "(linear), not the liquid's 0.17")
            entry["sigma_550_n_based_cm2"] = {"288.15": WATER_SIGMA_550_REPORT, "273.15": WATER_SIGMA_550_REPORT_273}
        out[species] = entry
    return out


# Species whose solubility the research says grows with pH or dissociation (part 01 section 8); the Henry constant is the physical one only.
HENRY_DISSOCIATES = ("CO2", "SO2", "NH3", "H2S", "HCl", "HCN")

HENRY_SPREAD_FACTOR = 3.0
"""float: A compiled range wider than this factor (high over low) makes a Henry entry partly verified."""
HENRY_DISAGREEMENT_FACTOR = 2.0
"""float: A NIST median and an independent thermo fit more than this factor apart make a Henry entry partly verified."""
VISCOSITY_CHECK_PERCENT = 5.0
"""float: A Chapman-Enskog viscosity within this of the NIST fluid value at 300 K counts as checked (the research: within 2 percent for non-polar gases)."""


def build_henry(rows):
    out = {}
    for species, nist, d_k, thermo_kh, thermo_d, status in rows:
        match = re.fullmatch(rf"({NUMBER}) \[({NUMBER})-({NUMBER}), n=(\d+)\]", nist)
        if not match:
            raise SystemExit(f"cannot read the Henry cell {nist!r}")
        kh, low, high, count = float(match.group(1)), float(match.group(2)), float(match.group(3)), int(match.group(4))
        reasons = []
        if count == 1:
            reasons.append("one compiled value")
        if high / low > HENRY_SPREAD_FACTOR:
            reasons.append(f"the compiled values span a factor {high / low:.0f}")
        ratio = max(kh, float(thermo_kh)) / min(kh, float(thermo_kh))
        if ratio > HENRY_DISAGREEMENT_FACTOR:
            reasons.append(f"the NIST median and the thermo fit differ by a factor {ratio:.1f}")
        out[species] = {
            "kh_298_mol_kg_bar": kh, "kh_range_mol_kg_bar": [low, high], "kh_n": count, "d_k": float(d_k),
            "thermo_kh_298_mol_kg_bar": float(thermo_kh), "thermo_d_k": float(thermo_d),
            "dissociates": species in HENRY_DISSOCIATES,
            "research_status": status,
            "state": "partly verified" if reasons else "verified",
            "limit": "; ".join(reasons),
        }
    return out


def build_viscosity(rows):
    out = {}
    for species, sigma, epsilon, mu_300, nist_300, diff, status in rows:
        nist = None if nist_300 == "-" else float(nist_300)
        difference = None if diff == "-" else float(diff)
        if species == "H2O":
            state, limit = "partly verified", ("polar: Chapman-Enskog is 10 percent high and the research gives the NIST value at 300 K "
                                              "only, no temperature law")
        elif nist is None:
            state, limit = "partly verified", "Chapman-Enskog only, not compared with a NIST fluid value"
        elif abs(difference) > VISCOSITY_CHECK_PERCENT:
            state, limit = "partly verified", f"Chapman-Enskog is {difference:+.1f} percent from the NIST fluid value at 300 K"
        else:
            state, limit = "verified", ""
        out[species] = {
            "sigma_angstrom": float(sigma), "epsilon_over_k_k": float(epsilon), "mu_300_research_upa_s": float(mu_300),
            "nist_300_upa_s": nist, "difference_percent": difference, "polar": species == "H2O",
            "research_status": status, "state": state, "limit": limit,
        }
    return out


# ---- Infrared bands (part 03 table 2.1) ----

def _band(label, cm1=None, low=None, high=None, about=False):
    """One band: its label as the research writes it, a centre in cm-1 (`about` where the research writes "about") or a range."""
    return {"label": label, "wavenumber_cm1": cm1, "range_cm1": None if low is None else [low, high], "approximate": about}


# The rows of table 2.1: the row label of the research table, then the band centres typed from its "Fundamentals" cell (the check in
# `build_infrared` fails when a number is not in that cell), and the rank of the gas in the saturated-absorption ordering of Byrne and Goldblatt 2014b
# where the strength cell gives one. A gas with no band has an empty list and the research's reason in the row.
INFRARED = {
    "H2O": ("H2O", [_band("nu2 bend", 1594.59), _band("nu1", 3656.65), _band("nu3", 3755.79),
                    _band("pure rotation edge", low=0.0, high=600.0, about=True)], 11),
    "CO2": ("CO2", [_band("nu2", 667, about=True), _band("nu3", 2349, about=True), _band("2nu2/nu1 Fermi dyad", 1388, about=True),
                    _band("2nu2/nu1 Fermi dyad", 1285, about=True), _band("hot band", 960, about=True),
                    _band("hot band", 1063, about=True)], 16),
    "CH4": ("CH4", [_band("nu4", 1306.2), _band("nu2", 1533.6), _band("nu1", 2917), _band("nu3", 3018.9)], 27),
    "NH3": ("NH3", [_band("nu2 umbrella", 932.5), _band("nu2 umbrella", 968.3), _band("nu4", 1626), _band("nu1", 3336),
                    _band("nu3", 3444)], None),
    "H2S": ("H2S", [_band("nu2", 1182.7), _band("nu1", 2614.6), _band("nu3", 2626)], None),
    "HCN": ("HCN", [_band("nu2 bend", 711.98), _band("nu3 CN", 2096.85), _band("nu1 CH", 3311.47)], None),
    "C2H6": ("C2H6", [_band("nu12 CH3 rock", 821.6), _band("CH3 deform", low=1379, high=1469),
                      _band("CH stretch", low=2896, high=2985), _band("torsion", 289)], None),
    "C2H2": ("C2H2", [_band("nu5 CH bend", 730.3), _band("nu2 CC", 1973.8), _band("nu3", 3289)], 6),
    "CH2O": ("H2CO", [_band("nu2 CO", 1746.1), _band("nu3 CH2 scissor", 1500.1), _band("nu5", 1249.1), _band("nu6", 1167.3),
                      _band("nu1", 2782.5), _band("nu4", 2843.1)], None),
    "CH3Cl": ("CH3Cl", [_band("nu3 CCl", 732.1), _band("nu6 rock", 1017.3), _band("nu2", 1354.9), _band("nu5", 1452.1)], None),
    "C2H4": ("C2H4", [_band("nu7 wag", 949.3), _band("nu10 rock", 826.0), _band("nu12", 1443.5), _band("nu11", 2988.66)], None),
    "PH3": ("PH3", [_band("nu2", 992.1), _band("nu4", 1118.3), _band("nu1", 2322.9), _band("nu3", 2327.7)], None),
    "CH3OH": ("CH3OH", [_band("nu8 CO stretch", 1033), _band("nu6 OH bend", 1345), _band("CH3 modes", low=1455, high=1477),
                        _band("CH3 modes", low=2844, high=3000), _band("OH stretch", 3681)], 1),
    "O3": ("O3", [_band("nu3", 1042, about=True), _band("nu1", 1110, about=True), _band("nu2", 705, about=True)], None),
    "N2O": ("N2O", [_band("nu2", 589, about=True), _band("nu1", 1285, about=True), _band("nu3", 2224, about=True)], None),
    "SO2": ("SO2", [_band("nu2", 518, about=True), _band("nu1", 1151, about=True), _band("nu3", 1362, about=True)], 8),
    "CO": ("CO", [_band("fundamental", 2143, about=True)], None),
    "OCS": ("OCS, NO2, HCl, HF, HNO3", [_band("OCS", 520, about=True), _band("OCS", 859, about=True), _band("OCS", 2062, about=True)], None),
    "NO2": ("OCS, NO2, HCl, HF, HNO3", [_band("NO2", 750, about=True), _band("NO2", 1318, about=True), _band("NO2", 1618, about=True)], None),
    "HCl": ("OCS, NO2, HCl, HF, HNO3", [_band("HCl", 2886, about=True)], None),
    "HF": ("OCS, NO2, HCl, HF, HNO3", [_band("HF", 3962, about=True)], None),
    "HNO3": ("OCS, NO2, HCl, HF, HNO3", [], None),
    "SiO": ("SiO", [_band("fundamental", 1242, about=True)], None),
    **{f: ("N2, O2, H2, He, Ar, Ne, Kr, Xe", [], None) for f in ("N2", "O2", "H2", "He", "Ne", "Ar", "Kr", "Xe")},
    **{f: ("TiO, VO, Na, K, Fe, Mg", [], None) for f in ("TiO", "VO", "Na", "K", "Fe", "Mg")},
}


def build_infrared(rows):
    """`{formula: infrared entry}` from the rows of part 03 table 2.1 and the typed bands of `INFRARED`."""
    by_label = {cells[0]: cells for cells in rows}
    out = {}
    for formula, (label, bands, rank) in INFRARED.items():
        if label not in by_label:
            raise SystemExit(f"part 03 table 2.1 has no row {label!r}")
        _gas, fundamentals, strength, source, status = by_label[label]
        for band in bands:
            for value in (band["wavenumber_cm1"], *(band["range_cm1"] or ())):
                if value is not None and value != 0.0 and f"{value:g}" not in fundamentals.replace(",", ""):
                    raise SystemExit(f"{formula}: {value:g} cm-1 is not in the research row {label!r}: {fundamentals!r}")
        if rank is not None and f"rank {rank}" not in strength:
            raise SystemExit(f"{formula}: the research row {label!r} does not say 'rank {rank}'")
        recalled = "RECALLED" in status
        out[formula] = {
            "bands": bands, "research_row": label, "fundamentals_text": fundamentals, "strength_text": strength,
            "saturated_rank": rank, "no_band": not bands, "source": source, "research_status": status,
            "state": "partly verified" if recalled else "verified",
            "limit": f"research status \"{status}\": what is RECALLED was recalled, not read from a table (part 03 section 12)" if recalled else "",
        }
    return out


# ---- Ultraviolet and visible absorption (part 02 section 2.3) ----

UV_SPECIES = {"H2CO": "CH2O"}
"""dict: The research's name of a species and the formula of its row where they differ."""

UV_TABLE_HEADERS = ("| bin nm | O2 |", "| bin nm | NH3 |", "| bin nm | H2CO |")

UV_LIMITS = {
    "O2": ("band means understate the transmission of thick columns of the structured Schumann-Runge bands (175 to 200 nm): "
           "the research says the generator needs a k-distribution or effective cross section there"),
    "CH4": ("the visible bins rest on Karkoschka 1994, which the research did not reach; giant-planet methane bands are one to two "
            "decades stronger than these bin means and need a band model (part 03)"),
}
"""dict: The limits the research names for the bin tables of a species (the other species' cells are DATA, COMPUTED)."""


def _uv_cell(cell):
    """`(value or None, partial)` of a table cell: `-` is outside the data, `0` a zero in the file, `*` less than half the bin covered."""
    if cell == "-":
        return None, False
    partial = cell.endswith("*")
    return float(cell.rstrip("*")), partial


def build_uv():
    """`{formula: uv entry}` from the three bin tables and the coverage sentence of part 02 section 2.3."""
    with open(os.path.join(PARTS, PART_02), encoding="utf-8") as handle:
        text = handle.read()
    lines = text.split("\n")
    edges, columns = None, {}
    for header in UV_TABLE_HEADERS:
        rows = read_table(PART_02, header)
        names = [cell.strip() for cell in next(line for line in lines if line.startswith(header)).strip().strip("|").split("|")][1:]
        bins = []
        for cells in rows:
            low, high = cells[0].split("-")
            bins.append([float(low), float(high)])
        if edges is None:
            edges = bins
        elif bins != edges:
            raise SystemExit("the three ultraviolet tables do not share their bins")
        for index, name in enumerate(names, 1):
            columns[name] = [_uv_cell(cells[index]) for cells in rows]
    sentence = re.search(r"coverage in nm is (.*?)\. Unlisted origins are (.*?)\.", text, re.S)
    if not sentence:
        raise SystemExit("part 02 section 2.3 has no coverage sentence")
    default_origin = sentence.group(2)
    coverage = {}
    for item in sentence.group(1).split("; "):
        match = re.fullmatch(r"(\S+) ([\d.]+) to ([\d.]+)(?: \((.*)\))?", item.strip(), re.S)
        if not match:
            raise SystemExit(f"cannot read the coverage item {item!r}")
        coverage[match.group(1)] = (float(match.group(2)), float(match.group(3)), (match.group(4) or "").replace("**", ""))
    out = {}
    for name, cells in columns.items():
        if name not in coverage:
            raise SystemExit(f"part 02 section 2.3 gives no coverage for {name}")
        low, high, origin = coverage[name]
        for (value, _partial), (bin_low, bin_high) in zip(cells, edges):
            overlaps = bin_high > low and bin_low < high
            if (value is None) == overlaps:
                raise SystemExit(f"{name}: the bin {bin_low}-{bin_high} nm disagrees with the coverage {low} to {high} nm")
        limit = UV_LIMITS.get(name, "")
        out[UV_SPECIES.get(name, name)] = {
            "research_name": name, "coverage_nm": [low, high], "origin": origin or default_origin,
            "cross_section_cm2": [value for value, _partial in cells], "partial_bin": [partial for _value, partial in cells],
            "research_status": "DATA, COMPUTED", "state": "partly verified" if limit else "verified", "limit": limit,
        }
    return {"bins_nm": edges, "species": out}


def state_of(blocks):
    """
    The source state of a row from the status words of its blocks (`{block name: status text}`), in the vocabulary of GEN.472:
    `(state, limit)`. Every block read or computed is verified; a block not reached or second-hand makes the row partly verified
    and the limit names it; a row with no block reached at all is open.
    """
    weak = {name: ("not reached" if "NOT REACHED" in text else "second-hand, no measurement reached")
            for name, text in blocks.items() if "NOT REACHED" in text or "2ND-HAND" in text}
    if len(weak) == len(blocks) and all("NOT REACHED" in blocks[name] for name in blocks):
        return "open", "no block reached by the research"
    if weak:
        return "partly verified", "; ".join(f"{name} {why}" for name, why in weak.items())
    return "verified", ""


def rayleigh_state(entry):
    """The source state of a Rayleigh entry: `(state, limit)`."""
    status, label = entry["research_status"], entry["research_row"]
    if "2ND-HAND" in status:
        return "partly verified", "second-hand form, no measurement reached"
    if label.startswith("H2O"):
        return "partly verified", "reference temperature of the refractivity unknown: the generic value is 9 to 21 percent above the n-based 3.11e-27 or 2.80e-27 at 550 nm"
    if entry["king_factor_assumed"]:
        return "partly verified", "King factor assumed 1.00 (not found)"
    if entry["reference_temperature_k"] is None and "alpha" in entry["research_note"]:
        return "partly verified", "polarizability route good to -9 to +6 percent on 13 verified gases"
    return "verified", ""


def build():
    identity = build_identity(read_table(PART_01, "| Species | CAS |"))
    phase = build_phase(read_table(PART_01, "| Species | Tt K |"))
    thermal = build_thermal(read_table(PART_01, "| Species | Cp 100 | Cp 300 | Cp 1000 | gamma 300"))
    extras = build_extras(read_table(PART_01, "| Species | M g/mol |"))
    henry = build_henry(read_table(PART_01, "| Species | kH NIST median"))
    viscosity = build_viscosity(read_table(PART_01, "| Species | sigma A |"))
    for name, table in (("Henry", henry), ("viscosity", viscosity)):
        unknown = [f for f in table if f not in identity and f not in extras]
        if unknown:
            raise SystemExit(f"the {name} table names species that are not rows: {unknown}")
    rayleigh = build_rayleigh(read_table(PART_02, "| Gas | sigma 450 nm |"), read_table(PART_02, "| Gas | alpha (A^3) | rho_n | F_K"))
    infrared = build_infrared(read_table(PART_03, "| Gas | Fundamentals"))
    uv = build_uv()
    for name, table in (("infrared", infrared), ("ultraviolet", uv["species"])):
        unknown = [f for f in table if f not in identity and f not in extras]
        if unknown:
            raise SystemExit(f"the {name} table names species that are not rows: {unknown}")

    gases = []
    for formula in list(identity) + [f for f in extras if f not in identity]:
        row = {"formula": formula, "name": NAMES[formula]}
        if formula in identity:
            row.update(identity[formula])
        else:
            row.update({"cas": None, "molar_mass_research_g_mol": extras[formula]["molar_mass_research_g_mol"], "ir": None,
                        "volatility": None, "photolysis_threshold_nm": None, "photolysis_channel": None,
                        "absorbs_near_uv": None, "identity_status": "NOT REACHED"})
        if formula in ISOTOPES:
            row["isotopes"] = ISOTOPES[formula]
        row.update(phase.get(formula, {}))
        if formula in thermal:
            row.update(thermal[formula])
        if formula in extras:
            # Part 01 section 5 carries Cp of the species it adds; TiO is in T3 as well and keeps the T3 numbers.
            row.setdefault("cp_j_mol_k", extras[formula]["cp_j_mol_k"])
            row["extras_note"] = extras[formula]["note"]
            row.setdefault("thermal_status", extras[formula]["thermal_status"])
        if formula in rayleigh:
            row["rayleigh"] = rayleigh[formula]
            row["rayleigh"]["state"], row["rayleigh"]["limit"] = rayleigh_state(rayleigh[formula])
        if formula in henry:
            row["henry"] = henry[formula]
        if formula in viscosity:
            row["viscosity"] = viscosity[formula]
        if formula in infrared:
            row["infrared"] = infrared[formula]
        if formula in uv["species"]:
            row["uv"] = uv["species"][formula]
        blocks = {name: row[key] for name, key in (("identity", "identity_status"), ("phase points", "phase_status"),
                                                   ("heat capacity", "thermal_status")) if row.get(key)}
        state, limit = state_of(blocks)
        row["state"] = state
        row["limit"] = limit
        row["source"] = (f"docs/design/composition/atmosphere-composition-parts/{PART_01} T1, T2, T3 (section 5 for the species it adds)"
                         + (f" and {PART_02} table 2.1" if formula in rayleigh else "")
                         + (f"; {PART_01} section 8" if formula in henry else "")
                         + (f"; {PART_01} section 9" if formula in viscosity else "")
                         + (f"; {PART_03} table 2.1" if formula in infrared else "")
                         + (f"; {PART_02} section 2.3" if formula in uv["species"] else ""))
        gases.append(row)

    return {
        "retrieved": RETRIEVED,
        "title": "The gases and compounds that can be in an atmosphere: identity, phase points, heat capacity, Rayleigh scattering, Henry constants, viscosity, infrared bands, ultraviolet and visible absorption (GEN.410)",
        "sources": [
            "chemicals 1.5.2 (heos_constants, critical, triple and vapour-pressure tables), thermo 0.6.1, periodictable 2.1.0, "
            "cantera nasa_gas.yaml, NIST WebBook and NIST fluid saturation tables, Yaws and CRC values as tagged in part 01",
            "He et al. 2021, ACP 21, 14927 (N2O, CH4, CO2, O2 refractivity); Peck and Huang 1977; Mansfield and Peck; "
            "Bideau-Mehu et al.; Cuthbertson; Allen/Cox; NIST CCCBDB polarizabilities as read in part 02",
            "Shimanouchi 1972 gas-phase fundamentals at NIST and Byrne and Goldblatt 2014b (arXiv:1409.1880) as read in part 03 table 2.1; "
            "photochem_clima_data 0.3.4 (Heays 2017, Huebner 2015, Schmidt 2013, Karkoschka 1994, Keller-Rudek 2013) "
            "binned to 35 bins in part 02 section 2.3",
        ],
        "notes": ("Each block of a row has its own status word from the research part (DATA, COMPUTED, READ, 2ND-HAND, "
                  "NOT REACHED); `state` and `limit` are the row's source state in the vocabulary of GEN.472. The molar mass "
                  "in the file is the research part's, kept for the audit; the loader computes the working value from "
                  "`physics.elements` (GEN.438). Rayleigh sigma is per molecule in cm2; `reference_temperature_k` is the "
                  "temperature the refractivity fit was made at (GEN.408). Henry constants are mol/(kg bar) at 298.15 K "
                  "with d = -dln kH/d(1/T) in K; viscosity needs the Lennard-Jones sigma (Angstrom) and epsilon/k (K). `infrared` lists "
                  "the band centres (cm-1) the research gives, `approximate` where it writes 'about'; `uv` holds each gas's mean "
                  "photoabsorption cross section (cm2 per molecule) in the bins of `uv_bins_nm`, `null` outside the data and `partial_bin` "
                  "where under half of the bin is covered."),
        "uv_bins_nm": uv["bins_nm"],
        "gases": gases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="fail if the committed file differs from a fresh build")
    args = parser.parse_args()
    text = json.dumps(build(), indent=1, ensure_ascii=False) + "\n"
    if args.check:
        with open(OUT, encoding="utf-8") as handle:
            if handle.read() != text:
                print("gases.json differs from a fresh build; run scripts/build_gas_table.py", file=sys.stderr)
                return 1
        return 0
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(text)
    print(f"wrote {os.path.relpath(OUT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
