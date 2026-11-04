#!/usr/bin/env python3
"""
Builds the shared materials table (GEN.427) from the research tables it is typed in from and the values read on 2026-10-22.

    python scripts/build_materials_table.py            # write src/planetgen/physics/data/materials.json
    python scripts/build_materials_table.py --check    # build in memory and fail if the committed file differs

The file is `src/planetgen/physics/data/materials.json`. The base is part E's materials table (section 2.10 of
`docs/design/composition/body-composition-parts/E-ground-radiogenic-shielding.md`, 44 rows, label column kept); part B's [R4233] phases
(table 2.1a) and iron rows (table 2.1b) and part C's ice phases (section 2.3) add the equation-of-state rows; the values in
`materials_data.py` ([R3992], the Handbook of Mineralogy, PubChem) replace every value the parts labelled
RECALLED and fill the minerals the shared table needs (every name of `ASTEROID_COMPONENTS`, `COMET_COMPOSITION`, the field
families and the class compositions of `tuning`). A value with no source is not kept: its cell is `null` with the code `N`.
Nothing is fitted or invented here.
"""

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "src"))

import materials_data as data  # noqa: E402
from planetgen.physics import elements  # noqa: E402

PARTS = os.path.join(HERE, "..", "docs", "design", "composition", "body-composition-parts")
OUT = os.path.join(HERE, "..", "src", "planetgen", "physics", "data", "materials.json")
PART_E = os.path.join(PARTS, "E-ground-radiogenic-shielding.md")
PART_B = os.path.join(PARTS, "B-rocky.md")

RETRIEVED = "2026-10-22"

ANY = "any"
PLANETARY = "planetary"

SOURCES = {
    "research-e": "Part E of the body composition research, section 2.10 (docs/design/composition/body-composition-parts/E-ground-radiogenic-shielding.md): the "
                  "values its label column marks R (read), H (second-hand) or C (computed); a value marked Q (recalled) is not kept",
    "nist-webbook": "NIST Chemistry WebBook phase-change and critical data, as read in part E (label R)",
    "crc-chemicals": "CRC Handbook tables as carried in the chemicals 1.5.2 package, as read in part E (label H, second-hand)",
    "pdg-nuclear": "Particle Data Group atomic and nuclear properties, nuclear interaction length, as read in part E (label R)",
    "slb2011-burnman": "Stixrude and Lithgow-Bertelloni mineral parameters read in the burnman package, as read in part E (label H)",
    "e-table-2.6": "Band reflectance of part E table 2.6, computed from the spectral library values it lists (label C)",
    "computed-e": "Computed in part E from its other columns (molar mass, latent heat per kilogram, interaction length from oxides)",
    "handbook-mineralogy": "Anthony, Bideaux, Bladh and Nichols, Handbook of Mineralogy, Mineralogical Society of America (Mineral Data "
                           "Publishing, 2001-2005), the page of the mineral: D(calc.) and D(meas.)",
    "webmineral": "Webmineral Mineralogy Database (Barthelmy), the data page of the mineral: the density listed there",
    "crc-density-wiki": "Densities of the elements (data page), Wikipedia, the row marked 'use' (CRC and WebElements values); the handbook itself was not read",
    "pubchem": "PubChem compound pages, section Density (experimental values from the sources named there: CRC, USCG, NTP, EPA, HSDB)",
    "seager2007": "[R3992], Mass-radius relationships for solid exoplanets, ApJ 669, 1279 "
                  "[R3992: Tables 1], 2 and 3 (read)",
    "hakim2018": "[R1752] (iron and light-element equations of state), Tables D.5 and F.7, as read in part B",
    "part-b-fit": "Part B's single Vinet fit to hcp iron (this work, computed, docs/design/composition/body-composition-parts/B-rocky.md section 2.1b)",
    "slb2024-burnman": "[R4233], Thermodynamics of mantle minerals III, GJI 237, 1699 "
                       "[R4233], as parameters read in burnman 2.1.0 SLB_2024.py (second-hand; part B table 2.1a)",
    "exoplex": "ExoPlex 1.0.1a1 source distribution ([R0111] liquid iron; water BM4 constants), second-hand (part B table 2.1b and 2.1c)",
    "hp2011-burnman": "[R1981] data set as read in burnman HP_2011_ds62 (diamond), second-hand (part B table 2.1b)",
    "journaux2020": "[R2204], JGR Planets, Gibbs-energy splines as distributed in SeaFreeze 1.1.3 (spline data only), evaluated in part C "
                    "section 2.3 (computed from second-hand data)",
    "blum-bischoff": "Blum and Bischoff, [R0430], as read in part C section 2.3 (comet densities)",
    "cod-part-d": "Part D of the body composition research, section 2.1: densities computed from the Crystallography Open Database cell of the same formula "
                  "(rho = Z M / (N_A V)), label COMPUTED; not a measurement of the solid at the body's temperature",
    "ideal-gas": "The ideal gas law rho = P M / (R T); exact, no constants of its own",
}

# A property: key -> (unit, the column key of the part E label string)
PROPERTIES = {
    "rho_grain_g_cm3": ("g/cm3", "rho"),
    "t_melt_k": ("K", "Tm"),
    "t_boil_or_sub_k": ("K", "Tb"),
    "t_triple_k": ("K", "Tt"),
    "p_triple_bar": ("bar", "Pt"),
    "t_crit_k": ("K", "Tc"),
    "p_crit_bar": ("bar", "Pc"),
    "h_fus_kj_mol": ("kJ/mol", "Hfus"),
    "h_vap_kj_mol": ("kJ/mol", "Hvap"),
    "h_sub_kj_mol": ("kJ/mol", "Hsub"),
    "h_fus_kj_kg": ("kJ/kg", "Hfus_kg"),
    "lambda_i_g_cm2": ("g/cm2", "lamI"),
    "debye_k": ("K", "debye"),
    "alb_040_050": ("1", "alb"),
    "alb_050_060": ("1", "alb"),
    "alb_060_070": ("1", "alb"),
    "rho_gas_stp_g_l": ("g/L", "rhog"),
}

CSV_COLUMNS = {  # part E csv header -> property key
    "rho_grain_g_cm3": "rho_grain_g_cm3", "T_melt_K": "t_melt_k", "T_boil_or_sub_K": "t_boil_or_sub_k", "T_triple_K": "t_triple_k",
    "P_triple_bar": "p_triple_bar", "T_crit_K": "t_crit_k", "P_crit_bar": "p_crit_bar", "H_fus_kJ_mol": "h_fus_kj_mol",
    "H_vap_kJ_mol": "h_vap_kj_mol", "H_sub_kJ_mol": "h_sub_kj_mol", "H_fus_kJ_kg": "h_fus_kj_kg", "lambda_I_g_cm2": "lambda_i_g_cm2",
    "debye_K": "debye_k", "alb_040_050": "alb_040_050", "alb_050_060": "alb_050_060", "alb_060_070": "alb_060_070",
}
LABEL_TO_PROPS = {  # label key of part E -> property keys
    "rho": ["rho_grain_g_cm3"], "Tm": ["t_melt_k"], "Tb": ["t_boil_or_sub_k"], "Tt": ["t_triple_k"], "Pt": ["p_triple_bar"],
    "Tc": ["t_crit_k"], "Pc": ["p_crit_bar"], "Hfus": ["h_fus_kj_mol"], "Hvap": ["h_vap_kj_mol"], "Hsub": ["h_sub_kj_mol"],
    "Hfus_kg": ["h_fus_kj_kg"], "lamI": ["lambda_i_g_cm2"], "debye": ["debye_k"],
    "alb": ["alb_040_050", "alb_050_060", "alb_060_070"], "rhog": ["rho_gas_stp_g_l"],
}
E_SOURCE = {  # (label key, code) -> source key
    ("Tm", "R"): "nist-webbook", ("Tb", "R"): "nist-webbook", ("Tt", "R"): "nist-webbook", ("Pt", "R"): "nist-webbook",
    ("Tc", "R"): "nist-webbook", ("Pc", "R"): "nist-webbook", ("Hvap", "R"): "nist-webbook", ("Hsub", "R"): "nist-webbook",
    ("Tm", "H"): "crc-chemicals", ("Tb", "H"): "crc-chemicals", ("Tt", "H"): "crc-chemicals", ("Tc", "H"): "crc-chemicals",
    ("Hfus", "H"): "crc-chemicals", ("Hvap", "H"): "crc-chemicals", ("Hsub", "H"): "crc-chemicals",
    ("lamI", "R"): "pdg-nuclear", ("lamI", "C"): "computed-e", ("rho", "R"): "research-e", ("rho", "H"): "slb2011-burnman",
    ("debye", "H"): "slb2011-burnman", ("alb", "C"): "e-table-2.6", ("Hfus_kg", "C"): "computed-e",
}
CODE_STATE = {"R": "verified", "C": "verified", "H": "partly verified", "S": "partly verified", "N": "open"}
HYDROGEN_FROM = "gas_H2"
CO_ICE_FROM = "gas_CO"

DENSITY_LIMIT = {
    "handbook-mineralogy": "a compilation read, not the original measurement",
    "webmineral": "a web compilation: a listed average or single value",
    "crc-density-wiki": "the CRC values as listed on a Wikipedia data page; the handbook was not read",
    "pubchem": "an annotation quoted with its source; the temperature is the annotation's",
}


def read(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


def parse_part_e():
    text = read(PART_E)
    block = re.search(r"```csv\n(id,family.*?)\n```", text, re.S).group(1)
    lines = block.split("\n")
    header = lines[0].split(",")
    rows = []
    for line in lines[1:]:
        cells = line.split(",")
        if len(cells) != len(header):
            raise SystemExit(f"part E csv row {cells[0]!r}: {len(cells)} cells for {len(header)} columns")
        rows.append(dict(zip(header, cells)))
    return rows


def parse_slb():
    text = read(PART_B)
    start = text.index("**Table 2.1a")
    section = text[start:text.index("**Table 2.1b")]
    out = []
    for line in section.split("\n"):
        if not line.startswith("|") or line.startswith("|---") or line.startswith("| Phase"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        phase, formula, mass, v0, rho0, k0, k0p, theta, gamma, q, stable = cells[:11]
        out.append({"phase": phase, "formula": formula, "m": float(mass), "v0": float(v0), "rho0": float(rho0), "k0": float(k0),
                    "k0p": float(k0p), "theta": float(theta), "gamma": float(gamma), "q": float(q), "stable": stable})
    return out


def prop(value, code, source, unit, **extra):
    entry = {"value": value, "unit": unit, "source": source, "code": code}
    entry.update({k: v for k, v in extra.items() if v is not None})
    return entry


def null_prop(key):
    return {"value": None, "unit": PROPERTIES[key][0], "source": None, "code": "N"}


def new_row(row_id, name, family, formula=None, aliases=(), gate=ANY, gas=None, members=None, note="", composition=""):
    return {"id": row_id, "name": name, "family": family, "formula": formula, "composition": composition,
            "aliases": list(aliases), "members": list(members) if members else None, "gas_formula": gas, "gate": gate,
            "props": {}, "eos": [], "labels": "", "note": note}


def number(cell):
    return float(cell) if cell.strip() else None


def e_rows():
    """Part E's solid rows with their labels; the gas rows are pointers and only supply hydrogen and carbon monoxide."""
    result, gas_rows = {}, {}
    for item in parse_part_e():
        labels = {}
        for pair in item["labels"].split(";"):
            if pair:
                key, _, code = pair.partition("=")
                labels[key] = code
        row = new_row(item["id"], item["id"].replace("_", " "), item["family"], note="")
        row["composition"] = item["composition"]
        if "all" in labels:
            row["props"] = {}
            result[item["id"]] = (row, item, labels)
            continue
        dropped = []
        for header, key in CSV_COLUMNS.items():
            label_key = PROPERTIES[key][1]
            code = labels.get(label_key)
            cell = item[header].strip()
            if code is None:
                continue
            if cell == "" or code in ("N", "Q"):
                if cell and code == "Q":
                    dropped.append(f"{key} {cell} (recalled, no source)")
                row["props"][key] = null_prop(key)
                continue
            row["props"][key] = prop(float(cell), code, E_SOURCE[(label_key, code)], PROPERTIES[key][0])
        if dropped:
            row["note"] = "Dropped as recalled with no source: " + "; ".join(dropped) + "."
        target = gas_rows if item["owner"] == "atmosphere" else result
        target[item["id"]] = (row, item, labels)
    return result, gas_rows


def formula_for(item):
    text = item["composition"]
    try:
        elements.parse_formula(text)
        return text
    except (ValueError, KeyError):
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    solid, gas_rows = e_rows()
    rows = {}
    for row_id, (row, item, labels) in solid.items():
        row["formula"] = formula_for(item)
        rows[row_id] = row

    # --- aliases, family fixes and the densities that replace recalled ones on part E's rows
    alias = {
        "water_ice": ["water ice", "ice"], "co2_ice": ["carbon dioxide ice"], "ch4_ice": ["methane ice", "methane"],
        "nh3_ice": ["ammonia ice", "ammonia"], "basalt_powder": ["basalt"], "serpentine": [], "iron": [], "quartz": ["silicon dioxide", "igneous silica"],
    }
    for row_id, names in alias.items():
        rows[row_id]["aliases"] = list(names)

    def set_density(row_id, value, source, code, lo=None, hi=None, basis=""):
        extra = {}
        if lo is not None and hi is not None:
            extra["measured_range"] = [lo, hi]
        old = rows[row_id]["props"].get("rho_grain_g_cm3")
        was = old["value"] if old and old.get("value") is not None else None
        rows[row_id]["props"]["rho_grain_g_cm3"] = prop(value, code, source, "g/cm3", basis=basis or None, **extra)
        if was is not None and abs(was - value) / value > 0.005:
            rows[row_id]["note"] = (rows[row_id].get("note", "") + f" Part E's {was} g/cm3 is replaced by the read value ({source}).").strip()

    # Handbook values on part E's mineral rows (verified wins over a second-hand parameter set).
    for row_id, handbook_id in (("halite", "halite"), ("gypsum", "gypsum"), ("calcite", "calcite"), ("kaolinite", "kaolinite"),
                                ("serpentine", "serpentine"), ("goethite", "goethite"), ("ilmenite", "ilmenite"),
                                ("hematite_fine", "hematite"), ("hematite_lt10um", "hematite"), ("forsterite", "forsterite"),
                                ("fayalite", "fayalite"), ("enstatite", "enstatite"), ("ferrosilite", "ferrosilite"),
                                ("diopside", "diopside"), ("albite", "albite"), ("quartz", "quartz")):
        value, lo, hi, source, basis = data.MINERAL_DENSITY[handbook_id]
        set_density(row_id, value, source, "R" if source == "handbook-mineralogy" else "H", lo, hi, basis)
    set_density("sulfur", 2.076, "handbook-mineralogy", "R", 2.07, 2.07, "D(calc.) of native sulfur (the page is filed as sulphur)")
    for row_id, (value, temperature) in (("co2_ice", data.PUBCHEM_DENSITY["co2_ice"]), ("n2_ice", data.PUBCHEM_DENSITY["n2_ice"])):
        rows[row_id]["props"]["rho_grain_g_cm3"] = prop(value, "R", "pubchem", "g/cm3", basis=temperature)

    # --- the pure elements the lists name
    element_rows = [
        ("carbon_graphite", "graphite", "C", ["carbon"]), ("carbon_diamond", "diamond", "C", []), ("silicon", None, "Si", []),
        ("magnesium", None, "Mg", []), ("aluminum", None, "Al", []), ("calcium", None, "Ca", []), ("phosphorus", None, "P", []),
        ("nickel", None, "Ni", []), ("iridium", None, "Ir", []), ("palladium", None, "Pd", []), ("platinum", None, "Pt", []),
        ("gold", None, "Au", []), ("osmium", None, "Os", []), ("ruthenium", None, "Ru", []), ("rhodium", None, "Rh", []),
        ("titanium", None, "Ti", []), ("potassium", None, "K", []), ("tungsten", None, "W", []),
    ]
    for row_id, name, formula, names in element_rows:
        row = new_row(row_id, name or row_id, "element", formula, names)
        row["props"]["rho_grain_g_cm3"] = prop(data.ELEMENT_DENSITY[row_id], "H", "crc-density-wiki", "g/cm3")
        rows[row_id] = row
    rows["phosphorus"]["note"] = "Red phosphorus; white is 1.823 and black 2.69 g/cm3 on the same page."
    rows["sulfur"]["aliases"] = ["sulfur"]
    rows["iron"]["aliases"] = ["iron"]

    # --- minerals the lists name, with Handbook or Webmineral densities
    minerals = [
        # id, name, family, formula, aliases, composition note
        ("troilite", "FeS", "sulfide", [], ""), ("schreibersite", None, "phosphide", [], "(Fe,Ni)3P"),
        ("cohenite", "Fe3C", "carbide", [], "(Fe,Ni,Co)3C; the formula is the iron end member"),
        ("magnetite", "Fe3O4", "oxide", [], ""), ("hematite", "Fe2O3", "oxide", [], ""), ("chromite", "FeCr2O4", "oxide", [], ""),
        ("moissanite", "SiC", "carbide", ["silicon carbide"], ""),
        ("bastnasite", "CeCO3F", "rare-earth", [], "bastnasite-(Ce)"), ("monazite", "CePO4", "rare-earth", [], "monazite-(Ce); (Ce,La,Nd,Th)PO4"),
        ("xenotime", "YPO4", "rare-earth", [], "xenotime-(Y)"), ("cerite", None, "rare-earth", [], "cerite-(Ce), a hydrous cerium silicate"),
        ("gadolinite", None, "rare-earth", [], "gadolinite-(Ce), (Ce,La,Nd,Y)2FeBe2Si2O10"),
        ("samarskite", None, "rare-earth", [], "samarskite-(Y), (Y,Fe,U)(Nb,Ta)5O4"), ("fergusonite", "YNbO4", "rare-earth", [], "fergusonite-(Y)"),
        ("euxenite", None, "rare-earth", [], "euxenite-(Y), (Y,Ca,Ce)(Nb,Ta,Ti)2O6"),
        ("osmiridium", None, "platinum-group", [], "an osmium-iridium alloy; the density is that of iridosmine (Os,Ir)"),
        ("sperrylite", "PtAs2", "platinum-group", [], ""), ("cooperite", "PtS", "platinum-group", [], "(Pt,Pd)S"),
        ("braggite", None, "platinum-group", [], "(Pt,Pd,Ni)S"), ("laurite", "RuS2", "platinum-group", [], ""),
        ("vysotskite", None, "platinum-group", [], "(Pd,Ni)S"),
        ("uraninite", "UO2", "radioactive", [], ""), ("thorianite", "ThO2", "radioactive", [], ""),
        ("carnotite", "K2(UO2)2(VO4)2.3H2O", "radioactive", [], ""), ("autunite", "Ca(UO2)2(PO4)2.11H2O", "radioactive", [], "10 to 12 H2O"),
        ("brannerite", None, "radioactive", [], "(U,Ca,Y,Ce)(Ti,Fe)2O6"), ("torbernite", "Cu(UO2)2(PO4)2.12H2O", "radioactive", [], "8 to 12 H2O"),
        ("coffinite", "USiO4", "radioactive", [], "U(SiO4)1-x(OH)4x"),
        ("rutile", "TiO2", "oxide", [], ""), ("titanite", "CaTiSiO5", "silicate", [], ""), ("perovskite", "CaTiO3", "oxide", [],
         "the calcium titanate mineral; not the magnesium silicate perovskite (bridgmanite)"),
        ("brookite", "TiO2", "oxide", [], ""), ("anatase", "TiO2", "oxide", [], ""),
        ("corundum", "Al2O3", "oxide", ["aluminum oxide"], ""), ("beryl", "Be3Al2Si6O18", "silicate", ["beryllium aluminum cyclosilicate"], ""),
        ("topaz", "Al2SiO4F2", "silicate", ["aluminum silicate fluoride hydroxide"], "Al2SiO4(F,OH)2; the formula is the fluorine end member"),
        ("spinel", "MgAl2O4", "oxide", ["magnesium aluminum oxide"], ""), ("zircon", "ZrSiO4", "silicate", ["zirconium silicate"], ""),
        ("schorl", "NaFe3Al6(BO3)3Si6O18(OH)4", "silicate", ["borosilicate"], "a tourmaline"),
        ("orthoclase", "KAlSi3O8", "silicate", ["potassium aluminum silicate"], ""),
        ("fluorite", "CaF2", "halide", ["calcium fluoride"], ""), ("fluorapatite", "Ca5(PO4)3F", "phosphate", ["calcium fluorophosphate"], ""),
        ("turquoise", "CuAl6(PO4)4(OH)8.4H2O", "phosphate", ["hydrated copper aluminum phosphate"], ""),
        ("malachite", "Cu2(CO3)(OH)2", "carbonate", ["copper carbonate hydroxide"], ""),
        ("almandine", "Fe3Al2(SiO4)3", "silicate", [], ""), ("pyrope", "Mg3Al2(SiO4)3", "silicate", [], ""),
        ("spessartine", "Mn3Al2(SiO4)3", "silicate", [], ""), ("grossular", "Ca3Al2(SiO4)3", "silicate", [], ""),
        ("muscovite", "KAl3Si3O10(OH)2", "silicate", [], ""), ("biotite", None, "silicate", [], "K(Mg,Fe)3(Al,Fe)Si3O10(OH,F)2"),
        ("chrysoberyl", "BeAl2O4", "oxide", [], ""),
        ("kamacite", None, "metal", [], "alpha iron-nickel, about 5 to 7 percent nickel by mass in the Handbook's analyses"),
        ("taenite", None, "metal", [], "gamma iron-nickel, about 25 to 31 percent nickel by mass in the Handbook's analyses"),
        ("periclase", "MgO", "oxide", [], ""), ("wustite", "FeO", "oxide", [], "Fe(1-x)O"), ("anorthite", "CaAl2Si2O8", "silicate", [], ""),
    ]
    for row_id, formula, family, aliases_here, composition in minerals:
        row = new_row(row_id, row_id, family, formula, aliases_here, composition=composition)
        value, lo, hi, source, basis = data.MINERAL_DENSITY[row_id]
        extra = {"measured_range": [lo, hi]} if lo is not None and hi is not None else {}
        code = "R" if source == "handbook-mineralogy" else "H"
        row["props"]["rho_grain_g_cm3"] = prop(value, code, source, "g/cm3", basis=basis, **extra)
        rows[row_id] = row
    rows["calcite"]["aliases"] = ["calcium carbonate"]
    rows["halite"]["aliases"] = ["sodium chloride"]
    # tourmaline's gem name and the mica gem names are not separate rows

    # --- the organics and the ices the lists name
    for row_id, formula, composition in (("naphthalene", "C10H8", ""), ("anthracene", "C14H10", ""), ("phenanthrene", "C14H10", ""),
                                        ("pyrene", "C16H10", ""), ("coronene", "C24H12", ""), ("fluoranthene", "C16H10", "")):
        row = new_row(row_id, row_id, "polycyclic aromatic hydrocarbon", formula)
        if row_id in data.PUBCHEM_DENSITY:
            value, temperature = data.PUBCHEM_DENSITY[row_id]
            row["props"]["rho_grain_g_cm3"] = prop(value, "R", "pubchem", "g/cm3", basis=temperature)
        else:
            row["props"]["rho_grain_g_cm3"] = null_prop("rho_grain_g_cm3")
            row["note"] = "No density reached (PubChem has no record of the compound's density)."
        rows[row_id] = row
    gas_h2, h2_item, h2_labels = gas_rows[HYDROGEN_FROM]
    gas_co, co_item, co_labels = gas_rows[CO_ICE_FROM]
    co = new_row("co_ice", "carbon monoxide ice", "ice", "CO", ["carbon monoxide ice"], gas="CO")
    co["props"] = {k: v for k, v in gas_co["props"].items() if k not in ("lambda_i_g_cm2", "h_sub_kj_mol", "debye_k")
                   and not k.startswith("alb_")}
    co["props"]["rho_grain_g_cm3"] = prop(data.PUBCHEM_DENSITY["co_ice"][0], "R", "pubchem", "g/cm3", basis=data.PUBCHEM_DENSITY["co_ice"][1])
    co["note"] = "Phase points are part E's carbon monoxide gas row; the solid density is PubChem's hexagonal solid at 65 K."
    rows["co_ice"] = co
    for row_id, gas_formula in (("ch4_ice", "CH4"), ("nh3_ice", "NH3"), ("n2_ice", "N2"), ("so2_ice", "SO2"), ("co2_ice", "CO2")):
        rows[row_id]["gas_formula"] = gas_formula
    for row_id in ("ch4_ice", "nh3_ice", "so2_ice"):
        rows[row_id]["props"]["rho_grain_g_cm3"] = null_prop("rho_grain_g_cm3")
        rows[row_id]["note"] = ("No primary solid density was reached (part C: PubChem and the handbooks list the liquid only; the "
                                "crystallographic estimates are recalled and are not kept).")
    hcn = new_row("hcn_ice", "hydrogen cyanide ice", "ice", "HCN", ["hydrogen cyanide"], gas="HCN")
    hcn["props"]["rho_grain_g_cm3"] = null_prop("rho_grain_g_cm3")
    hcn["note"] = "PubChem lists the liquid (0.6875 g/cm3 at 20 C) only; no solid density was reached."
    rows["hcn_ice"] = hcn
    fa = new_row("formaldehyde", "formaldehyde", "organic ice", "CH2O", ["formaldehyde"], gas="CH2O")
    fa["props"]["rho_grain_g_cm3"] = null_prop("rho_grain_g_cm3")
    fa["note"] = "PubChem lists the anhydrous liquid (0.815 g/cm3 at -20 C) only; no solid density was reached."
    rows["formaldehyde"] = fa
    for row_id, name, formula, aliases_here in (("ch3oh_ice", "methanol ice", "CH3OH", ["methanol ice", "methanol"]),
                                                 ("c2h6_ice", "ethane ice", "C2H6", ["ethane ice", "ethane"]),
                                                 ("h2s_ice", "hydrogen sulfide ice", "H2S", ["hydrogen sulfide ice", "hydrogen sulfide"])):
        # The comet ices the volatile population (GEN.495) names that the lists did not: the gas table owns the molecule's gas figures.
        row = new_row(row_id, name, "ice", formula, aliases_here, gas=formula)
        if row_id in data.COD_DENSITY:
            value, basis = data.COD_DENSITY[row_id]
            row["props"]["rho_grain_g_cm3"] = prop(value, "C", "cod-part-d", "g/cm3", basis=basis)
        else:
            row["props"]["rho_grain_g_cm3"] = null_prop("rho_grain_g_cm3")
            row["note"] = "No solid density was reached (part D lists no cell for the solid)."
        rows[row_id] = row
    ac = new_row("amorphous_carbon", "amorphous carbon", "carbon", "C", ["amorphous carbon"])
    ac["props"]["rho_grain_g_cm3"] = null_prop("rho_grain_g_cm3")
    ac["note"] = "No density source reached; the value depends on the hydrogen content and the formation."
    rows["amorphous_carbon"] = ac
    co_ = new_row("complex_organics", "complex organic compounds", "organic", None, ["complex organic compounds"],
                  composition="refractory organics of comets and primitive bodies",
                  note="[R0776] (read, section 3.3, the CONSERT permittivity review) gives only a rough "
                       "'some organics (permittivity about 2, rho about 2)' with no unit written (g/cm3 by the ice value beside it) and no "
                       "organics density of its own elsewhere; Blum and Bischoff's 1.0 to 2.0 covers it, so the row keeps the middle, 1.5, "
                       "inside the widest range (two sources in question: GEN.660). The comet draw reads this row and holds no density of its "
                       "own for the organics (GEN.712).")
    co_["props"]["rho_grain_g_cm3"] = prop(1.5, "R", "blum-bischoff", "g/cm3", measured_range=[1.0, 2.0], basis="organics 1000 to 2000 kg/m3: the middle")
    rows["complex_organics"] = co_
    anth = new_row("anthracite", "anthracite", "carbon", None, ["anthracite"], composition="high-rank coal")
    anth["props"]["rho_grain_g_cm3"] = null_prop("rho_grain_g_cm3")
    anth["note"] = "No density source reached."
    rows["anthracite"] = anth

    # --- the groups: a name the lists use for several rows; no number of its own
    groups = [
        ("olivine", "olivine", "mineral group", ["forsterite", "fayalite"], ["olivine"]),
        ("pyroxene", "pyroxene", "mineral group", ["enstatite", "ferrosilite", "diopside"], ["pyroxene"]),
        ("plagioclase", "plagioclase feldspars", "mineral group", ["albite", "anorthite"], ["plagioclase feldspars"]),
        ("silicate_dust", "silicate dust", "dust", ["forsterite", "fayalite", "enstatite", "ferrosilite"], ["silicate dust"]),
        ("hydrocarbons", "hydrocarbons", "organic group", ["naphthalene", "anthracene", "phenanthrene", "pyrene", "coronene", "fluoranthene",
                                                            "complex_organics"], ["hydrocarbons"]),
        ("frozen_hydrocarbons", "frozen hydrocarbons", "organic ice group", ["ch4_ice", "complex_organics"], ["frozen hydrocarbons"]),
        ("rock", "rock", "rock group", ["basalt_powder", "felsic_granite_powder_a", "ultramafic_dunite_powder", "anorthosite_powder",
                                         "forsterite", "enstatite"], ["rock", "silicate rock"]),
    ]
    for row_id, name, family, members, names in groups:
        rows[row_id] = new_row(row_id, name, family, None, names, members=members)

    # --- hydrogen and helium
    hydrogen = new_row("hydrogen", "hydrogen", "envelope", "H2", ["hydrogen"], gate=PLANETARY, gas="H2")
    hydrogen["props"] = {k: v for k, v in gas_h2["props"].items() if k not in ("lambda_i_g_cm2", "h_sub_kj_mol", "debye_k", "rho_grain_g_cm3")
                         and not k.startswith("alb_") and k != "h_fus_kj_kg"}
    hydrogen["props"]["rho_gas_stp_g_l"] = prop(0.08988, "H", "crc-density-wiki", "g/L", basis="0 C, 1 atm")
    helium = new_row("helium", "helium", "envelope", "He", ["helium"], gate=PLANETARY, gas="He")
    helium["props"]["rho_gas_stp_g_l"] = prop(0.1786, "H", "crc-density-wiki", "g/L", basis="0 C, 1 atm")
    for row in (hydrogen, helium):
        row["eos"].append({"model": "ideal-gas", "params": {}, "valid": {}, "source": "ideal-gas", "code": "C",
                           "note": "rho = P M / (R T) with the molar mass of the row; the envelope radius comes from the "
                                   "[R1429] tables (part C), which GEN.428 reads, not from a constants set; no hydrogen-helium EOS "
                                   "constants are held here"})
    rows["hydrogen"], rows["helium"] = hydrogen, helium

    # --- the equation-of-state rows
    def eos_entry(model, params, source, code, valid=None, uncertainty=None, note=""):
        entry = {"model": model, "params": params, "valid": valid or {}, "source": source, "code": code, "note": note}
        if uncertainty:
            entry["uncertainty"] = uncertainty
        return entry

    iron_epsilon = new_row("iron_epsilon", "epsilon iron (hcp)", "metal", "Fe", [], note="The high-pressure iron phase of cores; the ambient row is iron.")
    iron_epsilon["props"]["rho_grain_g_cm3"] = prop(8.30, "R", "seager2007", "g/cm3", basis="rho0 of the Vinet fit")
    rows["iron_epsilon"] = iron_epsilon
    rows["iron_liquid"] = new_row("iron_liquid", "liquid iron", "metal", "Fe", [], gate=PLANETARY)
    rows["iron_liquid"]["eos"].append(eos_entry(
        "bme", {"order": 4, "rho0_kg_m3": 7019.0, "k0_gpa": 109.7, "k0p": 4.66, "k0pp_per_gpa": -0.043}, "exoplex", "H",
        note="[R0111] as carried in the ExoPlex source; the original was not reached"))
    rows["iron_liquid"]["props"]["rho_grain_g_cm3"] = prop(7.019, "H", "exoplex", "g/cm3", basis="rho0 of the BM4 fit")
    rows["fesi"] = new_row("fesi", "iron silicide", "light-element compound", "FeSi", [], gate=PLANETARY)
    rows["fesi"]["props"]["rho_grain_g_cm3"] = prop(6.522, "R", "hakim2018", "g/cm3", basis="rho0 of the [R3923] refit")

    for row_id, model, rho0, k0, k0p, source, note in data.PART_B_IRON:
        rows[row_id]["eos"].append(eos_entry(
            model, ({"order": 3} if model == "bme" else {}) | {"rho0_kg_m3": rho0, "k0_gpa": k0, "k0p": k0p}, source,
            "C" if source == "part-b-fit" else "R", note=note))

    rows["mgfe_bridgmanite"] = new_row("mgfe_bridgmanite", "iron-bearing bridgmanite (Mg0.88,Fe0.12)SiO3", "silicate", None, [],
                                       gate=PLANETARY, composition="(Mg0.88,Fe0.12)SiO3 perovskite")
    rows["mgfe_bridgmanite"]["props"]["rho_grain_g_cm3"] = prop(4.26, "R", "seager2007", "g/cm3", basis="rho0 of Table 1")

    # the SLB phases
    slb_ids = {
        "Forsterite (olivine)": "forsterite", "Fayalite": "fayalite", "Mg-wadsleyite": "mg_wadsleyite", "Fe-wadsleyite": "fe_wadsleyite",
        "Mg-ringwoodite": "mg_ringwoodite", "Fe-ringwoodite": "fe_ringwoodite", "Mg-bridgmanite": "mg_bridgmanite",
        "Fe-bridgmanite": "fe_bridgmanite", "Mg-post-perovskite": "mg_post_perovskite", "Fe-post-perovskite": "fe_post_perovskite",
        "Ca-perovskite": "ca_perovskite", "Enstatite (orthopyroxene)": "enstatite", "Ferrosilite": "ferrosilite",
        "HP-clinoenstatite": "hp_clinoenstatite", "Clinoenstatite": "clinoenstatite", "Diopside": "diopside",
        "Mg-akimotoite": "mg_akimotoite", "Periclase": "periclase", "Wuestite": "wustite", "Wuestite, low spin": "wustite",
        "Quartz": "quartz", "Coesite": "coesite", "Stishovite": "stishovite", "alpha-PbO2 SiO2 (seifertite)": "seifertite",
        "Corundum": "corundum", "Al-perovskite": "al_perovskite", "Al-post-perovskite": "al_post_perovskite",
        "Anorthite": "anorthite", "Pyrope": "pyrope", "Grossular": "grossular", "Mg-majorite": "mg_majorite", "Spinel": "spinel",
        "Hematite": "hematite", "Magnetite": "magnetite",
    }
    row_formula = {"Mg-wadsleyite": "Mg2SiO4", "Fe-wadsleyite": "Fe2SiO4", "Mg-ringwoodite": "Mg2SiO4", "Fe-ringwoodite": "Fe2SiO4",
                   "Mg-bridgmanite": "MgSiO3", "Fe-bridgmanite": "FeSiO3", "Mg-post-perovskite": "MgSiO3", "Fe-post-perovskite": "FeSiO3",
                   "Ca-perovskite": "CaSiO3", "HP-clinoenstatite": "Mg2Si2O6", "Clinoenstatite": "Mg2Si2O6", "Mg-akimotoite": "MgSiO3",
                   "Coesite": "SiO2", "Stishovite": "SiO2", "alpha-PbO2 SiO2 (seifertite)": "SiO2", "Al-perovskite": "Al2O3",
                   "Al-post-perovskite": "Al2O3", "Mg-majorite": "Mg4Si4O12"}
    for phase in parse_slb():
        row_id = slb_ids[phase["phase"]]
        burnman_unit = {"Periclase": "Mg4O4", "Wuestite": "Fe4O4", "Wuestite, low spin": "Fe4O4"}.get(phase["phase"], row_formula.get(
            phase["phase"], None))
        if row_id not in rows:
            rows[row_id] = new_row(row_id, phase["phase"].lower().replace("alpha-pbo2 sio2 ", ""), "silicate", row_formula[phase["phase"]],
                                   [], gate=PLANETARY)
            rows[row_id]["props"]["rho_grain_g_cm3"] = prop(phase["rho0"], "H", "slb2024-burnman", "g/cm3", basis="rho0 at 1 bar and 300 K")
            rows[row_id]["note"] = "Stable in the mantle: " + (phase["stable"] or "no range given in part B") + "."
        formula = burnman_unit or rows[row_id]["formula"]
        params = {"v0_cm3_mol": phase["v0"], "rho0_g_cm3": phase["rho0"], "k0_gpa": phase["k0"], "k0p": phase["k0p"],
                  "theta0_k": phase["theta"], "gamma0": phase["gamma"], "q0": phase["q"], "unit_formula": formula,
                  "unit_molar_mass_g_mol": phase["m"]}
        note = {"Wuestite, low spin": "the low-spin state", "Periclase": "K0 161.14 GPa here; some older tables give 177"}.get(phase["phase"], "")
        rows[row_id]["eos"].append(eos_entry("mgd-slb3", params, "slb2024-burnman", "H", note=note))

    # [R3992]
    for row_id, (model, order, rho0, k0, k0_sd, k0p, k0p_sd, note) in data.SEAGER_TABLE_1.items():
        if row_id == "ice_vii":
            continue
        if row_id not in rows:
            raise SystemExit(f"Seager row {row_id} has no material")
        params = {"order": order} if model == "bme" else {}
        params.update({"rho0_kg_m3": rho0, "k0_gpa": k0, "k0p": k0p})
        uncertainty = {k: v for k, v in (("k0_gpa", k0_sd), ("k0p", k0p_sd)) if v is not None}
        if row_id == "mg_bridgmanite":
            params["k0pp_per_gpa"] = data.SEAGER_K0PP_PEROVSKITE
        valid = {"p_gpa": [0.0, 200.0]}
        rows[row_id]["eos"].append(eos_entry(model, params, "seager2007", "R", valid=valid, uncertainty=uncertainty, note=note))
    # --- water: the liquid and the ice phases of part C
    ice_rows = {
        "ice_ii": ("ice II", "II"), "ice_iii": ("ice III", "III"), "ice_v": ("ice V", "V"), "ice_vi": ("ice VI", "VI"),
        "ice_vii": ("ice VII", "VII"), "ice_x": ("ice X (and superionic water above about 100 GPa and 2000 K)", "X"),
    }
    for row_id, (name, phase) in ice_rows.items():
        row = new_row(row_id, name, "ice", "H2O", [], gate=PLANETARY if row_id in ("ice_vii", "ice_x") else ANY, gas="H2O")
        rows[row_id] = row
    rows["water_liquid"] = new_row("water_liquid", "liquid water", "liquid", "H2O", [], gas="H2O")
    for phase_id in ["water_ice"] + list(ice_rows):
        points = [p for p in data.ICE_POINTS if p[0] == phase_id]
        rows[phase_id]["eos"].append(eos_entry(
            "reference-points",
            {"columns": ["p_mpa", "t_k", "rho_kg_m3", "k_t_gpa", "dk_dp", "alpha_per_k"], "points": [list(p[1:]) for p in points]},
            "journaux2020", "H", valid=data.ICE_VALID[phase_id],
            note="density, isothermal bulk modulus, its pressure derivative and the expansivity evaluated from the Gibbs spline "
                 "(rho = 1/(dG/dP), K = -(dG/dP)/(d2G/dP2)); the paper was not read"))
    rows["water_liquid"]["props"]["rho_grain_g_cm3"] = prop(data.LIQUID_WATER_POINT[2] / 1000.0, "C", "journaux2020", "g/cm3",
                                                            basis="1 bar and 300 K (Bollengier spline)")
    rows["water_liquid"]["eos"].append(eos_entry(
        "bme", {"order": 4, "v0_cm3_mol": 18.797, "k0_gpa": 2.06, "k0p": 6.29, "k0pp_per_gpa": None}, "exoplex", "H",
        note="the ExoPlex water BM4 constants; K'' is not given there"))
    rows["ice_vii"]["eos"].append(eos_entry(
        "bme", {"order": 3, "rho0_kg_m3": 1460.0, "k0_gpa": 23.7, "k0p": 4.15}, "seager2007", "R",
        uncertainty={"k0_gpa": 0.9, "k0p": 0.07}, valid={"p_gpa": [3.0, 60.0]}, note="Table 1, reference 6 (Hemley et al.); zero-pressure form"))
    rows["ice_x"]["eos"].append(eos_entry(
        "vinet", dict(zip(("rho0_kg_m3", "k0_gpa", "k0p"), data.SEAGER_ICE_X_VINET)), "seager2007", "R", valid={"p_gpa": [2.0, 7700.0]},
        note="Vinet fit to the DFT ices VIII and X, deviation under 2.5 percent in density (section 3.2)"))
    rows["ice_x"]["eos"].append(eos_entry(
        "tabulated", {"columns": ["p_gpa", "rho_kg_m3"], "points": [list(p) for p in data.SEAGER_TABLE_2]}, "seager2007", "R",
        valid={"p_gpa": [2.32, 7686.171]}, note="Table 2, DFT ices VIII and X; the table prints Mg/m3 and is multiplied by 1000 here"))
    rows["water_ice"]["note"] = ("Ice Ih at ambient pressure; the pressure phases are ice_ii to ice_x. Comet ice is 934 kg/m3 at low "
                                 "temperature (Blum and Bischoff).")

    for row_id, (rho0, c, n, note) in data.SEAGER_TABLE_3.items():
        target = row_id
        rows[target]["eos"].append(eos_entry(
            "modified-polytrope", {"rho0_kg_m3": rho0, "c_kg_m3_pa_n": c, "n": n}, "seager2007", "R", valid={"p_pa": [0.0, 1e16]},
            note="Table 3, rho = rho0 + c P^n, P in Pa, the Vinet or BME fit merged with Thomas-Fermi-Dirac; " + note))

    # --- gates, labels, state
    for row in rows.values():
        if row["id"] in ("mg_wadsleyite", "fe_wadsleyite", "mg_ringwoodite", "fe_ringwoodite", "mg_bridgmanite", "fe_bridgmanite",
                         "mg_post_perovskite", "fe_post_perovskite", "ca_perovskite", "hp_clinoenstatite", "mg_akimotoite", "coesite",
                         "stishovite", "seifertite", "al_perovskite", "al_post_perovskite", "mg_majorite", "mgfe_bridgmanite", "iron_epsilon",
                         "iron_liquid", "fesi"):
            row["gate"] = PLANETARY
        labels = []
        for key, entry in row["props"].items():
            labels.append(f"{PROPERTIES[key][1]}={entry['code']}")
        row["labels"] = ";".join(dict.fromkeys(labels))
        weak = [f"{key}: {DENSITY_LIMIT.get(e['source'], 'second-hand')}" if key == "rho_grain_g_cm3" and e["code"] in ("H", "S") else
                f"{key}: second-hand" for key, e in row["props"].items() if e["code"] in ("H", "S")]
        weak += [f"{key}: no source" for key, e in row["props"].items() if e["code"] == "N"]
        weak += [f"eos {e['model']}: second-hand" for e in row["eos"] if e["code"] == "H"]
        sources = []
        for entry in row["props"].values():
            if entry["source"] and entry["source"] not in sources:
                sources.append(entry["source"])
        for entry in row["eos"]:
            if entry["source"] not in sources:
                sources.append(entry["source"])
        if row["members"]:
            row["source"] = "a group of the rows it names; no number of its own"
            row["state"], row["limit"] = "verified", ""
        else:
            row["source"] = "; ".join(sources) if sources else "part E of the body composition research (row listed with no value yet)"
            row["state"] = "partly verified" if weak or not sources else "verified"
            row["limit"] = "; ".join(weak) if row["state"] == "partly verified" else ""
            if row["state"] == "partly verified" and not row["limit"]:
                row["limit"] = "no property value yet"

    # molar-mass check against the element table
    for row in rows.values():
        if row["formula"]:
            elements.molar_mass(row["formula"])

    result = {
        "retrieved": RETRIEVED,
        "title": "The shared materials and components table: every element, mineral, ice, organic and alloy a body can contain, with the "
                 "properties and equation-of-state constants of each (GEN.427)",
        "notes": "Properties: a number, its unit, the source key and the code of part E's label column (R read, C computed, H second-hand, "
                 "S snippet, N not reached); a value with no source is null with the code N. Equations of state: the model and its "
                 "constants as the source gives them (pressures in GPa, densities in kg/m3 unless the name says otherwise), the range "
                 "they are valid over, the source and the code. Gate groups: `any` (every planet and moon class and every small-body "
                 "kind) and `planetary` (planet and moon classes only: the high-pressure phases and the envelopes). A group row lists "
                 "member rows and carries no number. The molar mass is computed from the formula with `physics.elements`.",
        "sources": SOURCES, "legacy_names": {"class_compositions": data.CLASS_COMPOSITIONS},
        "materials": list(rows.values()),
    }
    text = json.dumps(result, indent=1, ensure_ascii=False, sort_keys=False) + "\n"
    if args.check:
        with open(OUT, encoding="utf-8") as handle:
            if handle.read() != text:
                raise SystemExit("materials.json differs from what scripts/build_materials_table.py builds; run it and commit the file")
        print("materials.json is current")
        return
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(text)
    print(f"wrote {OUT}: {len(rows)} rows")


if __name__ == "__main__":
    main()
