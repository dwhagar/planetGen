#!/usr/bin/env python3
"""
Builds the species thermodynamics file (GEN.451, GEN.454, GEN.499) from the NASA Glenn `thermo.inp` it is read from.

    python scripts/build_species_thermo.py PATH_TO_THERMO_INP           # write src/planetgen/physics/data/nasa_glenn_species.json
    python scripts/build_species_thermo.py PATH_TO_THERMO_INP --check   # build in memory and fail if the committed file differs

`thermo.inp` is the NASA Glenn thermodynamic database of the CEA program ([R2968], NASA/TP-2002-211556,
9-coefficient fits; the file at https://github.com/nasa/cea/blob/main/data/thermo.inp, a work of the United States
Government). The raw file (1.2 MB, 2111 species) is never in the repository: this script reads the rows of the species
`WANTED` lists, by their exact names in that file, and writes only those. A name that is not in the file stops the build.

The file holds, per species: its elements, its phase (`g` gas, `c` condensed), its molar mass, the enthalpy of formation at
298.15 K and the temperature intervals with their seven coefficients `a` and two integration constants `b`:

    H/RT = -a1/T^2 + a2 ln(T)/T + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4 + a7 T^4/5 + b1/T
    S/R  = -a1/(2 T^2) - a2/T + a3 ln T + a4 T + a5 T^2/2 + a6 T^3/3 + a7 T^4/4 + b2
"""

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "src", "planetgen", "physics", "data", "nasa_glenn_species.json")

SOURCE = (
    "[R2968], NASA/TP-2002-211556 (NASA Glenn 9-coefficient fits), as released in the thermo.inp "
    "of NASA CEA (github.com/nasa/cea, data/thermo.inp)"
)

# The names are the species' names in thermo.inp. Gases first, then the condensed phases the equilibrium may form. The
# reason for each group is the research part that needs it.
WANTED = {
    # Hydrogen, helium and the noble gases (parts 01, 05 section 2.2)
    "H": "g", "H2": "g", "He": "g", "Ne": "g", "Ar": "g", "Kr": "g", "Xe": "g",
    # Oxygen, water and its relatives
    "O": "g", "O2": "g", "O3": "g", "OH": "g", "H2O": "g", "H2O2": "g",
    # Carbon (part 05 section 3.3, 7.3)
    "C": "g", "CO": "g", "CO2": "g", "CH4": "g", "C2H2,acetylene": "g", "C2H4": "g", "C2H6": "g", "HCN": "g",
    "CH3OH": "g", "HCOOH": "g",
    # Nitrogen
    "N": "g", "N2": "g", "NH3": "g", "NO": "g", "NO2": "g", "N2O": "g",
    # Sulfur
    "S": "g", "S2": "g", "H2S": "g", "SO": "g", "SO2": "g", "SO3": "g", "COS": "g", "CS2": "g",
    # Phosphorus and the halogens
    "P2": "g", "PH3": "g", "HF": "g",
    # Alkalis, silicon, iron, magnesium, titanium, vanadium, aluminium, calcium (lava vapour, hot giants, part 05 sections 4, 7, 8)
    "Na": "g", "NaOH": "g", "Na2": "g", "NaO": "g", "K": "g", "KOH": "g",
    "AL": "g", "K2": "g", "KO": "g", "K2O": "g", "Na2O": "g",     # the rest of the lava vapour's alkali and aluminium species (GEN.751; "AL" is the file's spelling)
    "ALO": "g", "AL2O": "g", "ALO2": "g", "AL2O2": "g",     # the aluminium oxides over a melt with an Al2O3 activity (GEN.818)
    "Fe": "g", "FeO": "g", "Fe(OH)2": "g", "Mg": "g", "MgO": "g", "MgOH": "g", "Mg(OH)2": "g",
    "Si": "g", "SiO": "g", "SiO2": "g", "SiH4": "g", "SiS": "g",
    "Ti": "g", "TiO": "g", "TiO2": "g", "V": "g", "VO": "g", "Cr": "g", "Ca": "g", "CaO": "g", "CaOH": "g",
    "Li": "g", "LiOH": "g", "Mn": "g",
    # Condensed phases (the rainout of part 05 section 7.2 and 12, the melt oxides of the lava worlds)
    "H2O(cr)": "c", "H2O(L)": "c", "C(gr)": "c",
    "Fe(a)": "c", "Fe(L)": "c", "MgO(cr)": "c", "MgO(L)": "c", "SiO2(a-qz)": "c", "SiO2(L)": "c",
    "Na(cr)": "c", "Na(L)": "c", "K(cr)": "c", "K(L)": "c",
    "CaO(cr)": "c", "CaO(L)": "c", "TiO2(cr)": "c", "TiO2(L)": "c",
    "Na2O(L)": "c", "K2O(L)": "c", "AL2O3(L)": "c",     # the melt's alkali oxides (the Na and K vapour over a lava world, part 05 section 4)
}

ATOMIC_MASS_TOLERANCE = 0.01  # g/mol between the file's molar mass and the sum over its elements (physics.elements)


def parse_float(text):
    return float(text.strip().replace("D", "E").replace("d", "e"))


def read_species(path):
    """Every species of `thermo.inp` as a dict name -> row (elements, phase, molar mass, formation enthalpy, intervals)."""
    with open(path, encoding="utf-8") as handle:
        lines = handle.read().split("\n")
    species = {}
    index = 0
    while index < len(lines) - 1:
        line = lines[index]
        header = lines[index + 1] if index + 1 < len(lines) else ""
        if line and line[0] not in " !-" and re.match(r"^\s*\d\s", header) and len(header) >= 65 and not line.startswith(("thermo", "END")):
            name = line[:18].strip()
            count = int(header[0:2])
            elements = {}
            for slot in range(5):
                symbol = header[10 + 8 * slot:12 + 8 * slot].strip()
                number = header[12 + 8 * slot:18 + 8 * slot].strip()
                if symbol and number and float(number) != 0.0:
                    elements[symbol.capitalize()] = float(number)
            phase = int(header[51:52] or 0)
            molar_mass = parse_float(header[52:65])
            hf = parse_float(header[65:80]) if header[65:80].strip() else 0.0
            intervals = []
            cursor = index + 2
            for _ in range(count):
                top = lines[cursor]
                t_min, t_max = parse_float(top[0:11]), parse_float(top[11:22])
                if int(top[22:23]) != 7:
                    raise ValueError(f"{name}: expected 7 coefficients")
                first, second = lines[cursor + 1], lines[cursor + 2]
                a = [parse_float(first[16 * k:16 * k + 16]) for k in range(5)]
                a += [parse_float(second[0:16]), parse_float(second[16:32])]
                b = [parse_float(second[48:64]), parse_float(second[64:80])]
                intervals.append({"t_min_k": t_min, "t_max_k": t_max, "a": a, "b": b})
                cursor += 3
            if count == 0:
                cursor = index + 2
            species.setdefault(name, {"elements": elements, "condensed": phase != 0, "molar_mass_g_mol": molar_mass,
                                      "formation_enthalpy_j_mol": hf, "intervals": intervals})
            index = cursor
        else:
            index += 1
    return species


def build(path):
    species = read_species(path)
    rows = {}
    for name, kind in WANTED.items():
        if name not in species:
            raise SystemExit(f"{name!r} is not in {path}")
        row = species[name]
        if row["condensed"] != (kind == "c"):
            raise SystemExit(f"{name!r}: the file says condensed={row['condensed']}")
        if not row["intervals"]:
            raise SystemExit(f"{name!r} has no fitted interval")
        rows[name] = row
    return {
        "title": "NASA Glenn 9-coefficient thermodynamic fits of the gases and condensed phases of the atmosphere chemistry (ideal gas, 1 bar standard state)",
        "source": SOURCE,
        "form": ("H/RT = -a1/T^2 + a2 ln(T)/T + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4 + a7 T^4/5 + b1/T; "
                 "S/R = -a1/(2 T^2) - a2/T + a3 ln T + a4 T + a5 T^2/2 + a6 T^3/3 + a7 T^4/4 + b2; G/RT = H/RT - S/R"),
        "species": rows,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("thermo_inp", help="path to the NASA CEA thermo.inp")
    parser.add_argument("--check", action="store_true", help="fail if the committed file differs")
    args = parser.parse_args(argv)
    text = json.dumps(build(args.thermo_inp), indent=1, ensure_ascii=False) + "\n"
    if args.check:
        with open(OUT, encoding="utf-8") as handle:
            if handle.read() != text:
                print("src/planetgen/physics/data/nasa_glenn_species.json differs from the build", file=sys.stderr)
                return 1
        return 0
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(text)
    print(f"wrote {os.path.relpath(OUT)} ({len(WANTED)} species)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
