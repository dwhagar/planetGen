#!/usr/bin/env python3
"""
Builds `src/planetgen/physics/data/thermo_nasa7.json`, the NASA seven-coefficient thermodynamic polynomials of the gas species the
equilibrium chemistry (GEN.518, `planetgen/physics/equilibrium.py`) reads.

The coefficients are those of McBride, Gordon and Reno, NASA Technical Memorandum 4513 (1993), a United States government
publication. They are read here from the `nasa_gas.yaml` file of the Cantera project (the same coefficients in YAML), which is
downloaded for the build and is not part of the repository; only the fitted coefficients of the species below are kept.

Usage: python3 scripts/build_thermo_table.py --source /path/to/nasa_gas.yaml
"""

import argparse
import json
import os
import sys

import yaml

SPECIES = ("H2 H O O2 OH H2O HO2 H2O2 O3 C CH CH2 CH3 CH4 C2 C2H4 C2H6 CO CO2 HCO CH3OH HCN CN N NH NH2 NH3 N2 N2O NO2 "
           "S SO SO2 SO3 H2S COS CS2 He Ar F HF CL HCL Na K Ti TiO V VO Si SiO SiH4 P PH3").split()
"""list: The species kept (NASA names)."""

ELEMENTS = {"H": "H", "C": "C", "N": "N", "O": "O", "S": "S", "HE": "He", "AR": "Ar", "F": "F", "CL": "Cl", "NA": "Na", "K": "K",
            "TI": "Ti", "V": "V", "SI": "Si", "P": "P", "E": "E"}
"""dict: The upper-case element names of the source file and the symbols they stand for."""

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "planetgen", "physics", "data", "thermo_nasa7.json")


def build(source):
    with open(source, encoding="utf-8") as handle:
        document = yaml.safe_load(handle)
    by_name = {str(s["name"]): s for s in document["species"]}
    missing = [n for n in SPECIES if n not in by_name]
    if missing:
        raise SystemExit(f"the source has no species {', '.join(missing)}")
    rows = []
    for name in SPECIES:
        item = by_name[name]
        thermo = item["thermo"]
        ranges, data = thermo["temperature-ranges"], thermo["data"]
        if thermo["model"] != "NASA7" or len(data) not in (1, 2) or len(ranges) != len(data) + 1:
            raise SystemExit(f"{name}: not a one- or two-range NASA7 polynomial")
        if len(data) == 1:  # one polynomial for the whole range: the same coefficients on both sides of the midpoint
            ranges, data = [ranges[0], ranges[1], ranges[1]], [data[0], data[0]]
        composition = {}
        for key, count in item["composition"].items():
            symbol = ELEMENTS.get(str(key).upper())
            if symbol is None:
                raise SystemExit(f"{name}: the element {key} is not in the table")
            composition[symbol] = int(count)
        rows.append({"name": name, "composition": composition, "t_low": ranges[0],
                     "t_mid": ranges[1], "t_high": ranges[2], "low": data[0], "high": data[1], "note": str(thermo.get("note", "")).strip()})
    return {"title": "NASA seven-coefficient thermodynamic polynomials of the gas species the equilibrium chemistry reads (GEN.518)",
            "source": "McBride, Gordon and Reno 1993, NASA Technical Memorandum 4513, 'Coefficients for Calculating Thermodynamic "
                      "and Transport Properties of Individual Species' (a United States government publication)",
            "retrieved": "2026-10-24", "standard_state_bar": 1.0, "species": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--source", required=True, help="the nasa_gas.yaml file")
    args = parser.parse_args()
    table = build(args.source)
    with open(OUT, "w", encoding="utf-8") as handle:
        json.dump(table, handle, indent=1)
        handle.write("\n")
    print(f"wrote {os.path.normpath(OUT)} ({len(table['species'])} species)")


if __name__ == "__main__":
    sys.exit(main())
