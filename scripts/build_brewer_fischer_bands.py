#!/usr/bin/env python3
"""
Reduces the [R0535] table of 1,615 FGK dwarfs (VizieR J/ApJS/225/32, SPOCS) to the bands the generated stellar
abundances are held to (TEST.456): per [Fe/H] bin, the number of stars, the share older than 8 Gyr (the thick-disc age of
`tuning.STELLAR_THICK_START_GY`), and for [C/Fe], [O/Fe], [Mg/Fe], [Si/Fe] and [Ca/Fe] the median, the standard deviation and the 16th
and 84th percentiles.

    python scripts/build_brewer_fischer_bands.py brewer2016_vizier_J_ApJS_225_32.tsv

The TSV is the VizieR `asu-tsv` answer (`-source=J/ApJS/225/32/stars`, `-out.all=1`); it is not in the repository (a raw dataset)
and is held in the project files under reports/gap-filling-plan-2026-10-17-files/. The reduced table is written to
`src/tests/data/brewer_fischer_xfe_bands.json`; the stars themselves are not kept.
"""

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "src", "tests", "data", "brewer_fischer_xfe_bands.json")
EDGES = (-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6)
ELEMENTS = ("C", "O", "Mg", "Si", "Ca")
MIN_STARS = 30
THICK_AGE_GY = 8.0


def read(path):
    """The table as `{column name: float array}` (blank cells are NaN)."""
    with open(path, encoding="utf-8") as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if not line.startswith("#") and line.strip()]
    header = rows[0]
    data = [row for row in rows[2:] if len(row) == len(header)]
    table = {}
    for i, name in enumerate(header):
        values = []
        for row in data:
            try:
                values.append(float(row[i]))
            except ValueError:
                values.append(float("nan"))
        table[name] = np.asarray(values)
    return table


def bands(table):
    fe_h, age = table["[Fe/H]"], table["Age"]
    out = []
    for low, high in zip(EDGES, EDGES[1:]):
        inside = (fe_h >= low) & (fe_h < high)
        if inside.sum() < MIN_STARS:
            continue
        entry = {"fe_h_low": low, "fe_h_high": high, "stars": int(inside.sum()),
                 "older_than_thick_start": round(float(np.nanmean(age[inside] > THICK_AGE_GY)), 3)}
        for element in ELEMENTS:
            offset = (table[f"[{element}/H]"] - fe_h)[inside]
            offset = offset[~np.isnan(offset)]
            p16, p50, p84 = np.percentile(offset, [16, 50, 84])
            entry[element] = {"n": int(offset.size), "median": round(float(p50), 3), "sd": round(float(offset.std()), 3),
                              "p16": round(float(p16), 3), "p84": round(float(p84), 3)}
        out.append(entry)
    return out


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__.strip().splitlines()[0], file=sys.stderr)
        return 2
    table = read(argv[0])
    document = {"source": "[R0535] (VizieR J/ApJS/225/32), 1,615 SPOCS dwarfs; [X/Fe] = [X/H] - [Fe/H]",
                "thick_start_gy": THICK_AGE_GY, "bins": bands(table)}
    with open(OUT, "w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=1)
        handle.write("\n")
    print(f"wrote {os.path.relpath(OUT)}: {len(document['bins'])} bins")
    return 0


if __name__ == "__main__":
    sys.exit(main())
