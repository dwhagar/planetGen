#!/usr/bin/env python3
# scripts/fit_comet_admission.py

"""
Fits src/planetgen/physics/data/comet_admission.json (GEN.507)
==============================================================

Which comets have a published coma composition, and at what heliocentric distance it was measured. Two published tables are joined:

* the compositional measurements of comets, one row per species and epoch with the comet's heliocentric distance `Rh(au)`
  ([R3782], database [R4972], CC BY 4.0), 909 rows of 96
  comets;
* the orbits of the known comets (JPL Small-Body Database query API: `q` and the orbit class of every comet).

The comets are grouped as the generator groups them (`jupiter_family` = JFc, JFC, ETc; `halley_type` = HTC; `long_period` = COM, HYP, PAR) and
the perihelion distance is binned (under 1, 1 to 2, 2 AU and over). Two things are written:

1. `admission`: the share of known comets in each group and perihelion bin that have a measurement in the database, with the counts
   (comets with `0.3 <= q <= 8` AU, so the sun-grazers that nobody measures do not set the base rate);
2. `distance`: `ln(Rh) = intercept + slope ln(q) + scatter eps` for the rows of two species kinds, parent molecules (CO, CO2, HCN,
   CH3OH, ...) and radicals and dust (C2, C3, CN, CH, NH, NH2, OH, Afrho, ...), with the number of rows and comets.

Only the fitted numbers are written; the database's rows are not kept.

    python3 scripts/fit_comet_admission.py                        query and write the file
    python3 scripts/fit_comet_admission.py --csv db.csv --sbdb sbdb.json     read saved copies instead
    python3 scripts/fit_comet_admission.py ... --check            fail if the file is not what the saved copies would write

Needs network access for the queries.
"""

import csv
import io
import json
import math
import os
import statistics
import sys
import urllib.parse
import urllib.request

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "src", "planetgen", "physics", "data", "comet_admission.json")
DATABASE_URL = "https://datashare.ed.ac.uk/bitstreams/6850db26-82a6-4bbe-adf9-487e99bf2cca/download"
SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb_query.api?" + urllib.parse.urlencode({"fields": "pdes,class,q", "sb-kind": "c"})

GROUPS = {"JFc": "jupiter_family", "JFC": "jupiter_family", "ETc": "jupiter_family", "HTC": "halley_type",
          "COM": "long_period", "HYP": "long_period", "PAR": "long_period"}
PERIHELION_BINS = ((0.3, 1.0), (1.0, 2.0), (2.0, 8.0))
RADICALS = {"C2", "C3", "CN", "CH", "NH", "NH2", "Afrho", "OH", "OI", "CS", "CS2"}


def load(path_or_url, binary=False):
    if path_or_url.startswith("http"):
        with urllib.request.urlopen(path_or_url, timeout=120) as response:
            return response.read().decode("utf-8")
    with open(path_or_url, encoding="utf-8") as handle:
        return handle.read()


def comets(sbdb_text):
    document = json.loads(sbdb_text)
    fields = document["fields"]
    out = {}
    for raw in document["data"]:
        row = dict(zip(fields, raw))
        try:
            out[row["pdes"]] = (row["class"], float(row["q"]))
        except (TypeError, ValueError):
            continue
    return out


def comet_key(row, known):
    """The designation the database row has in the orbit table, or None."""
    if row["type"] == "P" and row["number"]:
        number = int(float(row["number"]))
        for suffix in ("P", "D"):
            if f"{number}{suffix}" in known:
                return f"{number}{suffix}"
        return None
    return row["desig"] if row["desig"] in known else None


def regression(points):
    """`(intercept, slope, scatter)` of ln Rh on ln q."""
    x, y = [p[1] for p in points], [p[0] for p in points]
    mean_x, mean_y = statistics.fmean(x), statistics.fmean(y)
    slope = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y)) / sum((a - mean_x) ** 2 for a in x)
    intercept = mean_y - slope * mean_x
    residual = [b - (intercept + slope * a) for a, b in zip(x, y)]
    return intercept, slope, math.sqrt(sum(r * r for r in residual) / (len(x) - 2))


def build(database_text, sbdb_text, retrieved):
    known = comets(sbdb_text)
    measured, points = set(), {"parent": [], "radical": []}
    for row in csv.DictReader(io.StringIO(database_text)):
        key = comet_key(row, known)
        if key is None:
            continue
        measured.add(key)
        try:
            rh = float(row["Rh(au)"])
        except (KeyError, ValueError):
            continue
        if rh > 0:
            points["radical" if row["X"] in RADICALS else "parent"].append((math.log(rh), math.log(known[key][1]), key))
    admission = []
    for group in ("jupiter_family", "halley_type", "long_period"):
        for low, high in PERIHELION_BINS:
            members = [k for k, (cls, q) in known.items() if GROUPS.get(cls) == group and low <= q < high]
            admission.append({"group": group, "q_min_au": low, "q_max_au": high, "comets": len(members),
                              "measured": sum(1 for k in members if k in measured)})
    distance = {}
    for kind, pts in points.items():
        intercept, slope, scatter = regression(pts)
        distance[kind] = {"intercept": round(intercept, 4), "slope": round(slope, 4), "scatter": round(scatter, 4),
                          "rows": len(pts), "comets": len({p[2] for p in pts})}
    return {
        "title": "Which comets have a measured coma composition and at what heliocentric distance (GEN.507)",
        "source": "[R3782] database ([R4972], CC BY 4.0) joined to the JPL Small-Body Database orbits; "
                  "queries in scripts/fit_comet_admission.py",
        "retrieved": retrieved,
        "label": "DATA, COMPUTED",
        "groups": GROUPS,
        "known_comets": sum(e["comets"] for e in admission),
        "measured_comets": sum(e["measured"] for e in admission),
        "admission": admission,
        "distance": distance,
    }


def main(argv):
    csv_path = argv[argv.index("--csv") + 1] if "--csv" in argv else DATABASE_URL
    sbdb_path = argv[argv.index("--sbdb") + 1] if "--sbdb" in argv else SBDB_URL
    retrieved = os.environ.get("COMET_RETRIEVED", "2026-10-23")
    rendered = json.dumps(build(load(csv_path), load(sbdb_path), retrieved), indent=1) + "\n"
    if "--check" in argv:
        with open(OUT, encoding="utf-8") as handle:
            if handle.read() != rendered:
                raise SystemExit(f"{OUT} is not what this writes; run the script")
        print("ok")
        return
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(rendered)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main(sys.argv[1:])
