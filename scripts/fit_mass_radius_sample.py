#!/usr/bin/env python3
# scripts/fit_mass_radius_sample.py

"""
Fits src/planetgen/physics/data/mass_radius_sample.json (GEN.505)
=================================================================

The error model of the sample of precisely measured planets (mass better than 25 percent, radius better than 8 percent;
[R2774]), fitted to the NASA Exoplanet Archive's table of composite planetary parameters (`pscomppars`): the transiting
planets with a radius, a mass from a measurement (`pl_bmassprov == "Mass"`) and both errors. Only the fitted numbers are written:
the archive's rows are not kept.

For each planet the radial-velocity semi-amplitude K (circular orbit, edge-on) and the transit depth (R/R*)^2 follow from its
stored period, masses and radii, and the 1-sigma errors that were published become sigma_K = K sigma_M/M and
sigma_depth = 2 depth sigma_R/R. The model regresses ln sigma_K on ln K and ln sigma_depth on ln depth, with the scatter of each and
their correlation, in two sub-samples: planets smaller than 4 Earth radii (`small`) and those around hosts cooler than 4000 K (`m-dwarf`).

    python3 scripts/fit_mass_radius_sample.py                     query the archive and write the file
    python3 scripts/fit_mass_radius_sample.py --csv path.csv      read a saved query result instead
    python3 scripts/fit_mass_radius_sample.py --check             fail if the file is not what the saved query would write

Needs network access for the query (the archive is public; its acknowledgment is in DATA-LICENSE.md).
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
OUT = os.path.join(ROOT, "src", "planetgen", "physics", "data", "mass_radius_sample.json")
COLUMNS = ("pl_name,hostname,pl_orbper,pl_rade,pl_radeerr1,pl_radeerr2,pl_bmasse,pl_bmasseerr1,pl_bmasseerr2,pl_bmassprov,"
           "st_teff,st_mass,st_rad,tran_flag")
QUERY = f"select {COLUMNS} from pscomppars"
URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode({"query": QUERY, "format": "csv"})

K_PER_EARTH_MASS = 28.4329 / 317.83
"""float: K in m/s of 1 Earth mass at 1 yr around 1 solar mass (28.4329 m/s per Jupiter mass, 317.83 Earth masses per Jupiter mass)."""
R_EARTH_KM, R_SUN_KM = 6371.0, 695700.0
MASS_CUT, RADIUS_CUT = 0.25, 0.08
SUBSAMPLES = {
    "small": {"description": "transiting planets smaller than 4 Earth radii", "radius_max_earth": 4.0, "teff_max_k": None},
    "m-dwarf": {"description": "transiting planets smaller than 4 Earth radii around hosts cooler than 4000 K",
                "radius_max_earth": 4.0, "teff_max_k": 4000.0},
}


def number(row, key):
    try:
        value = float(row[key])
    except (KeyError, ValueError):
        return None
    return value if math.isfinite(value) else None


def planets(text):
    """Yields `(ln K, ln depth, ln sigma_K, ln sigma_depth, mass error share, radius error share, teff)` per usable planet."""
    for row in csv.DictReader(io.StringIO(text)):
        if row.get("tran_flag") != "1" or row.get("pl_bmassprov") != "Mass":
            continue
        values = [number(row, k) for k in ("pl_bmasse", "pl_bmasseerr1", "pl_bmasseerr2", "pl_rade", "pl_radeerr1", "pl_radeerr2",
                                           "pl_orbper", "st_mass", "st_rad")]
        if None in values:
            continue
        mass, e1, e2, radius, r1, r2, period, host_mass, host_radius = values
        sigma_m, sigma_r = (e1 - e2) / 2.0, (r1 - r2) / 2.0
        if mass <= 0 or sigma_m <= 0 or sigma_r <= 0 or radius <= 0 or period <= 0 or host_mass <= 0 or host_radius <= 0:
            continue
        k = K_PER_EARTH_MASS * mass * (period / 365.25) ** (-1.0 / 3.0) * host_mass ** (-2.0 / 3.0)
        depth = (radius * R_EARTH_KM / (host_radius * R_SUN_KM)) ** 2 * 1e6
        yield (math.log(k), math.log(depth), math.log(k * sigma_m / mass), math.log(2.0 * depth * sigma_r / radius),
               sigma_m / mass, sigma_r / radius, radius, number(row, "st_teff"))


def regression(x, y):
    """`(intercept, slope, residual standard deviation, residuals)` of the least-squares line."""
    mean_x, mean_y = statistics.fmean(x), statistics.fmean(y)
    slope = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y)) / sum((a - mean_x) ** 2 for a in x)
    intercept = mean_y - slope * mean_x
    residuals = [b - (intercept + slope * a) for a, b in zip(x, y)]
    return intercept, slope, math.sqrt(sum(r * r for r in residuals) / (len(x) - 2)), residuals


def correlation(x, y):
    mx, my = statistics.fmean(x), statistics.fmean(y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))


def fit(rows, spec):
    chosen = [r for r in rows if r[6] < spec["radius_max_earth"] and (spec["teff_max_k"] is None
                                                                       or (r[7] is not None and r[7] < spec["teff_max_k"]))]
    ln_k, ln_d = [r[0] for r in chosen], [r[1] for r in chosen]
    a_k, b_k, s_k, res_k = regression(ln_k, [r[2] for r in chosen])
    a_d, b_d, s_d, res_d = regression(ln_d, [r[3] for r in chosen])
    passing = sum(1 for r in chosen if r[4] < MASS_CUT and r[5] < RADIUS_CUT)
    return {
        "description": spec["description"], "n": len(chosen), "passing": passing,
        "ln_sigma_k": {"intercept": round(a_k, 4), "slope": round(b_k, 4), "scatter": round(s_k, 4)},
        "ln_sigma_depth": {"intercept": round(a_d, 4), "slope": round(b_d, 4), "scatter": round(s_d, 4)},
        "residual_correlation": round(correlation(res_k, res_d), 4),
        "population": {"mean_ln_k": round(statistics.fmean(ln_k), 4), "sd_ln_k": round(statistics.stdev(ln_k), 4),
                       "mean_ln_depth": round(statistics.fmean(ln_d), 4), "sd_ln_depth": round(statistics.stdev(ln_d), 4),
                       "correlation": round(correlation(ln_k, ln_d), 4)},
    }


def build(text, retrieved):
    rows = list(planets(text))
    return {
        "title": "The error model of the sample of precisely measured planets (GEN.505)",
        "source": "NASA Exoplanet Archive, table pscomppars (composite planetary parameters), transiting planets with a measured mass; "
                  "query in scripts/fit_mass_radius_sample.py",
        "retrieved": retrieved,
        "label": "DATA, COMPUTED",
        "cut": {"mass_error_share": MASS_CUT, "radius_error_share": RADIUS_CUT,
                "source": "[R2774] (the planets with mass better than 25 percent and radius better than 8 percent)"},
        "k_per_earth_mass_ms": K_PER_EARTH_MASS,
        "subsamples": {name: fit(rows, spec) for name, spec in SUBSAMPLES.items()},
    }


def main(argv):
    csv_path = None
    if "--csv" in argv:
        csv_path = argv[argv.index("--csv") + 1]
    retrieved = os.environ.get("MASS_RADIUS_RETRIEVED", "2026-10-23")
    if csv_path:
        with open(csv_path, encoding="utf-8") as handle:
            text = handle.read()
    else:
        with urllib.request.urlopen(URL, timeout=120) as response:
            text = response.read().decode("utf-8")
    document = build(text, retrieved)
    rendered = json.dumps(document, indent=1) + "\n"
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
