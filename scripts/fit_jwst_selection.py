#!/usr/bin/env python3
# scripts/fit_jwst_selection.py

"""
Fits src/planetgen/physics/data/jwst_selection.json (GEN.506)
=============================================================

Which transiting planets have a published JWST spectrum. The sample is the NASA Exoplanet Archive's table of composite planetary
parameters (`pscomppars`): the transiting planets with a radius, a mass (measured or from the archive's mass-radius relation), an
equilibrium temperature, the host radius and effective temperature and the host's J and K magnitudes. The label is a published
transmission or emission spectrum from JWST in the archive's `spectra` table (direct images are not included). A logistic regression
fits the label on log10 of the transmission spectroscopy metric of [R2315] (equation 1), log10 of their emission metric
(equation 4) and an offset for each radius class (under 1.5, 1.5 to 4, 4 to 10, over 10 Earth radii: JWST programs are chosen by class,
and the fit then reproduces the number observed in each class). The host's magnitudes are modelled too (J as a normal, J-K as a line in log10 Teff) so a
generated planet without a magnitude can be given one. Only fitted numbers are written; the archive's rows are not kept.

    python3 scripts/fit_jwst_selection.py                       query the archive and write the file
    python3 scripts/fit_jwst_selection.py --csv planets.csv --spectra spectra.csv      read saved query results instead
    python3 scripts/fit_jwst_selection.py ... --check            fail if the file is not what the saved queries would write

Needs network access for the queries (the archive is public; its acknowledgment is in DATA-LICENSE.md).
"""

import bisect
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
OUT = os.path.join(ROOT, "src", "planetgen", "physics", "data", "jwst_selection.json")
PLANET_QUERY = ("select pl_name,pl_rade,pl_bmasse,pl_eqt,st_teff,st_rad,sy_jmag,sy_kmag,tran_flag from pscomppars")
SPECTRA_QUERY = "select pl_name,spec_type,facility from spectra"
BASE = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"

TSM_SCALE = ((1.5, 0.190), (2.75, 1.26), (4.0, 1.28), (float("inf"), 1.15))
ESM_PREFACTOR = 4.29e6
RADIUS_RATIO = 0.009168
"""float: Earth radii per solar radius (6371 km / 695700 km)."""
RADIUS_BINS = ((0.0, 1.5), (1.5, 4.0), (4.0, 10.0), (10.0, float("inf")))


def fetch(query):
    url = BASE + urllib.parse.urlencode({"query": query, "format": "csv"})
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read().decode("utf-8")


def planck_7_5(temperature):
    h, c, k = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    return 1.0 / math.expm1(h * c / (7.5e-6 * k * temperature))


def tsm(radius, mass, teq, host_radius, j_mag):
    scale = next(s for limit, s in TSM_SCALE if radius < limit)
    return scale * radius ** 3 * teq / (mass * host_radius ** 2) * 10.0 ** (-j_mag / 5.0)


def esm(radius, teq, host_teff, host_radius, k_mag):
    return (ESM_PREFACTOR * planck_7_5(1.10 * teq) / planck_7_5(host_teff) * (radius / host_radius * RADIUS_RATIO) ** 2
            * 10.0 ** (-k_mag / 5.0))


def number(row, key):
    try:
        value = float(row[key])
    except (KeyError, ValueError):
        return None
    return value if math.isfinite(value) and value > 0 else None


def sample(planets_text, spectra_text):
    """`[(label, log10 TSM, log10 ESM, log10 R, J, K, Teff, R, in class 2, in class 3, in class 4)]` over the usable transiting planets."""
    observed = {row["pl_name"] for row in csv.DictReader(io.StringIO(spectra_text))
                if "James Webb" in row["facility"] and row["spec_type"] in ("Transmission", "Eclipse")}
    out = []
    for row in csv.DictReader(io.StringIO(planets_text)):
        if row.get("tran_flag") != "1":
            continue
        values = [number(row, k) for k in ("pl_rade", "pl_bmasse", "pl_eqt", "st_teff", "st_rad", "sy_jmag", "sy_kmag")]
        if None in values:
            continue
        radius, mass, teq, teff, host_radius, j_mag, k_mag = values
        out.append((row["pl_name"] in observed, math.log10(tsm(radius, mass, teq, host_radius, j_mag)),
                    math.log10(esm(radius, teq, teff, host_radius, k_mag)), math.log10(radius), j_mag, k_mag, teff, radius)
                   + tuple(1.0 if low <= radius < high else 0.0 for low, high in RADIUS_BINS[1:]))
    return out


def logistic(rows, columns):
    """Newton fit of `logit P = c0 + sum(c_i x_i)` on standardised columns; returns the raw coefficients and the log-likelihood."""
    count = len(columns) + 1
    x = [[1.0] + [r[c] for c in columns] for r in rows]
    mean = [0.0] + [statistics.fmean(row[j] for row in x) for j in range(1, count)]
    sd = [1.0] + [statistics.stdev(row[j] for row in x) for j in range(1, count)]
    z = [[(row[j] - mean[j]) / sd[j] for j in range(count)] for row in x]
    y = [1.0 if r[0] else 0.0 for r in rows]
    w = [0.0] * count
    for _ in range(60):
        g = [0.0] * count
        h = [[0.0] * count for _ in range(count)]
        for zi, yi in zip(z, y):
            p = 1.0 / (1.0 + math.exp(-sum(a * b for a, b in zip(w, zi))))
            for i in range(count):
                g[i] += (yi - p) * zi[i]
                for j in range(count):
                    h[i][j] += p * (1.0 - p) * zi[i] * zi[j]
        a = [h[i][:] + [g[i]] for i in range(count)]
        for i in range(count):
            pivot = a[i][i]
            a[i] = [v / pivot for v in a[i]]
            for r in range(count):
                if r != i:
                    factor = a[r][i]
                    a[r] = [u - factor * v for u, v in zip(a[r], a[i])]
        step = [a[i][count] for i in range(count)]
        w = [u + v for u, v in zip(w, step)]
        if max(abs(v) for v in step) < 1e-10:
            break
    raw = [0.0] * count
    for j in range(1, count):
        raw[j] = w[j] / sd[j]
    raw[0] = w[0] - sum(w[j] * mean[j] / sd[j] for j in range(1, count))
    log_likelihood = 0.0
    probabilities = []
    for row, yi in zip(x, y):
        p = 1.0 / (1.0 + math.exp(-sum(a * b for a, b in zip(raw, row))))
        probabilities.append(p)
        log_likelihood += math.log(p if yi else 1.0 - p)
    return raw, log_likelihood, probabilities


def auc(scores, labels):
    negatives = sorted(s for s, y in zip(scores, labels) if not y)
    total = 0.0
    positives = [s for s, y in zip(scores, labels) if y]
    for s in positives:
        total += bisect.bisect_left(negatives, s) + 0.5 * (bisect.bisect_right(negatives, s) - bisect.bisect_left(negatives, s))
    return total / (len(positives) * len(negatives))


def build(planets_text, spectra_text, retrieved):
    rows = sample(planets_text, spectra_text)
    raw, log_likelihood, probabilities = logistic(rows, (1, 2, 8, 9, 10))
    labels = [r[0] for r in rows]
    j = [r[4] for r in rows]
    log_teff = [math.log10(r[6]) for r in rows]
    colour = [r[4] - r[5] for r in rows]
    mean_x, mean_y = statistics.fmean(log_teff), statistics.fmean(colour)
    slope = sum((a - mean_x) * (b - mean_y) for a, b in zip(log_teff, colour)) / sum((a - mean_x) ** 2 for a in log_teff)
    intercept = mean_y - slope * mean_x
    residual = statistics.stdev(b - (intercept + slope * a) for a, b in zip(log_teff, colour))
    calibration = []
    for low, high in RADIUS_BINS:
        members = [i for i, r in enumerate(rows) if low <= r[7] < high]
        calibration.append({"radius_min_earth": low, "radius_max_earth": None if math.isinf(high) else high, "planets": len(members),
                            "observed": sum(1 for i in members if labels[i]), "predicted": round(sum(probabilities[i] for i in members), 2)})
    return {
        "title": "Which transiting planets have a published JWST spectrum (GEN.506)",
        "source": "NASA Exoplanet Archive, tables pscomppars and spectra; queries in scripts/fit_jwst_selection.py",
        "retrieved": retrieved,
        "label": "DATA, COMPUTED",
        "metrics": "[R2315]: transmission metric equation 1 with the scale factors of Table 1, "
                   "emission metric equation 4 with Tday = 1.10 Teq",
        "planets": len(rows), "observed": sum(labels), "log_likelihood": round(log_likelihood, 3),
        "auc_logistic": round(auc(probabilities, labels), 4),
        "coefficients": {"intercept": round(raw[0], 4), "log10_tsm": round(raw[1], 4), "log10_esm": round(raw[2], 4),
                         "class_offsets": [0.0, round(raw[3], 4), round(raw[4], 4), round(raw[5], 4)]},
        "radius_classes_earth": [[low, None if math.isinf(high) else high] for low, high in RADIUS_BINS],
        "host_magnitudes": {"mean_j": round(statistics.fmean(j), 4), "sd_j": round(statistics.stdev(j), 4),
                            "min_j": round(min(j), 3), "max_j": round(max(j), 3),
                            "colour_intercept": round(intercept, 4), "colour_slope_per_log10_teff": round(slope, 4),
                            "colour_scatter": round(residual, 4)},
        "calibration": calibration,
    }


def main(argv):
    options = {}
    for flag in ("--csv", "--spectra"):
        if flag in argv:
            options[flag] = argv[argv.index(flag) + 1]
    retrieved = os.environ.get("JWST_RETRIEVED", "2026-10-23")
    planets = open(options["--csv"], encoding="utf-8").read() if "--csv" in options else fetch(PLANET_QUERY)
    spectra = open(options["--spectra"], encoding="utf-8").read() if "--spectra" in options else fetch(SPECTRA_QUERY)
    rendered = json.dumps(build(planets, spectra, retrieved), indent=1) + "\n"
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
