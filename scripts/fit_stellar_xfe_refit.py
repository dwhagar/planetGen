#!/usr/bin/env python3
"""
Refits the generated stellar [X/Fe] to the [R0535] table of 1,615 FGK dwarfs (GEN.988, GEN.989, GEN.990, GEN.991).

    python scripts/fit_stellar_xfe_refit.py brewer2016_vizier_J_ApJS_225_32.tsv

Prints the constants `tuning.py` takes (they are copied by hand; nothing is written):

* the thin-disc [Mg/Fe] trend (GEN.989), with it the [Si/Fe] one and the [O/Fe] one (GEN.991; the first fit's line fell to -0.14 dex at
  [Fe/H] 0.3 where the table's stars stay at -0.03): the `("min", c, b, k)` form of `tuning.STELLAR_XFE_FITS`, robust
  least squares to the stars younger than 8 Gyr (the thin disc of `tuning.STELLAR_THICK_START_GY`) at [Fe/H] -0.6 and up. Magnesium
  alone would pull the Mg/Si ratio of every rock 0.07 dex below the table's, where the generated silicon sits 0.03 to 0.05 dex above
  it, so both come to the table's scale together;
* the residual correlation of C, O, Mg, Si and Ca (the same table): the residuals are the [X/Fe] less the median of their 0.1 dex
  [Fe/H] bin, so the trend does not enter;
* the residual scatter of C, O, Mg, Si and Ca by branch (GEN.988, GEN.990): one scale per element on the ratio the branches already
  have, chosen so that the variance in each of the bins of `src/tests/data/brewer_fischer_xfe_bands.json` is the table's robust one,
  `((p84 - p16) / 2)^2`, the variance of the generated means across the bin and the branch mixture included. Ca had none; it takes
  Mg's thick-to-thin ratio.

The robust spread is used and not the standard deviation because the table has a few stars far from the rest ([Ca/Fe] of the
0.2 to 0.4 bin has a standard deviation of 0.105 and a robust one of 0.025); a generator with no outliers is held to the body of the
table. The table's spread includes its own measurement error, which is not taken off. The TSV is a raw dataset, held in the project
files under reports/gap-filling-plan-2026-10-17-files/ and not in the repository.
"""

import json
import os
import sys

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "src"))

import build_brewer_fischer_bands as bands  # noqa: E402
from planetgen import tuning  # noqa: E402
from planetgen.physics import stellar_abundances as sa  # noqa: E402

ELEMENTS = ("C", "O", "Mg", "Si", "Ca")
BAND_FILE = os.path.join(HERE, "..", "src", "tests", "data", "brewer_fischer_xfe_bands.json")
DRAWS = 20000
MIN_STARS = 100
FIRST_SIGMA = {"thin": {"C": 0.044, "O": 0.052, "Mg": 0.066, "Si": 0.042, "Ca": 0.040},
               "high-alpha": {"C": 0.050, "O": 0.065, "Mg": 0.073, "Si": 0.073, "Ca": 0.040 * 0.073 / 0.066}}
"""dict: The scatter the branches had before this refit ([R0349, R3294]'s residuals) and calcium's start, Mg's thick-to-thin ratio on
an arbitrary 0.04: only the ratio between the branches of an element is kept, the scale is solved for."""


def fit_thin(table, element, start):
    """`(c, b, k)` of `c + b * min(F - k, 0)` to the [X/Fe] of the thin-disc stars (younger than the thick-disc start) at [Fe/H] >= -0.6."""
    fe_h, age = table["[Fe/H]"], table["Age"]
    xfe = table[f"[{element}/H]"] - fe_h
    keep = np.isfinite(age) & (age < tuning.STELLAR_THICK_START_GY) & (fe_h >= -0.6) & np.isfinite(xfe)
    x, y = fe_h[keep], xfe[keep]

    def residual(p):
        return p[0] + p[1] * np.minimum(x - p[2], 0.0) - y

    best = least_squares(residual, start, loss="soft_l1", f_scale=0.05)
    return tuple(float(v) for v in best.x), int(keep.sum())


def correlation(table):
    """The correlation of the five residuals: [X/Fe] less its median in the 0.1 dex bin of [Fe/H], over the stars with all five."""
    fe_h = table["[Fe/H]"]
    values = {e: table[f"[{e}/H]"] - fe_h for e in ELEMENTS}
    ok = np.all([np.isfinite(values[e]) for e in ELEMENTS], axis=0) & np.isfinite(fe_h)
    residual = {e: np.full(fe_h.size, np.nan) for e in ELEMENTS}
    edges = np.arange(-0.7, 0.65, 0.1)
    for low, high in zip(edges, edges[1:]):
        inside = ok & (fe_h >= low) & (fe_h < high)
        if inside.sum() >= 20:
            for e in ELEMENTS:
                residual[e][inside] = values[e][inside] - np.median(values[e][inside])
    used = np.isfinite(residual["C"])
    return np.array([[np.corrcoef(residual[a][used], residual[b][used])[0, 1] for b in ELEMENTS] for a in ELEMENTS]), int(used.sum())


def scatter_scales(thin_fits):
    """
    The scale of each element's residual scatter (thin, thick) that puts the generated variance of every populous bin on the table's.

    The generated variance of a bin is the variance of the means across the bin and the branches, plus the residual variance
    averaged over the branch mixture; the means are taken with the new thin magnesium and silicon trends.
    """
    document = json.load(open(BAND_FILE, encoding="utf-8"))
    chosen = [b for b in document["bins"] if b["stars"] >= MIN_STARS]
    rng = np.random.default_rng(988)
    fits = {branch: dict(fit) for branch, fit in tuning.STELLAR_XFE_FITS.items()}
    for element, fit in thin_fits.items():
        fits["thin"][element] = ("min",) + tuple(fit)
    result = {}
    for element in ELEMENTS:
        thin0, thick0 = FIRST_SIGMA["thin"][element], FIRST_SIGMA["high-alpha"][element]
        weights, means_var, mean_sigma2, target = [], [], [], []
        for entry in chosen:
            fe_h = rng.uniform(entry["fe_h_low"], entry["fe_h_high"], DRAWS)
            old = rng.random(DRAWS) < entry["older_than_thick_start"]
            u = rng.random(DRAWS)
            branches = [sa.branch_from_uniform(f, 10.0 if o else 4.0, v) for f, o, v in zip(fe_h, old, u)]
            means = np.array([sa._fit_value(fits[b][element], min(max(f, tuning.STELLAR_FE_H_FIT_RANGE[b][0]),
                                                                  tuning.STELLAR_FE_H_FIT_RANGE[b][1])) for f, b in zip(fe_h, branches)])
            thick = np.array([b == sa.BRANCH_HIGH_ALPHA for b in branches])
            weights.append(entry["stars"])
            means_var.append(float(means.var()))
            mean_sigma2.append(float(np.mean(np.where(thick, thick0, thin0) ** 2)))
            target.append(((entry[element]["p84"] - entry[element]["p16"]) / 2.0) ** 2)
        w, v, s2, t = (np.asarray(a) for a in (weights, means_var, mean_sigma2, target))
        scale_squared = max(float(np.sum(w * (t - v)) / np.sum(w * s2)), 0.0)
        scale = scale_squared ** 0.5
        result[element] = (round(thin0 * scale, 3), round(thick0 * scale, 3))
    return result


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__.strip().splitlines()[0], file=sys.stderr)
        return 2
    table = bands.read(argv[0])
    thin = {}
    for element, start in (("Mg", [-0.05, -0.3, 0.1]), ("Si", [0.0, -0.17, -0.07]), ("O", [-0.03, -0.4, 0.07])):
        (c, b, k), n = fit_thin(table, element, start)
        thin[element] = (round(c, 3), round(b, 3), round(k, 2))
        print(f"thin [{element}/Fe] (GEN.989), {n} stars: (\"min\", {thin[element][0]}, {thin[element][1]}, {thin[element][2]})")
    matrix, used = correlation(table)
    print(f"residual correlation of {', '.join(ELEMENTS)} (GEN.988), {used} stars:")
    for row in matrix:
        print("    (" + ", ".join(f"{value:.2f}" for value in row) + "),")
    scales = scatter_scales(thin)
    print("residual scatter (thin, high-alpha) (GEN.988, GEN.990):")
    for element, (thin, thick) in scales.items():
        print(f"    {element}: {thin}, {thick}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
