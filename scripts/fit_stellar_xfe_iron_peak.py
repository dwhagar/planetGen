#!/usr/bin/env python3
"""
Fits the generated stellar [X/Fe] of Ti, V, Cr and Mn to the [R0535] table of 1,615 FGK dwarfs (GEN.1360).

    python scripts/fit_stellar_xfe_iron_peak.py brewer2016_vizier_J_ApJS_225_32.tsv

Prints the constants `tuning.py` takes (they are copied by hand; nothing is written):

* the mean [X/Fe] of each element on each branch as `("linear", a, b)`, `a + b [Fe/H]`: a weighted least-squares line to the dwarfs, each
  weighted by the probability that the generator would put a star of its [Fe/H] and age on that branch (`stellar_abundances.branch_probability`:
  zero for a star under 8 Gyr, the logistic of [Fe/H] above), over the [Fe/H] range of the branch's fits (`tuning.STELLAR_FE_H_FIT_RANGE`). The
  branches are the repository's own, not a split in [Mg/Fe], so the fit is the mean the generator should give. The high-alpha line rests on about 40
  stars' weight;
* the residual scatter of each element on each branch: the weighted robust spread `(p84 - p16) / 2` of the residuals about that line, the
  table's measurement error included (the table's spread is what `src/tests/test_stellar_xfe_brewer_fischer.py` holds the generator to);
* the residual correlation of C, O, Mg, Si, Ca, Ti, V, Cr and Mn (the residual being [X/Fe] less its median in the 0.1 dex bin of [Fe/H],
  over the stars with all nine), whose first five rows are those of `scripts/fit_stellar_xfe_refit.py`.

The TSV is a raw dataset, held in the project files under reports/gap-filling-plan-2026-10-17-files/ and not in the repository.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "src"))

import build_brewer_fischer_bands as bands  # noqa: E402
from planetgen import tuning  # noqa: E402
from planetgen.physics import stellar_abundances as sa  # noqa: E402

FITTED = ("Ti", "V", "Cr", "Mn")
CORRELATED = ("C", "O", "Mg", "Si", "Ca") + FITTED
BRANCH_NAMES = {"thin": sa.BRANCH_THIN, "high-alpha": sa.BRANCH_HIGH_ALPHA}


def branch_fit(table, element, branch):
    """`(a, b, robust sigma, effective stars)` of the weighted line of [X/Fe] against [Fe/H] for `branch`."""
    fe_h, age = table["[Fe/H]"], table["Age"]
    xfe = table[f"[{element}/H]"] - fe_h
    high = np.array([sa.branch_probability(f, a) if np.isfinite(f) and np.isfinite(a) else np.nan for f, a in zip(fe_h, age)])
    weight = high if branch == sa.BRANCH_HIGH_ALPHA else 1.0 - high
    low, top = tuning.STELLAR_FE_H_FIT_RANGE[branch]
    keep = np.isfinite(weight) & np.isfinite(xfe) & (weight > 0.01) & (fe_h >= low) & (fe_h <= top)
    design = np.vstack([np.ones(keep.sum()), fe_h[keep]]).T
    w = weight[keep]
    a, b = np.linalg.lstsq(design * np.sqrt(w)[:, None], xfe[keep] * np.sqrt(w), rcond=None)[0]
    residual = xfe[keep] - design @ np.array([a, b])
    order = np.argsort(residual)
    cumulative = np.cumsum(w[order]) / w.sum()
    p16, p84 = residual[order][np.searchsorted(cumulative, 0.16)], residual[order][np.searchsorted(cumulative, 0.84)]
    return float(a), float(b), float((p84 - p16) / 2.0), float(w.sum())


def correlation(table):
    """The correlation of the nine residuals."""
    fe_h = table["[Fe/H]"]
    values = {e: table[f"[{e}/H]"] - fe_h for e in CORRELATED}
    ok = np.all([np.isfinite(values[e]) for e in CORRELATED], axis=0) & np.isfinite(fe_h)
    residual = {e: np.full(fe_h.size, np.nan) for e in CORRELATED}
    edges = np.arange(-0.7, 0.65, 0.1)
    for low, high in zip(edges, edges[1:]):
        inside = ok & (fe_h >= low) & (fe_h < high)
        if inside.sum() >= 20:
            for e in CORRELATED:
                residual[e][inside] = values[e][inside] - np.median(values[e][inside])
    used = np.isfinite(residual["C"])
    return np.array([[np.corrcoef(residual[a][used], residual[b][used])[0, 1] for b in CORRELATED] for a in CORRELATED]), int(used.sum())


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__.strip().splitlines()[0], file=sys.stderr)
        return 2
    table = bands.read(argv[0])
    for branch in sa.BRANCHES:
        for element in FITTED:
            a, b, sigma, stars = branch_fit(table, element, branch)
            print(f'{branch} [{element}/Fe] (GEN.1360), {stars:.0f} stars of weight: ("linear", {a:.3f}, {b:.3f}), scatter {sigma:.3f}')
    matrix, used = correlation(table)
    print(f"residual correlation of {', '.join(CORRELATED)} (GEN.1360), {used} stars:")
    for row in matrix:
        print("    (" + ", ".join(f"{value:.2f}" for value in row) + "),")
    return 0


if __name__ == "__main__":
    sys.exit(main())
