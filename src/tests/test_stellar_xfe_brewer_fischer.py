"""
The generated stellar [X/Fe] against the [R0535] table, not only its trend and scatter (TEST.456).

The fits of `tuning.STELLAR_XFE_FITS` were made to [R0349, R3294] and only the residual correlation of C, O, Mg and Si comes
from Brewer and Fischer (1,615 FGK dwarfs, VizieR J/ApJS/225/32). `tests/data/brewer_fischer_xfe_bands.json` holds that table reduced to
bands by [Fe/H] (built by `scripts/build_brewer_fischer_bands.py`; the stars are not kept). The generator is asked, bin by bin, for the
[X/Fe] of stars drawn across the bin with the table's share of thick-disc ages, and its median and scatter are held to the table's.

Brewer and Fischer's scatter includes their own measurement error (about 0.03 to 0.06 dex), and their [X/H] sit on their own solar
scale, so the tolerances are those of a cross-check between two calibrations: the medians within 0.10 dex, the scatter within a factor
of three either way (the tolerance of the first tests). The scatter and the thin-disc magnesium, silicon and oxygen trend were then refitted to
the table (GEN.988 to GEN.991, `scripts/fit_stellar_xfe_refit.py`): their tests hold the scatter to half again the table's robust spread
(`(p84 - p16) / 2`, which leaves out the handful of far stars the table has) and the medians to 0.05 dex.
"""

import json
import pathlib

import numpy as np
import pytest

from planetgen.physics import stellar_abundances as sa

ROOT = pathlib.Path(__file__).resolve().parent
DOCUMENT = json.loads((ROOT / "data" / "brewer_fischer_xfe_bands.json").read_text(encoding="utf-8"))
BINS = [entry for entry in DOCUMENT["bins"] if entry["stars"] >= 100]
DRAWS = 4000
THIN_AGE_GY, THICK_AGE_GY = 4.0, 10.0
MEDIAN_TOLERANCE_DEX = 0.10
SCATTER_FACTOR = 3.0
SCATTER_ROBUST = 1.5


def generated(entry, element, seed=456):
    """`DRAWS` generated [X/Fe] of `element` for stars spread across the bin, a share `older_than_thick_start` of them thick-disc old."""
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(DRAWS):
        fe_h = rng.uniform(entry["fe_h_low"], entry["fe_h_high"])
        age = THICK_AGE_GY if rng.random() < entry["older_than_thick_start"] else THIN_AGE_GY
        branch = sa.branch_from_uniform(fe_h, age, rng.random())
        values.append(sa.abundance_offsets(fe_h, branch, list(rng.standard_normal(5)))[element])
    return np.asarray(values)


def test_the_table_is_reduced_to_bins_with_enough_stars():
    assert len(BINS) >= 3 and all(entry["stars"] >= 100 for entry in BINS)
    assert DOCUMENT["thick_start_gy"] == sa.tuning.STELLAR_THICK_START_GY


@pytest.mark.parametrize("element", ["C", "O", "Mg", "Si", "Ca"])
def test_the_generated_median_is_within_a_tenth_of_a_dex_of_the_table_in_every_bin(element):
    for entry in BINS:
        shift = float(np.median(generated(entry, element))) - entry[element]["median"]
        assert abs(shift) <= MEDIAN_TOLERANCE_DEX, (
            f"[{element}/Fe] at [Fe/H] {entry['fe_h_low']} to {entry['fe_h_high']}: generated median is {shift:+.3f} dex from the table's")


@pytest.mark.parametrize("element", ["C", "O", "Mg", "Si", "Ca"])
def test_the_generated_scatter_is_within_half_again_of_the_tables_robust_scatter_in_every_bin(element):
    """GEN.988 and GEN.990: calcium gained a scatter and C, O and Si doubled theirs, each fitted to the table's `(p84 - p16) / 2`."""
    for entry in BINS:
        table = (entry[element]["p84"] - entry[element]["p16"]) / 2.0
        ratio = float(np.std(generated(entry, element))) / table
        assert 1 / SCATTER_ROBUST <= ratio <= SCATTER_ROBUST, (
            f"[{element}/Fe] at [Fe/H] {entry['fe_h_low']} to {entry['fe_h_high']}: generated scatter is {ratio:.2f} of the table's")


@pytest.mark.parametrize("element", ["O", "Mg", "Si"])
def test_the_thin_disc_oxygen_magnesium_and_silicon_medians_are_within_a_twentieth_of_a_dex_of_the_table(element):
    """GEN.989 and GEN.991: the generated [Mg/Fe] sat 0.05 to 0.07 dex above the table from [Fe/H] -0.2 up (silicon 0.03 to 0.05, oxygen 0.03 to 0.09 below)."""
    for entry in BINS:
        shift = float(np.median(generated(entry, element))) - entry[element]["median"]
        assert abs(shift) <= 0.05, f"[{element}/Fe] at [Fe/H] {entry['fe_h_low']} to {entry['fe_h_high']}: {shift:+.3f} dex"


def test_the_magnesium_to_silicon_ratio_stays_on_the_tables_scale():
    """The two thin-disc trends moved together, so the rock's Mg/Si did not drop by the 0.07 dex of magnesium alone."""
    for entry in BINS:
        table = entry["Mg"]["median"] - entry["Si"]["median"]
        shift = float(np.median(generated(entry, "Mg") - generated(entry, "Si"))) - table
        assert abs(shift) <= 0.04, f"[Mg/Si] at [Fe/H] {entry['fe_h_low']} to {entry['fe_h_high']}: {shift:+.3f} dex"


def test_calcium_sulfur_and_argon_move_together_and_calcium_is_correlated_with_the_others():
    rng = np.random.default_rng(988)
    draws = [sa.abundance_offsets(0.0, sa.BRANCH_THIN, list(rng.standard_normal(5))) for _ in range(DRAWS)]
    calcium, sulfur, argon = (np.array([d[element] for d in draws]) for element in ("Ca", "S", "Ar"))
    assert np.array_equal(calcium, sulfur) and np.array_equal(calcium, argon)
    magnesium = np.array([d["Mg"] for d in draws])
    assert abs(np.corrcoef(calcium, magnesium)[0, 1] - sa.tuning.STELLAR_XFE_CORRELATION[4][2]) < 0.05
