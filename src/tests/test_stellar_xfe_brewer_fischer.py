"""
The generated stellar [X/Fe] against the Brewer and Fischer 2016 table, not only its trend and scatter (TEST.456).

The fits of `tuning.STELLAR_XFE_FITS` were made to Bensby 2014 and Nissen 2014 and only the residual correlation of C, O, Mg and Si comes
from Brewer and Fischer (1,615 FGK dwarfs, VizieR J/ApJS/225/32). `tests/data/brewer_fischer_xfe_bands.json` holds that table reduced to
bands by [Fe/H] (built by `scripts/build_brewer_fischer_bands.py`; the stars are not kept). The generator is asked, bin by bin, for the
[X/Fe] of stars drawn across the bin with the table's share of thick-disc ages, and its median and scatter are held to the table's.

Brewer and Fischer's scatter includes their own measurement error (about 0.03 to 0.06 dex), and their [X/H] sit on their own solar
scale, so the tolerances are those of a cross-check between two calibrations: the medians within 0.10 dex, the scatter within a factor
of three either way. The strict expected failures name what the table shows and the generator lacks.
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


def generated(entry, element, seed=456):
    """`DRAWS` generated [X/Fe] of `element` for stars spread across the bin, a share `older_than_thick_start` of them thick-disc old."""
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(DRAWS):
        fe_h = rng.uniform(entry["fe_h_low"], entry["fe_h_high"])
        age = THICK_AGE_GY if rng.random() < entry["older_than_thick_start"] else THIN_AGE_GY
        branch = sa.branch_from_uniform(fe_h, age, rng.random())
        values.append(sa.abundance_offsets(fe_h, branch, list(rng.standard_normal(4)))[element])
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


@pytest.mark.parametrize("element", ["C", "O", "Mg", "Si"])
def test_the_generated_scatter_is_within_a_factor_of_three_of_the_tables_in_every_bin(element):
    for entry in BINS:
        ratio = float(np.std(generated(entry, element))) / entry[element]["sd"]
        assert 1 / SCATTER_FACTOR <= ratio <= SCATTER_FACTOR, (
            f"[{element}/Fe] at [Fe/H] {entry['fe_h_low']} to {entry['fe_h_high']}: generated scatter is {ratio:.2f} of the table's")


@pytest.mark.xfail(strict=True, reason="TEST.456 finding: calcium, sulfur and argon carry no scatter, the table's [Ca/Fe] has 0.04 to 0.10 dex")
def test_calcium_has_a_scatter_like_the_tables():
    for entry in BINS:
        assert float(np.std(generated(entry, "Ca"))) / entry["Ca"]["sd"] >= 0.5


@pytest.mark.xfail(strict=True, reason="TEST.456 finding: the generated [Mg/Fe] sits 0.05 to 0.07 dex above the table's from [Fe/H] -0.2 up")
def test_magnesium_median_is_within_a_twentieth_of_a_dex_of_the_table():
    for entry in BINS:
        shift = float(np.median(generated(entry, "Mg"))) - entry["Mg"]["median"]
        assert abs(shift) <= 0.05
