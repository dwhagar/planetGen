# tests/test_sector_separation_fit.py

"""GEN.1534: the Hill share of the galaxy's random placement (`tuning.SECTOR_HILL_SEPARATION_SHARE`, SEP1) is held to the nearest-neighbour distribution of the Gaia sample within 25 pc."""

import importlib.util
import pathlib

import pytest

from planetgen import tuning
from planetgen.datasets import fetch, gaia_nearby

SPEC = importlib.util.spec_from_file_location("fit_sector_separation", pathlib.Path(__file__).resolve().parents[2] / "scripts" / "fit_sector_separation.py")
fit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fit)


def test_a_placement_with_no_exclusion_is_a_random_set_and_a_big_share_pushes_the_neighbours_out():
    free, kept_free = fit.simulate(0.0, 0.07, box_pc=30.0)
    wide, kept_wide = fit.simulate(0.5, 0.07, box_pc=30.0)
    assert kept_free == 1.0 and float(wide[len(wide) // 2]) > float(free[len(free) // 2]) > 1.0
    assert float(wide[0]) >= 0.5 * 2.0 * fit.HILL_MEDIAN_PC * 0.7        # nothing nearer than half the two smallest radii summed


def test_the_ks_distance_of_a_distribution_from_itself_is_zero_and_the_band_narrows_with_the_sample():
    sample = fit.simulate(0.0, 0.07, box_pc=30.0)[0]
    assert fit.ks_distance(sample, sample) == 0.0
    assert fit.critical_ks(1000) < fit.critical_ks(100)


def test_the_tuned_share_is_inside_the_band_of_the_data_and_the_old_half_is_not():
    if fetch.locate(gaia_nearby.DATASET).missing:
        pytest.skip("the Gaia sample is not downloaded (python3 -m planetgen.cli.datasets fetch)")
    rows, count, _density = fit.fit(shares=(tuning.SECTOR_HILL_SEPARATION_SHARE, 0.5), box_pc=50.0)
    (tuned, tuned_ks, _), (old, old_ks, _) = rows
    band = fit.critical_ks(count) + 0.02            # the 50 pc box is noisier than the 70 pc one of the script
    assert tuned == tuning.SECTOR_HILL_SEPARATION_SHARE and tuned_ks <= band
    assert old == 0.5 and old_ks > 0.25
