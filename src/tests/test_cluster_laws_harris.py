"""
GEN.1525 (row 6 of the GEN.945 audit): the globular cluster mass laws of `tuning` held against the Harris catalog of 157 Milky Way clusters ([R1827], 2010 edition,
`datasets.harris_clusters`). The catalog prints each cluster's total absolute magnitude `M_V`, so a mass is `(M/L_V) * 10 ** (0.4 * (4.83 - M_V))` for a mass-to-light
ratio of 1.5 to 2.6 (the range the cluster dynamics papers give for the old populations; the registry rows AUD10 and AUD11 say the same). The tests skip where the catalog
is not downloaded (setup fetches it).
"""

import math

import pytest

from planetgen import tuning
from planetgen.datasets import fetch, harris_clusters

SUN_ABSOLUTE_V = 4.83
MASS_TO_LIGHT_RANGE = (1.5, 2.6)
BRIGHT_CUT_MV = -6.0       # the completeness of the catalog's faint end is poor; the law is fitted to the clusters brighter than this


@pytest.fixture(scope="module")
def magnitudes():
    try:
        clusters = harris_clusters.clusters()
    except fetch.DataMissing:
        pytest.skip("the Harris catalog is not downloaded")
    values = [c.m_v_total for c in clusters if c.m_v_total is not None]
    assert len(values) > 150
    return values


def _log_masses(magnitudes, mass_to_light, cut=None):
    return [math.log10(mass_to_light) + 0.4 * (SUN_ABSOLUTE_V - m) for m in magnitudes if cut is None or m < cut]


def test_the_peak_of_the_cluster_mass_function_is_the_catalogs_turnover_for_a_mass_to_light_ratio_in_range(magnitudes):
    peaks = []
    for ratio in MASS_TO_LIGHT_RANGE:
        logs = sorted(_log_masses(magnitudes, ratio, BRIGHT_CUT_MV))
        peaks.append(10.0 ** logs[len(logs) // 2])        # the median of the bright clusters, the middle of the log-normal
    assert peaks[0] <= tuning.CLUSTER_MASS_PEAK_SOL <= peaks[1] * 1.15, (peaks, tuning.CLUSTER_MASS_PEAK_SOL)


def test_the_width_of_the_cluster_mass_function_is_between_the_bright_clusters_and_all_of_them(magnitudes):
    def sigma(logs):
        mean = sum(logs) / len(logs)
        return math.sqrt(sum((x - mean) ** 2 for x in logs) / len(logs))

    bright = sigma(_log_masses(magnitudes, 2.0, BRIGHT_CUT_MV))      # 0.37 dex
    everything = sigma(_log_masses(magnitudes, 2.0))                  # 0.67 dex, the faint tail included
    assert bright <= tuning.CLUSTER_MASS_SIGMA_DEX <= everything, (bright, everything)


def test_the_mass_range_holds_all_but_the_faint_tail_and_no_cluster_is_heavier(magnitudes):
    low, high = tuning.CLUSTER_MASS_RANGE_SOL
    for ratio in MASS_TO_LIGHT_RANGE:
        masses = [10.0 ** x for x in _log_masses(magnitudes, ratio)]
        assert max(masses) <= high, ratio
        assert sum(1 for m in masses if m < low) / len(masses) < 0.12, ratio       # 8 to 10 percent of the catalog is under the floor: evolved, nearly dissolved clusters
