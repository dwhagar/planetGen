# tests/test_cluster_count_law.py

"""
GEN.403: how many globular clusters a galaxy holds follows its halo mass (the system's mass is `eta` times the halo's, over
the mean cluster mass), with the stated scatter, from the seed; a halo-less galaxy holds none; the Milky Way's 157 is only a
plausibility target.
"""

import math
import statistics

import pytest

from planetgen import tuning
from planetgen.generation import clusters

MILKY_WAY_HALO = 1.5e12


def _seeds(count):
    return [bytes([i % 256, i // 256]) + bytes(14) for i in range(count)]


def _counts(halo, seeds=400, family_factor=1.0):
    return [clusters.cluster_count(seed, halo, family_factor) for seed in _seeds(seeds)]


def test_the_mean_cluster_mass_is_the_truncated_lognormals():
    sample = []
    stream = clusters._stream(b"mean")
    for _ in range(40000):
        sample.append(clusters._mass(stream))
    assert clusters.mean_cluster_mass_solar() == pytest.approx(statistics.fmean(sample), rel=0.02)


def test_a_milky_way_mass_halo_gives_a_count_inside_the_observed_range():
    counts = _counts(MILKY_WAY_HALO)
    assert 100 <= statistics.median(counts) <= 300         # 157 in the Harris file, 140 to 190 by the three etas
    assert min(counts) > 15 and max(counts) < 1000


def test_the_count_follows_the_halo_mass_and_the_stated_scatter():
    low_eta, high_eta = tuning.CLUSTER_SYSTEM_FRACTION_RANGE
    mean = clusters.mean_cluster_mass_solar()
    for halo in (1.0e11, 1.0e12, 1.0e13):
        counts = _counts(halo)
        floor, ceiling = low_eta * halo / mean, high_eta * halo / mean
        assert floor / 3.0 < statistics.median(counts) < ceiling * 3.0
        spread = statistics.pstdev([math.log10(c) for c in counts if c])
        assert 0.1 < spread < 0.6                           # eta's factor of 2 (0.3 dex) and the 0.2 dex scatter
    assert statistics.median(_counts(1.0e13)) > 5.0 * statistics.median(_counts(1.0e12)) > 25.0 * statistics.median(_counts(1.0e11)) / 2.0


def test_a_small_halo_holds_few_and_the_mean_count_is_exact_under_one_cluster():
    counts = _counts(3.0e9, seeds=2000)
    mean = clusters.mean_cluster_mass_solar()
    low_eta, high_eta = tuning.CLUSTER_SYSTEM_FRACTION_RANGE
    assert max(counts) <= 4 and 0 in counts
    assert low_eta * 3.0e9 / mean / 2.0 < statistics.fmean(counts) < high_eta * 3.0e9 / mean * 2.0


def test_no_halo_no_clusters_and_the_family_factor_multiplies_the_system():
    assert clusters.cluster_count(bytes(16), None) == 0 and clusters.cluster_count(bytes(16), 0.0) == 0
    one, three = _counts(MILKY_WAY_HALO, seeds=100), _counts(MILKY_WAY_HALO, seeds=100, family_factor=3.0)
    assert statistics.fmean(three) / statistics.fmean(one) == pytest.approx(3.0, rel=0.05)


def test_the_same_seed_gives_the_same_count_and_the_count_does_not_move_the_clusters():
    seed = bytes(range(16))
    assert clusters.cluster_count(seed, MILKY_WAY_HALO) == clusters.cluster_count(seed, MILKY_WAY_HALO)
    assert clusters.generate_clusters(seed, 10) == clusters.generate_clusters(seed, 30)[:10]


def test_the_old_fixed_count_is_gone():
    assert not hasattr(tuning, "CLUSTER_COUNT_MILKY_WAY")
