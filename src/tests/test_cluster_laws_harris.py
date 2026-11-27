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


# --- GEN.297: the generated population against the file's distributions (statistically, not row for row) ---

from planetgen.generation import clusters as generator  # noqa: E402
from planetgen.physics import king  # noqa: E402

GENERATED = 3000
"""int: Clusters drawn for the comparison (one seed): the percentiles of 3,000 are good to a few percent, far inside the bands below."""
METAL_RICH_SPLIT = -1.0
"""float: The [Fe/H] the file's two groups are told apart at (the dip of the bimodal distribution)."""
PROJECTED_TO_3D = 1.33
"""float: The factor from a half-light radius seen in projection to the half-mass radius in 3D for a cluster-like profile (the generator draws the 3D one)."""


@pytest.fixture(scope="module")
def population():
    try:
        file = harris_clusters.clusters()
    except fetch.DataMissing:
        pytest.skip("the Harris catalog is not downloaded")
    return file, generator.generate_clusters(bytes(range(16)), GENERATED)


def _percentile(values, fraction):
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(fraction * len(ordered)))]


def _within_dex(generated, observed, dex, label):
    for fraction in (0.16, 0.5, 0.84):
        got, want = _percentile(generated, fraction), _percentile(observed, fraction)
        assert abs(math.log10(got / want)) <= dex, (label, fraction, got, want)


def _distance_kpc(cluster):
    return math.sqrt(sum(c * c for c in cluster.center_pc)) / 1000.0


def test_the_halo_clusters_lie_where_the_files_metal_poor_clusters_do(population):
    """GEN.297: the metal-poor group's distance from the center follows the cored power law of [R5363] (row GC1 of 6.335), whose 16th, 50th and 84th percentiles fall
    within 0.2 dex of the file's metal-poor ones (2.1, 7.1 and 21 kpc), and all clusters within 0.2 dex of the file's (1.9, 5.1 and 19)."""
    file, made = population
    poor = lambda f: f is not None and f <= METAL_RICH_SPLIT
    _within_dex([_distance_kpc(c) for c in made if poor(c.metallicity_fe_h)],
                [c.r_gc_kpc for c in file if poor(c.fe_h) and c.r_gc_kpc is not None], 0.2, "metal-poor distance")
    _within_dex([_distance_kpc(c) for c in made], [c.r_gc_kpc for c in file if c.r_gc_kpc is not None], 0.2, "distance")


def test_the_metal_rich_clusters_are_nearer_the_center_than_the_metal_poor_as_in_the_file(population):
    file, made = population
    rich = [_distance_kpc(c) for c in made if c.metallicity_fe_h > METAL_RICH_SPLIT]
    poor = [_distance_kpc(c) for c in made if c.metallicity_fe_h <= METAL_RICH_SPLIT]
    assert _percentile(rich, 0.5) < _percentile(poor, 0.5)
    _within_dex(rich, [c.r_gc_kpc for c in file if c.fe_h is not None and c.fe_h > METAL_RICH_SPLIT and c.r_gc_kpc is not None], 0.2, "metal-rich distance")


def test_the_metallicity_has_the_files_two_groups_and_a_metal_rich_share_inside_the_band_between_galaxies(population):
    """The two peaks sit within 0.15 dex of the file's group medians; the metal-rich share is held to 25 to 45 percent, the band of the published red shares (30 to 40 percent
    of the clusters of early-type hosts, [R5362]; the file's 30 percent and the Milky Way's 28 percent at -0.75 [R5361]), not set to the file's one galaxy."""
    file, made = population
    values = [c.fe_h for c in file if c.fe_h is not None]
    generated = [c.metallicity_fe_h for c in made]
    for want in (lambda f: f <= METAL_RICH_SPLIT, lambda f: f > METAL_RICH_SPLIT):
        theirs = sorted(f for f in values if want(f))
        ours = sorted(f for f in generated if want(f))
        assert abs(ours[len(ours) // 2] - theirs[len(theirs) // 2]) <= 0.15
    share = sum(1 for f in generated if f > METAL_RICH_SPLIT) / len(generated)
    assert 0.25 <= share <= 0.45, share
    assert 0.25 <= sum(1 for f in values if f > METAL_RICH_SPLIT) / len(values) <= 0.45


def test_the_king_concentration_and_the_collapsed_share_match_the_files(population):
    """The concentration `c = log10(rt / rc)` has its median within 0.15 and its 16th percentile within 0.25 of the file's (1.5 and 0.85); the core-collapsed share is within 6 points of the file's 18 percent."""
    file, made = population
    ours = [math.log10(c.tidal_radius_pc / king.scale_radius_pc(c.king_w0, c.half_mass_radius_pc)) for c in made]
    theirs = [c.concentration for c in file if c.concentration is not None]
    assert abs(_percentile(ours, 0.5) - _percentile(theirs, 0.5)) <= 0.15
    assert abs(_percentile(ours, 0.16) - _percentile(theirs, 0.16)) <= 0.25
    assert abs(sum(c.core_collapsed for c in made) / len(made) - sum(c.core_collapsed for c in file) / len(file)) <= 0.06


def test_the_half_mass_radius_matches_the_files_half_light_radius_seen_in_projection(population):
    file, made = population
    theirs = [c.half_light_radius_arcmin * c.r_sun_kpc * 1000.0 * math.pi / (180.0 * 60.0) * PROJECTED_TO_3D
              for c in file if c.half_light_radius_arcmin and c.r_sun_kpc]
    _within_dex([c.half_mass_radius_pc for c in made], theirs, 0.2, "half-mass radius")
