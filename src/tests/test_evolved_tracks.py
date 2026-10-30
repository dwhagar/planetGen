"""
GEN.935: the post-main-sequence tracks built from the MIST isochrone rows (`datasets/evolved_tracks.py`), on a synthetic
set of isochrones whose tracks are known exactly, so the construction (the age a mass reaches each equivalent evolutionary
point, found across the isochrone ages) is checked without the real archive.

The synthetic star of initial mass m reaches EEP e at the age `T0 m^-2.5 (1 + (e - 454) / 1000)`, so the share of its
main-sequence lifetime it has spent after the main sequence at e is exactly `(e - 454) / 1000`, whatever its mass.
"""

import io
import math

import numpy as np
import pytest

from planetgen.datasets import evolved_tracks as et
from planetgen.datasets import isochrones

T0 = 1.0e10
LOG_AGES = np.round(np.arange(6.0, 10.3001, 0.05), 2)


def _age(eep, mass):
    return T0 * mass ** -2.5 * (1.0 + (eep - 454) / 1000.0)


def _rows(eep_step=1):
    """Post-main-sequence rows in the order of `et.ROW_COLUMNS` of the synthetic isochrones."""
    rows = []
    for log_age in LOG_AGES:
        for eep in range(454, 1711, eep_step):
            mass = (10.0 ** log_age / (T0 * (1.0 + (eep - 454) / 1000.0))) ** (-1.0 / 2.5)
            if not 0.5 <= mass <= 400.0:
                continue
            rows.append([eep, log_age, mass, mass * (1.0 - 1.0e-4 * (eep - 454)), 2.0 + 0.5 * math.log10(mass) + 1.0e-3 * (eep - 454),
                         3.6 - 1.0e-4 * (eep - 454), 1.0 + 0.1 * math.log10(mass), 0.7, 0.003, 2.0 if eep < 605 else 3.0])
    return np.array(rows)


@pytest.fixture(scope="module")
def built():
    tracks, ms_age = et.build(_rows())
    return tracks, ms_age


def test_the_share_of_the_main_sequence_lifetime_after_it_ended_is_found_exactly_across_the_isochrone_ages(built):
    tracks, ms_age = built
    node = 70
    mass = et.MASS_NODES_SOL[node]
    assert ms_age[node] == pytest.approx(T0 * mass ** -2.5, rel=1e-6)
    for j, eep in enumerate(et.SAMPLE_EEPS):
        if not np.isnan(tracks[node, j, 0]):
            assert tracks[node, j, 0] == pytest.approx((eep - 454) / 1000.0, abs=2e-5), eep
            assert tracks[node, j, 1] == pytest.approx(2.0 + 0.5 * math.log10(mass) + 1.0e-3 * (eep - 454), abs=2e-5)


def test_the_present_mass_and_the_other_columns_are_read_at_the_same_point(built):
    tracks, _ = built
    node, j = 60, et.SAMPLE_EEPS.index(1000)
    mass = et.MASS_NODES_SOL[node]
    assert tracks[node, j, 2] == pytest.approx(mass * (1.0 - 1.0e-4 * 546), rel=5e-4)
    assert tracks[node, j, 3] == pytest.approx(3.6 - 1.0e-4 * 546, abs=1e-4)
    assert tracks[node, j, 4] == pytest.approx(1.0 + 0.1 * math.log10(mass), abs=1e-4)
    assert tracks[node, j, 5] == pytest.approx(0.7) and tracks[node, j, 6] == pytest.approx(0.003)
    assert tracks[node, j, 7] == 3.0 and tracks[node, et.SAMPLE_EEPS.index(502), 7] == 2.0


def test_a_mass_the_grid_does_not_reach_is_missing_where_it_is_not_held(built):
    tracks, ms_age = built
    # The lightest nodes do not reach late equivalent points inside the oldest isochrone (20 Gy) and are cut there.
    assert np.isnan(ms_age[0]) or np.isnan(tracks[0, -1, 0])
    lightest = int(np.nonzero(~np.isnan(tracks[:, et.SAMPLE_EEPS.index(808), 0]))[0][0])
    assert 0.85 < et.MASS_NODES_SOL[lightest] < 1.0
    # Every track is a run of samples with no gap, and its time never goes back.
    for node in range(len(et.MASS_NODES_SOL)):
        valid = ~np.isnan(tracks[node, :, 1])
        if valid.any():
            count = int(valid.sum())
            assert valid[:count].all() and not valid[count:].any()
            assert (np.diff(tracks[node, :count, 0]) >= 0.0).all()


def test_the_grid_blends_the_nearest_tracks_and_a_lighter_mass_takes_the_lightest_track(built):
    tracks, ms_age = built
    grid = et.Grid([-1.0, 0.0], [tracks, tracks], [ms_age / 1e9, ms_age / 1e9])
    track = grid.track(3.0, -0.4)
    assert track.mass_sol == 3.0 and track.eep[0] == 454 and track.eep[-1] == 1710
    j = int(np.searchsorted(track.eep, 1000))
    assert track.frac[j] == pytest.approx(0.546, abs=1e-4)
    assert track.log_l[j] == pytest.approx(2.0 + 0.5 * math.log10(3.0) + 0.546, abs=2e-3)
    assert track.ms_age_gy == pytest.approx(T0 * 3.0 ** -2.5 / 1e9, rel=2e-3)
    lighter = grid.track(0.8, 0.0)
    lightest = int(np.nonzero(~np.isnan(tracks[:, et.SAMPLE_EEPS.index(808), 0]))[0][0])
    assert lighter is not None and lighter.eep[-1] >= 808
    assert lighter.ms_age_gy == pytest.approx(ms_age[lightest] / 1e9)
    # A [Fe/H] past the grid takes the last one; the landmarks are the blend's samples.
    tip, helium_end, asymptotic_end, last = grid.landmark_fracs(3.0, 5.0)
    assert (tip, helium_end, asymptotic_end, last) == pytest.approx((0.151, 0.253, 0.955, 1.256), abs=1e-4)


def test_a_track_is_read_between_its_samples_inside_a_window(built):
    tracks, ms_age = built
    grid = et.Grid([0.0], [tracks], [ms_age / 1e9])
    track = grid.track(5.0, 0.0)
    star_mass, log_l, log_teff, log_r, h1, c12, phase = track.at_frac(0.1, 454, 605)
    assert log_l == pytest.approx(2.0 + 0.5 * math.log10(5.0) + 0.1, abs=2e-3)
    assert star_mass == pytest.approx(5.0 * (1.0 - 1.0e-4 * 100), rel=1e-3) and phase == 2
    # A share past the window is held to its last sample, one before it to its first.
    assert track.at_frac(9.0, 454, 605)[1] == pytest.approx(track.log_l[track.index_of(605)])
    assert track.at_frac(-1.0, 606, 707)[1] == pytest.approx(track.log_l[track.index_of(606)])
    assert track.reaches(1409) and not track.reaches(1711)


def test_an_isochrone_set_without_rows_after_the_main_sequence_gives_empty_tracks():
    tracks, ms_age = et.build(np.zeros((0, 10)))
    assert np.isnan(tracks).all() and np.isnan(ms_age).all()
    grid = et.Grid([0.0], [tracks], [ms_age])
    assert grid.track(1.0, 0.0) is None and not grid.has_tracks()


def test_the_derived_file_of_a_real_archive_layout_holds_the_tracks(tmp_path, monkeypatch):
    """`isochrones.derive` writes the tracks of every [Fe/H] beside the main sequence, and its reader serves them."""
    lines = _iso_text(0.0)
    y_init, z_init, fe_h, rows, post, pre_edge = isochrones.parse_isochrone_file(io.BytesIO(lines).readlines())
    assert fe_h == 0.0 and rows.shape[1] == 6 and post.shape[1] == 10 and post[:, 0].min() >= 454
    assert pre_edge.shape == (0, 2)
    assert set(np.unique(post[:, 9])) == {2.0, 3.0}


def _iso_text(fe_h):
    header = ["# MIST version number  = 1.2", "# ---", "#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit",
              f"# 0.27  1.4E-02     {fe_h:.2f}     0.00     0.40", "# ---", "# number of isochrones =   2"]
    body = []
    for log_age in (8.0, 9.0):
        for eep, mass, phase in ((10, 1.0, 0.0), (11, 1.1, 0.0), (454, 1.2, 2.0), (700, 1.2, 3.0)):
            columns = [0.0] * 25
            columns[0], columns[1], columns[2], columns[3], columns[24] = eep, log_age, mass, mass * 0.99, phase
            body.append(" ".join(f"{value:.6E}" for value in columns))
    return ("\n".join(header + body) + "\n").encode()
