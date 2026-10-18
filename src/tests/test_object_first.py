# tests/test_object_first.py

"""PERF.58: the object-first star sampler matches the density it samples from, and its majorant holds."""

import collections
import math
import random
import statistics

import pytest

from planetgen import tuning
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import object_first
from planetgen.galaxy import spherical
from planetgen.galaxy.spherical import wedge_count
from tests.test_bright_star_scatter import E_VALUE, EDGE_PC, EXTENTS, SHAPE, THRESHOLD

MILKY_WAY = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
    arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4,
)


def _shell_means(shape, shell, high_band, fractions, e=E_VALUE, edge=EDGE_PC, points=600):
    """The expected stars of a shell, integrating the density over each band's cells (Monte Carlo, fixed seed):
    a band's cells hold `e` times their volume in edge cubes times the mean density of the pass's populations."""
    rng = random.Random(0)
    total = 0.0
    for band in range(-high_band, high_band + 1):
        wedges = wedge_count(shell, band)
        used, summed = 0, 0.0
        for _ in range(points):
            point = brightStars._point_in_cell(rng, shell, band, rng.randrange(wedges), edge)
            if point is None:
                continue
            densities = brightStars._densities(tuple(v / brightStars.MPC_PER_PC for v in point), shape)
            summed += sum(densities[p] * fractions[p] for p in brightStars.SCATTERED_POPULATIONS)
            used += 1
        total += e * wedges * brightStars.volume_ratio(shell, band, edge) * summed / used
    return total


def _scatter(shape, shell, high_band, seed, **kwargs):
    return list(brightStars.scatter_shell(shape, shell, high_band, EDGE_PC, E_VALUE, THRESHOLD, seed, **kwargs))


@pytest.mark.parametrize("shape, shell", [(SHAPE, 0), (SHAPE, 1), (SHAPE, 12), (MILKY_WAY, 0), (MILKY_WAY, 400),
                                          (MILKY_WAY, 3000), (MILKY_WAY, 9000)])
def test_the_majorant_is_never_below_the_true_density(shape, shell):
    """Random points in every band of the shell, for fractions of different populations."""
    rng = random.Random(shell + 7)
    high = spherical.band_range(shell)[1]
    for fractions in ({"young": 1.0, "intermediate": 1.0, "old": 1.0, "bulge": 1.0, "halo": 1.0, "cluster": 1.0,
                       "satellite": 1.0, "nuclear_cluster": 1.0, "nuclear_disc": 1.0},
                      {"young": 0.9, "intermediate": 0.1, "old": 0.0, "bulge": 0.02, "halo": 0.0, "cluster": 0.0,
                       "satellite": 0.0, "nuclear_cluster": 0.0, "nuclear_disc": 0.5},
                      {"young": 0.0, "intermediate": 0.0, "old": 0.7, "bulge": 0.0, "halo": 0.3, "cluster": 0.0,
                       "satellite": 0.0, "nuclear_cluster": 0.6, "nuclear_disc": 0.0},
                      {"young": 0.01, "intermediate": 0.3, "old": 0.05, "bulge": 0.9, "halo": 0.5, "cluster": 0.2,
                       "satellite": 0.2, "nuclear_cluster": 0.9, "nuclear_disc": 0.2},
                      {"young": 0.0, "intermediate": 0.0, "old": 0.0, "bulge": 0.0, "halo": 1.0, "cluster": 0.0,
                       "satellite": 0.0, "nuclear_cluster": 0.0, "nuclear_disc": 0.0}):
        fractions = {"satellite_young": fractions["satellite"], "satellite_intermediate": fractions["satellite"],
                     **fractions}           # a satellite's three populations take its share (GEN.526)
        majorants = object_first.band_majorants(shape, shell, high, EDGE_PC, fractions)
        assert len(majorants) == 2 * high + 1
        for band in range(-high, high + 1, max(1, high // 20)):
            for _ in range(40):
                wedge = rng.randrange(wedge_count(shell, band))
                point = brightStars._point_in_cell(random.Random(rng.random()), shell, band, wedge, EDGE_PC)
                if point is None:
                    continue
                densities = brightStars._densities(tuple(v / brightStars.MPC_PER_PC for v in point), shape)
                wanted = sum(densities[p] * fractions[p] for p in brightStars.SCATTERED_POPULATIONS)
                assert wanted <= majorants[band + high] * (1 + 1e-9), (shell, band, fractions)


@pytest.mark.parametrize("shell", [100, 400, 800, 900, 1000, 1100, 1250, 2000, 3000])
def test_the_majorant_holds_across_the_bar_gap_and_the_outer_disk(shell):
    """GEN.134: the recent star-formation factor moves young and intermediate stars by radius (a gap, a ring at its
    edge, a fall outward), and the majorant widens their weights by the factor's true range over each band's radii.
    Bands near the plane and near the axis of shells from the bulge to the outer disk, for passes that take only the
    recent populations, the old ones, or a mix."""
    rng = random.Random(shell)
    high = spherical.band_range(shell)[1]
    zero = {name: 0.0 for name in brightStars.POPULATIONS}
    passes = (dict(zero, young=1.0), dict(zero, intermediate=1.0), dict(zero, young=0.7, intermediate=0.3),
              dict(zero, old=1.0), dict(zero, young=0.01, intermediate=0.2, old=0.6, bulge=1.0))
    bands = sorted({0, 1, -1, 2, high // 8, -high // 8, high // 3, -high // 3, high // 2, high - 1, 1 - high})
    for fractions in passes:
        majorants = object_first.band_majorants(MILKY_WAY, shell, high, EDGE_PC, fractions)
        for band in bands:
            for _ in range(120):
                wedge = rng.randrange(wedge_count(shell, band))
                point = brightStars._point_in_cell(random.Random(rng.random()), shell, band, wedge, EDGE_PC)
                if point is None:
                    continue
                densities = brightStars._densities(tuple(v / brightStars.MPC_PER_PC for v in point), MILKY_WAY)
                wanted = sum(densities[p] * fractions[p] for p in brightStars.POPULATIONS)
                assert wanted <= majorants[band + high] * (1 + 1e-9), (shell, band, fractions)


def test_check_mode_runs_whole_shells_without_the_majorant_being_exceeded():
    for shell, high in EXTENTS:
        list(brightStars.scatter_shell(SHAPE, shell, high, EDGE_PC, E_VALUE, THRESHOLD, 5, check=True))


def test_the_mean_count_matches_the_expected_count():
    fractions = brightStars.band_fractions(THRESHOLD)
    expected = _shell_means(SHAPE, 4, 3, fractions)
    counts = [len(_scatter(SHAPE, 4, 3, seed)) for seed in range(60)]
    assert expected > 5
    assert statistics.mean(counts) == pytest.approx(expected, abs=5 * math.sqrt(expected / len(counts)) + 0.03 * expected)
    # A Poisson count: its variance is about its mean.
    assert statistics.pvariance(counts) == pytest.approx(expected, rel=0.5)


def test_the_stars_fall_in_the_shells_and_populations_the_density_gives():
    fractions = brightStars.band_fractions(THRESHOLD)
    per_shell = {shell: _shell_means(SHAPE, shell, high, fractions) for shell, high in EXTENTS}
    rows = [row for seed in range(120) for row in brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, E_VALUE, THRESHOLD, seed)]
    seen = collections.Counter(row[0] for row in rows)
    total = sum(per_shell.values()) * 120
    assert len(rows) == pytest.approx(total, rel=0.06)
    for shell, mean in per_shell.items():
        expected = mean * 120
        assert seen[shell] == pytest.approx(expected, abs=5 * math.sqrt(expected) + 2 + 0.03 * expected)


def test_stars_in_skipped_sectors_are_thinned_not_moved():
    full = [row for seed in range(40) for row in _scatter(SHAPE, 4, 3, seed)]
    # A third of the addresses that hold stars (not the addresses of a third of the rows: the first rows can fill
    # nearly every sector, which left 18 stars to compare and failed on a different draw stream).
    skip = set(sorted({(row[0], row[1], row[2]) for row in full})[::3])
    thinned = [row for seed in range(40) for row in _scatter(SHAPE, 4, 3, seed, skip_addresses=skip)]
    assert not any((row[0], row[1], row[2]) in skip for row in thinned)
    in_skip = sum(1 for row in full if (row[0], row[1], row[2]) in skip)
    assert len(thinned) == pytest.approx(len(full) - in_skip, rel=0.1)


def test_a_sector_takes_several_stars_up_to_its_capacity_where_the_density_is_high():
    e_value = E_VALUE * 60.0     # a few stars per sector in the middle
    rows = list(brightStars.scatter_shell(SHAPE, 3, 3, EDGE_PC, e_value, THRESHOLD, 11))
    per_sector = collections.Counter((row[0], row[1], row[2]) for row in rows)
    assert max(per_sector.values()) > 1
    fractions = brightStars.band_fractions(THRESHOLD)
    for (shell, band, wedge), held in per_sector.items():
        centre = brightStars._densities(brightStars.sector_position_pc(shell, band, wedge, EDGE_PC), SHAPE)
        expected = (e_value * brightStars.volume_ratio(shell, band, EDGE_PC)
                    * sum(centre[p] * fractions[p] for p in brightStars.POPULATIONS))
        assert held <= max(object_first.sector_capacity(expected), 1)


def test_capacity_grows_with_the_expected_count():
    capacities = [object_first.sector_capacity(lam) for lam in (0.0, 1e-4, 0.001, 0.01, 0.1, 1.0, 5.0, 40.0, 1e6)]
    assert capacities == sorted(capacities)
    assert capacities[0] == capacities[1] == 1 and capacities[-1] == 16
    assert object_first.sector_capacity(0.0004) == 1 and object_first.sector_capacity(0.05) >= 2


def test_an_empty_shell_costs_only_its_majorant():
    # Far outside the disk: a few band bounds, one Poisson count, no walk.
    rows = list(brightStars.scatter_shell(SHAPE, 300, 3, EDGE_PC, E_VALUE, THRESHOLD, 3, check=True))
    assert len(rows) <= 2


def test_the_same_seed_gives_the_same_stars():
    first = _scatter(SHAPE, 4, 3, 4)
    assert first == _scatter(SHAPE, 4, 3, 4) and first
