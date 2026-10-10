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
from planetgen.galaxy.geometry import ring_sector_count
from tests.test_bright_star_scatter import E_VALUE, EDGE_PC, SHAPE, THRESHOLD

MILKY_WAY = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
    arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4,
)


def _ring_means(shape, layer_index, outer_ring, fractions, e=E_VALUE, edge=EDGE_PC, points=4000):
    """The expected stars of each ring of a layer, integrating the density over the ring's cells (Monte Carlo,
    fixed seed). The old ring-by-ring walk read the density at bin centres on the ring's centreline, which
    overstates it where the density curves across the ring."""
    rng = random.Random(0)
    means = []
    for ring_index in range(outer_ring + 1):
        slots = ring_sector_count(ring_index)
        total, used = 0.0, 0
        for _ in range(points):
            point = brightStars._point_in_cell(rng, ring_index, layer_index, rng.randrange(slots), slots, edge)
            if point is None:
                continue
            densities = brightStars._densities(tuple(v / brightStars.MPC_PER_PC for v in point), shape)
            total += sum(densities[p] * fractions[p] for p in brightStars.POPULATIONS)
            used += 1
        means.append(e * slots * total / used)
    return means


@pytest.mark.parametrize("shape, layer", [(SHAPE, 0), (SHAPE, 1), (SHAPE, -2), (MILKY_WAY, 0), (MILKY_WAY, 40)])
def test_the_majorant_is_never_below_the_true_density(shape, layer):
    """Random points in every ring of the layer, for fractions of different populations."""
    rng = random.Random(layer + 7)
    outer = 60
    for fractions in ({"young": 1.0, "intermediate": 1.0, "old": 1.0, "bulge": 1.0},
                      {"young": 0.9, "intermediate": 0.1, "old": 0.0, "bulge": 0.02},
                      {"young": 0.0, "intermediate": 0.0, "old": 0.7, "bulge": 0.0},
                      {"young": 0.01, "intermediate": 0.3, "old": 0.05, "bulge": 0.9}):
        majorants = object_first.ring_majorants(shape, layer, outer, EDGE_PC, fractions)
        for ring in range(0, outer + 1, 3):
            for _ in range(60):
                slots = ring_sector_count(ring)
                slot = rng.randrange(slots)
                point = brightStars._point_in_cell(random.Random(rng.random()), ring, layer, slot, slots, EDGE_PC)
                if point is None:
                    continue
                densities = brightStars._densities(tuple(v / brightStars.MPC_PER_PC for v in point), shape)
                wanted = sum(densities[p] * fractions[p] for p in brightStars.POPULATIONS)
                assert wanted <= majorants[ring] * (1 + 1e-9), (ring, layer, fractions)


def test_check_mode_runs_a_whole_layer_without_the_majorant_being_exceeded():
    for layer in (-1, 0, 1):
        list(brightStars.scatter_layer(SHAPE, layer, 8, EDGE_PC, E_VALUE, THRESHOLD, 5, check=True))


def test_the_mean_count_matches_the_expected_count():
    fractions = brightStars.band_fractions(THRESHOLD)
    expected = sum(_ring_means(SHAPE, 0, 8, fractions))
    counts = [len(list(brightStars.scatter_layer(SHAPE, 0, 8, EDGE_PC, E_VALUE, THRESHOLD, seed)))
              for seed in range(60)]
    assert expected > 5
    assert statistics.mean(counts) == pytest.approx(expected, abs=5 * math.sqrt(expected / len(counts)) + 0.03 * expected)
    # A Poisson count: its variance is about its mean.
    assert statistics.pvariance(counts) == pytest.approx(expected, rel=0.5)


def test_the_stars_fall_in_the_rings_and_populations_the_density_gives():
    fractions = brightStars.band_fractions(THRESHOLD)
    per_ring = _ring_means(SHAPE, 0, 8, fractions)
    rows = [row for seed in range(120) for row in brightStars.scatter_layer(SHAPE, 0, 8, EDGE_PC, E_VALUE, THRESHOLD, seed)]
    seen = collections.Counter(row[0] for row in rows)
    total = sum(per_ring) * 120
    assert len(rows) == pytest.approx(total, rel=0.06)
    for ring_index, mean in enumerate(per_ring):
        expected = mean * 120
        assert seen[ring_index] == pytest.approx(expected, abs=5 * math.sqrt(expected) + 2)


def test_stars_in_skipped_sectors_are_thinned_not_moved():
    full = [row for seed in range(40) for row in brightStars.scatter_layer(SHAPE, 0, 8, EDGE_PC, E_VALUE, THRESHOLD, seed)]
    skip = {(row[0], 0, row[2]) for row in full[: len(full) // 3]}
    thinned = [row for seed in range(40) for row in
               brightStars.scatter_layer(SHAPE, 0, 8, EDGE_PC, E_VALUE, THRESHOLD, seed, skip_addresses=skip)]
    assert not any((row[0], row[1], row[2]) in skip for row in thinned)
    in_skip = sum(1 for row in full if (row[0], row[1], row[2]) in skip)
    assert len(thinned) == pytest.approx(len(full) - in_skip, rel=0.1)


def test_a_sector_takes_several_stars_up_to_its_capacity_where_the_density_is_high():
    e_value = E_VALUE * 60.0     # a few stars per sector in the middle
    rows = list(brightStars.scatter_layer(SHAPE, 0, 2, EDGE_PC, e_value, THRESHOLD, 11))
    per_sector = collections.Counter((row[0], row[1], row[2]) for row in rows)
    assert max(per_sector.values()) > 1
    fractions = brightStars.band_fractions(THRESHOLD)
    for (ring, layer, slot), held in per_sector.items():
        centre = brightStars._densities(brightStars.sector_position_pc(ring, layer, slot, EDGE_PC), SHAPE)
        expected = e_value * sum(centre[p] * fractions[p] for p in brightStars.POPULATIONS)
        assert held <= max(object_first.sector_capacity(expected), 1)


def test_capacity_grows_with_the_expected_count():
    capacities = [object_first.sector_capacity(lam) for lam in (0.0, 1e-4, 0.001, 0.01, 0.1, 1.0, 5.0, 40.0, 1e6)]
    assert capacities == sorted(capacities)
    assert capacities[0] == capacities[1] == 1 and capacities[-1] == 16
    assert object_first.sector_capacity(0.0004) == 1 and object_first.sector_capacity(0.05) >= 2


def test_an_empty_layer_costs_only_its_majorant():
    # Far above the disk: a few ring bounds, one Poisson count, no walk.
    rows = list(brightStars.scatter_layer(SHAPE, 60, 8, EDGE_PC, E_VALUE, THRESHOLD, 3, check=True))
    assert len(rows) <= 2


def test_the_same_seed_gives_the_same_stars():
    first = list(brightStars.scatter_layer(SHAPE, 0, 8, EDGE_PC, E_VALUE, THRESHOLD, 4))
    assert first == list(brightStars.scatter_layer(SHAPE, 0, 8, EDGE_PC, E_VALUE, THRESHOLD, 4)) and first
