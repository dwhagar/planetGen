"""PERF.89: the array form of `population_densities` and the bright-star stages that use it."""
import math
import random

import numpy as np

from planetgen import tuning
from planetgen.galaxy.density import build_galaxy_shape, population_densities, population_densities_array
from planetgen.generation import bright_stars
from planetgen.generation.bright_stars import _densities

SHAPES = [
    build_galaxy_shape(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
                       bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4),
    build_galaxy_shape(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=4000.0,
                       bulge_amplitude=20.0, arm_count=4, pitch_angle_rad=math.radians(12), arm_amplitude=0.6),
]


def test_array_form_matches_one_point_form():
    rng = random.Random(89)
    for shape in SHAPES:
        for z in (0.0, 150.0, -900.0, 6000.0, 60000.0):
            xs = [rng.uniform(-20000, 20000) for _ in range(60)] + [0.0]
            ys = [rng.uniform(-20000, 20000) for _ in range(60)] + [0.0]
            arrays = population_densities_array(np.array(xs), np.array(ys), z, shape)
            for i, (x, y) in enumerate(zip(xs, ys)):
                scalar = population_densities((x, y, z), shape)
                for population, value in scalar.items():
                    assert math.isclose(arrays[population][i], value, rel_tol=1e-12, abs_tol=1e-300), (
                        population, x, y, z)


def test_bright_star_bins_match_one_point_form():
    shape = SHAPES[0]
    edge = tuning.DEFAULT_SECTOR_EDGE_PC
    shell, high = 40, 3
    bands = bright_stars._band_bins(shell, 0, high, shape, edge)
    assert len(bands) == 2 * high + 1
    for band, (wedges, bins) in zip(range(-high, high + 1), bands):
        assert wedges == bright_stars.wedge_count(shell, band) and len(bins) == min(wedges, bright_stars.ANGLE_BINS)
        plane, z = bright_stars.band_center_pc(shell, band, edge)
        for k, got in enumerate(bins):
            theta = (k + 0.5) * 2 * math.pi / len(bins)
            want = _densities((plane * math.cos(theta), plane * math.sin(theta), z), shape)
            for population, value in want.items():
                assert math.isclose(got[population], value, rel_tol=1e-12, abs_tol=1e-300)


def test_shell_expected_stars_matches_one_point_form():
    shape = SHAPES[0]
    edge = tuning.DEFAULT_SECTOR_EDGE_PC
    fractions = {population: 0.3 for population in bright_stars.POPULATIONS}
    for shell, high in ((0, 0), (4, 3), (30, 12)):
        got = bright_stars.shell_expected_stars(shape, shell, high, edge, 2.0, fractions)
        bands = 2 * high + 1
        count = min(bands, bright_stars.WEIGHT_BAND_SAMPLES)
        step = bands / count
        want = 0.0
        for k in range(count):
            band = -high + min(int((k + 0.5) * step), bands - 1)
            wedges = bright_stars.wedge_count(shell, band)
            bins = min(wedges, bright_stars.WEIGHT_ANGLE_BINS)
            plane, z = bright_stars.band_center_pc(shell, band, edge)
            band_expected = 0.0
            for j in range(bins):
                theta = (j + 0.5) * 2 * math.pi / bins
                densities = _densities((plane * math.cos(theta), plane * math.sin(theta), z), shape)
                band_expected += sum(densities[p] * fractions[p] for p in bright_stars.POPULATIONS)
            want += 2.0 * wedges * bright_stars.volume_ratio(shell, band, edge) / bins * band_expected * step
        assert math.isclose(got, want, rel_tol=1e-12)
