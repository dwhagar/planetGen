# tests/test_gen_density_floor.py

"""
GEN.78 Some regions have a star probability of zero, and GEN.79 Bright
stars only land between layers -121 and 121 (docs/TODO.md).

GEN.78: the density model never falls below the halo floor
(`tuning.MIN_RELATIVE_DENSITY`), so the sector draw, the bright-star
scatter and the backfill give every cell inside the galaxy's outline some
chance of a star, at the far edge of the disk and in the halo above it.

GEN.79: old disk and bulge giants reach the 1000 Lsun bright-star threshold
(their short bright tip, `tuning.GIANT_TIP_FRACTION`), so a galaxy with a
large bulge gets bright stars far above layer 121.
"""

import math

import pytest

from planetgen import tuning
from planetgen.galaxy.density import _raw_density, build_galaxy_shape, population_densities, relative_density
from planetgen.galaxy.geometry import ring_sector_count, sector_position_pc
from planetgen.galaxy.skeleton import build_layer_extents, expected_system_count_at_density_1
from planetgen.generation import bright_stars as brightStars
from planetgen.generation.star_population import bright_star_fraction

MILKY_WAY = build_galaxy_shape(
    disk_scale_length_pc=2800.0, disk_scale_height_pc=350.0, bulge_scale_radius_pc=200.0,
    bulge_amplitude=1.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
LARGE_BULGE = build_galaxy_shape(
    disk_scale_length_pc=2800.0, disk_scale_height_pc=350.0, bulge_scale_radius_pc=1500.0,
    bulge_amplitude=20.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
EDGE_PC = tuning.DEFAULT_SECTOR_EDGE_PC
E_VALUE = expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)

FAR_POINTS = [
    (14000.0, 0.0, 0.0),        # the far edge of the disk
    (0.0, -16000.0, 300.0),
    (9000.0, 9000.0, 4000.0),   # the halo above the disk
    (0.0, 0.0, 30000.0),        # far above the centre
    (0.0, 0.0, -100000.0),      # where every profile underflows
]


@pytest.mark.parametrize("point", FAR_POINTS)
def test_the_density_never_falls_below_the_halo_floor(point):
    assert relative_density(point, MILKY_WAY) >= tuning.MIN_RELATIVE_DENSITY
    densities = population_densities(point, MILKY_WAY)
    assert sum(densities.values()) == pytest.approx(relative_density(point, MILKY_WAY))
    # The halo's stars are old.
    assert densities["old"] >= 0.99 * tuning.MIN_RELATIVE_DENSITY


def test_the_floor_leaves_the_dense_disk_alone():
    # Near the Sun (density about 1) the model is untouched.
    point = (8000.0, 0.0, 0.0)
    assert relative_density(point, MILKY_WAY) == MILKY_WAY.k_norm * _raw_density(point, MILKY_WAY)


def test_the_floor_does_not_widen_the_outline():
    # A sector at the floor expects fewer than one star, so the skeleton's
    # one-star outline (which reads the model, not the floor) is unchanged.
    assert tuning.MIN_RELATIVE_DENSITY * E_VALUE < 1.0


def _halo_and_edge_cells(shape):
    """Cells inside the outline where the model alone predicts under one
    star: the top layer's rings and the plane's outermost ring."""
    extents, outer, _confirmed = build_layer_extents(shape, EDGE_PC, 1.0 / E_VALUE)
    top_layer, top_ring = extents[0]
    cells = [(ring, top_layer, slot) for ring in range(top_ring + 1) for slot in range(ring_sector_count(ring))]
    cells += [(outer, 0, slot) for slot in range(0, ring_sector_count(outer), 97)]
    return cells


def test_every_cell_at_the_far_edge_and_in_the_halo_can_get_a_bright_star():
    cells = _halo_and_edge_cells(MILKY_WAY)
    fractions = brightStars.band_fractions(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL)
    for address in cells:
        densities = population_densities(sector_position_pc(*address, EDGE_PC), MILKY_WAY)
        mean = E_VALUE * sum(densities[population] * fractions[population] for population in densities)
        assert mean > 0.0, address
    # The layer scatter's own per-ring means (no bin is zeroed).
    top_layer, top_ring = build_layer_extents(MILKY_WAY, EDGE_PC, 1.0 / E_VALUE)[0][0]
    assert brightStars.layer_expected_stars(MILKY_WAY, top_layer, top_ring, EDGE_PC, E_VALUE, fractions) > 0.0


def test_the_backfill_draws_stars_in_sparse_cells():
    # At a dense enough setting that a few hundred sparse cells should get
    # some: before GEN.78 the backfill skipped every one of them.
    cells = _halo_and_edge_cells(MILKY_WAY)
    rows = list(brightStars.backfill_cells(MILKY_WAY, cells, EDGE_PC, E_VALUE * 1000.0, 200.0, None, 5))
    assert rows
    assert {(row[0], row[1], row[2]) for row in rows} <= set(cells)


# --- GEN.79 ---------------------------------------------------------------

@pytest.mark.parametrize("population", ["old", "bulge"])
def test_old_and_bulge_stars_reach_the_bright_star_threshold(population):
    assert bright_star_fraction(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL, population) > 0.0
    # Still far rarer than among young stars.
    assert (bright_star_fraction(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL, population)
            < bright_star_fraction(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL, "young") / 5)


@pytest.mark.parametrize("layer_index", [400, -800])
def test_a_large_bulge_gets_bright_stars_far_above_layer_121(layer_index):
    extents = dict(build_layer_extents(LARGE_BULGE, EDGE_PC, 1.0 / E_VALUE)[0])
    assert extents[layer_index] > 60
    # Its inner rings, where the bulge is.
    rows = list(brightStars.scatter_layer(LARGE_BULGE, layer_index, 60, EDGE_PC, E_VALUE,
                                          tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL, 3))
    assert len(rows) > 20
    assert {row[1] for row in rows} == {layer_index}
    assert {row[6] for row in rows} <= {"bulge", "old"}


# --- GEN.117 --------------------------------------------------------------

INNER_RINGS = 242
"""int: The inner disk out to about 970 pc (Block 0-2 on the Galaxy Map)."""


def _scatter(layer_index):
    return list(brightStars.scatter_layer(MILKY_WAY, layer_index, INNER_RINGS, EDGE_PC, E_VALUE,
                                          tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL, 7))


@pytest.mark.parametrize("layer_index", [121, -121])
def test_bright_stars_rise_far_above_the_thin_young_band(layer_index):
    # About 485 pc from the plane, where young stars are all but gone,
    # old disk and bulge giants still give hundreds of bright stars.
    rows = _scatter(layer_index)
    populations = [row[6] for row in rows]
    assert len(rows) > 200
    assert populations.count("old") + populations.count("bulge") > 0.8 * len(rows)


def test_the_bright_stars_z_spread_by_population():
    plane = [row[6] for row in _scatter(0)]
    high = [row[6] for row in _scatter(60)]       # about 240 pc up
    # The plane is mostly young stars; 240 pc up young ones are rare,
    # but every older population is still there in strength.
    assert plane.count("young") > 0.5 * len(plane)
    assert high.count("young") < 0.1 * len(high)
    for population in ("intermediate", "old", "bulge"):
        assert high.count(population) > 0.1 * plane.count(population), population
    # Altogether, 240 pc up holds more than a tenth of the plane's count.
    assert len(high) > 0.1 * len(plane)
