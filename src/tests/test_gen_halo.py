# tests/test_gen_halo.py

"""
GEN.234 A measured stellar halo replaces the flat density floor, and GEN.79 Bright
stars only land between layers -121 and 121 (docs/TODO.md); on the shell grid
"far above the plane" is a band near a pole and a height of more than 485 pc.

GEN.234: the density model has a halo term (a broken power law on a spheroid, 1.5e9 solar masses, 3.75e9
systems for the Milky Way) in place of GEN.78's flat 0.1% floor. It is half the floor's density at the Sun, a
fifth to 340,000th of it at 20 to 262 kpc, and it changes no cell's chance of a star inside the outline.

GEN.79: old disk and bulge giants reach the 1000 Lsun bright-star threshold
(the top of their branch in the MIST tracks), so a galaxy with a
large bulge gets bright stars far above the plane (more than 485 pc).
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import density, spherical, structure as galaxyStructure
from planetgen.galaxy.density import (
    _halo, build_galaxy_shape, integrated_masses, halo_density, halo_mass, halo_terms,
    population_densities, population_densities_array, relative_density,
)
from planetgen.galaxy.skeleton import build_shell_extents, expected_system_count_at_density_1
from planetgen.galaxy.spherical import sector_position_pc
from planetgen.generation import bright_stars as brightStars
from planetgen.generation.star_population import bright_star_fraction
from planetgen.physics import mathcheck

def tuning_local_density_per_pc3():
    from planetgen.physics import constants as physicalConstants
    from planetgen.physics.units import ly_to_pc
    return physicalConstants.LOCAL_STELLAR_DENSITY_LY3 / ly_to_pc(1.0) ** 3


MILKY_WAY = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
    bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
LARGE_BULGE = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=4000.0,
    bulge_amplitude=20.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
RECALLED_OVER_DERIVED = 0.4 / density.mean_system_mass_solar(MILKY_WAY)
"""float: The recalled mean system mass the design note's numbers were worked out at, over the derived one (GEN.545)."""
EDGE_PC = tuning.DEFAULT_SECTOR_EDGE_PC
E_VALUE = expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)
LOW_FLOOR = 1000.0
"""float: The old 1000 Lsun floor; these GEN.79/117 tests need old-giant populations that the 2500 minimum (GEN.184) no longer allows."""

FAR_POINTS = [
    (14000.0, 0.0, 0.0),        # the far edge of the disk
    (0.0, -16000.0, 300.0),
    (9000.0, 9000.0, 4000.0),   # the halo above the disk
    (0.0, 0.0, 30000.0),        # far above the centre
    (0.0, 0.0, -100000.0),      # where every disk profile underflows
]


@pytest.mark.parametrize("point", FAR_POINTS)
def test_the_populations_still_sum_to_the_density_and_the_halo_is_the_far_stars(point):
    densities = population_densities(point, MILKY_WAY)
    assert set(densities) == set(tuning.STELLAR_POPULATION_AGE_RANGES_GY)
    assert sum(densities.values()) == pytest.approx(relative_density(point, MILKY_WAY))
    assert densities["halo"] == pytest.approx(halo_density(point, MILKY_WAY))
    assert densities["cluster"] == 0.0
    x, y, z = point
    array = population_densities_array(np.array([x]), np.array([y]), z, MILKY_WAY)
    assert {name: float(value[0]) for name, value in array.items()} == pytest.approx(densities)


def test_far_above_the_disk_only_the_halo_is_left():
    point = (0.0, 0.0, -100000.0)
    densities = population_densities(point, MILKY_WAY)
    assert densities["halo"] == pytest.approx(relative_density(point, MILKY_WAY))
    assert densities["halo"] < 1e-6    # nowhere near the old floor's 1e-3


# The design note's table (docs/design/galaxy/intergalactic-and-halo-contents.md section 4), relative density in the plane.
@pytest.mark.parametrize("radius_kpc, relative", [(8.2, 1.8e-3), (50.0, 5.0e-6), (100.0, 2.2e-7), (262.0, 2.9e-9)])
def test_the_halo_matches_the_design_notes_table(radius_kpc, relative):
    # The table was worked out at the recalled 0.4 solar masses a system; the derived mean (GEN.545) scales the density.
    assert halo_density((radius_kpc * 1000.0, 0.0, 0.0), MILKY_WAY) == pytest.approx(
        relative * RECALLED_OVER_DERIVED, rel=0.15)


def test_the_halo_is_far_below_the_old_floor_in_the_distance_and_near_it_at_the_sun():
    old_floor = 1.0e-3
    assert halo_density((8200.0, 0.0, 0.0), MILKY_WAY) == pytest.approx(old_floor * 1.7 * RECALLED_OVER_DERIVED, rel=0.1)
    assert halo_density((50000.0, 0.0, 0.0), MILKY_WAY) < old_floor / 100
    assert halo_density((262000.0, 0.0, 0.0), MILKY_WAY) < old_floor / 1e5


def test_the_halo_integrates_to_its_share_of_the_visible_systems():
    visible = density.model_terms(MILKY_WAY)["visible_mass"]
    assert visible == pytest.approx(sum(integrated_masses(MILKY_WAY, steps=60).values()), rel=2e-3)
    assert halo_mass(MILKY_WAY) == pytest.approx(galaxyStructure.stellar_halo_system_fraction(MILKY_WAY.structure) * visible, rel=1e-12)
    # And the closed form agrees with adding the density up: shells of the spheroid, 4 pi q m^2 dm, m log-spaced.
    terms = halo_terms(MILKY_WAY)
    m = np.logspace(0.0, 8.0, 400001)
    profile = _halo(m * terms["halo_axis_ratio"] ** 0, 0.0, terms)
    total = float(np.trapezoid(profile * 4.0 * math.pi * terms["halo_axis_ratio"] * m ** 2 * m, np.log(m)))
    assert total == pytest.approx(halo_mass(MILKY_WAY), rel=2e-3)


def test_the_halo_is_its_share_of_the_stars_over_the_rest():
    # 1.4e9 solar masses [R1026] of the 5e10 of stars: against the other 4.86e10, the same ratio of systems.
    structure = galaxyStructure.REFERENCE_STRUCTURE
    assert structure.halo_stellar_ratio * 5.0e10 == pytest.approx(1.4e9)
    assert galaxyStructure.stellar_halo_system_fraction(structure) == pytest.approx(1.4e9 / (5.0e10 - 1.4e9))
    # and the whole model holds the systems the stellar mass makes
    held = density.total_mass(MILKY_WAY) * MILKY_WAY.k_norm * tuning_local_density_per_pc3() * density.mean_system_mass_solar(MILKY_WAY)
    assert held == pytest.approx(5.0e10, rel=1e-9)


def test_the_halo_is_continuous_at_the_break_and_flat_in_the_core():
    terms = halo_terms(MILKY_WAY)
    below = _halo(terms["halo_break_pc"] * (1 - 1e-9), 0.0, terms)
    above = _halo(terms["halo_break_pc"] * (1 + 1e-9), 0.0, terms)
    assert below == pytest.approx(above, rel=1e-6)
    assert _halo(10.0, 0.0, terms) == _halo(terms["halo_core_pc"], 0.0, terms)


def test_the_halo_is_flattened_toward_the_poles():
    # q = 0.7: a point q times as high is as far out in the spheroid's own terms.
    terms = halo_terms(MILKY_WAY)
    assert _halo(0.0, 20000.0 * terms["halo_axis_ratio"], terms) == pytest.approx(_halo(20000.0, 0.0, terms))


def test_the_halo_scales_with_the_galaxy():
    small = build_galaxy_shape(
        disk_scale_length_pc=1300.0, disk_scale_height_pc=150.0, bulge_scale_radius_pc=790.0,
        bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    assert halo_terms(small)["halo_break_pc"] == pytest.approx(halo_terms(MILKY_WAY)["halo_break_pc"] / 2)
    assert halo_mass(small) == pytest.approx(galaxyStructure.stellar_halo_system_fraction(small.structure) * density.model_terms(small)["visible_mass"],
                                             rel=1e-12)


def test_no_cell_beyond_the_halos_dense_core_expects_a_whole_system():
    # At most 0.02 systems in a cell from 8.7 kpc out, at any angle (the note: about 0.01 there).
    for radius in (8700.0, 12000.0, 25000.0, 100000.0):
        for colatitude in np.linspace(0.0, math.pi, 13):
            point = (radius * math.sin(colatitude), 0.0, radius * math.cos(colatitude))
            assert halo_density(point, MILKY_WAY) * E_VALUE < 0.02


def test_the_halo_does_not_widen_the_outline():
    # The outline reads the disk and bulge only (`skeleton._bound_raw_density_at`), and the halo is under 0.3% of the
    # one-star threshold where it ends. The galaxy holds the systems its stellar mass makes (GEN.404) over an
    # exponential disk (GEN.544) and the derived mean system mass (GEN.545), so the reference galaxy's edge sits at shell
    # 3,421 (3,454 before the primary mass function of GEN.929 gave 5% fewer systems for the same stars; 3,576 at the
    # recalled 0.4 a system, 3,762 for the density tuned to the real count at the Sun, 3,127 with the isothermal disk).
    extents, outer, _confirmed = build_shell_extents(mathcheck._reference_galaxy_shape(), EDGE_PC, 1.0 / E_VALUE)
    assert outer == 3421
    assert len(extents) == 3422
    assert sum(high for _shell, high in extents) == 848394


def _halo_and_edge_cells(shape):
    """Cells inside the outline where the model alone predicts under one
    star: the highest band of every hundredth shell and the outermost shell's
    equatorial band."""
    extents, outer, _confirmed = build_shell_extents(shape, EDGE_PC, 1.0 / E_VALUE)
    cells = [(shell, high, wedge) for shell, high in extents[::100]
             for wedge in range(0, spherical.wedge_count(shell, high), 7)]
    cells += [(outer, 0, wedge) for wedge in range(0, spherical.wedge_count(outer, 0), 5)]
    return cells


def test_every_cell_at_the_far_edge_and_in_the_halo_can_get_a_bright_star():
    cells = _halo_and_edge_cells(MILKY_WAY)
    fractions = brightStars.band_fractions(LOW_FLOOR)
    for address in cells:
        densities = population_densities(sector_position_pc(*address, EDGE_PC), MILKY_WAY)
        mean = E_VALUE * sum(densities[population] * fractions[population] for population in densities)
        assert mean > 0.0, address
    # The shell scatter's own per-shell means (no bin is zeroed).
    extents = build_shell_extents(MILKY_WAY, EDGE_PC, 1.0 / E_VALUE)[0]
    shell, high = extents[-1]
    assert brightStars.shell_expected_stars(MILKY_WAY, shell, high, EDGE_PC, E_VALUE, fractions) > 0.0


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
    assert bright_star_fraction(LOW_FLOOR, population) > 0.0
    # Still far rarer than among young stars.
    assert (bright_star_fraction(LOW_FLOOR, population)
            < bright_star_fraction(LOW_FLOOR, "young") / 5)


def _scatter(shell, high_band=None, shape=MILKY_WAY):
    high = spherical.band_range(shell)[1] if high_band is None else high_band
    return list(brightStars.scatter_shell(shape, shell, high, EDGE_PC, E_VALUE, LOW_FLOOR, 7))


def _height_pc(row):
    return row[5] / brightStars.MPC_PER_PC


def test_a_large_bulge_gets_bright_stars_far_above_the_plane():
    extents = dict(build_shell_extents(LARGE_BULGE, EDGE_PC, 1.0 / E_VALUE)[0])
    assert extents[400] == spherical.band_range(400)[1]       # the outline reaches the axis there
    rows = _scatter(400, shape=LARGE_BULGE)
    assert len(rows) > 20
    assert max(abs(_height_pc(row)) for row in rows) > 1000.0
    high = [row for row in rows if abs(_height_pc(row)) > 485.0]
    assert len(high) > 20
    populations = [row[6] for row in high]
    assert populations.count("bulge") + populations.count("old") > 0.8 * len(high)


# --- GEN.117 --------------------------------------------------------------

INNER_SHELL = 240
"""int: A shell about 960 pc out (Block 0-2 on the Galaxy Map's inner disk)."""


def test_bright_stars_rise_far_above_the_thin_young_band():
    # More than 485 pc from the plane, where young stars are all but gone,
    # old disk and bulge giants still give hundreds of bright stars.
    rows = [row for row in _scatter(INNER_SHELL) if abs(_height_pc(row)) > 485.0]
    populations = [row[6] for row in rows]
    assert len(rows) > 200
    assert populations.count("old") + populations.count("bulge") > 0.8 * len(rows)


def test_the_bright_stars_z_spread_by_population():
    # Three shells pooled: the bright intermediate-age stars are about 1% of a shell's rows, so one shell holds 0 to 2 of them
    # 240 pc up and the comparison with the plane's 8 or 9 is noise (GEN.931 moved which rows the seed draws there).
    rows = [row for shell in (INNER_SHELL, INNER_SHELL + 2, INNER_SHELL + 5) for row in _scatter(shell)]
    plane = [row[6] for row in rows if abs(_height_pc(row)) < 100.0]
    high = [row[6] for row in rows if 200.0 < abs(_height_pc(row)) < 300.0]       # about 240 pc up
    # Young stars are a clear share of the plane's, though the bar gap holds them down here (GEN.134; the rest is
    # mostly the bar bulge, which fills the inner kiloparsec) and a galaxy this heavy formed its stars early
    # (GEN.469); 240 pc up young ones are rare, but every older population is still there in strength.
    # (A share of 3% before GEN.545; the sparser galaxy draws other rows from the same seed, 2.5% of 320.)
    assert plane.count("young") > 0.02 * len(plane)
    assert high.count("young") < 0.1 * len(high)
    for population in ("intermediate", "old", "bulge"):
        assert high.count(population) > 0.1 * plane.count(population) * (100.0 / 200.0), population
