# tests/test_unresolved_crossings.py
"""NAV.74: the mean number of unresolved keep-out spheres a straight course crosses (`galaxy.unresolved_crossings`)."""

import math

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import keepout, resolution, shape_draw, spherical, unresolved_crossings as uc
from planetgen.physics import unit_registry as units
from planetgen.physics.units import ly_to_pc

SEED = bytes(range(16))
SHAPE = shape_draw.default_shape(SEED)
EDGE = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
CEILING = 5_000.0


def skeleton(per_density_one=6.0, ceiling=CEILING, cap=tuning.FILL_LUMINOSITY_CAP_SOL):
    return store.GalaxySkeletonInfo(
        shape=SHAPE, edge_pc=EDGE, outer_shell_index=200, expected_system_count_at_density_1=per_density_one, galaxy_seed=SEED,
        reach_shell_index=500, cell_star_ceiling=ceiling, fill_luminosity_cap_sol=cap)


@pytest.fixture(scope="module")
def field():
    """One field for the module, so the cells it works out are shared (the first cells of the core are the slow ones)."""
    return uc.UnresolvedField(skeleton())


def test_a_course_far_from_the_capped_region_reads_no_cell_and_crosses_nothing():
    quiet = uc.UnresolvedField(skeleton())
    far = quiet.reach_pc + 500.0
    found = quiet.crossings((far, -300.0, 0.0), (far, 300.0, 0.0))
    assert found == uc.NONE_CROSSED
    assert quiet._cells == {}


def test_a_galaxy_with_no_capped_cell_has_no_reach():
    sparse = uc.UnresolvedField(skeleton(per_density_one=1.0e-6))
    assert sparse.reach_pc == 0.0
    assert sparse.crossings((-300.0, 0.0, 0.0), (300.0, 0.0, 0.0)) == uc.NONE_CROSSED


def test_a_zero_length_segment_crosses_nothing(field):
    assert field.crossings((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)) == uc.NONE_CROSSED


def test_inside_one_cell_the_crossings_are_the_formula(field):
    point = (60.0, 10.0, 0.0)
    address = spherical.sector_address_at(point, EDGE)
    cell = field.cell(address)
    assert cell is not None, "the test point must lie in a capped cell"
    systems, volume, mean_kg = cell
    start, end = (point[0] - 0.4, point[1], point[2]), (point[0] + 0.4, point[1], point[2])
    assert spherical.sector_address_at(start, EDGE) == address == spherical.sector_address_at(end, EDGE)
    found = field.crossings(start, end, step_pc=1.0)
    hill = field.hill_radius_pc(mean_kg, math.dist(point, (0.0, 0.0, 0.0)))
    assert found.crossings == pytest.approx(systems / volume * math.pi * hill ** 2 * 0.8, rel=1e-9)
    assert found.cells == 1 and found.unresolved_systems == pytest.approx(systems) and found.length_pc == pytest.approx(0.8)


def test_the_hill_radius_is_the_keep_out_rule_of_a_star(field):
    systems, volume, mean_kg = field.cell(spherical.sector_address_at((60.0, 10.0, 0.0), EDGE))
    radius = 60.8
    from planetgen.galaxy import enclosed_mass
    central = enclosed_mass.enclosed_mass_kg(SHAPE, SEED, radius)
    expected_km = keepout.galactic_hill_radius_km(mean_kg, radius, central)
    assert field.hill_radius_pc(mean_kg, radius) == pytest.approx(units.converter("km", "pc")(expected_km), rel=1e-9)


def test_a_finer_step_gives_the_same_answer_and_the_other_direction_too(field):
    start, end = (-90.0, 7.0, 3.0), (90.0, -7.0, -3.0)
    coarse = field.crossings(start, end)
    fine = field.crossings(start, end, step_pc=EDGE / 32)
    assert coarse.crossings == pytest.approx(fine.crossings, rel=0.03)
    assert field.crossings(end, start).crossings == pytest.approx(coarse.crossings, rel=0.02)
    assert coarse.cells > 10 and coarse.crossings > 0.0


def test_the_crossings_add_up_along_a_segment(field):
    whole = field.crossings((-90.0, 7.0, 0.0), (90.0, 7.0, 0.0), step_pc=1.0)
    first = field.crossings((-90.0, 7.0, 0.0), (0.0, 7.0, 0.0), step_pc=1.0)
    second = field.crossings((0.0, 7.0, 0.0), (90.0, 7.0, 0.0), step_pc=1.0)
    assert whole.crossings == pytest.approx(first.crossings + second.crossings, rel=0.02)
    assert whole.length_pc == pytest.approx(first.length_pc + second.length_pc)


def test_a_higher_ceiling_leaves_less_unresolved():
    low = uc.UnresolvedField(skeleton(ceiling=2_000.0))
    high = uc.UnresolvedField(skeleton(ceiling=50_000.0))
    segment = ((-60.0, 7.0, 0.0), (60.0, 7.0, 0.0))
    assert high.crossings(*segment).crossings < low.crossings(*segment).crossings


def test_a_cell_filled_whole_adds_nothing():
    whole = uc.UnresolvedField(skeleton(ceiling=1.0e12), reach_pc=200.0)
    assert whole.crossings((-60.0, 7.0, 0.0), (60.0, 7.0, 0.0)) == uc.NONE_CROSSED


def test_the_unresolved_systems_of_a_cell_are_the_cell_resolutions(field):
    address = spherical.sector_address_at((60.0, 10.0, 0.0), EDGE)
    systems, volume, mean_kg = field.cell(address)
    position = spherical.sector_position_pc(*address, EDGE)
    from planetgen.galaxy import density
    from types import SimpleNamespace
    mix = {name: max(v, 0.0) for name, v in density.population_densities(position, SHAPE).items()}
    light = SimpleNamespace(expected_systems=6.0 * max(density.relative_density(position, SHAPE), 0.0), mix=mix, volume_pc3=volume,
                            edge_pc=EDGE, shape=SHAPE)
    assert systems == pytest.approx(resolution.cell_resolution(light, CEILING, tuning.FILL_LUMINOSITY_CAP_SOL).unresolved_systems)
