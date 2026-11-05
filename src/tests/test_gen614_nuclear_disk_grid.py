# tests/test_gen614_nuclear_disk_grid.py

"""
GEN.614: the nuclear disk is switched off when its length (`NUCLEAR_DISK_LENGTH_RATIO` times the bulge radius) is under the
grid's core ball (`NUCLEAR_RESOLUTION_PC`, 2 pc), because the cluster is softened to the grid and the disk was not (the
density-limits report, section 2.3). Galaxies with a bulge radius from 0.7 pc to 100 pc, and the Milky Way's, are built.
"""

import math

import pytest

from planetgen import tuning
from planetgen.galaxy import density

CUT_PC = tuning.NUCLEAR_RESOLUTION_PC / tuning.NUCLEAR_DISK_LENGTH_RATIO      # the bulge radius where the disk's length is 2 pc


def _shape(bulge_radius_pc):
    shape = density.build_galaxy_shape(
        disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=bulge_radius_pc, bulge_amplitude=3.11,
        arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4, nuclear_amplitude=1.0)
    structure = shape.structure._replace(nuclear_disk_mass_solar=1.5e9, nuclear_cluster_mass_solar=2.5e7)
    return shape._replace(structure=structure)


@pytest.mark.parametrize("radius", [0.7, 1.0, 2.0, 5.0, 10.0, 20.0, 34.0, 36.0, 50.0, 100.0, 1580.0])
def test_the_disk_is_on_only_when_it_is_at_least_as_long_as_the_grids_core_ball(radius):
    terms = density.model_terms(_shape(radius))
    on = terms["nuclear_disk_length_pc"] >= tuning.NUCLEAR_RESOLUTION_PC
    assert on == (radius >= CUT_PC)
    assert (terms["nuclear_disk_mass"] > 0) == on and (terms["nuclear_disk_density"] > 0) == on
    assert terms["nuclear_cluster_mass"] > 0          # the cluster stays, softened to the grid


@pytest.mark.parametrize("radius", [0.7, 5.0, 34.0])
def test_a_switched_off_disk_adds_no_density_and_no_mass_anywhere(radius):
    shape = _shape(radius)
    terms = density.model_terms(shape)
    assert density._nuclear_disk(0.0, 0.0, terms) == 0.0 and density._nuclear_disk(1.0, 0.5, terms) == 0.0
    assert float(density._nuclear_disk_array([0.0, 3.0], 0.0, terms).sum()) == 0.0
    by_name = {c.name: c for c in density.density_components(shape)}
    assert by_name["nuclear_disk"].mass == 0.0 and by_name["nuclear_disk"].mass_solar == 0.0


def test_the_disk_of_a_wide_enough_bulge_is_unchanged_and_the_milky_way_keeps_its_disk():
    terms = density.model_terms(_shape(1580.0))
    assert terms["nuclear_disk_length_pc"] == pytest.approx(90.0)
    assert terms["nuclear_disk_density"] == pytest.approx(
        terms["nuclear_disk_mass"] / (8.0 * math.pi * 90.0 ** 2 * terms["nuclear_disk_height_pc"]))


def test_no_bulge_radius_from_0_7_pc_to_100_pc_divides_by_zero_or_leaves_a_disk_denser_than_the_cut_allows():
    peak_at_cut = None
    for radius in [0.7 * 1.2 ** i for i in range(40) if 0.7 * 1.2 ** i <= 100.0] + [CUT_PC, CUT_PC * 1.0001]:
        terms = density.model_terms(_shape(radius))            # no ZeroDivisionError
        assert math.isfinite(terms["nuclear_disk_density"])
        if terms["nuclear_disk_mass"]:
            peak = terms["nuclear_disk_density"] / terms["nuclear_disk_mass"]       # per unit of the disk's mass: 1 / (8 pi L^2 h)
            peak_at_cut = max(peak_at_cut or 0.0, peak)
    assert peak_at_cut is None or peak_at_cut <= 1.0 / (8.0 * math.pi * tuning.NUCLEAR_RESOLUTION_PC ** 2
                                                        * tuning.NUCLEAR_DISK_HEIGHT_RATIO * CUT_PC) * 1.001
