# tests/test_plan_figures.py

"""
GEN.626: the plan's crest and planned-share figures come from the built shape, not from a table and a fixed system mass.
The crest is the thin disc's peak expected systems per sector on an arm crest (the number the outline walks on), it
follows the Census dial, and the planned share equals an independent integral of the model's density.
"""

import math

import numpy as np
import pytest

from planetgen.galaxy import density, plan_figures
from planetgen.generation import plan_outline
from planetgen.physics import constants as physical_constants
from planetgen.physics.units import ly_to_pc

SEED = bytes.fromhex("2139DD5841E680B009A3D17A6319B2AA")
SMALL = ["--mass", "1e8", "--disk-scale-length", "2000pc", "--flatness", "20", "--bulge-share", "0.001"]


def _outline(argv):
    return plan_outline.build(plan_outline._parser().parse_args(argv), SEED, 4.0)


def _thin_disc_crest_from_the_density(shape, e_value):
    """The thin disc's expected systems per sector on an arm crest, read from `component_values` a little off the
    centre (the arms are undefined on the axis) and carried back along the exponential."""
    radius = 100.0
    arm_angle = (shape.spiral_reference_angle_rad
                 + math.log(radius / shape.spiral_reference_radius_pc) / math.tan(shape.pitch_angle_rad))
    theta = arm_angle                                   # cos(arm_count * 0) = 1: an arm crest
    point = (radius * math.cos(theta), radius * math.sin(theta), 0.0)
    values = density.component_values(point, shape)
    assert values.arm_cosine == pytest.approx(1.0)
    thin = shape.k_norm * values.thin_profile * shape.arm_level * (1.0 + shape.arm_amplitude * values.arm_cosine)
    return e_value * thin * math.exp(radius / shape.disk_scale_length_pc)


@pytest.mark.parametrize("extra", [[], ["--richness", "1.0"], ["--richness", "3.0"]])
def test_the_crest_is_the_thin_discs_peak_expected_systems_per_sector(extra):
    built = _outline(extra)
    assert built.shape.nuclear_amplitude > 0
    assert plan_outline.report(built).crest == pytest.approx(
        _thin_disc_crest_from_the_density(built.shape, built.e_value), rel=1e-9)


def test_the_crest_follows_the_census_dial():
    low = plan_outline.report(_outline(["--richness", "1.0"])).crest
    high = plan_outline.report(_outline(["--richness", "3.0"])).crest
    assert high / low == pytest.approx(3.0, rel=1e-9)


def _independent_share(shape, e_value, radius_pc):
    """The planned share on a spherical grid (radius, cosine of the polar angle, azimuth), a different grid from
    `plan_figures`' cylindrical one."""
    radial, polar, around = 300, 500, 36
    edges = np.geomspace(0.25, radius_pc, radial + 1)
    edges[0] = 0.0
    middles = 0.5 * (edges[:-1] + edges[1:])
    volume = (edges[1:] ** 3 - edges[:-1] ** 3) / 3.0
    steps = np.sinh(np.linspace(0.0, 1.0, polar + 1) * math.asinh(40.0)) / 40.0
    cosines = 0.5 * (steps[:-1] + steps[1:])
    widths = 2.0 * np.diff(steps)
    azimuths = (np.arange(around) + 0.5) * (2.0 * math.pi / around)
    r, u = middles[:, None, None], cosines[None, :, None]
    sine = np.sqrt(1.0 - u * u)
    rho = density.total_density_array(r * sine * np.cos(azimuths)[None, None, :],
                                      r * sine * np.sin(azimuths)[None, None, :], r * u, shape)
    weight = volume[:, None, None] * widths[None, :, None] * (2.0 * math.pi / around)
    local = physical_constants.LOCAL_STELLAR_DENSITY_LY3 / ly_to_pc(1.0) ** 3
    total = shape.structure.stellar_mass_solar / density.mean_system_mass_solar(shape)
    return (rho * weight)[rho * e_value >= 1.0].sum() * local / total


@pytest.mark.parametrize("extra", [[], SMALL, ["--richness", "1.0"]])
def test_the_planned_share_equals_an_independent_integral(extra):
    built = _outline(extra)
    made = plan_outline.report(built)
    independent = _independent_share(built.shape, built.e_value, made.radius_pc)
    assert 0.0 < made.share <= 1.0
    assert made.share == pytest.approx(independent, rel=0.02)


@pytest.mark.parametrize("extra", [[], SMALL])
def test_one_integral_gives_the_share_and_the_reaching_sectors_the_two_calls_gave(extra):
    """PERF.169: the report integrates the density once and hands the figures to `reaching_sectors`; the answers are the
    ones a second integral gave, bit for bit."""
    built = _outline(extra)
    made = plan_outline.report(built)
    again = plan_figures.planned(built.shape, built.e_value, made.radius_pc)
    assert made.share == again.share
    assert made.reaching == round(plan_figures.reaching_sectors(built.shape, built.e_value, made.radius_pc, built.edge_pc))
    assert made.reaching == round(plan_figures.reaching_sectors(
        built.shape, built.e_value, made.radius_pc, built.edge_pc, figures=again))


def test_the_share_falls_with_the_census_dial_and_is_not_a_table_lookup():
    shares = [plan_outline.report(_outline(["--richness", value])).share for value in ("1.0", "3.0")]
    assert shares[0] < shares[1] and all(0.0 < share < 1.0 for share in shares)


def test_the_fixed_system_mass_and_the_lookup_table_are_gone():
    from planetgen.galaxy import dials
    for name in ("OUTLINE_DENSITY_PER_HEIGHT", "PLANNED_SHARE_BY_CREST", "planned_share", "disc_crest"):
        assert not hasattr(dials, name)
