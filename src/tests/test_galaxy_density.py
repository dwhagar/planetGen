# tests/test_galaxy_density.py

"""
Tests for `planetgen.galaxy.density` -- the exponential-disk-plus-
bulge-plus-spiral-arm density model from
`docs/design/galaxy/galaxy-disk-density.md` revision 2.
"""

import math

import pytest
from planetgen.util import draw
from planetgen import tuning

from planetgen.galaxy import density
from planetgen.physics import constants as physicalConstants
from planetgen.physics.units import ly_to_pc
from planetgen.galaxy.density import (
    build_galaxy_shape,
    predicted_star_count,
    relative_density,
)

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0,
    disk_scale_height_pc=12.0,
    bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0,
    arm_count=2,
    pitch_angle_rad=math.radians(15),
    arm_amplitude=0.4,
)


def _systems_per_raw_unit(shape):
    return shape.k_norm * physicalConstants.LOCAL_STELLAR_DENSITY_LY3 / ly_to_pc(1.0) ** 3


def test_the_galaxy_holds_the_systems_its_stellar_mass_makes():
    """GEN.404: no point of the galaxy is pinned to 1.0; the whole holds its stellar mass at the mean system mass."""
    systems = density.total_mass(SHAPE) * _systems_per_raw_unit(SHAPE)
    assert systems == pytest.approx(SHAPE.structure.stellar_mass_solar / density.mean_system_mass_solar(SHAPE), rel=1e-12)


def test_a_more_massive_galaxy_of_the_same_shape_is_proportionally_denser():
    heavy = build_galaxy_shape(
        disk_scale_length_pc=SHAPE.disk_scale_length_pc, disk_scale_height_pc=SHAPE.disk_scale_height_pc,
        bulge_scale_radius_pc=SHAPE.bulge_scale_radius_pc, bulge_amplitude=SHAPE.bulge_amplitude,
        arm_count=SHAPE.arm_count, pitch_angle_rad=SHAPE.pitch_angle_rad, arm_amplitude=SHAPE.arm_amplitude,
        structure=SHAPE.structure._replace(stellar_mass_solar=2.0 * SHAPE.structure.stellar_mass_solar))
    assert heavy.k_norm / SHAPE.k_norm == pytest.approx(2.0, rel=1e-6)
    assert relative_density((9000.0, 0.0, 0.0), heavy) == pytest.approx(2.0 * relative_density((9000.0, 0.0, 0.0), SHAPE), rel=1e-6)


def test_relative_density_peaks_at_galactic_center():
    center = relative_density((0.0, 0.0, 0.0), SHAPE)
    nearby = relative_density((5.0, 0.0, 0.0), SHAPE)
    far = relative_density((500.0, 0.0, 0.0), SHAPE)
    assert center > nearby > far


def test_relative_density_falls_off_with_height_above_plane():
    in_plane = relative_density((60.0, 0.0, 0.0), SHAPE)
    above_plane = relative_density((60.0, 0.0, 40.0), SHAPE)
    assert in_plane > above_plane > 0


def test_relative_density_arm_crest_exceeds_inter_arm_at_same_radius():
    r_cyl = 60.0
    theta_arm = math.log(r_cyl / SHAPE.disk_scale_length_pc) / math.tan(SHAPE.pitch_angle_rad)
    on_arm = relative_density((r_cyl * math.cos(theta_arm), r_cyl * math.sin(theta_arm), 0.0), SHAPE)
    theta_inter = theta_arm + math.pi / SHAPE.arm_count
    inter_arm = relative_density((r_cyl * math.cos(theta_inter), r_cyl * math.sin(theta_inter), 0.0), SHAPE)
    assert on_arm > inter_arm


def test_relative_density_is_never_negative():
    draw.set_run_seed(0)
    for _ in range(200):
        position = (draw.uniform(-200, 200), draw.uniform(-200, 200), draw.uniform(-100, 100))
        assert relative_density(position, SHAPE) >= 0.0


def test_predicted_star_count_scales_linearly_with_expected_count():
    position = (10.0, 0.0, 0.0)
    a = predicted_star_count(position, SHAPE, expected_system_count_at_density_1=1.0)
    b = predicted_star_count(position, SHAPE, expected_system_count_at_density_1=3.0)
    assert b == pytest.approx(3.0 * a)


def test_relative_density_does_not_overflow_far_off_the_galactic_plane():
    # Regression test: 1/cosh(x)**2 computed as `1.0 / math.cosh(x)**2`
    # raises OverflowError once `|z| / disk_scale_height_pc` exceeds
    # ~710 (cosh(x) itself overflows there) -- easily reached with a
    # small disk_scale_height_pc relative to the position being
    # evaluated (this SHAPE's disk_scale_height_pc=12 overflows at
    # |z| > ~8520 pc under the old formula). The true value is 0.0 in
    # this limit, not a crash.
    position = (0.0, 0.0, 100000.0)
    assert relative_density(position, SHAPE) >= 0.0
    assert relative_density(position, SHAPE) == pytest.approx(
        relative_density((0.0, 0.0, -100000.0), SHAPE)
    )


def test_a_toy_galaxy_is_normalized_to_its_own_mass_too():
    toy = build_galaxy_shape(
        disk_scale_length_pc=40.0, disk_scale_height_pc=12.0,
        bulge_scale_radius_pc=10.0, bulge_amplitude=2.0,
        arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
    )
    assert density.total_mass(toy) * _systems_per_raw_unit(toy) == pytest.approx(
        toy.structure.stellar_mass_solar / density.mean_system_mass_solar(toy), rel=1e-12)
