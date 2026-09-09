# tests/test_galaxy_density.py

"""
Tests for `stellarObjects.galaxyDensity` -- the exponential-disk-plus-
bulge-plus-spiral-arm density model from
`docs/design/galaxy-disk-density.md` revision 2.
"""

import math

import pytest

from stellarObjects.galaxyDensity import (
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


def test_relative_density_is_exactly_one_at_calibration_point():
    calibration_radius_pc = 2.82 * SHAPE.disk_scale_length_pc
    theta_arm = math.log(calibration_radius_pc / SHAPE.disk_scale_length_pc) / math.tan(SHAPE.pitch_angle_rad)
    theta_interarm = theta_arm + math.pi / SHAPE.arm_count
    position = (
        calibration_radius_pc * math.cos(theta_interarm),
        calibration_radius_pc * math.sin(theta_interarm),
        0.0,
    )
    assert relative_density(position, SHAPE) == pytest.approx(1.0)


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
    import random
    random.seed(0)
    for _ in range(200):
        position = (random.uniform(-200, 200), random.uniform(-200, 200), random.uniform(-100, 100))
        assert relative_density(position, SHAPE) >= 0.0


def test_predicted_star_count_scales_linearly_with_expected_count():
    position = (10.0, 0.0, 0.0)
    a = predicted_star_count(position, SHAPE, expected_system_count_at_density_1=1.0)
    b = predicted_star_count(position, SHAPE, expected_system_count_at_density_1=3.0)
    assert b == pytest.approx(3.0 * a)


def test_build_galaxy_shape_default_calibration_radius():
    shape_default = build_galaxy_shape(
        disk_scale_length_pc=40.0, disk_scale_height_pc=12.0,
        bulge_scale_radius_pc=10.0, bulge_amplitude=2.0,
        arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
    )
    shape_explicit = build_galaxy_shape(
        disk_scale_length_pc=40.0, disk_scale_height_pc=12.0,
        bulge_scale_radius_pc=10.0, bulge_amplitude=2.0,
        arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
        calibration_radius_pc=2.82 * 40.0,
    )
    assert shape_default.k_norm == pytest.approx(shape_explicit.k_norm)
