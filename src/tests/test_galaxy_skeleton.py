# tests/test_galaxy_skeleton.py

"""
Tests for the density bound `planetgen.galaxy.skeleton` builds the outline
from. See that module's docstring for the reasoning: an exact upper bound over
`theta` falls monotonically with both `R` and `|z|`, so a shell's qualifying
bands are one run about the plane. The outline itself (`build_shell_extents`,
`GalaxyBounds`) is tested in `test_galaxy_shell_outline.py`.
"""

import math

import pytest

from planetgen.galaxy.density import halo_density, relative_density
from planetgen.galaxy.skeleton import bound_relative_density_at, expected_system_count_at_density_1
from tests.toy_shape import toy_shape

SHAPE = toy_shape(
    disk_scale_length_pc=40.0,
    disk_scale_height_pc=12.0,
    bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0,
    arm_count=2,
    pitch_angle_rad=math.radians(15),
    arm_amplitude=0.4,
)
EDGE_PC = 4.0


def test_expected_system_count_at_density_1_matches_spacesector():
    from planetgen.galaxy.sector import SpaceSector
    expected = SpaceSector(name="x", edge_ly=11.5).expected_system_count()
    assert expected_system_count_at_density_1(11.5) == pytest.approx(expected)


def _radius(shell):
    return shell * EDGE_PC


def test_bound_is_symmetric_in_z_and_falls_with_height():
    r = _radius(20)
    previous = bound_relative_density_at(SHAPE, r, 0.0)
    for z in (1.0, 5.0, 20.0, 60.0):
        assert bound_relative_density_at(SHAPE, r, z) == pytest.approx(bound_relative_density_at(SHAPE, r, -z))
        current = bound_relative_density_at(SHAPE, r, z)
        assert current < previous
        previous = current


def test_bound_falls_with_radius():
    previous = bound_relative_density_at(SHAPE, _radius(0), 8.0)
    for shell in range(1, 60):
        current = bound_relative_density_at(SHAPE, _radius(shell), 8.0)
        assert current < previous
        previous = current


def test_bound_is_never_below_the_exact_density():
    for shell in (0, 3, 15, 40):
        r = _radius(shell)
        for z in (0.0, 4.0, -9.0):
            bound = bound_relative_density_at(SHAPE, r, z)
            for k in range(0, 360, 15):
                theta = math.radians(k)
                point = (r * math.cos(theta), r * math.sin(theta), z)
                # The stellar halo (GEN.234) is left out of the bound on purpose, so it can't widen the outline.
                exact = relative_density(point, SHAPE) - halo_density(point, SHAPE)
                assert exact <= bound * (1 + 1e-12)
