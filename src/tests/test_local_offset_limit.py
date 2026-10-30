"""
TEST.200: the nightly update's local offset (docs/design/orbits/orbital-updates.md, "The local offset after a rotation"; docs/design/architecture/
numeric-precision-audit.md sections 3.1 and 5.1). The update turns an object's absolute galaxy point about the galactic axis and stores the offset
from its sector's centre as `point - centre` along the sector's axes, in milliparsecs. The offset is therefore only as good as the absolute point it
came from, the 64-bit spacing of the coordinate (28 km at the Sun, 898 km at the 262 kpc reach), 1e7 and 3e5 times below the smallest move the galaxy
stores (0.01 milliparsec). These tests pin that: the offset put back on its centre lands within one unit in the last place of the absolute point.
"""

import math
import random

import pytest

from planetgen.db import store
from planetgen.galaxy.frame import galaxy_to_local_pc, local_to_galaxy_pc

EDGE_PC = 8.0
"""float: The sector edge the offsets are drawn inside (the default edge is of this size)."""

RADII_PC = (8_000.0, 26_000.0, 262_144.0)
"""tuple: The Sun's distance from the centre, a far disc, and the 262 kpc reach."""

SMALLEST_STORED_MOVE_PC = 1e-5
"""float: 0.01 milliparsec, the smallest move the galaxy stores."""

TRIALS = 400


def _ulp_of(point):
    return max(math.ulp(abs(value)) for value in point)


def _centre(radius_pc, rng):
    colatitude = math.acos(rng.uniform(-1.0, 1.0))
    longitude = rng.uniform(0.0, 2.0 * math.pi)
    return (radius_pc * math.sin(colatitude) * math.cos(longitude), radius_pc * math.sin(colatitude) * math.sin(longitude),
            radius_pc * math.cos(colatitude))


@pytest.mark.parametrize("radius_pc", RADII_PC)
def test_local_offset_after_a_rotation_is_within_one_ulp_of_the_absolute_point(radius_pc):
    rng = random.Random(int(radius_pc))
    worst = 0.0
    angles = [0.0, 1e-9, math.pi / 2, math.pi, 2.0 * math.pi * (1.0 - 1e-12), store._galactic_turn(249_999_999.0, 0.25)]
    for trial in range(TRIALS):
        centre = _centre(radius_pc, rng)
        offset = tuple(rng.uniform(-EDGE_PC / 2, EDGE_PC / 2) for _ in range(3))
        angle = angles[trial % len(angles)] if trial < 3 * len(angles) else rng.uniform(0.0, 2.0 * math.pi)
        # The sector the object ends in may be another cell than it started in: its centre is the turned centre moved by up to a sector edge.
        landing = tuple(value + rng.uniform(-EDGE_PC, EDGE_PC) for value in store._rotate_about_axis(centre, angle))
        point = store._rotate_about_axis(local_to_galaxy_pc(centre, offset), angle)
        local_mpc = tuple(value * 1000.0 for value in galaxy_to_local_pc(landing, point))
        # The stored offset, put back on its centre, must give back the stored absolute point to the coordinate's last place.
        back = local_to_galaxy_pc(landing, tuple(value / 1000.0 for value in local_mpc))
        error = max(abs(back[i] - point[i]) for i in range(3))
        worst = max(worst, error / _ulp_of(point))
        assert error <= _ulp_of(point), (radius_pc, trial, angle, error, _ulp_of(point))
    assert worst <= 1.0


@pytest.mark.parametrize("radius_pc", RADII_PC)
def test_the_offset_error_is_far_below_the_smallest_move_the_galaxy_stores(radius_pc):
    """The limit in numbers: one unit in the last place at this distance, against the smallest stored move."""
    rounding = math.ulp(radius_pc) / 2.0   # the most a coordinate this size is wrong by when it is rounded (898 km at the reach)
    assert SMALLEST_STORED_MOVE_PC / rounding > 3e5, (radius_pc, rounding)
