# tests/test_gen375_field_origin.py

"""
GEN.375: the field kernel measures from a local origin, so a sector 262 kpc from the galaxy's center gets the same path and the same
pull as one beside it.

A float64 near the 262 kpc reach is 1,049 km apart (`docs/design/architecture/numeric-precision-audit.md` section 3.4). The pull of
a mass is checked against exact decimal arithmetic on the same numbers, and a whole path is checked against the path of the same
bodies set down beside the origin: the far knots, less the origin, must agree to the spacing of the far numbers.
"""

from decimal import Decimal, getcontext

import numpy as np
import pytest

from planetgen.physics import constants, field_kernel
from planetgen.physics import sector_path as sp

PC = constants.PARSEC_M
MSUN = constants.SOLAR_MASS_TO_KG
EDGE = 20.0 * PC
FAR = (262_000.0 * PC, -91_000.0 * PC, 41_000.0 * PC)
"""tuple: A sector close to the 262,144 pc reach, metres."""

OFFSETS = (((1.5 * PC, -2.0 * PC, 0.4 * PC), 1.2), ((-3.0 * PC, 0.7 * PC, -1.1 * PC), 0.6), ((0.2 * PC, 4.5 * PC, 2.0 * PC), 3.0),
           ((0.8 * PC, 0.35 * PC, 0.0), 40.0))
"""tuple: `(offset from the origin, mass in suns)` of four masses, the last one heavy and close to the paths below."""


def _far_masses():
    return [sp.PointMass(tuple(float(f + o) for f, o in zip(FAR, offset)), mass * MSUN) for offset, mass in OFFSETS]


def _local(masses):
    """The same masses seen from `FAR`: the difference of two nearby numbers, which is exact."""
    return [sp.PointMass(tuple(a - b for a, b in zip(m.position, FAR)), m.mass_kg) for m in masses]


def _far_box():
    low = tuple(f - EDGE / 2.0 for f in FAR)
    high = tuple(f + EDGE / 2.0 for f in FAR)
    return sp.box_cell(low, high), sp.box_cell(tuple(a - b for a, b in zip(low, FAR)), tuple(a - b for a, b in zip(high, FAR)))


def _decimal_pull(point, masses, softening):
    getcontext().prec = 60
    px, py, pz = (Decimal(c) for c in point)
    pull, potential = [Decimal(0)] * 3, Decimal(0)
    for m in masses:
        dx, dy, dz = (Decimal(c) - p for c, p in zip(m.position, (px, py, pz)))
        r2 = dx * dx + dy * dy + dz * dz + Decimal(softening) ** 2
        r = r2.sqrt()
        k = Decimal(m.mu) / (r2 * r)
        pull = [pull[0] + k * dx, pull[1] + k * dy, pull[2] + k * dz]
        potential -= Decimal(m.mu) / r
    return [float(c) for c in pull], float(potential)


def test_the_pull_and_potential_far_from_the_origin_match_exact_arithmetic():
    masses = _far_masses()
    table = sp.MassTable(masses)
    point = tuple(float(f + o) for f, o in zip(FAR, (0.9 * PC, 0.3 * PC, -0.6 * PC)))
    exact_pull, exact_potential = _decimal_pull(point, masses, 0.01 * PC)
    pull = field_kernel.accelerations(point, table, 0.01 * PC)[0]
    potential = field_kernel.potentials(point, table, 0.01 * PC)[0]
    assert pull == pytest.approx(exact_pull, rel=1e-13)
    assert potential == pytest.approx(exact_potential, rel=1e-13)


def test_the_pull_does_not_depend_on_which_point_is_first():
    table = sp.MassTable(_far_masses())
    points = np.array([[f + o for f, o in zip(FAR, offset)] for offset in ((0.9 * PC, 0.3 * PC, -0.6 * PC), (-2.0 * PC, 1.0 * PC, 0.1 * PC))])
    forward = field_kernel.accelerations(points, table, 0.01 * PC)
    backward = field_kernel.accelerations(points[::-1], table, 0.01 * PC)[::-1]
    assert forward == pytest.approx(backward, rel=1e-12)


def test_the_smooth_term_still_gets_absolute_positions():
    table = sp.MassTable(_far_masses())
    seen = []

    def smooth(points):
        seen.append(points.copy())
        return np.zeros((points.shape[0], 3))

    point = (FAR[0] + PC, FAR[1], FAR[2])
    field_kernel.accelerations(point, table, 0.01 * PC, smooth=smooth)
    assert seen[0][0, 0] == point[0]


def _integrate(masses, cell, start, velocity):
    return sp.integrate_path(start, velocity, masses, cell, EDGE, softening_m=0.05 * PC, exclude_key=None)


def test_a_path_far_from_the_origin_is_the_path_beside_it_to_the_spacing_of_the_far_numbers():
    far_masses = _far_masses()
    far_cell, near_cell = _far_box()
    velocity = (2.0e3, 0.1e3, -0.03e3)
    far_start = tuple(float(f + o) for f, o in zip(FAR, (-9.0 * PC, 0.0, 0.0)))
    near_start = tuple(a - b for a, b in zip(far_start, FAR))
    far = _integrate(far_masses, far_cell, far_start, velocity)
    near = _integrate(_local(far_masses), near_cell, near_start, velocity)
    assert far.exited and near.exited
    assert far.steps == near.steps and len(far.knots) == len(near.knots) > 2, (far.steps, len(far.knots))
    spacing = float(np.spacing(max(abs(c) for c in FAR))) * 2.0
    for far_knot, near_knot in zip(far.knots, near.knots):
        # The exit knot is found by asking whether a point is in the cell, and that question is asked of absolute metres, so it
        # can differ by half a spacing of the far numbers; every other knot is the same to the last bits.
        assert far_knot.t_s == pytest.approx(near_knot.t_s, rel=1e-11)
        assert far_knot.velocity == pytest.approx(near_knot.velocity, rel=1e-9)
        for axis in range(3):
            assert far_knot.position[axis] - FAR[axis] == pytest.approx(near_knot.position[axis], abs=spacing, rel=0.0)


def test_the_path_starts_exactly_where_it_was_given_and_bends_the_same_way_at_the_reach():
    far_masses = _far_masses()
    far_cell, _near = _far_box()
    start = tuple(float(f + o) for f, o in zip(FAR, (-9.0 * PC, 0.0, 0.0)))
    path = _integrate(far_masses, far_cell, start, (2.0e3, 0.1e3, -0.03e3))
    assert path.start.position == start
    assert abs(path.end.velocity[1]) > 100.0 or abs(path.end.velocity[0] - 2.0e3) > 100.0   # the heavy mass bent it
