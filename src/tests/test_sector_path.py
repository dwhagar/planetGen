# tests/test_sector_path.py

"""GEN.123: the path of a body through a sector, as cubic Hermite spline knots bent by the sector's masses."""

import math

import pytest

from planetgen.galaxy import spherical
from planetgen.physics import constants
from planetgen.physics import sector_path as sp

PC = constants.PARSEC_M
MSUN = constants.SOLAR_MASS_TO_KG
EDGE = 20.0 * PC


def _box(half=EDGE / 2.0, center=(0.0, 0.0, 0.0)):
    return sp.box_cell(tuple(c - half for c in center), tuple(c + half for c in center))


def _angle_between(a, b):
    cos = sum(x * y for x, y in zip(a, b)) / (math.sqrt(sum(x * x for x in a)) * math.sqrt(sum(y * y for y in b)))
    return math.acos(max(-1.0, min(1.0, cos)))


def test_hermite_reproduces_a_cubic_exactly():
    def at(t):
        return ((1.0 + 2.0 * t - 3.0 * t ** 2 + t ** 3), 0.5 * t ** 2, -t ** 3 + 4.0)

    def slope(t):
        return ((2.0 - 6.0 * t + 3.0 * t ** 2), t, -3.0 * t ** 2)

    a, b = sp.PathKnot(0.5, at(0.5), slope(0.5)), sp.PathKnot(2.5, at(2.5), slope(2.5))
    for t in (0.5, 1.0, 1.7, 2.5):
        position, velocity = sp.hermite(a, b, t)
        assert position == pytest.approx(at(t), abs=1e-12)
        assert velocity == pytest.approx(slope(t), abs=1e-12)


def test_with_no_masses_the_path_is_two_knots_on_a_straight_line_to_the_edge():
    velocity = (100.0e3, 0.0, 0.0)
    path = sp.integrate_path((-9.0 * PC, 1.0 * PC, 0.0), velocity, [], _box(), EDGE)
    assert path.exited and len(path.knots) == 2
    assert path.end.position[0] == pytest.approx(EDGE / 2.0, abs=1e-6 * PC)
    assert path.end.position[1] == pytest.approx(1.0 * PC)
    assert path.end.velocity == pytest.approx(velocity)
    assert path.duration_s == pytest.approx(19.0 * PC / 100.0e3, rel=1e-9)
    midpoint = path.position_at(path.duration_s / 2.0)
    assert midpoint[0] == pytest.approx(0.5 * (-9.0 + 10.0) * PC, abs=1e-6 * PC)


def test_a_heavy_mass_bends_the_path_like_a_hyperbola():
    mu_mass = 1.0e6 * MSUN
    b, v = 3.0 * PC, 200.0e3
    mass = sp.PointMass((0.0, 0.0, 0.0), mu_mass)
    edge = 400.0 * PC
    path = sp.integrate_path((-edge / 2.0 + PC, b, 0.0), (v, 0.0, 0.0), [mass], _box(edge / 2.0), edge,
                             softening_m=1.0e-3 * PC)
    assert path.exited and 2 < len(path.knots) <= sp.DEFAULT_MAX_KNOTS
    mu = constants.G * mu_mass
    eccentricity = math.sqrt(1.0 + (b * v * v / mu) ** 2)
    turn = 2.0 * math.asin(1.0 / eccentricity)
    assert _angle_between((v, 0.0, 0.0), path.end.velocity) == pytest.approx(turn, rel=0.02)
    assert path.end.velocity[1] < 0.0  # pulled toward the mass
    speed_out = math.sqrt(sum(c * c for c in path.end.velocity))
    assert speed_out == pytest.approx(v, rel=2e-3)  # energy is kept


def test_the_spline_stays_within_the_tolerance_of_a_fine_integration():
    mass = sp.PointMass((0.0, 0.0, 0.0), 1.0e6 * MSUN)
    edge = 400.0 * PC
    args = ((-edge / 2.0 + PC, 3.0 * PC, 0.0), (200.0e3, 0.0, 0.0), [mass], _box(edge / 2.0), edge)
    coarse = sp.integrate_path(*args, softening_m=1.0e-3 * PC)
    fine = sp.integrate_path(*args, softening_m=1.0e-3 * PC, tolerance_m=1.0e-4 * PC, max_knots=400)
    assert len(fine.knots) > len(coarse.knots)
    for k in range(0, 41):
        t = fine.duration_s * k / 40.0
        error = math.dist(coarse.position_at(t * coarse.duration_s / fine.duration_s),
                          fine.position_at(t))
        assert error < 4.0 * sp.DEFAULT_TOLERANCE_FRACTION * edge  # same path, time slightly rescaled by the cut


def test_the_knot_count_is_capped():
    mass = sp.PointMass((0.0, 0.0, 0.0), 1.0e6 * MSUN)
    edge = 400.0 * PC
    path = sp.integrate_path((-edge / 2.0 + PC, 3.0 * PC, 0.0), (200.0e3, 0.0, 0.0), [mass], _box(edge / 2.0), edge,
                             softening_m=1.0e-3 * PC, tolerance_m=1.0, max_knots=7)
    assert len(path.knots) == 7


def test_masses_too_weak_or_far_to_matter_are_skipped_and_the_strongest_come_first():
    edge = 4.0 * PC
    here, heading = (0.0, 0.0, 0.0), (100.0e3, 0.0, 0.0)
    faint = sp.PointMass((0.0, 3.0 * edge, 0.0), 1.0 * MSUN)
    near = sp.PointMass((edge, 0.3 * PC, 0.0), 1.0e4 * MSUN)
    heavy = sp.PointMass((edge, -1.0 * PC, 0.0), 1.0e6 * MSUN)
    kept = sp.relevant_masses(here, heading, [faint, near, heavy], edge)
    assert faint not in kept
    assert kept[0] is heavy and near in kept
    assert sp.relevant_masses(here, heading, [near, heavy], edge, limit=1) == [heavy]
    assert sp.relevant_masses(here, (0.0, 0.0, 0.0), [heavy], edge) == []


def test_a_body_at_rest_has_one_knot_and_does_not_exit():
    path = sp.integrate_path((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), [], _box(), EDGE)
    assert len(path.knots) == 1 and not path.exited
    assert path.position_at(123.0) == (0.0, 0.0, 0.0)


def test_a_path_must_start_inside_and_needs_two_knots_at_least():
    with pytest.raises(ValueError):
        sp.integrate_path((EDGE, 0.0, 0.0), (1.0, 0.0, 0.0), [], _box(), EDGE)
    with pytest.raises(ValueError):
        sp.integrate_path((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), [], _box(), EDGE, max_knots=1)


def test_a_body_held_by_a_mass_is_cut_off_instead_of_running_forever():
    mass = sp.PointMass((0.0, 0.0, 0.0), 1.0e6 * MSUN)
    edge = 400.0 * PC
    path = sp.integrate_path((3.0 * PC, 0.0, 0.0), (1.0e3, 0.0, 0.0), [mass], _box(edge / 2.0), edge, max_steps=300)
    assert not path.exited and 2 <= len(path.knots) <= sp.DEFAULT_MAX_KNOTS


def test_a_body_held_by_a_heavy_mass_costs_one_revolution_not_the_step_limit():
    """GEN.232: thousands of bodies orbit the core's black hole; each used to run every step of the limit."""
    mass = sp.PointMass((0.0, 0.0, 0.0), 1.0e8 * MSUN)
    edge = 4.0 * PC
    speed = 0.4 * math.sqrt(constants.G * mass.mass_kg / PC)
    path = sp.integrate_path((PC, 0.0, 0.0), (0.0, speed, 0.0), [mass], sp.ball_cell((0.0, 0.0, 0.0), 2.5 * PC), edge)
    assert not path.exited and 2 < len(path.knots) <= sp.DEFAULT_MAX_KNOTS
    assert path.steps < sp.DEFAULT_MAX_STEPS // 2, path.steps


def test_paths_chain_across_sectors_from_the_exit_state():
    first = sp.integrate_path((-9.0 * PC, 0.0, 0.0), (50.0e3, 5.0e3, 0.0), [], _box(), EDGE)
    entry = first.restart()
    assert entry.t_s == 0.0 and entry.position == first.end.position and entry.velocity == first.end.velocity
    shifted = tuple(c - EDGE for c in entry.position)  # the neighbour cell's own box, centred one edge along x
    second = sp.integrate_path(entry.position, entry.velocity, [], _box(center=(EDGE, 0.0, 0.0)), EDGE)
    assert second.exited and shifted[0] == pytest.approx(-EDGE / 2.0, abs=1e-3 * PC)
    assert second.end.position[0] == pytest.approx(1.5 * EDGE, abs=1e-6 * PC)


def test_sector_inside_follows_the_grid_cell():
    edge_pc = 4.0
    address = (1500, 2, 700)
    centre_pc = spherical.sector_position_pc(*address, edge_pc)
    inside = sp.sector_inside(address, edge_pc)
    assert inside(tuple(c * PC for c in centre_pc))
    far = tuple(c * PC for c in spherical.sector_position_pc(1500, 2, 703, edge_pc))
    assert not inside(far)


def test_sector_inside_agrees_with_the_grid_everywhere_including_the_faces():
    """The compiled cell test is `spherical.sector_address_at`, for the core, a cap and an interior band, with points
    scattered over the cell and just either side of each face."""
    import random

    rnd = random.Random(11)
    edge_pc = 4.0
    for address in ((0, 0, 0), (1, 1, 0), (1, -1, 0), (2, 0, 1), (5, 2, 7), (40, 3, 11), (1500, 2, 700), (1500, 2, 0)):
        inside = sp.sector_inside(address, edge_pc)
        centre = spherical.sector_position_pc(*address, edge_pc)
        points = [tuple(c + rnd.uniform(-1.5, 1.5) * edge_pc for c in centre) for _ in range(400)]
        for shell_delta in (-1e-9, 1e-9):   # a hair inside and outside the cell's faces, along its own radius
            radius = math.sqrt(sum(c * c for c in centre)) or edge_pc
            points.append(tuple(c * (1 + shell_delta) for c in centre))
            points.append(tuple(c * (1 + shell_delta) * (radius + edge_pc / 2) / radius for c in centre))
        for point in points:
            wanted = spherical.sector_address_at(point, edge_pc) == address
            assert inside(tuple(c * PC for c in point)) == wanted, (address, point)


def test_the_field_kernel_equals_a_plain_numpy_sum():
    import numpy as np

    from planetgen.physics import field_kernel

    rnd = np.random.default_rng(3)
    masses = sp.MassTable([sp.PointMass(tuple(rnd.uniform(-5, 5, 3) * PC), 10 ** rnd.uniform(1, 8) * MSUN) for _ in range(40)])
    points = rnd.uniform(-6, 6, (25, 3)) * PC
    softening = 0.3 * PC
    offsets = masses.positions[None, :, :] - points[:, None, :]
    r2 = (offsets ** 2).sum(axis=2) + softening ** 2
    expected_a = ((masses.mus[None, :] / r2 ** 1.5)[:, :, None] * offsets).sum(axis=1)
    expected_phi = -(masses.mus[None, :] / np.sqrt(r2)).sum(axis=1)
    assert np.allclose(field_kernel.accelerations(points, masses, softening), expected_a, rtol=1e-12, atol=0)
    assert np.allclose(field_kernel.potentials(points, masses, softening), expected_phi, rtol=1e-12, atol=0)
    # The galaxy's smooth term (GEN.115) is added on top.
    bump = field_kernel.accelerations(points, masses, softening, smooth=lambda p: np.full((len(p), 3), 2.0))
    assert np.allclose(bump - expected_a, 2.0, rtol=1e-9)
    # A single point is one row.
    assert field_kernel.accelerations(points[0], masses, softening).shape == (1, 3)
    with pytest.raises(ValueError):
        field_kernel.accelerations(np.zeros((3, 2)), masses, softening)


def test_the_body_own_mass_does_not_pull_it_and_masses_far_off_the_path_are_left_out():
    near = sp.PointMass((0.0, 5.0 * PC, 0.0), 1.0e6 * MSUN, ("star_systems", 1))
    far = sp.PointMass((0.0, 1.0e6 * PC, 0.0), 1.0 * MSUN, ("star_systems", 2))
    with_own = sp.integrate_path((-9.0 * PC, 0.0, 0.0), (1.0e5, 0.0, 0.0), [near, far], _box(), EDGE)
    without = sp.integrate_path((-9.0 * PC, 0.0, 0.0), (1.0e5, 0.0, 0.0), [near, far], _box(), EDGE,
                                exclude_key=("star_systems", 1))
    assert with_own.end.velocity[1] > 0.0
    assert without.end.velocity[1] == 0.0 and len(without.knots) == 2
