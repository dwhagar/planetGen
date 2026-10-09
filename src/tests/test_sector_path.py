# tests/test_sector_path.py

"""GEN.123: the path of a body through a sector, as cubic Hermite spline knots bent by the sector's masses."""

import math

import pytest

from planetgen.galaxy import geometry
from planetgen.physics import constants
from planetgen.physics import sector_path as sp

PC = constants.PARSEC_M
MSUN = constants.SOLAR_MASS_TO_KG
EDGE = 20.0 * PC


def _box(half=EDGE / 2.0):
    return lambda p: all(abs(c) <= half for c in p)


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


def test_paths_chain_across_sectors_from_the_exit_state():
    first = sp.integrate_path((-9.0 * PC, 0.0, 0.0), (50.0e3, 5.0e3, 0.0), [], _box(), EDGE)
    entry = first.restart()
    assert entry.t_s == 0.0 and entry.position == first.end.position and entry.velocity == first.end.velocity
    shifted = tuple(c - EDGE for c in entry.position)  # the neighbour cell's own box, centred one edge along x
    second = sp.integrate_path(entry.position, entry.velocity, [], lambda p: _box()((p[0] - EDGE, p[1], p[2])), EDGE)
    assert second.exited and shifted[0] == pytest.approx(-EDGE / 2.0, abs=1e-3 * PC)
    assert second.end.position[0] == pytest.approx(1.5 * EDGE, abs=1e-6 * PC)


def test_sector_inside_follows_the_grid_cell():
    edge_pc = 4.0
    address = (1500, 2, 700)
    centre_pc = geometry.sector_position_pc(*address, edge_pc)
    inside = sp.sector_inside(address, edge_pc)
    assert inside(tuple(c * PC for c in centre_pc))
    far = tuple(c * PC for c in geometry.sector_position_pc(1500, 2, 703, edge_pc))
    assert not inside(far)
