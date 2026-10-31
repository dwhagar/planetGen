# tests/test_galaxy_acceleration.py

"""
GEN.342: the acceleration kernel. The pull of the influence set plus the galactic gradient equals a direct sum for random sets, momentum is conserved for an
isolated pair, an empty set leaves a star's orbit to the potential, and the pull is finite everywhere.
"""

import math

import numpy as np
import pytest

from planetgen.galaxy import acceleration as A
from planetgen.galaxy import influence
from planetgen.galaxy import potential as P
from planetgen.physics import constants

G = P.G_PC_KMS2_PER_MSUN
MSUN = constants.SOLAR_MASS_TO_KG


def _set(positions, masses_msun, monopoles=()):
    positions = np.asarray(positions, float).reshape(-1, 3)
    return influence.InfluenceSet(10.0, tuple(range(len(positions))), positions, np.asarray(masses_msun, float) * MSUN, tuple(monopoles), 1)


def _galaxy():
    return P.Potential(spheres=(P.Hernquist("bulge", 1.0e10, 600.0),), disks=(P.Disk("thin_disk", 5.0e10, 3000.0, 300.0),),
                       halo=P.Halo(8.0e11, 12.0, 220000.0), black_hole_solar=4.0e6)


def test_the_set_pull_is_a_direct_sum_for_random_sets():
    rng = np.random.default_rng(3)
    for count in (1, 5, 40):
        positions = rng.normal(0.0, 3.0, (count, 3))
        masses = rng.uniform(0.1, 30.0, count)
        here = rng.normal(0.0, 3.0, 3)
        expected = np.zeros(3)
        for p, m in zip(positions, masses):
            d = p - here
            expected += G * m * d / (d @ d + A.MEMBER_SOFTENING_PC ** 2) ** 1.5
        assert A.set_pull(here, _set(positions, masses)) == pytest.approx(tuple(expected), rel=1e-12)


def test_a_monopole_pulls_like_a_member_of_its_mass_at_its_centre_of_mass():
    monopole = influence.Monopole((1, 2, 3), (4.0, 0.0, 0.0), 12.0 * MSUN, 7)
    assert A.set_pull((0.0, 0.0, 0.0), _set([], [], [monopole])) == pytest.approx(A.set_pull((0.0, 0.0, 0.0), _set([(4.0, 0.0, 0.0)], [12.0])), rel=1e-12)


def test_the_kernel_is_the_potentials_pull_plus_the_sets():
    pot, here = _galaxy(), (6000.0, 1000.0, 20.0)
    members = _set([(6001.0, 1000.0, 20.0), (5990.0, 1004.0, 22.0)], [3.0, 0.5])
    total = A.pull(here, pot, members)
    expected = tuple(a + b for a, b in zip(P.acceleration(pot, *here), A.set_pull(here, members)))
    assert total == pytest.approx(expected, rel=1e-12)
    assert A.pull(here, pot) == P.acceleration(pot, *here) == A.pull(here, pot, _set([], []))


def test_an_isolated_pair_conserves_momentum():
    m1, m2 = 3.0, 11.0
    p1, p2 = np.array([0.0, 0.0, 0.0]), np.array([0.7, -0.2, 0.4])
    a1 = np.array(A.set_pull(p1, _set([p2], [m2])))
    a2 = np.array(A.set_pull(p2, _set([p1], [m1])))
    assert m1 * a1 + m2 * a2 == pytest.approx(np.zeros(3), abs=1e-15)


def test_an_empty_set_leaves_a_stars_orbit_to_the_potential():
    pot = _galaxy()
    speed = P.circular_speed_kms(pot, 8000.0)

    def run(kicker):
        position, velocity = (8000.0, 0.0, 0.0), (0.0, speed, 0.0)
        dt = 1.0
        for _ in range(500):
            velocity = kicker(position, velocity, 0.5 * dt)                       # kick-drift-kick
            position = tuple(p + v * dt * P.KMS_TO_PC_PER_MYR for p, v in zip(position, velocity))
            velocity = kicker(position, velocity, 0.5 * dt)
        return position, velocity

    empty = _set([], [])
    with_kernel = run(lambda x, v, h: A.velocity_update(x, v, pot, empty, h))
    alone = run(lambda x, v, h: A.kick(v, P.acceleration(pot, *x), h))
    assert with_kernel == alone
    assert math.hypot(*with_kernel[0][:2]) == pytest.approx(8000.0, abs=1.0)


def test_a_planets_pull_is_the_perturbation_in_its_hosts_frame():
    pot = _galaxy()
    host, planet = (6000.0, 0.0, 0.0), (6000.0, 0.0001, 0.0)
    far = _set([(6004.0, 3.0, 1.0)], [2.0])
    direct = A.host_frame_pull(planet, host, pot, far)
    assert direct == pytest.approx(tuple(a - b for a, b in zip(A.pull(planet, pot, far), A.pull(host, pot, far))), rel=1e-12)
    assert max(abs(c) for c in A.host_frame_pull(host, host, pot, far)) == 0.0
    # a planet 0.0001 pc from its host feels only the tide: far smaller than the pull at either place
    assert max(abs(c) for c in direct) < 1e-3 * max(abs(c) for c in A.pull(host, pot, far))


def test_the_kick_is_acceleration_times_the_step_in_the_codes_units():
    assert A.kick((10.0, 0.0, -5.0), (1.0e-3, 0.0, 2.0e-3), 0.0) == (10.0, 0.0, -5.0)
    one = A.kick((0.0, 0.0, 0.0), (1.0e-3, 0.0, 0.0), 1.0)[0]
    assert one == pytest.approx(1.0e-3 * P.KMS_TO_PC_PER_MYR)
    assert A.kick((0.0, 0.0, 0.0), (1.0e-3, 0.0, 0.0), 2.0)[0] == pytest.approx(2.0 * one)


@pytest.mark.parametrize("here", [(0.0, 0.0, 0.0), (5.0, 5.0, 5.0), (1.0e7, -2.0e7, 3.0e6), (0.001, 0.0, 0.0)])
def test_the_pull_is_finite_everywhere_even_at_a_members_own_point(here):
    pot = _galaxy()
    members = _set([here, (5.0, 5.0, 5.0)], [1.0, 1.0])
    assert all(math.isfinite(c) for c in A.pull(here, pot, members))
    assert A.set_pull(here, _set([here], [1.0])) == (0.0, 0.0, 0.0)
