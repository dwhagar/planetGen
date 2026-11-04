# tests/test_galaxy_stepper.py

"""
GEN.340: the star stepper. Kick-drift-kick in the potential keeps a circular orbit's radius and angular momentum over 5 Gyr, a century in one step agrees with a
hundred one-year steps, the stars of a batch step as if alone, the step cap sets a flag, and the extra-pull hook is added to the same kick.
"""

import math

import numpy as np
import pytest

from planetgen.galaxy import potential as P
from planetgen.galaxy import stepper as S


def _galaxy():
    return P.Potential(spheres=(P.Hernquist("bulge", 1.0e10, 600.0),), disks=(P.Disk("thin_disk", 5.0e10, 3000.0, 300.0),),
                       halo=P.Halo(8.0e11, 12.0, 220000.0), black_hole_solar=4.0e6)


def _circular(pot, radius, count=1):
    speed = P.circular_speed_kms(pot, radius)
    return np.tile([radius, 0.0, 0.0], (count, 1)), np.tile([0.0, speed, 0.0], (count, 1))


def test_a_star_at_the_suns_kind_of_radius_keeps_its_radius_and_angular_momentum_over_five_gyr():
    pot = _galaxy()
    x, v = _circular(pot, 8000.0)
    lz0 = x[0, 0] * v[0, 1] - x[0, 1] * v[0, 0]
    x, v, capped = S.advance(x, v, 5000.0, pot)
    assert not capped.any()
    assert abs(math.hypot(x[0, 0], x[0, 1]) - 8000.0) < 1.0
    assert abs(x[0, 2]) < 1e-6
    assert abs((x[0, 0] * v[0, 1] - x[0, 1] * v[0, 0]) / lz0 - 1.0) < 1e-10


def test_a_century_in_one_step_agrees_with_a_hundred_one_year_steps():
    pot = _galaxy()
    x0, v0 = np.array([[5000.0, 2000.0, 300.0], [800.0, -400.0, 50.0]]), np.array([[-40.0, 180.0, 10.0], [90.0, 120.0, -20.0]])
    year = 1.0e-6
    one, one_v, _ = S.advance(x0, v0, 100 * year, pot)
    x, v = x0, v0
    for _ in range(100):
        x, v, _ = S.advance(x, v, year, pot)
    assert np.abs(one - x).max() < 1e-6 and np.abs(one_v - v).max() < 1e-9


def test_a_batch_steps_each_star_as_if_it_were_alone():
    pot = _galaxy()
    rng = np.random.default_rng(1)
    x0 = rng.normal(0.0, 4000.0, (6, 3))
    v0 = rng.normal(0.0, 120.0, (6, 3))
    both, both_v, _ = S.advance(x0, v0, 50.0, pot)
    for i in range(6):
        alone, alone_v, _ = S.advance(x0[i:i + 1], v0[i:i + 1], 50.0, pot)
        assert alone[0] == pytest.approx(both[i], abs=1e-9) and alone_v[0] == pytest.approx(both_v[i], abs=1e-12)


def test_the_step_is_a_twentieth_of_a_crossing_and_never_over_the_cap():
    pot = _galaxy()
    x, v = _circular(pot, 8000.0, count=2)
    x[1] *= 0.01
    crossing = S.crossing_time_myr(x, v, pot)
    wanted = S.step_myr(x, v, pot)
    assert wanted[0] == pytest.approx(min(S.STEP_ETA * crossing[0], S.STEP_CAP_MYR)) and wanted[1] < wanted[0]
    faraway = S.step_myr(np.array([[2.0e6, 0.0, 0.0]]), np.array([[1.0, 0.0, 0.0]]), pot)
    assert faraway[0] == S.STEP_CAP_MYR


def test_block_levels_cut_the_span_in_powers_of_two_and_flag_the_cap():
    levels, capped = S.block_levels(np.array([10.0, 1.0, 0.3, 1.0e-9]), 8.0, max_level=6)
    assert levels.tolist() == [0, 3, 5, 6] and capped.tolist() == [False, False, False, True]
    pot = _galaxy()
    x, v = _circular(pot, 8000.0)
    _, _, capped = S.advance(x, v, 5000.0, pot, max_level=3)
    assert capped.all()


def test_the_extra_pull_is_added_to_the_same_kick():
    pot = P.Potential()
    x0, v0 = np.zeros((1, 3)), np.zeros((1, 3))
    push = 1.0e-4                                                    # (km/s)^2 per pc along x
    x, v, _ = S.advance(x0, v0, 2.0, pot, extra=lambda pos, vel, t: np.tile([push, 0.0, 0.0], (len(pos), 1)))
    t = 2.0 * P.KMS_TO_PC_PER_MYR
    assert v[0, 0] == pytest.approx(push * t, rel=1e-9) and x[0, 0] == pytest.approx(0.5 * push * t * t, rel=1e-9)
    still, _, _ = S.advance(x0, v0, 2.0, pot, extra=None)
    assert (still == 0.0).all()


def test_a_star_at_rest_at_the_centre_and_a_zero_span_are_safe():
    pot = _galaxy()
    x, v, capped = S.advance(np.zeros((1, 3)), np.zeros((1, 3)), 10.0, pot)
    assert np.isfinite(x).all() and np.isfinite(v).all() and not capped.any()
    keep, keep_v, _ = S.advance(np.array([[10.0, 0.0, 0.0]]), np.array([[1.0, 2.0, 3.0]]), 0.0, pot)
    assert keep.tolist() == [[10.0, 0.0, 0.0]] and keep_v.tolist() == [[1.0, 2.0, 3.0]]
    with pytest.raises(ValueError):
        S.advance(keep, keep_v, -1.0, pot)
    empty = S.advance(np.zeros((0, 3)), np.zeros((0, 3)), 5.0, pot)
    assert empty[0].shape == (0, 3)


def test_an_eccentric_orbit_keeps_its_energy_to_a_part_in_a_thousand_a_gyr():
    pot = _galaxy()
    x0, v0 = np.array([[9000.0, 0.0, 400.0]]), np.array([[0.0, 90.0, 15.0]])
    energy = lambda x, v: 0.5 * float(v[0] @ v[0]) + P.potential(pot, *x[0])
    e0 = energy(x0, v0)
    x, v, capped = S.advance(x0, v0, 1000.0, pot)
    assert not capped.any() and abs(energy(x, v) - e0) / abs(e0) < 1e-3


def test_each_star_of_a_batch_takes_its_own_span():
    """The rows of an orbital update have each their own epoch: an array of spans moves every star as if it were stepped alone over its span."""
    pot = _galaxy()
    rng = np.random.default_rng(2)
    x0 = rng.normal(0.0, 4000.0, (5, 3))
    v0 = rng.normal(0.0, 120.0, (5, 3))
    spans = np.array([0.0, 0.3, 7.0, 50.0, 400.0])
    both, both_v, capped = S.advance(x0, v0, spans, pot)
    assert not capped.any()
    for i, span in enumerate(spans):
        alone, alone_v, _ = S.advance(x0[i:i + 1], v0[i:i + 1], float(span), pot)
        assert alone[0] == pytest.approx(both[i], abs=1e-9) and alone_v[0] == pytest.approx(both_v[i], abs=1e-12)
    assert both[0] == pytest.approx(x0[0]) and both_v[0] == pytest.approx(v0[0])      # a span of 0 leaves the star where it is
    with pytest.raises(ValueError):
        S.advance(x0, v0, np.array([1.0, -1.0, 1.0, 1.0, 1.0]), pot)
