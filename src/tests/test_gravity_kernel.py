# tests/test_gravity_kernel.py

"""
MAP.384: the gravity map's field kernel (`galaxy.gravity_kernel`) against a plain numpy sum, against the closed forms of one point mass, against finite
differences of itself, and for its determinism.
"""

import numpy as np
import pytest

from planetgen.galaxy import gravity_kernel as kernel
from planetgen.galaxy import potential

G = potential.G_PC_KMS2_PER_MSUN


def _reference(points, positions, masses, eps):
    """The same sum in numpy: points x sources x 3, the formulas of the module doc."""
    d = positions[None, :, :] - points[:, None, :]
    r2 = (d * d).sum(-1) + eps * eps
    w = G * masses[None, :]
    a = w / r2 ** 1.5
    b = 3.0 * w / r2 ** 2.5
    out = np.empty((len(points), 10))
    out[:, 0:3] = (a[..., None] * d).sum(1)
    out[:, 3] = -(w / np.sqrt(r2)).sum(1)
    for column, (i, j) in zip(range(4, 10), ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))):
        out[:, column] = (b * d[..., i] * d[..., j] - (a if i == j else 0.0)).sum(1)
    return out


def _random(sources=60, points=40, seed=3):
    rng = np.random.default_rng(seed)
    return rng.random((points, 3)) * 4.0 - 2.0, rng.random((sources, 3)) * 4.0 - 2.0, 0.08 + rng.random(sources) * 5.0


def test_the_kernel_matches_a_numpy_sum_of_the_same_terms():
    points, positions, masses = _random()
    got = kernel.fields(points, positions, masses, 0.01)
    np.testing.assert_allclose(got, _reference(points, positions, masses, 0.01), rtol=1e-9, atol=1e-12)


def test_one_point_mass_gives_the_closed_forms():
    out = kernel.fields(np.array([[3.0, 0.0, 0.0], [0.0, 2.0, 0.0]]), np.array([[0.0, 0.0, 0.0]]), np.array([2.0]), 0.0)
    # at (3,0,0): the pull points back at the source, -x
    assert out[0, kernel.PULL_X] == pytest.approx(-G * 2.0 / 9.0, rel=1e-12)
    assert out[0, kernel.POTENTIAL] == pytest.approx(-G * 2.0 / 3.0, rel=1e-12)
    assert out[0, kernel.T_XX] == pytest.approx(2.0 * G * 2.0 / 27.0, rel=1e-12)       # a stretch along the line to the source
    assert out[0, kernel.T_YY] == pytest.approx(-G * 2.0 / 27.0, rel=1e-12)
    assert out[0, kernel.T_ZZ] == pytest.approx(-G * 2.0 / 27.0, rel=1e-12)
    assert out[1, kernel.PULL_Y] == pytest.approx(-G * 2.0 / 4.0, rel=1e-12)
    assert kernel.tidal_strength(out)[0] == pytest.approx(2.0 * G * 2.0 / 27.0, rel=1e-12)


def test_the_pull_is_the_gradient_of_the_potential_and_the_tensor_that_of_the_pull():
    points, positions, masses = _random(sources=20, points=10)
    positions = positions + 6.0                      # keep every point well away from the sources
    h = 1.0e-4
    base = kernel.fields(points, positions, masses, 0.01)
    for axis in range(3):
        step = np.zeros(3)
        step[axis] = h
        up, down = kernel.fields(points + step, positions, masses, 0.01), kernel.fields(points - step, positions, masses, 0.01)
        np.testing.assert_allclose(-(up[:, kernel.POTENTIAL] - down[:, kernel.POTENTIAL]) / (2 * h), base[:, axis], rtol=1e-6, atol=1e-9)
        gradient = (up[:, 0:3] - down[:, 0:3]) / (2 * h)       # d a_i / d x_axis
        np.testing.assert_allclose(gradient, kernel.tensors(base)[:, :, axis], rtol=1e-5, atol=1e-8)


def test_the_tidal_tensor_has_no_trace_away_from_the_sources():
    points, positions, masses = _random()
    positions = positions + 8.0
    out = kernel.fields(points, positions, masses, 0.0)
    trace = out[:, kernel.T_XX] + out[:, kernel.T_YY] + out[:, kernel.T_ZZ]
    assert np.abs(trace).max() < 1e-12 * np.abs(out[:, kernel.T_XX]).max()


def test_a_point_at_a_source_stays_finite_with_the_softening():
    out = kernel.fields(np.zeros((1, 3)), np.zeros((1, 3)), np.array([1.0]), kernel.DEFAULT_SOFTENING_PC)
    assert np.isfinite(out).all()
    assert out[0, kernel.PULL_X] == 0.0 and out[0, kernel.POTENTIAL] == pytest.approx(-G / kernel.DEFAULT_SOFTENING_PC)


def test_the_same_inputs_give_the_same_bits_at_any_thread_count():
    points, positions, masses = _random(sources=300, points=500)
    one = kernel.fields(points, positions, masses, threads=1)
    several = kernel.fields(points, positions, masses, threads=4)
    assert one.tobytes() == several.tobytes()
    assert kernel.fields(points, positions, masses).tobytes() == kernel.fields(points, positions, masses).tobytes()


def test_another_order_of_the_sources_changes_only_the_last_bits():
    points, positions, masses = _random()
    order = np.random.default_rng(9).permutation(len(masses))
    np.testing.assert_allclose(kernel.fields(points, positions[order], masses[order]), kernel.fields(points, positions, masses), rtol=1e-12, atol=1e-14)


def test_no_sources_or_no_points_is_zeros_or_empty_and_bad_input_is_refused():
    points, positions, masses = _random(sources=3, points=4)
    assert not kernel.fields(points, positions[:0], masses[:0]).any() and kernel.fields(points, positions[:0], masses[:0]).shape == (4, 10)
    assert kernel.fields(points[:0], positions, masses).shape == (0, 10)
    with pytest.raises(ValueError):
        kernel.fields(points, positions, masses[:2])
    with pytest.raises(ValueError):
        kernel.fields(points, positions, -masses)
    with pytest.raises(ValueError):
        kernel.fields(points, positions, masses, softening_pc=-1.0)
    with pytest.raises(ValueError):
        kernel.fields(points[:, :2], positions, masses)


def test_the_softening_is_the_orbital_updates_member_softening():
    from planetgen.galaxy import acceleration
    assert kernel.DEFAULT_SOFTENING_PC == acceleration.MEMBER_SOFTENING_PC
    positions = np.array([[1.0, 0.5, -0.3], [-2.0, 1.0, 0.7]])
    masses = np.array([1.5, 0.4])
    point = np.array([[0.1, 0.2, 0.3]])
    d = positions - point
    expected = (G * masses[:, None] * d / ((d * d).sum(1)[:, None] + kernel.DEFAULT_SOFTENING_PC ** 2) ** 1.5).sum(0)
    np.testing.assert_allclose(kernel.fields(point, positions, masses)[0, 0:3], expected, rtol=1e-12)
