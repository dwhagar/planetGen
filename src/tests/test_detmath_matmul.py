# tests/test_detmath_matmul.py

"""
PERF.201: `detmath.matmul` adds a long inner index with `numpy.add.accumulate` instead of a loop, and the bits stay those of
the loop (a stored value must not depend on how the sum is made).
"""

import numpy as np
import pytest

from planetgen.util import detmath


def loop(a, b):
    """The definition: the products added left to right, one inner index at a time."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    total = a[..., 0:1] * b[0]
    for k in range(1, b.shape[0]):
        total = total + a[..., k:k + 1] * b[k]
    return total


@pytest.mark.parametrize("seed", range(40))
def test_the_sum_is_bit_for_bit_the_loops(seed):
    rng = np.random.default_rng(seed)
    inner = int(rng.integers(1, 150))
    columns = int(rng.integers(1, 20))
    leading = tuple(int(x) for x in rng.integers(1, 12, size=int(rng.integers(0, 3))))
    a = rng.standard_normal(leading + (inner,)) * 10.0 ** rng.integers(-30, 30, size=inner)
    b = rng.standard_normal((inner, columns)) * rng.choice([1e-9, 1.0, 1e9])
    a[..., int(rng.integers(0, inner))] = 0.0
    b[int(rng.integers(0, inner))] = -0.0
    got, want = detmath.matmul(a, b), loop(a, b)
    assert got.shape == want.shape and got.tobytes() == want.tobytes()


def test_a_sum_that_pairwise_addition_would_change_stays_left_to_right():
    a = np.array([[1.0] + [2.0 ** -53] * 40])
    b = np.ones((41, 1))
    assert detmath.matmul(a, b).tobytes() == loop(a, b).tobytes()
    assert detmath.matmul(a, b)[0, 0] == 1.0                    # each tiny term is lost to the 1.0 before it, as in the loop
    assert np.sum(a[0]) != 1.0                                  # numpy's own sum is pairwise and does not give it


def test_a_product_too_large_to_form_at_once_goes_through_the_loop(monkeypatch):
    monkeypatch.setattr(detmath, "_ACCUMULATE_LIMIT", 10)
    a, b = np.random.default_rng(3).standard_normal((4, 30)), np.random.default_rng(4).standard_normal((30, 5))
    assert detmath.matmul(a, b).tobytes() == loop(a, b).tobytes()


@pytest.mark.parametrize("limit", [200, 330, 1000])
@pytest.mark.parametrize("seed", range(12))
def test_a_product_added_an_output_at_a_time_is_bit_for_bit_the_loops(monkeypatch, seed, limit):
    """PERF.246: when the products do not fit at once and there are fewer outputs than terms, each output is added alone; the bits are the loop's."""
    monkeypatch.setattr(detmath, "_ACCUMULATE_LIMIT", limit)
    rng = np.random.default_rng(100 + seed)
    inner = int(rng.integers(40, 400))
    columns = int(rng.integers(1, 6))
    leading = tuple(int(x) for x in rng.integers(1, 5, size=int(rng.integers(0, 3))))
    a = rng.standard_normal(leading + (inner,)) * 10.0 ** rng.integers(-30, 30, size=inner)
    b = rng.standard_normal((inner, columns)) * rng.choice([1e-9, 1.0, 1e9])
    got, want = detmath.matmul(a, b), loop(a, b)
    assert got.shape == want.shape and got.tobytes() == want.tobytes()


def test_an_output_at_a_time_keeps_a_tiny_term_lost_to_the_sum_before_it(monkeypatch):
    monkeypatch.setattr(detmath, "_ACCUMULATE_LIMIT", 20)
    a = np.array([[1.0] + [2.0 ** -53] * 99])
    b = np.ones((100, 1))
    assert detmath.matmul(a, b)[0, 0] == 1.0


def test_a_long_gram_matrix_is_bit_for_bit_the_loops():
    """The seasonal swing's case: 6 series of 138,240 values, whose products do not fit at once."""
    a = np.random.default_rng(5).standard_normal((6, 138240)) * 10.0
    assert detmath.matmul(a, a.T).tobytes() == loop(a, a.T).tobytes()


def complex_loop(a, b):
    """The definition for complex numbers: each product from real ones (`ar*br - ai*bi`, `ar*bi + ai*br`), added left to right, one inner index at a time."""
    a = np.asarray(a, dtype=complex)
    b = np.asarray(b, dtype=complex)
    real, imaginary = None, None
    for k in range(b.shape[0]):
        re = a[..., k:k + 1].real * b[k].real - a[..., k:k + 1].imag * b[k].imag
        im = a[..., k:k + 1].real * b[k].imag + a[..., k:k + 1].imag * b[k].real
        real, imaginary = (re, im) if real is None else (real + re, imaginary + im)
    return _join(real, imaginary)


def _join(real, imaginary):
    out = np.empty(real.shape, dtype=complex)
    out.real, out.imag = real, imaginary
    return out


@pytest.mark.parametrize("seed", range(30))
def test_a_complex_sum_is_bit_for_bit_the_complex_loops(seed):
    """GEN.1787: the secular chain's modes are complex (`V @ (exp(i w t) * (V^-1 @ z))`); `matmul` and `dot` take them, with the real and imaginary parts formed as numpy forms them."""
    rng = np.random.default_rng(500 + seed)
    inner = int(rng.integers(1, 30))
    columns = int(rng.integers(1, 8))
    a = rng.standard_normal((int(rng.integers(1, 9)), inner)) + 1j * rng.standard_normal((1, inner)) * 10.0 ** rng.integers(-5, 5)
    b = rng.standard_normal((inner, columns)) * (1.0 + 2.0j)
    got, want = detmath.matmul(a, b), complex_loop(a, b)
    assert got.dtype == complex and got.tobytes() == want.tobytes()
    vector = b[:, 0]
    assert detmath.dot(a, vector).tobytes() == complex_loop(a, vector[:, None])[..., 0].tobytes()


def test_a_real_matrix_times_a_complex_vector_keeps_the_imaginary_part():
    a, z = np.array([[1.0, 2.0], [3.0, 4.0]]), np.array([1.0 + 1.0j, -2.0j])
    assert np.array_equal(detmath.dot(a, z), np.array([1.0 + 1.0j - 4.0j, 3.0 + 3.0j - 8.0j]))
    long_a, long_z = np.arange(1.0, 41.0).reshape(1, 40), (np.arange(40) + 1j)[:, None]
    assert detmath.matmul(long_a, long_z).tobytes() == complex_loop(long_a, long_z).tobytes()


def test_a_dot_of_two_vectors_is_a_python_number():
    assert detmath.dot(np.array([1.0, 2.0]), np.array([3.0, 4.0])) == 11.0 and isinstance(detmath.dot(np.array([1.0, 2.0]), np.array([3.0, 4.0])), float)
    assert detmath.dot(np.array([1j, 2.0]), np.array([1j, 1.0])) == -1.0 + 2.0 + 0.0j
