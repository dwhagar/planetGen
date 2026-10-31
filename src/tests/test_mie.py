# tests/test_mie.py

"""GEN.1020: Mie spheres (`planetgen.physics.mie`) against their closed-form limits."""

import numpy as np
import pytest

from planetgen.physics import mie


def test_a_large_clear_sphere_extinguishes_twice_its_cross_section():
    q_ext, q_sca, _g = mie.efficiencies(np.array([100.0, 200.0]), 1.33 + 0j)
    assert q_ext == pytest.approx(q_sca) and np.all(np.abs(q_ext - 2.0) < 0.15)


def test_a_small_sphere_scatters_as_rayleigh():
    x, m = 0.01, 1.5 + 0j
    assert mie.efficiencies(x, m)[1] == pytest.approx(8.0 / 3.0 * x ** 4 * abs((m ** 2 - 1) / (m ** 2 + 2)) ** 2, rel=1e-3)


def test_a_small_absorbing_sphere_absorbs_by_the_dipole_formula():
    x, m = 0.01, 1.5 + 0.1j
    q_ext, q_sca, _g = mie.efficiencies(x, m)
    assert q_ext - q_sca == pytest.approx(4.0 * x * np.imag((m ** 2 - 1) / (m ** 2 + 2)), rel=1e-3)


def test_the_bohren_and_huffman_example():
    """x = 5.213, m = 1.55: the extinction of a clear sphere there is 3.105 (their table)."""
    assert mie.efficiencies(5.213, 1.55 + 0j)[0] == pytest.approx(3.105, abs=0.002)


def test_the_asymmetry_grows_with_the_size_parameter():
    g = mie.efficiencies(np.array([0.1, 1.0, 10.0, 50.0]), 1.43 + 0j)[2]
    assert np.all(np.diff(g) > 0.0) and g[0] < 0.05 and g[-1] > 0.8


def test_an_absorbing_sphere_has_less_scattering_than_extinction():
    q_ext, q_sca, _g = mie.efficiencies(20.0, 1.4 + 0.2j)
    assert q_ext > q_sca > 0.0


def test_the_gamma_average_of_big_droplets_is_two():
    q_ext, q_sca, g = mie.gamma_average(np.array([0.55]), np.array([1.43 + 0j]), 5.0, 0.2)
    assert q_ext[0] == pytest.approx(2.0, abs=0.25) and q_sca[0] == pytest.approx(q_ext[0]) and 0.75 < g[0] < 0.9
