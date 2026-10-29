# tests/test_planck_fraction.py

"""
PERF.201: `infrared_opacity.planck_fraction_below` forms the terms of its series at once and adds only as many as each `x` needs, and
the bits stay those of the loop that adds all 40 (a stored optical depth must not depend on how the sum is made).
"""

import math
import warnings

import numpy as np
import pytest

from planetgen.physics import infrared_opacity as io


def loop(x):
    """The definition: 40 terms added to a running total, one `n` at a time."""
    x = np.asarray(x, dtype=float)
    small = np.minimum(x, 2.0)
    series = 15.0 / math.pi ** 4 * (small ** 3 / 3.0 - small ** 4 / 8.0 + small ** 5 / 60.0 - small ** 7 / 5040.0
                                    + small ** 9 / 272160.0 - small ** 11 / 13305600.0)
    large = np.maximum(x, 2.0)
    total = np.zeros_like(large)
    for n in range(1, 41):
        total += np.exp(-n * large) / n * (large ** 3 + 3.0 * large ** 2 / n + 6.0 * large / n ** 2 + 6.0 / n ** 3)
    return np.where(x < 2.0, series, 1.0 - 15.0 / math.pi ** 4 * total)


def same_bits(x):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        want, got = loop(x), io.planck_fraction_below(x)
    return want.shape == got.shape and want.tobytes() == got.tobytes()


@pytest.mark.parametrize("temperature_k", [3.0, 20.0, 150.0, 255.0, 288.0, 400.0, 1200.0, 6000.0])
def test_the_grid_of_a_planet_has_the_bits_of_the_loop(temperature_k):
    assert same_bits(io.EDGES_CM1 * io.SECOND_RADIATION_CM_K / temperature_k)


@pytest.mark.parametrize("seed", range(25))
def test_random_arguments_have_the_bits_of_the_loop(seed):
    rng = np.random.default_rng(seed)
    x = rng.uniform(0.0, rng.choice([3.0, 10.0, 50.0, 500.0]), size=int(rng.integers(1, 700)))
    assert same_bits(x)


def test_the_edges_of_the_ranges_have_the_bits_of_the_loop():
    edges = [1.9999999, 2.0, 2.0000001, 5.999999, 6.0, 19.999999, 20.0, 49.999999, 50.0, 700.0, 1e300, np.inf, np.nan, 0.0, 1e-300]
    assert same_bits(np.array(edges))
    assert same_bits(np.array(3.0)) and same_bits(np.array([[1.0, 2.5], [7.0, 30.0]]))
