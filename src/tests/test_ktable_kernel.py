"""
The compiled cut of a k-distribution (PERF.291, `planetgen.physics.ktable_kernel`): the same bits as the numpy arithmetic it replaced, on rows with ties, on float32 rows,
on any leading shape, and on the weights and edges the infrared tables use.
"""

import numpy as np
import pytest

from planetgen.datasets import ir_ktables
from planetgen.physics import line_opacity


def _numpy_reduce(k, weights, edges):
    """The function as it was before the kernel (the numpy arithmetic, kept here as the reference)."""
    order = np.argsort(k, axis=-1)
    k_sorted = np.take_along_axis(k, order, axis=-1)
    w_sorted = weights[order]
    cumulative_w = np.concatenate([np.zeros(w_sorted.shape[:-1] + (1,)), np.cumsum(w_sorted, axis=-1)], axis=-1)
    cumulative_kw = np.concatenate([np.zeros(w_sorted.shape[:-1] + (1,)), np.cumsum(k_sorted * w_sorted, axis=-1)], axis=-1)
    size = cumulative_w.shape[-1]
    at_edges = np.empty(k.shape[:-1] + (len(edges),))
    for j, edge in enumerate(edges):
        index = np.clip((cumulative_w < edge - 1.0e-12).sum(axis=-1), 1, size - 1)
        w0 = np.take_along_axis(cumulative_w, (index - 1)[..., None], -1)[..., 0]
        w1 = np.take_along_axis(cumulative_w, index[..., None], -1)[..., 0]
        k0 = np.take_along_axis(cumulative_kw, (index - 1)[..., None], -1)[..., 0]
        k1 = np.take_along_axis(cumulative_kw, index[..., None], -1)[..., 0]
        share = np.where(w1 > w0, (edge - w0) / np.maximum(w1 - w0, 1.0e-300), 0.0)
        at_edges[..., j] = k0 + share * (k1 - k0)
    return np.diff(at_edges, axis=-1) / np.diff(edges)


G_EDGES = np.array(line_opacity.G_EDGES, dtype=float)


def _weights(count, rng):
    weights = rng.uniform(0.2, 1.0, count)
    return weights / weights.sum()


@pytest.mark.parametrize("shape", [(100,), (7, 100), (3, 5, 25), (40, 130, 25), (1, 10)])
def test_the_cut_equals_the_numpy_arithmetic_bit_for_bit(shape):
    rng = np.random.default_rng(shape[-1])
    k = np.exp(rng.normal(0.0, 4.0, shape))
    weights = _weights(shape[-1], rng)
    assert np.array_equal(ir_ktables.reduce_distribution(k, weights, G_EDGES), _numpy_reduce(k, weights, G_EDGES))


def test_rows_with_ties_and_zeros_come_out_the_same():
    rng = np.random.default_rng(5)
    k = np.round(np.exp(rng.normal(0.0, 2.0, (300, 100))), 1)       # many equal values
    k[::3, :40] = 0.0                                               # and rows with a run of zeros
    weights = _weights(100, rng)
    assert np.array_equal(ir_ktables.reduce_distribution(k, weights, G_EDGES), _numpy_reduce(k, weights, G_EDGES))


def test_float32_rows_come_out_the_same():
    rng = np.random.default_rng(6)
    k = np.exp(rng.normal(0.0, 3.0, (50, 25))).astype(np.float32)
    weights = _weights(25, rng)
    assert np.array_equal(ir_ktables.reduce_distribution(k, weights, G_EDGES), _numpy_reduce(k, weights, G_EDGES))


def test_an_overlap_of_two_distributions_goes_through_the_cut():
    """`line_opacity.resort`: the 100 sums of two sorted distributions, cut at the g edges."""
    rng = np.random.default_rng(7)
    first = np.sort(np.exp(rng.normal(0.0, 3.0, (60, 10))), axis=-1)
    second = np.sort(np.exp(rng.normal(0.0, 3.0, (60, 10))), axis=-1)
    combined = (first[..., :, None] + second[..., None, :]).reshape(60, 100)
    weights = (line_opacity.G_WEIGHTS[:, None] * line_opacity.G_WEIGHTS[None, :]).reshape(-1)
    assert np.array_equal(line_opacity.resort(first, second), _numpy_reduce(combined, weights, line_opacity.G_EDGES))


def test_the_kernel_is_not_on_the_web_apps_import_path():
    import subprocess
    import sys
    code = "import sys\nfrom planetgen.datasets import ir_ktables\nassert 'planetgen.physics.ktable_kernel' not in sys.modules\nassert 'numba' not in sys.modules\n"
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-800:]
