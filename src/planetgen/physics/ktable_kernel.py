# planetgen/physics/ktable_kernel.py

"""
The compiled cut of a k-distribution at the g edges (PERF.291)
==============================================================

`datasets.ir_ktables.reduce_distribution` sorts each row of opacities, accumulates the weight and the weighted opacity along the sorted row and reads both
at the g edges. The numpy version made the sorted copies, two cumulative sums and, for each edge, four gathers and a comparison over the whole array;
in an air history that is one fifth of the time of a wide binary's planets (0.012 s a call, 1,416 calls of 106 s, PERF.291).

The arithmetic is the numpy function's, in the same order: the order of each row is numpy's own `argsort` (passed in, since a compiled sort breaks ties in another
order and moves the last bits), the sums are sequential in float64 with no fast-math, so the answer is the same bit for bit
(`tests/test_ktable_kernel.py` holds it to `array_equal` on random rows with ties and on real rows). The pattern is `radiative_kernel`'s; like it, the module is
**kept out of the web app's import path** (numba's cache needs a writable code tree, `test_the_web_app_never_imports_numba`): `reduce_distribution` imports it when it is called.
"""

import numpy as np
from numba import njit

@njit(cache=True)
def cut(k, order, weights, edges, out):
    """
    Fills `out` (row, len(edges) - 1) with the mean of each g interval of each row of `k` (row, n), taken in `order` (row, n; the ascending order of each row),
    with `weights` (n) summing to 1: the cumulative integral of `k` over the weight is read at `edges`, linear between sorted points, and differenced.
    """
    rows, count = k.shape
    edge_count = edges.shape[0]
    cumulative_w = np.empty(count + 1)
    cumulative_kw = np.empty(count + 1)
    for row in range(rows):
        cumulative_w[0] = 0.0
        cumulative_kw[0] = 0.0
        total_w = 0.0
        total_kw = 0.0
        for position in range(count):
            source = order[row, position]
            weight = weights[source]
            total_w = total_w + weight
            total_kw = total_kw + k[row, source] * weight
            cumulative_w[position + 1] = total_w
            cumulative_kw[position + 1] = total_kw
        previous = 0.0
        for j in range(edge_count):
            edge = edges[j]
            threshold = edge - 1.0e-12                     # the numpy version's tolerance: a cumulative weight under this lies below the edge
            index = 0
            for position in range(count + 1):
                if cumulative_w[position] < threshold:
                    index += 1
            if index < 1:
                index = 1
            if index > count:
                index = count
            w0 = cumulative_w[index - 1]
            w1 = cumulative_w[index]
            k0 = cumulative_kw[index - 1]
            k1 = cumulative_kw[index]
            if w1 > w0:
                gap = w1 - w0
                if gap < 1.0e-300:                          # keeps the division finite where the interval has no weight
                    gap = 1.0e-300
                share = (edge - w0) / gap
            else:
                share = 0.0
            value = k0 + share * (k1 - k0)
            if j > 0:
                out[row, j - 1] = (value - previous) / (edge - edges[j - 1])
            previous = value
