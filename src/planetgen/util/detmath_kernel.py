# planetgen/util/detmath_kernel.py

"""
The long sums of `detmath.matmul`, compiled (PERF.279)
======================================================

`detmath.matmul` adds the products of a long inner index left to right, so the bits are the loop's and not BLAS's. For a Gram matrix of the seasonal
swing's series (6 series of up to 190,080 values) the numpy way forms every product in an array and adds with `numpy.add.accumulate`: each sum is one chain
of dependent additions, and the arrays are written and read again. This is the same sum, written out: every output starts from the product at index 0 and
adds the product at index 1, 2 and so on in that order, each product a single correctly rounded multiplication added by a single correctly rounded addition (no
fused multiply-add, no reordering: numba gives LLVM no licence for either). The chains of different outputs are interleaved in one pass over the data, which
a processor runs at the speed of its adders instead of their latency, so a sum comes out in a quarter of the time with the same bits
(`tests/test_detmath_matmul.py` compares them with the loop on random shapes, signed zeros and strided views).

Like the other kernels the module is imported inside the function that needs it (`detmath.matmul`, from `_kernel`), not at the top of `detmath`: the web
app imports `detmath` and must not import numba (`test_the_web_app_never_imports_numba`).
"""

import numpy as np
from numba import njit


@njit(cache=True)
def chains(a, b):
    """`a @ b` for a (rows, n) `a` and (n, m) `b` whose rows * m outputs are few: the index runs outside, so every output's chain advances at each step."""
    rows, n, columns = a.shape[0], a.shape[1], b.shape[1]
    out = np.empty((rows, columns))
    for r in range(rows):
        for j in range(columns):
            out[r, j] = a[r, 0] * b[0, j]
    for k in range(1, n):
        for r in range(rows):
            product = a[r, k]
            for j in range(columns):
                out[r, j] += product * b[k, j]
    return out


@njit(cache=True)
def wide(a, b):
    """`a @ b` for many outputs: each row of `a` adds its terms one index after another over the whole row of outputs, which is one vector operation each."""
    rows, n, columns = a.shape[0], a.shape[1], b.shape[1]
    out = np.empty((rows, columns))
    for r in range(rows):
        for j in range(columns):
            out[r, j] = a[r, 0] * b[0, j]
        for k in range(1, n):
            product = a[r, k]
            for j in range(columns):
                out[r, j] += product * b[k, j]
    return out
