# planetgen/util/detmath.py

"""
Maths that gives the same bits on every Python and every machine (GEN.135).

Python's `math.hypot` and `math.dist` have changed their algorithm between releases (3.9 differs from 3.10 to 3.14 on
about a third of random inputs), and `sum()` of three or more floats changed in 3.12 (it sums with compensation now);
a value that goes into the database through either can differ from one interpreter to the next. These helpers use
only operations IEEE 754 fixes the result of (`+`, `-`, `*` and `sqrt`) in a fixed order, or `math.fsum`, which is
correctly rounded and the same in every version, so a galaxy does not depend on which Python made it.

    from planetgen.util import detmath

    r = detmath.hypot(x, y)            # sqrt(x*x + y*y), summed left to right
    d = detmath.dist(a, b)             # the same over the differences of two points
    s = detmath.fsum(masses)           # the exactly rounded sum

`matmul` replaces `@` and `numpy.dot`, whose results depend on the BLAS library installed.

`tests/test_reproducible_draws.py` fails on `math.hypot` or `math.dist` anywhere in the package outside `web/` (a
drawing's last bit does not matter) and the migrations (a revision that has run is never edited), and on `@` and
`numpy.dot`.
"""

import math

import numpy as np

fsum = math.fsum
"""The correctly rounded sum of an iterable of floats; the same result in every Python version."""

_sqrt = math.sqrt


def hypot(*components):
    """The square root of the squares of `components` added left to right: `hypot(x, y)` is `sqrt(x*x + y*y)`."""
    if len(components) == 2:
        x, y = components
        return _sqrt(x * x + y * y)
    if len(components) == 3:
        x, y, z = components
        return _sqrt(x * x + y * y + z * z)
    total = 0.0
    for component in components:
        total += component * component
    return _sqrt(total)


def dist(a, b):
    """The distance between points `a` and `b` (sequences of the same length): `hypot` of the differences."""
    count = len(a)
    if count != len(b):
        raise ValueError("both points must have the same number of coordinates")
    if count == 3:
        dx = a[0] - b[0]
        dy = a[1] - b[1]
        dz = a[2] - b[2]
        return _sqrt(dx * dx + dy * dy + dz * dz)
    if count == 2:
        dx = a[0] - b[0]
        dy = a[1] - b[1]
        return _sqrt(dx * dx + dy * dy)
    total = 0.0
    for first, second in zip(a, b):
        delta = first - second
        total += delta * delta
    return _sqrt(total)


def matmul(a, b):
    """
    `a @ b` for an array `a` of shape (..., n) and a matrix `b` of shape (n, m), summed over n left to right with
    element-wise operations only (no BLAS, whose blocking and fused multiply-adds differ between libraries).
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    total = a[..., 0:1] * b[0]
    for k in range(1, b.shape[0]):
        total = total + a[..., k:k + 1] * b[k]
    return total
