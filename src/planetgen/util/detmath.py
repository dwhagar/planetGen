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

`matmul` replaces `@` and `numpy.dot`, whose results depend on the BLAS library installed; `dot` is its form for a vector.

`tests/test_reproducible_draws.py` fails on `math.hypot` or `math.dist` anywhere in the package outside `web/` (a
drawing's last bit does not matter) and the migrations (a revision that has run is never edited), and on `@` and
`numpy.dot`.
"""

import math

import numpy as np

from planetgen.util import log

fsum = math.fsum
"""The correctly rounded sum of an iterable of floats; the same result in every Python version."""

_sqrt = math.sqrt


def float_or(value, empty=math.nan):
    """`value` as a float, or `empty` for `None` (a NULL column read into an array)."""
    return empty if value is None else float(value)


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


_ACCUMULATE_FROM = 6
"""int: The inner length from which `matmul` forms the products at once and adds them with `numpy.add.accumulate`; below it
the loop over the inner index is cheaper than the larger array."""

_ACCUMULATE_LIMIT = 4_000_000
"""int: The most products `matmul` forms at once; a larger product is added one output at a time (PERF.246) or goes through the loop."""

_KERNEL_FROM = 100_000
"""int: The number of products (rows of `a` times the inner index times the columns of `b`) from which `matmul` adds with the compiled kernel (`detmath_kernel`, PERF.279); fewer go through numpy."""

_KERNEL_CHAINS = 64
"""int: Up to this many outputs the kernel advances every output's chain at each step of the index (a Gram row, a few dot products); more, it adds a whole row of outputs at each step (a few terms over a long vector)."""

_kernel_module = []
"""list: The kernel module once imported (one entry), or `[None]` where numba cannot be imported and `matmul` keeps the numpy way."""


def _kernel():
    """The compiled kernels, imported on first use; `None` if numba is not there or cannot cache in this code tree (the numpy way gives the same bits)."""
    if not _kernel_module:
        try:
            from planetgen.util import detmath_kernel
            _kernel_module.append(detmath_kernel)
        except Exception as exc:  # noqa: BLE001 -- numba missing, or no writable cache folder: the slower sum has the same bits
            log.swallowed("util/detmath", "importing the compiled kernels", exc)
            _kernel_module.append(None)
    return _kernel_module[0]


def dot(a, b):
    """
    `a @ b` for a vector `b` of length n: two vectors give a float, an array `a` of shape (..., n) gives the array of its dot products with `b`. Summed over
    n left to right by `matmul`, so it has the same bits.
    """
    out = matmul(a, np.asarray(b)[:, None])[..., 0]
    return out.item() if out.ndim == 0 else out


def _complex_matmul(a, b):
    """
    `matmul` for complex arrays. numpy's own complex product is not fixed (its vector loops fuse the multiplications and additions on some shapes and processors), so each product is formed from real
    ones, as the definition: the real part `ar*br - ai*bi` and the imaginary part `ar*bi + ai*br`, each a plain IEEE operation, and the products are added left to right (the real and imaginary parts separately).
    """
    ar, ai, br, bi = a.real, a.imag, b.real, b.imag
    real = ar[..., 0:1] * br[0] - ai[..., 0:1] * bi[0]
    imaginary = ar[..., 0:1] * bi[0] + ai[..., 0:1] * br[0]
    for k in range(1, b.shape[0]):
        real = real + (ar[..., k:k + 1] * br[k] - ai[..., k:k + 1] * bi[k])
        imaginary = imaginary + (ar[..., k:k + 1] * bi[k] + ai[..., k:k + 1] * br[k])
    out = np.empty(real.shape, dtype=complex)
    out.real, out.imag = real, imaginary
    return out


def matmul(a, b):
    """
    `a @ b` for an array `a` of shape (..., n) and a matrix `b` of shape (n, m), summed over n left to right with
    element-wise operations only (no BLAS, whose blocking and fused multiply-adds differ between libraries).

    A long inner index forms every product at once (each one a single correctly rounded multiplication, as in the loop) and
    adds them with `numpy.add.accumulate`, which adds in order and never pairwise (`numpy.add.reduce` does), so the bits
    are those of the loop (`tests/test_detmath.py` compares them on random shapes).
    """
    a, b = np.asarray(a), np.asarray(b)
    if np.iscomplexobj(a) or np.iscomplexobj(b):
        return _complex_matmul(a.astype(complex, copy=False), b.astype(complex, copy=False))
    a, b = a.astype(float, copy=False), b.astype(float, copy=False)
    inner = b.shape[0]
    if inner and a.size * b.shape[1] >= _KERNEL_FROM:
        kernel = _kernel()
        if kernel is not None:
            # PERF.279: the same sums, every output's products added left to right, in one compiled pass (no array of products, the chains of the outputs advance together)
            flat = np.ascontiguousarray(a.reshape(-1, inner))
            outputs = flat.shape[0] * b.shape[1]
            if outputs > _KERNEL_CHAINS:
                out = kernel.wide(flat, b)
            elif flat.shape[0] > b.shape[1]:
                out = kernel.chains(b.T, flat.T).T          # a product commutes exactly, so the output is the same and its few chains run along the longer side
            else:
                out = kernel.chains(flat, b)
            return out.reshape(a.shape[:-1] + (b.shape[1],))
    if inner >= _ACCUMULATE_FROM and a.size * b.shape[1] <= _ACCUMULATE_LIMIT:
        return np.add.accumulate(a[..., :, None] * b, axis=-2)[..., -1, :]
    rows = a.size // inner if inner else 0
    if inner >= _ACCUMULATE_FROM and rows * b.shape[1] <= inner:
        # PERF.246: a product too large to form at once, with fewer outputs than inner terms (a Gram matrix of 6 series of 190,080 values in the seasonal swing), adds each
        # output's terms alone: its products in one buffer and `numpy.add.accumulate` over it, left to right as the loop does, so the bits are the loop's; the loop
        # below is `inner` array operations and took 0.26 s a call, this is rows * columns of them on long contiguous vectors and takes 0.05 s.
        flat, columns = a.reshape(rows, inner), np.ascontiguousarray(b.T)
        out = np.empty((rows, b.shape[1]))
        buffer = np.empty(inner)
        for row in range(rows):
            for column in range(b.shape[1]):
                np.multiply(flat[row], columns[column], out=buffer)
                np.add.accumulate(buffer, out=buffer)
                out[row, column] = buffer[-1]
        return out.reshape(a.shape[:-1] + (b.shape[1],))
    total = a[..., 0:1] * b[0]
    for k in range(1, inner):
        total = total + a[..., k:k + 1] * b[k]
    return total
