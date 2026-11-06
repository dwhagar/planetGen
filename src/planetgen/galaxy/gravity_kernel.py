# planetgen/galaxy/gravity_kernel.py

"""
The gravity map's field kernel (MAP.384; docs/design/orbits/gravity-map.md, PERF.72's report `perf72-gravity-grid-cost-2026-11-05.md`)
=======================================================================================================================================

One compiled function, `fields`, sums the pull, the potential and the tidal tensor of a set of point masses at a set of points. It is
the inner loop of the gravity map's sector grid: the caller (the interaction-list builder of MAP.385) hands it, for each target
block, the exact stars near the block and the point masses that stand for the far cells, and takes back the ten numbers per point.

- **Term.** Plummer-softened Newtonian, with `d` the source less the point and `r2 = |d|^2 + eps^2`:
  pull `a = G m d / r2^(3/2)`, potential `phi = -G m / r2^(1/2)`, tidal tensor `T_ij = d a_i / d x_j = G m (3 d_i d_j / r2^(5/2) - delta_ij / r2^(3/2))`.
  The tensor is the exact derivative of the softened pull, not a difference over a step, so the map's gradient mode needs no step and its trace is
  zero away from the sources (Laplace). A positive eigenvalue is a stretch, along the line to a source it is `2 G m / r^3`.
- **Output.** A float64 array `(points, 10)`, columns `PULL_X, PULL_Y, PULL_Z, POTENTIAL, T_XX, T_YY, T_ZZ, T_XY, T_XZ, T_YZ`.
- **Determinism.** Each point sums its sources in the order given, in float64, with no fast-math, and the points are shared between the threads,
  not the sources, so the same inputs give the same bits at any thread count. A different source order differs in the last bits.
- **Cost** (PERF.72, 4 cores at 2.1 GHz): about 0.55 billion point-source pairs a second, so 4096 zones against 1,000 sources take about 8 ms;
  the numpy sum of the same pairs is 54 to 65 times slower and is not offered.

Units: positions in parsecs, masses in solar masses, the pull in (km/s)^2 per pc, the potential in (km/s)^2 and the tensor in (km/s)^2 per pc^2
(`galaxy.potential.G_PC_KMS2_PER_MSUN`). The module imports numba, so the web app must not import it (OPS.42, `test_the_web_app_never_imports_numba`);
the sector grid is a job that stores its result.
"""

import math

import numpy as np
from numba import njit, prange, get_num_threads, set_num_threads

from planetgen.galaxy import potential as galactic

PULL_X, PULL_Y, PULL_Z, POTENTIAL, T_XX, T_YY, T_ZZ, T_XY, T_XZ, T_YZ = range(10)
"""int: The columns of `fields`' result."""
COLUMNS = 10
"""int: The width of `fields`' result."""
DEFAULT_SOFTENING_PC = 1.0e-3
"""float: The Plummer softening of a source, parsecs; the same 0.001 pc as `galaxy.acceleration.MEMBER_SOFTENING_PC`."""


@njit(parallel=True, cache=True)
def _sum(points, positions, masses, softening2, big_g, out):
    """Fills `out` (points, 10) with the sums over every source, point by point (see `fields`)."""
    for i in prange(points.shape[0]):
        gx = gy = gz = phi = txx = tyy = tzz = txy = txz = tyz = 0.0
        px, py, pz = points[i, 0], points[i, 1], points[i, 2]
        for j in range(positions.shape[0]):
            dx = positions[j, 0] - px
            dy = positions[j, 1] - py
            dz = positions[j, 2] - pz
            r2 = dx * dx + dy * dy + dz * dz + softening2
            inverse = 1.0 / math.sqrt(r2)
            inverse3 = inverse / r2
            inverse5 = inverse3 / r2
            weight = big_g * masses[j]
            pull = weight * inverse3
            spread = 3.0 * weight * inverse5
            gx += pull * dx
            gy += pull * dy
            gz += pull * dz
            phi -= weight * inverse
            txx += spread * dx * dx - pull
            tyy += spread * dy * dy - pull
            tzz += spread * dz * dz - pull
            txy += spread * dx * dy
            txz += spread * dx * dz
            tyz += spread * dy * dz
        out[i, 0] = gx
        out[i, 1] = gy
        out[i, 2] = gz
        out[i, 3] = phi
        out[i, 4] = txx
        out[i, 5] = tyy
        out[i, 6] = tzz
        out[i, 7] = txy
        out[i, 8] = txz
        out[i, 9] = tyz


def fields(points, positions, masses_solar, softening_pc=DEFAULT_SOFTENING_PC, threads=None):
    """
    The pull, potential and tidal tensor of the point masses `positions` (sources, parsecs, shape `(n, 3)`) of `masses_solar` (shape `(n,)`) at `points`
    (parsecs, shape `(p, 3)`): a float64 array `(p, 10)` with the columns `PULL_X` to `T_YZ`. No sources, or no points, gives zeros, or an empty array.

    Args:
        softening_pc (float): The Plummer softening, parsecs, zero or more (zero is legal when no point is at a source).
        threads (int, optional): The number of threads to use for this call; `None` leaves numba's setting as it is.

    Raises:
        ValueError: If a shape is wrong, or a mass or the softening is negative or not finite.
    """
    points = np.ascontiguousarray(points, dtype=np.float64)
    positions = np.ascontiguousarray(positions, dtype=np.float64)
    masses = np.ascontiguousarray(masses_solar, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or positions.ndim != 2 or positions.shape[1] != 3 or masses.shape != (positions.shape[0],):
        raise ValueError("points and positions are (n, 3) arrays and masses_solar has one mass per position")
    if not (math.isfinite(softening_pc) and softening_pc >= 0.0):
        raise ValueError(f"the softening must be zero or more and finite, not {softening_pc!r}")
    if masses.size and not (np.all(np.isfinite(masses)) and np.all(masses >= 0.0)):
        raise ValueError("every mass must be zero or more and finite")
    out = np.zeros((points.shape[0], COLUMNS))
    if points.shape[0] == 0 or positions.shape[0] == 0:
        return out
    previous = get_num_threads()
    if threads is not None:
        set_num_threads(max(1, int(threads)))
    try:
        _sum(points, positions, masses, softening_pc * softening_pc, galactic.G_PC_KMS2_PER_MSUN, out)
    finally:
        set_num_threads(previous)
    return out
