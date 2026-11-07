# planetgen/galaxy/stellar_table.py

"""
The stars of the density model as a table on top of the analytic potential (GEN.1285; reports/gen1285-potential-fit-2026-11-05.md)
===================================================================================================================================

The potential (`galaxy.potential`) stands the density model's disks and spheres in as Miyamoto-Nagai and Hernquist profiles of the same mass (stand-ins GP1 and
GP2), and those hold 25 to 40 percent too little stellar density through the disk, a 40-fold cusp at the center and far too much beyond 16 kpc. This module
adds the difference: the potential of the **residual density** (the density model less the analytic profiles' own density, signed) so that analytic part and
table together are the density model, to the accuracy of the expansion. The analytic part keeps the large, steep part of the field in closed form; the table
carries only a smooth residual, which a short expansion holds well (the method of Dehnen and Binney's GalPot, applied to the residual rather than the whole).

- **Expansion.** The residual `d(r, cos(theta))`, averaged round the axis (the arms and the bar are averaged, as `enclosed_mass` averages them) and even in `z`,
  is cut into Legendre terms `d_l(r)` for even `l` up to `LMAX`, on a log radial grid and Gauss-Legendre nodes in `cos(theta)` packed towards the plane
  (`cos(theta) = t^2`). The potential of each term follows from two radial integrals, `A_l(r) = integral of d_l r' (r'/r)^(l+1) dr'` over `r' < r`
  and `B_l(r) = integral of d_l r' (r/r')^l dr'` over `r' > r`: `phi_l = -4 pi G (A_l + B_l) / (2 l + 1)` and the radial pull
  `4 pi G (l B_l - (l + 1) A_l) / ((2 l + 1) r)`. The integrals are exact for the density taken as piecewise linear in `ln r`, and every ratio is at most 1, so no
  power of the radius overflows.
- **Table.** The expansion is evaluated once on a grid of `ln r` and `t = sqrt(|cos theta|)` and kept: the residual's potential and its radial and polar
  pull (and, from the model itself, not the series, its density). A query is a bilinear interpolation, a few microseconds a thousand points, so the orbit stepper's million-star calls stay in closed-form time.
  Outside the grid the residual falls as a quadrupole (its mass is nil by construction) towards the outer edge and as `r` towards the center.
- **Cost.** About `RADIAL * ANGLES * AZIMUTHS` evaluations of the density model, once per galaxy and process (`potential.for_galaxy` keeps the result).
- **Limits.** A razor-thin disk is not held by any short expansion; the density model's thinnest disk is a tenth as thin as it is long, which is, to the
  accuracy the tests state.

Units: parsecs, solar masses, `(km/s)^2` for the potential, `(km/s)^2 / pc` for the pull.
"""

import hashlib
import math

import numpy as np
from scipy.interpolate import CubicSpline

from planetgen.util import detmath

G = 4.30091727e-3
"""float: Overridden by `galaxy.potential.G_PC_KMS2_PER_MSUN` at build time (see `build`)."""
LMAX = 24
"""int: The highest Legendre degree kept (even degrees only)."""
RADIAL = 280
"""int: Radial grid points of the expansion, log spaced from `INNER_PC` to `OUTER_PC`."""
ANGLES = 96
"""int: Gauss-Legendre nodes in `t` (`cos theta = t^2`) over the upper half sphere."""
AZIMUTHS = 12
"""int: Points round the axis the density is averaged over."""
INNER_PC = 0.25
OUTER_PC = 3.0e5
TABLE_U = 360
"""int: Points of the table along `ln r`."""
TABLE_T = 120
"""int: Points of the table along `t = sqrt(|cos theta|)`."""


def _legendre(x, lmax):
    """`(P, dP)`: Legendre polynomials and their derivatives of every degree 0 to `lmax` at the array `x`, shape `(lmax + 1,) + x.shape`."""
    p = np.empty((lmax + 1,) + x.shape)
    dp = np.empty_like(p)
    p[0], dp[0] = 1.0, 0.0
    if lmax >= 1:
        p[1], dp[1] = x, 1.0
    for n in range(2, lmax + 1):
        p[n] = ((2 * n - 1) * x * p[n - 1] - (n - 1) * p[n - 2]) / n
        dp[n] = dp[n - 2] + (2 * n - 1) * p[n - 1]
    return p, dp


def _segment_sums(q, du, degrees):
    """
    `(A, B)`, each `(RADIAL, len(degrees))`, from `q = d_l r^2` on the log grid (spacing `du`) by the exact integrals of the piecewise linear `q` against the
    exponentials (see the module doc), with the power-law tails beyond the grid's two ends.
    """
    n = q.shape[0]
    lam_a = degrees + 1.0
    lam_b = degrees.astype(float)
    ea = np.exp(-lam_a * du)
    eb = np.exp(-lam_b * du)
    wa0 = (1.0 - ea) / lam_a                                                    # integral of e^(-lam w) over the segment
    wa1 = (1.0 - ea * (1.0 + lam_a * du)) / (lam_a ** 2 * du)                  # integral of (w / h) e^(-lam w)
    safe_b = np.where(lam_b > 0.0, lam_b, 1.0)
    wb0 = np.where(lam_b > 0.0, (1.0 - eb) / safe_b, du)
    wb1 = np.where(lam_b > 0.0, (1.0 - eb * (1.0 + lam_b * du)) / (safe_b ** 2 * du), du / 2.0)
    a = np.empty_like(q)
    b = np.empty_like(q)
    a[0] = q[0] / (degrees + 3.0)
    for i in range(n - 1):
        a[i + 1] = ea * a[i] + q[i + 1] * wa0 + (q[i] - q[i + 1]) * wa1
    b[-1] = q[-1] / (degrees + 2.0)
    for i in range(n - 2, -1, -1):
        b[i] = eb * b[i + 1] + q[i] * wb0 + (q[i + 1] - q[i]) * wb1
    return a, b


class StellarTable:
    """The tabulated residual of the stars: build with `build`, query with `potential`, `pull` and `density` (arrays of galaxy-frame parsecs in, arrays out)."""

    def __init__(self, u, t, potential, radial_pull, polar_pull, density, digest):
        self._u, self._t = u, t
        self._fields = (potential, radial_pull, polar_pull, density)
        self._digest = digest

    def __repr__(self):
        return f"StellarTable({self._digest})"

    @property
    def digest(self):
        """str: 16 hex characters naming the table's numbers (the stamp of the gravity grid reads it)."""
        return self._digest

    def _locate(self, x, y, z):
        r = np.sqrt(x * x + y * y + z * z)
        safe = np.where(r > 0.0, r, 1.0)
        cos_t = np.abs(z) / safe
        u = np.log(np.clip(safe, INNER_PC, OUTER_PC))
        fu = (u - self._u[0]) / (self._u[1] - self._u[0])
        ft = np.sqrt(np.clip(cos_t, 0.0, 1.0)) / (self._t[1])
        iu = np.clip(fu.astype(np.int64), 0, len(self._u) - 2)
        it = np.clip(ft.astype(np.int64), 0, len(self._t) - 2)
        return r, safe, iu, fu - iu, it, ft - it

    def _interpolate(self, field, iu, wu, it, wt):
        f = field
        return ((1 - wu) * (1 - wt) * f[iu, it] + wu * (1 - wt) * f[iu + 1, it] + (1 - wu) * wt * f[iu, it + 1] + wu * wt * f[iu + 1, it + 1])

    def _scale(self, r, inner_power, outer_power):
        """The factor taking a value at the grid's edge to radius `r` beyond it: `(r / INNER)^inner_power` inside, `(r / OUTER)^-outer_power` outside."""
        return np.minimum(r / INNER_PC, 1.0) ** inner_power * (OUTER_PC / np.maximum(r, OUTER_PC)) ** outer_power

    def potential(self, x, y, z):
        """The residual's potential, `(km/s)^2`."""
        r, safe, iu, wu, it, wt = self._locate(x, y, z)
        return self._interpolate(self._fields[0], iu, wu, it, wt) * self._scale(r, 0.0, 3.0)

    def density(self, x, y, z):
        """The residual density, solar masses per pc^3 (the truncated series; the caller adds it to the analytic density)."""
        r, safe, iu, wu, it, wt = self._locate(x, y, z)
        return self._interpolate(self._fields[3], iu, wu, it, wt) * self._scale(r, 0.0, 4.0)

    def pull(self, x, y, z):
        """The residual's pull `(ax, ay, az)`, `(km/s)^2` per pc."""
        r, safe, iu, wu, it, wt = self._locate(x, y, z)
        a_r = self._interpolate(self._fields[1], iu, wu, it, wt) * self._scale(r, 1.0, 4.0)
        a_t = self._interpolate(self._fields[2], iu, wu, it, wt) * self._scale(r, 1.0, 4.0)
        cos_t = z / safe
        sin_t = np.sqrt(np.maximum(1.0 - cos_t * cos_t, 0.0))
        cyl = np.sqrt(x * x + y * y)
        cos_p = np.where(cyl > 0.0, x / np.where(cyl > 0.0, cyl, 1.0), 1.0)
        sin_p = np.where(cyl > 0.0, y / np.where(cyl > 0.0, cyl, 1.0), 0.0)
        # a_t is stored for z >= 0 along increasing theta; a mirror point at z < 0 has the mirrored polar component
        a_theta = np.where(z >= 0.0, a_t, -a_t)
        ax = a_r * sin_t * cos_p + a_theta * cos_t * cos_p
        ay = a_r * sin_t * sin_p + a_theta * cos_t * sin_p
        az = a_r * cos_t - a_theta * sin_t
        return ax, ay, az


def build(density_fn, analytic_density_fn, big_g):
    """
    The `StellarTable` of the residual `density_fn - analytic_density_fn` (each `f(x, y, z)` over arrays, solar masses per pc^3; the first is the density model,
    the second the analytic profiles' own density). `big_g` is G in pc (km/s)^2 per solar mass.
    """
    degrees = np.arange(0, LMAX + 1, 2)
    nodes, weights = np.polynomial.legendre.leggauss(ANGLES)
    t = (nodes + 1.0) / 2.0
    w = weights / 2.0 * 2.0 * t                    # integral over cos(theta) = t^2 from 0 to 1: f(t^2) 2 t dt
    x_nodes = t * t
    p_nodes, _ = _legendre(x_nodes, LMAX)
    p_nodes = p_nodes[degrees]
    r = np.geomspace(INNER_PC, OUTER_PC, RADIAL)
    du = math.log(OUTER_PC / INNER_PC) / (RADIAL - 1)
    phis = (np.arange(AZIMUTHS) + 0.5) * (2.0 * math.pi / AZIMUTHS)
    sin_t = np.sqrt(1.0 - x_nodes ** 2)
    rr = r[:, None, None]
    px = rr * sin_t[None, :, None] * np.cos(phis)[None, None, :]
    py = rr * sin_t[None, :, None] * np.sin(phis)[None, None, :]
    pz = rr * x_nodes[None, :, None] + 0.0 * px
    residual = (density_fn(px, py, pz) - analytic_density_fn(px, py, pz)).mean(axis=2)          # (RADIAL, ANGLES)
    residual = np.where(np.isfinite(residual), residual, 0.0)
    d_l = detmath.matmul(residual * w[None, :], p_nodes.T) * (2 * degrees + 1)                                # (RADIAL, n_degrees)
    a, b = _segment_sums(d_l * r[:, None] ** 2, du, degrees)
    coefficient = 4.0 * math.pi * big_g / (2 * degrees + 1)
    phi_l = -coefficient * (a + b)
    force_l = coefficient * (degrees * b - (degrees + 1) * a) / r[:, None]
    # tables on (ln r, t)
    u = np.linspace(math.log(INNER_PC), math.log(OUTER_PC), TABLE_U)
    tt = np.linspace(0.0, 1.0, TABLE_T)
    x_grid = tt * tt
    p_grid, dp_grid = _legendre(x_grid, LMAX)
    p_grid, dp_grid = p_grid[degrees], dp_grid[degrees]
    sin_grid = np.sqrt(np.maximum(1.0 - x_grid ** 2, 0.0))
    spline = CubicSpline(np.log(r), np.stack([phi_l, force_l, d_l], axis=0), axis=1)
    s_phi, s_force, _ = spline(u)                                                                 # each (TABLE_U, n_degrees)
    potential = detmath.matmul(s_phi, p_grid)
    radial = detmath.matmul(s_force, p_grid)
    # a_theta = (sin theta / r) sum phi_l P_l'(x)
    polar = (np.exp(-u)[:, None] * detmath.matmul(s_phi, dp_grid)) * sin_grid[None, :]
    # the density table is the residual itself on the grid, not the truncated series (the series cannot follow a thin disk's edge or the model's fall beyond 16 kpc)
    grid_r = np.exp(u)[:, None, None]
    grid_sin = sin_grid[None, :, None]
    gx = grid_r * grid_sin * np.cos(phis)[None, None, :]
    gy = grid_r * grid_sin * np.sin(phis)[None, None, :]
    gz = grid_r * x_grid[None, :, None] + 0.0 * gx
    density = (density_fn(gx, gy, gz) - analytic_density_fn(gx, gy, gz)).mean(axis=2)
    density = np.where(np.isfinite(density), density, 0.0)
    digest = hashlib.sha256(np.concatenate([potential.ravel()[::97], radial.ravel()[::97], polar.ravel()[::97]]).astype("<f8").tobytes()).hexdigest()[:16]
    return StellarTable(u, tt, potential, radial, polar, density, digest)
