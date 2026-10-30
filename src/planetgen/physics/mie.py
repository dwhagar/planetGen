# planetgen/physics/mie.py

"""
Mie scattering of a homogeneous sphere (GEN.1020; Bohren and Huffman 1983, chapter 4), vectorised over spheres and size distributions
=======================================================================================================================================

`efficiencies(x, m)` is the extinction and scattering efficiencies and the asymmetry factor of spheres of size parameter `x = 2 pi r / lambda`
and complex refractive index `m = n + i k` (the sign of the imaginary part is the absorbing one). `gamma_average` averages them over a
gamma size distribution of effective radius and effective variance (Hansen and Travis 1974), weighted by the cross section.

The coefficients `a_n` and `b_n` come from the usual recurrences: the logarithmic derivative `D_n(m x)` downward, the Riccati-Bessel
functions of `x` upward, each sphere summed to its own `n_max = x + 4 x^(1/3) + 2`. Units: lengths in the same unit for `r` and `lambda`.
"""

import numpy as np


@np.errstate(all="ignore")
def efficiencies(x, m):
    """
    `(q_ext, q_sca, g)` of spheres: `x` and `m` broadcast to one shape; `g` is the mean cosine of the scattering angle (0 where nothing scatters).
    """
    x, m = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(m, dtype=complex))
    shape = x.shape
    x, m = x.ravel(), m.ravel()
    count = len(x)
    n_max = (x + 4.0 * x ** (1.0 / 3.0) + 2.0).astype(int)
    top = int(n_max.max())
    y = m * x
    start = int(max(top, np.abs(y).max())) + 16
    d = np.zeros((top + 1, count), dtype=complex)
    d_run = np.zeros(count, dtype=complex)
    for n in range(start, 0, -1):
        d_run = n / y - 1.0 / (d_run + n / y)
        if n - 1 <= top:
            d[n - 1] = d_run
    psi_prev, psi = np.cos(x), np.sin(x)
    xi_prev, xi = psi_prev - 1j * (-np.sin(x)), psi - 1j * np.cos(x)     # xi_n = psi_n - i chi_n, chi_0 = cos x, chi_-1 = -sin x
    q_ext = np.zeros(count)
    q_sca = np.zeros(count)
    asym = np.zeros(count)
    a_old = b_old = None
    for n in range(1, top + 1):
        factor = (2.0 * n - 1.0) / x
        psi_new = factor * psi - psi_prev
        xi_new = factor * xi - xi_prev
        psi_prev, psi = psi, psi_new
        xi_prev, xi = xi, xi_new
        a = ((d[n] / m + n / x) * psi - psi_prev) / ((d[n] / m + n / x) * xi - xi_prev)
        b = ((d[n] * m + n / x) * psi - psi_prev) / ((d[n] * m + n / x) * xi - xi_prev)
        live = n <= n_max
        a = np.where(live, a, 0.0)
        b = np.where(live, b, 0.0)
        q_ext += (2.0 * n + 1.0) * np.real(a + b)
        q_sca += (2.0 * n + 1.0) * (np.abs(a) ** 2 + np.abs(b) ** 2)
        asym += (2.0 * n + 1.0) / (n * (n + 1.0)) * np.real(a * np.conj(b))
        if a_old is not None:
            asym += (n - 1.0) * (n + 1.0) / n * np.real(a_old * np.conj(a) + b_old * np.conj(b))
        a_old, b_old = a, b
    q_ext = 2.0 / x ** 2 * q_ext
    q_sca = 2.0 / x ** 2 * q_sca
    g = np.where(q_sca > 0.0, 4.0 / (x ** 2 * np.maximum(q_sca, 1e-300)) * asym, 0.0)
    return q_ext.reshape(shape), q_sca.reshape(shape), g.reshape(shape)


def gamma_average(wavelength, m, effective_radius, effective_variance=0.2, points=48):
    """
    `(q_ext, q_sca, g)` averaged over a gamma size distribution `n(r) ~ r^((1 - 3 v) / v) exp(-r / (a v))` of effective radius `a` and effective
    variance `v`, weighted by the geometric cross section (so `q` is the efficiency of the cross-section-weighted sphere).
    `wavelength` and `effective_radius` share a length unit; `m` broadcasts with `wavelength`.
    """
    a, v = float(effective_radius), float(effective_variance)
    r = a * np.exp(np.linspace(np.log(0.03), np.log(12.0), points))
    k = (1.0 - 3.0 * v) / v
    density = r ** k * np.exp(-r / (a * v))
    weight = density * np.pi * r ** 2 * r          # per ln r
    weight = weight / weight.sum()
    wavelength = np.atleast_1d(np.asarray(wavelength, dtype=float))
    index = np.broadcast_to(np.asarray(m, dtype=complex), wavelength.shape)
    x = 2.0 * np.pi * r[None, :] / wavelength[:, None]
    q_ext, q_sca, g = efficiencies(x, index[:, None])
    # the cross-section weighting uses q_sca for the asymmetry mean
    mean_sca = (weight * q_sca).sum(axis=1)
    mean_g = np.where(mean_sca > 0.0, (weight * q_sca * g).sum(axis=1) / np.maximum(mean_sca, 1e-300), 0.0)
    return (weight * q_ext).sum(axis=1), mean_sca, mean_g
