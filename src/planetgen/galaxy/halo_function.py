# planetgen/galaxy/halo_function.py

"""
The halo mass function and the halo mass a galaxy of a given stellar mass has (GEN.1154; GEN.399 report
`reports/gen399-halo-mass-function-2026-10-31.md`).

`mass_chain` draws the stellar mass first, then the halo. A stellar-to-halo relation gives the stellar mass at a halo
mass with a scatter in `log M*` at fixed `Mh`, so the halos that hold a galaxy of a given `M*` are drawn from

    P(Mh | M*)  proportional to  n(Mh) x LogNormal(log M*; log m*(Mh), scatter),

where `n(Mh)` is the halo mass function. Halos fall steeply in number with mass, so the heavy side is thinner than the
light side. Inverting the median relation and adding a Gaussian in `log Mh` (the first build, GEN.368) is exact at
the light end and puts the halo up to 0.34 dex too heavy at 2e11 solar masses of stars.

* The halo mass function is Tinker et al. 2008 (ApJ 688, 709, arXiv:0803.2706) Eq. 3, `dn/dlnM = f(sigma) rho_m / M
  |dln sigma^-1 / dlnM|`, `f = A [(sigma / b)^-a + 1] exp(-c / sigma^2)`, its Table 2 interpolated in `log Delta`
  (Appendix B practice) at the overdensity the chain's halos are measured at (`mass_chain.OVERDENSITY` times the
  critical density, `Delta = OVERDENSITY / Omega_m` times the mean).
* `sigma(M)` is the top-hat variance of the linear power spectrum, normalised to `SIGMA_8`, with the Eisenstein and Hu
  1998 (ApJ 496, 605, arXiv:astro-ph/9709112) no-wiggle transfer function (Eqs. 26 and 28 to 31).
* The cosmology is Planck 2013 (Planck Collaboration XVI, arXiv:1303.5076, Table 2, Planck best fit: Omega_b h^2
  0.022068, n_s 0.9624, Omega_m 0.3175, sigma_8 0.8344, h 0.6711); `mass_chain` already carried Omega_m 0.3175 and h 0.671.

Nothing here reads the Milky Way or a dataset. The function is built once per process (about a second) on a grid
of `log Mh` and then every conditional is a cumulative sum over the same grid.
"""

import bisect
import functools
import math

from planetgen.physics import constants
from planetgen.physics import unit_registry as units

HUBBLE_H = 0.671
"""float: The Hubble constant in units of 100 km/s/Mpc."""
OMEGA_MATTER = 0.3175
"""float: The matter density parameter."""
OMEGA_BARYON_H2 = 0.022068
"""float: The baryon density parameter times h squared (Planck 2013 XVI Table 2). [Stand-in HF1, tier T1, tweak L0.]"""
SPECTRAL_INDEX = 0.9624
"""float: The primordial spectral index n_s (Planck 2013 XVI Table 2). [Stand-in HF1, tier T1, tweak L0.]"""
SIGMA_8 = 0.8344
"""float: The rms linear fluctuation in 8 Mpc/h spheres (Planck 2013 XVI Table 2). [Stand-in HF1, tier T1, tweak L0.]"""
CMB_TEMPERATURE_K = 2.728
"""float: The CMB temperature the transfer function's Theta is built on (Eisenstein and Hu 1998 use 2.728 K)."""

# Tinker et al. 2008 Table 2: Delta (times the mean density), then A, a, b, c.
TINKER_TABLE = (
    (200.0, 0.186, 1.47, 2.57, 1.19), (300.0, 0.200, 1.52, 2.25, 1.27), (400.0, 0.212, 1.56, 2.05, 1.34),
    (600.0, 0.218, 1.61, 1.87, 1.45), (800.0, 0.248, 1.87, 1.59, 1.58), (1200.0, 0.255, 2.13, 1.51, 1.80),
    (1600.0, 0.260, 2.30, 1.46, 1.97), (2400.0, 0.260, 2.53, 1.44, 2.24), (3200.0, 0.260, 2.66, 1.41, 2.44),
)
"""tuple: Rows of (Delta, A, a, b, c) of the Tinker et al. 2008 fit to the halo mass function, z = 0."""

LOG_MASS_LOW = 8.0
LOG_MASS_HIGH = 17.0
"""float: The halo masses the conditional is built over, log10 solar masses."""
GRID_STEP = 0.01
"""float: The step of the halo mass grid the conditional is summed on, dex."""
_COARSE_STEP = 0.05
"""float: The step the function itself is computed on; between its points `ln n` is interpolated linearly. [Stand-in
HF4, tier T9, tweak L0, with `GRID_STEP` and the wave number grid.]"""
_LN_K_LOW = math.log(1.0e-4)
_LN_K_HIGH = math.log(2.0e3)
_LN_K_STEP = 0.02

_msun_to_kg = units.converter("Msun", "kg")
_kpc_to_m = units.converter("kpc", "m")


def _mean_matter_density_msun_mpc3():
    """The mean matter density today, solar masses per cubic megaparsec."""
    mpc_m = _kpc_to_m(1.0) * 1000.0
    h0_per_s = 100.0 * HUBBLE_H * 1000.0 / mpc_m
    critical_kg_m3 = 3.0 * h0_per_s ** 2 / (8.0 * math.pi * constants.G)
    return OMEGA_MATTER * critical_kg_m3 * mpc_m ** 3 / _msun_to_kg(1.0)


def tinker_parameters(delta):
    """(A, a, b, c) at `delta` times the mean density, interpolated linearly in `log Delta` between the rows of Table 2.

    Raises:
        ValueError: If `delta` is outside the table (200 to 3200).
    """
    rows = TINKER_TABLE
    if not rows[0][0] <= delta <= rows[-1][0]:
        raise ValueError(f"Tinker et al. 2008 give no fit at Delta = {delta:g} (200 to 3200 times the mean density)")
    for low, high in zip(rows, rows[1:]):
        if low[0] <= delta <= high[0]:
            t = (math.log(delta) - math.log(low[0])) / (math.log(high[0]) - math.log(low[0]))
            return tuple(a + (b - a) * t for a, b in zip(low[1:], high[1:]))
    raise AssertionError("unreachable")


def _transfer_function(k_per_mpc):
    """Eisenstein and Hu 1998's no-wiggle fit to the CDM plus baryon transfer function (Eqs. 26, 28 to 31). [Stand-in
    HF2, tier T1, tweak L1.]"""
    om_h2 = OMEGA_MATTER * HUBBLE_H ** 2
    baryon_share = OMEGA_BARYON_H2 / om_h2
    sound_horizon = 44.5 * math.log(9.83 / om_h2) / math.sqrt(1.0 + 10.0 * OMEGA_BARYON_H2 ** 0.75)
    alpha = (1.0 - 0.328 * math.log(431.0 * om_h2) * baryon_share
             + 0.38 * math.log(22.3 * om_h2) * baryon_share ** 2)
    theta = CMB_TEMPERATURE_K / 2.7
    gamma = OMEGA_MATTER * HUBBLE_H * (alpha + (1.0 - alpha) / (1.0 + (0.43 * k_per_mpc * sound_horizon) ** 4))
    q = (k_per_mpc / HUBBLE_H) * theta ** 2 / gamma
    l0 = math.log(2.0 * math.e + 1.8 * q)
    c0 = 14.2 + 731.0 / (1.0 + 62.5 * q)
    return l0 / (l0 + c0 * q * q)


def _window(x):
    """The Fourier transform of a top-hat sphere, `3 (sin x - x cos x) / x^3`."""
    if x < 0.02:
        return 1.0 - x * x / 10.0
    return 3.0 * (math.sin(x) - x * math.cos(x)) / x ** 3


@functools.lru_cache(maxsize=1)
def _power_grid():
    """(k, k^3 P(k) / 2 pi^2 per unit `ln k`) on the log-spaced wave number grid, P unnormalised."""
    count = int(round((_LN_K_HIGH - _LN_K_LOW) / _LN_K_STEP)) + 1
    ks, weights = [], []
    for i in range(count):
        k = math.exp(_LN_K_LOW + i * _LN_K_STEP)
        ks.append(k)
        weights.append(k ** (3.0 + SPECTRAL_INDEX) * _transfer_function(k) ** 2 / (2.0 * math.pi ** 2))
    return tuple(ks), tuple(weights)


def _variance_unnormalised(radius_mpc):
    ks, weights = _power_grid()
    terms = [w * _window(k * radius_mpc) ** 2 for k, w in zip(ks, weights)]
    return _LN_K_STEP * (math.fsum(terms) - 0.5 * (terms[0] + terms[-1]))


@functools.lru_cache(maxsize=1)
def _normalisation():
    return SIGMA_8 ** 2 / _variance_unnormalised(8.0 / HUBBLE_H)


def sigma(halo_mass_solar):
    """The rms linear fluctuation of the density in a top-hat sphere that holds `halo_mass_solar` of mean matter."""
    radius = (3.0 * halo_mass_solar / (4.0 * math.pi * _mean_matter_density_msun_mpc3())) ** (1.0 / 3.0)
    return math.sqrt(_normalisation() * _variance_unnormalised(radius))


def mass_function_per_dex(log_halo_mass, delta):
    """Tinker et al. 2008's number density of halos, per cubic megaparsec per dex of mass, at `10 ** log_halo_mass`
    solar masses and `delta` times the mean density; the logarithmic slope of `sigma` is a central difference of
    0.025 dex."""
    step = 0.025
    s = sigma(10.0 ** log_halo_mass)
    dln_sigma_dlog10m = (math.log(sigma(10.0 ** (log_halo_mass + step)))
                         - math.log(sigma(10.0 ** (log_halo_mass - step)))) / (2.0 * step)
    a_norm, a, b, c = tinker_parameters(delta)
    f = a_norm * ((s / b) ** -a + 1.0) * math.exp(-c / s ** 2)
    return f * _mean_matter_density_msun_mpc3() / 10.0 ** log_halo_mass * -dln_sigma_dlog10m


@functools.lru_cache(maxsize=4)
def _ln_function_grid(delta):
    """`ln(dn/dlog10M)` on the fine grid, linearly interpolated between the coarse points."""
    coarse_count = int(round((LOG_MASS_HIGH - LOG_MASS_LOW) / _COARSE_STEP)) + 1
    coarse = [math.log(mass_function_per_dex(LOG_MASS_LOW + i * _COARSE_STEP, delta)) for i in range(coarse_count)]
    per_coarse = int(round(_COARSE_STEP / GRID_STEP))
    fine = []
    for i in range(coarse_count - 1):
        for j in range(per_coarse):
            fine.append(coarse[i] + (coarse[i + 1] - coarse[i]) * j / per_coarse)
    fine.append(coarse[-1])
    return tuple(fine)


def grid_log_masses():
    """The halo mass grid, log10 solar masses."""
    count = int(round((LOG_MASS_HIGH - LOG_MASS_LOW) / GRID_STEP)) + 1
    return tuple(LOG_MASS_LOW + i * GRID_STEP for i in range(count))


def halo_log_mass_at(log_stellar_mass, share, log_median_stellar_grid, scatter_dex, delta):
    """`log10` of the halo mass below which a fraction `share` of the halos that hold `10 ** log_stellar_mass` solar
    masses of stars lie.

    Args:
        log_stellar_mass (float): The galaxy's stellar mass, log10 solar masses.
        share (float): The cumulative share, 0 to 1 (0.5 is the median).
        log_median_stellar_grid (sequence): `log10` of the relation's median stellar mass at each point of
            `grid_log_masses()`.
        scatter_dex (float): The relation's scatter in `log M*` at fixed `Mh`.
        delta (float): The halo overdensity, times the mean matter density.

    Returns:
        float: `log10` of the halo mass, solar masses, interpolated on the grid.
    """
    ln_n = _ln_function_grid(delta)
    exponents = [n - 0.5 * ((log_stellar_mass - m) / scatter_dex) ** 2 for n, m in zip(ln_n, log_median_stellar_grid)]
    peak = max(exponents)
    weights = [math.exp(e - peak) for e in exponents]
    cumulative = [0.0]
    for i in range(1, len(weights)):
        cumulative.append(cumulative[-1] + 0.5 * (weights[i - 1] + weights[i]))
    target = min(max(share, 0.0), 1.0) * cumulative[-1]
    i = min(max(bisect.bisect_left(cumulative, target), 1), len(cumulative) - 1)
    span = cumulative[i] - cumulative[i - 1]
    fraction = 0.0 if span <= 0.0 else (target - cumulative[i - 1]) / span
    return LOG_MASS_LOW + (i - 1 + fraction) * GRID_STEP
