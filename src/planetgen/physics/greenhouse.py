# planetgen/physics/greenhouse.py

"""
Greenhouse closure: the surface temperature a grey infrared opacity gives, and the opacity a surface temperature needs
(GEN.443, docs/design/composition/radiation-model.md section 4; docs/design/composition/atmosphere-composition-parts/03-infrared-greenhouse-thermal.md 5.3 and 5.4)
==========================================================================================================================================================

The repo's greenhouse variable is `G = (T_s / T_e)^4 - 1`, with `T_e` the effective temperature of the absorbed flux.
Robinson and Catling (2012, 2014) give a closed analytic atmosphere for it: an optical depth that rises with pressure as
`tau(p) = tau0 (p / p0)^n`, a convective troposphere below the radiative-convective boundary `tau_rc`, a radiative
stratosphere above it, the diffusivity factor `D = 1.66` and the lapse-rate scaling `beta = alpha (gamma - 1) / gamma`
(`gamma` the heat-capacity ratio, `alpha` 0.6 to 0.9). With a purely radiative atmosphere `T_0^4 = T_e^4 (1 + D tau0 /
2)`. Otherwise the radiative-convective boundary is where temperature and upwelling flux are continuous (their Eqs. 27
and 30), a root found numerically (`radiative_convective_boundary`).

* `surface_temperature_k(tau0, absorbed_flux_w_m2, gamma, alpha, n)` returns `(T_0, tau_rc)` (`tau_rc` is `None` for a
  purely radiative atmosphere).
* `tau0_implied(T_s, absorbed_flux_w_m2, gamma, alpha, n)` inverts it (bisection; `T_0` rises with `tau0`).
* `check_stored_multiplier(...)` is the check of section 4: from the stored albedo, flux and surface temperature it
  returns `T_e`, `G`, the `tau0` the closure needs, and the `tau0` the stored multiplier would need if it were the
  surface temperature. The composition side of the check (the tau0 of the CIA floor, the water continuum and the lines)
  is not here: it waits for the infrared opacity data (R2, R10) and lane D's stored multiplier.

Absorbed flux is the globally averaged flux the planet absorbs, W/m2: `(1 - A) S / 4` for a planet that spreads the
heat. Limits (stated): a grey single-band closure; the pressure exponent `n` is 1 for Doppler- or line-wing-dominated
opacity and 2 for pressure-broadened or collision-induced; the continuity root is unique in tested ranges but a thin
atmosphere (tau0 under about 0.6) gives a radiative closure and a continuity root that differ by a factor of 2 in
`tau0` (R11): `surface_temperature_k` returns the continuity answer when it exists and the radiative one otherwise.
"""

import math

import numpy as np

from planetgen.physics import constants

DIFFUSIVITY = 1.66
"""float: The diffusivity factor D of Robinson and Catling (2012)."""
STEPS = 400
"""int: Simpson steps of the continuity integral (an even number)."""
INTEGRAL_REACH = 60.0
"""float: How many e-folds of `exp(-D s)` the continuity integral runs over (the rest is under 1e-26 of it)."""
TAU0_MIN = 1e-6
TAU0_MAX = 1e6
"""float: The range `tau0_implied` searches."""


def effective_temperature_k(absorbed_flux_w_m2):
    """`T_e`, K, the temperature that radiates the absorbed flux: `(F / sigma)^(1/4)`."""
    return (absorbed_flux_w_m2 / constants.STEFAN_BOLTZMANN_CONSTANT) ** 0.25


def absorbed_flux_w_m2(stellar_flux_w_m2, bond_albedo, redistribution=0.25):
    """The flux a planet absorbs averaged over its surface, W/m2: `(1 - A) S f`, `f` 1/4 for a planet that spreads the
    heat (a tidally locked world with a hot dayside uses a larger `f`)."""
    return (1.0 - bond_albedo) * stellar_flux_w_m2 * redistribution


def greenhouse_factor(surface_temperature_k, effective_k):
    """`G = (T_s / T_e)^4 - 1`, the variable the repo stores."""
    return (surface_temperature_k / effective_k) ** 4 - 1.0


def _exponent(gamma, alpha, n):
    return 4.0 * alpha * (gamma - 1.0) / gamma / n


def _boundary_residual(t, tau0, m):
    """The continuity condition of Robinson and Catling at optical depths `t` (an array; zero at `tau_rc`)."""
    t = np.asarray(t, dtype=float)
    # the integrand falls as exp(-D (s - t)), so past t + REACH / D it adds nothing; integrating the whole of a very thick atmosphere
    # with a fixed number of steps would step over the part that counts
    upper = np.minimum(tau0, t + INTEGRAL_REACH / DIFFUSIVITY)
    fractions = np.linspace(0.0, 1.0, STEPS + 1)
    s = t[..., None] + (upper - t)[..., None] * fractions
    weights = np.full(STEPS + 1, 2.0)
    weights[1::2] = 4.0
    weights[0] = weights[-1] = 1.0
    integrand = (s / tau0) ** m * np.exp(-DIFFUSIVITY * (s - t[..., None]))
    integral = (integrand * weights).sum(axis=-1) * (upper - t) / STEPS / 3.0
    return ((np.exp(-DIFFUSIVITY * (tau0 - t)) + DIFFUSIVITY * integral) / (t / tau0) ** m
            - (2.0 + DIFFUSIVITY * t) / (1.0 + DIFFUSIVITY * t))


def radiative_convective_boundary(tau0, gamma, alpha, n):
    """The optical depth `tau_rc` where the convective troposphere meets the radiative stratosphere, or `None` when the
    atmosphere is radiative all the way (no sign change of the continuity residual between 1e-6 and `tau0`)."""
    m = _exponent(gamma, alpha, n)
    steps = 200
    lo = math.log10(1e-6)
    hi = math.log10(tau0)
    if hi <= lo:
        return None
    xs = 10 ** (lo + (hi - lo) * np.arange(steps + 1) / steps)
    values = _boundary_residual(xs, tau0, m)
    changes = np.nonzero(values[:-1] * values[1:] < 0)[0]
    if not len(changes):
        return None
    index = int(changes[0])
    a, b, fa = float(xs[index]), float(xs[index + 1]), float(values[index])
    for _ in range(80):
        mid = 0.5 * (a + b)
        fm = float(_boundary_residual(np.array([mid]), tau0, m)[0])
        if fa * fm <= 0:
            b = mid
        else:
            a, fa = mid, fm
        if b - a < 1e-12 * a:
            break
    return 0.5 * (a + b)


def surface_temperature_k(tau0, absorbed_flux_w_m2_, gamma, alpha, n):
    """`(T_0, tau_rc)`: the surface temperature, K, of a grey atmosphere of surface optical depth `tau0` under an
    absorbed flux, and its radiative-convective boundary (`None` for a radiative atmosphere)."""
    if tau0 <= 0:
        return effective_temperature_k(absorbed_flux_w_m2_), None
    boundary = radiative_convective_boundary(tau0, gamma, alpha, n)
    if boundary is None:
        return (absorbed_flux_w_m2_ / constants.STEFAN_BOLTZMANN_CONSTANT * (1.0 + DIFFUSIVITY * tau0 / 2.0)) ** 0.25, None
    m = _exponent(gamma, alpha, n)
    flux = absorbed_flux_w_m2_ * (1.0 + DIFFUSIVITY * boundary) / 2.0
    return (flux / constants.STEFAN_BOLTZMANN_CONSTANT / (boundary / tau0) ** m) ** 0.25, boundary


def tau0_implied(surface_temperature_k_, absorbed_flux_w_m2_, gamma, alpha, n):
    """The `tau0` for which the closure gives `surface_temperature_k_`. Raises `ValueError` when the temperature is below
    the effective temperature (no greenhouse can cool a surface) or above what `TAU0_MAX` reaches."""
    if surface_temperature_k_ <= effective_temperature_k(absorbed_flux_w_m2_):
        raise ValueError("a surface at or below the effective temperature needs no greenhouse (tau0 = 0)")
    lo, hi = TAU0_MIN, TAU0_MAX
    if surface_temperature_k(hi, absorbed_flux_w_m2_, gamma, alpha, n)[0] < surface_temperature_k_:
        raise ValueError("the surface temperature is beyond what a tau0 of %g reaches" % TAU0_MAX)
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        if surface_temperature_k(mid, absorbed_flux_w_m2_, gamma, alpha, n)[0] < surface_temperature_k_:
            lo = mid
        else:
            hi = mid
        if hi / lo < 1 + 1e-9:
            break
    return math.sqrt(lo * hi)


def check_stored_multiplier(surface_temperature_k_, absorbed_flux_w_m2_, gamma, alpha, n, stored_multiplier=None):
    """The check of a stored greenhouse (radiation-model section 4, steps 1 and 2): `T_e`, `G` from the stored
    temperature, the `tau0` it implies and its `tau_rc`; with `stored_multiplier` (a `G`) also the temperature that
    multiplier gives and the `tau0` that temperature implies."""
    effective = effective_temperature_k(absorbed_flux_w_m2_)
    result = {"effective_temperature_k": effective,
              "greenhouse_factor": greenhouse_factor(surface_temperature_k_, effective)}
    tau0 = tau0_implied(surface_temperature_k_, absorbed_flux_w_m2_, gamma, alpha, n)
    result["tau0"] = tau0
    result["tau_rc"] = surface_temperature_k(tau0, absorbed_flux_w_m2_, gamma, alpha, n)[1]
    if stored_multiplier is not None:
        stored_k = effective * (1.0 + stored_multiplier) ** 0.25
        result["stored_surface_temperature_k"] = stored_k
        result["stored_tau0"] = tau0_implied(stored_k, absorbed_flux_w_m2_, gamma, alpha, n)
    return result
