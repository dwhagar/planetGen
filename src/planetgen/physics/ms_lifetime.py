# planetgen/physics/ms_lifetime.py

"""
The main-sequence lifetime of a star from its mass and metallicity (GEN.467).

It replaces the fixed `10 Gy * M^-2.5`, a law anchored on the Sun, with the analytic fit of Hurley, Pols and Tout
2000 (MNRAS 315, 543, "Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity",
[R2055]), their equations 4 to 7 and the Appendix coefficients: the time to the base of the giant
branch `t_BGB(M, Z)` (rms 1.9%, at most 4.8% off the detailed models over 0.1 to 100 solar masses and Z = 0.0001 to
0.03), and the main sequence ends at `t_MS = max(t_hook, x t_BGB)`, a share `x` of that time (0.95 to 0.99) or a share
`mu` where a hook feature makes it shorter. Above 100 solar masses the fit is carried on as it stands: it flattens to
about 3 Myr, the lifetime the very massive stars have.

Metallicity enters as `zeta = log10(Z / 0.02)`. A star's heavy-element mass fraction is `Z = Z_ref * 10^[Fe/H]`, with
`tuning.MS_LIFETIME_REFERENCE_Z` the proto-solar value of the composition scale ([R0168], `stellar-composition.md`
D1); every caller that has no [Fe/H] yet takes the default 0, a star of the reference composition, for which a solar
mass lasts 10.2 Gyr.
"""

import functools
import math

from planetgen import tuning

HPT_REFERENCE_Z = 0.02
"""float: The `Z` the formulae's `zeta` is measured from [R2055]."""

# The Appendix coefficients, `a_n = alpha + beta zeta + gamma zeta^2 + eta zeta^3`, lowest power first; `a5`, `a7` and
# `a10` do not depend on metallicity. Times are in Myr.
_COEFFICIENTS = {
    1: (1.593890e3, 2.053038e3, 1.231226e3, 2.327785e2),
    2: (2.706708e3, 1.483131e3, 5.772723e2, 7.411230e1),
    3: (1.466143e2, -1.048442e2, -6.795374e1, -1.391127e1),
    4: (4.141960e-2, 4.564888e-2, 2.958542e-2, 5.571483e-3),
    5: (3.426349e-1,),
    6: (1.949814e1, 1.758178, -6.008212, -4.470533),
    7: (4.903830,),
    8: (5.212154e-2, 3.166411e-2, -2.750074e-3, -2.271549e-3),
    9: (1.312179, -3.294936e-1, 9.231860e-2, 2.610989e-2),
    10: (8.073972e-1,),
}
MIN_MASS_SOL = 0.05
MAX_MASS_SOL = 300.0
"""float: The range `mass_for_lifetime_gy` searches; a lifetime outside what these two masses give returns the end."""
_MYR_PER_GYR = 1.0e3
_BISECTIONS = 50


@functools.lru_cache(maxsize=4096)
@functools.lru_cache(maxsize=1 << 14)
def _terms(zeta):
    """`a1` to `a10` at this `zeta`, and the share `x` of the base-of-giant-branch time the main sequence takes."""
    terms = tuple(sum(c * zeta ** power for power, c in enumerate(_COEFFICIENTS[n])) for n in range(1, 11))
    x = max(0.95, min(0.95 - 0.03 * (zeta + 0.30103), 0.99))
    return terms, x


@functools.lru_cache(maxsize=1 << 14)
def zeta(fe_h):
    """`log10(Z / 0.02)` of a star of this [Fe/H]: the metallicity argument of every formula of the paper, held to the
    range the paper's fits are made over (`tuning.MS_LIFETIME_Z_RANGE`): below Z = 1e-4 the polynomials in it run away
    (a 0.8 solar-mass star at [Fe/H] -3.25 would live 5 Gyr and a table of bright stars would be 50 times as big there
    as anywhere), so a more metal-poor star has the lifetimes and branches of Z = 1e-4."""
    low, high = (math.log10(z / HPT_REFERENCE_Z) for z in tuning.MS_LIFETIME_Z_RANGE)
    return min(max(math.log10(tuning.MS_LIFETIME_REFERENCE_Z / HPT_REFERENCE_Z) + fe_h, low), high)


def _base_of_giant_branch_gy(terms, mass_sol):
    """`base_of_giant_branch_gy` from the `terms` of its [Fe/H]."""
    a1, a2, a3, a4, a5 = terms[:5]
    m = mass_sol
    return (a1 + a2 * m ** 4 + a3 * m ** 5.5 + m ** 7) / (a4 * m ** 2 + a5 * m ** 7) / _MYR_PER_GYR


def base_of_giant_branch_gy(mass_sol, fe_h=0.0):
    """The age at which a star of this mass leaves the Hertzsprung gap for the giant branch (eq. 4), Gyr."""
    return _base_of_giant_branch_gy(_terms(zeta(fe_h))[0], mass_sol)


def _lifetime_gy(terms, x, mass_sol):
    """`lifetime_gy` from the `terms` and `x` of its [Fe/H]."""
    _a1, _a2, _a3, _a4, _a5, a6, a7, a8, a9, a10 = terms
    hook = max(0.5, 1.0 - 0.01 * max(a6 / mass_sol ** a7, a8 + a9 / mass_sol ** a10))
    return max(hook, x) * _base_of_giant_branch_gy(terms, mass_sol)


def lifetime_gy(mass_sol, fe_h=0.0):
    """The main-sequence lifetime of a star of this initial mass and [Fe/H], Gyr (eq. 5)."""
    terms, x = _terms(zeta(fe_h))
    return _lifetime_gy(terms, x, mass_sol)


def mass_for_lifetime_gy(lifetime, fe_h=0.0):
    """
    The initial mass whose main-sequence lifetime is `lifetime` Gyr: the lifetime falls with mass over the whole
    range, so a bisection in log mass finds it. A lifetime longer than the lightest mass (`MIN_MASS_SOL`) lives gives
    that mass, a shorter one than the heaviest (`MAX_MASS_SOL`) the heaviest.
    """
    low, high = math.log(MIN_MASS_SOL), math.log(MAX_MASS_SOL)
    terms, x = _terms(zeta(fe_h))
    if lifetime >= _lifetime_gy(terms, x, MIN_MASS_SOL):
        return MIN_MASS_SOL
    if lifetime <= _lifetime_gy(terms, x, MAX_MASS_SOL):
        return MAX_MASS_SOL
    for _ in range(_BISECTIONS):
        middle = 0.5 * (low + high)
        if _lifetime_gy(terms, x, math.exp(middle)) > lifetime:
            low = middle
        else:
            high = middle
    return math.exp(0.5 * (low + high))
