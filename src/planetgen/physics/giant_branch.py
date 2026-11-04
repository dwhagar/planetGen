# planetgen/physics/giant_branch.py

"""
The first giant branch of a star of mass `M` and metallicity (GEN.502): where it starts and ends, how long it takes and
how the star's luminosity is spread over that time.

It follows [R2055] (MNRAS 315, 543, [R2055], section 5.2 and the Appendix), whose
formulae `physics/ms_lifetime.py` already reads for the time to the base of the branch. On the branch the helium core
`Mc` grows by the hydrogen burning of the shell, `dMc/dt = A_H L`, and the luminosity is a power of the core mass,
`L = min(B Mc^q, D Mc^p)` (equation 37; `p = 6` and `q = 3` below the helium-flash mass, 5 and 2 above 2.5 solar
masses, `D`, `B` and `A_H` from the equations of section 5.2). Integrating gives the time to reach each luminosity
(equations 36, 40 to 43), so the branch's duration, its end at helium ignition `L_HeI` (equation 49, the tip: about
2,700 solar luminosities for a low-mass star of reference composition, 1,900 at Z = 1e-4 against the 1,650 to 1,780 of
[R4024]) and the time a giant spends at each luminosity (the slope flattening along the branch, from
-0.83 at the base to -0.56 at the tip) follow from the core-mass-luminosity relation instead of being set. The
duration reproduces the paper's Table 1: `(t_HeI - t_BGB) / t_BGB` is 6.4e-2 at one solar mass, 2.0e-2 at two and 2.4e-3
at five.

The base luminosity `L_BGB` is equation 10. Metallicity enters as `zeta = log10(Z / 0.02)`, with `Z` from the star's
[Fe/H] as in `ms_lifetime.zeta`; `fe_h` defaults to the reference composition, as in `ms_lifetime`.
"""

import functools
import math

from planetgen.physics import ms_lifetime

_BASE = {
    27: (9.511033e1, 6.819618e1, -1.045625e1, -1.474939e1),
    28: (3.113458e1, 1.012033e1, -4.650511, -2.463185),
    29: (1.413057, 4.578814e-1, -6.850581e-2, -5.588658e-2),  # a'29: a29 = a'29 ^ a32
    30: (3.910862e1, 5.196646e1, 2.264970e1, 2.873680),
    31: (4.597479, -2.855179e-1, 2.709724e-1),
    32: (6.682518, 2.827718e-1, -7.294429e-2),
}
"""dict: The Appendix coefficients of `L_BGB`, `alpha + beta zeta + gamma zeta^2 + eta zeta^3`, lowest power first."""
_BASE_C2 = 9.301992
_BASE_C3 = 4.637345
_TIP = {
    9: (2.751631e3, 3.557098e2),
    10: (-3.820831e-2, 5.872664e-2),
    11: (1.071738e2, -8.970339e1, -3.949739e1),  # b'11: b11 = b'11 ^ 2
    12: (7.348793e2, -1.531020e2, -3.793700e1),
    13: (9.219293, -2.005865, -5.561309e-1),  # b'13: b13 = b'13 ^ 2
}
"""dict: The Appendix coefficients of `L_HeI`, as `_BASE`."""
_TRANSITION_MASS_SOL = 2.5
"""float: Above it the intermediate-mass values of `p`, `q` and `D` hold; the paper interpolates linearly from `M_HeF`."""
_MYR_PER_GYR = 1.0e3


def _poly(coefficients, zeta):
    return sum(c * zeta ** power for power, c in enumerate(coefficients))


@functools.lru_cache(maxsize=256)
def _base_terms(zeta):
    """`a27 .. a32` of equation 10 at this `zeta` (cached: they depend on the metallicity alone, not on the mass)."""
    return tuple(_poly(_BASE[n], zeta) for n in (27, 28, 29, 30, 31, 32))


@functools.lru_cache(maxsize=256)
def _tip_terms(zeta):
    """`(b9, b10, b11, b12, b13)` of equation 49 at this `zeta`, `b11` and `b13` squared as the paper has them."""
    return (_poly(_TIP[9], zeta), _poly(_TIP[10], zeta), _poly(_TIP[11], zeta) ** 2, _poly(_TIP[12], zeta),
            _poly(_TIP[13], zeta) ** 2)


def helium_flash_mass_sol(fe_h=0.0):
    """The heaviest initial mass whose helium ignites degenerately, in a flash (equation 2), solar masses."""
    zeta = ms_lifetime.zeta(fe_h)
    return 1.995 + 0.25 * zeta + 0.087 * zeta * zeta


def base_luminosity_sol(mass_sol, fe_h=0.0):
    """The luminosity at the base of the giant branch (equation 10), solar luminosities."""
    zeta = ms_lifetime.zeta(fe_h)
    a27, a28, a29p, a30, a31, a32 = _base_terms(zeta)
    m = mass_sol
    return ((a27 * m ** a31 + a28 * m ** _BASE_C2)
            / (a29p ** a32 + a30 * m ** _BASE_C3 + m ** a32))


def tip_luminosity_sol(mass_sol, fe_h=0.0):
    """The luminosity at helium ignition, where the branch ends (equation 49), solar luminosities."""
    zeta = ms_lifetime.zeta(fe_h)
    flash_mass = helium_flash_mass_sol(fe_h)
    b9, b10, b11, b12, b13 = _tip_terms(zeta)

    def intermediate(m):
        return (b11 + b12 * m ** 3.8) / (b13 + m * m)

    if mass_sol >= flash_mass:
        return intermediate(mass_sol)
    alpha = (b9 * flash_mass ** b10 - intermediate(flash_mass)) / intermediate(flash_mass)
    return b9 * mass_sol ** b10 / (1.0 + alpha * math.exp(15.0 * (mass_sol - flash_mass)))


class Branch:
    """
    The giant branch of one initial mass and [Fe/H]: `base` and `tip` luminosities (solar), `start_myr` and `end_myr`
    (the ages the branch begins and ends), and the core-mass-luminosity law that spreads the time between them.
    """

    __slots__ = ("mass_sol", "base", "tip", "start_myr", "end_myr", "_p", "_q", "_d", "_b", "_rate", "_crossover",
                 "_inf_1", "_inf_2", "_cross_myr")

    def __init__(self, mass_sol, fe_h=0.0):
        zeta = ms_lifetime.zeta(fe_h)
        flash_mass = helium_flash_mass_sol(fe_h)
        d0 = 5.37 + 0.135 * zeta
        if mass_sol <= flash_mass:
            p, q, log_d = 6.0, 3.0, d0
        else:
            heavy_log_d = max(-1.0, 0.975 * d0 - 0.18 * max(mass_sol, _TRANSITION_MASS_SOL),
                              0.5 * d0 - 0.06 * max(mass_sol, _TRANSITION_MASS_SOL))
            weight = min((mass_sol - flash_mass) / (_TRANSITION_MASS_SOL - flash_mass), 1.0)
            p, q = 6.0 - weight, 3.0 - weight
            log_d = d0 + (heavy_log_d - d0) * weight
        d = 10.0 ** log_d
        b = max(3.0e4, 500.0 + 1.75e4 * mass_sol ** 0.6)
        rate = 10.0 ** max(-4.8, min(-5.7 + 0.8 * mass_sol, -4.1 + 0.14 * mass_sol))
        base, tip = base_luminosity_sol(mass_sol, fe_h), tip_luminosity_sol(mass_sol, fe_h)
        crossover = d * (b / d) ** (p / (p - q))
        start = ms_lifetime.base_of_giant_branch_gy(mass_sol, fe_h) * _MYR_PER_GYR
        inf_1 = start + (d / base) ** ((p - 1.0) / p) / ((p - 1.0) * rate * d)
        cross = inf_1 - (inf_1 - start) * (base / crossover) ** ((p - 1.0) / p)
        inf_2 = cross + (b / crossover) ** ((q - 1.0) / q) / ((q - 1.0) * rate * b)
        self.mass_sol, self.base, self.tip, self.start_myr = mass_sol, base, tip, start
        self._p, self._q, self._d, self._b, self._rate, self._crossover = p, q, d, b, rate, crossover
        self._inf_1, self._inf_2, self._cross_myr = inf_1, inf_2, cross
        self.end_myr = self._age_at(tip)

    def _age_at(self, luminosity):
        """The age, Myr, at which the star on the branch has this luminosity (equations 36 and 40 to 43)."""
        if luminosity <= self._crossover:
            p = self._p
            return self._inf_1 - (self._d / luminosity) ** ((p - 1.0) / p) / ((p - 1.0) * self._rate * self._d)
        q = self._q
        return self._inf_2 - (self._b / luminosity) ** ((q - 1.0) / q) / ((q - 1.0) * self._rate * self._b)

    def _luminosity_at(self, age_myr):
        """The luminosity at this age (equation 35 for each of the two parts of the law)."""
        if age_myr <= self._cross_myr:
            p = self._p
            return self._d * ((p - 1.0) * self._rate * self._d * (self._inf_1 - age_myr)) ** (p / (1.0 - p))
        q = self._q
        return self._b * ((q - 1.0) * self._rate * self._b * (self._inf_2 - age_myr)) ** (q / (1.0 - q))

    @property
    def duration_myr(self):
        return self.end_myr - self.start_myr

    def share_above(self, min_luminosity_sol):
        """The part of the branch's time spent at or above this luminosity; 1.0 for `None`."""
        if min_luminosity_sol is None or min_luminosity_sol <= self.base:
            return 1.0
        if min_luminosity_sol >= self.tip:
            return 0.0
        return (self.end_myr - self._age_at(min_luminosity_sol)) / self.duration_myr

    def quantile(self, fraction, min_luminosity_sol=None):
        """The luminosity below which `fraction` (0 to 1) of the time spent at or above `min_luminosity_sol` (`None`:
        from the base) is spent."""
        fraction = min(max(fraction, 0.0), 1.0)
        start = 0.0 if min_luminosity_sol is None or min_luminosity_sol <= self.base \
            else min(max((self._age_at(min_luminosity_sol) - self.start_myr) / self.duration_myr, 0.0), 1.0)
        progress = start + (1.0 - start) * fraction
        if progress >= 1.0:
            return self.tip
        return min(max(self._luminosity_at(self.start_myr + progress * self.duration_myr), self.base), self.tip)


@functools.lru_cache(maxsize=1 << 17)
def branch(mass_sol, fe_h=0.0):
    """The `Branch` of this mass and [Fe/H] (cached: the sampler asks for one mass again and again)."""
    return Branch(mass_sol, fe_h)


def subgiant_end_fraction(mass_sol, fe_h=0.0):
    """The age the star leaves the Hertzsprung gap for the giant branch, as a multiple of its main-sequence lifetime."""
    return ms_lifetime.base_of_giant_branch_gy(mass_sol, fe_h) / ms_lifetime.lifetime_gy(mass_sol, fe_h)


def giant_end_fraction(mass_sol, fe_h=0.0):
    """The age of helium ignition, the branch's end, as a multiple of the main-sequence lifetime."""
    return branch(mass_sol, fe_h).end_myr / _MYR_PER_GYR / ms_lifetime.lifetime_gy(mass_sol, fe_h)
