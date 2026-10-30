# planetgen/physics/resonant_chains.py

"""
Resonant chains: which planets lock their periods in a row, and the eccentricity that forces on them (GEN.1043)
==============================================================================================================

Three or more planets whose neighbouring periods sit near ratios of small whole numbers form a chain: TRAPPIST-1's seven planets (8:5, 5:3, 3:2, 3:2, 4:3, 3:2), the four
of Kepler-223 (4:3, 3:2, 4:3), Io, Europa and Ganymede (2:1, 2:1). A planet that a neighbour tugs at near a resonance does not keep a circular orbit: it is given a
**forced eccentricity**, and a planet on an eccentric orbit is flexed by its star's tide (`tidal_heating`, the Peale, Cassen and Reynolds 1979 power Io shows). This module
finds the chains in a set of periods and gives the forced eccentricity and the orbits for `tidal_heating` to flex. Nothing calls it yet (GEN.1049 wires it into the
architecture and the interior heat, `weird-worlds.md` W5).

**The pair.** Planets in order of period; for the period ratio `P'/P` of a neighbouring pair the nearest of `TARGET_RATIOS` (2:1, 3:2, 4:3, 5:3) is taken if the ratio is within
`RESONANCE_TOLERANCE` (5 percent) of it. The distance from exact resonance is `Delta = (P'/P) (q/p) - 1` for the ratio `p:q`; the Kepler multis cluster just wide of exact,
`Delta` of about +0.01 to +0.02 (Fabrycky et al. 2014, ApJ 790:146; Goldberg and Batygin 2022, AJ 163:201; recalled, not opened). A **chain** is a run of three or more planets all
of whose neighbouring pairs are resonant.

**The forced eccentricity** of a first-order resonance (p - q = 1: 2:1, 3:2, 4:3), the planets not yet in the resonance's own libration, from the secular-free resonant terms of the
disturbing function (Murray and Dermott 1999, *Solar System Dynamics*, section 8.8 and Appendix B; Lithwick, Xie and Wu 2012, ApJ 761:122; recalled, derived here):

    e_inner = mu' alpha |f_27| / ((p - 1) |Delta|),      e_outer = mu |f_31| / (p |Delta|)

with `mu`, `mu'` the inner and outer planets' masses over their star's, `alpha` the ratio of the semi-major axes, and `f_27`, `f_31` the resonance's coefficients in Laplace
coefficients of `alpha` (`inner_coefficient`, `outer_coefficient`: -1.1905 and +0.4284 for 2:1, -2.0252 and +2.4840 for 3:2, the published table). For Io pushed by Europa it gives
0.0052 where the observed forced value is 0.0041; the difference is Io's precession (Jupiter's flattening, which adds to the detuning and is not in the formula). A
second-order resonance (5:3) forces nothing at this order: its forced eccentricity is proportional to a free one, so it is zero here. A planet with two neighbours is forced by
both, the two added in quadrature (their phases are independent). Near exact resonance the denominator would vanish; the resonance has a width, so `|Delta|` is not taken
below `mu_max ** WIDTH_EXPONENT` (the near-resonance approximation of Lithwick et al. holds for `|Delta|` above `mu^(2/3)`), which bounds the result at about `mu^(1/3)`. The
eccentricity is also held to `CROSSING_SHARE` of the most that keeps two neighbouring orbits from crossing, `(a' - a) / (a' + a)`.

**What it does not do.** No libration (a planet inside the resonance's width keeps a fixed eccentricity set by its tidal history, not by the formula), no second-order or
three-body resonances, no stability integration (the crossing limit is the necessary condition), and the planets are taken as having no other perturber. The eccentricity
feeds `tidal_heating.forcing_of`-style `Forcing` objects (`forcing`) and from there the solid, melt-limited or icy power as for a moon. Stand-ins are in
`docs/design/research/gap-filling-plan.md` section 6.92.
"""

import itertools
import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import constants, tidal_heating

TARGET_RATIOS = ((2, 1), (3, 2), (4, 3), (5, 3))
"""tuple: The period ratios `p:q` a pair is taken to be near (the design note's list). Source tier T6, tweak level L0."""
RESONANCE_TOLERANCE = 0.05
"""float: How near a ratio must be to a target, as a share of it (5 percent, the design note's gate). Source tier T6, tweak level L0."""
WIDTH_EXPONENT = 2.0 / 3.0
"""float: The detuning is not taken below `mu ** WIDTH_EXPONENT` (Lithwick, Xie and Wu 2012's limit of the near-resonance formula; the half-width of a first-order resonance is of
that order). Source tier T4, tweak level L0."""
CROSSING_SHARE = 0.5
"""float: The forced eccentricity of a planet is held to this share of the most that keeps its orbit from crossing a neighbour's, `(a' - a) / (a' + a)` for the sum of the two.
Source tier T6, tweak level L1."""
_INTEGRAL_POINTS = 4096


class ResonantChainError(ValueError):
    """An input the model cannot use."""


@dataclass(frozen=True)
class Pair:
    """
    A neighbouring pair near a resonance.

    Attributes:
        inner, outer (int): Their places in the order of period.
        ratio (float): `P_outer / P_inner`.
        p, q (int): The target `p:q`.
        delta (float): The distance from exact resonance, `ratio * q / p - 1`.
    """

    inner: int
    outer: int
    ratio: float
    p: int
    q: int
    delta: float

    @property
    def order(self):
        """`p - q`: 1 for a first-order resonance."""
        return self.p - self.q


@dataclass(frozen=True)
class Chain:
    """
    Three or more planets whose neighbours are all near resonances.

    Attributes:
        members (tuple): Their places in the order of period.
        pairs (tuple): The `Pair` between each neighbour.
    """

    members: tuple
    pairs: tuple


def _laplace(order, alpha, derivative=False):
    """`b_{1/2}^{(order)}(alpha)` or its derivative in `alpha`, by the trapezoid rule on the full period (spectrally exact for this smooth integrand)."""
    psi = (np.arange(_INTEGRAL_POINTS) + 0.5) / _INTEGRAL_POINTS * 2.0 * math.pi
    base = 1.0 - 2.0 * alpha * np.cos(psi) + alpha * alpha
    if derivative:
        return float(2.0 / _INTEGRAL_POINTS * np.sum(np.cos(order * psi) * -0.5 * (2.0 * alpha - 2.0 * np.cos(psi)) / base ** 1.5))
    return float(2.0 / _INTEGRAL_POINTS * np.sum(np.cos(order * psi) / base ** 0.5))


def inner_coefficient(p, alpha):
    """`f_27`, the coefficient of the inner planet's eccentricity in the first-order `p:(p-1)` term (Murray and Dermott Appendix B): `-1/2 (2p + alpha D) b_{1/2}^{(p)}`."""
    return -0.5 * (2.0 * p * _laplace(p, alpha) + alpha * _laplace(p, alpha, True))


def outer_coefficient(p, alpha):
    """`f_31`, the same for the outer planet: `1/2 ((2p - 1) + alpha D) b_{1/2}^{(p-1)}`, less the indirect term `2 alpha` for 2:1."""
    value = 0.5 * ((2.0 * p - 1.0) * _laplace(p - 1, alpha) + alpha * _laplace(p - 1, alpha, True))
    return value - 2.0 * alpha if p == 2 else value


def classify_pair(period_inner, period_outer):
    """The `(p, q, delta)` of the target the ratio of the periods is nearest to within `RESONANCE_TOLERANCE`, or `None`."""
    ratio = period_outer / period_inner
    best = None
    for p, q in TARGET_RATIOS:
        miss = abs(ratio * q / p - 1.0)
        if miss <= RESONANCE_TOLERANCE and (best is None or miss < best[0]):
            best = (miss, p, q, ratio * q / p - 1.0)
    return None if best is None else best[1:]


def _order(periods_days):
    periods = [float(p) for p in periods_days]
    if any(not math.isfinite(p) or p <= 0.0 for p in periods):
        raise ResonantChainError("periods must be positive")
    return sorted(range(len(periods)), key=lambda i: periods[i]), periods


def find_pairs(periods_days):
    """The resonant neighbouring `Pair`s among the periods (days, in any order); `inner` and `outer` index the input list."""
    order, periods = _order(periods_days)
    out = []
    for a, b in itertools.pairwise(order):
        found = classify_pair(periods[a], periods[b])
        if found is not None:
            p, q, delta = found
            out.append(Pair(a, b, periods[b] / periods[a], p, q, delta))
    return tuple(out)


def find_chains(periods_days):
    """The `Chain`s, three planets or more, among the periods (days, any order); members index the input list, innermost first."""
    pairs = find_pairs(periods_days)
    chains, run = [], []
    for pair in pairs:
        if run and run[-1].outer == pair.inner:
            run.append(pair)
        else:
            if len(run) >= 2:
                chains.append(run)
            run = [pair]
    if len(run) >= 2:
        chains.append(run)
    return tuple(Chain(tuple([run[0].inner] + [pair.outer for pair in run]), tuple(run)) for run in chains)


def semi_major_axis_m(period_days, host_mass_kg, planet_mass_kg=0.0):
    """The semi-major axis of an orbit of the period round a host (Kepler's third law with both masses), m."""
    return (constants.G * (host_mass_kg + planet_mass_kg) * (period_days * 86400.0 / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)


def _pair_eccentricities(pair, periods, masses, host_mass_kg):
    """`(e_inner, e_outer)` forced by one pair (zero for a second-order resonance), before the crossing limit."""
    if pair.order != 1:
        return 0.0, 0.0
    a_in = semi_major_axis_m(periods[pair.inner], host_mass_kg, masses[pair.inner])
    a_out = semi_major_axis_m(periods[pair.outer], host_mass_kg, masses[pair.outer])
    alpha = a_in / a_out
    mu_in, mu_out = masses[pair.inner] / host_mass_kg, masses[pair.outer] / host_mass_kg
    width = max(mu_in, mu_out) ** WIDTH_EXPONENT
    delta = max(abs(pair.delta), width)
    inner = mu_out * alpha * abs(inner_coefficient(pair.p, alpha)) / ((pair.p - 1) * delta)
    outer = mu_in * abs(outer_coefficient(pair.p, alpha)) / (pair.p * delta)
    return inner, outer


def forced_eccentricities(periods_days, masses_kg, host_mass_kg):
    """
    The forced eccentricity of each planet in a chain, in the order of the input (0.0 for a planet in no chain).

    Args:
        periods_days (sequence): The orbital periods, days.
        masses_kg (sequence): The planets' masses, kg, in the same order.
        host_mass_kg (float): The star's mass.
    """
    if len(periods_days) != len(masses_kg):
        raise ResonantChainError("a mass for each period")
    if not host_mass_kg > 0.0 or any(not m > 0.0 for m in masses_kg):
        raise ResonantChainError("masses must be positive")
    periods = [float(p) for p in periods_days]
    squares = [0.0] * len(periods)
    for chain in find_chains(periods):
        for pair in chain.pairs:
            inner, outer = _pair_eccentricities(pair, periods, masses_kg, host_mass_kg)
            squares[pair.inner] += inner ** 2
            squares[pair.outer] += outer ** 2
        for pair in chain.pairs:
            a_in = semi_major_axis_m(periods[pair.inner], host_mass_kg, masses_kg[pair.inner])
            a_out = semi_major_axis_m(periods[pair.outer], host_mass_kg, masses_kg[pair.outer])
            limit = CROSSING_SHARE * (a_out - a_in) / (a_out + a_in)
            for member in (pair.inner, pair.outer):
                squares[member] = min(squares[member], limit ** 2)
    return tuple(math.sqrt(s) for s in squares)


def forcing(host_mass_kg, period_days, planet_mass_kg, eccentricity):
    """The `tidal_heating.Forcing` of a planet on the orbit of the period with the eccentricity, to give `solid_power_w`, `molten_power_w` or `icy_flux_w_m2` as for a moon."""
    return tidal_heating.Forcing(host_mass_kg, semi_major_axis_m(period_days, host_mass_kg, planet_mass_kg), eccentricity)
