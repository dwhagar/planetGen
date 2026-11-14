# planetgen/physics/resonant_chains.py

"""
Resonant chains: which planets lock their periods in a row, and the eccentricity that forces on them (GEN.1043)
==============================================================================================================

Three or more planets whose neighbouring periods sit near ratios of small whole numbers form a chain: TRAPPIST-1's seven planets (8:5, 5:3, 3:2, 3:2, 4:3, 3:2), the four
of Kepler-223 (4:3, 3:2, 4:3), Io, Europa and Ganymede (2:1, 2:1). A planet that a neighbour tugs at near a resonance does not keep a circular orbit: it is given a
**forced eccentricity**, and a planet on an eccentric orbit is flexed by its star's tide (`tidal_heating`, the [R3464] power Io shows). This module
finds the chains in a set of periods and gives the forced eccentricity and the orbits for `tidal_heating` to flex. Oddity W5 (GEN.1049) reads it: `give_forced_eccentricities`
hands each planet of a chain its forced eccentricity as a quantity of the body, and the `resonance_forced` effect raises the planet's eccentricity to it for the rows that roll
(about a fifth of the chains' planets, WR8), which the tide of the orbit then flexes into heat (`planets.finish_orbit_oddities` makes the planet's end stages again).

**The pair.** Planets in order of period; for the period ratio `P'/P` of a neighbouring pair the nearest of `TARGET_RATIOS` (2:1, 3:2, 4:3, 5:3) is taken if the ratio is within
`RESONANCE_TOLERANCE` (5 percent) of it. The distance from exact resonance is `Delta = (P'/P) (q/p) - 1` for the ratio `p:q`; the Kepler multis cluster just wide of exact,
`Delta` of about +0.01 to +0.02 ([R1311]:146; [R1605]:201; recalled, not opened). A **chain** is a run of three or more planets all
of whose neighbouring pairs are resonant.

**The forced eccentricity** of a first-order resonance (p - q = 1: 2:1, 3:2, 4:3), the planets not yet in the resonance's own libration, from the secular-free resonant terms of the
disturbing function ([R3220], *Solar System Dynamics*, section 8.8 and Appendix B; Lithwick, Xie and Wu 2012, ApJ 761:122; recalled, derived here):

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

from planetgen.physics import clearance_ledger, constants, oddities, tidal_heating

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
_PSI = (np.arange(_INTEGRAL_POINTS) + 0.5) / _INTEGRAL_POINTS * 2.0 * math.pi
"""ndarray: The midpoints of the full period the Laplace coefficients are integrated over."""


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


def _laplace(order, alpha, derivative=False, s=0.5):
    """`b_s^{(order)}(alpha)` (`s` is 1/2 for the resonant terms; the secular exchange of `secular_chain` reads 3/2) or its derivative in `alpha`, by the trapezoid rule on the full
    period (spectrally exact for this smooth integrand)."""
    psi = _PSI
    base = 1.0 - 2.0 * alpha * np.cos(psi) + alpha * alpha
    if derivative:
        return float(2.0 / _INTEGRAL_POINTS * np.sum(np.cos(order * psi) * -s * (2.0 * alpha - 2.0 * np.cos(psi)) / base ** (s + 1.0)))
    return float(2.0 / _INTEGRAL_POINTS * np.sum(np.cos(order * psi) / base ** s))


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


def _pair_eccentricities(pair, axes, masses, host_mass_kg):
    """`(e_inner, e_outer)` forced by one pair (zero for a second-order resonance), before the crossing limit; `axes` are the planets' semi-major axes, m."""
    if pair.order != 1:
        return 0.0, 0.0
    alpha = axes[pair.inner] / axes[pair.outer]
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
    chains = find_chains(periods)
    axes = [semi_major_axis_m(period, host_mass_kg, mass) for period, mass in zip(periods, masses_kg)] if chains else []
    for chain in chains:
        for pair in chain.pairs:
            inner, outer = _pair_eccentricities(pair, axes, masses_kg, host_mass_kg)
            squares[pair.inner] += inner ** 2
            squares[pair.outer] += outer ** 2
        for pair in chain.pairs:
            a_in, a_out = axes[pair.inner], axes[pair.outer]
            limit = CROSSING_SHARE * (a_out - a_in) / (a_out + a_in)
            for member in (pair.inner, pair.outer):
                squares[member] = min(squares[member], limit ** 2)
    return tuple(math.sqrt(s) for s in squares)


def forcing(host_mass_kg, period_days, planet_mass_kg, eccentricity):
    """The `tidal_heating.Forcing` of a planet on the orbit of the period with the eccentricity, to give `solid_power_w`, `molten_power_w` or `icy_flux_w_m2` as for a moon."""
    return tidal_heating.Forcing(host_mass_kg, semi_major_axis_m(period_days, host_mass_kg, planet_mass_kg), eccentricity)


def _clear_of_neighbours(bodies, forced, host_mass_kg):
    """
    `forced` (one eccentricity for each of `bodies`) cut to what still leaves each planet its spacing from the planets either side of it (the mutual Hill radii of
    `clearance_ledger.planet_exclusion`, which `validation.space_orbits` enforced when the orbits were placed): a resonance raises an eccentricity after the spacing was
    made, so an eccentricity that takes the apoapsis into the next planet's exclusion is held where the gap closes. A neighbour is taken at the higher of its own and its
    forced eccentricity, so two neighbours that both rise still clear each other.
    """
    ranked = sorted(range(len(bodies)), key=lambda i: bodies[i].distance)
    held = list(forced)

    def reach(i):
        return max(getattr(bodies[i], "eccentricity", None) or 0.0, held[i])

    for _ in range(len(ranked)):
        moved = False
        for place, i in enumerate(ranked):
            if not held[i] > 0.0:
                continue
            body, cap = bodies[i], held[i]
            if place > 0:
                inner = bodies[ranked[place - 1]]
                kappa = clearance_ledger.spacing_kappa(inner.mass, body.mass, host_mass_kg)
                cap = min(cap, 1.0 - kappa / 2.0 - inner.distance / body.distance * (1.0 + reach(ranked[place - 1]) + kappa / 2.0))
            if place + 1 < len(ranked):
                outer = bodies[ranked[place + 1]]
                kappa = clearance_ledger.spacing_kappa(body.mass, outer.mass, host_mass_kg)
                cap = min(cap, outer.distance / body.distance * (1.0 - reach(ranked[place + 1]) - kappa / 2.0) - 1.0 - kappa / 2.0)
            cap = max(cap, 0.0)
            if cap < held[i]:
                held[i], moved = cap, True
        if not moved:
            break
    return held


def give_forced_eccentricities(planets, host_mass_kg):
    """
    Finds the chains among the periods of `planets` (top-level planets of one star, any order) and keeps on each planet of a chain, as the quantity `forced_eccentricity`
    (`body.oddity_quantities`), the eccentricity its neighbours force on it; a planet in no chain has none. Reads and draws nothing else.

    Returns:
        int: How many planets were given one.
    """
    bodies = [p for p in planets if getattr(p, "body_type", None) in ("t", "g") and not getattr(p, "is_moon", False) and getattr(p, "period", None) and p.mass]
    for body in planets:
        quantities = dict(getattr(body, "oddity_quantities", None) or {})
        quantities.pop("forced_eccentricity", None)
        body.oddity_quantities = quantities
    if len(bodies) < 3 or not host_mass_kg > 0.0:
        return 0
    forced = _clear_of_neighbours(bodies, forced_eccentricities([p.period * constants.DAYS_PER_YEAR for p in bodies], [p.mass for p in bodies], host_mass_kg), host_mass_kg)
    given = 0
    for body, value in zip(bodies, forced):
        if value > 0.0:
            body.oddity_quantities = {**body.oddity_quantities, "forced_eccentricity": value}
            given += 1
    return given


@oddities.register("resonance_forced")
def resonance_forced(body, row):
    """
    The effect of oddity W5: the planet is held on the eccentricity its resonant neighbours force on it. The eccentricity rises to the forced one (held under what the star's Roche
    limit allows); a planet already at or above it is declined. The orbit is kept: the star's tide does not wear it down again (the resonance renews it), so `settle_orbit`
    leaves it. The caller makes the position, the moons' host eccentricity and the end stages again (`planets.finish_orbit_oddities`).
    """
    from planetgen.physics import orbit_elements

    forced = (getattr(body, "oddity_quantities", None) or {}).get("forced_eccentricity")
    if not forced:
        return False
    ceiling = orbit_elements.max_planet_eccentricity(body.distance, body.star, body)
    target = min(forced, ceiling)
    if getattr(body, "_resonance_e", None) == target:       # the orbit was already held on it (this step runs again): nothing more to do
        return None
    if target <= body.eccentricity:
        return False
    body.eccentricity = target
    body._orbit_birth = None
    body._resonance_e = target
    return None
