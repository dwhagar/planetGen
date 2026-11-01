# planetgen/physics/remnant_mass.py

"""
The mass of a remnant, in one place (GEN.562)
=============================================

A white dwarf's mass comes from its progenitor's initial mass through the initial-final relation
(`tuning.WD_IFMR_PIECES`, [R0961], held inside `WD_MASS_RANGE_SOL`; `progenitor_mass_sol` is its inverse); a neutron star's and a stellar black hole's
mass is a Gaussian of the observed slow pulsars held to its range (`tuning.STELLAR_NS_MASS_GAUSSIANS`, `NEUTRON_STAR_MASS_RANGE_SOLAR`) and a black hole's uniform over its range (`BLACK_HOLE_MASS_RANGE_SOLAR`). The
star model, the multi-star member draw, the scatter's companions, the compact objects and the placement terms each held
their own copy of these lines; they all call this module. The numbers stay in `tuning`. A derived mass function
(U4/U5 of docs/design/architecture/unified-star-like-objects.md) changes this module and nothing else.
"""

from statistics import NormalDist

import numpy as np

from planetgen import tuning
from planetgen.util import draw

KINDS = ("neutron-star", "black-hole")
"""tuple: The kinds whose mass is drawn over a range."""


def _ifmr_nodes():
    """`(initial, final)` at the joins of `tuning.WD_IFMR_PIECES`: the relation is continuous and linear between them."""
    pieces = tuning.WD_IFMR_PIECES
    initial = [pieces[0][0]] + [piece[1] for piece in pieces]
    final = [pieces[0][3] + pieces[0][2] * initial[0]] + [piece[3] + piece[2] * piece[1] for piece in pieces]
    return initial, final


def white_dwarf_mass_sol(initial_mass_sol):
    """The white dwarf mass in solar masses of a progenitor of this initial mass (a number, or a NumPy array of them):
    the initial-final relation of `tuning.WD_IFMR_PIECES`, its end pieces carried beyond their ranges, held inside
    `tuning.WD_MASS_RANGE_SOL`."""
    initial, final = _ifmr_nodes()
    pieces = tuning.WD_IFMR_PIECES
    mass = np.asarray(initial_mass_sol, dtype=float)
    inside = np.interp(mass, initial, final)
    below = pieces[0][3] + pieces[0][2] * mass
    above = pieces[-1][3] + pieces[-1][2] * mass
    found = np.where(mass < initial[0], below, np.where(mass > initial[-1], above, inside))
    found = np.clip(found, *tuning.WD_MASS_RANGE_SOL)
    return float(found) if np.ndim(initial_mass_sol) == 0 else found


def progenitor_mass_sol(final_mass_sol):
    """The initial mass in solar masses a white dwarf of this mass was born with: the initial-final relation inverted (a
    number, or a NumPy array), its end pieces carried beyond their ranges."""
    initial, final = _ifmr_nodes()
    pieces = tuning.WD_IFMR_PIECES
    mass = np.asarray(final_mass_sol, dtype=float)
    inside = np.interp(mass, final, initial)
    below = (mass - pieces[0][3]) / pieces[0][2]
    above = (mass - pieces[-1][3]) / pieces[-1][2]
    found = np.where(mass < final[0], below, np.where(mass > final[-1], above, inside))
    return float(found) if np.ndim(final_mass_sol) == 0 else found


def mass_range_sol(kind):
    """`(low, high)` solar masses a neutron star or a stellar black hole draws from."""
    if kind == "neutron-star":
        return tuning.NEUTRON_STAR_MASS_RANGE_SOLAR
    if kind == "black-hole":
        return tuning.BLACK_HOLE_MASS_RANGE_SOLAR
    raise ValueError(f"{kind!r} is not a kind drawn over a mass range ({', '.join(KINDS)})")


def within(law_range, mass_range):
    """`law_range` truncated to `mass_range` (both `(low, high)` solar masses; `None` keeps the whole law)."""
    if mass_range is None:
        return law_range
    low, high = max(law_range[0], mass_range[0]), min(law_range[1], mass_range[1])
    if low > high:
        raise ValueError(f"mass range {mass_range!r} lies outside {law_range!r}")
    return low, high


def _neutron_star_normal():
    """`(NormalDist, (low, high))` of the neutron star mass law: the observed Gaussian of the slow pulsars, held to the range
    `tuning.NEUTRON_STAR_MASS_RANGE_SOLAR` (the lightest and heaviest neutron stars the generator makes)."""
    _population, mean, sigma = next(row for row in tuning.STELLAR_NS_MASS_GAUSSIANS if row[0] == "slow")
    return NormalDist(mean, sigma), tuning.NEUTRON_STAR_MASS_RANGE_SOLAR


def draw_mass_sol(kind, mass_range=None, rng=draw):
    """A neutron star's or stellar black hole's mass in solar masses, truncated to `mass_range` when given: a neutron star's
    from the Gaussian of the slow pulsars ([R3396]: 1.49 +- 0.19) held to the kind's range, a black hole's
    uniform over its range (the widest of the published ranges stands, `docs/design/galaxy/stellar-mass-function.md`).
    One draw from `rng`."""
    low, high = within(mass_range_sol(kind), mass_range)
    if kind == "neutron-star":
        normal, _range = _neutron_star_normal()
        cdf_low, cdf_high = normal.cdf(low), normal.cdf(high)
        return normal.inv_cdf(cdf_low + rng.random() * (cdf_high - cdf_low))
    return rng.uniform(low, high)


def share_above(kind, mass_sol):
    """The share of the kind's mass law at or above `mass_sol` (1 below its range, 0 above it)."""
    low, high = mass_range_sol(kind)
    if mass_sol <= low:
        return 1.0
    if mass_sol >= high:
        return 0.0
    if kind == "neutron-star":
        normal, _range = _neutron_star_normal()
        return (normal.cdf(high) - normal.cdf(mass_sol)) / (normal.cdf(high) - normal.cdf(low))
    return (high - mass_sol) / (high - low)


def moments(kind):
    """`(mean, mean of the square)` of the kind's mass draw, for the placement terms."""
    low, high = mass_range_sol(kind)
    if kind == "neutron-star":
        normal, _range = _neutron_star_normal()
        grid = np.linspace(low, high, 20001)
        weight = np.array([normal.pdf(float(m)) for m in grid])
        weight /= np.trapezoid(weight, grid)
        return float(np.trapezoid(weight * grid, grid)), float(np.trapezoid(weight * grid * grid, grid))
    return 0.5 * (low + high), (low * low + low * high + high * high) / 3.0
