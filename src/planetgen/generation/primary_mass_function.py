# planetgen/generation/primary_mass_function.py

"""
The Primary Mass Function (GEN.928)
===================================

The mass function of the stars that are the primary of a system, derived and not set. The single-star function
(`physics/imf.py`, `imf.SINGLE`) counts every star born; a system's companions are stars too, and each is lighter than its
primary. So at each mass the stars born are the primaries of that mass plus the companions of the heavier primaries that fall
there:

    xi(m) = P(m) + sum over primaries M >= m of P(M) * n_c(M) * (share of M's companions that fall at m)

with `n_c(M)` the companions a primary of mass M has (the mean of the `star-count` table's count less one, by the class of the
star, `hosts.STAR_COUNT_PROBABILITIES`) and the companions' masses as `multistar` draws them: each uniform between a tenth of
the primary's mass (or the hydrogen-burning limit, where that is higher) and the primary's mass. Companions are never heavier
than the primary, so the sum runs over primaries at or above `m` and the equation is solved from the top mass down, no
iteration. A system's primary is then drawn from `P`, and the stars of its tree from the tree's own draws, and the stars of all
the systems together follow `xi` by construction: the single-star function stays the one set curve and the systems' function
(peaking higher than the single stars', as the published ones do) comes out of the multiplicity rows.

`imf.PRIMARY` is this function, as the same continuous piecewise power law `imf` keeps `xi` in (`segments`), so every consumer
of the function takes either by naming it. Systems per star born is the ratio of the two numbers (`systems_per_star`) and
systems per solar mass the ratio of `P`'s number to `xi`'s mass (`systems_per_solar_mass`, GEN.929).

What the derivation leaves out, and the reason: the companions past the seventh star of a system (`multistar.TAIL_RATIOS`, a
chance of about 1e-5 of a system); the brown dwarf companion (below the hydrogen-burning limit, not in `xi`); and the mergers and
stripped stars of `multistar`, which move a star between mass bins after its birth.
"""

import functools
import math

import numpy as np

from planetgen import tuning
from planetgen.generation import mass_ratio, multistar
from planetgen.util import detmath
from planetgen.generation.phenomena import hosts
from planetgen.physics import imf, mass_boundaries

GRID_STEP_DEX = 0.002
"""float: The width, in dex of mass, of the cells the equation is solved in. The answer is read off at the nodes
(`NODE_SPACING_DEX`) and does not move at this resolution (a test compares it with a grid twice as coarse)."""

NODE_SPACING_DEX = 0.05
"""float: The spacing of the nodes the primary function is cut into, as for `tuning.IMF_NODE_SPACING_DEX`."""


def companions_per_primary(mass_sol):
    """`n_c(M)`: the companions a primary of this birth mass has, the mean of the `star-count` table's count less one for its class."""
    weights = hosts.star_count_weights(multistar.companion_class(mass_sol), 1.0)
    return sum(weight * stars for stars, weight in enumerate(weights)) / sum(weights)


def _edges():
    low, high = imf.lowest_mass_sol(), imf.highest_mass_sol()
    cells = int(math.ceil(math.log10(high / low) / GRID_STEP_DEX))
    return low * (high / low) ** (np.arange(cells + 1) / cells)


@functools.lru_cache(maxsize=None)
def _solve():
    """`(edges, number)`: the cell edges (solar masses) and the number of primaries born in each, on the scale of `xi`."""
    edges = _edges()
    cells = len(edges) - 1
    centre = np.sqrt(edges[:-1] * edges[1:])
    born = np.array([imf.integral(0.0, a, b) for a, b in zip(edges[:-1], edges[1:])])
    n_c = np.array([companions_per_primary(m) for m in centre])
    ratios = mass_ratio.Ratios(centre)
    primaries = np.zeros(cells)
    for i in range(cells - 1, -1, -1):
        # A companion is no lighter than `tuning.BINARY_MASS_RATIO_RANGE`'s low end of its primary, so the primaries that put one
        # in this cell are within that factor above it.
        end = min(cells, int(np.searchsorted(centre, edges[i + 1] / tuning.BINARY_MASS_RATIO_RANGE[0])) + 1)
        rows = slice(i, end)
        share = (ratios.fraction_below(edges[i + 1] / centre[rows], rows) - ratios.fraction_below(edges[i] / centre[rows], rows))
        own = n_c[i] * share[0]
        heavier = n_c[i + 1:end] * share[1:]
        from_heavier = float(detmath.matmul(primaries[i + 1:end], heavier[:, None])[0]) if heavier.size else 0.0
        primaries[i] = max(0.0, (born[i] - from_heavier) / (1.0 + own))
    return edges, primaries


def _nodes():
    """The masses the primary function is cut at: every `NODE_SPACING_DEX`, and every mass where a rule changes."""
    low, high = imf.lowest_mass_sol(), imf.highest_mass_sol()
    cells = int(math.ceil(math.log10(high / low) / NODE_SPACING_DEX))
    nodes = {low * (high / low) ** (i / cells) for i in range(cells + 1)}
    nodes.update(edge for _letter, edge in mass_boundaries.SPECTRAL_CLASS_LOWER_MASS_SOL if low < edge < high)
    nodes.update(edge for _name, edge in mass_boundaries.B_SUBCLASS_LOWER_MASS_SOL)
    nodes.update((tuning.IMF_BELL_TOP_SOL, mass_boundaries.sol("core_collapse_limit")))
    nodes = sorted(node for node in nodes if low <= node <= high)
    kept = [nodes[0]]
    for node in nodes[1:]:
        if node > kept[-1] * (1.0 + 1.0e-9):
            kept.append(node)
    return kept


def density_at(masses):
    """`P(M)` per solar mass at each of `masses` (solar masses), on the scale of `xi`, read off the solved cells (log-log)."""
    edges, primaries = _solve()
    centre = np.sqrt(edges[:-1] * edges[1:])
    per_mass = primaries / (edges[1:] - edges[:-1])
    where, log_density, at = np.log(masses), np.log(np.maximum(per_mass, 1.0e-300)), np.log(centre)
    found = np.interp(where, at, log_density)
    # Past the first and last cell centre the line of the nearest two continues (not flat: the density falls on).
    found = np.where(where > at[-1], log_density[-1] + (where - at[-1]) * (log_density[-1] - log_density[-2]) / (at[-1] - at[-2]), found)
    found = np.where(where < at[0], log_density[0] + (where - at[0]) * (log_density[1] - log_density[0]) / (at[1] - at[0]), found)
    return np.exp(found)


@functools.lru_cache(maxsize=None)
def segments():
    """`((low, high, alpha, k), ...)` of `P`, the continuous piecewise power law through its value at each node."""
    nodes = _nodes()
    values = density_at(np.array(nodes))
    if not np.all(values > 0.0):
        raise ValueError("the primary mass function is not positive: the companions the rows give outnumber the stars born")
    out = []
    for a, b, fa, fb in zip(nodes, nodes[1:], values, values[1:]):
        alpha = -math.log(fb / fa) / math.log(b / a)
        out.append((a, b, alpha, float(fa) * a ** alpha))
    return tuple(out)


def systems_per_star():
    """The systems born for every star born: the number of primaries over the number of stars."""
    return imf.number(imf.lowest_mass_sol(), imf.highest_mass_sol(), imf.PRIMARY) / imf.number(
        imf.lowest_mass_sol(), imf.highest_mass_sol())


def systems_per_solar_mass():
    """The systems born for every solar mass of stars born (GEN.929)."""
    low, high = imf.lowest_mass_sol(), imf.highest_mass_sol()
    return imf.number(low, high, imf.PRIMARY) / imf.integral(1.0, low, high)
