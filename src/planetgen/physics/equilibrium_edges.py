# planetgen/physics/equilibrium_edges.py

"""
The temperatures at which a hot atmosphere's chemistry changes regime (GEN.518, GEN.672, docs/design/composition/equilibrium-chemistry.md)
===========================================================================================================================================

Three edges of the hot-Jupiter regimes (rows h1 to h3 of the recipe tables), found by bisection on `gas_equilibrium.solve`, the one
equilibrium solver of the atmosphere chain:

* `ch4_co_temperature`: carbon monoxide and methane equally abundant (above it CO holds the carbon);
* `h2o_survival_temperature`: the water has fallen to a share of the oxygen left once carbon monoxide has taken its carbon;
* `h2_dissociation_temperature`: half of the hydrogen nuclei in H2 and H are atomic.

The budget is that of the independent computation the edges were checked against (tunable-gaps-sources/group4-giants-bd-clouds.md
section 2.4): hydrogen, helium, carbon, nitrogen, oxygen, sulfur and argon at the solar ratios of `mixture_giant.SOLAR_PER_HYDROGEN`,
the metals scaled by the metallicity and the carbon set to `c_to_o` times the oxygen (`reference_budget`).

`regime_edge` is the one the recipe rows h1 to h3 read (GEN.994): the edge at a pressure for the reference gas of solar metals and a
C/O of 0.55, cached by pressure.

Use:

    budget = reference_budget(metallicity=1.0, c_to_o=0.55)
    ch4_co_temperature(0.1, budget)            # K at 0.1 bar
"""

import functools
import math

from planetgen.physics import gas_equilibrium
from planetgen.physics.gas_equilibrium import EquilibriumError
from planetgen.physics.mixture_giant import SOLAR_PER_HYDROGEN

BUDGET_METALS = ("C", "N", "O", "S", "Ar")
"""tuple: The elements scaled by the metallicity in `reference_budget`."""

REFERENCE_C_TO_O = 0.55
"""float: The C/O of the reference gas of the recipe rows h1 to h3 (the solar value, Asplund et al. 2009). [Stand-in HE1, tier T1, tweak L0.]"""

WATER_EDGE_SHARE = 0.9
"""float: The share of the possible water that survives at the edge row h3 and the hot side of row h2 call "water destroyed": a tenth is gone
(a definition of the regime, GEN.518). [Stand-in HE2, tier "from needed output", tweak L0.]"""

REGIMES = ("methane_carbon_monoxide", "water_destroyed")
"""tuple: The regime edges `regime_edge` gives, the lower and the upper edge of the hot-Jupiter rows."""

TOLERANCE_K = 0.5
"""float: The bisection stops when the bracket is this narrow."""


def reference_budget(metallicity=1.0, c_to_o=None):
    """
    The atoms of each element per 1 hydrogen atom for a gas of `metallicity` times the reference metals, helium at the reference
    ratio. `c_to_o` sets the carbon to that multiple of the oxygen (None keeps the reference carbon, scaled with the other metals).

    Raises:
        EquilibriumError: The metallicity or the C/O is not positive.
    """
    if not metallicity > 0:
        raise EquilibriumError(f"the metallicity {metallicity!r} is not positive")
    budget = {"H": 1.0, "He": SOLAR_PER_HYDROGEN["He"]}
    for symbol in BUDGET_METALS:
        budget[symbol] = metallicity * SOLAR_PER_HYDROGEN[symbol]
    if c_to_o is not None:
        if not c_to_o > 0:
            raise EquilibriumError(f"the C/O {c_to_o!r} is not positive")
        budget["C"] = c_to_o * budget["O"]
    return budget


def _bisect(function, low, high, budget, pressure_bar):
    """The temperature in [low, high] where `function(equilibrium)` crosses zero (negative at `low`, positive at `high`)."""
    f_low = function(gas_equilibrium.solve(budget, low, pressure_bar))
    f_high = function(gas_equilibrium.solve(budget, high, pressure_bar))
    if not f_low < 0 < f_high:
        raise EquilibriumError(f"no crossing between {low} and {high} K at {pressure_bar} bar")
    while high - low > TOLERANCE_K:
        mid = 0.5 * (low + high)
        if function(gas_equilibrium.solve(budget, mid, pressure_bar)) < 0:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def ch4_co_temperature(pressure_bar, budget):
    """The temperature, K, where carbon monoxide and methane are equally abundant (above it CO holds the carbon)."""
    return _bisect(lambda e: math.log(e.moles.get("CO", 1e-300) / e.moles.get("CH4", 1e-300)), 300.0, 4000.0, budget, pressure_bar)


def h2o_survival_temperature(pressure_bar, budget, share=0.5):
    """
    The temperature, K, where the water has fallen to `share` of the oxygen left once carbon monoxide has taken its carbon
    (O - C atoms), i.e. where thermal dissociation has destroyed `1 - share` of the water that could exist.

    Raises:
        EquilibriumError: The carbon is at least the oxygen (no oxygen is left for water).
    """
    available = budget["O"] - budget["C"]
    if not available > 0:
        raise EquilibriumError("no oxygen is left for water once carbon monoxide has taken its share (C/O of 1 or more)")
    return _bisect(lambda e: share - e.moles.get("H2O", 0.0) / available, 1500.0, 6000.0, budget, pressure_bar)


def h2_dissociation_temperature(pressure_bar, budget):
    """The temperature, K, where half of the hydrogen nuclei in H2 and H are atomic (`n(H) = 2 n(H2)`)."""
    return _bisect(lambda e: e.moles.get("H", 0.0) - 2.0 * e.moles.get("H2", 0.0), 1000.0, 6000.0, budget, pressure_bar)


@functools.lru_cache(maxsize=256)
def regime_edge(regime, pressure_bar):
    """
    The temperature, K, of the hot-Jupiter regime edge `regime` (one of `REGIMES`) at `pressure_bar`, for a gas of solar metals at
    `REFERENCE_C_TO_O`; the edges the recipe rows h1 to h3 are half open at.

    Raises:
        ValueError: `regime` is not one of `REGIMES`.
    """
    budget = reference_budget(1.0, REFERENCE_C_TO_O)
    if regime == "methane_carbon_monoxide":
        return ch4_co_temperature(pressure_bar, budget)
    if regime == "water_destroyed":
        return h2o_survival_temperature(pressure_bar, budget, WATER_EDGE_SHARE)
    raise ValueError(f"not a regime edge: {regime!r}")
