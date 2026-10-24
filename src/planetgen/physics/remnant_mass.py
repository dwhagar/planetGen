# planetgen/physics/remnant_mass.py

"""
The mass of a remnant, in one place (GEN.562)
=============================================

A white dwarf's mass comes from its progenitor's initial mass through the initial-final relation
(`tuning.WD_IFMR_SLOPE`, `WD_IFMR_INTERCEPT_SOL`, held inside `WD_MASS_RANGE_SOL`); a neutron star's and a stellar black hole's
mass is drawn uniformly over their ranges (`tuning.NEUTRON_STAR_MASS_RANGE_SOLAR`, `BLACK_HOLE_MASS_RANGE_SOLAR`). The
star model, the multi-star member draw, the scatter's companions, the compact objects and the placement terms each held
their own copy of these lines; they all call this module. The numbers stay in `tuning`. A derived mass function
(U4/U5 of docs/design/architecture/unified-star-like-objects.md) changes this module and nothing else.
"""

import numpy as np

from planetgen import tuning
from planetgen.util import draw

KINDS = ("neutron-star", "black-hole")
"""tuple: The kinds whose mass is drawn over a range."""


def white_dwarf_mass_sol(initial_mass_sol):
    """The white dwarf mass in solar masses of a progenitor of this initial mass (a number, or a NumPy array of them):
    the initial-final relation, held inside `tuning.WD_MASS_RANGE_SOL`."""
    mass = tuning.WD_IFMR_SLOPE * initial_mass_sol + tuning.WD_IFMR_INTERCEPT_SOL
    low, high = tuning.WD_MASS_RANGE_SOL
    if isinstance(initial_mass_sol, np.ndarray):
        return np.clip(mass, low, high)
    return min(max(mass, low), high)


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


def draw_mass_sol(kind, mass_range=None, rng=draw):
    """A neutron star's or stellar black hole's mass in solar masses: one uniform draw over the kind's range, truncated to
    `mass_range` when given."""
    return rng.uniform(*within(mass_range_sol(kind), mass_range))


def uniform_moments(kind):
    """`(mean, mean of the square)` of the kind's uniform mass draw, for the placement terms."""
    low, high = mass_range_sol(kind)
    return 0.5 * (low + high), (low * low + low * high + high * high) / 3.0
