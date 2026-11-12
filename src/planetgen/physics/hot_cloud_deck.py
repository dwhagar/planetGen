# planetgen/physics/hot_cloud_deck.py

"""
Where the clouds of a hot atmosphere sit, and the most each trace gas can hold (GEN.994, built ahead by Lane 07 for Lane 14, stacked on PR 5976)
====================================================================================================================================================

`hot_condensates` gives the temperature at which a cloud species starts to condense at a pressure (`cloud_temperature_k`) and the most molecules of a trace gas the host
leaves room for (`gas_ceiling`). The hot-giant regime wiring needs the two answers it will store: the *deck*, where each cloud species meets a temperature-pressure
profile, and the *ceiling* of each trace gas. This module is that step and nothing more: pure functions of the host and the profile, no draw, no call into
`mixture_giant`, no stored column yet.

**The deck.** A profile is a function `T(P)` in K of the pressure in bar that rises with depth. A species condenses where the profile is cooler than its condensation
curve, so the cloud base is the deepest pressure where the two meet and the cloud lies above it. `cloud_base` finds that meeting on a log-spaced grid of pressures
(`GRID_BAR`) and refines it by bisection. It returns None for a species the profile never cools enough to condense, and flags `reaches_bottom` when the profile is cooler than
the curve down to the deepest grid pressure (the cloud has no base inside the grid).

**Column plan (option (a), Lane 14).** One stored quantity per species, not one row per layer: for each of `hot_condensates.CLOUD_SPECIES` the base pressure in bar (NULL for
no cloud), and for each of `hot_condensates.HOT_TRACE_GASES` the ceiling in molecules per hydrogen atom. `column_names` gives the column names in a fixed order
(`cloud_base_bar_<species lowercased>`, `trace_ceiling_<formula lowercased>`) so the migration, `db/column_units.py` (unit "bar" and "mol/mol H") and the describe page
are written from one list. The migration, the unit registry entry and the epoch belong to the item that wires this (GEN.994); none is claimed here.
"""

import math
from dataclasses import dataclass

from planetgen.physics import hot_condensates

GRID_BAR = tuple(10.0 ** (exponent / 4.0) for exponent in range(-24, 13))
"""tuple: The pressures, bar, the profile is searched on: 1e-6 to 1e3, four to the decade."""

TOLERANCE = 1.0e-3
"""float: The bisection stops when the bracket is this narrow in log10 of the pressure."""


@dataclass(frozen=True)
class CloudBase:
    """One species' cloud: the deepest pressure where the profile meets the condensation curve, and the temperature there (`reaches_bottom`: no meeting inside the grid)."""

    species: str
    base_bar: float
    base_k: float
    reaches_bottom: bool = False


def _gap(species, host, profile, pressure_bar):
    """Profile minus condensation temperature, K: positive where the gas is too hot to condense."""
    return profile(pressure_bar) - hot_condensates.cloud_temperature_k(species, pressure_bar, host)


def cloud_base(species, host, profile):
    """The `CloudBase` of `species` in `host`'s gas along `profile` (a function of pressure in bar giving K), or None when the profile is nowhere cool enough."""
    gaps = [_gap(species, host, profile, pressure) for pressure in GRID_BAR]
    if gaps[-1] < 0.0:
        return CloudBase(species, GRID_BAR[-1], profile(GRID_BAR[-1]), reaches_bottom=True)
    for index in range(len(GRID_BAR) - 1, 0, -1):
        if gaps[index] >= 0.0 > gaps[index - 1]:
            low, high = math.log10(GRID_BAR[index - 1]), math.log10(GRID_BAR[index])
            while high - low > TOLERANCE:
                mid = 0.5 * (low + high)
                if _gap(species, host, profile, 10.0 ** mid) < 0.0:
                    low = mid
                else:
                    high = mid
            base = 10.0 ** (0.5 * (low + high))
            return CloudBase(species, base, profile(base))
    return None


def deck(host, profile, species=None):
    """The clouds that form along `profile`, deepest base first, as a tuple of `CloudBase`; `species` defaults to every species of `hot_condensates.CLOUD_SPECIES`."""
    found = (cloud_base(name, host, profile) for name in (species or hot_condensates.CLOUD_SPECIES))
    return tuple(sorted((base for base in found if base is not None), key=lambda base: -base.base_bar))


def trace_gas_ceilings(host):
    """The most molecules per hydrogen atom each trace gas of `hot_condensates.HOT_TRACE_GASES` can hold in `host`'s gas, by formula."""
    return {formula: hot_condensates.gas_ceiling(formula, host) for formula in hot_condensates.HOT_TRACE_GASES}


def column_names():
    """The stored columns of the plan, in order: a base pressure per cloud species, then a ceiling per trace gas."""
    return (tuple(f"cloud_base_bar_{name.lower()}" for name in hot_condensates.CLOUD_SPECIES)
            + tuple(f"trace_ceiling_{formula.lower()}" for formula in hot_condensates.HOT_TRACE_GASES))


def stored_quantities(host, profile):
    """`column_names` mapped to the numbers to store: the base pressure of each cloud (None where none forms) and the ceiling of each trace gas."""
    bases = {base.species: base.base_bar for base in deck(host, profile)}
    values = [bases.get(name) for name in hot_condensates.CLOUD_SPECIES]
    values += [ceiling for ceiling in trace_gas_ceilings(host).values()]
    return dict(zip(column_names(), values))
