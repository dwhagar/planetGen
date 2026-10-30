# planetgen/physics/ice_species.py

"""
The non-water ices of a body, one amount each (GEN.834, docs/design/composition/bulk-composition.md)
===================================================================================================

A body's `ice_mass_fraction` is every ice but water as one share of its whole mass. This module splits it by species into the columns
`ice_<formula>_mass_fraction` (`co2`, `ch4`, `nh3`, `n2`, `so2`, `co`, `hcn`), and the rest, the ices the surface rows have no entry for
(ethane, methanol, formaldehyde, hydrogen sulfide), into `ice_other_mass_fraction`. The eight add up to `ice_mass_fraction`, so the object page,
the sheet and the ground draw read the same amounts the bulk draw made.

1. `shares_of_percent(percent)` turns the drawn comet volatiles (percent of the water, by formula) into shares of the ice.
2. `median_shares()` is the stand-in for a body whose ice was not drawn species by species (a moon's makeup, a giant's ice mantle, a small
   body): the comet volatiles at the medians of the population (registry row IS1).
3. `attributes(ice, shares)` gives the eight column values, each `ice * share`; `shares_of(body)` reads the stored amounts back as shares of
   the ice, or `None` for a body saved before (the caller then uses `median_shares()`).
"""

import math

from planetgen.physics import comet_population

SPECIES = ("CO2", "CH4", "NH3", "N2", "SO2", "CO", "HCN")
"""tuple: The formulas of the ices that have a column each."""

OTHER = "other"
"""str: The key of every ice that has no column of its own."""

FIELD_OF = {**{formula: f"ice_{formula.lower()}_mass_fraction" for formula in SPECIES}, OTHER: f"ice_{OTHER}_mass_fraction"}
"""dict: The body attribute (and column) that holds each species' share of the whole body's mass."""

BODY_FIELDS = tuple(FIELD_OF.values())
"""tuple: Every attribute (and column) a body holds for its non-water ices."""

COMET_CLASS = "jupiter_family"
"""str: The comet class whose volatiles a body's ice follows (`bulk_draw.RANGES["other_ice_comet_class"]`)."""


def shares_of_percent(percent):
    """`{formula or OTHER: share of the ice}` from comet volatiles (percent of H2O by formula), adding up to 1; empty if there are none."""
    amounts = {}
    for formula, value in percent.items():
        if value > 0.0:
            key = formula if formula in SPECIES else OTHER
            amounts[key] = amounts.get(key, 0.0) + value
    total = math.fsum(amounts.values())
    return {key: value / total for key, value in amounts.items()} if total > 0.0 else {}


def median_shares(period_class=COMET_CLASS):
    """The stand-in shares: the comet volatiles of `period_class` at the medians of the population, with no draw."""
    return shares_of_percent(comet_population.median_volatiles(period_class))


def attributes(ice, shares=None):
    """The eight column values of a body with the ice mass fraction `ice`, split by `shares` (`median_shares()` if none are given)."""
    if not ice or ice <= 0.0:
        return {field: 0.0 for field in BODY_FIELDS}
    shares = median_shares() if not shares else shares
    return {FIELD_OF[key]: (ice * shares.get(key, 0.0)) for key in FIELD_OF}


def shares_of(body):
    """The stored ice of `body` as `{formula or OTHER: share of its ice}`, or `None` when it holds no species amounts or has no ice."""
    amounts = {key: getattr(body, field, None) for key, field in FIELD_OF.items()}
    if any(value is None for value in amounts.values()):
        return None
    total = math.fsum(amounts.values())
    return {key: value / total for key, value in amounts.items() if value > 0.0} if total > 0.0 else None
