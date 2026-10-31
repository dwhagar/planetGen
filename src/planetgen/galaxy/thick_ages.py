# planetgen/galaxy/thick_ages.py

"""
The thick disc's ages in the old population (GEN.1025, after GEN.944 part 1)
============================================================================

The old population's window ended at 10 Gy, the thin disc's age, but the thick disc (which `galaxy.density.population_units` holds in `old`) formed
its stars from 13 to 8 Gy ago (Haywood et al. 2013, `tuning.THICK_DISC_AGE_RANGE_GY`), so no star of 10 to 13 Gy existed. This module builds the
`AgeLaw` of the old population as a mixture: the thin disc's own law (`disk_history`, 3 to 10 Gy) with the weight the thin disc has among the
old systems, and the thick disc's, flat over its window (stand-in TD2), with the thick disc's share of them (`thick_disk_units` over the old
population's units).

`mixed_old_law` is the pure mixture and `old_age_law(shape, thin)` the galaxy's (`disk_history.population_age_laws` calls it). Only the shape of the law matters (`AgeLaw`), so the counts of the
old population per sector, which come from the density model, do not change: the stars are spread over a longer range of ages.
"""

from planetgen import tuning
from planetgen.galaxy.age_law import AgeLaw

EDGE_GY = 0.01
"""float: The width, Gy, of the ramp that stands for a step in the law (an `AgeLaw` is piecewise linear through increasing ages)."""


def thick_share_of_old(units, thick_units):
    """The thick disc's share of the old population's systems, `units` being `density.population_units(shape)["old"]` and `thick_units`
    `density.thick_disk_units(shape)`."""
    if units <= 0.0:
        return 0.0
    return min(max(thick_units / units, 0.0), 1.0)


def _rate_at(law, age):
    """The law's relative density at `age` (held to its window), linear between its points."""
    low, high = law.window
    age = min(max(age, low), high)
    for (a0, a1), (r0, r1) in zip(zip(law.ages, law.ages[1:]), zip(law.rates, law.rates[1:])):
        if age <= a1:
            return r0 + (r1 - r0) * (age - a0) / (a1 - a0)
    return law.rates[-1]


def mixed_old_law(thin, thick_share, thick_window=None):
    """
    The `AgeLaw` of the old population: `thin` (the thin disc's law over its window) with weight `1 - thick_share` and the thick disc's flat law
    over `thick_window` (default `tuning.THICK_DISC_AGE_RANGE_GY`) with weight `thick_share`, both taken to unit area. With no thick disc
    (`thick_share` 0) it is `thin` itself.
    """
    if thick_share <= 0.0:
        return thin
    low, high = tuple(tuning.THICK_DISC_AGE_RANGE_GY if thick_window is None else thick_window)
    thick_density = thick_share / (high - low)
    thin_weight = 1.0 - thick_share
    thin_low, thin_high = thin.window
    ages = {float(age) for age in thin.ages} | {low - EDGE_GY, low, thin_high, thin_high + EDGE_GY, high}
    ages = sorted(age for age in ages if age >= thin_low and age <= max(high, thin_high + EDGE_GY))

    def density(age):
        in_thick = 0.0
        if low <= age <= high:
            in_thick = thick_density
        elif low - EDGE_GY < age < low:
            in_thick = thick_density * (age - (low - EDGE_GY)) / EDGE_GY
        thin_part = 0.0
        if age <= thin_high:
            thin_part = thin_weight * _rate_at(thin, age) / thin._cumulative[-1]
        elif age < thin_high + EDGE_GY:
            thin_part = thin_weight * _rate_at(thin, thin_high) / thin._cumulative[-1] * (1.0 - (age - thin_high) / EDGE_GY)
        return thin_part + in_thick

    rates = [density(age) for age in ages]
    top = max(rates)
    return AgeLaw(tuple(ages), tuple(round(rate / top, 6) for rate in rates))


def old_age_law(shape, thin):
    """The old population's `AgeLaw` for a galaxy `shape`, given the thin disc's law `thin`: `thin` mixed with the thick disc's by the thick
    disc's share of the old systems."""
    from planetgen.galaxy import density
    share = thick_share_of_old(density.population_units(shape)["old"], density.thick_disk_units(shape))
    return mixed_old_law(thin, share)
