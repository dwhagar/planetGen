# stellarObjects/stellarPopulation.py

"""
Bright and Dim Star Sampling
============================

The star population model (`stellarEvolution`) split at a luminosity
threshold, so a whole galaxy's bright stars can be drawn and placed at
plan time and each sector later filled with only the dimmer rest (see the
project's bright-star pre-placement plan):

- `bright_star_fraction(min_luminosity_sol, population)`: the share of a
  population's living stars at or above the threshold.
- `sample_bright_stars(n, min_luminosity_sol, population, rng)`: `n` stars
  drawn from the model conditional on being at least that bright, each as
  `stellarEvolution.star_params` (what `Star.from_params` rebuilds).
- `sample_dim_star(max_luminosity_sol, population, rng)`: one star from the
  model conditional on being dimmer than the threshold.

`population` is one of `stellarEvolution.population_age_range_gy`'s names
("young", "intermediate", "old", "bulge") or None for the whole disk.

How the bright draw works
-------------------------
Bright stars are rare (well under 1 in 1,000), so drawing whole-model
stars and keeping the bright ones would take thousands of draws each.
Instead, for a star of initial mass m, every phase it can be bright in is
an age window (its main sequence if its main-sequence luminosity is
already enough, the late part of its subgiant phase, its giant or
supergiant phase) times the chance its (random) luminosity in that phase
clears the threshold. Their overlap with the population's age range is
the star's "bright measure" (`_bright_windows`). The mass range is cut
into log-spaced cells (`BRIGHT_STAR_MASS_GRID_CELLS`, split at every mass
where a phase's rule changes), each weighted by its IMF share times an
upper bound on the measure inside it. A draw picks a cell by weight, a
mass in it from the IMF, keeps it with chance measure / bound (so the
result is exact, not a grid approximation), then picks a phase by its
share of the measure, an age inside that phase's bright window and a
luminosity from the part of the phase's range above the threshold.

White dwarfs are never pre-placed (`WD_LUMINOSITY_RANGE_SOL` tops out at
100 Lsun), so the threshold must be at least that range's top, and a
sector's dim draw keeps every white dwarf whatever the cap.
"""

import bisect
import functools
import math
import random

from . import program_constants
from .stellarEvolution import (
    _power_law_integral, _sample_power_law, evolve_star, main_sequence_lifetime_gy,
    main_sequence_luminosity_sol, population_age_range_gy, sample_living_star, star_params,
)

POPULATIONS = tuple(program_constants.STELLAR_POPULATION_AGE_RANGES_GY)
"""The named stellar populations, in `population_densities` order."""


def _overlap(low, high, window):
    return max(0.0, min(high, window[1]) - max(low, window[0]))


def _giant_luminosity_range(mass_sol):
    if mass_sol >= program_constants.BRIGHT_GIANT_MIN_MASS_SOL:
        return program_constants.BRIGHT_GIANT_LUMINOSITY_RANGE_SOL
    return program_constants.GIANT_LUMINOSITY_RANGE_SOL


def _giant_bright_chance(mass_sol, min_luminosity_sol):
    """The chance a giant's log-uniform luminosity is >= the threshold."""
    low, high = _giant_luminosity_range(mass_sol)
    if min_luminosity_sol <= low:
        return 1.0
    if min_luminosity_sol >= high:
        return 0.0
    return math.log(high / min_luminosity_sol) / math.log(high / low)


def _supergiant_bright_chance(ms_luminosity_sol, min_luminosity_sol):
    """The chance a supergiant's uniform luminosity growth takes it >= the threshold."""
    low, high = program_constants.SUPERGIANT_LUMINOSITY_GROWTH_RANGE
    need = min_luminosity_sol / ms_luminosity_sol
    if need <= low:
        return 1.0
    if need >= high:
        return 0.0
    return (high - need) / (high - low)


def _bright_windows(mass_sol, min_luminosity_sol):
    """
    Every age window in which a star of initial mass `mass_sol` can be at
    least `min_luminosity_sol` bright, as `(start_gy, end_gy, chance)`:
    within the window its age puts it in a phase whose luminosity clears
    the threshold with probability `chance`. Mirrors `evolve_star`.
    """
    pc = program_constants
    t_ms = main_sequence_lifetime_gy(mass_sol)
    t_sub = t_ms * pc.SUBGIANT_PHASE_END_MS_FRACTION
    t_giant = t_ms * pc.GIANT_PHASE_END_MS_FRACTION
    l_ms = main_sequence_luminosity_sol(mass_sol)
    windows = []
    if l_ms >= min_luminosity_sol:
        windows.append((0.0, t_ms, 1.0))
    if mass_sol >= pc.SUPERGIANT_MIN_MASS_SOL:
        chance = _supergiant_bright_chance(l_ms, min_luminosity_sol)
        if chance > 0:
            windows.append((t_ms, t_giant, chance))
        return windows
    # A subgiant brightens linearly from l_ms to SUBGIANT_MAX_BRIGHTENING * l_ms.
    progress = (min_luminosity_sol / l_ms - 1) / (pc.SUBGIANT_MAX_BRIGHTENING - 1)
    if progress < 1:
        windows.append((t_ms + (t_sub - t_ms) * max(progress, 0.0), t_sub, 1.0))
    chance = _giant_bright_chance(mass_sol, min_luminosity_sol)
    if chance > 0:
        windows.append((t_sub, t_giant, chance))
    return windows


def _bright_measure(mass_sol, min_luminosity_sol, window):
    """The chance a star of this mass, of an age uniform in `window`, is
    at least `min_luminosity_sol` bright (times the window's length)."""
    return sum(_overlap(start, end, window) * chance
               for start, end, chance in _bright_windows(mass_sol, min_luminosity_sol))


def _living_measure(mass_sol, window):
    """The chance a star of this mass, of an age uniform in `window`, has
    not collapsed (times the window's length); see `evolve_star`."""
    span = window[1] - window[0]
    if mass_sol < program_constants.SUPERGIANT_MIN_MASS_SOL:
        return span
    dead_from = main_sequence_lifetime_gy(mass_sol) * program_constants.GIANT_PHASE_END_MS_FRACTION
    return span - _overlap(dead_from, math.inf, window)


def _bright_measure_bound(low_mass, high_mass, min_luminosity_sol, window):
    """
    An upper bound on `_bright_measure` for every mass in `[low_mass,
    high_mass]`, a cell that no phase-rule mass (`_mass_grid`'s breaks)
    splits: each phase's window, stretched over every mass in the cell
    (lifetimes fall and luminosities rise with mass), at its best chance.
    """
    pc = program_constants
    t_long = main_sequence_lifetime_gy(low_mass)
    t_short = main_sequence_lifetime_gy(high_mass)
    l_max = main_sequence_luminosity_sol(high_mass)
    bound = 0.0
    if l_max >= min_luminosity_sol:
        bound += _overlap(0.0, t_long, window)
    if low_mass >= pc.SUPERGIANT_MIN_MASS_SOL:
        chance = _supergiant_bright_chance(l_max, min_luminosity_sol)
        return bound + chance * _overlap(t_short, t_long * pc.GIANT_PHASE_END_MS_FRACTION, window)
    if l_max * pc.SUBGIANT_MAX_BRIGHTENING > min_luminosity_sol:
        bound += _overlap(t_short, t_long * pc.SUBGIANT_PHASE_END_MS_FRACTION, window)
    chance = _giant_bright_chance(low_mass, min_luminosity_sol)
    bound += chance * _overlap(t_short * pc.SUBGIANT_PHASE_END_MS_FRACTION,
                               t_long * pc.GIANT_PHASE_END_MS_FRACTION, window)
    return bound


def _mass_grid():
    """Log-spaced cell edges over the IMF's range, with every mass where
    the IMF slope or a phase rule changes added as an edge, and each
    cell's IMF slope and continuity factor."""
    pc = program_constants
    low, high = pc.IMF_BREAKS_SOL[0], pc.IMF_BREAKS_SOL[-1]
    cells = pc.BRIGHT_STAR_MASS_GRID_CELLS
    edges = {low * (high / low) ** (i / cells) for i in range(cells + 1)}
    edges.update(pc.IMF_BREAKS_SOL)
    edges.update((pc.BRIGHT_GIANT_MIN_MASS_SOL, pc.SUPERGIANT_MIN_MASS_SOL))
    edges = sorted(edge for edge in edges if low <= edge <= high)

    norms = [1.0]
    for i in range(1, len(pc.IMF_SLOPES)):
        norms.append(norms[-1] * pc.IMF_BREAKS_SOL[i] ** (pc.IMF_SLOPES[i] - pc.IMF_SLOPES[i - 1]))
    grid = []
    for a, b in zip(edges, edges[1:]):
        segment = bisect.bisect_right(pc.IMF_BREAKS_SOL, a) - 1
        grid.append((a, b, pc.IMF_SLOPES[segment], norms[segment]))
    return grid


def _check_threshold(min_luminosity_sol):
    if not (isinstance(min_luminosity_sol, (int, float)) and math.isfinite(min_luminosity_sol)):
        raise ValueError(f"luminosity threshold must be a finite number, got {min_luminosity_sol!r}")
    if min_luminosity_sol < program_constants.WD_LUMINOSITY_RANGE_SOL[1]:
        raise ValueError(f"luminosity threshold {min_luminosity_sol} Lsun must be at least the brightest white dwarf "
                         f"({program_constants.WD_LUMINOSITY_RANGE_SOL[1]} Lsun)")


# Quadrature points per cell for `bright_star_fraction`: the cell's IMF
# quantiles at (j + 1/2) / n.
_FRACTION_POINTS_PER_CELL = 4


@functools.lru_cache(maxsize=None)
def _bright_table(min_luminosity_sol, population):
    """
    The cached per-(threshold, population) sampling table: cell edges,
    slopes and continuity factors, each cell's measure bound, the running
    total of cell weights (IMF share times bound), and the bright fraction.
    """
    _check_threshold(min_luminosity_sol)
    window = population_age_range_gy(population)
    cells, cumulative = [], []
    running = bright = living = 0.0
    for a, b, alpha, norm in _mass_grid():
        count = norm * _power_law_integral(a, b, alpha)
        # E[measure] over the cell by stratified IMF quantiles.
        points = [_power_law_quantile(a, b, alpha, (j + 0.5) / _FRACTION_POINTS_PER_CELL)
                  for j in range(_FRACTION_POINTS_PER_CELL)]
        bright += count * sum(_bright_measure(m, min_luminosity_sol, window) for m in points) / len(points)
        living += count * sum(_living_measure(m, window) for m in points) / len(points)
        bound = _bright_measure_bound(a, b, min_luminosity_sol, window)
        if bound > 0:
            running += count * bound
            cells.append((a, b, alpha, bound))
            cumulative.append(running)
    return {"window": window, "cells": cells, "cumulative": cumulative,
            "fraction": bright / living if living > 0 else 0.0}


def _power_law_quantile(a, b, alpha, u):
    """The `u` quantile of a power law of slope `alpha` on `[a, b]`."""
    if alpha == 1.0:
        return a * (b / a) ** u
    p = 1 - alpha
    return (a ** p + u * (b ** p - a ** p)) ** (1 / p)


def bright_star_fraction(min_luminosity_sol, population=None):
    """
    The share of a population's living stars (those `sample_living_star`
    draws) whose luminosity is at least `min_luminosity_sol`. Computed once
    per threshold and population, then cached.

    Args:
        min_luminosity_sol (float): The threshold (Lsun), above the
            brightest white dwarf.
        population (str or None): "young", "intermediate", "old", "bulge",
            or None for the whole disk.

    Returns:
        float: The fraction, in `[0, 1]`.
    """
    return _bright_table(float(min_luminosity_sol), population)["fraction"]


def _sample_one_bright(table, min_luminosity_sol, rng):
    cells, cumulative, window = table["cells"], table["cumulative"], table["window"]
    if not cells:
        raise ValueError(f"no star in population window {window} Gy reaches {min_luminosity_sol} Lsun")
    for _ in range(program_constants.STAR_MODEL_MAX_REDRAWS):
        index = min(bisect.bisect_right(cumulative, rng.random() * cumulative[-1]), len(cells) - 1)
        a, b, alpha, bound = cells[index]
        mass = _sample_power_law(a, b, alpha, rng)
        pieces = [(max(start, window[0]), min(end, window[1]), chance)
                  for start, end, chance in _bright_windows(mass, min_luminosity_sol)]
        weights = [max(0.0, end - start) * chance for start, end, chance in pieces]
        total = sum(weights)
        if rng.random() * bound >= total:
            continue
        pick = rng.random() * total
        for (start, end, _), weight in zip(pieces, weights):
            if pick < weight:
                break
            pick -= weight
        age = rng.uniform(start, end)
        state = evolve_star(mass, age, rng, min_luminosity_sol=min_luminosity_sol)
        # Guards the window edges against rounding (an age landing exactly
        # on a phase boundary).
        if state is not None and state["luminosity_sol"] >= min_luminosity_sol:
            return star_params(mass, age, state)
    raise ValueError(f"no bright star drawn in {program_constants.STAR_MODEL_MAX_REDRAWS} tries")


def sample_bright_stars(n, min_luminosity_sol, population=None, rng=random):
    """
    Draws `n` stars from the population model conditional on luminosity
    `>= min_luminosity_sol` (see the module docstring for how).

    Args:
        n (int): How many to draw.
        min_luminosity_sol (float): The threshold (Lsun), above the
            brightest white dwarf.
        population (str or None): As in `bright_star_fraction`.
        rng (random.Random): The random source (the module-level one by
            default); a seeded one gives the same stars every time.

    Returns:
        list: `n` dicts from `stellarEvolution.star_params` (type, Yerkes
              class, mass, radius, temperature, luminosity, age, lifespan,
              initial mass and phase end age), each for `Star.from_params`.
    """
    table = _bright_table(float(min_luminosity_sol), population)
    return [_sample_one_bright(table, min_luminosity_sol, rng) for _ in range(n)]


def sample_dim_star(max_luminosity_sol, population=None, rng=random):
    """
    Draws one star from the population model conditional on luminosity
    `< max_luminosity_sol` (a sector's own stars once its bright ones were
    pre-placed). Nearly every star is that dim, so it simply redraws the
    rare bright one.

    Returns:
        dict: As one of `sample_bright_stars`'s.
    """
    mass, age, state = sample_living_star(rng=rng, population=population, max_luminosity_sol=max_luminosity_sol)
    return star_params(mass, age, state)


def pick_population(densities, rng=random):
    """
    Picks a population for one star at a position, in proportion to
    `galaxyDensity.population_densities` there.

    Args:
        densities (dict): Population name -> density, as
            `population_densities` returns.

    Returns:
        str or None: The population's name; None when every density is
                     zero (the whole disk's ages).
    """
    total = sum(densities.values())
    if total <= 0:
        return None
    pick = rng.random() * total
    for name, density in densities.items():
        if pick < density:
            return name
        pick -= density
    return name
