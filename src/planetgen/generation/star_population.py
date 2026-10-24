# planetgen/generation/star_population.py

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

`age_laws` (GEN.546, optional on every function that takes a `population`) is a dict of population name to
`galaxy.age_law.AgeLaw`, `galaxy.disk_history.population_age_laws(shape)` for a galaxy: the density of ages inside the
population's window. A population without an entry draws its ages evenly across the window, as every population did
before. The bright measure, the tables and the age draw all use the same law, so a scatter and a sector fill agree.

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
50 Lsun), so the threshold must be at least that range's top, and a
sector's dim draw keeps every white dwarf whatever the cap.

Mass bands (GEN.187)
--------------------
The backfill around generated sectors fills by initial mass, not
luminosity: `mass_band_fraction` and `sample_mass_band_stars` take every
living star born in `[low, high)` solar masses (any luminosity, a white
dwarf included) below an optional luminosity ceiling (what the luminosity
pass already placed). Their sampling table is the bright table's
counterpart: each cell is weighted by its IMF share times the longest
age window a star in it is still alive in.
"""

import bisect
import functools
import math

from planetgen.galaxy.age_law import AgeLaw
from planetgen.galaxy.component_field import POPULATIONS_IN_THE_GATE
from planetgen.physics import constants, wolf_rayet
from planetgen import tuning
from planetgen.physics.stellar_evolution import (
    _power_law_integral, _sample_power_law, agb_bright_chance, birth_fe_h, evolve_star, giant_bright_chance,
    giant_end_fraction,
    has_core_helium_phases, horizontal_branch_bright_chance, main_sequence_lifetime_gy, main_sequence_luminosity_sol,
    population_age_range_gy, sample_living_star, star_params, subgiant_end_fraction, subgiant_end_luminosity_sol,
)
from planetgen.util import draw

POPULATIONS = tuple(tuning.STELLAR_POPULATION_AGE_RANGES_GY)
"""The named stellar populations, in `population_densities` order."""


def _overlap(low, high, law):
    """The window's length times the chance an age falls in `[low, high]` (the overlap itself, for an even spread)."""
    return law.weight(low, high)


def _law(population, age_laws):
    """The `AgeLaw` of `population`: its entry in `age_laws`, else an even spread across its window."""
    found = None if age_laws is None else age_laws.get(population)
    return found if found is not None else AgeLaw.uniform(population_age_range_gy(population))


def _giant_bright_chance(mass_sol, min_luminosity_sol, fe_h=0.0):
    """The chance a giant's luminosity is >= the threshold
    (`stellarEvolution.giant_bright_chance`)."""
    return giant_bright_chance(mass_sol, min_luminosity_sol, fe_h)


def _supergiant_bright_chance(ms_luminosity_sol, min_luminosity_sol):
    """The chance a supergiant's uniform luminosity growth takes it >= the threshold."""
    low, high = tuning.SUPERGIANT_LUMINOSITY_GROWTH_RANGE
    need = min_luminosity_sol / ms_luminosity_sol
    if need <= low:
        return 1.0
    if need >= high:
        return 0.0
    return (high - need) / (high - low)


def _wolf_rayet_bright_chance(min_luminosity_sol):
    """The chance a Wolf-Rayet star's log-uniform luminosity is >= the threshold (GEN.243)."""
    low, high = tuning.WR_LUMINOSITY_RANGE_SOL
    if min_luminosity_sol <= low:
        return 1.0
    if min_luminosity_sol >= high:
        return 0.0
    return math.log(high / min_luminosity_sol) / math.log(high / low)


def _bright_windows(mass_sol, min_luminosity_sol, fe_h=0.0):
    """
    Every age window in which a star of initial mass `mass_sol` and birth [Fe/H] `fe_h` can be at
    least `min_luminosity_sol` bright, as `(start_gy, end_gy, chance)`:
    within the window its age puts it in a phase whose luminosity clears
    the threshold with probability `chance`. Mirrors `evolve_star`.
    """
    pc = tuning
    t_ms = main_sequence_lifetime_gy(mass_sol, fe_h)
    t_sub = t_ms * subgiant_end_fraction(mass_sol, fe_h)
    t_giant = t_ms * giant_end_fraction(mass_sol, fe_h)
    l_ms = main_sequence_luminosity_sol(mass_sol)
    windows = []
    if l_ms >= min_luminosity_sol:
        windows.append((0.0, t_ms, 1.0))
    if mass_sol >= pc.SUPERGIANT_MIN_MASS_SOL:
        chance = _supergiant_bright_chance(l_ms, min_luminosity_sol)
        wr_start = wolf_rayet.start_gy(mass_sol, t_ms, t_giant, fe_h)
        if wr_start is None:
            wr_start = t_giant
        elif _wolf_rayet_bright_chance(min_luminosity_sol) > 0:
            windows.append((wr_start, t_giant, _wolf_rayet_bright_chance(min_luminosity_sol)))
        if chance > 0:
            windows.append((t_ms, wr_start, chance))
        return windows
    # A subgiant brightens linearly from l_ms to the base of its giant branch.
    progress = (min_luminosity_sol - l_ms) / (subgiant_end_luminosity_sol(mass_sol, fe_h) - l_ms)
    if progress < 1:
        windows.append((t_ms + (t_sub - t_ms) * max(progress, 0.0), t_sub, 1.0))
    chance = _giant_bright_chance(mass_sol, min_luminosity_sol, fe_h)
    if chance > 0:
        windows.append((t_sub, t_giant, chance))
    if has_core_helium_phases(mass_sol):
        t_horizontal = t_giant + pc.HORIZONTAL_BRANCH_DURATION_GY
        t_agb = t_horizontal + t_ms * pc.AGB_MS_FRACTION
        chance = horizontal_branch_bright_chance(min_luminosity_sol)
        if chance > 0:
            windows.append((t_giant, t_horizontal, chance))
        chance = agb_bright_chance(min_luminosity_sol)
        if chance > 0:
            windows.append((t_horizontal, t_agb, chance))
    return windows


def _bright_measure(mass_sol, min_luminosity_sol, window, fe_h=0.0):
    """The chance a star of this mass and birth [Fe/H], of an age spread over the age law `window`, is
    at least `min_luminosity_sol` bright (times the window's length)."""
    return sum(_overlap(start, end, window) * chance
               for start, end, chance in _bright_windows(mass_sol, min_luminosity_sol, fe_h))


def _living_measure(mass_sol, window, fe_h=0.0):
    """The chance a star of this mass and birth [Fe/H], of an age spread over the age law `window`, has
    not collapsed (times the window's length); see `evolve_star`."""
    span = window.span
    if mass_sol < tuning.SUPERGIANT_MIN_MASS_SOL:
        return span
    dead_from = main_sequence_lifetime_gy(mass_sol, fe_h) * giant_end_fraction(mass_sol, fe_h)
    return span - _overlap(dead_from, math.inf, window)


@functools.lru_cache(maxsize=None)
def _light_star_phase_envelope(low_mass, high_mass, fe_h=0.0):
    """
    `(sub_low, sub_high, giant_low, giant_high)`: the least and greatest `subgiant_end_fraction` and
    `giant_end_fraction` of a star of birth [Fe/H] `fe_h` lighter than `HELIUM_FLASH_MASS_SOL` between `low_mass` and
    `high_mass` (GEN.502),
    over a grid of masses across the cell and a percent either side. The fractions move with the mass, so a cell's
    bound takes them over its own masses; a bound over all light stars would be a hundred times too loose.
    """
    top = min(high_mass, tuning.HELIUM_FLASH_MASS_SOL)
    masses = [low_mass * (top / low_mass) ** (i / _ENVELOPE_STEPS) for i in range(_ENVELOPE_STEPS + 1)]
    subgiant = [subgiant_end_fraction(m, fe_h) for m in masses]
    giant = [giant_end_fraction(m, fe_h) for m in masses]
    return min(subgiant) * 0.99, max(subgiant) * 1.01, min(giant) * 0.99, max(giant) * 1.01


_ENVELOPE_STEPS = 16
"""int: How many steps across a cell the phase fractions are sampled for its bound."""


def _phase_fraction_bounds(low_mass, high_mass, fe_h=0.0):
    """`(sub_low, sub_high, giant_low, giant_high)` over the masses `low_mass` to `high_mass` of one cell, which no
    phase-rule mass splits: the light stars' envelope, or the constants of the heavier ones."""
    pc = tuning
    heavy = (pc.SUBGIANT_PHASE_END_MS_FRACTION, pc.SUBGIANT_PHASE_END_MS_FRACTION,
             pc.GIANT_PHASE_END_MS_FRACTION, pc.GIANT_PHASE_END_MS_FRACTION)
    if low_mass >= pc.HELIUM_FLASH_MASS_SOL:
        return heavy
    light = _light_star_phase_envelope(low_mass, high_mass, fe_h)
    if high_mass <= pc.HELIUM_FLASH_MASS_SOL:
        return light
    return tuple(min(a, b) if i % 2 == 0 else max(a, b) for i, (a, b) in enumerate(zip(light, heavy)))


def _bright_measure_terms(low_mass, high_mass, min_luminosity_sol, fe_h):
    """
    The pieces of `_bright_measure_bound` for every mass in `[low_mass, high_mass]`, a cell that no phase-rule mass
    (`_mass_grid`'s breaks) splits, at birth [Fe/H] `fe_h`: `{phase: (start_gy, end_gy, chance)}`, each phase's window
    stretched over every mass in the cell (lifetimes fall and luminosities rise with mass) at its best chance.
    """
    pc = tuning
    t_long = main_sequence_lifetime_gy(low_mass, fe_h)
    t_short = main_sequence_lifetime_gy(high_mass, fe_h)
    sub_low, sub_high, giant_low, giant_high = _phase_fraction_bounds(low_mass, high_mass, fe_h)
    l_max = main_sequence_luminosity_sol(high_mass)
    terms = {}
    if l_max >= min_luminosity_sol:
        terms["main sequence"] = (0.0, t_long, 1.0)
    if low_mass >= pc.SUPERGIANT_MIN_MASS_SOL:
        chance = _supergiant_bright_chance(l_max, min_luminosity_sol)
        terms["supergiant"] = (t_short, t_long * pc.GIANT_PHASE_END_MS_FRACTION, chance)
        if high_mass >= wolf_rayet.threshold_mass_sol(fe_h):
            # GEN.243: the Wolf-Rayet window at the end of the supergiant window, over every mass of the cell.
            terms["Wolf-Rayet"] = (
                t_short * pc.GIANT_PHASE_END_MS_FRACTION - tuning.WR_PHASE_MYR / wolf_rayet.MYR_PER_GY,
                t_long * pc.GIANT_PHASE_END_MS_FRACTION, _wolf_rayet_bright_chance(min_luminosity_sol))
        return terms
    if subgiant_end_luminosity_sol(high_mass, fe_h) > min_luminosity_sol:
        terms["subgiant"] = (t_short, t_long * sub_high, 1.0)
    # A heavier giant starts its branch higher, so the chance at the cell's top mass bounds the cell's, except right
    # at a light star's tip, which falls a little with the mass: the greatest of three masses covers that.
    chance = max(_giant_bright_chance(mass, min_luminosity_sol, fe_h)
                 for mass in (low_mass, math.sqrt(low_mass * high_mass), high_mass))
    terms["giant"] = (t_short * sub_low, t_long * giant_high, chance)
    if has_core_helium_phases(high_mass):
        horizontal = pc.HORIZONTAL_BRANCH_DURATION_GY
        terms["horizontal branch"] = (t_short * giant_low, t_long * giant_high + horizontal,
                                      horizontal_branch_bright_chance(min_luminosity_sol))
        terms["asymptotic branch"] = (t_short * giant_low + horizontal,
                                      t_long * (giant_high + pc.AGB_MS_FRACTION) + horizontal,
                                      agb_bright_chance(min_luminosity_sol))
    return terms


def _bright_measure_bound(low_mass, high_mass, min_luminosity_sol, window, fe_hs=(0.0,)):
    """
    An upper bound on `_bright_measure` for every mass in `[low_mass, high_mass]` at every birth [Fe/H] in `fe_hs`
    (a table's bin: its edges and centre): each phase's window stretched over every mass of the cell and, as a
    lifetime moves with the composition, over every composition, at its best chance. Taking the bound at each
    composition and the greatest would miss a window that slides across the edge of the age law between them.
    """
    stretched = {}
    for fe_h in fe_hs:
        for phase, (start, end, chance) in _bright_measure_terms(low_mass, high_mass, min_luminosity_sol,
                                                                 fe_h).items():
            low, high, best = stretched.get(phase, (math.inf, -math.inf, 0.0))
            stretched[phase] = (min(low, start), max(high, end), max(best, chance))
    return sum(chance * _overlap(start, end, window) for start, end, chance in stretched.values())


@functools.lru_cache(maxsize=None)
def _mass_grid(extra_edges=(), cells=None):
    """Log-spaced cell edges over the IMF's range (`cells` of them, `tuning.BRIGHT_STAR_MASS_GRID_CELLS` without),
    with every mass where the IMF slope or a phase rule changes added as an edge (and `extra_edges`, a mass cut a
    draw splits at), and each cell's IMF slope and continuity factor."""
    pc = tuning
    low, high = pc.IMF_BREAKS_SOL[0], pc.IMF_BREAKS_SOL[-1]
    cells = pc.BRIGHT_STAR_MASS_GRID_CELLS if cells is None else cells
    edges = {low * (high / low) ** (i / cells) for i in range(cells + 1)}
    edges.update(pc.IMF_BREAKS_SOL)
    edges.update((pc.BRIGHT_GIANT_MIN_MASS_SOL, pc.SUPERGIANT_MIN_MASS_SOL, pc.HELIUM_FLASH_MASS_SOL))
    edges.update(extra_edges)
    edges = sorted(edge for edge in edges if low <= edge <= high)

    norms = [1.0]
    for i in range(1, len(pc.IMF_SLOPES)):
        norms.append(norms[-1] * pc.IMF_BREAKS_SOL[i] ** (pc.IMF_SLOPES[i] - pc.IMF_SLOPES[i - 1]))
    grid = []
    for a, b in zip(edges, edges[1:]):
        segment = bisect.bisect_right(pc.IMF_BREAKS_SOL, a) - 1
        grid.append((a, b, pc.IMF_SLOPES[segment], norms[segment]))
    return grid


def check_threshold(min_luminosity_sol):
    if not (isinstance(min_luminosity_sol, (int, float)) and math.isfinite(min_luminosity_sol)):
        raise ValueError(f"luminosity threshold must be a finite number, got {min_luminosity_sol!r}")
    if min_luminosity_sol < tuning.WD_LUMINOSITY_RANGE_SOL[1]:
        raise ValueError(f"luminosity threshold {min_luminosity_sol} Lsun must be at least the brightest white dwarf "
                         f"({tuning.WD_LUMINOSITY_RANGE_SOL[1]} Lsun)")


# Quadrature points per cell for `bright_star_fraction`: the cell's IMF
# quantiles at (j + 1/2) / n.
_FRACTION_POINTS_PER_CELL = 4


def _check_mass_range(mass_range):
    """A mass range as the `(low, high)` pair of solar masses `_bright_table`
    caches (a missing end `None`), or `None` for every mass."""
    if mass_range is None:
        return None
    low, high = mass_range
    if low is not None and high is not None and not low < high:
        raise ValueError(f"empty mass range [{low:g}, {high:g}) solar masses")
    if low is None and high is None:
        return None
    return (None if low is None else float(low), None if high is None else float(high))


def fe_h_bin(fe_h):
    """The centre of the `tuning.BRIGHT_TABLE_FE_H_STEP_DEX` bin of a birth [Fe/H]: the composition a sampling table is
    built at (GEN.651)."""
    step = tuning.BRIGHT_TABLE_FE_H_STEP_DEX
    return round(fe_h / step) * step + 0.0


FE_H_BINS = tuple(sorted({fe_h_bin(tuning.METALLICITY_RANGE_FE_H[0] + i * tuning.BRIGHT_TABLE_FE_H_STEP_DEX / 2.0)
                          for i in range(int(2 * (tuning.METALLICITY_RANGE_FE_H[1] - tuning.METALLICITY_RANGE_FE_H[0])
                                             / tuning.BRIGHT_TABLE_FE_H_STEP_DEX) + 2)}
                         | {fe_h_bin(tuning.METALLICITY_RANGE_FE_H[1])}))
"""tuple: Every bin a star's birth [Fe/H] can fall in (`tuning.METALLICITY_RANGE_FE_H`, clamped), ascending."""


def warm_bins(min_luminosity_sol, population, mass_range=None, age_laws=None):
    """Builds the sampling table of every [Fe/H] bin for one threshold and population (GEN.651), so a worker that
    forks a horse per job holds them all; a star draws from the table of its own bin."""
    window = _law(population, age_laws)
    for fe_bin in FE_H_BINS:
        _bright_table(float(min_luminosity_sol), population, _check_mass_range(mass_range), window, fe_bin)


def _grid_cells(fe_bin):
    """How many mass cells a table of this [Fe/H] bin has: the full grid for the reference composition's (it works out
    the share of stars the counts use), a coarser one for the others, whose bounds already stretch over the bin's
    lifetimes by more than a coarse cell's width, and which are eight times cheaper to build."""
    return tuning.BRIGHT_STAR_MASS_GRID_CELLS if fe_bin == 0.0 else tuning.BRIGHT_TABLE_BIN_GRID_CELLS


def _bin_fe_hs(fe_bin):
    """The composition a sampling table's bounds are taken over: the bin's two edges and its centre. A star in the bin
    draws from the table with its own [Fe/H], so a cell's bound has to hold anywhere in the bin."""
    half = tuning.BRIGHT_TABLE_FE_H_STEP_DEX / 2.0
    return (fe_bin - half, fe_bin, fe_bin + half)


@functools.lru_cache(maxsize=None)
def _bright_table(min_luminosity_sol, population, mass_range=None, window=None, fe_bin=0.0):
    """
    The cached per-(threshold, population, mass range, [Fe/H] bin) sampling table: cell
    edges, slopes and continuity factors, each cell's measure bound, the
    running total of cell weights (IMF share times bound), and the bright
    fraction. `mass_range` (`(low, high)` initial masses in solar masses,
    either end `None`) keeps only the stars born in it: the fraction is
    still the share of ALL living stars, so a mass pass and a luminosity
    pass below it add up to the share that one threshold gives. `window` is the population's `AgeLaw` (an even spread
    across its window when `None`). `fe_bin` is the birth [Fe/H] the table is built at (`fe_h_bin`; a star's lifetimes
    and giant branch follow its own, GEN.651), so a star is drawn from the table of its own bin and the bounds are
    valid for it. The `fraction` is worked out for the reference composition's table (bin 0) only, the one the counts
    use (a stand-in for the mixture of a place's compositions, GEN.661); another bin's is `None`.
    """
    check_threshold(min_luminosity_sol)
    if window is None:
        window = _law(population, None)
    cells, cumulative = [], []
    running = bright = living = 0.0
    mass_low, mass_high = mass_range if mass_range is not None else (None, None)
    with_fraction = fe_bin == 0.0
    for a, b, alpha, norm in _mass_grid(tuple(edge for edge in (mass_low, mass_high) if edge is not None),
                                        _grid_cells(fe_bin)):
        count = norm * _power_law_integral(a, b, alpha)
        if with_fraction:
            # E[measure] over the cell by stratified IMF quantiles.
            points = [_power_law_quantile(a, b, alpha, (j + 0.5) / _FRACTION_POINTS_PER_CELL)
                      for j in range(_FRACTION_POINTS_PER_CELL)]
            living += count * sum(_living_measure(m, window, fe_bin) for m in points) / len(points)
        if (mass_low is not None and b <= mass_low) or (mass_high is not None and a >= mass_high):
            continue
        if with_fraction:
            bright += count * sum(_bright_measure(m, min_luminosity_sol, window, fe_bin) for m in points) / len(points)
        bound = _bright_measure_bound(a, b, min_luminosity_sol, window, _bin_fe_hs(fe_bin)) * tuning.BRIGHT_TABLE_BOUND_MARGIN
        if bound > 0:
            running += count * bound
            cells.append((a, b, alpha, bound))
            cumulative.append(running)
    return {"window": window.window, "law": window, "cells": cells, "cumulative": cumulative, "population": population,
            "fraction": (bright / living if living > 0 else 0.0) if with_fraction else None}


def _power_law_quantile(a, b, alpha, u):
    """The `u` quantile of a power law of slope `alpha` on `[a, b]`."""
    if alpha == 1.0:
        return a * (b / a) ** u
    p = 1 - alpha
    return (a ** p + u * (b ** p - a ** p)) ** (1 / p)


def bright_star_fraction(min_luminosity_sol, population=None, mass_range=None, age_laws=None):
    """
    The share of a population's living stars (those `sample_living_star`
    draws) whose luminosity is at least `min_luminosity_sol`, among those
    born in `mass_range` (initial masses `(low, high)` in solar masses, an
    end `None` for no limit; every mass by default). Computed once per
    threshold, population and range, then cached.

    Args:
        min_luminosity_sol (float): The threshold (Lsun), above the
            brightest white dwarf.
        population (str or None): "young", "intermediate", "old", "bulge",
            or None for the whole disk.
        mass_range (tuple, optional): As above.
        age_laws (dict, optional): Population name -> `AgeLaw`, as in the module docstring.

    Returns:
        float: The fraction, in `[0, 1]`.
    """
    return _bright_table(float(min_luminosity_sol), population, _check_mass_range(mass_range),
                         _law(population, age_laws))["fraction"]


def massive_star_fraction(min_mass_sol, population=None, age_laws=None):
    """
    The share of a population's living stars born with at least
    `min_mass_sol` solar masses, at any luminosity: what the mass pass of
    the scatter places (a star born with 8 or more is brighter than
    `MASS_PASS_MIN_LUMINOSITY_SOL` all its life), and what a mass backfill
    (GEN.187) takes a sector down to.
    """
    return _mass_band_table((float(min_mass_sol), None), population, _law(population, age_laws))["fraction"]


def _prior_source(population, rng, fe_h_source, index):
    """The draw of a birth [Fe/H] for the `index`-th star: the place's (`fe_h_source(index)`, `sample_bright_stars`)
    or, with no place, one for `population` at the Sun's distance. A star of a population in the gate takes its
    component's [Fe/H] later and the row stores none (GEN.159, GEN.237), so it is drawn at the reference composition
    its tree reads it back at (`multistar.Placed.tree_fe_h`): a stand-in until the component's value reaches the draw
    (GEN.662)."""
    if fe_h_source is None:
        if population in POPULATIONS_IN_THE_GATE:
            return lambda: 0.0
        return lambda: birth_fe_h(rng, population)
    return lambda: fe_h_source(index)


def _table_total(table):
    """What a table draws a try's cell from: its cells' running total of weight (IMF share times bound), 0 for a table
    with none."""
    return table["cumulative"][-1] if table["cells"] else 0.0


FE_H_REDRAWS = 100 * tuning.STAR_MODEL_MAX_REDRAWS
"""int: How many compositions `_draw_fe_h` tries for one that a table keeps (a metal-poor place against a metal-rich
table total, rare in the luminosity window it asks for)."""


def _draw_fe_h(source, totals, rng):
    """
    The birth [Fe/H] for a try at a star the tables condition on (it is bright, or in a mass band): `source()` draws one
    from the place's distribution and it is kept with the chance its [Fe/H] bin's table total (`totals`, a dict by bin) is
    of the greatest bin's. A try then draws a cell and a mass from the table of that bin and keeps the star with its
    measure over its cell's bound, so a composition is kept as often as its stars are there: the composition a place's
    bright stars have is its prior times each composition's bright share, not the prior (a metal-poor halo has no 2500
    Lsun giants in a 3 to 10 Gy window, and its metal-rich tail makes all of them). A try that fails draws a new
    composition with the next.
    """
    best = max(totals.values())
    if best <= 0:
        raise ValueError("no composition makes a star of this kind in this population window")
    for _ in range(FE_H_REDRAWS):
        fe_h = source()
        if rng.random() * best < totals.get(fe_h_bin(fe_h), 0.0):
            return fe_h
    raise ValueError(f"no birth [Fe/H] kept in {FE_H_REDRAWS} tries")


def _bin_totals(table_at):
    """What each [Fe/H] bin's table draws a try's cell from (`_table_total`), by bin: one lookup of each of the tables
    (built on first use), made once for a star, not once for each composition tried."""
    return {fe_bin: _table_total(table_at(fe_bin)) for fe_bin in FE_H_BINS}


def _sample_one_bright(min_luminosity_sol, population, mass_range, window, rng, source):
    def table_at(fe_bin):
        return _bright_table(float(min_luminosity_sol), population, mass_range, window, fe_bin)

    for _ in range(tuning.STAR_MODEL_MAX_REDRAWS):
        fe_h = _draw_fe_h(source, _bin_totals(table_at), rng)
        table = table_at(fe_h_bin(fe_h))
        cells, cumulative, law = table["cells"], table["cumulative"], table["law"]
        span = law.window
        index = min(bisect.bisect_right(cumulative, rng.random() * cumulative[-1]), len(cells) - 1)
        a, b, alpha, bound = cells[index]
        mass = _sample_power_law(a, b, alpha, rng)
        pieces = [(max(start, span[0]), min(end, span[1]), chance)
                  for start, end, chance in _bright_windows(mass, min_luminosity_sol, fe_h)]
        weights = [max(0.0, end - start) * chance for start, end, chance in pieces]
        total = sum(weights)
        if rng.random() * bound >= total:
            continue
        pick = rng.random() * total
        for (start, end, _), weight in zip(pieces, weights):
            if pick < weight:
                break
            pick -= weight
        age = law.sample(start, end, rng)
        state = evolve_star(mass, age, rng, min_luminosity_sol=min_luminosity_sol, fe_h=fe_h,
                            population=table["population"])
        # Guards the window edges against rounding (an age landing exactly
        # on a phase boundary).
        if state is not None and state["luminosity_sol"] >= min_luminosity_sol:
            return star_params(mass, age, state)
    raise ValueError(f"no bright star drawn in {tuning.STAR_MODEL_MAX_REDRAWS} tries")


BAND_REDRAWS = 10000
"""int: How many stars `sample_bright_stars` draws above a band's lower
limit before giving up on one landing below its upper limit (a band so
thin almost every star overshoots it)."""


MASS_PASS_MIN_LUMINOSITY_SOL = tuning.WD_LUMINOSITY_RANGE_SOL[1]
"""float: The luminosity floor of the mass pass: the lowest `check_threshold`
allows (the brightest white dwarf). A star born with 8 solar masses or more
is brighter than that through all of its life, so the floor cuts none of
them."""


def _sample_one_in_band(min_luminosity_sol, max_luminosity_sol, population, mass_range, window, rng, source):
    """One star with `min <= luminosity < max`: a star drawn above `min`,
    redrawn (its composition with it) while it reaches `max` (the brighter band already placed)."""
    max_w = max_luminosity_sol * constants.SOLAR_LUMINOSITY
    for _ in range(BAND_REDRAWS):
        params = _sample_one_bright(min_luminosity_sol, population, mass_range, window, rng, source)
        if params["luminosity_w"] < max_w:
            return params
    raise ValueError(f"no star between {min_luminosity_sol:g} and {max_luminosity_sol:g} Lsun drawn "
                     f"in {BAND_REDRAWS} tries")


def sample_bright_stars(n, min_luminosity_sol, population=None, rng=draw, max_luminosity_sol=None,
                        mass_range=None, age_laws=None, fe_h_source=None):
    """
    Draws `n` stars from the population model conditional on luminosity
    `>= min_luminosity_sol` (see the module docstring for how), and below
    `max_luminosity_sol` when given (one band of a staged scatter), and
    born in `mass_range` (`(low, high)` initial solar masses, an end
    `None` for no limit) when given.

    Args:
        n (int): How many to draw.
        min_luminosity_sol (float): The threshold (Lsun), above the
            brightest white dwarf.
        population (str or None): As in `bright_star_fraction`.
        rng (draw.Stream): The random source (the module-level one by
            default); a seeded one gives the same stars every time.
        max_luminosity_sol (float, optional): The band's upper limit
            (Lsun, exclusive); `None` for no limit.
        mass_range (tuple, optional): As in `bright_star_fraction`.
        age_laws (dict, optional): As in `bright_star_fraction`.
        fe_h_source (callable, optional): `fe_h_source(index)` draws a birth [Fe/H] from the distribution of the
            place of the `index`-th star (`birth_fe_h` there, GEN.651); each star's lifetimes and giant branch follow
            its own, and it is kept with the chance that composition makes such a star (`_draw_fe_h`). Without it
            each star draws from `population` at the Sun's distance. The sampling tables are built per [Fe/H] bin
            (`fe_h_bin`), so each star is drawn from the table of its own bin.

    Returns:
        list: `n` dicts from `stellarEvolution.star_params` (type, Yerkes
              class, mass, radius, temperature, luminosity, age, lifespan,
              initial mass, phase end age and birth [Fe/H]), each for `Star.from_params`.
    """
    if max_luminosity_sol is not None and max_luminosity_sol <= min_luminosity_sol:
        raise ValueError(f"empty luminosity band [{min_luminosity_sol:g}, {max_luminosity_sol:g}) Lsun")
    mass_range, window = _check_mass_range(mass_range), _law(population, age_laws)
    if max_luminosity_sol is None:
        return [_sample_one_bright(min_luminosity_sol, population, mass_range, window, rng,
                                   _prior_source(population, rng, fe_h_source, i)) for i in range(n)]
    return [_sample_one_in_band(min_luminosity_sol, max_luminosity_sol, population, mass_range, window, rng,
                                _prior_source(population, rng, fe_h_source, i)) for i in range(n)]


def bright_band_fraction(min_luminosity_sol, max_luminosity_sol=None, population=None, mass_range=None,
                         age_laws=None):
    """
    The share of a population's living stars with `min <= luminosity <
    max` (`bright_star_fraction(min)` when `max` is `None`) among those
    born in `mass_range`: what a staged scatter's new band adds.
    """
    fraction = bright_star_fraction(min_luminosity_sol, population, mass_range, age_laws)
    if max_luminosity_sol is not None:
        fraction -= bright_star_fraction(max_luminosity_sol, population, mass_range, age_laws)
    return max(fraction, 0.0)


@functools.lru_cache(maxsize=None)
def _mass_band_table(mass_range, population, window=None, fe_bin=0.0):
    """
    The cached sampling table of the stars born in `mass_range` (`(low,
    high)` initial solar masses, `high` `None` for no limit) that are still
    alive, at any luminosity: the cells, their running weights (IMF share
    times the longest age window in the cell), the population's age
    window and the band's `fraction` of ALL living stars (the same
    denominator `bright_star_fraction` uses, so the two subtract). `window` is the population's `AgeLaw` (an even
    spread across its window when `None`).
    """
    if window is None:
        window = _law(population, None)
    mass_low, mass_high = mass_range
    cells, cumulative = [], []
    running = born = living = 0.0
    with_fraction = fe_bin == 0.0
    for a, b, alpha, norm in _mass_grid(tuple(edge for edge in mass_range if edge is not None), _grid_cells(fe_bin)):
        count = norm * _power_law_integral(a, b, alpha)
        measure = 0.0
        if with_fraction:
            points = [_power_law_quantile(a, b, alpha, (j + 0.5) / _FRACTION_POINTS_PER_CELL)
                      for j in range(_FRACTION_POINTS_PER_CELL)]
            measure = count * sum(_living_measure(m, window, fe_bin) for m in points) / len(points)
        living += measure
        if b <= mass_low or (mass_high is not None and a >= mass_high):
            continue
        born += measure
        # Falls with mass, and no cell straddles the supergiant mass.
        bound = max(_living_measure(a, window, fe_h) for fe_h in _bin_fe_hs(fe_bin)) * tuning.BRIGHT_TABLE_BOUND_MARGIN
        if bound > 0:
            running += count * bound
            cells.append((a, b, alpha, bound))
            cumulative.append(running)
    return {"window": window.window, "law": window, "cells": cells, "cumulative": cumulative, "population": population,
            "fraction": (born / living if living > 0 else 0.0) if with_fraction else None}


def _check_mass_band(min_mass_sol, max_mass_sol):
    if min_mass_sol is None or not min_mass_sol > 0:
        raise ValueError(f"a mass band needs a positive lower mass, got {min_mass_sol!r}")
    if max_mass_sol is not None and not min_mass_sol < max_mass_sol:
        raise ValueError(f"empty mass band [{min_mass_sol:g}, {max_mass_sol:g}) solar masses")
    return (float(min_mass_sol), None if max_mass_sol is None else float(max_mass_sol))


def mass_band_fraction(min_mass_sol, max_mass_sol=None, max_luminosity_sol=None, population=None, age_laws=None):
    """
    The share of a population's living stars born with `min_mass_sol <=
    mass < max_mass_sol` (`None`: no upper mass) that are dimmer than
    `max_luminosity_sol` (`None`: any luminosity; a white dwarf is always
    dimmer): what a mass backfill adds to a sector whose luminosity pass
    already placed the brighter ones.
    """
    band = _check_mass_band(min_mass_sol, max_mass_sol)
    born = _mass_band_table(band, population, _law(population, age_laws))["fraction"]
    if max_luminosity_sol is not None:
        born -= bright_star_fraction(max_luminosity_sol, population, band, age_laws)
    return max(born, 0.0)


def _sample_one_in_mass_band(band, population, law, max_luminosity_sol, rng, source):
    def table_at(fe_bin):
        return _mass_band_table(band, population, law, fe_bin)

    for _ in range(tuning.STAR_MODEL_MAX_REDRAWS):
        fe_h = _draw_fe_h(source, _bin_totals(table_at), rng)
        table = table_at(fe_h_bin(fe_h))
        cells, cumulative = table["cells"], table["cumulative"]
        window = table["law"].window
        index = min(bisect.bisect_right(cumulative, rng.random() * cumulative[-1]), len(cells) - 1)
        a, b, alpha, bound = cells[index]
        mass = _sample_power_law(a, b, alpha, rng)
        if rng.random() * bound >= _living_measure(mass, law, fe_h):
            continue
        oldest = window[1]
        if mass >= tuning.SUPERGIANT_MIN_MASS_SOL:
            oldest = min(oldest, main_sequence_lifetime_gy(mass, fe_h) * giant_end_fraction(mass, fe_h))
        age = law.sample(window[0], oldest, rng)
        state = evolve_star(mass, age, rng, fe_h=fe_h, population=table["population"])
        if state is None:
            continue
        if (max_luminosity_sol is not None and state["yerkes_class"] != "VII"
                and state["luminosity_sol"] >= max_luminosity_sol):
            continue
        return star_params(mass, age, state)
    raise ValueError(f"no star drawn in {tuning.STAR_MODEL_MAX_REDRAWS} tries")


def sample_mass_band_stars(n, min_mass_sol, max_mass_sol=None, population=None, rng=draw, max_luminosity_sol=None,
                           age_laws=None, fe_h_source=None):
    """
    Draws `n` living stars born with `min_mass_sol <= mass < max_mass_sol`
    and, when given, dimmer than `max_luminosity_sol` (the bright ones
    were placed already), as `stellarEvolution.star_params` dicts. `fe_h_source`: as in `sample_bright_stars`.
    """
    band = _check_mass_band(min_mass_sol, max_mass_sol)
    law = _law(population, age_laws)
    return [_sample_one_in_mass_band(band, population, law, max_luminosity_sol, rng,
                                     _prior_source(population, rng, fe_h_source, i)) for i in range(n)]


def placed_star_fraction(min_luminosity_sol, min_mass_sol, population=None, age_laws=None):
    """
    The share of a population's living stars the galaxy scatter places in
    advance: those born with at least `min_mass_sol` solar masses (the mass
    pass; `None`: no mass pass) and, of the lighter ones, those at least
    `min_luminosity_sol` bright (`None`: no luminosity pass). A sector's own
    draw is the rest.
    """
    if min_mass_sol is None:
        return 0.0 if min_luminosity_sol is None else bright_star_fraction(min_luminosity_sol, population,
                                                                           age_laws=age_laws)
    if min_luminosity_sol is None:
        return massive_star_fraction(min_mass_sol, population, age_laws)
    return (massive_star_fraction(min_mass_sol, population, age_laws)
            + bright_star_fraction(min_luminosity_sol, population, (None, min_mass_sol), age_laws))


def sample_dim_star(max_luminosity_sol, population=None, rng=draw, max_mass_sol=None, age_laws=None,
                    fe_h_peaks=None, galactic_center_dist_ly=None):
    """
    Draws one star from the population model conditional on luminosity
    `< max_luminosity_sol` (a sector's own stars once its bright ones were
    pre-placed). Nearly every star is that dim, so it simply redraws the
    rare bright one.

    Returns:
        dict: As one of `sample_bright_stars`'s.
    """
    mass, age, state = sample_living_star(rng=rng, population=population, max_luminosity_sol=max_luminosity_sol,
                                          max_mass_sol=max_mass_sol, fe_h_peaks=fe_h_peaks,
                                          galactic_center_dist_ly=galactic_center_dist_ly,
                                          age_law=None if age_laws is None else age_laws.get(population))
    return star_params(mass, age, state)


def pick_population(densities, rng=draw):
    """
    Picks a population for one star at a position, in proportion to
    `galaxy.density.population_densities` there.

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
