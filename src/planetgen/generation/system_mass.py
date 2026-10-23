# planetgen/generation/system_mass.py

"""
The mass of a system, from the generator's own distributions (GEN.545).

A galaxy's stellar mass sets its number of systems once the density is normalized to it (`galaxy.density`, GEN.404), so
the mean mass of a system is the divisor of every galaxy's star count. It used to be a recalled constant (0.4 solar
masses, "within a factor of 1.5"); here it is worked out from the draws the generator makes, for stars formed at a given
age:

1. **The placed star.** An initial mass from the initial mass function (`tuning.IMF_BREAKS_SOL`, `IMF_SLOPES`), taken to
   its state at the age: still on the main sequence at its birth mass, a white dwarf by the initial-final relation, a
   neutron star or a black hole at the mean of the kind's mass range (`present_day_moments`). The remnants are systems
   too: a remnant is the placed object of the system it left.
2. **The other stars.** The number of stars of its system is the `star-count` table's, by the placed star's class and
   with the 8 to 12 tail (`multistar.expected_star_count`); each further star has a birth mass of the placed star's
   times a uniform 0.1 to 1 (`multistar.companion_ratio_range`) and is taken to its own state at the same age.

Not counted: a compact companion added to a massive star's tree (1% of the stars of 2 to 14 solar masses a white dwarf,
0.03% a neutron star, 0.02% a black hole: under 0.2% of a system's mass, `COMPACT_COMPANION_SHARE_BOUND`), the companions
of the remnants the scatter places (a few percent of systems whose remnant has a companion at all), planets, and the
mass the giant branch's wind takes (a living star is at its birth mass, which `present_day_moments` has always done).

The mean over a population's ages weights each age by the star formation there (`cohort_mean_system_mass_solar`): the
disk's by its history (`galaxy.age_law`), the others uniformly. The mean over a galaxy or a cell weights each population
by its systems (`galaxy.system_mass`): the number of systems is the density's, so the mean is the average of the
populations' by their densities, not the mass-weighted one.
"""

import functools
import math

import numpy as np

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.physics import ms_lifetime, remnant_mass, remnant_progenitors, stellar_evolution

COHORT_AGES = 24
"""int: Ages of the midpoint rule that averages a population's cohorts, which formed at a steady rate over its range."""

REMNANT_GRID_POINTS = 400
"""int: Points of the log grid the white dwarf masses are integrated on (the initial-final relation is clipped)."""

MASS_BINS = 1200
"""int: Bins of the log grid of placed-star masses the companions are summed over (about 0.6% of the mass each)."""

COMPANION_MASS_STEPS = 32
"""int: Midpoints of the uniform draw of a further star's birth mass over the placed star's."""

COMPACT_COMPANION_SHARE_BOUND = 0.002
"""float: The most the compact companions of the massive stars' trees add to the mean mass, as a share of it (see the
module docstring); a test holds the figure."""


def imf_segments():
    """`((low, high, alpha, k), ...)`: the generator's initial mass function (`tuning.IMF_BREAKS_SOL`, `IMF_SLOPES`),
    each segment `k M^-alpha`, joined continuously."""
    breaks, slopes = tuning.IMF_BREAKS_SOL, tuning.IMF_SLOPES
    norms = [1.0]
    for i in range(1, len(slopes)):
        norms.append(norms[-1] * breaks[i] ** (slopes[i] - slopes[i - 1]))
    return tuple((breaks[i], breaks[i + 1], slopes[i], norms[i]) for i in range(len(slopes)))


def imf_integral(power, low, high):
    """`integral of M^power xi(M) dM` from `low` to `high` (initial masses, solar masses; `xi` the IMF)."""
    total = 0.0
    for seg_low, seg_high, alpha, k in imf_segments():
        a, b = max(low, seg_low), min(high, seg_high)
        if a >= b:
            continue
        exponent = power + 1.0 - alpha
        total += k * (math.log(b / a) if exponent == 0.0 else (b ** exponent - a ** exponent) / exponent)
    return total


def present_day_moments(age_gy, remnants=True):
    """
    `(sum of M, sum of M squared)` over everything a cohort of one age holds now, per unit of the IMF's number: the
    stars still on the main sequence (initial mass up to the turn-off, where the main-sequence lifetime reaches the
    age) and the remnants of the rest, by the progenitor ranges and windows the generator draws remnants from
    (`tuning.WHITE_DWARF_PROGENITOR_SOL` and the others) and their mass rules: white dwarfs by the initial-final
    relation, neutron stars and black holes uniform over their mass ranges. Neutron stars and black holes take the
    progenitors of their windows above the white dwarf range's top (`remnant_progenitors`).
    """
    top = tuning.IMF_BREAKS_SOL[-1]
    turn_off = min(top, ms_lifetime.mass_for_lifetime_gy(age_gy))
    low_edge = tuning.IMF_BREAKS_SOL[0]
    mass = imf_integral(1.0, low_edge, turn_off)
    square = imf_integral(2.0, low_edge, turn_off)
    if not remnants:
        return mass, square

    wd_low, wd_high = tuning.WHITE_DWARF_PROGENITOR_SOL
    wd_low = max(wd_low, turn_off)
    if wd_low < wd_high:
        grid = np.geomspace(wd_low, wd_high, REMNANT_GRID_POINTS)
        middle = 0.5 * (grid[:-1] + grid[1:])
        weight = np.array([imf_integral(0.0, a, b) for a, b in zip(grid[:-1], grid[1:])])
        remnant = remnant_mass.white_dwarf_mass_sol(middle)
        mass += float(np.sum(weight * remnant))
        square += float(np.sum(weight * remnant * remnant))

    for kind in remnant_progenitors.KINDS:
        number = sum(imf_integral(0.0, max(low, wd_high, turn_off), high)
                     for low, high in remnant_progenitors.kind_windows(kind) if max(low, wd_high, turn_off) < high)
        first, second = remnant_mass.uniform_moments(kind)
        mass += number * first
        square += number * second
    return mass, square


@functools.lru_cache(maxsize=256)
def cohort_moments(age_low_gy, age_high_gy, remnants=True):
    """`present_day_moments` averaged over a population whose stars formed at a steady rate between two ages."""
    low, high = sorted((max(age_low_gy, 1.0e-4), max(age_high_gy, 1.0e-4)))
    ages = low + (high - low) * (np.arange(COHORT_AGES) + 0.5) / COHORT_AGES
    moments = [present_day_moments(float(age), remnants) for age in ages]
    return (sum(m for m, _s in moments) / COHORT_AGES, sum(s for _m, s in moments) / COHORT_AGES)


def _remnant_states(initial):
    """
    What `present_mass_sol` takes from the birth masses alone, not the age: `(initial, remnant)`, `initial` an array of
    birth masses and `remnant` the mass each has once its star has ended: the white dwarf of the initial-final relation
    below the white dwarf progenitors' top, and above it a black hole of the kind's mean mass where the progenitor lies in
    a failure window (`remnant_progenitors`) and a neutron star of the kind's mean mass elsewhere. A cohort's age only
    moves the turn-off, so the means over many ages (a galaxy's thousands of clusters) work these out once.
    """
    wd_high = tuning.WHITE_DWARF_PROGENITOR_SOL[1]
    fails = np.zeros(initial.shape, dtype=bool)
    for low, high in remnant_progenitors.kind_windows(remnant_progenitors.BLACK_HOLE):
        fails |= (initial >= low) & (initial <= high)
    remnant = np.where(initial < wd_high, remnant_mass.white_dwarf_mass_sol(np.minimum(initial, wd_high)),
                       np.where(~fails, remnant_mass.uniform_moments("neutron-star")[0],
                                remnant_mass.uniform_moments("black-hole")[0]))
    return initial, remnant


def _turn_off_sol(age_gy):
    """The heaviest birth mass a star can have and still be alive at `age_gy`: the lightest progenitor whose giant phases
    have ended, capped at the IMF's top."""
    return min(tuning.IMF_BREAKS_SOL[-1], stellar_evolution.youngest_progenitor_sol(age_gy * 1.0e9))


def _present_mass(states, turn_off):
    """`present_mass_sol` for birth masses whose `_remnant_states` are `states`, for a population whose turn-off is
    `turn_off` solar masses (`_turn_off_sol`): the birth mass up to the turn-off, the remnant above it."""
    initial, remnant = states
    return np.where(initial <= turn_off, initial, remnant)


def present_mass_sol(initial_mass_sol, age_gy):
    """
    The mass at `age_gy` of a star born with `initial_mass_sol` (an array of masses), by the rules
    `present_day_moments` sums except for the turn-off: its birth mass until its last giant phase ends
    (`stellar_evolution.youngest_progenitor_sol`, the lifetime the star model gives it, not the main sequence's alone),
    then the white dwarf of the initial-final relation below the white dwarf progenitors' top, and above it a black hole
    of the kind's mean mass where the progenitor lies in a failure window (`remnant_progenitors`) and a neutron star of
    the kind's mean mass elsewhere.
    """
    return _present_mass(_remnant_states(np.asarray(initial_mass_sol, dtype=float)), _turn_off_sol(age_gy))


@functools.lru_cache(maxsize=1)
def _placed_bins():
    """
    The log grid of placed-star masses: `(number, mean birth mass, other stars)` arrays over `MASS_BINS` bins of the
    IMF between its edges (and the edges of the star-count classes, so a class's count never straddles a bin): the IMF's
    number in each, the mean birth mass of those stars, and the stars a system of that mass holds besides the placed one.
    """
    class_edges = [lower for _letter, lower in multistar.SPECTRAL_CLASS_LOWER_MASS_SOL]
    class_edges += [lower for _name, lower in tuning.B_SUBCLASS_LOWER_MASS_SOL]
    low, top = tuning.IMF_BREAKS_SOL[0], tuning.IMF_BREAKS_SOL[-1]
    edges = np.concatenate([np.geomspace(low, top, MASS_BINS + 1), tuning.IMF_BREAKS_SOL, class_edges])
    edges = np.unique(edges[(edges >= low) & (edges <= top)])
    number = np.array([imf_integral(0.0, a, b) for a, b in zip(edges[:-1], edges[1:])])
    first = np.array([imf_integral(1.0, a, b) for a, b in zip(edges[:-1], edges[1:])])
    mean_mass = first / number
    others = np.array([multistar.expected_star_count(m) - 1.0 for m in mean_mass])
    return number, mean_mass, others


def companion_fractions(birth_mass_sol):
    """`(len(birth_mass_sol), COMPANION_MASS_STEPS)` midpoints of the uniform draw of a further star's birth mass over the
    placed star's, for each placed star of `birth_mass_sol` (`multistar.companion_ratio_range`)."""
    ranges = np.array([multistar.companion_ratio_range(float(m)) for m in np.atleast_1d(birth_mass_sol)])
    steps = (np.arange(COMPANION_MASS_STEPS) + 0.5) / COMPANION_MASS_STEPS
    return ranges[:, :1] + (ranges[:, 1:] - ranges[:, :1]) * steps[None, :]


@functools.lru_cache(maxsize=1)
def _placed_states():
    """`(number, others, placed states, companion states)`: `_placed_bins`' number and other stars, and the
    `_remnant_states` of the placed stars' birth masses and of their further stars' (`companion_fractions`), which no age
    changes."""
    number, mean_mass, others = _placed_bins()
    companions = mean_mass[:, None] * companion_fractions(mean_mass)
    return number, others, _remnant_states(mean_mass), _remnant_states(companions)


@functools.lru_cache(maxsize=65536)
def mean_system_mass_at_age_solar(age_gy):
    """The mean mass of a system whose stars formed `age_gy` ago, solar masses: the placed star's (remnants included)
    and the other stars of its system (see the module docstring)."""
    age_gy = max(float(age_gy), 1.0e-4)
    number, others, placed_states, companion_states = _placed_states()
    turn_off = _turn_off_sol(age_gy)
    placed = float(np.sum(number * _present_mass(placed_states, turn_off)) / np.sum(number))
    further = _present_mass(companion_states, turn_off).mean(axis=1)
    return placed + float(np.sum(number * others * further) / np.sum(number))


def cohort_mean_system_mass_solar(ages_gy, weights=None):
    """
    The mean mass of a system of a population whose stars formed at the given ages, solar masses, each age weighted by
    its share of the star formation (equal weights when none are given).

    Args:
        ages_gy (sequence): Ages, Gy.
        weights (sequence, optional): The number of systems formed at each age, relative.
    """
    ages = [float(age) for age in ages_gy]
    weights = [1.0] * len(ages) if weights is None else [float(w) for w in weights]
    total = sum(weights)
    if total <= 0.0:
        raise ValueError("a cohort needs at least one age with a weight")
    return sum(w * mean_system_mass_at_age_solar(round(age, 6)) for age, w in zip(ages, weights)) / total


@functools.lru_cache(maxsize=256)
def window_mean_system_mass_solar(age_low_gy, age_high_gy):
    """`mean_system_mass_at_age_solar` averaged over a window of ages in which the stars formed at a steady rate."""
    low, high = sorted((max(age_low_gy, 1.0e-4), max(age_high_gy, 1.0e-4)))
    return cohort_mean_system_mass_solar(low + (high - low) * (np.arange(COHORT_AGES) + 0.5) / COHORT_AGES)
