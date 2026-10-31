# generation/companion_classes.py

"""
Which kind of companion a star gets and how far away it is (GEN.448; docs/design/systems/multistar-and-compact-systems.md
section 4 and placement-by-population.md 3.1): one distribution of companions over orbital period for the whole range,
in place of the coin flip between a close pair (0.05 to 0.25 AU) and a wide pair (50 AU and out) with nothing between.

The distribution is Moe and Di Stefano's (2017, Table 13: `tuning.COMPANIONS_PER_DECADE_BY_PRIMARY_MASS`), the
companions per decade of period at log P 1, 3, 5 and 7 for five primary masses, joined by straight lines in log P and
log mass, held flat past the ends. Its share between three separation ranges is the share of a companion's chances:

- the close range, `CLOSE_RANGE_AU` (0.05 to 0.25 AU), where the pair is a P-type pair (`generation.binary`);
- the inner range, 0.25 AU to `INNER_LIMIT_AU` (10 AU), an S-type pair that the metallicity trend of GEN.103 reaches
  (Moe, Kratter and Badenes 2019 measured it for periods under 10^4 days, which is 10 AU for a pair of two solar
  masses: the trend belongs to the close and inner ranges together);
- the outer range, from `INNER_LIMIT_AU` to the widest axis the place keeps (GEN.103), where the trend does not reach.

The odds factor of GEN.103 multiplies the odds of the close and inner ranges together; the outer range keeps its chance
unless the two would pass 1. A wide pair then draws its axis from the same distribution, cut to 0.25 AU and the place's
limit, with the inner part of the range carrying the factor's weight, so the trend shapes the separations too.
"""

import functools
import math
from dataclasses import dataclass

import numpy as np

from planetgen import tuning
from planetgen.generation import close_binary_trend, mass_ratio, multistar, place_terms
from planetgen.util import draw
from planetgen.physics import unit_registry as units

_yr_to_d = units.converter("yr", "d")

CLOSE_RANGE_AU = (0.05, 0.25)
"""tuple: The separations of a close (P-type) pair, `generation.binary.BinaryStarProxy`."""

INNER_LIMIT_AU = 10.0
"""float: Where the metallicity trend stops (Moe, Kratter and Badenes 2019: periods under 10^4 days, 10 AU)."""

_GRID_POINTS = 801
"""int: Points of the grid over log P that the distribution is integrated and inverted on."""

_DAYS_PER_YEAR = _yr_to_d(1.0)


def mean_mass_ratio(primary_born_mass_solar):
    """The mean mass ratio of a companion of a primary of this birth mass (`mass_ratio.mean_ratio`): the secondary the shares
    are worked out for, before the real one is drawn."""
    return mass_ratio.mean_ratio(primary_born_mass_solar)


def log_period_days(separation_au, total_mass_solar):
    """log10 of the period in days of a pair of this axis and total mass (Kepler's third law)."""
    return math.log10(math.sqrt(separation_au ** 3 / total_mass_solar) * _DAYS_PER_YEAR)


def separation_au(log_period, total_mass_solar):
    """The axis in AU of a pair of this total mass and log10 period in days."""
    years = 10.0 ** log_period / _DAYS_PER_YEAR
    return (years * years * total_mass_solar) ** (1.0 / 3.0)


@functools.lru_cache(maxsize=512)
def _anchors(born_mass_solar):
    """The four companions-per-decade values for a primary of this mass: straight in log mass between the table's
    bins, flat past the first and last."""
    table = tuning.COMPANIONS_PER_DECADE_BY_PRIMARY_MASS
    masses = [math.log10(mass) for mass, _values in table]
    columns = list(zip(*[values for _mass, values in table]))
    where = math.log10(max(born_mass_solar, 1e-6))
    return tuple(float(np.interp(where, masses, column)) for column in columns)


def companions_per_decade(log_period, born_mass_solar):
    """The companions per decade of period at `log_period` (array or number) for this primary mass."""
    return np.interp(log_period, tuning.COMPANION_LOG_PERIODS_DAYS, _anchors(round(born_mass_solar, 4)))


def _grid(low_au, high_au, total_mass_solar):
    """The log P grid between two axes."""
    return np.linspace(log_period_days(low_au, total_mass_solar), log_period_days(high_au, total_mass_solar),
                       _GRID_POINTS)


def _weight(born_mass_solar, total_mass_solar, low_au, high_au):
    """The expected companions in an axis range: the density integrated over log P (0 for an empty range)."""
    if high_au <= low_au:
        return 0.0
    grid = _grid(low_au, high_au, total_mass_solar)
    return float(np.trapezoid(companions_per_decade(grid, born_mass_solar), grid))


@dataclass(frozen=True)
class CompanionOdds:
    """
    What a star's companion roll uses.

    Attributes:
        close (float): Chance of a close (P-type) companion.
        wide (float): Chance of a wide (S-type) one, 0.25 AU and out.
        inner_boost (float): How much the metallicity factor multiplies the weight of the inner range (0.25 to 10 AU)
            among the wide pair's separations (1 with no [Fe/H] or at the line's zero).
        wide_share (float): The wide kinds' share of the companions before the trend.
        within_10_au (float): The chance of a companion inside `INNER_LIMIT_AU` (close or wide), what the trend of
            Moe, Kratter and Badenes 2019 measured.
    """
    close: float
    wide: float
    inner_boost: float
    wide_share: float
    within_10_au: float = 0.0


def odds(primary_born_mass_solar, chance, fe_h, limit_au, wide_share_scale=None):
    """
    The chances of a close and a wide companion for a star that has a `chance` of any, at its own [Fe/H], where the
    widest axis its place keeps is `limit_au`.

    Args:
        primary_born_mass_solar (float): The primary's birth mass.
        chance (float): The chance of a companion at all.
        fe_h (float or None): The primary's [Fe/H].
        limit_au (float): The widest axis the place keeps for this pair (`place_terms.wide_axis_limit_au`).
        wide_share_scale (callable or None): Maps the wide kinds' share of the companions (before the trend) to the
            share to use; the run's `wide_binary` prevalence (GEN.52). None leaves it.

    Returns:
        CompanionOdds: The chance of each kind and the inner range's weight among the wide separations.
    """
    total = primary_born_mass_solar * (1.0 + mean_mass_ratio(primary_born_mass_solar))
    birth_top = multistar.ROOT_MAX_AXIS_AU
    top = min(limit_au, tuning.WIDE_BINARY_SEPARATION_CEILING_AU, birth_top)
    outer_low = max(INNER_LIMIT_AU, CLOSE_RANGE_AU[1])
    close = _weight(primary_born_mass_solar, total, *CLOSE_RANGE_AU)
    inner = _weight(primary_born_mass_solar, total, CLOSE_RANGE_AU[1], INNER_LIMIT_AU)
    outer = _weight(primary_born_mass_solar, total, outer_low, birth_top)
    # What the place keeps of each wide range: the part of the birth distribution under its widest axis.
    kept_inner = _weight(primary_born_mass_solar, total, CLOSE_RANGE_AU[1], min(top, INNER_LIMIT_AU)) / inner
    kept_outer = _weight(primary_born_mass_solar, total, outer_low, top) / outer
    born = close + inner + outer
    close_share, inner_share, outer_share = close / born, inner / born, outer / born
    wide_share = inner_share + outer_share
    if wide_share_scale is not None and wide_share > 0.0:
        scaled = min(max(wide_share_scale(wide_share), 0.0), 1.0)
        inner_share, outer_share = inner_share * scaled / wide_share, outer_share * scaled / wide_share
        close_share, wide_share = 1.0 - scaled, scaled
    inside_share = close_share + inner_share
    p_inside = close_binary_trend.scaled_close_chance(chance * inside_share, fe_h) if inside_share > 0.0 else 0.0
    p_outside = min(chance * outer_share, 1.0 - p_inside)
    close_part = p_inside * close_share / inside_share if inside_share > 0.0 else 0.0
    inner_born = p_inside * inner_share / inside_share if inside_share > 0.0 else 0.0
    boost = inner_born / (chance * inner_share) if chance * inner_share > 0.0 else 1.0
    return CompanionOdds(close=close_part, wide=inner_born * kept_inner + p_outside * kept_outer, inner_boost=boost,
                         wide_share=wide_share, within_10_au=close_part + inner_born * kept_inner)


def sample_separation_au(born_mass_solar, total_mass_solar, limit_au, inner_boost=1.0, rng=draw):
    """
    A wide pair's axis in AU from the companion distribution cut to `[WIDE_BINARY_SEPARATION_MIN_AU, limit_au]`, the
    inner range (up to `INNER_LIMIT_AU`) carrying `inner_boost` times its weight. One uniform draw.

    Args:
        born_mass_solar (float): The primary's birth mass (it picks the distribution).
        total_mass_solar (float): The pair's total mass (it sets the period of an axis).
        limit_au (float): The widest axis the place keeps; a limit under the lower bound gives the lower bound.
        inner_boost (float): The weight factor of the inner range.
        rng: The draw source (default: the unit's stream).
    """
    pick = rng.random()
    low = tuning.WIDE_BINARY_SEPARATION_MIN_AU
    high = min(max(limit_au, low), tuning.WIDE_BINARY_SEPARATION_CEILING_AU)
    if high <= low * (1.0 + 1e-12):
        return low
    grid = _grid(low, high, total_mass_solar)
    density = companions_per_decade(grid, born_mass_solar)
    middle = 0.5 * (grid[:-1] + grid[1:])
    step = 0.5 * (density[:-1] + density[1:]) * np.diff(grid)
    knot = log_period_days(INNER_LIMIT_AU, total_mass_solar)
    step = np.where(middle < knot, step * inner_boost, step)
    cumulative = np.concatenate(([0.0], np.cumsum(step)))
    total = float(cumulative[-1])
    if total <= 0.0:
        return math.exp(math.log(low) + (math.log(high) - math.log(low)) * pick)
    return separation_au(float(np.interp(pick * total, cumulative, grid)), total_mass_solar)
