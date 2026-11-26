"""
GEN.708: the double white dwarf merger products follow the share of their mass (`planetgen/physics/white_dwarf_merger.py`, rows WM1 to WM3
of docs/design/research/gap-filling-plan.md 6.330), tested against the cited samples: [R5355] by mass bin and [R5354] for the whole
old population. No database or Redis.
"""

import numpy as np
import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.physics import imf, ms_lifetime, remnant_mass, white_dwarf_merger


@pytest.mark.parametrize("mass, share", [(0.5, 0.0), (0.69, 0.0), (0.7, 0.25), (0.89, 0.25), (0.9, 0.25), (1.0, 0.36), (1.1, 0.78),
                                         (1.2, 0.49), (1.3, 0.38), (1.35, 0.38), (1.4, 0.38)])
def test_the_share_of_a_mass_is_the_bin_it_is_in(mass, share):
    assert white_dwarf_merger.merger_fraction(mass) == share


def test_the_bins_of_the_sample_are_contiguous_and_end_at_the_top_of_the_mass_range():
    bins = tuning.WHITE_DWARF_MERGER_FRACTION_BY_MASS
    assert [high for _low, high, _share in bins[:-1]] == [low for low, _high, _share in bins[1:]]
    assert bins[0][0] == 0.9 and bins[-1][1] == tuning.WD_MASS_RANGE_SOL[1]
    assert tuning.WHITE_DWARF_MERGER_FLOOR_SOL < bins[0][0]


def _old_population_share(age_gy):
    """The share of the white dwarfs of an old population that the mass law flags: progenitors from the turn-off up, by the initial
    mass function, each with the mass the initial-final relation gives."""
    low, high = tuning.WHITE_DWARF_PROGENITOR_SOL
    grid = np.geomspace(max(low, ms_lifetime.mass_for_lifetime_gy(age_gy)), high, 600)
    middle = 0.5 * (grid[:-1] + grid[1:])
    weight = np.array([imf.integral(0.0, a, b) for a, b in zip(grid[:-1], grid[1:])])
    flagged = np.array([white_dwarf_merger.merger_fraction(remnant_mass.white_dwarf_mass_sol(m)) for m in middle])
    return float(np.sum(weight * flagged) / np.sum(weight))


@pytest.mark.parametrize("age_gy", [3.0, 8.0, 12.0])
def test_an_old_population_is_a_merger_product_share_inside_the_factor_of_two_of_the_100_pc_sample(age_gy):
    """[R5354]: about 14 percent of the single white dwarfs of the 100 pc sample formed by a merger, a factor of 2 either way."""
    share = _old_population_share(age_gy)
    assert 0.07 <= share <= 0.28


def test_the_heavy_white_dwarfs_of_a_population_are_mostly_mergers_as_the_100_pc_sample_has_them():
    """[R5355]: from 0.9 to 1.3 solar masses the merger products are 49 percent of the 204 white dwarfs (99 of 204); the generator's
    heavy white dwarfs (the same bins weighted by what the relation gives an old population) come out near it."""
    low, high = tuning.WHITE_DWARF_PROGENITOR_SOL
    grid = np.geomspace(max(low, ms_lifetime.mass_for_lifetime_gy(8.0)), high, 600)
    middle = 0.5 * (grid[:-1] + grid[1:])
    weight = np.array([imf.integral(0.0, a, b) for a, b in zip(grid[:-1], grid[1:])])
    mass = np.array([remnant_mass.white_dwarf_mass_sol(m) for m in middle])
    heavy = mass >= 0.9
    flagged = np.array([white_dwarf_merger.merger_fraction(m) for m in mass])
    share = float(np.sum(weight[heavy] * flagged[heavy]) / np.sum(weight[heavy]))
    assert 0.3 <= share <= 0.65


def test_a_white_dwarf_primary_of_a_system_is_flagged_at_the_share_of_its_mass_not_at_the_blue_stragglers():
    assert multistar.draw_white_dwarf_merger(1.15, 8.0e9, 0.77)
    assert not multistar.draw_white_dwarf_merger(0.6, 8.0e9, 0.0)
