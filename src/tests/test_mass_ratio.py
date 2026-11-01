# tests/test_mass_ratio.py

"""
GEN.951: the mass ratio of a companion is drawn from the law Moe and Di Stefano (2017, ApJS 230:15, Eqn. 2 and Table 13) measured:
a power law by primary mass and period with a twin excess, not a flat one.
"""

import numpy as np
import pytest

from planetgen import tuning
from planetgen.generation import mass_ratio
from planetgen.util import draw

MASSES = (0.1, 0.3, 0.8, 1.0, 2.0, 3.16, 7.0, 12.0, 30.0, 100.0)


@pytest.mark.parametrize("mass", MASSES)
def test_the_law_is_one_distribution_over_the_ratios_a_companion_can_have(mass):
    ratios = mass_ratio.Ratios(mass)
    floor = float(ratios.floor[0])
    assert floor == pytest.approx(float(mass_ratio.floor_ratio(mass)))
    below = ratios.fraction_below(np.array([floor, 0.5 * (floor + 1.0), 1.0]))
    assert below[0] == pytest.approx(0.0, abs=1.0e-9) and below[2] == pytest.approx(1.0, abs=1.0e-9)
    grid = np.linspace(floor, 1.0, 200)
    assert np.all(np.diff(ratios.fraction_below(grid)) >= -1.0e-12)


@pytest.mark.parametrize("mass", MASSES)
def test_the_draw_follows_the_law_and_its_mean_is_the_means(mass):
    ratios = mass_ratio.Ratios(mass)
    rng = draw.Stream(f"gen951-law-{mass}")
    drawn = np.sort([mass_ratio.draw(mass, rng.random()) for _ in range(8000)])
    assert drawn[0] >= float(ratios.floor[0]) - 1.0e-9 and drawn[-1] <= 1.0
    grid = np.linspace(float(ratios.floor[0]), 1.0, 120)
    assert np.max(np.abs(np.searchsorted(drawn, grid, side="right") / len(drawn) - ratios.fraction_below(grid))) < 0.02
    assert float(np.mean(drawn)) == pytest.approx(mass_ratio.mean_ratio(mass), abs=0.012)


def test_the_quantiles_invert_the_fraction_below():
    ratios = mass_ratio.Ratios(np.array([0.5, 1.0, 5.0, 40.0]))
    units = (np.arange(16) + 0.5) / 16
    found = ratios.quantiles(units)
    for row in range(4):
        assert ratios.fraction_below(found[row], np.full(16, row)) == pytest.approx(units, abs=1.0e-6)


def test_a_solar_type_primary_is_near_flat_with_a_twin_excess_and_a_massive_one_leans_to_light_companions():
    """Table 13: solar-type slopes about -0.5 to +0.3 and a twin fraction 0.1 to 0.3 out to log P = 5; an O star's -1.7 to -2
    past log P = 3 and no twins there."""
    solar, massive = mass_ratio.Ratios(1.0), mass_ratio.Ratios(40.0)
    twins = lambda ratios: 1.0 - float(ratios.fraction_below(np.array([tuning.MASS_RATIO_TWIN_FROM]))[0])
    assert 0.10 < twins(solar) < 0.20 and twins(massive) < 0.06
    assert solar.mean() > massive.mean() + 0.1
    assert 0.45 < solar.mean() < 0.6 and 0.3 < massive.mean() < 0.45


def test_a_primary_under_the_hydrogen_burning_limit_has_a_companion_as_heavy_as_it_is():
    assert mass_ratio.draw(0.05, 0.37) == 1.0 and mass_ratio.mean_ratio(0.05) == 1.0


def test_the_mass_ratio_of_a_light_primary_is_cut_at_the_hydrogen_burning_limit_and_renormalised():
    ratios = mass_ratio.Ratios(0.2)
    assert float(ratios.floor[0]) == pytest.approx(0.4)
    assert ratios.fraction_below(np.array([0.4]))[0] == pytest.approx(0.0, abs=1.0e-9)
    assert 0.0 < ratios.fraction_below(np.array([0.7]))[0] < 1.0


def test_the_tables_have_a_row_per_primary_mass_and_a_column_per_period():
    rows, columns = len(tuning.COMPANIONS_PER_DECADE_BY_PRIMARY_MASS), len(tuning.COMPANION_LOG_PERIODS_DAYS)
    for table in (tuning.MASS_RATIO_TWIN_EXCESS, tuning.MASS_RATIO_LARGE_Q_SLOPE, tuning.MASS_RATIO_SMALL_Q_SLOPE):
        assert len(table) == rows and all(len(row) == columns for row in table)
    assert all(0.0 <= value < 1.0 for row in tuning.MASS_RATIO_TWIN_EXCESS for value in row)
