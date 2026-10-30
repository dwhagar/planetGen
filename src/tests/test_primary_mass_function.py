# tests/test_primary_mass_function.py

"""
GEN.928, GEN.929: the primaries' mass function is the single-star function less the companions of heavier primaries, and the
systems per solar mass are one number the placed systems' mean mass agrees with.
"""

import numpy as np
import pytest

from planetgen.generation import multistar, primary_mass_function as pmf, system_mass
from planetgen.physics import imf

MASSES = (0.1, 0.2, 0.5, 0.8, 1.0, 2.0, 5.0, 20.0, 100.0)


def test_the_primary_function_is_positive_and_below_the_single_star_function_at_every_node():
    nodes = np.array(pmf._nodes())
    primary = imf.density(nodes, curve=imf.PRIMARY)
    single = imf.density(nodes, curve=imf.SINGLE)
    assert np.all(primary > 0.0) and np.all(primary <= single * (1.0 + 1.0e-9))


def test_primaries_and_the_companions_of_heavier_primaries_add_up_to_the_stars_born():
    """The defining equation, checked on a window not used to solve it: stars born in the window = primaries in it + companions."""
    for low, high in ((0.1, 0.3), (0.3, 0.8), (0.8, 2.0), (2.0, 8.0), (8.0, 40.0)):
        born = imf.number(low, high, imf.SINGLE)
        primaries = imf.number(low, high, imf.PRIMARY)
        companions = 0.0
        edges = pmf._edges()
        mid = np.sqrt(edges[:-1] * edges[1:])
        _, cells = pmf._solve()
        for m, count in zip(mid, cells):
            window_low, window_high = pmf.companion_window(m)
            overlap = max(0.0, min(high, window_high) - max(low, window_low))
            width = window_high - window_low
            if width > 1.0e-12:
                companions += count * pmf.companions_per_primary(m) * overlap / width
        assert primaries + companions == pytest.approx(born, rel=0.01), (low, high)


def test_a_companion_is_never_heavier_than_its_primary_so_the_top_mass_is_all_primaries():
    top = imf.highest_mass_sol()
    assert imf.density(np.array([top * 0.999]), curve=imf.PRIMARY)[0] == pytest.approx(
        imf.density(np.array([top * 0.999]), curve=imf.SINGLE)[0], rel=0.02)


def test_more_systems_than_the_stars_in_them_is_impossible_and_one_per_primary_is_the_count():
    per_star = pmf.systems_per_star()
    assert 0.5 < per_star < 1.0
    stars_per_system = 1.0 / per_star
    assert stars_per_system == pytest.approx(sum(
        pmf.companions_per_primary(m) * c for m, c in zip(
            np.sqrt(pmf._edges()[:-1] * pmf._edges()[1:]), pmf._solve()[1])) / pmf._solve()[1].sum() + 1.0, rel=0.02)


def test_the_solution_does_not_move_when_the_cells_are_twice_as_wide(monkeypatch):
    fine = pmf.density_at(np.array(MASSES))
    monkeypatch.setattr(pmf, "GRID_STEP_DEX", pmf.GRID_STEP_DEX * 2.0)
    pmf._solve.cache_clear()
    try:
        coarse = pmf.density_at(np.array(MASSES))
    finally:
        monkeypatch.undo()
        pmf._solve.cache_clear()
    assert coarse == pytest.approx(fine, rel=0.02)


def test_the_systems_per_solar_mass_is_the_mean_system_mass_of_the_construction_at_birth():
    """GEN.929: one number, not a second anchor. The placed systems' mean mass at age 0 is 1 over the systems per solar mass."""
    assert system_mass.mean_system_mass_at_age_solar(0.0) == pytest.approx(1.0 / pmf.systems_per_solar_mass(), rel=0.005)


def test_the_companion_window_is_the_one_multistar_draws_from():
    for mass in MASSES:
        low, high = pmf.companion_window(mass)
        ratio_low, ratio_high = multistar.companion_ratio_range(mass)
        assert (low, high) == pytest.approx((ratio_low * mass, ratio_high * mass))
