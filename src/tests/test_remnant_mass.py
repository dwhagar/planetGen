# tests/test_remnant_mass.py

"""GEN.562: the remnant mass rules live in `physics/remnant_mass.py`. Each is checked against the line it replaced."""

import numpy as np
import pytest

from planetgen import tuning
from planetgen.physics import remnant_mass
from planetgen.util import draw


def _cummings_2018_mist(initial):
    """[R0961], equations 4, 5 and 6, written out here: the end pieces carried beyond their ranges, held inside the range."""
    if initial < 2.85:
        mass = 0.080 * initial + 0.489
    elif initial < 3.60:
        mass = 0.187 * initial + 0.184
    else:
        mass = 0.107 * initial + 0.471
    return min(max(mass, tuning.WD_MASS_RANGE_SOL[0]), tuning.WD_MASS_RANGE_SOL[1])


def test_the_white_dwarf_mass_is_the_cummings_relation_held_in_its_range():
    for initial in np.geomspace(0.5, 12.0, 400):
        assert remnant_mass.white_dwarf_mass_sol(float(initial)) == pytest.approx(_cummings_2018_mist(float(initial)), abs=2.0e-3)
    assert remnant_mass.white_dwarf_mass_sol(0.0) == tuning.WD_MASS_RANGE_SOL[0]
    assert remnant_mass.white_dwarf_mass_sol(100.0) == tuning.WD_MASS_RANGE_SOL[1]


def test_the_relation_is_continuous_at_its_joins_and_rises_with_the_initial_mass():
    for join in (2.85, 3.60):
        assert remnant_mass.white_dwarf_mass_sol(join - 1.0e-9) == pytest.approx(remnant_mass.white_dwarf_mass_sol(join + 1.0e-9), abs=1.0e-6)
    masses = remnant_mass.white_dwarf_mass_sol(np.linspace(0.8, 8.0, 200))
    assert np.all(np.diff(masses) > 0.0)
    assert remnant_mass.white_dwarf_mass_sol(8.0) < 1.44    # under the Chandrasekhar limit at the top of the progenitors


def test_the_progenitor_is_the_relation_inverted():
    for initial in np.linspace(0.8, 8.0, 60):
        assert remnant_mass.progenitor_mass_sol(remnant_mass.white_dwarf_mass_sol(float(initial))) == pytest.approx(initial, abs=1.0e-9)
    grid = np.linspace(0.8, 8.0, 50)
    assert remnant_mass.progenitor_mass_sol(remnant_mass.white_dwarf_mass_sol(grid)) == pytest.approx(grid, abs=1.0e-9)


def test_the_white_dwarf_mass_of_an_array_is_the_mass_of_each_member():
    grid = np.geomspace(0.5, 12.0, 50)
    masses = remnant_mass.white_dwarf_mass_sol(grid)
    assert isinstance(masses, np.ndarray)
    assert masses == pytest.approx([remnant_mass.white_dwarf_mass_sol(float(m)) for m in grid], rel=0, abs=0)


def test_a_black_hole_mass_is_one_uniform_draw_over_the_range():
    law = tuning.BLACK_HOLE_MASS_RANGE_SOLAR
    assert remnant_mass.mass_range_sol("black-hole") == law
    for seed in range(200):
        assert remnant_mass.draw_mass_sol("black-hole", rng=draw.Stream(seed)) == draw.Stream(seed).uniform(*law)
    low, high = law[0] + 0.1 * (law[1] - law[0]), law[0] + 0.2 * (law[1] - law[0])
    with draw.bound(5):
        masses = [remnant_mass.draw_mass_sol("black-hole", (low, high)) for _ in range(100)]
    assert all(low <= m <= high for m in masses)


def test_a_neutron_star_mass_is_the_slow_pulsar_gaussian_held_to_its_range():
    """[R3396]: slow pulsars 1.49 +- 0.19 (`tuning.STELLAR_NS_MASS_GAUSSIANS`); few are above 1.9 solar masses."""
    assert remnant_mass.mass_range_sol("neutron-star") == tuning.NEUTRON_STAR_MASS_RANGE_SOLAR
    low, high = tuning.NEUTRON_STAR_MASS_RANGE_SOLAR
    with draw.bound(3):
        masses = np.array([remnant_mass.draw_mass_sol("neutron-star") for _ in range(20000)])
    assert masses.min() >= low and masses.max() <= high
    assert masses.mean() == pytest.approx(1.49, abs=0.01)
    assert masses.std() == pytest.approx(0.18, abs=0.01)
    assert (masses > 1.9).mean() == pytest.approx(remnant_mass.share_above("neutron-star", 1.9), abs=0.004)
    assert 0.01 < remnant_mass.share_above("neutron-star", 1.9) < 0.025      # the 1.5% of the published Gaussian, not the uniform's 27%
    with draw.bound(5):
        cut = [remnant_mass.draw_mass_sol("neutron-star", (1.6, 1.8)) for _ in range(200)]
    assert all(1.6 <= m <= 1.8 for m in cut)


def test_the_share_above_a_mass_runs_from_all_to_none_across_the_range():
    for kind in remnant_mass.KINDS:
        low, high = remnant_mass.mass_range_sol(kind)
        assert remnant_mass.share_above(kind, low - 1.0) == 1.0 and remnant_mass.share_above(kind, high + 1.0) == 0.0
        shares = [remnant_mass.share_above(kind, low + (high - low) * f) for f in (0.1, 0.3, 0.5, 0.7, 0.9)]
        assert shares == sorted(shares, reverse=True)


def test_a_range_outside_the_law_is_refused_and_an_unknown_kind_too():
    with pytest.raises(ValueError):
        remnant_mass.within((1.0, 2.0), (3.0, 4.0))
    with pytest.raises(ValueError):
        remnant_mass.mass_range_sol("white-dwarf")
    assert remnant_mass.within((1.0, 2.0), None) == (1.0, 2.0)


def test_the_moments_are_those_of_the_draw():
    low, high = remnant_mass.mass_range_sol("black-hole")
    mean, square = remnant_mass.moments("black-hole")
    assert mean == pytest.approx((low + high) / 2)
    assert square == pytest.approx((high ** 3 - low ** 3) / (3 * (high - low)))
    with draw.bound(8):
        masses = np.array([remnant_mass.draw_mass_sol("neutron-star") for _ in range(40000)])
    mean, square = remnant_mass.moments("neutron-star")
    assert mean == pytest.approx(masses.mean(), abs=0.01)
    assert square == pytest.approx((masses ** 2).mean(), abs=0.02)
