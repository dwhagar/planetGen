# tests/test_remnant_mass.py

"""GEN.562: the remnant mass rules live in `physics/remnant_mass.py`. Each is checked against the line it replaced."""

import numpy as np
import pytest

from planetgen import tuning
from planetgen.physics import remnant_mass
from planetgen.util import draw


def _original_white_dwarf_mass(initial):
    return min(max(tuning.WD_IFMR_SLOPE * initial + tuning.WD_IFMR_INTERCEPT_SOL, tuning.WD_MASS_RANGE_SOL[0]),
               tuning.WD_MASS_RANGE_SOL[1])


def test_the_white_dwarf_mass_is_the_initial_final_relation_held_in_its_range():
    for initial in np.geomspace(0.5, 12.0, 400):
        assert remnant_mass.white_dwarf_mass_sol(float(initial)) == _original_white_dwarf_mass(float(initial))
    assert remnant_mass.white_dwarf_mass_sol(0.0) == tuning.WD_MASS_RANGE_SOL[0]
    assert remnant_mass.white_dwarf_mass_sol(100.0) == tuning.WD_MASS_RANGE_SOL[1]


def test_the_white_dwarf_mass_of_an_array_is_the_mass_of_each_member():
    grid = np.geomspace(0.5, 12.0, 50)
    masses = remnant_mass.white_dwarf_mass_sol(grid)
    assert isinstance(masses, np.ndarray)
    assert masses == pytest.approx([_original_white_dwarf_mass(float(m)) for m in grid], rel=0, abs=0)


@pytest.mark.parametrize("kind, law", [("neutron-star", tuning.NEUTRON_STAR_MASS_RANGE_SOLAR),
                                       ("black-hole", tuning.BLACK_HOLE_MASS_RANGE_SOLAR)])
def test_a_remnant_mass_is_one_uniform_draw_over_the_kinds_range(kind, law):
    assert remnant_mass.mass_range_sol(kind) == law
    for seed in range(200):
        assert remnant_mass.draw_mass_sol(kind, rng=draw.Stream(seed)) == draw.Stream(seed).uniform(*law)
    low, high = law[0] + 0.1 * (law[1] - law[0]), law[0] + 0.2 * (law[1] - law[0])
    with draw.bound(5):
        masses = [remnant_mass.draw_mass_sol(kind, (low, high)) for _ in range(100)]
    assert all(low <= m <= high for m in masses)


def test_a_range_outside_the_law_is_refused_and_an_unknown_kind_too():
    with pytest.raises(ValueError):
        remnant_mass.within((1.0, 2.0), (3.0, 4.0))
    with pytest.raises(ValueError):
        remnant_mass.mass_range_sol("white-dwarf")
    assert remnant_mass.within((1.0, 2.0), None) == (1.0, 2.0)


def test_the_moments_are_those_of_the_uniform_draw():
    for kind in remnant_mass.KINDS:
        low, high = remnant_mass.mass_range_sol(kind)
        mean, square = remnant_mass.uniform_moments(kind)
        assert mean == pytest.approx((low + high) / 2)
        assert square == pytest.approx((high ** 3 - low ** 3) / (3 * (high - low)))
