"""GEN.684 part A, GEN.462: the Bond albedo of a dry body's ground and its space weathering (`planetgen/physics/ground_albedo.py`)."""

import types

import pytest

from planetgen.physics import atmosphere_colour, ground_albedo as ga, surface_composition
from planetgen.util import draw


def _body(carbon_rock=False, metal=0.0, age=4.5):
    return types.SimpleNamespace(carbon_rock=carbon_rock, metal_mass_fraction=metal, star=types.SimpleNamespace(age=age))


def test_the_dry_ground_is_carbon_metal_or_basalt():
    assert ga.dry_ground_type(_body()) == ga.BASALTIC
    assert ga.dry_ground_type(_body(carbon_rock=True)) == ga.CARBON
    assert ga.dry_ground_type(_body(metal=surface_composition.RANGES["metal_stripped"])) == ga.METAL
    assert ga.dry_ground_type(_body(metal=0.3)) == ga.BASALTIC


def test_a_basaltic_powder_has_the_mean_of_its_table_row_and_carbon_is_dark():
    """The table's basalt powder is (0.207, 0.214, 0.224) in the three visible bands."""
    assert ga.fresh_bond_albedo(_body()) == pytest.approx(ga.VISIBLE_TO_BOND * (0.207 + 0.214 + 0.224) / 3.0)
    assert ga.fresh_bond_albedo(_body(carbon_rock=True)) < 0.05


def test_a_young_surface_is_a_slab_and_the_metal_takes_the_stand_in():
    old, young = ga.fresh_bond_albedo(_body(age=4.5)), ga.fresh_bond_albedo(_body(age=0.01))
    assert old != young and 0.0 < young < 1.0
    assert ga.fresh_bond_albedo(_body(metal=0.9)) == atmosphere_colour.STAND_IN_GROUND_ALBEDO


def test_weathering_darkens_with_age_and_keeps_a_fresh_surface():
    assert ga.weathering_factor(0.0, 1.5, 0.5) == pytest.approx(1.0)
    assert ga.weathering_factor(1e6, 1.5, 0.5) == pytest.approx(0.5)
    assert ga.weathering_factor(4.0, 1.5, 0.5) < ga.weathering_factor(1.0, 1.5, 0.5) < 1.0
    with pytest.raises(ValueError):
        ga.weathering_factor(-1.0, 1.0, 0.5)
    with pytest.raises(ValueError):
        ga.weathering_factor(1.0, 0.0, 0.5)


def test_an_old_basaltic_surface_is_as_dark_as_the_moon_and_mercury():
    """Moon 0.11, Mercury 0.088 (the measured Bond albedos): the weathered fresh basalt over the ranges of the draw sits between 0.06 and 0.17 at 4 Gyr."""
    fresh = ga.fresh_bond_albedo(_body())
    with draw.bound(5):
        old = [ga.weathered_albedo(fresh, 4.0) for _ in range(300)]
    assert 0.06 < min(old) and max(old) < 0.17
    assert 0.09 < sum(old) / len(old) < 0.14


def test_the_weathering_draw_stays_in_its_ranges_and_repeats():
    with draw.bound(1):
        values = [ga.draw_weathering() for _ in range(100)]
    assert all(1.0 <= tau <= 2.0 and 0.3 <= d <= 0.7 for tau, d in values)
    with draw.bound(1):
        assert [ga.draw_weathering() for _ in range(100)] == values
