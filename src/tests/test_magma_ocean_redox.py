# tests/test_magma_ocean_redox.py

"""GEN.1005: the surface redox of a magma ocean from the pressure at its base (`physics/magma_ocean_redox.py`). The Moon, Mars and the Earth are checks, never inputs."""

import math

import pytest

from planetgen.physics import magma_ocean_redox as mor

MOON, MARS, EARTH = 0.0123, 0.107, 1.0

# What the sources print ([R1951] section 6.3 and Table 1; [R1075]), the bodies as test points only.
PUBLISHED_SURFACE_DIW = {MOON: (-2.2, -1.5), MARS: (-1.4, 0.2), EARTH: (-1.2, 1.8)}
PUBLISHED_BASE_GPA = {MOON: 5.0, MARS: 14.0, EARTH: 55.0}


def test_the_base_pressure_of_three_bodies_is_near_the_ones_deng_et_al_compute():
    for mass, published in PUBLISHED_BASE_GPA.items():
        assert mor.base_pressure_gpa(mass) == pytest.approx(published, rel=0.2)


def test_the_base_pressure_stops_at_the_end_of_the_calculation():
    assert mor.base_pressure_gpa(5.0) < mor.CALCULATION_LIMIT_GPA
    assert mor.base_pressure_gpa(6.0) == mor.CALCULATION_LIMIT_GPA == mor.base_pressure_gpa(50.0)


def test_the_relation_reproduces_the_pairs_of_hirschmann_2022_for_the_earth():
    assert mor.diw_from_fe3_fraction(0.034) == pytest.approx(-0.2, abs=0.1)
    assert mor.diw_from_fe3_fraction(0.103) == pytest.approx(1.8, abs=0.1)


def test_four_decades_of_oxygen_fugacity_per_decade_of_ferric_ratio():
    ratio = 0.05
    fraction = ratio / (1.0 + ratio)
    assert mor.diw_from_fe3_fraction(fraction * 10.0 / (1.0 + 9.0 * fraction)) - mor.diw_from_fe3_fraction(fraction) == pytest.approx(4.0)


@pytest.mark.parametrize("mass", sorted(PUBLISHED_SURFACE_DIW))
def test_the_model_lands_inside_the_published_range_of_each_body_within_one_sigma(mass):
    low, high = PUBLISHED_SURFACE_DIW[mass]
    mean, sigma = mor.surface_diw(mass), mor.diw_sigma(mass)
    assert low - sigma <= mean <= high + sigma
    assert abs(mean - (low + high) / 2) <= (high - low) / 2 + 0.2


def test_a_heavier_body_has_a_more_oxidised_surface_up_to_the_end_of_the_calculation_and_not_beyond():
    values = [mor.surface_diw(mass) for mass in (0.01, 0.1, 1.0, 3.0, 5.0)]
    assert values == sorted(values) and len(set(values)) == len(values)
    assert mor.surface_diw(6.0) == mor.surface_diw(30.0)


def test_the_exponent_is_the_middle_of_its_bracket_and_both_ends_of_the_bracket_move_the_surface_the_same_way_either_side():
    low, high = mor.FE3_EXPONENT_BRACKET
    assert low < mor.FE3_PRESSURE_EXPONENT < high
    # a small body: a steeper power means less ferric iron, so a lower dIW; a heavy body the other way round
    assert mor.diw_from_fe3_fraction(mor.fe3_fraction(MOON, exponent=high)) < mor.diw_from_fe3_fraction(mor.fe3_fraction(MOON, exponent=low))
    assert mor.diw_from_fe3_fraction(mor.fe3_fraction(3.0, exponent=high)) > mor.diw_from_fe3_fraction(mor.fe3_fraction(3.0, exponent=low))


def test_the_width_is_the_quadrature_of_the_two_log_spreads_and_least_where_the_pressure_is_the_reference():
    low, high = mor.FE3_FRACTION_RANGE
    share = (math.log10(high / (1 - high)) - math.log10(low / (1 - low))) / math.sqrt(12.0)
    assert mor.diw_sigma(1.0 * (55.0 / mor.BASE_PRESSURE_1_EARTH_GPA) ** (1 / mor.BASE_PRESSURE_EXPONENT)) == pytest.approx(4.0 * share)
    assert mor.diw_sigma(MOON) > mor.diw_sigma(1.0)
    assert mor.diw_sigma(0.001) > mor.diw_sigma(MOON)


def test_the_later_mantle_is_above_the_ocean_surface_and_held_at_one_earth_mass():
    for mass in (0.001, MOON, MARS, 1.0):
        assert mor.later_mantle_diw(mass) > mor.surface_diw(mass)
    assert mor.later_mantle_diw(1.0) - mor.surface_diw(1.0) == pytest.approx(mor.OXIDATION_INCREMENT[1])
    assert mor.later_mantle_diw(0.05) - mor.surface_diw(0.05) == pytest.approx(mor.OXIDATION_INCREMENT[0])
    assert mor.later_mantle_diw(4.0) == mor.later_mantle_diw(1.0)


@pytest.mark.parametrize("mass,published", [(MOON, (-2.0, 0.0)), (MARS, (0.0, 1.0)), (EARTH, (2.5, 4.5))])
def test_the_later_mantle_is_inside_the_published_ranges_within_one_sigma(mass, published):
    low, high = published
    assert low - mor.LATER_MANTLE_SIGMA <= mor.later_mantle_diw(mass) <= high + mor.LATER_MANTLE_SIGMA


def test_a_body_without_a_mass_or_a_ferric_share_is_a_refusal():
    with pytest.raises(ValueError):
        mor.base_pressure_gpa(0.0)
    with pytest.raises(ValueError):
        mor.diw_from_fe3_fraction(1.0)


def test_the_recipe_expressions_can_call_the_model():
    """The rows a1 and a2 read the model through the expression language (GEN.1005 part 2): the same numbers as the module."""
    from planetgen.physics import atmosphere_recipes as recipes

    class Env:
        def get(self, name):
            return {"mass_earth": 0.3}[name]

    for text, expected in (("magma_ocean_diw(mass_earth)", mor.surface_diw(0.3)), ("magma_ocean_diw_sigma(mass_earth)", mor.diw_sigma(0.3)),
                           ("later_mantle_diw(mass_earth)", mor.later_mantle_diw(0.3))):
        assert recipes.Expression(text, "test")(Env()) == pytest.approx(expected)
