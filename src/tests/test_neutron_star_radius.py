# tests/test_neutron_star_radius.py

"""
GEN.463, GEN.969: the radius of a neutron star from its mass (`physics/neutron_star.py`): a weighted quantile over the radii of the
21 equations of state of `neutron_star_structure`, weighted by the [R3075] constraints. The measurements it must
contain, that a star's radius lies between the lowest and highest equation of state at its mass, monotonicity in the draw,
and determinism.
"""

import statistics

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import NeutronStar
from planetgen.physics import constants, neutron_star, neutron_star_structure
from planetgen.util import draw


def test_the_weights_favour_the_equations_that_fit_the_constraints():
    weights = {eos: neutron_star.eos_weight(eos) for eos in range(1, 22)}
    assert all(0.0 < weight <= 1.0 for weight in weights.values())
    assert weights[1] < 0.05  # R(1.4) 13.8 km, 2 sigma above Miller et al.
    assert weights[8] > 0.8  # R(1.4) 12.7 km, inside both
    assert weights[6] < weights[8]  # holds only 2.01 solar masses, under PSR J0740+6620


def test_the_centre_is_the_weighted_median_and_lies_in_the_miller_range_with_the_tables_stiffness():
    centre = neutron_star.radius_km(1.4)
    assert 12.45 - 0.65 < centre < 12.45 + 0.65 + 0.5  # the 21 tables are stiffer than the 68% interval's centre
    assert neutron_star.radius_km(1.4, 0.0) == centre


@pytest.mark.parametrize("mass, radius", [
    (1.34, 12.71),  # PSR J0030+0451, [R3766]
    (2.08, 13.7),  # PSR J0740+6620, [R3075]
])
def test_the_measured_radii_lie_inside_the_limits(mass, radius):
    lowest, highest = neutron_star.radius_limits_km(mass)
    assert lowest < radius < highest


def test_the_limits_over_the_mass_range_are_the_equations_of_state():
    for mass in tuning.NEUTRON_STAR_MASS_RANGE_SOLAR + (1.4, 1.8):
        lowest, highest = neutron_star.radius_limits_km(mass)
        band = neutron_star_structure.radius_band_km(mass)
        assert (lowest, highest) == (band[0], band[2])
        assert 11.0 < lowest < 12.0 and 13.5 < highest < 14.2


def test_the_radius_grows_with_the_draw_and_stops_at_the_limits():
    for mass in (1.2, 1.4, 2.0):
        radii = [neutron_star.radius_km(mass, deviation / 2.0) for deviation in range(-8, 9)]
        assert radii == sorted(radii)
        assert neutron_star.radius_km(mass, 50.0) == neutron_star.radius_limits_km(mass)[1]
        assert neutron_star.radius_km(mass, -50.0) == neutron_star.radius_limits_km(mass)[0]


def test_a_mass_no_table_holds_is_refused():
    with pytest.raises(ValueError):
        neutron_star.radius_km(3.5)


def test_a_drawn_star_is_between_the_limits_at_its_mass_and_its_radii_spread():
    draw.set_run_seed(463)
    stars = [NeutronStar(SystemConfig()) for _ in range(300)]
    for star in stars:
        low, high = neutron_star.radius_limits_km(star.mass_solar)
        assert low - 1e-9 <= star.radius <= high + 1e-9
    radii = [star.radius for star in stars]
    assert 11.8 < statistics.mean(radii) < 13.2
    assert 0.2 < statistics.pstdev(radii) < 0.9
    assert min(radii) < 12.3 and max(radii) > 13.0


def test_a_neutron_star_companion_has_the_centre_radius_of_its_mass():
    leaf = multistar._remnant_leaf("neutron-star", 15.0, 1.4)
    assert leaf.radius_au * constants.AU_M / 1000.0 == pytest.approx(neutron_star.radius_km(1.4))


def test_the_same_seed_gives_the_same_radius():
    def radii():
        draw.set_run_seed(7)
        return [NeutronStar(SystemConfig()).radius for _ in range(5)]
    assert radii() == radii()
