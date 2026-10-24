# tests/test_neutron_star_radius.py

"""
GEN.463: the radius of a neutron star from its mass (`physics/neutron_star.py`): the band of Miller et al. 2021 at the
cited masses, the measurements it must contain, that a star's radius lies in the band at its mass, that no code reads
the old 10 to 13 km range, and determinism.
"""

import statistics

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import NeutronStar
from planetgen.physics import constants, neutron_star
from planetgen.util import draw


def test_the_band_is_the_papers_at_the_cited_masses():
    assert neutron_star.radius_band_km(1.4) == pytest.approx((12.45, 0.65))
    assert neutron_star.radius_band_km(2.08) == pytest.approx((12.35, 0.75))


def test_the_band_between_the_masses_is_the_line():
    centre, sigma = neutron_star.radius_band_km(1.74)
    assert centre == pytest.approx(12.40, abs=0.01) and sigma == pytest.approx(0.70, abs=0.01)


@pytest.mark.parametrize("mass, radius", [
    (1.34, 12.71),  # PSR J0030+0451, Riley et al. 2019
    (2.08, 13.7),  # PSR J0740+6620, Miller et al. 2021
])
def test_the_measured_radii_lie_inside_the_clipped_band(mass, radius):
    lowest, highest = neutron_star.radius_limits_km(mass)
    assert lowest < radius < highest


def test_the_limits_over_the_mass_range_are_about_11_to_14_km():
    """The old range reached 10 km, a floor no cited measurement supports."""
    for mass in tuning.NEUTRON_STAR_MASS_RANGE_SOLAR + (1.4, 1.8):
        lowest, highest = neutron_star.radius_limits_km(mass)
        assert 10.7 < lowest < 11.4 and 13.6 < highest < 14.0


def test_the_deviation_is_clipped():
    assert neutron_star.radius_km(1.4, 50.0) == neutron_star.radius_km(1.4, 2.0)
    assert neutron_star.radius_km(1.4, -50.0) == neutron_star.radius_km(1.4, -2.0)
    assert neutron_star.radius_km(1.4) == pytest.approx(12.45)


def test_a_drawn_star_is_in_the_band_at_its_mass_and_its_radii_spread():
    draw.set_run_seed(463)
    stars = [NeutronStar(SystemConfig()) for _ in range(300)]
    for star in stars:
        low, high = neutron_star.radius_limits_km(star.mass_solar)
        assert low - 1e-9 <= star.radius <= high + 1e-9
    radii = [star.radius for star in stars]
    assert 12.0 < statistics.mean(radii) < 12.9
    assert 0.4 < statistics.pstdev(radii) < 0.9
    assert min(radii) < 11.9 and max(radii) > 13.0


def test_a_neutron_star_companion_has_the_centre_radius_of_its_mass():
    leaf = multistar._remnant_leaf("neutron-star", 15.0, 1.4)
    assert leaf.radius_au * constants.AU_M / 1000.0 == pytest.approx(12.45)


def test_the_same_seed_gives_the_same_radius():
    def radii():
        draw.set_run_seed(7)
        return [NeutronStar(SystemConfig()).radius for _ in range(5)]
    assert radii() == radii()


