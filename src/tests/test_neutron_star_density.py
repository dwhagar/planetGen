# tests/test_neutron_star_density.py

"""
GEN.554: the neutron stars are what the black holes leave of the massive stars (docs/design/systems/neutron-star-density.md):
every star of 8 to 150 solar masses ends as one or the other, so a place's neutron stars per star are the same sum as its
black holes' with `1 - F_p` for `F_p`, less the share the galaxy ejects. Boss's 4.0e-3 per star is a test target.
"""

import math

import pytest

from planetgen import tuning
from planetgen.galaxy import ejecta_field, profile, remnant_density
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from planetgen.physics import black_hole_formation as formation
from planetgen.physics.units import ly_to_pc

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
    arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
E_VALUE = expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)
SUN = (8128.0, 0.0, 0.0)
NS, BH = remnant_density.NEUTRON_STAR, remnant_density.BLACK_HOLE


def _shape(solar, free, unbound=0.3):
    return SHAPE._replace(structure=SHAPE.structure._replace(
        black_hole_solar_share=solar, black_hole_metal_free_share=free, neutron_star_unbound_share=unbound))


def test_every_massive_star_leaves_a_neutron_star_or_a_black_hole():
    """At every distance, each population's formed neutron stars and black holes add up to the remnants it has made."""
    for distance in (1.0e3, 2.6e4, 1.0e5):
        neutron = remnant_density.population_factors(SHAPE.structure, distance, NS, bound=False)
        hole = remnant_density.population_factors(SHAPE.structure, distance, BH)
        for name in neutron:
            assert neutron[name] + hole[name] == pytest.approx(remnant_density._formed(name) * 1.0)


def test_the_galaxys_formed_remnants_never_exceed_its_massive_stars():
    for solar, free in ((0.0, 0.15), (0.2, 0.3), (0.247, 0.247)):
        shape = _shape(solar, free)
        total = remnant_density.formed_per_star(shape, NS) + remnant_density.formed_per_star(shape, BH)
        assert total <= formation.imf_massive_share() / (1.0 - formation.imf_massive_share()) * 1.0001


def test_bound_plus_ejected_neutron_stars_are_what_the_galaxy_formed():
    for unbound in (0.2, 0.3, 0.4):
        shape = _shape(0.09, 0.18, unbound)
        systems = ejecta_field.galaxy_systems(shape, E_VALUE, EDGE_PC)
        bound = remnant_density.base_rate_per_star(NS) * remnant_density.mean_factor(shape, NS) * systems
        ejected = ejecta_field.unbound_count(NS, shape, E_VALUE, EDGE_PC)
        formed = remnant_density.formed_per_star(shape, NS) * systems
        assert bound + ejected == pytest.approx(formed, rel=1e-9)
        assert ejected / formed == pytest.approx(unbound)


def test_a_galaxy_that_fails_more_massive_stars_makes_fewer_neutron_stars():
    few, many = _shape(0.0, 0.15), _shape(0.2, 0.3)
    assert remnant_density.formed_per_star(few, NS) > remnant_density.formed_per_star(many, NS)
    assert remnant_density.formed_per_star(few, BH) < remnant_density.formed_per_star(many, BH)


def test_a_metal_rich_place_makes_more_neutron_stars_than_a_metal_poor_one():
    rich = remnant_density.place_factor(SUN, SHAPE, NS)
    poor = remnant_density.place_factor((SUN[0] * 4.0, 0.0, 6000.0), SHAPE, NS)
    assert rich > poor > 0.0


def test_the_neutron_stars_per_star_at_the_solar_circle_are_boss_research_density():
    """The 4.0e-3 per star the generator was tuned to (1e9 neutron stars in 2.5e11, GEN.100) is a target, not an input:
    the derived rate lands within a quarter of it."""
    target = tuning.phenomenon_rate_per_star("neutron-star")
    derived = remnant_density.rate_per_star(kind=NS)
    assert derived == pytest.approx(target, rel=0.25)
    here = remnant_density.base_rate_per_star(NS) * remnant_density.place_factor(SUN, SHAPE, NS)
    assert here == pytest.approx(target, rel=0.25)


def test_a_profiles_neutron_star_factor_is_the_places_and_the_bound_is_never_exceeded():
    factor = profile.profile_at(SUN, SHAPE, EDGE_PC, E_VALUE).kind_factors["neutron-star"]
    assert factor == pytest.approx(remnant_density.place_factor(SUN, SHAPE, NS), rel=0.1)
    assert 0.0 < factor <= remnant_density.bound_factor(SHAPE.structure, NS)


def test_only_the_stellar_black_holes_and_plain_neutron_stars_follow_the_metals():
    follows = remnant_density.follows_metals
    assert follows("black-hole", "stellar") and follows("neutron-star", None)
    assert not follows("black-hole", "intermediate") and not follows("neutron-star", "magnetar")
    assert not follows("planetary-nebula", None)
    with pytest.raises(ValueError):
        remnant_density.base_rate_per_star("planetary-nebula")
