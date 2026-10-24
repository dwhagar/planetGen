# tests/test_galaxy_scales.py

"""GEN.404: no Sun or home radius is read; the laws that were quoted at the Sun read the galaxy's own scales
(`galaxy.scales`), and a system with no place is taken at the radius that holds half of its galaxy's disc stars."""

import math
from pathlib import Path

import pytest

from planetgen import tuning
from planetgen.galaxy import remnant_distribution, scales, seed as galaxySeed
from planetgen.generation import clusters
from planetgen.physics import metallicity
from planetgen.physics.units import ly_to_pc

SOURCE = Path(__file__).resolve().parent.parent / "planetgen"
SEED = bytes(range(16))


def test_the_scales_follow_the_galaxy_the_unit_runs_in():
    with galaxySeed.seeded(SEED, "sector", (1, 1, 1), disk_scale_length_pc=2000.0):
        assert scales.disk_scale_length_pc() == 2000.0
        assert scales.reference_radius_pc() == pytest.approx(tuning.REFERENCE_RADIUS_SCALE_LENGTHS * 2000.0)
        assert scales.hypervelocity_reference_radius_pc() == pytest.approx(tuning.HVS_REFERENCE_SCALE_LENGTHS * 2000.0)
        assert scales.default_galactocentric_pc() == pytest.approx(tuning.DISK_HALF_MASS_SCALE_LENGTHS * 2000.0)
    with galaxySeed.seeded(SEED, "sector", (1, 1, 1), disk_scale_length_pc=4000.0):
        assert ly_to_pc(scales.default_galactocentric_ly()) == pytest.approx(scales.default_galactocentric_pc())


def test_with_no_bound_length_the_scales_use_the_length_the_seed_draws():
    with galaxySeed.seeded(SEED, "sector", (1, 1, 1)):
        drawn = scales.disk_scale_length_pc()
    assert drawn > 0
    assert scales.disk_scale_length_pc() > 0      # outside any unit: the fallback seed's galaxy


def test_half_of_an_exponential_disc_lies_inside_the_default_radius():
    x = tuning.DISK_HALF_MASS_SCALE_LENGTHS
    assert 1.0 - (1.0 + x) * math.exp(-x) == pytest.approx(0.5, abs=1e-12)


def test_the_pulsar_profile_scales_with_the_galaxy():
    with galaxySeed.seeded(SEED, "sector", (1, 1, 1), disk_scale_length_pc=1300.0):
        small = remnant_distribution.pulsar_radial_factor(0.5 * scales.reference_radius_pc())
        assert remnant_distribution.pulsar_radial_factor(scales.reference_radius_pc()) == pytest.approx(1.0)
    with galaxySeed.seeded(SEED, "sector", (1, 1, 1), disk_scale_length_pc=2600.0):
        large = remnant_distribution.pulsar_radial_factor(0.5 * scales.reference_radius_pc())
    assert small == pytest.approx(large, rel=1e-9)      # a function of R / R0 alone once the disc's length scales too


def test_a_galaxy_twice_as_large_has_clusters_of_the_same_size_at_twice_the_distance():
    import random
    for distance in (3000.0, 8000.0, 20000.0):
        one = clusters._structure(random.Random(3), distance, 1.0)
        two = clusters._structure(random.Random(3), 2.0 * distance, 2.0)
        assert one == two


def test_the_disc_gradient_is_zero_at_the_galaxys_reference_radius():
    with galaxySeed.seeded(SEED, "sector", (1, 1, 1), disk_scale_length_pc=2600.0):
        reference_ly = scales.reference_radius_pc() * 3.2615637771674337
        assert metallicity.disc_gradient_fe_h(reference_ly) == pytest.approx(0.0, abs=1e-9)
    with galaxySeed.seeded(SEED, "sector", (1, 1, 1), disk_scale_length_pc=1300.0):
        reference_ly = scales.reference_radius_pc() * 3.2615637771674337
        assert metallicity.disc_gradient_fe_h(reference_ly) == pytest.approx(0.0, abs=1e-9)
