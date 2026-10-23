# tests/test_nucleus_planet_limit.py -- GEN.284: the nucleus' density and dispersion reach the one planet survival
# function (GEN.325), with no second cull, and the profile's dispersion is one thing: the three-axis total.

import math

import pytest

from planetgen import tuning
from planetgen.galaxy import density, profile as galaxy_profile
from planetgen.galaxy import structure as galaxy_structure
from planetgen.galaxy.profile import profile_at
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from planetgen.galaxy.spherical import cell_volume
from planetgen.generation import place_terms
from planetgen.physics import constants, planet_survival
from planetgen.physics.units import ly_to_pc

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
BASE = dict(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
            bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
SHAPE = density.build_galaxy_shape(**BASE, nuclear_amplitude=1.0)


def _terms_at(radius_pc, shape=SHAPE):
    profile = profile_at((float(radius_pc), 0.0, 0.0), shape, EDGE_PC,
                         expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY), address=(0, 0, 0),
                         volume_pc3=cell_volume(0, 0, 0, EDGE_PC))
    return profile, place_terms.from_profile(profile)


def test_the_nucleus_cell_feeds_its_density_and_dispersion_to_the_one_survival_function():
    profile, terms = _terms_at(1.0)
    assert profile.region == "nuclear"
    assert terms.systems_per_pc3 > 1.0e5
    assert 100.0 < terms.dispersion_kms < 200.0           # the cluster's 150 km/s, one-axis
    limit = place_terms.planet_limit_au(terms, tuning.PLANET_SURVIVAL_AGE_GY)
    assert limit == planet_survival.planet_survival_limit_au(terms.systems_per_pc3, terms.dispersion_kms)
    assert limit == pytest.approx(0.06, rel=0.25)         # the core note's row at 1 pc


def test_the_limit_follows_the_note_s_table_from_the_cluster_out_through_the_nuclear_disc():
    limits = [place_terms.planet_limit_au(_terms_at(r)[1], tuning.PLANET_SURVIVAL_AGE_GY)
              for r in (0.0, 3.0, 10.0, 30.0, 100.0, 200.0)]
    assert all(a <= b for a, b in zip(limits, limits[1:]))
    assert limits[0] < 0.15                                # inside a few parsecs the cluster holds no planets
    assert 0.5 < _limit(10.0) < 3.0                        # the note: about 2 AU at 10 pc
    assert 3.0 < _limit(100.0) < 11.0                      # the nuclear disc keeps terrestrial-zone orbits
    assert _limit(8000.0) > 100.0                          # the Sun's neighbourhood keeps orbits out to hundreds of AU


def _limit(radius_pc):
    return place_terms.planet_limit_au(_terms_at(radius_pc)[1], tuning.PLANET_SURVIVAL_AGE_GY)


def test_every_population_is_stated_as_a_three_axis_total():
    for name in ("cluster", "satellite", "nuclear_cluster", "nuclear_disc"):
        assert galaxy_profile._dispersion_kms(name, SHAPE) == pytest.approx(
            tuning.STELLAR_POPULATION_DISPERSION_KMS[name] * math.sqrt(3.0))
    assert galaxy_profile._dispersion_kms("old", SHAPE) == pytest.approx(
        galaxy_structure.disc_dispersion_kms(SHAPE.structure, "old"))
    assert constants.ISOTROPIC_TOTAL_OVER_ONE_AXIS == pytest.approx(math.sqrt(3.0))


def test_a_disc_systems_encounter_limit_uses_the_three_axis_total_over_root_three():
    profile, terms = _terms_at(8000.0)
    assert set(profile.populations) >= {"young", "old"}
    weighted = sum(p.share * p.dispersion_kms ** 2 for p in profile.populations.values())
    assert terms.dispersion_kms == pytest.approx(math.sqrt(weighted / 3.0))
    one_axis = place_terms.PlaceTerms(terms.mass_density_solar_pc3, terms.dispersion_kms, terms.perturber_mass_solar,
                                      terms.mean_system_mass_solar)
    three_axis = place_terms.PlaceTerms(terms.mass_density_solar_pc3, terms.dispersion_kms * math.sqrt(3.0),
                                        terms.perturber_mass_solar, terms.mean_system_mass_solar)
    assert (place_terms.encounter_limit_au(three_axis, 2.0, 5.0) / place_terms.encounter_limit_au(one_axis, 2.0, 5.0)
            == pytest.approx(math.sqrt(3.0)))
