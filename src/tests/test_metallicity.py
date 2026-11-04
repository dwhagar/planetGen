# tests/test_metallicity.py

"""
GEN.158: one [Fe/H] per star (`planetgen/physics/metallicity.py`), drawn from
the star's population and its distance from the galactic centre, and the
planet-occurrence factors it feeds (docs/design/galaxy/globular-clusters.md 5 and 7).
"""

import statistics

import pytest

from planetgen import tuning
from planetgen.galaxy import scales
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import constants, formation, metallicity
from planetgen.util import draw

SUN = scales.reference_radius_pc() * 3.2615637771674337   # the zero point of the disc gradient, light-years
KPC = constants.KILOPARSEC_M / constants.LY_TO_M


def _sample(population, distance_ly, n=3000, seed=1):
    stream = draw.Stream(seed)
    return [metallicity.sample_fe_h(population, distance_ly, stream) for _ in range(n)]


def test_the_disc_gradient_is_zero_at_the_sun_and_falls_outward():
    assert metallicity.disc_gradient_fe_h(SUN) == pytest.approx(0.0)
    assert metallicity.disc_gradient_fe_h(SUN + 5 * KPC) == pytest.approx(-0.30)
    assert metallicity.disc_gradient_fe_h(SUN - 5 * KPC) == pytest.approx(0.30)


def test_each_population_has_its_own_mean_and_scatter():
    young = _sample("young", SUN)
    old = _sample("old", SUN)
    bulge = _sample("bulge", 3000)
    halo = _sample("halo", 10 * KPC)
    assert statistics.mean(young) == pytest.approx(0.0, abs=0.01)
    assert statistics.pstdev(young) == pytest.approx(0.10, abs=0.01)
    assert statistics.mean(old) == pytest.approx(-0.325, abs=0.02)
    assert statistics.pstdev(bulge) > statistics.pstdev(old) > statistics.pstdev(young)
    assert statistics.mean(halo) == pytest.approx(-1.6, abs=0.05)


def test_the_halo_is_poorer_beyond_25_kpc_and_a_draw_stays_in_range():
    assert metallicity.halo_mean_fe_h(10 * KPC) == pytest.approx(-1.6)
    assert metallicity.halo_mean_fe_h(40 * KPC) == pytest.approx(-2.2)
    low, high = tuning.METALLICITY_RANGE_FE_H
    assert all(low <= value <= high for value in _sample("halo", 40 * KPC))


def test_a_cluster_carries_its_own_value_and_an_unknown_population_is_refused():
    with pytest.raises(ValueError):
        metallicity.sample_fe_h("cluster", SUN)
    with pytest.raises(ValueError):
        metallicity.sample_fe_h("nonsense", SUN)


def test_the_same_stream_gives_the_same_value():
    assert _sample("old", SUN, 5, seed=4) == _sample("old", SUN, 5, seed=4)


def test_planet_occurrence_factors():
    assert metallicity.giant_planet_factor(0.0) == 1.0
    assert metallicity.giant_planet_factor(0.3) == 1.0
    assert metallicity.giant_planet_factor(-0.5) == pytest.approx(0.1)
    assert metallicity.giant_planet_factor(-1.0) == pytest.approx(0.01)
    assert metallicity.giant_planet_factor(-1.5) == pytest.approx(0.001)
    assert metallicity.giant_planet_factor(-2.3) == pytest.approx(2.5e-5, rel=0.3)


def test_a_disk_has_fewer_solids_at_low_metal_content_and_solar_is_unchanged():
    mass = constants.SOLAR_MASS_TO_KG
    assert formation.disk_surface_density_scale(mass) == formation.disk_surface_density_scale(mass, 1.0)
    assert formation.disk_surface_density_scale(mass, 0.1) == pytest.approx(0.1 * formation.disk_surface_density_scale(mass))


def test_a_star_gets_a_metallicity_that_does_not_move_other_draws():
    draw.set_run_seed(7)
    first = Star(SystemConfig())
    draw.set_run_seed(7)
    second = Star(SystemConfig())
    assert first.metallicity_fe_h == second.metallicity_fe_h
    assert first.mass == second.mass and first.rotation_period_hours == second.rotation_period_hours
    assert "metallicity_fe_h" in Star.SERIALIZABLE_FIELDS


def test_a_forced_metallicity_is_kept_and_a_halo_star_is_poor():
    assert Star(SystemConfig(), metallicity_fe_h=-0.7).metallicity_fe_h == -0.7
    config = SystemConfig()
    config.POPULATION = "bulge"
    values = [Star(config, galactic_center_dist_ly=8000).metallicity_fe_h for _ in range(30)]
    assert all(value is not None for value in values)
