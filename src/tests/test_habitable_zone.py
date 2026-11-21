# tests/test_habitable_zone.py

"""
The habitable zone from the [R2447] fits (GEN.146): the runaway greenhouse limit is the inner edge and
the maximum greenhouse limit the outer edge, both moving with the star's temperature.
"""

import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import constants as physical_constants
from planetgen.physics import orbits

SUN = physical_constants.SOLAR_LUMINOSITY


def test_the_sun_has_the_published_edges():
    inner, outer = orbits.calculate_habitable_zone(SUN, 5780.0)
    assert inner == pytest.approx(0.950, abs=0.001)    # sqrt(1 / 1.107)
    assert outer == pytest.approx(1.676, abs=0.001)    # sqrt(1 / 0.356)


@pytest.mark.parametrize("teff, runaway", [(2800, 0.918), (4000, 0.947), (5780, 1.107), (7200, 1.296)])
def test_runaway_flux_at_the_checked_temperatures(teff, runaway):
    # The values the design note (atmospheres-retention-and-classes.md, 4.1) worked out from the paper's Table 1.
    assert orbits.kopparapu_flux(orbits.KOPPARAPU_RUNAWAY_GREENHOUSE, teff) == pytest.approx(runaway, abs=0.0006)


def test_a_cool_star_has_its_zone_nearer_in_than_the_fixed_fluxes_put_it():
    # The old fixed fluxes (1.1 and 0.53) put the outer edge of a 3,200 K star a third too near.
    luminosity = 0.0016 * SUN
    inner, outer = orbits.calculate_habitable_zone(luminosity, 3200.0)
    old_outer = math.sqrt(0.0016 / 0.53)
    assert outer == pytest.approx(math.sqrt(0.0016 / orbits.kopparapu_flux(
        orbits.KOPPARAPU_MAXIMUM_GREENHOUSE, 3200.0)))
    assert old_outer / outer == pytest.approx(0.67, abs=0.02)    # a third too near
    assert inner < outer


def test_the_zone_is_ordered_and_positive_across_the_fitted_range():
    low, high = orbits.KOPPARAPU_TEFF_RANGE_K
    for teff in range(int(low), int(high) + 1, 100):
        for coefficients in (orbits.KOPPARAPU_RUNAWAY_GREENHOUSE, orbits.KOPPARAPU_MAXIMUM_GREENHOUSE):
            assert orbits.kopparapu_flux(coefficients, teff) > 0
        inner, outer = orbits.calculate_habitable_zone(SUN, float(teff))
        assert 0 < inner < outer


def test_a_temperature_outside_the_fitted_range_takes_the_nearest_end():
    low, high = orbits.KOPPARAPU_TEFF_RANGE_K
    assert orbits.calculate_habitable_zone(SUN, 1500.0) == orbits.calculate_habitable_zone(SUN, low)
    assert orbits.calculate_habitable_zone(SUN, 0.0) == orbits.calculate_habitable_zone(SUN, low)
    assert orbits.calculate_habitable_zone(SUN, 5.0e6) == orbits.calculate_habitable_zone(SUN, high)


def test_the_zone_scales_with_the_square_root_of_the_luminosity():
    inner, outer = orbits.calculate_habitable_zone(SUN, 4500.0)
    inner4, outer4 = orbits.calculate_habitable_zone(4 * SUN, 4500.0)
    assert inner4 == pytest.approx(2 * inner) and outer4 == pytest.approx(2 * outer)


def test_a_generated_stars_zone_follows_its_temperature():
    for star_type in ("G2V", "M3V", "K5V", "F5V"):
        cfg = SystemConfig()
        cfg.STAR_TYPE = star_type
        star = Star(cfg)
        assert star.habitable_zone == orbits.calculate_habitable_zone(star.luminosity, star.temperature)


def test_a_binarys_zone_uses_the_pairs_light_weighted_temperature():
    for _ in range(40):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = False
        cfg.PLANETS = False
        proxy = StarSystem(system_config=cfg).star
        if not hasattr(proxy, "_primary"):
            continue
        first, second = proxy._primary, proxy._secondary
        weighted = ((first.luminosity * first.temperature + second.luminosity * second.temperature)
                    / (first.luminosity + second.luminosity))
        assert proxy.habitable_zone == orbits.calculate_habitable_zone(proxy._effective_luminosity, weighted)
        return
    pytest.fail("no binary system in 40 tries")


def test_the_effort_zone_is_the_habitable_zone_widened_by_half_its_width_each_side():
    # GEN.1807 (PA20): the span doubles, the habitable zone itself keeps the [R2447] edges
    inner, outer = orbits.calculate_habitable_zone(SUN, 5780.0)
    wide_inner, wide_outer = orbits.candidate_zone((inner, outer))
    assert wide_outer - wide_inner == pytest.approx(2.0 * (outer - inner))
    assert wide_inner == pytest.approx(inner - 0.5 * (outer - inner)) and wide_outer == pytest.approx(outer + 0.5 * (outer - inner))
    assert (inner, outer) == pytest.approx((0.950, 1.676), abs=0.001)
    assert (wide_inner, wide_outer) == pytest.approx((0.587, 2.039), abs=0.002)
    assert orbits.candidate_zone((1.0, 1.0)) == (1.0, 1.0)


def test_a_zone_letter_stays_the_habitable_zones_while_the_effort_zone_is_wider():
    from planetgen.physics import planets
    zone = (0.95, 1.676)
    assert planets.zone_for(zone, 0.7) == "h" and planets.in_candidate_zone(zone, 0.7)
    assert planets.zone_for(zone, 1.9) == "c" and planets.in_candidate_zone(zone, 1.9)
    assert not planets.in_candidate_zone(zone, 0.5) and not planets.in_candidate_zone(zone, 2.1)
