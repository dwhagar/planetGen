# tests/test_ring_origin.py

"""
GEN.1054: where a ring comes from (`planetgen.physics.ring_origin`): a moon that fell inside the Roche limit, a comet that passed inside it, none close to the star, Saturn and Phobos
as cases, and a giant with nothing in its past has no ring.
"""

import pytest

from planetgen.physics import ring_origin, rings
from planetgen.physics.ring_origin import CometPassage, FallenMoon, Planet
from planetgen.util import draw

YEAR = 1.0
SATURN = Planet(**ring_origin.SATURN)
JUPITER = Planet(1.898e27, 69911.0, 5.2, 1.989e30)
MARS = Planet(6.417e23, 3389.5, 1.524, 1.989e30)
URANUS = Planet(8.681e25, 25362.0, 19.2, 1.989e30)
AGE = 4.5e9


def test_a_giant_with_nothing_in_its_past_has_no_ring():
    for giant in (SATURN, JUPITER, URANUS):
        result = ring_origin.ring_origin(giant, AGE)
        assert not result.exists and result.origin is None and result.mass_kg == 0.0


def test_saturns_rings_come_out_of_a_mimas_sized_moon_a_hundred_million_years_ago():
    """Saturn's rings hold 1.5e19 kg ([R2077], recalled) and are put at 1e8 to 4e8 years old; a moon of Mimas's 3.7e19 kg that fell then leaves that to a factor of two."""
    result = ring_origin.ring_origin(SATURN, AGE, [FallenMoon(3.7e19, AGE - 1.5e8)])
    assert result.exists and result.origin == "disrupted_moon"
    assert result.mass_kg == pytest.approx(1.5e19, rel=0.5)
    assert result.age_yr == pytest.approx(1.5e8)
    assert result.mass_kg <= rings.MAX_MASS_OVER_PLANET * SATURN.mass_kg


def test_phobos_would_leave_mars_a_ring_that_lasts_about_a_billion_years():
    phobos = 1.066e16
    result = ring_origin.ring_origin(MARS, AGE, [FallenMoon(phobos, AGE - 4.0e7)])
    assert result.exists and result.origin == "disrupted_moon"
    assert result.mass_kg == pytest.approx(0.3 * phobos * 0.923, rel=0.2)
    assert not ring_origin.ring_origin(MARS, 1.0e10, [FallenMoon(phobos, 1.0)]).exists      # worn away in ten billion years


def test_a_ring_loses_mass_by_spreading_and_by_meteoroid_erosion_and_never_gains_any():
    assert ring_origin.spread(1.0e18, 0.0, SATURN) == 1.0e18
    assert ring_origin.spread(0.0, 1.0e9, SATURN) == 0.0
    last = 1.0e18
    for years in (1.0e6, 1.0e7, 1.0e8, 1.0e9, 4.5e9):
        now = ring_origin.spread(1.0e18, years, SATURN)
        assert 0.0 <= now < last
        last = now


def test_saturns_rings_lose_about_two_thirds_of_their_mass_in_the_age_of_the_solar_system():
    """GEN.1171: the spreading time of `RING_SPREAD_TIME_YR` at Saturn's own ring mass, with the erosion of [R5001]'s meteoroid flux added."""
    left = ring_origin.spread(ring_origin.SATURN_RING_MASS_KG, AGE, SATURN)
    assert 0.2 < left / ring_origin.SATURN_RING_MASS_KG < 0.45
    flux = ring_origin.meteoroid_flux_kg_m2_s(SATURN)
    assert flux == pytest.approx(ring_origin.METEOROID_FLUX_KG_M2_S)
    erosion = ring_origin.EROSION_YIELD * flux * ring_origin.annulus_area_m2(SATURN) * 3.15576e7
    assert 1.0e19 / erosion > 1.0e10          # the meteoroids alone would take Saturn's ring in more than ten billion years


def test_a_light_ring_is_worn_away_by_meteoroids_and_a_heavy_one_spreads_faster_than_a_light_one():
    """GEN.1171: 1e13 kg rings of fallen moons no longer persist for ten billion years, and a ring ten times as dense loses its share faster."""
    for planet in (SATURN, MARS, Planet(5.972e24, 6371.0, 1.0, 1.989e30)):
        assert ring_origin.spread(1.0e13, 1.0e9, planet) == 0.0
    heavy = ring_origin.spread(1.0e21, 1.0e8, SATURN) / 1.0e21
    light = ring_origin.spread(1.0e19, 1.0e8, SATURN) / 1.0e19
    assert heavy < light < 1.0


def test_meteoroids_fall_less_thickly_on_a_ring_further_from_the_star():
    near = ring_origin.meteoroid_flux_kg_m2_s(Planet(5.683e26, 58232.0, 5.0, 1.989e30))
    far = ring_origin.meteoroid_flux_kg_m2_s(Planet(5.683e26, 58232.0, 20.0, 1.989e30))
    assert near > far


def test_a_comet_that_passed_inside_the_limit_leaves_a_captured_body_ring():
    result = ring_origin.ring_origin(SATURN, AGE, comet_passages=[CometPassage(1.0e19, AGE - 1.0e8)])
    assert result.exists and result.origin == "captured_body" and result.from_comets_kg > 0.0 and result.from_moons_kg == 0.0
    small = ring_origin.ring_origin(SATURN, AGE, comet_passages=[CometPassage(1.0e7, AGE - 1.0e8)])
    assert not small.exists


def test_the_largest_source_names_the_origin_and_the_ring_is_capped_at_a_ten_thousandth_of_the_planet():
    result = ring_origin.ring_origin(SATURN, AGE, [FallenMoon(1.0e17, AGE - 1.0e7)], [CometPassage(1.0e19, AGE - 1.0e7)])
    assert result.origin == "captured_body"
    heavy = ring_origin.ring_origin(MARS, AGE, [FallenMoon(1.0e23, AGE - 1.0e3)])
    assert heavy.mass_kg == pytest.approx(rings.MAX_MASS_OVER_PLANET * MARS.mass_kg)


def test_no_ring_close_to_the_star():
    hot = Planet(SATURN.mass_kg, SATURN.radius_km, 0.05, SATURN.star_mass_kg)
    assert not ring_origin.ring_origin(hot, AGE, [FallenMoon(3.7e19, AGE - 1.0e7)]).exists
    assert ring_origin.comet_passage_rate_per_gyr(hot) == 0.0


def test_the_comet_rate_is_the_reference_at_saturn_and_scales_with_the_roche_sphere_and_focusing():
    assert ring_origin.comet_passage_rate_per_gyr(SATURN) == pytest.approx(ring_origin.COMET_PASSAGES_PER_GYR)
    jupiter, uranus = ring_origin.comet_passage_rate_per_gyr(JUPITER), ring_origin.comet_passage_rate_per_gyr(URANUS)
    assert jupiter > ring_origin.COMET_PASSAGES_PER_GYR > uranus > 0.0
    assert ring_origin.roche_limit_km(SATURN) == pytest.approx(1.3e5, rel=0.1)       # Saturn's A ring ends at 1.37e5 km


def test_drawn_comet_passages_follow_the_rate_and_come_from_the_units_stream():
    draw.set_run_seed(5)
    counts = [len(ring_origin.draw_comet_passages(JUPITER, 1.0e11)) for _ in range(40)]
    expected = ring_origin.comet_passage_rate_per_gyr(JUPITER) * 100.0
    assert sum(counts) / len(counts) == pytest.approx(expected, rel=0.25)
    draw.set_run_seed(5)
    again = ring_origin.draw_comet_passages(JUPITER, 1.0e11)
    draw.set_run_seed(5)
    assert again == ring_origin.draw_comet_passages(JUPITER, 1.0e11)
    low, high = ring_origin.COMET_MASS_KG_RANGE
    assert all(low <= c.mass_kg <= high and 0.0 <= c.time_yr <= 1.0e11 for c in again)
