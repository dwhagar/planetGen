"""
GEN.460 and GEN.642 wired: a giant's heavy elements, radius and interior drawn together (`planetgen.physics.giant_draw`).
"""

from types import SimpleNamespace

import pytest

from planetgen.physics import bulk_composition as bc, constants, giant_draw, giant_radius
from planetgen.util import draw

JUPITER_MASS_KG = constants.JUPITER_MASS_TO_KG


def _star(**extra):
    fields = dict(metallicity_fe_h=0.0, age=4.5, mass=2e30, galactic_center_dist_ly=26000.0, luminosity=3.828e26,
                  snow_line_au=2.7)
    fields.update(extra)
    return SimpleNamespace(**fields)


def _giant(planet_class, mass_kg, distance=5.2, **extra):
    fields = dict(planet_class=planet_class, radius=None, mass=mass_kg, star=_star(), is_moon=False, distance=distance)
    fields.update(extra)
    return SimpleNamespace(**fields)


def _fallback(planet_class, mass_kg):
    return 70000.0


def test_a_jupiter_mass_giant_takes_its_radius_from_the_tables():
    draw.set_run_seed(1)
    body = _giant("J", JUPITER_MASS_KG)
    giant_draw.draw_giant(body, False, _fallback)
    assert body.radius != 70000.0
    assert 0.8 < body.radius / constants.JUPITER_RADIUS_KM < 1.2
    assert body.density == pytest.approx(body.mass / (4 / 3 * 3.141592653589793 * (body.radius * 1000) ** 3) / 1000)


def test_a_hot_giant_is_inflated_and_a_cold_one_is_not():
    draw.set_run_seed(2)
    hot, cold = _giant("J", JUPITER_MASS_KG, distance=0.03), _giant("J", JUPITER_MASS_KG, distance=5.2)
    giant_draw.draw_giant(hot, False, _fallback)
    giant_draw.draw_giant(cold, False, _fallback)
    assert hot.radius > 1.1 * cold.radius


def test_an_ice_giant_is_outside_the_tables_and_uses_the_fallback_relation():
    draw.set_run_seed(3)
    body = _giant("I", 17.0 * constants.EARTH_MASS_TO_KG, distance=19.0)
    giant_draw.draw_giant(body, False, _fallback)
    assert body.radius == 70000.0
    bulk = bc.bulk_fractions(body)
    assert bulk.envelope < 0.3 and bulk.water + bulk.ice > 0.0


def test_a_given_radius_is_kept():
    draw.set_run_seed(4)
    body = _giant("J", JUPITER_MASS_KG, radius=65000.0)
    giant_draw.draw_giant(body, True, _fallback)
    assert body.radius == 65000.0


def test_the_composition_adds_up_and_holds_the_heavy_element_law():
    for seed, planet_class, masses in ((5, "J", (0.3, 1.0, 5.0)), (6, "T", (10.0, 30.0)), (7, "I", (15.0, 50.0))):
        draw.set_run_seed(seed)
        for mass in masses:
            kg = mass * (JUPITER_MASS_KG if planet_class == "J" else constants.EARTH_MASS_TO_KG)
            body = _giant(planet_class, kg)
            giant_draw.draw_giant(body, False, _fallback)
            bulk = bc.bulk_fractions(body)
            assert sum(bulk.shares().values()) == pytest.approx(1.0, abs=1e-9)
            assert 0.2 < bulk.envelope_helium < 0.35 and 0.0 <= bulk.envelope_z < 0.5
            if planet_class == "J":
                heavy = (1.0 - bulk.envelope + bulk.envelope * bulk.envelope_z) * body.mass / constants.EARTH_MASS_TO_KG
                assert heavy >= giant_radius.MIN_HEAVY_FRACTION * body.mass / constants.EARTH_MASS_TO_KG - 1e-9


def test_ices_come_only_from_beyond_the_snow_line():
    draw.set_run_seed(8)
    inside, outside = _giant("J", JUPITER_MASS_KG, distance=1.0), _giant("J", JUPITER_MASS_KG, distance=6.0)
    giant_draw.draw_giant(inside, False, _fallback)
    giant_draw.draw_giant(outside, False, _fallback)
    assert bc.bulk_fractions(inside).water == 0.0 and bc.bulk_fractions(outside).water > 0.0


def test_a_giant_moon_has_no_star_distance_and_takes_the_fallback():
    draw.set_run_seed(9)
    body = _giant("J", JUPITER_MASS_KG, is_moon=True)
    giant_draw.draw_giant(body, False, _fallback)
    assert body.radius == 70000.0 and bc.bulk_fractions(body) is not None


def test_the_same_seed_makes_the_same_giant():
    results = []
    for _ in range(2):
        draw.set_run_seed(77)
        body = _giant("J", 2.0 * JUPITER_MASS_KG)
        giant_draw.draw_giant(body, False, _fallback)
        results.append((body.radius, body.metal_mass_fraction, body.envelope_mass_fraction))
    assert results[0] == results[1]


def test_jupiter_gets_helium_rain_and_a_hot_giant_does_not():
    cold = giant_draw.ONE_BAR_OVER_EFFECTIVE * giant_draw.effective_temperature_k(50.0, 317.8)
    hot = giant_draw.ONE_BAR_OVER_EFFECTIVE * giant_draw.effective_temperature_k(5e5, 317.8)
    assert 140.0 < cold < 185.0 and hot > 1000.0
