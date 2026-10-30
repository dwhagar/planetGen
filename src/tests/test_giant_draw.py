"""
GEN.460 and GEN.642 wired: a giant's heavy elements, radius and interior drawn together (`planetgen.physics.giant_draw`).
"""

from types import SimpleNamespace

import pytest

from planetgen.physics import bulk_composition as bc, constants, giant_draw, giant_radius, orbit_elements
from planetgen.util import draw

JUPITER_MASS_KG = constants.JUPITER_MASS_TO_KG


def _star(**extra):
    fields = dict(metallicity_fe_h=0.0, age=4.5, mass=2e30, galactic_center_dist_ly=26000.0, luminosity=3.828e26,
                  snow_line_au=2.7)
    fields.update(extra)
    return SimpleNamespace(**fields)


def _giant(mass_kg, distance=5.2, **extra):
    fields = dict(radius=None, mass=mass_kg, star=_star(), is_moon=False, distance=distance, eccentricity=0.0)
    fields.update(extra)
    return SimpleNamespace(**fields)


def _fallback(mass_kg):
    return 70000.0


def test_a_jupiter_mass_giant_takes_its_radius_from_the_tables():
    draw.set_run_seed(1)
    body = _giant(JUPITER_MASS_KG)
    giant_draw.draw_giant(body, False, _fallback)
    assert body.radius != 70000.0
    assert 0.8 < body.radius / constants.JUPITER_RADIUS_KM < 1.2
    assert body.density == pytest.approx(body.mass / (4 / 3 * 3.141592653589793 * (body.radius * 1000) ** 3) / 1000)


def test_a_hot_giant_is_inflated_and_a_cold_one_is_not():
    draw.set_run_seed(2)
    hot, cold = _giant(JUPITER_MASS_KG, distance=0.03), _giant(JUPITER_MASS_KG, distance=5.2)
    giant_draw.draw_giant(hot, False, _fallback)
    giant_draw.draw_giant(cold, False, _fallback)
    assert hot.radius > 1.1 * cold.radius


def test_an_ice_giant_is_outside_the_tables_and_uses_the_fallback_relation():
    draw.set_run_seed(3)
    body = _giant(17.0 * constants.EARTH_MASS_TO_KG, distance=19.0)
    giant_draw.draw_giant(body, False, _fallback)
    assert body.radius == 70000.0
    bulk = bc.bulk_fractions(body)
    assert bulk.envelope < 0.3 and bulk.water + bulk.ice > 0.0


def test_a_given_radius_is_kept():
    draw.set_run_seed(4)
    body = _giant(JUPITER_MASS_KG, radius=65000.0)
    giant_draw.draw_giant(body, True, _fallback)
    assert body.radius == 65000.0


def test_the_composition_adds_up_and_holds_the_heavy_element_law():
    for seed, masses in ((5, (0.3, 1.0, 5.0)), (6, (10.0, 30.0)), (7, (15.0, 50.0))):
        draw.set_run_seed(seed)
        for mass in masses:
            kg = mass * (JUPITER_MASS_KG if seed == 5 else constants.EARTH_MASS_TO_KG)
            body = _giant(kg)
            giant_draw.draw_giant(body, False, _fallback)
            bulk = bc.bulk_fractions(body)
            assert sum(bulk.shares().values()) == pytest.approx(1.0, abs=1e-9)
            assert 0.0 <= bulk.envelope_z < 0.5 and 0.2 < bulk.envelope_helium / (1.0 - bulk.envelope_z) < 0.35
            if giant_draw.regime_of(body.mass / constants.EARTH_MASS_TO_KG) == giant_draw.GAS:
                heavy = (1.0 - bulk.envelope + bulk.envelope * bulk.envelope_z) * body.mass / constants.EARTH_MASS_TO_KG
                assert heavy >= giant_radius.MIN_HEAVY_FRACTION * body.mass / constants.EARTH_MASS_TO_KG - 1e-9


def test_ices_come_only_from_beyond_the_snow_line():
    draw.set_run_seed(8)
    inside, outside = _giant(JUPITER_MASS_KG, distance=1.0), _giant(JUPITER_MASS_KG, distance=6.0)
    giant_draw.draw_giant(inside, False, _fallback)
    giant_draw.draw_giant(outside, False, _fallback)
    assert bc.bulk_fractions(inside).water == 0.0 and bc.bulk_fractions(outside).water > 0.0


def test_a_giant_moon_has_no_star_distance_and_takes_the_fallback():
    draw.set_run_seed(9)
    body = _giant(JUPITER_MASS_KG, is_moon=True)
    giant_draw.draw_giant(body, False, _fallback)
    assert body.radius == 70000.0 and bc.bulk_fractions(body) is not None


def test_the_same_seed_makes_the_same_giant():
    results = []
    for _ in range(2):
        draw.set_run_seed(77)
        body = _giant(2.0 * JUPITER_MASS_KG)
        giant_draw.draw_giant(body, False, _fallback)
        results.append((body.radius, body.metal_mass_fraction, body.envelope_mass_fraction))
    assert results[0] == results[1]


def test_jupiter_gets_helium_rain_and_a_hot_giant_does_not():
    own = giant_draw.own_flux_w_m2(_giant(JUPITER_MASS_KG, radius=constants.JUPITER_RADIUS_KM))
    cold = giant_draw.ONE_BAR_OVER_EFFECTIVE * giant_draw.effective_temperature_k(50.0, own)
    hot = giant_draw.ONE_BAR_OVER_EFFECTIVE * giant_draw.effective_temperature_k(5e5, own)
    assert 140.0 < cold < 185.0 and hot > 1000.0


# GEN.816: the one-bar temperature carries the giant's own heat. (mass in Earth masses, radius km, distance AU, measured 1 bar K, tolerance)
SOLAR_SYSTEM_GIANTS = (("Jupiter", 317.8, 69911.0, 5.2, 166.0, 0.08), ("Uranus", 14.5, 25362.0, 19.2, 76.0, 0.10),
                       ("Neptune", 17.1, 24622.0, 30.07, 72.0, 0.10),
                       # Saturn's own heat is twice what the cooling law gives (helium rain, GEN.810 adds it): 134 K measured, 114 K here
                       ("Saturn", 95.2, 58232.0, 9.58, 134.0, 0.17))


@pytest.mark.parametrize("name, mass_earth, radius_km, distance_au, measured_k, tolerance", SOLAR_SYSTEM_GIANTS)
def test_the_one_bar_temperature_of_the_solar_system_giants_is_near_the_measured_one(name, mass_earth, radius_km, distance_au,
                                                                                    measured_k, tolerance):
    body = _giant(mass_earth * constants.EARTH_MASS_TO_KG, distance=distance_au, radius=radius_km)
    flux = body.star.luminosity / (4 * 3.141592653589793 * (distance_au * constants.AU_TO_KM) ** 2) / 1e6
    one_bar = giant_draw.ONE_BAR_OVER_EFFECTIVE * giant_draw.effective_temperature_k(flux, giant_draw.own_flux_w_m2(body))
    assert one_bar == pytest.approx(measured_k, rel=tolerance), name


def test_a_giant_without_its_own_heat_is_far_colder_at_one_bar():
    flux = 50.34
    with_heat = giant_draw.effective_temperature_k(flux, 5.6)
    assert with_heat > 1.15 * giant_draw.effective_temperature_k(flux, 0.0)


def test_a_giants_own_heat_fades_with_age_and_grows_with_mass():
    young = giant_draw.own_flux_w_m2(_giant(JUPITER_MASS_KG, radius=constants.JUPITER_RADIUS_KM, star=_star(age=1.0)))
    old = giant_draw.own_flux_w_m2(_giant(JUPITER_MASS_KG, radius=constants.JUPITER_RADIUS_KM, star=_star(age=9.0)))
    light = giant_draw.own_flux_w_m2(_giant(0.3 * JUPITER_MASS_KG, radius=constants.JUPITER_RADIUS_KM))
    assert young > old
    assert giant_draw.own_flux_w_m2(_giant(JUPITER_MASS_KG, radius=constants.JUPITER_RADIUS_KM)) > light


# --- GEN.961: a table giant is sized again once the star's age and the orbit are final ------------------------------------------------------

def _drawn_giant(seed, mass_earth=300.0, distance=5.2, age=1.0, eccentricity=0.0):
    draw.set_run_seed(seed)
    body = _giant(mass_earth * constants.EARTH_MASS_TO_KG, distance=distance, star=_star(age=age), eccentricity=eccentricity)
    giant_draw.draw_giant(body, False, _fallback)
    assert getattr(body, giant_draw.SETTLE_ATTRIBUTE) is not None, "this giant was not sized from the tables"
    return body


def _table_radius_at(body, age, eccentricity=0.0):
    drawn = getattr(body, giant_draw.SETTLE_ATTRIBUTE)
    mass_earth = body.mass / constants.EARTH_MASS_TO_KG
    core = giant_radius.table_core_mass_earth(mass_earth, drawn.heavy_earth, drawn.envelope_z)
    distance = orbit_elements.flux_distance_au(body.distance, eccentricity)
    light = body.star.luminosity / (4.0 * 3.141592653589793 * (distance * constants.AU_TO_KM) ** 2) / 1e6
    return giant_radius.radius_km(mass_earth, core, light / giant_radius.flux_w_m2(1.0), age)


def test_a_giant_drawn_young_is_settled_to_the_tables_radius_at_the_final_age():
    body = _drawn_giant(21, mass_earth=100.0, distance=0.3)
    drawn_radius = body.radius
    body.star.age = 9.0
    assert giant_draw.settle(body) is True
    assert body.radius < drawn_radius
    assert body.radius == pytest.approx(_table_radius_at(body, 9.0), rel=1e-9)
    assert getattr(body, giant_draw.SETTLE_ATTRIBUTE).radius_km == body.radius
    assert body.density == pytest.approx(body.mass / (4 / 3 * 3.141592653589793 * (body.radius * 1000) ** 3) / 1000)
    assert sum(bc.bulk_fractions(body).shares().values()) == pytest.approx(1.0, abs=1e-9)


def test_a_giant_whose_age_did_not_move_is_left_as_drawn():
    body = _drawn_giant(22)
    before = (body.radius, body.density, body.composition, body.metal_mass_fraction)
    assert giant_draw.settle(body) is False
    assert (body.radius, body.density, body.composition, body.metal_mass_fraction) == before


def test_a_giant_on_an_eccentric_orbit_is_sized_at_the_mean_light():
    body = _drawn_giant(23, mass_earth=300.0, distance=0.05)
    body.eccentricity = 0.6
    assert giant_draw.settle(body) is True
    assert body.radius == pytest.approx(_table_radius_at(body, body.star.age, 0.6), rel=1e-9)


def test_a_change_under_the_tolerance_keeps_the_interior():
    body = _drawn_giant(24, mass_earth=100.0, distance=5.2)
    before = (body.radius, body.composition)
    body.star.age = 1.0 + 1e-6
    assert giant_draw.settle(body) is False
    assert (body.radius, body.composition) == before


def test_giants_the_tables_do_not_size_are_not_settled():
    draw.set_run_seed(25)
    for body in (_giant(17.0 * constants.EARTH_MASS_TO_KG, distance=19.0, eccentricity=0.0),
                 _giant(JUPITER_MASS_KG, radius=65000.0, eccentricity=0.0),
                 _giant(JUPITER_MASS_KG, is_moon=True, eccentricity=0.0)):
        giant_draw.draw_giant(body, body.radius is not None, _fallback)
        radius = body.radius
        body.star.age = 9.0
        assert giant_draw.settle(body) is False and body.radius == radius


def test_settling_is_the_same_for_the_same_giant_and_leaves_the_run_stream_where_it_was():
    after = []
    for settle in (True, False):
        body = _drawn_giant(26, mass_earth=100.0, distance=0.3)
        body.star.age = 9.0
        if settle:
            giant_draw.settle(body)
            first = (body.radius, body.composition, body.metal_mass_fraction)
        after.append(draw.random())
    assert after[0] == after[1]
    again = _drawn_giant(26, mass_earth=100.0, distance=0.3)
    again.star.age = 9.0
    giant_draw.settle(again)
    assert (again.radius, again.composition, again.metal_mass_fraction) == first
