# tests/test_irregular_zone.py

"""
GEN.1076: where an irregular moon can stay (`planetgen.physics.irregular_zone`): the Hill radii of the four giants, the stable limits of both senses, the Kozai band and the inner edge from the
outermost regular moon, Jewitt and Haghighipour 2007 and Vokrouhlicky, Nesvorny and Levison 2008.
"""

import pytest

from planetgen.physics import irregular_zone as zone

SUN = 1.989e30
JUPITER = (5.203, 1.898e27)
SATURN = (9.537, 5.683e26)
URANUS = (19.19, 8.681e25)
NEPTUNE = (30.07, 1.024e26)
CALLISTO_AU, IAPETUS_AU, OBERON_AU = 1.8827e6 / 1.495979e8, 3.5613e6 / 1.495979e8, 5.8352e5 / 1.495979e8


def test_the_hill_radii_of_the_four_giants_are_those_of_jewitt_and_haghighipour_table_1():
    for (distance, mass), expected in zip((JUPITER, SATURN, URANUS, NEPTUNE), (0.35, 0.43, 0.47, 0.77)):
        assert zone.hill_radius_au(distance, mass, SUN) == pytest.approx(expected, abs=0.012)


def test_the_stable_limits_are_a_share_of_the_hill_radius_and_the_retrograde_one_reaches_farther():
    jupiter = zone.zone(*JUPITER, SUN, CALLISTO_AU)
    assert jupiter.prograde_limit_au == pytest.approx(0.4 * jupiter.hill_radius_au)
    assert jupiter.retrograde_limit_au == pytest.approx(0.55 * jupiter.hill_radius_au)
    assert jupiter.retrograde_limit_au > jupiter.prograde_limit_au
    assert zone.stable_limit_au(jupiter, 150.0) == jupiter.retrograde_limit_au
    assert zone.stable_limit_au(jupiter, 30.0) == jupiter.prograde_limit_au
    # 26 to 36 million km at Jupiter (design 4.5), 0.17 to 0.24 AU
    assert 0.14 < jupiter.prograde_limit_au < 0.15 and 0.19 < jupiter.retrograde_limit_au < 0.20


def test_the_inner_edge_is_the_outermost_regular_moon_of_each_giant():
    for planet, moon, edge in ((JUPITER, CALLISTO_AU, 0.013), (SATURN, IAPETUS_AU, 0.023), (URANUS, OBERON_AU, 0.004)):
        assert zone.zone(*planet, SUN, moon).inner_edge_au == pytest.approx(edge, abs=0.001)
    assert zone.zone(*JUPITER, SUN).inner_edge_au == 0.0


def test_no_body_stays_in_the_kozai_band_and_the_band_is_the_moon_codes():
    from planetgen.physics import moon_composition
    assert (55.0, 130.0) == moon_composition.KOZAI_FORBIDDEN_DEG
    jupiter = zone.zone(*JUPITER, SUN, CALLISTO_AU)
    for inclination in (55.0, 90.0, 130.0):
        assert zone.verdict(jupiter, 0.1, 0.3, inclination) == "kozai"
    assert zone.verdict(jupiter, 0.1, 0.3, 54.0) is None and zone.verdict(jupiter, 0.1, 0.3, 131.0) is None


def test_a_retrograde_orbit_stays_where_a_prograde_one_does_not():
    jupiter = zone.zone(*JUPITER, SUN, CALLISTO_AU)
    between = 0.5 * (jupiter.prograde_limit_au + jupiter.retrograde_limit_au)
    assert zone.verdict(jupiter, between, 0.2, 40.0) == "unstable"
    assert zone.verdict(jupiter, between, 0.2, 150.0) is None
    assert zone.verdict(jupiter, jupiter.retrograde_limit_au * 1.01, 0.2, 150.0) == "unstable"


def test_a_pericentre_inside_the_outermost_regular_moon_collides():
    jupiter = zone.zone(*JUPITER, SUN, CALLISTO_AU)
    assert zone.verdict(jupiter, 0.05, 0.7, 150.0) is None                                  # q = 0.015 is past Callisto's 0.0126
    assert zone.verdict(jupiter, 0.05, 0.8, 150.0) == "collides with regular moons"      # q = 0.010


def test_observed_medians_of_a_over_hill_radius_lie_inside_the_stable_limits():
    for (distance, mass), median in zip((JUPITER, SATURN, URANUS, NEPTUNE), (0.44, 0.29, 0.17, 0.19)):
        z = zone.zone(distance, mass, SUN)
        assert zone.verdict(z, median * z.hill_radius_au, 0.3, 150.0) is None


def test_a_hot_planets_zone_is_empty_and_an_earths_is_tiny():
    hot_jupiter = zone.zone(0.05, 1.898e27, SUN, 0.003)
    assert hot_jupiter.hill_radius_au == pytest.approx(0.0034, rel=0.1) and hot_jupiter.empty
    assert not zone.zone(*JUPITER, SUN, CALLISTO_AU).empty
    earth = zone.zone(1.0, 5.972e24, SUN)
    assert earth.hill_radius_au * 1.495979e8 == pytest.approx(1.5e6, rel=0.05)
