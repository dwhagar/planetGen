# tests/test_trojans.py

"""The Lagrange trojan rules (physics/trojans.py): the Gascheau-Routh limit, the finite-trojan rule, libration and the labels."""

import math

import pytest

from planetgen.physics import trojans

SUN, JUPITER, EARTH = 1.0, 9.546e-4, 3.003e-6
BROWN_DWARF = 0.03


def test_the_critical_mass_parameter_is_the_published_0385():
    assert trojans.MU_CRIT == pytest.approx(0.0385, abs=5e-5)
    assert trojans.NUMBERS["mu_crit"]["value"] == pytest.approx(trojans.MU_CRIT, abs=1e-7)


def test_a_sun_jupiter_pair_holds_a_class_m_and_a_class_j_trojan():
    assert trojans.stable(SUN, JUPITER, EARTH)
    assert trojans.stable(SUN, JUPITER, 35 * EARTH)
    assert trojans.stable(SUN, JUPITER, JUPITER)


def test_two_stars_hold_nothing_and_a_brown_dwarf_companion_does():
    assert not trojans.stable(1.0, 0.5, 0.0)
    assert not trojans.stable(1.0, 0.1, EARTH)
    assert trojans.stable(1.0, BROWN_DWARF, EARTH)
    assert not trojans.stable(1.0, 0.04, 0.01)


def test_the_room_for_a_trojan_closes_exactly_at_the_limit():
    room = trojans.max_trojan_mass(SUN, JUPITER)
    assert room == pytest.approx(0.0395, abs=1e-3)
    assert trojans.stable(SUN, JUPITER, room * 0.99)
    assert not trojans.stable(SUN, JUPITER, room * 1.01)
    assert trojans.max_trojan_mass(1.0, 0.5) == 0.0


def test_a_low_mass_host_limits_a_class_j_trojan_to_a_few_jupiters():
    host = 0.1
    room_jupiters = trojans.max_trojan_mass(host, JUPITER) / JUPITER
    assert 2.5 < room_jupiters < 4.0
    assert trojans.stable(host, JUPITER, EARTH)


def test_a_giants_moon_pair_scales_with_the_giants_mass():
    giant_earth = 318.0
    room_earth = trojans.max_trojan_mass(giant_earth, 1.0)
    assert 10.0 < room_earth < 12.0


def test_jupiters_tadpole_librates_in_about_twelve_orbits():
    mu = trojans.mass_parameter(SUN, JUPITER)
    assert trojans.libration_period(11.86, mu) / 11.86 == pytest.approx(12.5, abs=0.3)
    assert trojans.libration_period(1.0, trojans.MU_CRIT * 1.01) is None


def test_the_libration_slows_as_the_secondary_gets_lighter():
    assert trojans.libration_period(1.0, 1e-5) > trojans.libration_period(1.0, 1e-3)


def test_small_libration_is_a_tadpole_and_large_is_a_horseshoe():
    mu = 1e-6
    assert trojans.motion(10.0, mu) == "tadpole"
    assert trojans.motion(100.0, mu) == "horseshoe"
    assert trojans.motion(10.0, 0.05) is None
    assert trojans.motion(200.0, mu) is None


def test_only_l4_and_l5_hold_bodies():
    assert trojans.longitude_offset_deg("L4") == 60.0 and trojans.longitude_offset_deg("L5") == -60.0
    for point in ("L1", "L2", "L3"):
        with pytest.raises(trojans.TrojanError):
            trojans.longitude_offset_deg(point)
        with pytest.raises(trojans.TrojanError):
            trojans.point_label("Jupiter", point)


def test_the_map_label_names_the_secondary_and_the_point():
    assert trojans.point_label("Jupiter", "L4") == "Jupiter-L4"
    assert trojans.pair_text("Sun", "Jupiter", "L5") == "L5 of Sun and Jupiter"


def test_a_swarm_split_keeps_the_mass_and_the_draw_stays_in_its_range():
    low, high = trojans.LEADING_SHARE_RANGE
    for u in (0.0, 0.5, 1.0):
        lead = trojans.leading_share(u)
        assert low <= lead <= high
        l4, l5 = trojans.swarm_split(10.0, lead)
        assert math.isclose(l4 + l5, 10.0)
    with pytest.raises(trojans.TrojanError):
        trojans.swarm_split(1.0, 1.5)


def test_the_orbit_limits_and_allows():
    assert trojans.allows(SUN, JUPITER, EARTH, 0.05, 10.0)
    assert not trojans.allows(SUN, JUPITER, EARTH, 0.5, 10.0)
    assert not trojans.allows(SUN, JUPITER, EARTH, 0.05, 60.0)
    assert not trojans.allows(1.0, 0.5, EARTH)


def test_every_number_has_a_tier_and_a_tweak():
    for name, entry in trojans.NUMBERS.items():
        assert entry["tier"].startswith("T") and entry["tweak"].startswith("L"), name
