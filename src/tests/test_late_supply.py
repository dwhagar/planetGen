"""GEN.1493 slice 0: the late supply a planet takes from the reservoirs of its system (`planetgen/physics/late_supply.py`).
The Sun's bodies below are test targets only: the law must be able to reach the published veneer of an Earth, and nothing in the module is centered on them."""

import pytest

from planetgen.physics import influx, late_supply

EARTH_KG = 5.9722e24
# (radius_km, distance_au, mass_kg): Mercury, Venus, Earth, Mars, Jupiter, Saturn
SYSTEM = [(2440, 0.387, 3.3e23), (6052, 0.723, 4.87e24), (6371, 1.0, 5.97e24), (3390, 1.524, 6.42e23), (69911, 5.2, 1.898e27), (58232, 9.5, 5.683e26)]
EARTH = SYSTEM[2]
BELT = influx.Reservoir("belt", 3.0 * EARTH_KG, 2.1, 3.3)


def test_a_planet_with_no_belt_takes_nothing():
    assert late_supply.late_supply_kg(EARTH, SYSTEM, 1.0, []) == 0.0
    assert late_supply.veneer_fraction(EARTH, SYSTEM, 1.0, []) == 0.0


def test_the_supply_is_linear_in_the_belt_mass_and_the_inward_share(monkeypatch):
    base = late_supply.late_supply_kg(EARTH, SYSTEM, 1.0, [BELT])
    heavier = influx.Reservoir("belt", 6.0 * EARTH_KG, 2.1, 3.3)
    assert late_supply.late_supply_kg(EARTH, SYSTEM, 1.0, [heavier]) == pytest.approx(2.0 * base)
    monkeypatch.setattr(late_supply, "INWARD_SHARE", 2.0 * late_supply.INWARD_SHARE)
    assert late_supply.late_supply_kg(EARTH, SYSTEM, 1.0, [BELT]) == pytest.approx(2.0 * base)


def test_a_giant_next_to_the_belt_takes_most_of_what_would_reach_the_planets_inside_it():
    with_giant = late_supply.late_supply_kg(EARTH, SYSTEM, 1.0, [BELT])
    without = late_supply.late_supply_kg(EARTH, [b for b in SYSTEM if b[2] < 1.0e26], 1.0, [BELT])
    assert with_giant < 0.05 * without


def test_a_giant_far_beyond_the_belt_does_not_share_it():
    far = [*SYSTEM[:4], (69911, 40.0, 1.898e27)]
    near_free = late_supply.late_supply_kg(EARTH, SYSTEM[:4], 1.0, [BELT])
    assert late_supply.late_supply_kg(EARTH, far, 1.0, [BELT]) == pytest.approx(near_free)


def test_a_planet_beyond_the_belt_is_not_reached():
    outside = (6371, 10.0, 5.97e24)
    assert late_supply.late_supply_kg(outside, [*SYSTEM, outside], 1.0, [BELT]) == 0.0


def test_a_lone_planet_takes_the_whole_inward_share_wherever_it_lies_in_or_below_the_ring():
    for au in (1.0, 2.0, 2.7):
        lone = (6371, au, 5.97e24)
        assert late_supply.late_supply_kg(lone, [lone], 1.0, [BELT]) == pytest.approx(late_supply.INWARD_SHARE * BELT.mass_kg)


def test_the_shares_of_all_bodies_never_exceed_what_was_scattered():
    total = sum(late_supply.late_supply_kg(body, SYSTEM, 1.0, [BELT]) for body in SYSTEM)
    assert total <= late_supply.INWARD_SHARE * BELT.mass_kg * (1.0 + 1e-9)


def test_a_heavier_belt_can_give_an_earth_the_published_veneer():
    """The law can reach 0.2 to 0.7 percent of an Earth's mass with an inward share and a belt mass in the ranges a belt can have (not a calibration to the Sun)."""
    fractions = [late_supply.veneer_fraction(EARTH, SYSTEM, 1.0, [influx.Reservoir("belt", mass * EARTH_KG, 2.1, 3.3)]) for mass in (3.0, 10.0, 30.0)]
    assert any(0.002 <= f <= 0.007 for f in fractions)
