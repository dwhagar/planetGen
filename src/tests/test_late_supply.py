"""GEN.1493 slice 0: the late supply a planet takes from the reservoirs of its system (`planetgen/physics/late_supply.py`).
The Sun's bodies below are test targets only: the law must be able to reach the published veneer of an Earth, and nothing in the module is centered on them."""

import pytest

from planetgen.physics import late_supply, scatter

EARTH_KG = 5.9722e24
# (radius_km, distance_au, mass_kg): Mercury, Venus, Earth, Mars, Jupiter, Saturn
SYSTEM = [(2440, 0.387, 3.3e23), (6052, 0.723, 4.87e24), (6371, 1.0, 5.97e24), (3390, 1.524, 6.42e23), (69911, 5.2, 1.898e27), (58232, 9.5, 5.683e26)]
EARTH = SYSTEM[2]
BELT = late_supply.Early("belt", 3.0 * EARTH_KG, 2.1, 3.3)


def test_a_planet_with_no_belt_takes_nothing():
    assert late_supply.late_supply_kg(EARTH, SYSTEM, 1.0, []) == 0.0
    assert late_supply.veneer_fraction(EARTH, SYSTEM, 1.0, []) == 0.0


def test_the_supply_is_linear_in_the_belt_mass_and_the_inward_share():
    base = late_supply.late_supply_kg(EARTH, SYSTEM, 1.0, [BELT])
    heavier = late_supply.Early("belt", 6.0 * EARTH_KG, 2.1, 3.3)
    assert late_supply.late_supply_kg(EARTH, SYSTEM, 1.0, [heavier]) == pytest.approx(2.0 * base)
    richer = late_supply.Early("belt", 3.0 * EARTH_KG, 2.1, 3.3, inward_share=2.0 * late_supply.INWARD_SHARE)
    assert late_supply.late_supply_kg(EARTH, SYSTEM, 1.0, [richer]) == pytest.approx(2.0 * base)


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
    fractions = [late_supply.veneer_fraction(EARTH, SYSTEM, 1.0, [late_supply.Early("belt", mass * EARTH_KG, 2.1, 3.3)]) for mass in (3.0, 10.0, 30.0)]
    assert any(0.002 <= f <= 0.007 for f in fractions)


def _fates(**shares):
    """The populations of an outer disk (`scatter.POPULATIONS`) with the shares given; the rest is the cold classical belt."""
    fates = dict.fromkeys(scatter.POPULATIONS, 0.0)
    fates.update(shares)
    fates["cold_classical"] = 1.0 - sum(share for key, share in fates.items() if key != "cold_classical")
    return fates


def test_a_ring_no_giant_reached_sends_the_base_share_and_a_stirred_one_more():
    assert late_supply.inward_share_of(None, None) == late_supply.INWARD_SHARE
    assert late_supply.inward_share_of(_fates(), False) == pytest.approx(late_supply.INWARD_SHARE)
    assert late_supply.inward_share_of(_fates(), True) == pytest.approx(late_supply.INWARD_SHARE + late_supply.STIRRED_GAIN)


def test_a_disk_a_giant_threw_onto_crossing_orbits_sends_more_inward_than_a_quiet_one():
    quiet = late_supply.inward_share_of(_fates(), False)
    crossing = late_supply.inward_share_of(_fates(scattered_disk=0.4, centaurs=0.1), False)
    assert crossing > quiet
    assert crossing == pytest.approx(late_supply.INWARD_SHARE * 0.5 + late_supply.CROSSING_GAIN * 0.5)


def test_a_giant_that_threw_the_disk_out_of_the_system_sends_less_inward_than_one_that_kept_it_bound():
    kept_bound = late_supply.inward_share_of(_fates(scattered_disk=0.4, centaurs=0.1), False)
    thrown_out = late_supply.inward_share_of(_fates(ejected=0.5), False)
    assert thrown_out < kept_bound


def test_a_disk_swept_clean_still_feeds_the_planets_inside_it_through_what_it_threw_inward():
    fates = _fates(scattered_disk=0.6, centaurs=0.3, ejected=0.1)
    assert late_supply.inward_share_of(fates, False) == pytest.approx(late_supply.CROSSING_GAIN * 0.9)
    gone = late_supply.Early("belt", 5.0 * EARTH_KG, 30.0, 50.0, inward_share=late_supply.inward_share_of(fates, False))
    lone_giant_free = [(6371, 1.0, 5.97e24)]
    assert late_supply.late_supply_kg(EARTH, lone_giant_free, 1.0, [gone]) > 0.0


def test_the_share_never_passes_one():
    assert late_supply.inward_share_of(_fates(scattered_disk=1.0), True) <= 1.0


def test_a_lost_stretch_sends_only_what_the_giants_threw_onto_crossing_orbits():
    assert late_supply.inward_share_of(_fates(), True, gone=True) == 0.0
    assert late_supply.inward_share_of(None, False, gone=True) == 0.0
    assert late_supply.inward_share_of(_fates(scattered_disk=0.4, centaurs=0.1, ejected=0.2), True, gone=True) == pytest.approx(late_supply.CROSSING_GAIN * 0.5)


def test_a_solar_system_like_architecture_with_a_swept_outer_disk_stays_inside_the_published_veneer():
    """The two gains are fitted (rows LS1a, LS1b) so a synthetic Sun-like system, an Earth inside the four giants with an outer disk of 10 to 20 Earth masses the outermost giant swept, keeps 0.2 to 0.7 percent of its mass."""
    giants = [(69911, 5.2, 1.898e27), (58232, 9.5, 5.683e26), (25362, 19.2, 8.68e25), (24622, 30.1, 1.02e26)]
    system = [*SYSTEM[:4], *giants]
    fates = scatter.scatter_shares(1.0, 17.1, 30.1, 3.88, 0.0, 30.0, 50.0, 0.0)
    share = late_supply.inward_share_of(fates, False)
    for mass in (10.0, 20.0):
        disk = late_supply.Early("belt", mass * EARTH_KG, 30.0, 50.0, inward_share=share)
        assert 0.002 <= late_supply.veneer_fraction(EARTH, system, 1.0, [disk]) <= 0.007
