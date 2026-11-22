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


def test_a_giant_between_the_ring_and_the_planets_takes_what_would_reach_them():
    """Jupiter at 5.2 AU sits between a ring at 20 to 30 AU and the inner planets: it takes the planetesimals that cross its orbit within the residence time, and throws a share of the rest out."""
    ring = late_supply.Early("belt", 10.0 * EARTH_KG, 20.0, 30.0)
    with_giant = late_supply.late_supply_kg(EARTH, SYSTEM, 1.0, [ring])
    without = late_supply.late_supply_kg(EARTH, [b for b in SYSTEM if b[2] < 1.0e26], 1.0, [ring])
    assert without > 0.0 and with_giant < 1.0e-3 * without


def test_a_giant_beyond_the_ring_does_not_take_part_in_it():
    far = [*SYSTEM[:4], (69911, 40.0, 1.898e27)]
    near_free = late_supply.late_supply_kg(EARTH, SYSTEM[:4], 1.0, [BELT])
    assert late_supply.late_supply_kg(EARTH, far, 1.0, [BELT]) == pytest.approx(near_free)


def test_a_planet_beyond_the_belt_is_not_reached():
    outside = (6371, 10.0, 5.97e24)
    assert late_supply.late_supply_kg(outside, [*SYSTEM, outside], 1.0, [BELT]) == 0.0


def test_a_lone_planet_takes_the_share_its_cross_section_gives_of_the_part_of_the_ring_beyond_it():
    for au in (1.0, 2.0, 2.7):
        lone = (6371, au, 5.97e24)
        reached = late_supply.influx.reach(BELT, au)
        expected = late_supply.INWARD_SHARE * BELT.mass_kg * reached * late_supply.hit_share(*lone, 1.0)
        assert late_supply.late_supply_kg(lone, [lone], 1.0, [BELT]) == pytest.approx(expected)


def test_everything_scattered_is_taken_ejected_or_left_to_the_star():
    taken, ejected, star = late_supply.fates(BELT, SYSTEM, 1.0)
    assert sum(taken) + ejected + star == pytest.approx(late_supply.INWARD_SHARE * BELT.mass_kg)
    assert all(share >= 0.0 for share in taken) and ejected >= 0.0 and star >= 0.0


def test_a_planetesimal_that_is_not_taken_goes_on_to_the_planets_inside():
    """Two planets in a row: the inner one takes from what the outer one let past, so a heavier outer planet leaves it less."""
    inner = (6371, 1.0, 5.97e24)
    light, heavy = (3390, 1.5, 6.42e23), (25362, 1.5, 8.68e25)
    ring = late_supply.Early("belt", 3.0 * EARTH_KG, 2.1, 3.3)
    assert late_supply.late_supply_kg(inner, [inner, light], 1.0, [ring]) > late_supply.late_supply_kg(inner, [inner, heavy], 1.0, [ring])


def test_the_share_a_body_takes_grows_with_the_years_and_the_share_a_giant_throws_out_with_its_safronov_number():
    earth, jupiter = (6371, 1.0, 5.97e24), (69911, 5.2, 1.898e27)
    short = late_supply.hit_share(*earth, 1.0)
    original = late_supply.RESIDENCE_YR
    try:
        late_supply.RESIDENCE_YR = original / 10.0
        assert late_supply.hit_share(*earth, 1.0) < short
    finally:
        late_supply.RESIDENCE_YR = original
    assert late_supply.eject_share(*jupiter, 1.0) > 0.5 > 0.05 > late_supply.eject_share(*earth, 1.0)
    assert late_supply.hit_share(*jupiter, 1.0) > 0.99
    assert late_supply.hit_share(0.0, 1.0, 1.0, 1.0) == 0.0 and late_supply.eject_share(0.0, 1.0, 1.0, 1.0) == 0.0


def test_the_leftovers_of_the_zone_give_an_earth_the_published_veneer_with_the_giants_beyond_it():
    """The Sun's terrestrial planets grew 1.8 Earth masses between 0.58 and 1.83 AU; the leftover share of that (row LA3) puts 0.2 to 0.7 percent of an Earth's mass on the Earth, giants or not."""
    from planetgen.generation import late_accretion

    zone = late_supply.Early("leftover", late_accretion.LEFTOVER_SHARE * 1.8 * EARTH_KG, 0.58, 1.83, inward_share=1.0)
    assert 0.002 <= late_supply.veneer_fraction(EARTH, SYSTEM, 1.0, [zone]) <= 0.007
    assert 0.002 <= late_supply.veneer_fraction(EARTH, SYSTEM[:4], 1.0, [zone]) <= 0.007


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


def test_a_solar_system_like_architecture_with_a_swept_outer_disk_gives_the_planets_inside_the_giants_almost_nothing():
    """An outer disk the outermost giant swept sends planetesimals on crossing orbits, but the four giants between it and the Earth take or throw out nearly all of them: the Earth's veneer comes from the leftovers of its own zone (row LA3)."""
    giants = [(69911, 5.2, 1.898e27), (58232, 9.5, 5.683e26), (25362, 19.2, 8.68e25), (24622, 30.1, 1.02e26)]
    system = [*SYSTEM[:4], *giants]
    fates = scatter.scatter_shares(1.0, 17.1, 30.1, 3.88, 0.0, 30.0, 50.0, 0.0)
    share = late_supply.inward_share_of(fates, False)
    disk = late_supply.Early("belt", 20.0 * EARTH_KG, 30.0, 50.0, inward_share=share)
    assert late_supply.veneer_fraction(EARTH, system, 1.0, [disk]) < 1.0e-4
