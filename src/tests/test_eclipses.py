# tests/test_eclipses.py

"""GEN.1462: the orbit-averaged eclipses of a star by a body of the system (`planetgen.physics.eclipses`), against textbook cases."""

import math

import pytest

from planetgen.physics import constants, eclipses

KM = constants.AU_TO_KM
SUN_R = 695_700.0


def test_the_overlap_of_two_disks():
    assert eclipses.lens_area(1.0, 1.0, 2.0) == 0.0
    assert eclipses.lens_area(1.0, 1.0, 0.0) == pytest.approx(math.pi)
    assert eclipses.lens_area(0.5, 2.0, 0.5) == pytest.approx(math.pi * 0.25)                 # a small disk inside a large one
    assert eclipses.lens_area(1.0, 1.0, 1.0) == pytest.approx(2.0 * math.pi / 3.0 - math.sqrt(3.0) / 2.0)


def test_the_moon_covers_a_small_share_of_the_suns_light_and_touches_it_rarely():
    """The Moon touches the Sun's disk for about 1.5e-4 of the time (the track within 0.54 degrees of the Sun in 6.7 percent of the passes) and covers 3.6e-5 of its light."""
    covered, touched = eclipses.orbiting_the_observer(SUN_R, 1.0, 1.0, 1737.4, 384_400.0, 27.3216 / 365.25, 5.145, KM)
    assert 2.0e-5 < covered < 6.0e-5
    assert 1.0e-4 < touched < 2.0e-4


def test_venus_crosses_the_sun_for_earth_four_times_in_243_years():
    covered, touched = eclipses.orbiting_the_star(SUN_R, 1.0, 1.0, 6051.8, 0.723, 0.6152, 3.39, KM)
    assert touched == pytest.approx(4 * 6.0 / (243.0 * 8766.0), rel=0.5)      # about six hours each
    assert 0.0 < covered < touched


def test_a_blocker_outside_the_observer_or_with_no_size_covers_nothing():
    assert eclipses.orbiting_the_star(SUN_R, 1.0, 1.0, 6371.0, 1.5, 1.8, 1.0, KM) == (0.0, 0.0)      # farther from the star than the observer
    assert eclipses.orbiting_the_star(SUN_R, 1.0, 1.0, 0.0, 0.5, 0.35, 1.0, KM) == (0.0, 0.0)
    assert eclipses.orbiting_the_observer(SUN_R, 1.0, 1.0, 100.0, 0.0, 0.1, 1.0, KM) == (0.0, 0.0)


def test_a_small_blocker_covers_in_proportion_to_its_area():
    base, _ = eclipses.orbiting_the_star(SUN_R, 5.2, 11.86, 3000.0, 1.0, 1.0, 0.0, KM)
    double, _ = eclipses.orbiting_the_star(SUN_R, 5.2, 11.86, 6000.0, 1.0, 1.0, 0.0, KM)
    assert double / base == pytest.approx(4.0, rel=0.03)


def test_a_tilted_orbit_misses_the_star_more_often():
    flat, flat_touched = eclipses.orbiting_the_observer(SUN_R, 1.0, 1.0, 1737.4, 384_400.0, 0.0748, 0.0, KM)
    tilted, tilted_touched = eclipses.orbiting_the_observer(SUN_R, 1.0, 1.0, 1737.4, 384_400.0, 0.0748, 5.0, KM)
    assert flat > 20.0 * tilted and flat_touched > 10.0 * tilted_touched


def test_two_stars_of_a_close_pair_eclipse_each_other_for_a_planet_in_their_plane():
    """Two Suns 0.2 AU apart (a 15.6-day pair), seen from 1 AU in their plane: each crosses the other twice a period, about 3 percent of the time touching."""
    covered, touched = eclipses.pair(SUN_R, SUN_R, 0.2, 0.0426, 1.0, 0.0, KM)
    assert 0.02 < touched < 0.04
    assert 0.005 < covered < 0.02
    tilted, tilted_touched = eclipses.pair(SUN_R, SUN_R, 0.2, 0.0426, 1.0, 3.0, KM)
    assert tilted < covered and tilted_touched < touched
    far, _ = eclipses.pair(SUN_R, SUN_R, 0.2, 0.0426, 1.0, 40.0, KM)
    assert 0.0 < far < covered / 20.0                         # a planet far out of their plane sees a crossing only when it stands near a node


def _pair_planet(distance_au=1.0, inclination_deg=0.0, moon=False, plane=(0.0, 0.0), stars=2):
    """A planet of a close pair of two Suns 0.2 AU apart (a 15.6-day pair), the host a stand-in for `BinaryStarProxy`."""
    from types import SimpleNamespace
    sun = SimpleNamespace(radius=SUN_R, luminosity=3.828e26)
    host = SimpleNamespace(stars=[sun, sun][:stars], binary_mutual_orbital_period_years=0.0426, binary_separation_au=0.2,
                           binary_mutual_orbital_inclination_deg=plane[0], binary_mutual_orbital_ascending_node_deg=plane[1])
    return SimpleNamespace(star=host, distance=distance_au, is_moon=moon, orbital_inclination_deg=inclination_deg, orbital_ascending_node_deg=0.0)


def test_a_planet_in_the_plane_of_a_close_pair_loses_the_share_its_stars_eclipse_of_the_light():
    """GEN.1484: two Suns 0.2 AU apart, seen from 1 AU in their plane, eclipse for about 3 percent of the time and cover 0.5 to 2 percent of the light of the star behind, so 0.25 to 1 percent of the pair's light."""
    factor = eclipses.light_factor_of(_pair_planet())
    covered, _ = eclipses.pair(SUN_R, SUN_R, 0.2, 0.0426, 1.0, 0.0, KM)
    assert factor == pytest.approx(1.0 - covered / 2.0, rel=1e-9)         # the star behind carries half the light of two equal stars
    assert 0.99 < factor < 0.998


def test_a_tilted_planet_loses_less_and_every_other_body_keeps_all_its_light():
    flat = eclipses.light_factor_of(_pair_planet())
    tilted = eclipses.light_factor_of(_pair_planet(inclination_deg=20.0))
    assert flat < tilted <= 1.0
    assert eclipses.light_factor_of(_pair_planet(moon=True)) is None
    assert eclipses.light_factor_of(_pair_planet(stars=1)) is None            # a host of one star
    lone = _pair_planet()
    lone.star = type(lone.star)(radius=SUN_R, luminosity=3.8e26)               # a Star holds no `stars` list
    assert eclipses.light_factor_of(lone) is None


def test_two_unequal_stars_lose_the_light_of_the_star_behind_in_proportion_to_its_luminosity():
    from types import SimpleNamespace
    big = SimpleNamespace(radius=SUN_R, luminosity=3.828e26)
    small = SimpleNamespace(radius=0.5 * SUN_R, luminosity=0.06 * 3.828e26)
    same = eclipses.pair_light_lost(big, big, 0.2, 0.0426, 1.0, 0.0, KM)
    unequal = eclipses.pair_light_lost(big, small, 0.2, 0.0426, 1.0, 0.0, KM)
    assert unequal > 0.0
    assert unequal < same                         # the faint star has little light to lose and the big one's disk is covered only by a small disk
    assert eclipses.pair_light_lost(big, small, 0.2, 0.0426, 1.0, 0.0, KM) == eclipses.pair_light_lost(small, big, 0.2, 0.0426, 1.0, 0.0, KM)


def test_gen1606_a_factor_of_zero_keeps_no_light_and_none_keeps_all():
    from planetgen.physics import eclipses
    assert eclipses.kept_share(None) == 1.0 and eclipses.kept_share(0.0) == 0.0 and eclipses.kept_share(0.8) == 0.8
