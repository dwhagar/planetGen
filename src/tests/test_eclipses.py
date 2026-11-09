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


def test_venus_crosses_the_sun_for_earth_eight_times_in_243_years():
    covered, touched = eclipses.orbiting_the_star(SUN_R, 1.0, 1.0, 6051.8, 0.723, 0.6152, 3.39, KM)
    assert touched == pytest.approx(8 * 6.0 / (243.0 * 8766.0), rel=0.5)      # about six hours each
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
