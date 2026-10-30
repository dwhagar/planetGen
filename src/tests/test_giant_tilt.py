# tests/test_giant_tilt.py

"""The giant-impact tilt of an ice giant (GEN.1053): Uranus-like and Neptune-like cases, Kegerreis et al. 2018 angular momenta."""

import math

from planetgen.physics import giant_tilt as gt

URANUS_M, URANUS_R = 14.5, 4.0         # Earth masses and radii


def test_the_kegerreis_impact_brings_an_angular_momentum_in_their_range():
    # 2 Earth-mass impactor on 14.536 - 2 Earth masses at 5 km/s, grazing; their L runs from 1 to 10e36 kg m^2/s
    grazing = gt.impact_angular_momentum(12.536, 3.8, 2.0, 1.0)
    assert 1e36 < grazing < 1e37
    # and a head-on hit brings none; the angular momentum grows with the impact parameter
    assert gt.impact_angular_momentum(12.536, 3.8, 2.0, 0.0) == 0.0
    assert gt.impact_angular_momentum(12.536, 3.8, 2.0, 0.5) < grazing
    # the impactor's share of the mass raises it
    assert gt.impact_angular_momentum(13.536, 3.9, 1.0, 1.0) < grazing < gt.impact_angular_momentum(11.536, 3.7, 3.0, 1.0)


def test_uranus_is_a_two_earth_mass_grazing_hit_with_a_tilt_over_ninety_possible():
    ratio = gt.angular_momentum_ratio(URANUS_M, URANUS_R, 2.0, 0.8)
    assert ratio > 1.0
    angles = [gt.tilt_deg(ratio, -1.0 + 0.1 * i) for i in range(21)]
    assert max(angles) > 90.0
    # the 98 degree is reached at some direction
    assert any(abs(angle - 98.0) < 3.0 for angle in angles)
    # the spin is faster than the present 17.24 h in the Kegerreis case (L over 2e36), and the planet stays under breakup
    assert gt.final_period_h(URANUS_M, URANUS_R, ratio, 0.0) < 17.24
    assert gt.final_period_h(URANUS_M, URANUS_R, 1e6, 1.0) > 2.0


def test_neptune_is_a_small_impactor_a_tilt_under_asin_of_the_ratio():
    # a hit that brings under the planet's spin turns the axis by at most asin(ratio)
    ratio = gt.angular_momentum_ratio(17.1, 3.9, 0.3, 0.5)
    assert ratio < 1.0
    biggest = max(gt.tilt_deg(ratio, -1.0 + 0.01 * i) for i in range(201))
    assert math.isclose(biggest, math.degrees(math.asin(ratio)), abs_tol=0.5)
    assert biggest < 60.0
    assert gt.probability_of_tilt_over(ratio) == 0.0


def test_the_tilt_geometry():
    assert gt.tilt_deg(0.0, 0.3) == 0.0
    assert gt.tilt_deg(1.0, 1.0) == 0.0                   # along the axis: no turn
    assert math.isclose(gt.tilt_deg(1.0, 0.0), 45.0, abs_tol=1e-9)
    assert math.isclose(gt.tilt_deg(1e9, 0.0), 90.0, abs_tol=1e-3)     # a huge perpendicular hit: the axis follows it
    assert gt.tilt_deg(3.0, -1.0) > 179.0                 # a big impactor against the spin reverses it
    # more angular momentum, a larger chance of a large tilt
    chances = [gt.probability_of_tilt_over(r) for r in (0.5, 1.0, 2.0, 5.0, 20.0)]
    assert chances == sorted(chances) and chances[0] == 0.0 and chances[-1] > 0.6


def test_the_share_of_mergers_that_tilt_over_sixty_is_the_stand_in_and_the_chance_one_in_five():
    share = gt.share_tilted_over(URANUS_M, URANUS_R)
    assert 0.3 < share < 0.42
    chance = gt.large_tilt_chance(URANUS_M, URANUS_R)
    assert math.isclose(chance, 0.2, abs_tol=0.02)         # the Solar System: one ice giant in five
    assert gt.large_tilt_chance(URANUS_M, URANUS_R, has_merger=False) == 0.0
    # about the same across the 5 to 30 Earth mass range (the design table)
    for mass, radius in ((6.0, 2.6), (25.0, 4.8)):
        assert 0.25 < gt.share_tilted_over(mass, radius) < 0.45


def test_the_draw_maps_unit_numbers_to_a_tilt():
    assert 0.0 <= gt.draw_tilt_deg(URANUS_M, URANUS_R, 0.5, 0.5, 0.5) <= 180.0
    # a head-on hit (u of the parameter 0) leaves the axis alone
    assert gt.draw_tilt_deg(URANUS_M, URANUS_R, 0.5, 0.0, 0.2) == 0.0
    # the impactor mass runs over its share of the planet's
    assert math.isclose(gt.impactor_mass_earth(URANUS_M, 0.0), 0.02 * URANUS_M)
    assert math.isclose(gt.impactor_mass_earth(URANUS_M, 1.0), 0.25 * URANUS_M)


def test_the_rings_and_moons_lie_in_the_equator():
    assert gt.ring_plane_tilt_deg(97.8) == 97.8
