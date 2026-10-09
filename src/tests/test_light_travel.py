# tests/test_light_travel.py

"""VIEW.5: light-travel positions, against hand cases."""

import math

import pytest

from planetgen.physics import constants
from planetgen.physics.light_travel import apparent_position, apparent_position_of, uniform_motion
from planetgen.physics.position import SpatialPosition3D

C = constants.SPEED_OF_LIGHT_M_S
LY = constants.LY_TO_M
YEAR = LY / C  # a light year is a year of light, by the definition


def test_a_star_one_light_year_away_at_rest_is_seen_a_year_ago_where_it_is():
    position, light_time = apparent_position(uniform_motion((LY, 0.0, 0.0), (0.0, 0.0, 0.0)), (0.0, 0.0, 0.0))
    assert light_time == pytest.approx(YEAR)
    assert position == pytest.approx((LY, 0.0, 0.0))


def test_a_star_moving_away_is_seen_nearer_than_it_is():
    # 1 ly out and receding at 0.1 c: the light left when it was
    # x0 - v*tau, with c*tau = x0 - v*tau, so tau = x0 / (c + v) = 1 / 1.1 years.
    star = uniform_motion((LY, 0.0, 0.0), (0.1 * C, 0.0, 0.0))
    position, light_time = apparent_position(star, (0.0, 0.0, 0.0))
    assert light_time == pytest.approx(YEAR / 1.1)
    assert position[0] == pytest.approx(LY / 1.1)
    # Now it is at 1 ly, but it appears at 0.909 ly.
    assert star(0.0)[0] == pytest.approx(LY)


def test_a_star_moving_toward_us_is_seen_farther_than_it_is():
    position, light_time = apparent_position(uniform_motion((LY, 0.0, 0.0), (-0.1 * C, 0.0, 0.0)), (0.0, 0.0, 0.0))
    assert light_time == pytest.approx(YEAR / 0.9)
    assert position[0] == pytest.approx(LY / 0.9)


def test_sideways_motion_and_a_moving_observer_point():
    # A star 3 ly east of the observer and 4 ly north moving north at 0.5 c.
    star = uniform_motion((3 * LY, 4 * LY, 0.0), (0.0, 0.5 * C, 0.0))
    position, light_time = apparent_position(star, (0.0, 0.0, 0.0))
    # |(3, 4 - 0.5 t, 0)| = t (years, ly): 9 + (4 - t/2)^2 = t^2.
    t = light_time / YEAR
    assert 9 + (4 - 0.5 * t) ** 2 == pytest.approx(t ** 2)
    assert math.dist(position, (0.0, 0.0, 0.0)) == pytest.approx(light_time * C)
    assert position[1] == pytest.approx((4 - 0.5 * t) * LY)


def test_a_body_here_is_seen_where_it_is_now():
    position, light_time = apparent_position(uniform_motion((0.0, 0.0, 0.0), (1e4, 0.0, 0.0)), (0.0, 0.0, 0.0))
    assert light_time == 0.0 and position == (0.0, 0.0, 0.0)


def test_faster_than_light_toward_the_observer_never_settles():
    with pytest.raises(ValueError):
        apparent_position(uniform_motion((LY, 0.0, 0.0), (-2 * C, 0.0, 0.0)), (0.0, 0.0, 0.0))


def test_a_spatial_position_in_light_years_is_seen_a_year_ago():
    star = SpatialPosition3D((1.0, 0.0, 0.0), (0.0, 0.0, 0.0), velocity_vector_cartesian=(0.1 * C, 0.0, 0.0),
                             length_unit_m=LY, epoch_unix=1000.0)
    position, light_time = apparent_position_of(star, (0.0, 0.0, 0.0))
    assert light_time == pytest.approx(YEAR / 1.1)
    assert position[0] == pytest.approx(LY / 1.1)
