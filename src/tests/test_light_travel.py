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
                             length_unit="ly", epoch_unix=1000.0)
    position, light_time = apparent_position_of(star, (0.0, 0.0, 0.0))
    assert light_time == pytest.approx(YEAR / 1.1)
    assert position[0] == pytest.approx(LY / 1.1)


# --- VIEW.21: the same law on arrays -----------------------------------------------------------------------------------------------------

import numpy as np

from planetgen.physics import light_travel as lt

PC = constants.PARSEC_M


def _stars(offsets, velocities):
    return np.array(offsets, dtype=float), np.array(velocities, dtype=float)


def test_arrays_agree_with_the_one_body_function():
    rng = np.random.default_rng(7)
    offsets = rng.normal(size=(40, 3)) * 50 * LY
    velocities = rng.normal(size=(40, 3)) * 3e5
    elapsed = 1e9
    observer_velocity = (2e4, -1e4, 5e3)
    apparent, tau = lt.retarded_offsets(offsets, velocities, observer_velocity, elapsed)
    for index in range(40):
        position0 = tuple(offsets[index])
        velocity = tuple(velocities[index])
        # the one-body function with the observer fixed: shift the star by the observer's motion over the elapsed time
        star = uniform_motion(tuple(p - o * elapsed for p, o in zip(position0, observer_velocity)), velocity, 0.0)
        position, light_time = apparent_position(star, (0.0, 0.0, 0.0), now_s=elapsed)
        assert tau[index] == pytest.approx(light_time, rel=1e-9)
        assert apparent[index] == pytest.approx(position, rel=1e-9, abs=1e-3)


def test_a_star_at_rest_one_light_year_away_is_seen_a_year_ago_where_it_is():
    apparent, tau = lt.retarded_offsets([[LY, 0.0, 0.0]], [[0.0, 0.0, 0.0]], (0.0, 0.0, 0.0), 0.0)
    assert tau[0] == pytest.approx(YEAR) and apparent[0] == pytest.approx((LY, 0.0, 0.0))


def test_a_receding_star_is_seen_nearer_and_an_approaching_one_farther():
    offsets, velocities = _stars([[LY, 0.0, 0.0], [LY, 0.0, 0.0]], [[0.1 * C, 0.0, 0.0], [-0.1 * C, 0.0, 0.0]])
    apparent, tau = lt.retarded_offsets(offsets, velocities, (0.0, 0.0, 0.0), 0.0)
    assert apparent[0, 0] == pytest.approx(LY / 1.1) and tau[0] == pytest.approx(YEAR / 1.1)
    assert apparent[1, 0] == pytest.approx(LY / 0.9) and tau[1] == pytest.approx(YEAR / 0.9)


def test_a_hypervelocity_star_is_shifted_by_v_over_c_in_direction():
    """1,000 km/s across the line of sight at 100 pc: the plan's 11.5 arcminutes."""
    offsets, velocities = _stars([[100 * PC, 0.0, 0.0]], [[0.0, 1.0e6, 0.0]])
    apparent, _tau = lt.retarded_offsets(offsets, velocities, (0.0, 0.0, 0.0), 0.0)
    shift = math.atan2(abs(apparent[0, 1]), apparent[0, 0])
    assert math.degrees(shift) * 60 == pytest.approx(11.5, abs=0.1)
    assert apparent[0, 1] < 0.0         # seen behind its motion: where it was


def test_a_typical_star_is_shifted_by_about_twenty_arcseconds():
    offsets, velocities = _stars([[10 * PC, 0.0, 0.0]], [[0.0, 3.0e4, 0.0]])
    apparent, _tau = lt.retarded_offsets(offsets, velocities, (0.0, 0.0, 0.0), 0.0)
    assert math.degrees(math.atan2(abs(apparent[0, 1]), apparent[0, 0])) * 3600 == pytest.approx(20.6, abs=0.3)


def test_an_observer_moving_with_a_star_sees_it_where_it_is():
    """Co-moving: the star's offset never changes, whatever the elapsed time, and the light left it from where it is."""
    offsets, velocities = _stars([[3 * LY, 4 * LY, 0.0]], [[1.0e5, -2.0e4, 0.0]])
    apparent, tau = lt.retarded_offsets(offsets, velocities, (1.0e5, -2.0e4, 0.0), 5.0e11)
    # moving with the observer, the offset at the sky's time is the stored one less what the star moved during the light time
    assert apparent[0] == pytest.approx(offsets[0] - velocities[0] * tau[0], rel=1e-12)
    assert tau[0] == pytest.approx(math.dist(apparent[0], (0.0, 0.0, 0.0)) / C, rel=1e-9)


def test_two_passes_are_enough_for_slow_stars_and_more_run_for_fast_ones():
    offsets, velocities = _stars([[100 * PC, 0.0, 0.0]], [[0.0, 3.0e5, 0.0]])
    slow, _ = lt.retarded_offsets(offsets, velocities, (0.0, 0.0, 0.0), 0.0)
    assert slow[0, 0] == pytest.approx(100 * PC, rel=1e-3) and lt.MIN_PASSES == 2 and lt.MAX_PASSES > 3
    fast_offsets, fast_velocities = _stars([[LY, 0.0, 0.0]], [[0.5 * C, 0.0, 0.0]])
    fast, tau = lt.retarded_offsets(fast_offsets, fast_velocities, (0.0, 0.0, 0.0), 0.0)
    assert tau[0] == pytest.approx(YEAR / 1.5) and fast[0, 0] == pytest.approx(LY / 1.5)


def test_arrays_refuse_a_star_that_outruns_its_own_light():
    with pytest.raises(ValueError):
        lt.retarded_offsets([[LY, 0.0, 0.0]], [[-2 * C, 0.0, 0.0]], (0.0, 0.0, 0.0), 0.0)


def test_no_stars_give_empty_arrays():
    apparent, tau = lt.retarded_offsets(np.zeros((0, 3)), np.zeros((0, 3)), (0.0, 0.0, 0.0), 0.0)
    assert apparent.shape == (0, 3) and tau.shape == (0,)


def test_the_age_at_the_retarded_time_leaves_out_stars_not_yet_formed():
    age = np.array([5.0, 20.0, 40.0]) * YEAR
    tau = np.array([10.0, 10.0, 10.0]) * YEAR
    seen, formed, other = lt.seen_age(age, tau)
    assert list(formed) == [False, True, True] and seen[1] == pytest.approx(10 * YEAR) and not other.any()


def test_a_star_that_left_its_phase_within_the_light_time_is_seen_in_the_phase_before():
    age = np.array([100.0, 100.0, 100.0, 100.0]) * YEAR
    tau = np.array([30.0, 30.0, 30.0, 30.0]) * YEAR
    end = np.array([80.0, 60.0, np.nan, 120.0]) * YEAR       # ended 20 years ago, ended 40 years ago, none stored, ends later
    _seen, formed, other = lt.seen_age(age, tau, end)
    assert list(other) == [True, False, False, False] and formed.all()


def test_an_explosion_younger_than_the_light_time_has_not_been_seen():
    seen = lt.remnant_seen(np.array([500.0, 5000.0, 5000.0]) * YEAR, np.array([1000.0, 1000.0, 5000.0]) * YEAR)
    assert list(seen) == [False, True, True]


def test_straight_lines_are_exact_enough_for_a_star_ten_kiloparsecs_away():
    """The curvature bound: 0.003 pc (3e-7 radians, 0.07 arcsecond) for 220 km/s, 32,600 years of light and an 8.2 kpc orbit."""
    bound = lt.curvature_bound_m(2.2e5, 32600 * YEAR, 8200 * PC)
    assert bound / PC == pytest.approx(0.003, rel=0.2)
    assert bound / (10000 * PC) == pytest.approx(3.2e-7, rel=0.1)
    assert math.degrees(bound / (10000 * PC)) * 3600 < 0.1
