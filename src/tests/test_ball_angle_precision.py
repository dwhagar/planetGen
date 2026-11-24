# tests/test_ball_angle_precision.py

"""
GEN.376: the angles of the sphere walks (`galaxy/region.py` through `spherical.sphere_shells` and `band_window`, and `window.cells_in_ball`) are taken by the haversine
form below `spherical.SMALL_ANGLE` (1e-3 rad): the arccos of a cosine within a few ulps of 1 keeps hardly a digit of the angle, which at 262 kpc is the angle of a ball of a few
parsecs.
"""

import math

import pytest

from planetgen.galaxy import spherical, window


def _cosine(radius, reach, ball):
    return (radius * radius + reach * reach - ball * ball) / (2.0 * radius * reach)


def test_two_points_a_few_parsecs_apart_at_262_kpc_give_the_angle_to_a_tight_bound():
    radius = 262000.0
    ball = 3.0
    exact = 2.0 * math.asin(ball / (2.0 * radius))          # the chord of the ball on the shell, both at the same distance
    cosine = _cosine(radius, radius, ball)
    assert abs(math.acos(cosine) - exact) / exact > 1e-8    # the arccos form the walks used keeps about six digits
    assert spherical.ball_angle(radius, radius, ball, cosine) == pytest.approx(exact, rel=1e-12)


def test_the_haversine_form_holds_with_the_two_distances_apart():
    radius, reach, ball = 262000.0, 262002.0, 3.0
    gap = reach - radius
    exact = 2.0 * math.asin(math.sqrt((ball - gap) * (ball + gap) / (4.0 * radius * reach)))
    assert spherical.ball_angle(radius, reach, ball, _cosine(radius, reach, ball)) == pytest.approx(exact, rel=1e-12)


def test_a_ball_that_only_grazes_the_shell_gives_a_zero_angle_not_an_error():
    radius, reach, ball = 1000.0, 1003.0, 3.0
    assert spherical.ball_angle(radius, reach, ball + 1e-12, _cosine(radius, reach, ball + 1e-12)) < 1e-6
    assert spherical.ball_angle(radius, reach, ball, 1.0 + 1e-12) == 0.0


@pytest.mark.parametrize("angle", [0.9e-3, 1.1e-3])
def test_the_two_forms_agree_at_the_switch(angle):
    radius = reach = 1000.0
    ball = 2.0 * radius * math.sin(angle / 2.0)
    assert spherical.ball_angle(radius, reach, ball, _cosine(radius, reach, ball)) == pytest.approx(angle, rel=1e-9)


def test_a_large_angle_is_still_the_arccos_and_a_whole_shell_is_pi():
    assert spherical.ball_angle(10.0, 12.0, 8.0, _cosine(10.0, 12.0, 8.0)) == math.acos(_cosine(10.0, 12.0, 8.0))
    assert spherical.ball_angle(10.0, 12.0, 8.0, -1.5) == math.pi


def test_the_longitude_half_width_on_the_centers_own_parallel_follows_the_haversine():
    delta = 1.1450381679451866e-05
    colat = 0.9
    cosine = (math.cos(delta) - math.cos(colat) ** 2) / math.sin(colat) ** 2     # the law of cosines on the same parallel
    exact = 2.0 * math.asin(math.sin(delta / 2.0) / math.sin(colat))
    assert abs(math.acos(cosine) - exact) / exact > 1e-8
    assert spherical.longitude_half_width(delta, colat, colat, cosine) == pytest.approx(exact, rel=1e-12)
    assert spherical.longitude_half_width(2.0, 1.0, 0.5, -2.0) == math.pi


def test_a_ball_of_three_parsecs_at_262_kpc_lists_the_cells_within_it():
    center = (262000.0, 0.0, 0.0)
    cells = window.cells_in_ball(1, center, 3.0, 1.0)
    for cell in cells:
        x, y, z = window.cell_position_pc(cell, 1.0)
        assert math.dist((x, y, z), center) <= 3.0 + 1e-6
