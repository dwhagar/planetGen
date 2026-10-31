# tests/test_lagrange_geometry.py

"""MAP.354 geometry: the five Lagrange points of a pair in the system frame, and the labels of the occupied ones."""

import math

import pytest

from planetgen.physics import lagrange_geometry as lg
from planetgen.physics import trojans


def test_earth_moon_collinear_points_match_the_published_values():
    x1, x2, x3 = lg.collinear_offsets(0.01215)
    assert x1 == pytest.approx(0.8369, abs=2e-4)
    assert x2 == pytest.approx(1.1557, abs=2e-4)
    assert x3 == pytest.approx(-1.0051, abs=2e-4)


def test_l4_and_l5_make_equilateral_triangles_with_the_pair():
    pts = lg.points([0, 0, 0], 1.0, [5.2, 0, 0], 9.5e-4, [0, 13.0, 0])
    for name in ("L4", "L5"):
        p = pts[name]
        assert math.dist(p, [0, 0, 0]) == pytest.approx(5.2, rel=1e-9)
        assert math.dist(p, [5.2, 0, 0]) == pytest.approx(5.2, rel=1e-9)
    assert pts["L4"][1] > 0 > pts["L5"][1]


def test_l4_leads_the_secondary_for_either_sense_of_motion():
    ahead = lg.points([0, 0, 0], 1.0, [1, 0, 0], 1e-3, [0, 1, 0])
    behind = lg.points([0, 0, 0], 1.0, [1, 0, 0], 1e-3, [0, -1, 0])
    assert ahead["L4"][1] > 0 and behind["L4"][1] < 0


def test_the_collinear_points_lie_on_the_line_in_order():
    pts = lg.points([0, 0, 0], 1.0, [5.2, 0, 0], 9.5e-4, [0, 13.0, 0])
    assert pts["L3"][0] < 0 < pts["L1"][0] < 5.2 < pts["L2"][0]
    assert 5.2 - pts["L1"][0] == pytest.approx(5.2 * (9.5e-4 / 3) ** (1 / 3), rel=0.03)
    assert all(abs(pts[k][1]) < 1e-9 for k in ("L1", "L2", "L3"))


def test_a_tilted_orbit_puts_l4_and_l5_in_its_plane():
    pts = lg.points([0, 0, 0], 1.0, [1, 0, 0], 1e-3, [0, 0.7, 0.7])
    for k in ("L4", "L5"):
        assert pts[k][1] == pytest.approx(pts[k][2], abs=1e-9) or pts[k][1] == pytest.approx(-pts[k][2], abs=1e-9)
        assert math.dist(pts[k], [0, 0, 0]) == pytest.approx(1.0, rel=1e-9)


def test_occupied_points_carry_the_map_label():
    pts = lg.points([0, 0, 0], 1.0, [5.2, 0, 0], 9.5e-4, [0, 13.0, 0])
    found = lg.occupied(pts, "Jupiter", {"L4", "L5"})
    assert [f["label"] for f in found] == ["Jupiter-L4", "Jupiter-L5"]
    assert lg.occupied(pts, "Jupiter", set()) == []
    assert found[0]["label"] == trojans.point_label("Jupiter", "L4")


def test_a_bad_mass_parameter_or_a_zero_motion_is_refused():
    with pytest.raises(trojans.TrojanError):
        lg.collinear_offsets(0.7)
    with pytest.raises(trojans.TrojanError):
        lg.points([0, 0, 0], 1.0, [1, 0, 0], 1e-3, [1, 0, 0])
