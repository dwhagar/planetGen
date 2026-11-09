# tests/test_day_union.py

"""GEN.1461: the day of a planet with more than one sun is the union of the arcs every counted sun lights (`physics/daylight.py`)."""

import math
from types import SimpleNamespace

import pytest

from planetgen.physics import constants, daylight, radiation_field as rf

SUN_RADIUS_KM = 695_700.0
SUN_W = constants.SOLAR_LUMINOSITY


def test_arcs_that_overlap_are_counted_once():
    quarter = math.pi / 2.0
    assert daylight.lit_arcs_share([(0.0, quarter)]) == pytest.approx(0.5)
    assert daylight.lit_arcs_share([(0.0, quarter), (math.pi, quarter)]) == pytest.approx(1.0)       # opposite suns: never night
    assert daylight.lit_arcs_share([(0.0, quarter), (quarter, quarter)]) == pytest.approx(0.75)
    assert daylight.lit_arcs_share([(0.1, quarter), (0.0, quarter)]) == pytest.approx(0.5 + 0.1 / (2.0 * math.pi))
    assert daylight.lit_arcs_share([(6.2, 1.0), (0.1, 1.0)]) == pytest.approx((2.0 + (0.1 + 2.0 * math.pi - 6.2)) / (2.0 * math.pi))      # across the seam
    assert daylight.lit_arcs_share([(0.0, math.pi)]) == 1.0 and daylight.lit_arcs_share([]) == 0.0


def _two_suns(offset_au, distance_au=1.0, tilt_rad=0.0):
    sun = SimpleNamespace(radius=SUN_RADIUS_KM)
    return sun, [(sun, None), (sun, offset_au)]


def test_a_wide_pair_seen_opposite_lights_a_longer_day_than_one_sun():
    """Hand worked: with no tilt every sun stands on the equator, an arc of half a turn each; the host goes round the sky once a year against a companion fixed
    opposite its start, so the union is half a turn plus half the angle between them, which averages a quarter more over the year (a surface mean as well)."""
    _, counted = _two_suns((-1.0e5, 0.0, 0.0))
    equator, mean = daylight.sky_day_shares(0.0, 1.0, counted)
    assert equator == pytest.approx(0.75, abs=0.01)
    assert mean == pytest.approx(0.75, abs=0.01)
    assert equator > daylight.equator_day_share(0.0, daylight.angular_radius(SUN_RADIUS_KM, 1.0)) + 0.2


def test_a_close_pair_seen_from_far_away_is_nearly_one_sun():
    _, counted = _two_suns((0.001, 0.0, 0.0), distance_au=1.0)
    equator, mean = daylight.sky_day_shares(math.radians(23.4), 1.0, counted)
    one = daylight.equator_day_share(23.4, daylight.angular_radius(SUN_RADIUS_KM, 1.0))
    assert equator == pytest.approx(one, abs=0.002)
    assert mean == pytest.approx(daylight.mean_day_share(daylight.angular_radius(SUN_RADIUS_KM, 1.0)), abs=0.002)


def test_more_suns_never_shorten_the_day():
    _, counted = _two_suns((-30.0, 20.0, 5.0))
    equator, mean = daylight.sky_day_shares(math.radians(30.0), 1.0, counted)
    assert equator >= daylight.equator_day_share(30.0, daylight.angular_radius(SUN_RADIUS_KM, 1.0)) - 1e-9
    assert mean >= daylight.mean_day_share(daylight.angular_radius(SUN_RADIUS_KM, 1.0)) - 1e-9


def _planet_with(companion_w, offset, axis_au=1.0e4):
    host = SimpleNamespace(luminosity=SUN_W, radius=SUN_RADIUS_KM, temperature=5772.0, mass=1.989e30)
    other = SimpleNamespace(luminosity=companion_w, radius=SUN_RADIUS_KM, temperature=5772.0, mass=1.989e30)
    rf.attach(host, [rf.Companion(other, 1, axis_au, 0.0, offset)])
    e = math.radians(23.4)
    return SimpleNamespace(rotation_period_hours=24.0, period=1.0, axial_tilt_deg=23.4, distance=1.0, eccentricity=0.0, star=host, body_type="t",
                           spin_axis_x=0.0, spin_axis_y=math.sin(e), spin_axis_z=math.cos(e), orbital_inclination_deg=0.0, orbital_ascending_node_deg=0.0)


def test_a_lone_star_reads_as_before():
    planet = _planet_with(SUN_W, (1.0e4, 0.0, 0.0))
    planet.star = SimpleNamespace(luminosity=SUN_W, radius=SUN_RADIUS_KM)      # no companions attached
    fields = daylight.compute(planet)
    assert fields["day_share_mean"] == pytest.approx(daylight.mean_day_share(daylight.angular_radius(SUN_RADIUS_KM, 1.0)))
    assert fields["day_share_equator"] == pytest.approx(daylight.equator_day_share(23.4, daylight.angular_radius(SUN_RADIUS_KM, 1.0)))


def test_a_companion_too_faint_to_make_a_day_changes_nothing():
    """A companion that puts less than `DAY_SUN_W_M2` on the planet (here a 1e-8 Sun at 1 AU's worth of light) is no sun of the day."""
    faint = _planet_with(SUN_W * 1.0e-8, (-1.0e4, 0.0, 0.0))
    alone = _planet_with(SUN_W * 1.0e-8, (-1.0e4, 0.0, 0.0))
    alone.star = SimpleNamespace(luminosity=SUN_W, radius=SUN_RADIUS_KM)
    assert daylight.compute(faint)["day_share_equator"] == daylight.compute(alone)["day_share_equator"]
    assert daylight.compute(faint)["day_share_mean"] == daylight.compute(alone)["day_share_mean"]


def test_a_bright_wide_companion_lengthens_the_day_and_the_hours_of_daylight():
    planet = _planet_with(SUN_W, (-100.0, 0.0, 0.0), axis_au=100.0)       # 0.14 W/m2 on the planet, over the twilight light
    alone = _planet_with(SUN_W, (-100.0, 0.0, 0.0))
    alone.star = SimpleNamespace(luminosity=SUN_W, radius=SUN_RADIUS_KM)
    both, one = daylight.compute(planet), daylight.compute(alone)
    assert both["day_share_equator"] > one["day_share_equator"] + 0.1
    assert both["day_share_mean"] > one["day_share_mean"] + 0.1
    assert both["daylight_hours"] == pytest.approx(both["solar_day_hours"] * both["day_share_equator"])
    assert both["night_hours"] < one["night_hours"]
    assert both["solar_day_hours"] == one["solar_day_hours"]


def test_the_place_of_a_close_companion_is_taken_from_the_planet_not_the_host():
    """A companion half the planet's distance from the host is seen from the planet nearly where the host is (its direction is the planet's, not the host's)."""
    counted = [(SimpleNamespace(radius=SUN_RADIUS_KM), None), (SimpleNamespace(radius=SUN_RADIUS_KM), (0.0, 0.1, 0.0))]
    far = daylight._sky_sources(0.0, 1.0, counted, 0.0)
    assert far[0][1] == pytest.approx(0.0) and abs(far[1][1]) < math.radians(7.0)
