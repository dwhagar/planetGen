# tests/test_daylight.py

"""
GEN.767: the day and night time shares and the hours of daylight (`physics/daylight.py`), checked against hand-worked cases and the
bodies the generator makes.
"""

import math
from types import SimpleNamespace

import pytest

from planetgen.generation.system import StarSystem, SystemConfig
from planetgen.physics import constants, daylight, spin
from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy

SUN_RADIUS_KM = 695_700.0
EARTH_ALPHA = daylight.angular_radius(SUN_RADIUS_KM, 1.0)
EARTH_YEAR_HOURS = 8766.0


def _body(rotation_hours, period_years, tilt_deg, distance_au=1.0, **extra):
    star = SimpleNamespace(radius=SUN_RADIUS_KM)
    return SimpleNamespace(rotation_period_hours=rotation_hours, period=period_years, axial_tilt_deg=tilt_deg, distance=distance_au,
                           star=star, **extra)


def test_the_sun_from_earth_is_a_quarter_of_a_degree_in_radius():
    assert math.degrees(EARTH_ALPHA) == pytest.approx(0.2665, abs=0.001)
    assert daylight.angular_radius(None, 1.0) == 0.0 and daylight.angular_radius(SUN_RADIUS_KM, None) == 0.0


def test_earth_has_about_half_a_day_at_the_equator_and_twelve_hours_of_daylight():
    share = daylight.equator_day_share(23.44, EARTH_ALPHA)
    assert share == pytest.approx(0.5015, abs=0.0005)        # the upper limb adds a few minutes to the half
    assert daylight.equator_day_share(0.0, 0.0) == pytest.approx(0.5)
    earth = daylight.compute(_body(23.934, 1.0, 23.44))
    assert earth["solar_day_hours"] == pytest.approx(24.0, abs=0.01)
    assert earth["daylight_hours"] == pytest.approx(12.04, abs=0.02)
    assert earth["day_share_mean"] == pytest.approx(0.5 * (1.0 + math.sin(EARTH_ALPHA)))
    assert earth["parent_shadow_share"] is None


def test_the_polar_circle_follows_the_tilt():
    assert daylight.polar_circle_deg(23.44, EARTH_ALPHA) == pytest.approx(66.29, abs=0.01)
    assert daylight.polar_circle_deg(97.77, EARTH_ALPHA) == pytest.approx(7.77 - math.degrees(EARTH_ALPHA), abs=0.01)  # Uranus
    assert daylight.polar_circle_deg(180.0 - 23.44, EARTH_ALPHA) == daylight.polar_circle_deg(23.44, EARTH_ALPHA)
    assert daylight.polar_circle_deg(0.0, EARTH_ALPHA) > 89.0
    assert daylight.polar_circle_deg(90.0, EARTH_ALPHA) == 0.0


def test_a_day_at_a_latitude_gets_longer_in_summer_and_a_pole_has_a_polar_night():
    assert daylight.day_share_at(45.0, 23.44, EARTH_ALPHA) > 0.5 > daylight.day_share_at(45.0, -23.44, EARTH_ALPHA)
    assert daylight.day_share_at(80.0, 23.44, EARTH_ALPHA) == 1.0
    assert daylight.day_share_at(80.0, -23.44, EARTH_ALPHA) == 0.0


def test_a_fast_spinner_and_a_retrograde_one_get_their_solar_day_from_both_motions():
    fast = daylight.compute(_body(10.0, 1.0, 5.0))
    assert fast["solar_day_hours"] == pytest.approx(1.0 / (1.0 / 10.0 - 1.0 / EARTH_YEAR_HOURS), rel=1e-3)
    assert fast["daylight_hours"] == pytest.approx(fast["solar_day_hours"] * fast["day_share_equator"])
    assert fast["daylight_hours"] == pytest.approx(5.0, abs=0.1)
    venus = daylight.compute(_body(5832.5, 0.615, 177.4, distance_au=0.723))   # 243 d retrograde against a 224.7 d year
    assert venus["solar_day_hours"] == pytest.approx(116.75 * 24.0, rel=0.01)    # Venus's solar day is 116.75 Earth days


def test_io_and_the_moon_spend_the_shares_of_their_orbits_the_hand_work_gives_in_their_planets_shadow():
    io = daylight.parent_shadow_share(421_800.0, 71_492.0, 3.1)
    assert daylight.parent_shadow_share(421_800.0, 71_492.0, 0.0) == pytest.approx(
        math.asin(71_492.0 / 421_800.0) / math.pi, rel=1e-6)           # in the planet's plane: the disc's angular radius over pi
    assert io == pytest.approx(0.0528, abs=0.0005)                      # 2.2 hours of each 42.5 hour orbit
    moon = daylight.parent_shadow_share(384_400.0, 6_371.0, 5.145)
    assert moon == pytest.approx(4.9e-4, abs=0.3e-4)                    # eclipse seasons of about 43 days a year, hours each month
    assert moon < 0.1 * daylight.parent_shadow_share(384_400.0, 6_371.0, 0.0)


def test_the_planet_blocks_the_star_only_while_the_star_is_within_its_angular_radius_of_the_moons_plane():
    theta = math.degrees(math.asin(6_371.0 / 384_400.0))
    assert daylight.shadow_share_at(384_400.0, 6_371.0, 0.9 * theta) > 0.0
    assert daylight.shadow_share_at(384_400.0, 6_371.0, 1.1 * theta) == 0.0
    assert daylight.shadow_share_at(384_400.0, 6_371.0, 60.0) == 0.0
    assert daylight.shadow_share_at(5_000.0, 6_371.0, 0.0) == 0.0     # inside the planet: no orbit
    assert daylight.parent_shadow_share(384_400.0, 6_371.0, 60.0) < daylight.parent_shadow_share(384_400.0, 6_371.0, 5.0)


def test_a_moon_loses_the_shadow_share_of_its_day():
    star = SimpleNamespace(radius=SUN_RADIUS_KM)
    planet = SimpleNamespace(radius=71_492.0, distance=5.2, period=11.86, star=star, orbital_inclination_deg=1.3,
                             orbital_ascending_node_deg=100.0)
    normal = spin.orbit_normal(1.3, 100.0)
    moon = SimpleNamespace(rotation_period_hours=42.46, period=42.46 / EARTH_YEAR_HOURS, axial_tilt_deg=0.0,
                           distance=421_800.0 / constants.AU_TO_KM, star=star, orbital_inclination_deg=1.3 + 0.2,
                           orbital_ascending_node_deg=100.0, spin_axis_x=normal[0], spin_axis_y=normal[1], spin_axis_z=normal[2])
    alone = daylight.compute(moon, SimpleNamespace(**{**vars(planet), "radius": 1.0}))
    shadowed = daylight.compute(moon, planet)
    assert alone["parent_shadow_share"] < 1e-6
    assert shadowed["parent_shadow_share"] == pytest.approx(0.0542, abs=0.002)
    assert shadowed["day_share_mean"] == pytest.approx(alone["day_share_mean"] * (1.0 - shadowed["parent_shadow_share"]), rel=1e-6)
    assert shadowed["day_share_equator"] == pytest.approx(alone["day_share_equator"] * (1.0 - shadowed["parent_shadow_share"]), rel=1e-6)
    # a moon locked to its planet turns once in 42.46 h against the 11.86 y year of the planet: its noon is its orbit, a little over
    assert shadowed["solar_day_hours"] == pytest.approx(42.46, rel=0.01)


def test_a_body_with_no_rotation_has_nothing():
    assert all(value is None for value in daylight.compute(SimpleNamespace(rotation_period_hours=None, star=object())).values())


def test_every_generated_planet_and_moon_stores_shares_that_make_sense():
    seen = moons = 0
    for seed in range(6):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            for body in [planet, *getattr(planet, "moons", [])]:
                if not getattr(body, "rotation_period_hours", None):
                    continue
                seen += 1
                assert 0.0 < body.day_share_mean <= 1.0 and 0.0 <= body.polar_circle_deg <= 90.0
                if body.day_share_equator is None:
                    assert body.solar_day_hours is None and body.daylight_hours is None
                else:
                    assert 0.0 <= body.day_share_equator <= 1.0 and body.solar_day_hours > 0.0
                    assert body.daylight_hours == pytest.approx(body.solar_day_hours * body.day_share_equator)
                if body.is_moon:
                    moons += 1
                    assert 0.0 <= body.parent_shadow_share < 0.5
                else:
                    assert body.parent_shadow_share is None
    assert seen > 6 and moons > 0


def test_the_longest_eclipse_is_the_planets_angular_radius_over_pi_of_an_orbit():
    """UX.205: the pass through the middle of the shadow: Io (421,700 km from Jupiter, 71,492 km radius, 42.46 h orbit) is
    eclipsed for about 2.3 hours at the longest; a moon not outside its planet, or with no period, has none."""
    hours = daylight.longest_eclipse_hours(421700.0, 71492.0, 42.46)
    assert hours == pytest.approx(math.asin(71492.0 / 421700.0) / math.pi * 42.46, rel=1e-9)
    assert 2.2 < hours < 2.4
    assert daylight.longest_eclipse_hours(60000.0, 71492.0, 42.46) == 0.0
    assert daylight.longest_eclipse_hours(421700.0, 71492.0, 0.0) == 0.0


def test_a_world_that_turns_once_an_orbit_has_no_day_and_a_little_over_half_of_it_is_lit():
    locked = daylight.compute(_body(0.25 * EARTH_YEAR_HOURS, 0.25, 0.0, distance_au=0.2))
    assert locked["day_share_equator"] is None and locked["solar_day_hours"] is None and locked["daylight_hours"] is None
    alpha = daylight.angular_radius(SUN_RADIUS_KM, 0.2)
    assert locked["day_share_mean"] == pytest.approx(0.5 * (1.0 + math.sin(alpha)))
