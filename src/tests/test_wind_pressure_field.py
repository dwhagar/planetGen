# tests/test_wind_pressure_field.py

"""GEN.1452: the wind of every star on a body, and a moon inside its parent's magnetosphere (`planetgen.physics.magnetism`)."""

import math
from types import SimpleNamespace

import pytest

from planetgen.physics import constants, magnetism, radiation_belts, radiation_field as rf


def sun(**extra):
    base = dict(luminosity=constants.SOLAR_LUMINOSITY, temperature=5772.0, mass=1.989e30, radius=6.957e5, log_lx_lbol=-6.5)
    base.update(extra)
    return SimpleNamespace(**base)


def with_companion(host, companion, axis_au, eccentricity=0.0):
    rf.attach(host, [rf.Companion(companion, 1, axis_au, eccentricity)])
    return host


def test_a_lone_star_gives_the_hosts_wind_to_the_last_digit():
    host = sun()
    for distance in (0.1, 1.0, 5.2, 40.0):
        assert magnetism.wind_pressure_at(host, distance) == magnetism.wind_pressure_pa(host, distance)


def test_a_companion_adds_its_wind_by_the_mean_inverse_square_of_its_orbit():
    host = with_companion(sun(), sun(), 5.0)
    expected = magnetism.wind_pressure_pa(host, 1.5) + magnetism.wind_pressure_pa(sun(), 1.0) * rf.mean_inverse_square(1.5, 5.0, 0.0)
    assert magnetism.wind_pressure_at(host, 1.5) == pytest.approx(expected, rel=1e-12)
    share = magnetism.wind_pressure_at(host, 1.5) / magnetism.wind_pressure_pa(host, 1.5) - 1.0
    assert 0.09 <= share <= 0.11                                                    # the same tenth the light of an equal star at 5 AU adds to a planet at 1.5


def test_a_companion_with_no_wind_data_gives_none_and_a_pulsar_companion_gives_its_pulsar_wind(monkeypatch):
    bare = SimpleNamespace(luminosity=None, radius=None, log_lx_lbol=None, mass=1e30)
    host = with_companion(sun(), bare, 5.0)
    assert magnetism.wind_pressure_at(host, 1.0) == magnetism.wind_pressure_pa(host, 1.0)       # not the Earth's wind, which the host alone falls back to
    monkeypatch.setattr(magnetism.activity, "neutron_star_spin_down_w", lambda star: 1.0e26)
    pulsar = SimpleNamespace(pulsar_type="young", luminosity=None, radius=None, log_lx_lbol=None, mass=2.8e30)
    host2 = with_companion(sun(), pulsar, 5.0)
    expected = magnetism.wind_pressure_pa(host2, 1.0) + magnetism.pulsar_wind_pressure_pa(1.0e26, 1.0) * rf.mean_inverse_square(1.0, 5.0, 0.0)
    assert magnetism.wind_pressure_at(host2, 1.0) == pytest.approx(expected, rel=1e-12) and expected > magnetism.wind_pressure_pa(host2, 1.0)


def test_the_magnetopause_of_a_planet_shrinks_with_a_companion_and_is_unchanged_for_a_lone_star():
    planet = SimpleNamespace(magnetic_moment_a_m2=magnetism.EARTH_MOMENT_A_M2, radius=constants.EARTH_RADIUS_KM, star=sun())
    lone = magnetism.magnetopause_rp(planet, 1.0)
    assert lone == pytest.approx(magnetism.EARTH_STANDOFF_RP * (magnetism.EARTH_WIND_PA / magnetism.wind_pressure_pa(planet.star, 1.0)) ** (1 / 6), rel=1e-9)
    planet.star = with_companion(sun(), sun(), 3.0)
    assert magnetism.magnetopause_rp(planet, 1.0) < lone


def test_the_belts_of_a_planet_read_the_winds_of_every_star():
    host = with_companion(sun(), sun(), 3.0)
    lone = SimpleNamespace(magnetic_moment_a_m2=magnetism.EARTH_MOMENT_A_M2, radius=constants.EARTH_RADIUS_KM, mass=5.97e24, rotation_period_hours=24.0, star=sun(),
                           _star_distance_au=1.0)
    pair = SimpleNamespace(**{**vars(lone), "star": host})
    outer = lambda body: radiation_belts.belts_of_body(body)[1].outer_l
    assert outer(pair) < outer(lone)


def test_the_pressure_of_jupiters_field_at_ganymedes_orbit_gives_ganymede_a_standoff_near_two_radii():
    jupiter_moment, ganymede_orbit_m, ganymede_moment = 1.55e27, 1.0704e9, 1.3e20
    pressure = magnetism.magnetosphere_pressure_pa(jupiter_moment, ganymede_orbit_m)
    assert 2e-9 < pressure < 1e-8                                                  # B of about 120 nT, a few nPa
    standoff = magnetism.EARTH_STANDOFF_RP * (ganymede_moment / magnetism.EARTH_MOMENT_A_M2) ** (1 / 3) * (pressure / magnetism.EARTH_WIND_PA) ** (-1 / 6)
    assert 1.5 < standoff * constants.EARTH_RADIUS_KM / 2634.0 < 3.0               # observed: about 2 Ganymede radii


def test_a_moon_inside_its_parents_magnetopause_holds_off_the_parent_field_and_one_outside_keeps_the_wind():
    jupiter = SimpleNamespace(magnetic_moment_a_m2=1.55e27, magnetopause_rp=60.0, radius=71492.0)
    inside = SimpleNamespace(magnetic_moment_a_m2=1.3e20, distance=1.0704e9 / constants.AU_M, radius=2634.0, magnetopause_rp=999.0)
    outside = SimpleNamespace(magnetic_moment_a_m2=1.3e20, distance=61.0 * 71492.0e3 / constants.AU_M, radius=2634.0, magnetopause_rp=999.0)
    magnetism.apply_parent_magnetosphere(inside, jupiter)
    magnetism.apply_parent_magnetosphere(outside, jupiter)
    assert 1.5 < inside.magnetopause_rp < 3.0 and outside.magnetopause_rp == 999.0
    fieldless, parentless = SimpleNamespace(**vars(inside)), SimpleNamespace(**vars(inside))
    fieldless.magnetic_moment_a_m2, fieldless.magnetopause_rp = 0.0, 7.0
    parentless.magnetopause_rp = 7.0
    magnetism.apply_parent_magnetosphere(fieldless, jupiter)
    magnetism.apply_parent_magnetosphere(parentless, SimpleNamespace(magnetic_moment_a_m2=0.0, magnetopause_rp=None, radius=1.0))
    assert fieldless.magnetopause_rp == 7.0 and parentless.magnetopause_rp == 7.0
