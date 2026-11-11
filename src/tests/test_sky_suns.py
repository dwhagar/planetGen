# tests/test_sky_suns.py

"""GEN.1458: the other suns of a planet's sky (`planetgen.physics.sky_suns`) and what the light laws do with them."""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from planetgen.physics import constants, latitude_bands, radiation_field as rf, regional_temperature as rt, sky_suns, surface_cycle

SUN_W = constants.SOLAR_LUMINOSITY


def host_with(offset, luminosity_w=SUN_W, axis_au=100.0):
    host = SimpleNamespace(luminosity=SUN_W, temperature=5772.0, mass=1.989e30)
    other = SimpleNamespace(luminosity=luminosity_w, temperature=5772.0, mass=1.989e30)
    rf.attach(host, [rf.Companion(other, 1, axis_au, 0.0, offset)])
    return host


def tilted(tilt_deg):
    e = math.radians(tilt_deg)
    return (0.0, math.sin(e), math.cos(e)), (0.0, 0.0, 1.0)      # spin axis, orbit normal


def test_a_lone_star_has_no_other_suns():
    host = SimpleNamespace(luminosity=SUN_W)
    axis, normal = tilted(23.4)
    assert sky_suns.suns_of(host, axis, normal, 1.0) == ()
    assert sky_suns.encode(()) is None and sky_suns.decode(None) == ()


def test_the_direction_is_written_in_the_equinox_frame():
    axis, normal = tilted(30.0)
    x, y, n, tilt = sky_suns.frame(axis, normal)
    assert math.degrees(tilt) == pytest.approx(30.0)
    assert np.allclose(np.dot(x, axis), 0.0, atol=1e-12)           # the equinox line lies in the equator
    assert np.allclose(np.dot(x, n), 0.0, atol=1e-12)              # and in the orbit plane
    host = host_with((0.0, 10.0, 0.0), axis_au=10.0)                  # along the planet's +y
    (sun,) = sky_suns.suns_of(host, axis, normal, 1.0)
    assert sun.ux ** 2 + sun.uy ** 2 + sun.uz ** 2 == pytest.approx(1.0)
    assert sun.uz == pytest.approx(0.0, abs=1e-12)                 # the companion lies in the orbit plane


def test_the_share_is_the_companions_part_of_the_mean_light():
    axis, normal = tilted(20.0)
    host = host_with((10.0, 0.0, 0.0), axis_au=10.0)
    (sun,) = sky_suns.suns_of(host, axis, normal, 1.0)
    assert sun.share == pytest.approx(rf.companion_share(host, 1.0) / (1.0 + rf.companion_share(host, 1.0)), rel=1e-6)


def test_a_distant_companion_under_the_minimum_share_is_left_out():
    axis, normal = tilted(20.0)
    host = host_with((1.0e5, 0.0, 0.0), axis_au=1.0e5)
    assert sky_suns.suns_of(host, axis, normal, 1.0) == ()


def test_the_sphere_mean_of_the_year_is_still_a_quarter_with_other_suns():
    latitudes = np.radians(np.linspace(-89.5, 89.5, 359))
    weight = np.cos(latitudes)
    suns = (sky_suns.Sun(0.3, 0.6, 0.64, 0.48 ** 0.5),)
    tilt = math.radians(40.0)
    lone = float(4.0 * np.sum(weight * rt._annual_insolation(latitudes, tilt)) / np.sum(weight))
    both = float(4.0 * np.sum(weight * rt._annual_insolation(latitudes, tilt, suns)) / np.sum(weight))
    assert lone == pytest.approx(1.0, abs=2e-3) and both == pytest.approx(1.0, abs=2e-3)


def test_with_no_tilt_another_sun_changes_nothing():
    latitudes = np.radians(np.linspace(-85.0, 85.0, 35))
    suns = (sky_suns.Sun(0.4, 1.0, 0.0, 0.0),)
    assert np.allclose(rt._annual_insolation(latitudes, 0.0, suns), rt._annual_insolation(latitudes, 0.0))


def test_a_sun_over_the_pole_lights_the_polar_winter():
    latitudes = np.radians(np.array([80.0]))
    e = math.radians(60.0)
    over_the_pole = sky_suns.Sun(0.2, 0.0, math.sin(e), math.cos(e))     # along the spin axis
    plain = rt._annual_insolation(latitudes, e)
    lit = rt._annual_insolation(latitudes, e, (over_the_pole,))
    assert lit[0] > 0.8 * plain[0], "a pole sun at 20 percent of the light keeps the pole warm all year"
    assert rt._annual_insolation(np.radians(np.array([-80.0])), e, (over_the_pole,))[0] < rt._annual_insolation(np.radians(np.array([-80.0])), e)[0]


def test_the_cell_series_keeps_the_mean_light_with_other_suns():
    suns = (sky_suns.Sun(0.35, 0.8, 0.5, math.sqrt(1 - 0.64 - 0.25)),)
    for tilt in (0.0, 25.0, 70.0):
        lone = surface_cycle.series(24.0, 8766.0, tilt, rows=12, columns=12)
        both = surface_cycle.series(24.0, 8766.0, tilt, rows=12, columns=12, suns=suns)
        weight = np.repeat(lone.grid.row_weight / 12, 12).reshape(12, 12)
        mean = lambda s: float(np.sum(np.maximum(s.cosine, 0.0).mean(axis=2) * weight))
        assert mean(lone) == pytest.approx(0.25, rel=0.02)
        assert mean(both) == pytest.approx(0.25, rel=0.02), tilt


def test_the_cell_series_without_other_suns_is_what_it_was():
    a = surface_cycle.series(30.0, 5000.0, 20.0, rows=6, columns=6)
    b = surface_cycle.series(30.0, 5000.0, 20.0, rows=6, columns=6, suns=())
    assert np.array_equal(a.cosine, b.cosine)


def test_a_second_sun_lights_a_cell_the_host_leaves_in_the_dark():
    """A locked-like slow spin: where the host is down, a companion on the other side of the sky is up."""
    suns = (sky_suns.Sun(0.5, -1.0, 0.0, 0.0),)
    both = surface_cycle.series(24.0, 8766.0, 0.0, rows=6, columns=12, suns=suns)
    assert (np.maximum(both.cosine, 0.0) > 0.0).mean() > (np.maximum(surface_cycle.series(24.0, 8766.0, 0.0, rows=6, columns=12).cosine, 0.0) > 0.0).mean()


def test_the_stored_text_round_trips():
    suns = (sky_suns.Sun(0.123456, 0.1, -0.2, 0.9746794), sky_suns.Sun(0.01, 1.0, 0.0, 0.0))
    back = sky_suns.decode(sky_suns.encode(suns))
    assert [(round(s.share, 6), round(s.ux, 6)) for s in back] == [(0.123456, 0.1), (0.01, 1.0)]


def test_the_bands_read_the_stored_suns():
    values = {"surface_temperature_k": 288.0, "atmospheric_pressure_pa": 1.0e5, "gravity_g": 1.0, "rotation_period_hours": 24.0, "axial_tilt_deg": 50.0, "period_years": 1.0}
    lone = latitude_bands.compute(values, 0.27)
    e = math.radians(50.0)
    suns = sky_suns.encode((sky_suns.Sun(0.25, 0.0, math.sin(e), math.cos(e)),))
    both = latitude_bands.compute({**values, "sky_suns": suns}, 0.27)
    assert lone and both
    assert both[0].mean_k != lone[0].mean_k, "a sun over the north pole warms the north band"
