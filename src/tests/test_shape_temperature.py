# tests/test_shape_temperature.py

"""
The temperature a body took its shape at (GEN.1201): `physics/shape_temperature.py`.
"""

import math

import pytest

from planetgen.physics import constants, radiogenic, roundness, shape_temperature

CONTENT = radiogenic.Content(k0_ppm=260.0, th0_ppm=0.0834, u0_ppm=0.0218, al26_ratio=5.0e-5)
NO_AL = radiogenic.Content(k0_ppm=260.0, th0_ppm=0.0834, u0_ppm=0.0218, al26_ratio=0.0)


def _uniform_binding(mass, radius):
    return 0.6 * constants.G * mass ** 2 / radius


def _rock(radius_m, content=CONTENT, surface=150.0, density=3000.0, tidal_w=0.0):
    mass = 4.0 / 3.0 * math.pi * radius_m ** 3 * density
    return shape_temperature.shape_temperature_k(mass, radius_m, 0.0, 0.0, 1.0, surface, _uniform_binding(mass, radius_m), content, tidal_w=tidal_w)


def test_a_small_cold_body_keeps_its_surface_temperature():
    assert _rock(2.0e3, NO_AL) == pytest.approx(150.0, abs=0.5)


def test_26al_warms_a_body_the_more_the_larger_it_is_until_the_ceiling():
    radii = [1.0e3, 3.0e3, 1.0e4, 1.0e5, 1.0e6]
    temps = [_rock(r) for r in radii]
    assert temps == sorted(temps) and temps[0] < 200.0 and temps[-1] > 1500.0       # a kilometre-sized body keeps its cold; one of 10 km radius has melted


def test_no_aluminium_leaves_a_planetesimal_cold_and_a_moon_hot_from_its_assembly():
    assert _rock(1.0e5, NO_AL) < 200.0
    moon = _rock(1.737e6, NO_AL)
    assert 1000.0 < moon <= shape_temperature.ceiling_k(7.3e22, 1.737e6)


def test_the_ceiling_holds_a_planet_at_its_liquidus():
    mass = 5.972e24
    t = shape_temperature.shape_temperature_k(mass, 6.371e6, 0.0, 0.3, 0.7, 250.0, _uniform_binding(mass, 6.371e6), CONTENT)
    assert t == pytest.approx(shape_temperature.ceiling_k(mass, 6.371e6))


def test_the_tides_warm_a_small_moon():
    quiet = _rock(1.0e5, NO_AL)
    tidal = _rock(1.0e5, NO_AL, tidal_w=1.0e12)
    assert tidal > quiet + 10.0


def test_ice_has_a_bigger_heat_capacity_so_it_warms_less_than_rock():
    mass = 4.0 / 3.0 * math.pi * 2.0e5 ** 3 * 1500.0
    icy = shape_temperature.shape_temperature_k(mass, 2.0e5, 0.6, 0.0, 0.4, 100.0, _uniform_binding(mass, 2.0e5), None)
    rocky = shape_temperature.shape_temperature_k(mass, 2.0e5, 0.0, 0.0, 1.0, 100.0, _uniform_binding(mass, 2.0e5), None)
    assert 100.0 <= icy < rocky


def test_a_body_with_no_stored_content_makes_no_radiogenic_heat():
    assert shape_temperature.radiogenic_energy_j_kg(None, 1.0e13, 1.0) == 0.0
    assert shape_temperature.radiogenic_energy_j_kg(CONTENT, 1.0e13, 1.0) > 0.0


def test_the_conduction_time_grows_as_the_radius_squared():
    assert shape_temperature.conduction_time_s(2.0e4) == pytest.approx(4.0 * shape_temperature.conduction_time_s(1.0e4))
    assert shape_temperature.conduction_time_s(5.0e4) / 3.156e7 == pytest.approx(8.0e6, rel=0.1)       # about 8 Myr for a body 100 km across


def test_the_shape_temperature_feeds_the_strength_of_an_icy_body():
    cold = roundness.yield_strength_pa(0.5, 0.0, 0.0, None)
    warm = roundness.yield_strength_pa(0.5, 0.0, 0.0, 250.0)
    assert warm < cold


def test_bad_arguments_are_refused():
    with pytest.raises(shape_temperature.ShapeTemperatureError):
        shape_temperature.shape_temperature_k(0.0, 1.0e5, 0.0, 0.0, 1.0, 150.0)
    with pytest.raises(shape_temperature.ShapeTemperatureError):
        shape_temperature.shape_temperature_k(1.0e18, 1.0e5, 0.6, 0.6, 0.6, 150.0)
    with pytest.raises(shape_temperature.ShapeTemperatureError):
        shape_temperature.shape_temperature_k(1.0e18, 1.0e5, 0.0, 0.0, 1.0, 150.0, tidal_w=-1.0)
