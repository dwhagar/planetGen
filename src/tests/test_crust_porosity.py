# tests/test_crust_porosity.py

"""
GEN.780: the porosity of a world's upper crust from its gravity and its share of water and ice (physics/crust_porosity).
"""

import math

import pytest

from planetgen.physics import crust_porosity as cp

GRAVITY = {"earth": 9.81, "moon": 1.62, "mars": 3.71, "venus": 8.87}


def test_earth_moon_and_mars_come_out_where_they_are_measured():
    assert 0.02 <= cp.mean_porosity(GRAVITY["earth"]) <= 0.10
    assert 0.08 <= cp.mean_porosity(GRAVITY["moon"]) <= 0.20
    assert 0.10 <= cp.mean_porosity(GRAVITY["mars"]) <= 0.20


def test_a_lower_gravity_keeps_the_pores_open_deeper():
    values = [cp.mean_porosity(g) for g in (20.0, 9.81, 3.71, 1.62, 0.3)]
    assert values == sorted(values)
    assert cp.scale_depth_km(1.62) > cp.scale_depth_km(9.81) > cp.scale_depth_km(30.0)


def test_an_icy_crust_is_more_porous_than_a_rocky_one():
    assert cp.mean_porosity(1.35, 1.0) > cp.mean_porosity(1.35, 0.0)
    assert cp.mean_porosity(1.35, 0.5) == pytest.approx((cp.mean_porosity(1.35, 0.0) + cp.mean_porosity(1.35, 1.0)) / 2.0)
    assert cp.surface_porosity(5.0) == cp.surface_porosity(1.0) and cp.surface_porosity(-1.0) == cp.surface_porosity(0.0)


def test_the_mean_is_the_integral_of_the_exponential():
    depth = cp.scale_depth_km(9.81)
    steps = 20000
    total = sum(cp.surface_porosity(0.0) * math.exp(-(i + 0.5) * cp.REFERENCE_DEPTH_KM / steps / depth) for i in range(steps)) / steps
    assert cp.mean_porosity(9.81) == pytest.approx(total, rel=1e-4)


def test_the_gravity_of_a_body_is_gm_over_r_squared():
    assert cp.gravity_m_s2(5.972e24, 6371.0) == pytest.approx(9.82, rel=0.01)
