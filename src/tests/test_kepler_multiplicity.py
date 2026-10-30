# tests/test_kepler_multiplicity.py

"""
GEN.728: the Kepler transit multiplicity of a simulated sky (`planetgen.physics.kepler_multiplicity`): the dispersion law of the
mutual inclinations, one observation of a system, and the counts over many orientations against the geometry.
"""

import math

import pytest

from planetgen.physics import kepler_multiplicity as km
from planetgen.physics import kepler_recovery as kr
from planetgen.util import draw


def system(periods, radius=2.0, host_radius=1.0, host_mass=1.0, **extra):
    return {"host_mass_solar": host_mass, "host_radius_solar": host_radius,
            "planets": [{"period_d": p, "radius_earth": radius, **extra} for p in periods]}


def test_the_dispersion_falls_steeply_with_the_number_of_planets():
    assert km.inclination_sigma_deg(5) == pytest.approx(0.8)
    assert km.inclination_sigma_deg(10) == pytest.approx(0.8 / 8.0)
    assert km.inclination_sigma_deg(3) > km.inclination_sigma_deg(4) > km.inclination_sigma_deg(5)
    assert km.inclination_sigma_deg(1) == km.MAX_DISPERSION_DEG
    with pytest.raises(ValueError):
        km.inclination_sigma_deg(0)


def test_a_system_with_no_planets_shows_none():
    draw.set_run_seed(5)
    assert km.observe(system([])) == []
    assert km.transit_counts([system([])], orientations=5) == {}


def test_a_planet_inside_its_star_is_never_seen():
    draw.set_run_seed(5)
    assert km.transit_counts([system([0.05])], orientations=50) == {}


def test_the_transit_probability_of_a_single_planet_is_the_stellar_radius_over_the_orbit():
    draw.set_run_seed(11)
    # a 5-day orbit round the Sun: a/R* = 11.1 so about 9 percent of the orientations transit; a 3 Earth-radius planet is easily found
    n = 3000
    item = system([5.0], radius=3.0, inclination_deg=0.0)
    seen = sum(1 for _ in range(n) if km.observe(item))
    expected = 1.0 / kr.scaled_semi_major_axis({"period_d": 5.0, "host_mass_solar": 1.0, "host_radius_solar": 1.0})
    assert seen / n == pytest.approx(expected, rel=0.35)


def test_a_flat_compact_system_shows_several_planets_at_once_and_a_tilted_one_fewer():
    draw.set_run_seed(21)
    flat = system([4.0, 7.0, 12.0, 20.0], radius=4.0, inclination_deg=0.0)
    tilted = system([4.0, 7.0, 12.0, 20.0], radius=4.0, inclination_deg=8.0)
    flat_counts = km.transit_counts([flat], orientations=2500)
    tilted_counts = km.transit_counts([tilted], orientations=2500)
    assert max(flat_counts) >= 3
    assert sum(k * v for k, v in flat_counts.items()) / sum(flat_counts.values()) > \
        sum(k * v for k, v in tilted_counts.items()) / max(sum(tilted_counts.values()), 1)


def test_a_seeded_run_repeats():
    item = system([4.0, 7.0, 12.0], radius=3.0)
    draw.set_run_seed(3)
    first = km.transit_counts([item], orientations=200)
    draw.set_run_seed(3)
    assert km.transit_counts([item], orientations=200) == first


def test_a_found_planet_carries_the_impact_parameter_it_was_found_at():
    draw.set_run_seed(2)
    item = system([6.0], radius=5.0)
    found = []
    for _ in range(400):
        found.extend(km.observe(item))
    assert found and all(0.0 <= p["impact_parameter"] < 1.0 and math.isfinite(p["impact_parameter"]) for p in found)
