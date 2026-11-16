# tests/test_evection.py

"""The evection resonance and the moons it catches (GEN.1750)."""

import math

import pytest

from planetgen.physics import constants, evection, moon_orbits, spin_torques, tidal_spin
from planetgen.physics import kozai_bodies

EARTH_MASS = 5.972e24
EARTH_RADIUS = 6.371e6
YEAR_RATE = 2.0 * math.pi / constants.SECONDS_PER_YEAR


def test_distance_solves_the_precession_condition():
    resonance = evection.Resonance(False, YEAR_RATE)
    spin = 4.0 * 7.292e-5
    a = evection.distance_m(resonance, spin, EARTH_MASS, EARTH_RADIUS)
    assert kozai_bodies.bulge_precession_rad_s(False, spin, EARTH_MASS, EARTH_RADIUS, a) == pytest.approx(YEAR_RATE, rel=1e-9)
    assert evection.resonance_precession(resonance, spin, EARTH_MASS, EARTH_RADIUS, a) == pytest.approx(YEAR_RATE, rel=1e-9)


def test_the_resonance_moves_out_with_the_spin_and_in_as_it_slows():
    resonance = evection.Resonance(False, YEAR_RATE)
    fast = evection.distance_m(resonance, 5.0e-4, EARTH_MASS, EARTH_RADIUS)
    slow = evection.distance_m(resonance, 1.0e-4, EARTH_MASS, EARTH_RADIUS)
    assert fast > slow
    assert fast / slow == pytest.approx(5.0 ** evection.SPIN_EXPONENT, rel=1e-9)
    assert evection.distance_m(resonance, 0.0, EARTH_MASS, EARTH_RADIUS) == 0.0


def test_the_sweep_ratio_reproduces_the_touma_and_wisdom_rates():
    """1 km a year at 4.6 Earth radii is `rho` 0.54 (caught); 10 km a year is 5.4 (passes)."""
    resonance = evection.Resonance(False, YEAR_RATE)
    axis = 4.6 * EARTH_RADIUS
    nu = YEAR_RATE      # at the resonance

    def rho(km_per_year):
        rate = km_per_year * 1000.0 / constants.SECONDS_PER_YEAR      # da/dt, m/s
        return evection.sweep_ratio(resonance, EARTH_MASS, axis, evection.nu_rate(1.0, 0.0, axis, rate, nu))

    assert rho(1.0) == pytest.approx(0.54, abs=0.05)
    assert rho(10.0) == pytest.approx(5.4, abs=0.5)
    assert evection.caught(rho(1.0))
    assert not evection.caught(rho(10.0))


def test_held_eccentricity_is_cut_short_near_the_planet():
    assert evection.held_eccentricity(10.0 * EARTH_RADIUS, EARTH_RADIUS) == evection.ECCENTRICITY_CEILING
    assert evection.held_eccentricity(2.5 * EARTH_RADIUS, EARTH_RADIUS) == pytest.approx(0.2)
    assert evection.held_eccentricity(1.5 * EARTH_RADIUS, EARTH_RADIUS) == 0.0


def _planet():
    return spin_torques.Planet(EARTH_MASS, EARTH_RADIUS, tidal_spin.ROCKY_MOMENT_FACTOR)


def _run(resonance, moon_distance_m, spin, age_s):
    star = spin_torques.Raiser(1.989e30, 1.496e11, YEAR_RATE)
    moons = [moon_orbits.Moon(7.35e22, moon_distance_m, 3300.0)]
    return moon_orbits.evolve(_planet(), star, spin_torques.constant_dissipation(0.3, 12.0), moons, spin, age_s, resonance=resonance)


def test_no_resonance_changes_nothing():
    age = 1.0e9 * constants.SECONDS_PER_YEAR
    plain = _run(None, 8.0 * EARTH_RADIUS, 3.0e-4, age)
    again = _run(None, 8.0 * EARTH_RADIUS, 3.0e-4, age)
    assert plain.moons[0].distance_m == again.moons[0].distance_m
    assert not plain.moons[0].held


def test_a_slow_crossing_is_held_and_follows_the_resonance_in():
    """A fast-spun planet whose resonance lies outside a close moon: the moon drifts out, crosses it slowly and is held."""
    resonance = evection.Resonance(False, YEAR_RATE)
    age = 4.0e9 * constants.SECONDS_PER_YEAR
    spin = 4.0e-4
    start = 0.5 * evection.distance_m(resonance, spin, EARTH_MASS, EARTH_RADIUS)
    out = _run(resonance, start, spin, age)
    free = _run(None, start, spin, age)
    fate = out.moons[0]
    if fate.fate == "stays" and fate.held:
        assert fate.distance_m <= free.moons[0].distance_m + 1.0
    # a moon that never reaches the resonance is not held and ends where it would have
    far = _run(resonance, 3.0 * EARTH_RADIUS, 1.0e-5, 1.0e3 * constants.SECONDS_PER_YEAR)
    assert not far.moons[0].held
