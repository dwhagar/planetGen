"""GEN.1492 slice 0: a companion star's secular forcing of a belt's members (`physics/belt_forcing.py`)."""

import math

import pytest

from planetgen.physics import belt_forcing, constants, secular_chain

SUN = constants.SOLAR_MASS_TO_KG
JUPITER = 1.898e27


def test_a_far_companion_forces_what_the_small_alpha_law_says():
    a, axis, e = 3.0, 100.0, 0.4
    forced = belt_forcing.forced_vector(a, SUN, [(0.5 * SUN, axis, e, 0.0)], 0.0)
    assert forced.real == pytest.approx(1.25 * (a / axis) * e / (1.0 - e * e), rel=0.03)
    assert forced.imag == pytest.approx(0.0, abs=1e-12)


def test_the_forced_vector_points_along_the_companions_periapsis():
    forced = belt_forcing.forced_vector(3.0, SUN, [(0.5 * SUN, 100.0, 0.4, 1.0)], 0.0)
    assert math.atan2(forced.imag, forced.real) == pytest.approx(1.0)


def test_a_circular_companion_forces_nothing():
    assert belt_forcing.forced_vector(3.0, SUN, [(0.5 * SUN, 100.0, 0.0, 0.0)], 0.0) == 0j


def test_a_heavy_giant_beside_the_belt_dilutes_the_forcing():
    bare = abs(belt_forcing.forced_vector(3.0, SUN, [(0.5 * SUN, 100.0, 0.4, 0.0)], belt_forcing.planet_rate(3.0, [], [], SUN)))
    rate = belt_forcing.planet_rate(3.0, [5.2], [JUPITER], SUN)
    diluted = abs(belt_forcing.forced_vector(3.0, SUN, [(0.5 * SUN, 100.0, 0.4, 0.0)], rate))
    assert rate > belt_forcing.planet_rate(3.0, [], [], SUN) and diluted < 0.5 * bare


def test_two_companions_add_as_vectors():
    one = belt_forcing.forced_vector(3.0, SUN, [(0.5 * SUN, 100.0, 0.4, 0.0)], 0.0)
    both = belt_forcing.forced_vector(3.0, SUN, [(0.5 * SUN, 100.0, 0.4, 0.0), (0.5 * SUN, 100.0, 0.4, math.pi)], 0.0)
    assert abs(both) < 0.05 * abs(one)       # two equal and opposite pulls: the rates add, the forcing cancels


def test_outside_the_theory_the_member_keeps_its_orbit():
    assert belt_forcing.forced_vector(40.0, SUN, [(0.5 * SUN, 100.0, 0.4, 0.0)], 0.0) == 0j              # alpha 0.4
    assert belt_forcing.forced_vector(3.0, SUN, [], 0.0) == 0j                                           # a lone star
    assert belt_forcing.forced_vector(29.0, SUN, [(1.0 * SUN, 100.0, 0.9, 0.0)], 0.0) == 0j             # a forced eccentricity past the linear limit


def test_the_combined_orbit_is_the_vector_sum():
    e, varpi = belt_forcing.combined(0.1, 0.0, 0.05 + 0j)
    assert (e, varpi) == (pytest.approx(0.15), pytest.approx(0.0))
    e, varpi = belt_forcing.combined(0.1, math.pi, 0.05 + 0j)
    assert e == pytest.approx(0.05) and varpi == pytest.approx(math.pi)
    assert belt_forcing.combined(0.0, 1.0, 0j) == (0.0, 0.0) or belt_forcing.combined(0.0, 1.0, 0j)[0] == 0.0


def test_the_companions_periapsis_is_the_same_for_the_same_system_and_in_a_turn():
    first = belt_forcing.periapsis_rad(SUN, 1, 100.0)
    assert first == belt_forcing.periapsis_rad(SUN, 1, 100.0) and 0.0 <= first < 2.0 * math.pi
    assert first != belt_forcing.periapsis_rad(SUN, 2, 100.0)


def test_the_planet_rate_of_a_lone_star_is_relativity_alone():
    assert belt_forcing.planet_rate(3.0, [], [], SUN) == pytest.approx(secular_chain.relativistic_rates_rad_yr([3.0], SUN, [0.0])[0])
