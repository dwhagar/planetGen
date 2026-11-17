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


# --- the tilt (GEN.1492, second part) -----------------------------------------------------------------------------------------------------------

def _normal(inclination_deg, node_deg=0.0):
    from planetgen.physics import spin
    return spin.orbit_normal(inclination_deg, node_deg)


def _inclination_of(axis):
    return math.degrees(math.acos(axis[2]))


def test_a_bare_belt_precesses_about_the_companions_own_plane():
    axis = belt_forcing.forced_normal(3.0, SUN, [(0.5 * SUN, 100.0, 0.2, _normal(20.0, 40.0))], 0.0)
    assert axis == pytest.approx(_normal(20.0, 40.0), abs=1e-9)


def test_a_coplanar_companion_gives_no_tilt():
    axis = belt_forcing.forced_normal(3.0, SUN, [(0.5 * SUN, 100.0, 0.2, (0.0, 0.0, 1.0))], 0.0)
    assert axis == pytest.approx((0.0, 0.0, 1.0), abs=1e-12)


def test_a_heavy_giant_holds_the_belt_flat_and_a_light_one_less():
    companion = [(0.5 * SUN, 100.0, 0.2, _normal(20.0))]
    bare = _inclination_of(belt_forcing.forced_normal(3.0, SUN, companion, 0.0))
    giant = belt_forcing.planet_nodal_rate(3.0, [5.2], [JUPITER], SUN)
    held = _inclination_of(belt_forcing.forced_normal(3.0, SUN, companion, giant))
    heavy = _inclination_of(belt_forcing.forced_normal(3.0, SUN, companion, 100.0 * giant))
    assert heavy < held < bare
    assert heavy < 0.5


def test_equal_pulls_put_the_axis_half_way():
    rate, _ = belt_forcing.companion_rates(3.0, SUN, 0.5 * SUN, 100.0, 0.2)
    axis = belt_forcing.forced_normal(3.0, SUN, [(0.5 * SUN, 100.0, 0.2, _normal(20.0))], rate)
    assert _inclination_of(axis) == pytest.approx(10.0, abs=0.1)


def test_a_companion_out_of_the_theory_or_with_no_plane_changes_nothing():
    assert belt_forcing.forced_normal(3.0, SUN, [(0.5 * SUN, 5.0, 0.2, _normal(20.0))], 0.0) is None       # alpha beyond ALPHA_LIMIT
    assert belt_forcing.forced_normal(3.0, SUN, [(0.5 * SUN, 100.0, 0.2, None)], 0.0) is None
    assert belt_forcing.forced_normal(3.0, SUN, [], 1.0) is None
    assert belt_forcing.forced_normal(3.0, SUN, [(0.5 * SUN, 100.0, 0.2, _normal(60.0))], 0.0) is None     # a tilt in the Kozai-Lidov range keeps the drawn orbit


def test_a_lone_planet_gives_no_nodal_rate_without_planets_and_a_positive_one_with():
    assert belt_forcing.planet_nodal_rate(3.0, [], [], SUN) == 0.0
    assert belt_forcing.planet_nodal_rate(3.0, [5.2], [JUPITER], SUN) > 0.0


@pytest.mark.parametrize("phase", [0.0, 1.0, 2.5, 4.0])
def test_the_angle_to_the_forced_axis_is_kept_whatever_the_phase(phase):
    axis = _normal(20.0, 40.0)
    drawn = _normal(3.0, 200.0)
    inclination, node = belt_forcing.tilted_orbit(3.0, 200.0, axis, phase)
    turned = _normal(inclination, node)
    assert sum(a * b for a, b in zip(turned, axis)) == pytest.approx(sum(a * b for a, b in zip(drawn, axis)), abs=1e-9)


def test_a_flat_member_swings_between_flat_and_twice_the_forced_inclination():
    axis = _normal(20.0, 40.0)
    seen = [belt_forcing.tilted_orbit(0.0, 10.0, axis, 2.0 * math.pi * k / 72.0)[0] for k in range(72)]
    assert min(seen) == pytest.approx(0.0, abs=0.5)
    assert max(seen) == pytest.approx(40.0, abs=0.5)
    assert 20.0 < sum(seen) / len(seen) < 30.0         # the angle to the plane is not linear in the phase: about 25 degrees on average


def test_an_axis_on_the_planets_plane_keeps_the_inclination_and_turns_only_the_node():
    inclination, node = belt_forcing.tilted_orbit(5.0, 100.0, (0.0, 0.0, 1.0), 1.0)
    assert inclination == pytest.approx(5.0, abs=1e-9)
    assert node == pytest.approx((100.0 + math.degrees(1.0)) % 360.0, abs=1e-6)


def test_the_inclination_spreads_between_the_difference_and_the_sum_of_the_two_angles():
    axis = _normal(20.0, 0.0)
    drawn = _normal(5.0, 77.0)
    angle = math.degrees(math.acos(sum(a * b for a, b in zip(drawn, axis))))
    seen = [belt_forcing.tilted_orbit(5.0, 77.0, axis, 2.0 * math.pi * k / 72.0)[0] for k in range(72)]
    assert min(seen) >= abs(20.0 - angle) - 1e-6 and max(seen) <= 20.0 + angle + 1e-6
    assert min(seen) < abs(20.0 - angle) + 0.5 and max(seen) > 20.0 + angle - 0.5
