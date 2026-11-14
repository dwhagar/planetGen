"""
GEN.1454: how far the tilt of a planet wanders under the pull of the other bodies of its system (`physics/obliquity_variation.py`). The Solar System is the test point, never an input:
the eight planets' masses and orbits (J2000 values to the invariable plane, recalled), the Earth with its Moon and without, and Mars.
"""

import math

import pytest

from planetgen.physics import constants, obliquity_variation as ov, secular_chain

AU = constants.AU_TO_M
SUN = 1.98847e30
A_AU = [0.3871, 0.7233, 1.0, 1.5237, 5.2026, 9.5549, 19.2184, 30.1104]
MASS = [3.30e23, 4.87e24, 5.97e24, 6.42e23, 1.898e27, 5.683e26, 8.68e25, 1.024e26]
INCLINATION = [math.radians(x) for x in (6.34, 2.19, 1.57, 1.67, 0.32, 0.93, 1.02, 0.72)]
NODE = [math.radians(x) for x in (32.2, 52.3, 284.5, 352.9, 312.9, 118.8, 307.8, 192.0)]
EARTH, MARS = 2, 3
MOON = ov.Torquer(7.342e22, 3.844e8, 0.055)


def _earth(with_moon):
    ellipticity = ov.dynamical_ellipticity(5.972e24, 6.371e6, 7.292e-5, 0.3307)
    torquers = [ov.Torquer(SUN, 1.0 * AU, 0.0167)] + ([MOON] if with_moon else [])
    alpha = ov.precession_constant_rad_yr(ellipticity, 7.292e-5, torquers)
    return alpha, ov.all_sibling_modes(A_AU, MASS, SUN, INCLINATION, NODE)[EARTH]


def test_the_fluid_love_number_and_ellipticity_of_the_earth_are_the_measured_ones_to_a_few_percent():
    assert ov.fluid_love_number(0.3307) == pytest.approx(0.917, abs=0.005)
    assert ov.dynamical_ellipticity(5.972e24, 6.371e6, 7.292e-5, 0.3307) == pytest.approx(0.003274, rel=0.04)


def test_the_precession_constant_of_the_earth_is_the_sun_and_the_moon_together():
    ellipticity = ov.dynamical_ellipticity(5.972e24, 6.371e6, 7.292e-5, 0.3307)
    sun = ov.precession_constant_rad_yr(ellipticity, 7.292e-5, [ov.Torquer(SUN, AU, 0.0167)]) * ov.ARCSEC_PER_RAD
    both = ov.precession_constant_rad_yr(ellipticity, 7.292e-5, [ov.Torquer(SUN, AU, 0.0167), MOON]) * ov.ARCSEC_PER_RAD
    assert sun == pytest.approx(16.0, rel=0.1)             # the solar share of the 50.3 arcseconds a year
    assert both == pytest.approx(50.3, rel=0.1)
    assert both > 2.5 * sun                                # the Moon is what makes the constant large


def test_the_nodal_modes_of_the_solar_system_are_the_laplace_lagrange_ones():
    _a, B = secular_chain.matrices(A_AU, MASS, SUN)
    frequencies = sorted(f * ov.ARCSEC_PER_RAD for f in secular_chain._modes(B)[0])
    assert frequencies[0] == pytest.approx(-25.7, abs=1.0) and frequencies[1] == pytest.approx(-18.7, abs=1.0)
    assert frequencies[-2] == pytest.approx(-0.7, abs=0.3)                       # s1 to s8 of Laskar: -5.6, -7.1, -17.7, -18.8, -26.3 (first order, eight bodies)
    assert all(abs(m.frequency_rad_yr) > ov.MODE_FLOOR_RAD_YR for m in ov.all_sibling_modes(A_AU, MASS, SUN, INCLINATION, NODE)[EARTH])


def test_the_earth_with_its_moon_stays_inside_a_degree_or_so_of_its_measured_range():
    alpha, modes = _earth(True)
    variation = ov.analyse(alpha, 23.44, modes)
    assert not variation.chaotic and not variation.resonant
    assert 21.5 <= variation.minimum_deg <= 23.4 and 23.5 <= variation.maximum_deg <= 25.5          # measured 22.1 to 24.5
    assert variation.forced_deg < 2.0


def test_the_earth_without_its_moon_swings_several_times_as_far():
    with_moon, modes = _earth(True)
    without, _same = _earth(False)
    steady = ov.analyse(with_moon, 23.44, modes)
    wandering = ov.analyse(without, 23.44, modes)
    assert without < 0.4 * with_moon
    assert wandering.forced_deg > 5.0 * steady.forced_deg
    assert wandering.maximum_deg - wandering.minimum_deg > 5.0 * (steady.maximum_deg - steady.minimum_deg)


def test_mars_at_its_measured_ellipticity_sits_inside_a_resonance_with_a_wide_swing():
    alpha = ov.precession_constant_rad_yr(0.0048, 7.088e-5, [ov.Torquer(SUN, 1.5237 * AU, 0.0934)])
    variation = ov.analyse(alpha, 25.19, ov.all_sibling_modes(A_AU, MASS, SUN, INCLINATION, NODE)[MARS])
    assert variation.resonant and variation.maximum_deg - variation.minimum_deg > 15.0          # the full theory finds it chaotic, 0 to 60: this is the lower bound


def test_a_lone_planet_keeps_its_tilt_and_a_planet_with_no_torque_follows_the_nod_of_its_orbit():
    lone = ov.analyse(1.0e-4, 23.0, ())
    assert (lone.minimum_deg, lone.maximum_deg, lone.chaotic) == (23.0, 23.0, False)
    modes = (ov.Mode(-2.0e-5, math.radians(2.0)), ov.Mode(-5.0e-5, math.radians(1.0)))
    fixed = ov.analyse(0.0, 40.0, modes)
    assert fixed.forced_deg == pytest.approx(math.degrees(math.hypot(math.radians(2.0), math.radians(1.0))))


def test_two_resonances_close_enough_to_overlap_make_the_zone_chaotic_and_far_apart_ones_do_not():
    alpha = 1.0e-4
    near = (ov.Mode(-0.50e-4, math.radians(1.0)), ov.Mode(-0.56e-4, math.radians(1.0)))
    far = (ov.Mode(-0.20e-4, math.radians(0.2)), ov.Mode(-0.85e-4, math.radians(0.2)))
    center = math.degrees(math.acos(0.53))
    chaotic = ov.analyse(alpha, center, near)
    assert chaotic.chaotic and chaotic.resonant
    assert chaotic.minimum_deg < center < chaotic.maximum_deg
    lone = ov.analyse(alpha, math.degrees(math.acos(0.2)), far)
    assert lone.resonant and not lone.chaotic
    assert not ov.analyse(alpha, 5.0, far).resonant


def test_the_resonance_sits_where_the_axis_rate_matches_the_node_rate_and_its_width_follows_the_pendulum():
    alpha, mode = 1.0e-4, ov.Mode(-4.0e-5, math.radians(1.0))
    (resonance,) = ov.resonances_of(alpha, (mode,))
    assert resonance.cos_center == pytest.approx(0.4)
    sine = math.sqrt(1.0 - 0.16)
    assert resonance.cos_half_width == pytest.approx(2.0 * math.sqrt(2.0 * math.radians(1.0) * 4.0e-5 * sine / alpha))
    assert ov.resonances_of(alpha, (ov.Mode(-2.0e-4, 0.01),)) == ()          # |s| above alpha: no tilt meets it
    assert ov.resonances_of(0.0, (mode,)) == ()


def test_a_tilt_drawn_for_a_chaotic_planet_stays_in_its_range_and_is_even_in_the_cosine():
    variation = ov.Variation(10.0, 60.0, True, True, 0.0, (), 1.0e-4)
    assert ov.chaotic_tilt_deg(variation, 0.0) == pytest.approx(10.0)
    assert ov.chaotic_tilt_deg(variation, 1.0) == pytest.approx(60.0)
    assert math.cos(math.radians(ov.chaotic_tilt_deg(variation, 0.5))) == pytest.approx(0.5 * (math.cos(math.radians(10.0)) + math.cos(math.radians(60.0))))
    with pytest.raises(ValueError):
        ov.chaotic_tilt_deg(variation, 1.5)


def test_a_companion_star_drives_a_mode_only_when_its_orbit_is_inclined_and_far_enough():
    assert ov.companion_mode(1.0, SUN, 0.5 * SUN, 200.0, 0.3, 0.0) is None
    assert ov.companion_mode(1.0, SUN, 0.5 * SUN, 3.0, 0.3, 0.5) is None                       # inside five times the planet's distance: not a perturbation
    mode = ov.companion_mode(1.0, SUN, 0.5 * SUN, 200.0, 0.3, math.radians(20.0))
    assert mode.frequency_rad_yr < 0.0 and mode.amplitude_rad == pytest.approx(math.radians(20.0))
    n = math.sqrt(constants.G * SUN / AU ** 3) * constants.SECONDS_PER_YEAR
    assert mode.frequency_rad_yr == pytest.approx(-0.75 * n * 0.5 * (1.0 / 200.0) ** 3 / (1.0 - 0.09) ** 1.5 * math.cos(math.radians(20.0)))


def test_bad_inputs_are_refused():
    with pytest.raises(ValueError):
        ov.fluid_love_number(0.0)
    with pytest.raises(ValueError):
        ov.dynamical_ellipticity(-1.0, 1.0, 1.0, 0.33)
    with pytest.raises(ValueError):
        ov.analyse(1.0, 200.0, ())
    with pytest.raises(ValueError):
        ov.precession_constant_rad_yr(0.003, 7e-5, [ov.Torquer(1.0, 0.0)])
