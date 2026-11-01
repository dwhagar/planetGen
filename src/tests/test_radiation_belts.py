"""
Trapped radiation belts (GEN.749, `planetgen.physics.radiation_belts`): the belts of a dipole checked against the two bodies whose belts have been
flown through, Earth and Jupiter, and the rules that bound them.
"""

import math
import types

import pytest

from planetgen.physics import constants, magnetism, radiation_belts as rb

EARTH_RADIUS_M = constants.EARTH_RADIUS_KM * 1000.0
JUPITER_RADIUS_M = constants.JUPITER_RADIUS_KM * 1000.0
JUPITER_MOMENT = 1.5e27
EARTH_SIDEREAL_S = 86164.0
JUPITER_SPIN_S = 9.925 * 3600.0
JUPITER_WIND_PA = magnetism.EARTH_WIND_PA / 5.2 ** 2


def _earth():
    return rb.belts_for(magnetism.EARTH_MOMENT_A_M2, EARTH_RADIUS_M, constants.EARTH_MASS_TO_KG, EARTH_SIDEREAL_S,
                        magnetism.EARTH_WIND_PA, 8500.0)


def _jupiter():
    return rb.belts_for(JUPITER_MOMENT, JUPITER_RADIUS_M, constants.JUPITER_MASS_TO_KG, JUPITER_SPIN_S, JUPITER_WIND_PA, None)


def _belt(belts, name):
    return next(belt for belt in belts if belt.name == name)


def test_earth_standoff_and_outer_edge_are_the_magnetopause_scaled():
    assert rb.standoff_m(magnetism.EARTH_MOMENT_A_M2, magnetism.EARTH_WIND_PA) / EARTH_RADIUS_M == pytest.approx(9.75)
    assert _belt(_earth(), "outer").outer_l == pytest.approx(0.75 * 9.75)


def test_earth_belts_sit_where_they_were_flown_through():
    inner, outer = _belt(_earth(), "inner"), _belt(_earth(), "outer")
    assert 1.5 <= inner.peak_l <= 2.0
    assert 4.0 <= outer.peak_l <= 5.0
    assert 1.03 <= inner.inner_l <= 1.15
    assert 6.5 <= outer.outer_l <= 8.0
    assert inner.peak_l < outer.peak_l


def test_earth_has_a_slot_between_the_belts():
    belts = _earth()
    assert rb.dose_gy_yr(belts, 2.6) < 0.2 * rb.dose_gy_yr(belts, _belt(belts, "inner").peak_l)
    assert rb.dose_gy_yr(belts, 2.6) < 0.2 * rb.dose_gy_yr(belts, _belt(belts, "outer").peak_l)


@pytest.mark.parametrize("name, l_shell, low_gy, high_gy", [
    ("LEO 500 km", 1.08, 1.0, 100.0),
    ("MEO navigation 20 000 km", 4.14, 1000.0, 1e4),
    ("GEO", 6.6, 10.0, 500.0),
])
def test_earth_orbit_doses_are_the_published_orders_of_magnitude(name, l_shell, low_gy, high_gy):
    assert low_gy <= rb.dose_gy_yr(_earth(), l_shell) <= high_gy


def test_nothing_is_trapped_outside_the_outer_edge_or_inside_the_air():
    belts = _earth()
    assert rb.dose_gy_yr(belts, 8.0) == 0.0
    assert rb.dose_gy_yr(belts, 1.0) == 0.0


def test_jupiter_is_spin_limited_and_earth_is_wind_limited():
    assert _belt(_jupiter(), "outer").outer_l == pytest.approx(9.0 * rb.synchronous_l(constants.JUPITER_MASS_TO_KG, JUPITER_RADIUS_M, JUPITER_SPIN_S))
    assert 2.1 <= rb.synchronous_l(constants.JUPITER_MASS_TO_KG, JUPITER_RADIUS_M, JUPITER_SPIN_S) <= 2.4
    assert _belt(_jupiter(), "outer").outer_l < 0.75 * rb.standoff_m(JUPITER_MOMENT, JUPITER_WIND_PA) / JUPITER_RADIUS_M
    assert _belt(_earth(), "outer").outer_l == pytest.approx(0.75 * 9.75)


def test_jupiter_peak_is_inside_io_and_its_belts_reach_past_ganymede_not_callisto():
    outer = _belt(_jupiter(), "outer")
    assert 4.0 <= outer.peak_l <= 6.0
    assert 15.0 <= outer.outer_l <= 25.0
    assert rb.dose_gy_yr(_jupiter(), 26.3) == 0.0
    assert rb.dose_gy_yr(_jupiter(), 15.0) > 0.0


def test_jupiter_moon_doses_fall_with_distance_as_io_to_europa_near_seven():
    belts = _jupiter()
    io, europa, ganymede = (rb.dose_gy_yr(belts, l) for l in (5.9, 9.4, 15.0))
    assert io > europa > ganymede > 0.0
    assert 6.7 / 3.0 <= io / europa <= 6.7 * 3.0


def test_jupiter_belt_is_thousands_of_times_earths():
    ratio = max(b.peak_gy_yr for b in _jupiter()) / max(b.peak_gy_yr for b in _earth())
    assert 1e2 <= ratio <= 1e5


def test_no_field_or_no_wind_has_no_belts():
    assert rb.belts_for(0.0, EARTH_RADIUS_M, constants.EARTH_MASS_TO_KG, EARTH_SIDEREAL_S, magnetism.EARTH_WIND_PA) == ()
    assert rb.belts_for(None, EARTH_RADIUS_M, constants.EARTH_MASS_TO_KG, EARTH_SIDEREAL_S, magnetism.EARTH_WIND_PA) == ()
    assert rb.belts_for(magnetism.EARTH_MOMENT_A_M2, EARTH_RADIUS_M, constants.EARTH_MASS_TO_KG, EARTH_SIDEREAL_S, 0.0) == ()


def test_a_field_the_wind_presses_onto_the_air_has_no_belts():
    weak = 1e-3 * magnetism.EARTH_MOMENT_A_M2
    assert rb.belts_for(weak, EARTH_RADIUS_M, constants.EARTH_MASS_TO_KG, EARTH_SIDEREAL_S, 1e3 * magnetism.EARTH_WIND_PA, 8500.0) == ()


def test_a_stronger_wind_shrinks_the_belts_and_a_stronger_field_raises_them():
    calm = _belt(_earth(), "outer")
    gale = _belt(rb.belts_for(magnetism.EARTH_MOMENT_A_M2, EARTH_RADIUS_M, constants.EARTH_MASS_TO_KG, EARTH_SIDEREAL_S,
                              100.0 * magnetism.EARTH_WIND_PA, 8500.0), "outer")
    strong = _belt(rb.belts_for(10.0 * magnetism.EARTH_MOMENT_A_M2, EARTH_RADIUS_M, constants.EARTH_MASS_TO_KG, EARTH_SIDEREAL_S,
                                magnetism.EARTH_WIND_PA, 8500.0), "outer")
    assert gale.outer_l < calm.outer_l
    assert strong.peak_gy_yr == pytest.approx(calm.peak_gy_yr * 10.0 ** rb.INTENSITY_EXPONENT)


def test_a_slow_rotator_is_bounded_by_the_wind_not_the_spin():
    slow = rb.belts_for(magnetism.EARTH_MOMENT_A_M2, EARTH_RADIUS_M, constants.EARTH_MASS_TO_KG, 30.0 * 86400.0, magnetism.EARTH_WIND_PA, None)
    assert _belt(slow, "outer").outer_l == pytest.approx(0.75 * 9.75)


def test_profile_is_zero_at_the_edges_and_one_at_the_peak():
    assert rb.profile(1.0, 1.0, 3.0, 7.0) == 0.0
    assert rb.profile(7.0, 1.0, 3.0, 7.0) == 0.0
    assert rb.profile(3.0, 1.0, 3.0, 7.0) == pytest.approx(1.0)
    assert 0.0 < rb.profile(5.0, 1.0, 3.0, 7.0) < 1.0


def test_belts_of_body_reads_the_stored_values():
    star = types.SimpleNamespace(radius=None, luminosity=None, log_lx_lbol=None, pulsar_type=None)
    body = types.SimpleNamespace(
        magnetic_moment_a_m2=magnetism.EARTH_MOMENT_A_M2, body_type="t", is_moon=False, rotation_period_hours=23.93,
        radius=constants.EARTH_RADIUS_KM, mass=constants.EARTH_MASS_TO_KG, star=star, _star_distance_au=1.0,
        surface_temperature=None, gravity=None, atmospheric_pressure=0.0)
    belts = rb.belts_of_body(body)
    assert len(belts) == 2
    assert _belt(belts, "outer").outer_l == pytest.approx(0.75 * magnetism.magnetopause_rp(body, 1.0))
    body.magnetic_moment_a_m2 = None
    assert rb.belts_of_body(body) == ()
