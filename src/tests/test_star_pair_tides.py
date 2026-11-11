# tests/test_star_pair_tides.py

"""GEN.1467: the tide of a stellar partner on a star's spin, and the apsidal precession of the pair's orbit (`physics/star_pair_tides.py`)."""

import math
from types import SimpleNamespace

import pytest

from planetgen.physics import constants, stability
from planetgen.physics import star_pair_tides as st

SUN_KG = constants.SOLAR_MASS_TO_KG
SUN_KM = 695700.0


def _star(mass_sol=1.0, radius_sol=1.0, temperature=5772.0, cls="V", rotation_hours=25.4 * 24.0):
    return SimpleNamespace(mass=mass_sol * SUN_KG, radius=radius_sol * SUN_KM, temperature=temperature, yerkes_class=cls, rotation_period_hours=rotation_hours)


def _axis_au(period_days, total_mass_sol):
    return ((period_days / 365.25) ** 2 * total_mass_sol) ** (1.0 / 3.0)


def test_pseudo_synchronous_rate_is_the_mean_motion_for_a_circular_orbit_and_grows_with_eccentricity():
    assert st.pseudo_synchronous_ratio(0.0) == 1.0
    assert st.pseudo_synchronous_ratio(0.3) > st.pseudo_synchronous_ratio(0.1) > 1.0
    assert st.pseudo_synchronous_ratio(0.5) == pytest.approx(2.805, abs=0.01)       # equation 42 of the tidal theory (R5229) at e = 0.5


def test_the_synchronisation_time_falls_steeply_with_the_separation():
    star = _star()
    near, far = (st.sync_time_yr(star, SUN_KG, _axis_au(p, 2.0)) for p in (5.0, 10.0))
    assert near < far
    # convective tide: t ~ a^6 / f_conv at a given spin; doubling the period multiplies a^6 by 2^4, and f_conv (the tide slower than the eddy turnover is damped) rises with the period
    # about eightfold from 5 to 10 days for a Sun with a 25 day spin, so the time grows about twofold
    assert 1.5 < far / near < 40.0
    assert st.sync_time_yr(star, SUN_KG, _axis_au(100.0, 2.0)) > 1000.0 * far


def test_a_solar_type_pair_locks_inside_ten_days_and_not_at_a_month():
    star = _star()
    age = 4.5e9
    close = st.spin_state(star, SUN_KG, _axis_au(5.0, 2.0), 0.0, age)
    ten = st.spin_state(star, SUN_KG, _axis_au(10.0, 2.0), 0.0, age)
    wide = st.spin_state(star, SUN_KG, _axis_au(100.0, 2.0), 0.0, age)
    assert close.locked and close.rotation_period_hours == pytest.approx(5.0 * 24.0, rel=1e-3) and close.obliquity_factor < 1e-6
    assert ten.locked
    assert not wide.locked and not wide.changed and wide.rotation_period_hours == star.rotation_period_hours and wide.obliquity_factor == 1.0


def test_a_radiative_star_locks_only_much_closer():
    hot = _star(mass_sol=2.5, radius_sol=2.3, temperature=9500.0, rotation_hours=30.0)
    assert not st.has_convective_envelope(hot)
    age = 5.0e8
    two_days = st.spin_state(hot, 2.5 * SUN_KG, _axis_au(2.0, 5.0), 0.0, age)
    ten_days = st.spin_state(hot, 2.5 * SUN_KG, _axis_au(10.0, 5.0), 0.0, age)
    assert two_days.locked and not ten_days.locked and not ten_days.changed


def test_a_lone_compact_or_untouched_star_comes_back_unchanged():
    dwarf = _star(cls="DA", radius_sol=0.012, temperature=12000.0, rotation_hours=24.0)
    assert st.sync_time_yr(dwarf, SUN_KG, 0.05) == math.inf
    result = st.spin_state(dwarf, SUN_KG, 0.05, 0.0, 1e9)
    assert not result.changed and result.rotation_period_hours == 24.0
    assert st.sync_time_yr(_star(), 0.0, 1.0) == math.inf and st.sync_time_yr(_star(), SUN_KG, 0.0) == math.inf


def test_an_eccentric_orbit_settles_the_spin_above_the_mean_motion():
    star = _star()
    circular = st.spin_state(star, SUN_KG, _axis_au(4.0, 2.0), 0.0, 4.5e9)
    eccentric = st.spin_state(star, SUN_KG, _axis_au(4.0, 2.0), 0.4, 4.5e9)
    assert eccentric.rotation_period_hours < circular.rotation_period_hours


def test_the_spin_is_never_faster_than_the_breakup_cap():
    tiny = _star(mass_sol=0.1, radius_sol=0.12, temperature=2800.0, rotation_hours=100.0)
    result = st.spin_state(tiny, 0.1 * SUN_KG, 0.0008, 0.0, 1e9)
    cap = 2.0 * math.pi * tiny.radius / (0.85 * math.sqrt(2.0 / 3.0 * constants.G * tiny.mass / (tiny.radius * 1000.0)) / 1000.0) / 3600.0
    assert result.rotation_period_hours >= cap * (1.0 - 1e-9)


def test_the_apsidal_period_of_a_close_pair_is_the_size_observed_and_relativity_is_the_published_one():
    a = _axis_au(2.0, 2.0)
    star = _star(rotation_hours=48.0)
    pre = st.apsidal_precession(star, star, a, 0.0)
    assert pre.tide > 0.0 and pre.rotation > 0.0 and pre.relativity > 0.0
    assert 10.0 < pre.period_yr < 2000.0        # close main-sequence pairs turn their orbits in tens to thousands of years
    gr_period = stability.gr_precession_period_yr(a, 0.0, 2.0, stability.orbital_period_yr(a, 2.0))
    assert 2.0 * math.pi / pre.relativity == pytest.approx(gr_period, rel=1e-9)
    # a wide pair is dominated by relativity
    wide_a = _axis_au(1000.0, 2.0)
    wide = st.apsidal_precession(star, star, wide_a, 0.0)
    assert wide.relativity > 10.0 * (wide.tide + wide.rotation) and wide.relativity > 1e6 * wide.tide


def test_the_tidal_term_goes_as_the_fifth_power_of_the_radius_over_the_axis_and_grows_with_eccentricity():
    star = _star(rotation_hours=1e9)         # no rotational bulge
    near = st.apsidal_precession(star, star, 0.04, 0.0)
    far = st.apsidal_precession(star, star, 0.08, 0.0)
    # rate ~ n (R/a)^5 and n ~ a^-1.5: doubling the axis cuts the rate by 2^6.5
    assert near.tide / far.tide == pytest.approx(2.0 ** 6.5, rel=1e-6)
    assert st.apsidal_precession(star, star, 0.04, 0.3).tide > near.tide


def test_the_apsidal_constants_and_envelope_shares_follow_the_mass():
    assert st.apsidal_constant(_star(mass_sol=0.2)) == pytest.approx(st.APSIDAL_CONSTANT_FULLY_CONVECTIVE)
    assert st.apsidal_constant(_star(mass_sol=1.0)) == pytest.approx(st.APSIDAL_CONSTANT_SOLAR)
    assert st.apsidal_constant(_star(mass_sol=8.0)) == pytest.approx(st.APSIDAL_CONSTANT_MASSIVE)
    assert st.apsidal_constant(_star(cls="III")) == st.APSIDAL_CONSTANT_GIANT
    assert st.envelope_mass_share(_star(mass_sol=0.2)) == 1.0 and 0.02 < st.envelope_mass_share(_star(mass_sol=1.0)) < 0.04
