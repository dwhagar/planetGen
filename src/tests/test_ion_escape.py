# tests/test_ion_escape.py

"""The stellar wind's ion loss and the field that closes it (GEN.804): the open share of a dipole, the Venus and Mars rates it is tuned to, and the
retention of two otherwise equal worlds that differ only in their field."""

import math

import pytest

from planetgen.physics import activity, constants, ion_escape, thermosphere as ts

ME = constants.EARTH_MASS_TO_KG
MARS = dict(mass=0.107 * ME, radius=3.3895e6, teq=210.0)
VENUS = dict(mass=0.815 * ME, radius=6.0518e6, teq=232.0)
KG_PER_O_ION = 16.0 * 1.6605e-27


def _energy_limited(body, xuv_earth):
    return ts.energy_limited_mass_loss_kg_s(xuv_earth * constants.EARTH_XUV_FLUX_W_M2, 1.05 * body["radius"], body["mass"], efficiency=1.0)


def test_the_open_share_is_all_without_a_field_and_the_polar_caps_of_a_dipole_with_one():
    assert ion_escape.open_share(None) == 1.0 and ion_escape.open_share(1.0) == 1.0 and ion_escape.open_share(0.4) == 1.0
    assert ion_escape.open_share(3.0) == pytest.approx(1.0 - math.sqrt(2.0 / 3.0))
    assert ion_escape.open_share(9.75) == pytest.approx(0.052, abs=0.002)       # Earth
    assert ion_escape.open_share(30.0) < ion_escape.open_share(9.75) < ion_escape.open_share(3.0) < 1.0


def test_venus_and_mars_lose_the_ions_the_spacecraft_saw_to_within_a_factor_of_three():
    """Venus Express: about 1e25 ions/s (0.1 to 0.3 kg/s of H+ and O+ by number and mass); MAVEN: 1e24 to 1e25 at Mars (second-hand)."""
    venus = ion_escape.ion_loss_kg_s(_energy_limited(VENUS, 1.9))
    mars = ion_escape.ion_loss_kg_s(_energy_limited(MARS, 0.43))
    assert 0.2 / 3.0 < venus < 0.2 * 3.0
    assert 0.13 / 3.0 < mars < 0.13 * 3.0


def test_a_field_not_yet_drawn_leaves_the_ion_loss_out_and_a_drawn_one_with_no_field_puts_it_in():
    assert _thin_co2_world()["CO2"] == pytest.approx(_thin_co2_world(magnetopause_rp=1e9)["CO2"], rel=1e-3)
    assert _thin_co2_world(magnetopause_rp=0.0)["CO2"] < _thin_co2_world()["CO2"]


def test_no_energy_or_a_negative_one_loses_nothing():
    assert ion_escape.ion_loss_kg_s(0.0) == 0.0 and ion_escape.ion_loss_kg_s(-5.0, 4.0) == 0.0


def _thin_co2_world(**kwargs):
    air = {"CO2": 0.95, "N2": 0.03, "Ar": 0.02}
    return ts.kept_shares(air, MARS["mass"], MARS["radius"], MARS["teq"], 600.0, 4.5, 0.43, 1.0, constants.SOLAR_LUMINOSITY, **kwargs)[0]


def test_two_worlds_that_differ_only_in_their_field_keep_different_air():
    """GEN.804: a thin air under the wind with no field loses more than the same air behind a magnetopause ten radii out."""
    open_world, closed_world = _thin_co2_world(magnetopause_rp=0.0), _thin_co2_world(magnetopause_rp=10.0)
    assert open_world["CO2"] < closed_world["CO2"]
    assert closed_world["CO2"] > 0.5 > open_world["CO2"]
    assert _thin_co2_world(magnetopause_rp=3.0)["CO2"] < closed_world["CO2"]


def test_a_thick_air_hardly_notices_the_wind_and_a_field_makes_it_nothing():
    air = {"N2": 0.78, "O2": 0.21, "Ar": 0.01}
    args = (air, ME, 6.371e6, 255.0, 101325.0, 4.5, 1.0, 1.0, constants.SOLAR_LUMINOSITY)
    assert ts.kept_shares(*args, magnetopause_rp=0.0)[0]["N2"] > 0.97
    assert ts.kept_shares(*args, magnetopause_rp=9.75)[0]["N2"] > 0.99


def test_a_star_that_kept_its_saturation_longer_took_more_from_a_young_world():
    """The star's own `log_lx_lbol` at its age sets its decay: the same star now, with a history that was brighter earlier, costs the air more."""
    mass_sol, age, luminosity = 0.7, 4.0, 0.2 * constants.SOLAR_LUMINOSITY
    t_sat, decay = activity.saturation_age_gy(mass_sol), activity.decay_exponent(mass_sol)
    saturated = activity.saturated_lx_lbol(mass_sol)
    nominal = math.log10(activity.lx_lbol_at(age, saturated, t_sat, decay))
    slow = ts.xuv_ratio(mass_sol, luminosity, age, 0.1, nominal - 0.5)
    fast = ts.xuv_ratio(mass_sol, luminosity, age, 0.1, nominal + 0.3)
    assert ts.xuv_ratio(mass_sol, luminosity, age, 0.1) == pytest.approx(ts.xuv_ratio(mass_sol, luminosity, age, 0.1, nominal), rel=1e-9)
    assert fast < slow          # the star that is more active now than its table was less active before, relative to now


def test_the_effective_saturation_age_is_the_nominal_one_at_the_nominal_activity_and_moves_with_the_star():
    mass_sol, age = 0.7, 4.0
    t_sat, decay = activity.saturation_age_gy(mass_sol), activity.decay_exponent(mass_sol)
    nominal = math.log10(activity.lx_lbol_at(age, activity.saturated_lx_lbol(mass_sol), t_sat, decay))
    assert activity.effective_saturation_age_gy(mass_sol, age, None) == t_sat
    assert activity.effective_saturation_age_gy(mass_sol, age, nominal) == pytest.approx(t_sat)
    assert activity.effective_saturation_age_gy(mass_sol, age, nominal + 0.3) > t_sat
    assert activity.effective_saturation_age_gy(mass_sol, age, nominal - 0.3) < t_sat
    assert activity.effective_saturation_age_gy(mass_sol, 0.01, math.log10(activity.saturated_lx_lbol(mass_sol))) >= 0.01
