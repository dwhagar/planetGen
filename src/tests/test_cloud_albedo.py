"""
The cloud albedo of a world with an ocean (GEN.974, `physics/cloud_albedo.py`, design `docs/design/composition/cloud-albedo.md`): the amount term that thins the low cloud
as a world warms, the substellar deck of a world locked to its star, the estimate of the lock, and their wiring into the surface-air fixed point.
"""

import types

import pytest

from planetgen.physics import climate_branches as cb, cloud_albedo as ca, surface_air as sa


def test_the_cloud_feedback_of_the_earth_at_288_k_is_the_assessments_central_value():
    """`BETA` is sized so that `-BETA x S / 4` is +0.42 W/m2 per K (AR6 table 7.10, inside the very likely range -0.10 to +0.94), a pure-module test like the ice law's."""
    feedback = ca.BETA_PER_K * ca.REFERENCE_SOLAR_CONSTANT_W_M2 / 4.0       # the albedo lost per kelvin times the light a square metre of the sphere takes
    assert feedback == pytest.approx(0.42)
    assert -0.10 < feedback < 0.94


def test_a_body_at_the_reference_temperature_keeps_the_albedo_it_drew():
    assert ca.amount_shift(ca.REFERENCE_TEMPERATURE_K) == 0.0


def test_the_amount_falls_as_the_world_warms_and_is_bounded():
    temperatures = [150.0, 250.0, 288.0, 320.0, 400.0, 600.0]
    shifts = [ca.amount_shift(t) for t in temperatures]
    assert shifts == sorted(shifts, reverse=True)
    assert shifts[0] == ca.DELTA_MAX and shifts[-1] == -ca.DELTA_MAX
    assert ca.amount_shift(320.0) == pytest.approx(-ca.BETA_PER_K * 32.0)


def test_the_deck_reproduces_the_albedos_of_the_locked_ocean_planet_at_its_two_ends():
    """[R0914] table 1 (an M-star spectrum): 0.42 at 1,000 W/m2 and 0.53 at 2,200, so the deck adds 0.11 between them."""
    assert ca.deck_shift(1000.0) == 0.0 and ca.deck_shift(500.0) == 0.0
    assert ca.deck_shift(2200.0) == pytest.approx(0.53 - 0.42)
    assert ca.deck_shift(1.0e5) == ca.DECK_DELTA_MAX


def test_only_a_locked_body_has_a_deck_and_only_with_an_ocean():
    free, locked = ca.CloudLaw(2200.0, locked=False), ca.CloudLaw(2200.0, locked=True)
    assert locked.shift(288.0) - free.shift(288.0) == pytest.approx(ca.deck_shift(2200.0))
    assert free.shift(250.0, ocean=False) == 0.0 and locked.shift(250.0, ocean=False) == 0.0


def test_the_lock_estimate_follows_the_tide():
    """A hot Earth round a small star is locked in the age of the Galaxy, the Earth round the Sun is not (and the Earth under the Moon takes about 5e10 years by this form)."""
    gyr = 3.15576e16
    earth_kg, earth_m, sun_kg, au_m = 5.972e24, 6.371e6, 1.989e30, 1.496e11
    assert not ca.tidally_locked(earth_kg, earth_m, sun_kg, au_m, 4.5 * gyr)
    assert ca.tidally_locked(earth_kg, earth_m, 0.2 * sun_kg, 0.05 * au_m, 4.5 * gyr)
    assert ca.despin_time_s(earth_kg, earth_m, 7.35e22, 3.84e8) / 3.15576e7 == pytest.approx(5.0e10, rel=0.3)
    assert not ca.tidally_locked(0.0, earth_m, sun_kg, au_m, gyr)
    assert ca.despin_time_s(earth_kg, earth_m, sun_kg, 0.1 * au_m) < ca.despin_time_s(earth_kg, earth_m, sun_kg, 0.2 * au_m)


def _setting(cloud_law):
    return sa.Setting(gravity_m_s2=9.8, radius_m=6.371e6, mantle_kg=4.0e24, released_j=0.0, heat_flux_w_m2=0.087, ocean_cover=0.71, water_kg=1.4e21,
                      absorbed_flux_w_m2=240.0, released_at_seal_j=0.0, subject_of=lambda *_: None, starlight_w_m2=1361.0, albedo_law=cb.IceAlbedoLaw(),
                      cloud_law=cloud_law)


def test_the_albedo_of_a_step_is_the_ices_and_the_clouds_and_the_flux_follows_it():
    held, clouded = _setting(None), _setting(ca.CloudLaw(1361.0))
    assert held.albedo_at(300.0) == pytest.approx(cb.IceAlbedoLaw().albedo(300.0))
    assert clouded.albedo_at(300.0) == pytest.approx(cb.IceAlbedoLaw().albedo(300.0) + ca.amount_shift(300.0))
    assert clouded.albedo_at(300.0) < held.albedo_at(300.0) and clouded.absorbed_flux_at(300.0) > held.absorbed_flux_at(300.0)
    assert clouded.albedo_at(300.0, ocean=False) == pytest.approx(held.albedo_at(300.0))      # no liquid ocean at the step, no cloud
    assert sa.Setting.albedo_at(types.SimpleNamespace(albedo_law=None, cloud_law=None), 300.0) is None
    assert _setting(ca.CloudLaw(1361.0)).albedo_at(150.0) <= ca.MAX_ALBEDO


def test_a_setting_of_a_wet_body_has_a_cloud_law_and_a_dry_one_or_a_giant_has_none():
    body = types.SimpleNamespace(body_type="t", is_moon=False, mass=5.972e24, radius=6371.0,
                                 star=types.SimpleNamespace(age=4.5, mass=1.989e30, temperature=5772.0))
    wet = sa.cloud_law_of(body, 1361.0, 1.0)
    assert wet.starlight_w_m2 == 1361.0 and wet.locked is False
    close = sa.cloud_law_of(types.SimpleNamespace(**{**vars(body), "star": types.SimpleNamespace(age=4.5, mass=0.2 * 1.989e30, temperature=3300.0)}), 1361.0, 0.05)
    assert close.locked is True
    moon = sa.cloud_law_of(types.SimpleNamespace(**{**vars(body), "is_moon": True}), 1361.0, 1.0)
    assert moon.locked is True
    assert sa.ice_albedo_law(types.SimpleNamespace(body_type="g", mass=1.0e27, star=None), 0.3, 0.1) is None
