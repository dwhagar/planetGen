# tests/test_radiation_dose_field.py

"""GEN.1442: the XUV, flare, particle and far-ultraviolet doses and the loss of air and envelope read every radiating object of the system (the field of
`planetgen.physics.radiation_field`), and a lone star gives the numbers it always did."""

import math
from types import SimpleNamespace

import pytest

from planetgen.physics import aerosols, constants, envelope_escape, magnetism, radiation, radiation_field as rf, thermosphere

SUN_W = constants.SOLAR_LUMINOSITY
SUN_KG = constants.SOLAR_MASS_TO_KG


def star(l_xuv_w=2.0e-6 * SUN_W, flare=3.0, luminosity_w=SUN_W, mass_kg=SUN_KG, fluence=5.0e24):
    return SimpleNamespace(luminosity=luminosity_w, temperature=5772.0, mass=mass_kg, radius=6.957e8, l_xuv_w=l_xuv_w, xuv_fluence_j=fluence,
                           flare_n33_per_yr=flare, flare_alpha=None, log_lx_lbol=-5.7)


def pair(axis_au=5.0, **companion):
    host = star()
    rf.attach(host, [rf.Companion(star(**companion), 1, axis_au, 0.0)])
    return host


def body_of(host):
    return SimpleNamespace(star=host)


def test_a_lone_star_gives_the_exposure_it_always_did():
    host = star()
    body = body_of(host)
    magnetism.update_exposure(body, 1.3)
    assert body.xuv_flux_earth == pytest.approx(host.l_xuv_w / (4 * math.pi * (1.3 * constants.AU_M) ** 2) / constants.EARTH_XUV_FLUX_W_M2, rel=1e-12)
    assert body.flare_irradiation_index == host.flare_n33_per_yr / 1.3 ** 2


def test_an_equal_companion_raises_the_exposure_by_its_geometry():
    lone, together = body_of(star()), body_of(pair(5.0))
    magnetism.update_exposure(lone, 1.5)
    magnetism.update_exposure(together, 1.5)
    ratio = 1.0 + 1.0 / (5.0 ** 2 - 1.5 ** 2) * 1.5 ** 2
    for name in ("xuv_flux_earth", "xuv_exposure_index", "flare_irradiation_index"):
        assert getattr(together, name) == pytest.approx(getattr(lone, name) * ratio, rel=1e-5), name


def test_a_companion_at_a_hundred_au_barely_changes_the_particle_dose():
    host = pair(100.0)
    near = radiation.free_sep_msv_yr(host.flare_n33_per_yr, 1.0)
    assert radiation.companion_free_sep_msv_yr(host, 1.0) / near < 1.5e-4
    assert radiation.companion_free_sep_msv_yr(star(), 1.0) == 0.0


def test_the_particle_dose_of_a_companion_follows_the_square_root_of_its_flare_rate():
    quiet, active = pair(5.0, flare=1.0), pair(5.0, flare=9.0)
    assert radiation.companion_free_sep_msv_yr(active, 1.0) == pytest.approx(3.0 * radiation.companion_free_sep_msv_yr(quiet, 1.0), rel=1e-9)


def test_a_remnant_with_no_stored_xuv_adds_none():
    host = pair(5.0, l_xuv_w=None, fluence=None, flare=None)
    body = body_of(host)
    reference = body_of(star())
    magnetism.update_exposure(body, 1.5)
    magnetism.update_exposure(reference, 1.5)
    assert body.xuv_flux_earth == reference.xuv_flux_earth
    assert radiation.companion_free_sep_msv_yr(host, 1.5) == 0.0


def test_the_far_ultraviolet_adds_the_companion():
    lone, together = star(), pair(5.0)
    base = aerosols.far_uv_of(lone, 1.5)
    both = aerosols.far_uv_of(together, 1.5)
    assert both.ratio > base.ratio
    assert both.ratio == pytest.approx(base.ratio * (1.0 + 1.0 / (5.0 ** 2 - 1.5 ** 2) * 1.5 ** 2), rel=1e-5)
    assert both.chromospheric_share == pytest.approx(base.chromospheric_share, rel=1e-9)     # two equal stars have the same share


def test_the_envelope_a_planet_keeps_is_less_with_a_companions_xuv():
    host = pair(3.0)
    d = 0.4
    insolation = 1.0 / d ** 2
    lone = envelope_escape.envelope_left_earth(5.0, 0.1, insolation, 1.0, 1.0, 4.5)
    emitters = rf.xuv_emitters(host, d)
    total = insolation + rf.companion_insolation_sol(host, d)
    both = envelope_escape.envelope_left_earth(5.0, 0.1, total, 1.0, 1.0, 4.5, distance_au=d, emitters=emitters)
    assert both < lone
    assert envelope_escape.envelope_left_earth(5.0, 0.1, insolation, 1.0, 1.0, 4.5, distance_au=d, emitters=()) == lone


def test_the_xuv_history_mixes_the_decay_of_each_star():
    host = pair(5.0)
    companions = rf.xuv_emitters(host, 1.0)
    flux = 2.0 * rf.emitter_flux_earth(companions[0])
    mixed = thermosphere._mixed_history(0.5, flux, companions, 4.5, 1.0)
    own = thermosphere.xuv_ratio(1.0, SUN_W, 4.5, 1.0, -5.7)
    assert mixed == pytest.approx(0.5 * 0.5 + 0.5 * own, rel=1e-9)
    assert own > 0.5, "a younger moment is more active, so the companion's ratio is above the host's stand-in value here"


def test_the_subject_of_a_body_carries_the_hosts_distance_only_with_companions():
    host = pair(5.0)
    emitters = rf.xuv_emitters(host, 1.2)
    subject = {"host_mass_solar": 1.0, "age_gy": 4.5, "luminosity_solar": 1.0, "insolation": 1.0 / 1.2 ** 2 + rf.companion_insolation_sol(host, 1.2),
               "host_distance_au": 1.2, "xuv_companions": tuple(emitters)}
    lone_flux = thermosphere.xuv_flux_earth({k: v for k, v in subject.items() if k in ("host_mass_solar", "age_gy", "luminosity_solar")} | {"insolation": 1.0 / 1.2 ** 2})
    assert thermosphere.xuv_flux_earth(subject) > lone_flux
