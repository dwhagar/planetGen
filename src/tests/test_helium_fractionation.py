# tests/test_helium_fractionation.py

"""Helium enrichment of an envelope that loses hydrogen first (GEN.1045): Hu et al. 2015 and [R2839]."""

import math

from planetgen.physics import constants, envelope_escape, helium_fractionation as hf

MASS = 10.0 * constants.EARTH_MASS_TO_KG                     # [R2839: Figure 3]: 10 Earth masses, transit radius 2.14
RADIUS = 2.14 * constants.EARTH_RADIUS_KM * 1e3 * hf.HOMOPAUSE_RADIUS_OVER_TRANSIT
ENVELOPE_KG = 3.0e23                                         # 0.5 percent of 10 Earth masses, 3.0e26 g (their B.4)
YEAR = constants.SECONDS_PER_YEAR


def test_the_envelope_escape_model_keeps_a_solar_mix():
    # the audit: envelope_escape returns a mass only, so nothing in it can make a helium-rich air
    assert not any("helium" in name.lower() for name in dir(envelope_escape))
    assert isinstance(envelope_escape.envelope_left_earth(6.0, 0.01, 50.0, 1.0, 1.0, 10.0), float)


def test_the_diffusion_coefficient_is_the_papers():
    # Hu et al. 2015 Sec. 2.1: 10^4 K and an ionised share of 0.1 give b' = 8.0e20 cm^-1 s^-1
    assert math.isclose(hf.effective_diffusion_m(1.0e4, 0.1) / 100.0, 8.0e20, rel_tol=0.01)
    assert math.isclose(hf.binary_diffusion_m(2.0e4) / hf.binary_diffusion_m(1.0e4), 2 ** 0.732, rel_tol=1e-9)


def test_a_flux_under_the_critical_one_is_all_hydrogen():
    critical = hf.critical_mass_flux_kg_s(MASS, RADIUS, 0.245)
    assert 5e4 < critical < 1e6               # a few 1e8 g/s for a 10 Earth-mass planet, the range of the rates in Malsky et al. B.4
    assert hf.partition(0.9 * critical, MASS, RADIUS, 0.245) == (0.9 * critical, 0.0)
    hydrogen, helium = hf.partition(3.0 * critical, MASS, RADIUS, 0.245)
    assert helium > 0.0 and math.isclose(hydrogen + helium, 3.0 * critical, rel_tol=1e-12)
    # far over it the two leave in the envelope's proportions (no fractionation)
    hydrogen, helium = hf.partition(1.0e4 * critical, MASS, RADIUS, 0.245)
    assert math.isclose(helium / (hydrogen + helium), 0.245, abs_tol=0.01)


def test_fractionation_is_strongest_near_the_critical_flux():
    ratios = []
    for multiple in (1.5, 3.0, 10.0, 100.0):
        critical = hf.critical_mass_flux_kg_s(MASS, RADIUS, 0.245)
        _, y = hf.evolve(ENVELOPE_KG, 0.245, MASS, RADIUS, multiple * critical, 0.5 * ENVELOPE_KG / (multiple * critical))
        ratios.append(y)
    assert ratios == sorted(ratios, reverse=True)
    assert ratios[0] > 0.35 and ratios[-1] < 0.26


def test_over_ninety_percent_loss_leaves_a_helium_dominated_air_near_the_critical_flux():
    critical = hf.critical_mass_flux_kg_s(MASS, RADIUS, 0.245)
    left, y = hf.evolve(ENVELOPE_KG, 0.245, MASS, RADIUS, 1.5 * critical, 0.9 * ENVELOPE_KG / (1.5 * critical))
    assert math.isclose(left / ENVELOPE_KG, 0.1, rel_tol=0.02)
    assert y > 0.8 and hf.is_helium_enhanced(y)


def test_malsky_b4_example_loses_half_the_envelope_and_gains_helium():
    # 1e9 g/s (1e6 kg/s) for 5 Gyr takes just over half of a 0.5 percent envelope; the helium share rises, not to 0.4 at this rate
    left, y = hf.evolve(ENVELOPE_KG, 0.245, MASS, RADIUS, 1.0e6, 5.0e9 * YEAR)
    assert 0.4 < left / ENVELOPE_KG < 0.55
    assert 0.25 < y < 0.4


def test_nothing_is_lost_nothing_changes():
    assert hf.evolve(ENVELOPE_KG, 0.245, MASS, RADIUS, 0.0, 1.0e9 * YEAR) == (ENVELOPE_KG, 0.245)
    assert hf.partition(0.0, MASS, RADIUS, 0.245) == (0.0, 0.0)


def test_a_history_from_the_escape_model_gives_a_helium_rise_with_the_loss():
    ages = (1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0)
    results = {}
    for core, f_env, insolation in ((10.0, 0.01, 30.0), (6.0, 0.01, 50.0), (4.0, 0.02, 30.0)):
        history = [(age, envelope_escape.envelope_left_earth(core, f_env, insolation, 1.0, 1.0, age)) for age in ages]
        left, y = hf.helium_fraction_from_history(core, history, insolation)
        results[(core, f_env)] = (1.0 - left / history[0][1], y)
    lost_a, y_a = results[(10.0, 0.01)]
    lost_b, y_b = results[(6.0, 0.01)]
    lost_c, y_c = results[(4.0, 0.02)]
    assert lost_a < lost_b < lost_c
    assert 0.245 <= y_a < min(y_b, y_c)       # the more the planet loses, the richer in helium; the rate sets how much
    assert y_b > 0.3 and y_c > 0.3
    # a colder homopause fractionates more
    history = [(age, envelope_escape.envelope_left_earth(4.0, 0.02, 30.0, 1.0, 1.0, age)) for age in ages]
    assert hf.helium_fraction_from_history(4.0, history, 30.0, temperature_k=3000.0)[1] > y_c


def test_a_stripped_planet_is_a_bare_core():
    history = [(1.0, 0.04), (2.0, 0.02), (3.0, 0.0)]
    assert hf.helium_fraction_from_history(4.0, history, 100.0) == (0.0, hf.SOLAR_ENVELOPE_HELIUM)


def test_number_fractions():
    x_h, x_he = hf.number_fractions(0.245)
    assert math.isclose(x_h + x_he, 1.0)
    assert 0.9 < x_h < 0.93
    assert hf.number_fractions(0.0) == (1.0, 0.0)
