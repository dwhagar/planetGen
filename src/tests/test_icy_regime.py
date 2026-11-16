# tests/test_icy_regime.py

"""
GEN.1755 and GEN.1756: which tide model the spin of a water-rich body meets (`body_spin.icy_regime`): the ice shell over an ocean, ice solid through, or no ice at the surface (the rock under a sea).
"""

import types

import pytest

from planetgen.physics import body_spin, icy_tide, interior_heat, spin_torques, tidal_rheology

EUROPA = dict(mass=4.8e22, radius=1560.8, body_type="t", water_mass_fraction=0.5, surface_temperature=100.0, internal_heat_flux_w_m2=7.0e-4, tidal_flux_w_m2=0.0174)
SMALL_COLD = dict(mass=1.0e19, radius=100.0, body_type="t", water_mass_fraction=0.3, surface_temperature=60.0, internal_heat_flux_w_m2=1.0e-3, tidal_flux_w_m2=0.0)
WARM_WET = dict(mass=2.86e24, radius=5409.7, body_type="t", water_mass_fraction=0.0586, surface_temperature=346.0, internal_heat_flux_w_m2=0.087, tidal_flux_w_m2=0.0)
DRY = dict(EUROPA, water_mass_fraction=0.0)


def _regime(fields):
    body = types.SimpleNamespace(**fields)
    return body, interior_heat.icy_regime(body)


def test_a_shell_thinner_than_the_radius_is_the_shell_model():
    """Europa's 30 km shell over its ocean is what `icy_tide` is made for."""
    assert _regime(EUROPA)[1] == interior_heat.SHELL


def test_a_body_frozen_through_is_not_given_the_shell_models_k2_over_q():
    """GEN.1755: the shell of a cold 100 km body would be 860 km thick, more than the body is; the shell formula gave it a k2/Q of order one, where solid ice's is of order a thousandth."""
    body, regime = _regime(SMALL_COLD)
    assert regime == interior_heat.FROZEN
    planet = spin_torques.Planet(body.mass, body.radius * 1000.0)
    rate = 1.0e-4
    shell = icy_tide.k2_over_q(planet.mass_kg, planet.radius_m, rate, 1.0e-3, 60.0)
    frozen = body_spin.dissipation_of(body, planet)(rate)
    mean_k = interior_heat.frozen_mean_temperature_k(body)
    assert 60.0 < mean_k < 273.15
    assert frozen == pytest.approx(tidal_rheology.ice_response(planet.mass_kg, planet.radius_m, rate, mean_k).k2_over_q, rel=1e-3)
    assert shell > 100.0 * frozen


def test_a_water_rich_body_too_warm_for_ice_meets_the_rock_under_its_sea():
    """GEN.1756: a surface over the melting point of ice has no shell, the shell model gave 0 and the planet never despun; it is the rock's rheology (and the air's thermal tide) it meets."""
    assert _regime(WARM_WET)[1] == interior_heat.ROCK
    assert _regime(DRY)[1] == interior_heat.ROCK
