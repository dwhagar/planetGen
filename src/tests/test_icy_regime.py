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


# --- GEN.1759: the one chooser read by the circularisation, the heat and the sea's tide -----------------------------------------------------------------------------------

STAR_KG = 1.989e30


def _on_orbit(fields, distance_au=3.0, eccentricity=0.3):
    return types.SimpleNamespace(**fields, primary_mass_kg=STAR_KG, distance=distance_au, eccentricity=eccentricity, is_moon=False)


def test_the_choice_does_not_depend_on_the_tide_an_earlier_run_stored():
    """The shell's flux is its own heat plus the tide the shell makes on the orbit, worked out here, so a body stores the same regime whatever `tidal_flux_w_m2` holds."""
    body = _on_orbit(EUROPA, distance_au=5.2, eccentricity=0.01)
    first = interior_heat.icy_regime(body)
    body.tidal_flux_w_m2 = 5.0
    assert interior_heat.icy_regime(body) == first
    body.tidal_flux_w_m2 = 0.0
    assert interior_heat.icy_regime(body) == first


def test_a_body_frozen_through_is_flexed_as_solid_ice_by_the_circularisation_and_the_heat():
    """GEN.1759: `tidal_response` gives the ice sphere's `(k2, Q)` and the tide's power is the Peale power of its `k2/Q`, far under what the shell formula made of the same body."""
    body = _on_orbit(SMALL_COLD, distance_au=3.0, eccentricity=0.5)
    assert interior_heat.icy_regime(body) == interior_heat.FROZEN
    rate = 1.0e-7
    k2, q = interior_heat.tidal_response(body, rate)
    assert k2 / q == pytest.approx(interior_heat.frozen_response(body, rate).k2_over_q)
    forcing = interior_heat.tidal_heating.forcing_of(body)
    shell = interior_heat.tidal_heating.peale_power_w(icy_tide.k2_over_q(body.mass, body.radius * 1000.0, rate, 1.0e-3, 60.0), forcing, body.radius * 1000.0)
    assert interior_heat.frozen_tide_w(body, forcing) < shell / 10.0


def test_a_water_world_too_warm_for_ice_flexes_as_rock_in_the_sea_tide_too():
    from planetgen.physics import ocean_tides, tidal_heating
    body = _on_orbit(WARM_WET, distance_au=0.5, eccentricity=0.0)
    assert ocean_tides._love_k2(body, body.mass, body.radius * 1000.0, 1.0e-4) == tidal_heating.rigid_k2(body.mass, body.radius * 1000.0)
    cold = _on_orbit(EUROPA, distance_au=5.2, eccentricity=0.0)
    assert ocean_tides._love_k2(cold, cold.mass, cold.radius * 1000.0, 1.0e-4) == icy_tide.love_k2(cold.mass, cold.radius * 1000.0, 1.0e-4)


# --- GEN.1760: the orbit of a body too warm for ice is worn by the rock's tide, once its surface is known -----------------------------------------------------------------------


def test_a_surface_under_the_melting_point_of_ice_binds_the_body_to_the_ice_models_and_a_dry_or_warm_one_does_not():
    assert interior_heat.is_ice_bound(types.SimpleNamespace(**EUROPA))
    assert interior_heat.is_ice_bound(types.SimpleNamespace(**SMALL_COLD))
    assert not interior_heat.is_ice_bound(types.SimpleNamespace(**WARM_WET))
    assert not interior_heat.is_ice_bound(types.SimpleNamespace(**DRY))
    unknown = dict(EUROPA, surface_temperature=None)       # circularisation before the air: the cold fallback, a shell
    assert interior_heat.is_ice_bound(types.SimpleNamespace(**unknown))


def test_only_a_warm_water_rich_body_takes_the_tides_of_its_orbit_again(monkeypatch):
    from planetgen.physics import planets
    taken = []
    monkeypatch.setattr(planets, "settle_orbit", lambda body: taken.append(body) or True)
    for fields in (EUROPA, SMALL_COLD, DRY):
        assert not planets.settle_warm_icy_orbit(_on_orbit(fields))
    assert taken == []
    warm = _on_orbit(WARM_WET, distance_au=0.03, eccentricity=0.43)
    assert planets.settle_warm_icy_orbit(warm)
    assert taken == [warm]


def test_the_tide_of_a_warm_water_rich_body_is_in_the_box_not_left_out(monkeypatch):
    """GEN.1759 left it out because its orbit was never worn by the rock's tide; GEN.1760 wears it, so the heat model gives it the mantle's tide like any rocky body (and a body under ice still has none in the box)."""
    forcings = []
    monkeypatch.setattr(interior_heat.radiogenic, "content_of", lambda body: object())
    monkeypatch.setattr(interior_heat, "interior_of", lambda body: object())
    monkeypatch.setattr(interior_heat.radiogenic, "rock_of", lambda body: None)
    monkeypatch.setattr(interior_heat, "lid_of", lambda body: None)
    monkeypatch.setattr(interior_heat, "surface_of", lambda body: None)
    monkeypatch.setattr(interior_heat, "heat_flow", lambda content, interior, mass, radius, age, rock, forcing, lid, surface: forcings.append(forcing))
    interior_heat.heat_at(_on_orbit(WARM_WET, distance_au=0.5, eccentricity=0.2), 4.0)
    interior_heat.heat_at(_on_orbit(EUROPA, distance_au=5.2, eccentricity=0.2), 4.0)
    interior_heat.heat_at(_on_orbit(DRY, distance_au=5.2, eccentricity=0.2), 4.0)
    assert forcings[0] is not None and forcings[0].eccentricity == 0.2
    assert forcings[1] is None
    assert forcings[2] is not None
