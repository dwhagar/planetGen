# tests/test_induction_cores.py

"""
GEN.981: the induced magnetic moment of a metal core and of a rocky planet in its star's field (`magnetism.sphere_response`, `core_induced_moment`,
`generate_induced_field`, `generate_stellar_induction`).
"""

import math
import types

import pytest

from planetgen.physics import core_heat, magnetism


def _core(radius_m=3.48e6):
    return core_heat.Core(temperature_k=5000.0, boundary_k=4000.0, cmb_w=1e13, adiabatic_w=1e12, inner_fraction=0.3, dynamo_w=1e12,
                          dynamo_ended_gy=None, radius_m=radius_m, density_kg_m3=11000.0)


def test_a_sphere_much_thinner_than_a_skin_depth_answers_as_omega_mu_sigma_a2_over_15():
    omega, sigma, radius = 1e-3, 10.0, 1e3
    assert magnetism.sphere_response(omega, sigma, radius) == pytest.approx(omega * magnetism.VACUUM_PERMEABILITY * sigma * radius ** 2 / 15.0, rel=1e-3)


def test_a_sphere_many_skin_depths_across_answers_with_one():
    assert magnetism.sphere_response(7.3e-5, magnetism.CORE_CONDUCTIVITY_S_M, 3.48e6) == pytest.approx(1.0, abs=1e-3)
    assert magnetism.sphere_response(1e-9, magnetism.CORE_CONDUCTIVITY_S_M, 1e9) == 1.0


def test_the_response_grows_with_frequency_conductivity_and_size_and_never_passes_one():
    steps = [magnetism.sphere_response(1e-9 * 10 ** k, 1e6, 1e5) for k in range(0, 8)]
    assert steps == sorted(steps) and steps[-1] <= 1.0 and steps[0] < 0.6 * steps[-1]
    assert magnetism.sphere_response(1e-6, 1e6, 1e5) < magnetism.sphere_response(1e-6, 1e7, 1e5) < magnetism.sphere_response(1e-6, 1e7, 1e6)


def test_no_frequency_conductor_or_size_gives_no_response():
    assert magnetism.sphere_response(0.0, 1e6, 1e6) == 0.0
    assert magnetism.sphere_response(1e-4, 0.0, 1e6) == 0.0
    assert magnetism.sphere_response(1e-4, 1e6, 0.0) == 0.0
    assert magnetism.sphere_response(None, 1e6, 1e6) == 0.0


def test_a_perfect_conductor_core_has_the_moment_of_a_superconducting_sphere():
    excitation, period = 2e-7, 4e4
    moment = magnetism.core_induced_moment(_core(6e5), excitation, period)
    assert moment == pytest.approx(2.0 * math.pi * (6e5) ** 3 * excitation / magnetism.VACUUM_PERMEABILITY, rel=1e-3)
    assert magnetism.core_induced_moment(_core(6e5), excitation, period, transmission=0.5) == pytest.approx(0.5 * moment, rel=1e-3)
    assert magnetism.core_induced_moment(None, excitation, period) == 0.0
    assert magnetism.core_induced_moment(_core(), 0.0, period) == 0.0
    assert magnetism.core_induced_moment(_core(), excitation, None) == 0.0


def test_a_moons_core_adds_to_its_oceans_moment(monkeypatch):
    from planetgen.physics import interior_heat
    parent = types.SimpleNamespace(magnetic_moment_a_m2=1.56e27, rotation_period_hours=9.925, dipole_class="strong", body_type="g", is_moon=False)
    europa = types.SimpleNamespace(hydrosphere="ice-covered ocean", ocean_depth_km=100.0, ice_shell_km=20.0, ocean_class="neutral", body_type="t",
                                   period=85.228 / 24.0 / 365.25, distance=6.711e8 / 1.496e11, radius=1560.8)
    monkeypatch.setattr(interior_heat, "core_of", lambda body: None)
    magnetism.generate_induced_field(europa, parent)
    ocean_only = europa.induced_moment_a_m2
    monkeypatch.setattr(interior_heat, "core_of", lambda body: _core(6e5))
    magnetism.generate_induced_field(europa, parent)
    assert europa.induced_moment_a_m2 > ocean_only
    assert europa.induced_moment_a_m2 < 2.0 * ocean_only       # the ocean screens the core: it adds a fraction
    europa.hydrosphere = "dry"                                  # a dry moon: the core alone answers
    magnetism.generate_induced_field(europa, parent)
    assert 0 < europa.induced_moment_a_m2 < ocean_only


def _planet(hydrosphere="dry", period_hours=24.0):
    star = types.SimpleNamespace(rotation_period_hours=600.0, pulsar_type=None, radius=None, luminosity=None, log_lx_lbol=None)
    return types.SimpleNamespace(body_type="t", is_moon=False, star=star, _star_distance_au=1.0, period=1.0, rotation_period_hours=period_hours,
                                 radius=6371.0, hydrosphere=hydrosphere, ocean_depth_km=3.7, ice_shell_km=0.0, ocean_class="neutral",
                                 induced_moment_a_m2=None)


def test_a_rocky_planet_in_its_stars_field_gets_the_induced_moment_of_its_core(monkeypatch):
    from planetgen.physics import interior_heat
    monkeypatch.setattr(interior_heat, "core_of", lambda body: _core())
    planet = _planet()
    magnetism.generate_stellar_induction(planet)
    expected = 2.0 * math.pi * 3.48e6 ** 3 * magnetism.STELLAR_EXCITATION_SHARE * magnetism.IMF_EARTH_T / magnetism.VACUUM_PERMEABILITY
    assert planet.induced_moment_a_m2 == pytest.approx(expected, rel=0.01)
    assert planet.induced_moment_a_m2 < 1e-3 * magnetism.EARTH_MOMENT_A_M2    # the answer to the star's field is small beside a dynamo


def test_a_stronger_wind_gives_a_larger_moment(monkeypatch):
    from planetgen.physics import interior_heat
    monkeypatch.setattr(interior_heat, "core_of", lambda body: _core())
    near, far = _planet(), _planet()
    near._star_distance_au, far._star_distance_au = 0.1, 2.0
    for planet in (near, far):
        magnetism.generate_stellar_induction(planet)
    assert near.induced_moment_a_m2 == pytest.approx(far.induced_moment_a_m2 * (2.0 / 0.1), rel=1e-6)      # B goes as the root of the pressure: 1 / r


def test_a_tidally_locked_planet_feels_the_stars_rotation_and_a_giant_or_a_coreless_dry_one_none(monkeypatch):
    from planetgen.physics import interior_heat
    monkeypatch.setattr(interior_heat, "core_of", lambda body: _core())
    locked = _planet(period_hours=365.25 * 24.0)
    magnetism.generate_stellar_induction(locked)
    assert locked.induced_moment_a_m2 > 0
    giant = _planet()
    giant.body_type = "g"
    magnetism.generate_stellar_induction(giant)
    assert giant.induced_moment_a_m2 is None
    monkeypatch.setattr(interior_heat, "core_of", lambda body: None)
    coreless = _planet()
    magnetism.generate_stellar_induction(coreless)
    assert coreless.induced_moment_a_m2 == 0.0


def test_a_surface_ocean_adds_its_shell_to_the_planets_moment(monkeypatch):
    from planetgen.physics import interior_heat
    monkeypatch.setattr(interior_heat, "core_of", lambda body: None)
    dry, wet = _planet(), _planet(hydrosphere="surface ocean")
    for planet in (dry, wet):
        magnetism.generate_stellar_induction(planet)
    assert dry.induced_moment_a_m2 == 0.0
    assert 0 < wet.induced_moment_a_m2 < 1e-3 * magnetism.EARTH_MOMENT_A_M2
    wet.ocean_fraction = 0.1                 # an ocean that covers a tenth of the surface is a smaller shell
    magnetism.generate_stellar_induction(wet)
    assert 0 < wet.induced_moment_a_m2 < 0.6 * dry_to_full(wet)


def dry_to_full(planet):
    """The moment of the same planet with an ocean over all its surface."""
    full = types.SimpleNamespace(**vars(planet))
    full.ocean_fraction = 1.0
    magnetism.generate_stellar_induction(full)
    return full.induced_moment_a_m2
