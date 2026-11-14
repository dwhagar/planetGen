# tests/test_air_final.py

"""
GEN.1657: the air a body ends with (`physics/air_final.py`) is one value, whether the solver found an air, an oddity replaced it or the body kept none.

The airless cases are driven with the solver stubbed (a moon, so no oddity rolls): the body keeps no pressure, density, scale height or partials, its temperature is the starlight's alone,
or the ice-albedo law's own answer when the setting has one.
"""

import types

from planetgen.physics import air_branches, air_final, surface_air


def _moon():
    return types.SimpleNamespace(is_moon=True, body_type="t", mass=None, star=types.SimpleNamespace(temperature=None))


def test_a_body_with_no_air_keeps_the_temperature_starlight_gives(monkeypatch):
    monkeypatch.setattr(surface_air, "solve", lambda setting, start_k: None)
    final = air_final.make(_moon(), types.SimpleNamespace(albedo_law=None), 255.0, 1361.0, 1.0, 0.3, False, None)
    assert (final.air, final.temperature_k, final.albedo, final.greenhouse_factor) == (None, 255.0, 0.3, 0.0)
    assert final.pressure_pa == 0.0 and final.partials_pa == {}
    assert final.molar_density_kg_mol is final.density_kg_m3 is final.scale_height_km is None


def test_a_body_with_no_air_and_an_ice_law_keeps_the_laws_temperature_and_albedo(monkeypatch):
    monkeypatch.setattr(surface_air, "solve", lambda setting, start_k: None)
    seen = []
    monkeypatch.setattr(surface_air, "airless_state", lambda law, starlight, start: seen.append((law, starlight, start)) or (123.0, 0.55))
    final = air_final.make(_moon(), types.SimpleNamespace(albedo_law="law"), 255.0, 1361.0, 1.0, 0.3, False, None)
    assert seen == [("law", 1361.0, 255.0)]
    assert (final.temperature_k, final.albedo, final.greenhouse_factor, final.pressure_pa) == (123.0, 0.55, 0.0, 0.0)


def test_a_body_past_its_runaway_limit_with_water_starts_the_solver_hot(monkeypatch):
    starts = []
    monkeypatch.setattr(surface_air, "solve", lambda setting, start_k: starts.append(start_k))
    monkeypatch.setattr(air_final, "past_runaway_with_water", lambda planet, starlight, water: water is not None)
    air_final.make(_moon(), types.SimpleNamespace(albedo_law=None), 255.0, 1361.0, 1.0, 0.3, False, 1e-3)
    air_final.make(_moon(), types.SimpleNamespace(albedo_law=None), 255.0, 1361.0, 1.0, 0.3, False, None)
    assert starts[0] > 255.0 and starts[1] == 255.0


def test_a_body_with_the_ice_law_enters_the_solver_from_the_warm_side_and_a_giant_or_a_body_without_the_law_does_not(monkeypatch):
    """GEN.1696: a planet cooled from its formation, so a world that can be warm or frozen is solved from the warm side."""
    starts = []
    monkeypatch.setattr(surface_air, "solve", lambda setting, start_k: starts.append(start_k))
    monkeypatch.setattr(surface_air, "airless_state", lambda law, starlight, start: (start, 0.3))
    air_final.make(_moon(), types.SimpleNamespace(albedo_law="law"), 255.0, 1361.0, 1.0, 0.3, False, 1e-3)
    air_final.make(_moon(), types.SimpleNamespace(albedo_law="law"), 255.0, 1361.0, 1.0, 0.3, True, 1e-3)
    air_final.make(_moon(), types.SimpleNamespace(albedo_law=None), 255.0, 1361.0, 1.0, 0.3, False, 1e-3)
    assert starts == [air_branches.WARM_START_K, 255.0, 255.0]
    air_final.make(_moon(), types.SimpleNamespace(albedo_law="law"), 400.0, 1361.0, 1.0, 0.3, False, 1e-3)
    assert starts[-1] == 400.0      # a body already hotter than the warm side keeps its own start
