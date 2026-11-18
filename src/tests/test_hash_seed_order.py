# tests/test_hash_seed_order.py

"""
The order of the gases of a blended air is the order they were given, not the order a set gives them (the hash seed of the process orders strings): every sum over the air adds in that order, so a set made
the same sector different in the last digits under two PYTHONHASHSEED values (`test_reproducible_draws`).
"""

from types import SimpleNamespace

from planetgen.physics import nebula_composition as nc
from planetgen.physics import photochemical_loss, surface_air


def test_a_blend_of_two_mixtures_keeps_the_order_the_first_gives_then_the_second_adds():
    first = {"N2": 0.7, "CO2": 0.2, "Ar": 0.1}
    second = {"CH4": 0.5, "CO2": 0.4, "N2": 0.1}
    mixed = nc._mix(first, second, 0.25)
    assert list(mixed) == ["N2", "CO2", "Ar", "CH4"]
    assert mixed["CH4"] == 0.25 * 0.5 and mixed["Ar"] == 0.75 * 0.1


def test_the_air_after_photolysis_keeps_the_order_of_the_air_before_it(monkeypatch):
    partials = {"N2": 8.0e4, "NH3": 5.0e2, "CO2": 1.0e3}
    monkeypatch.setattr(photochemical_loss, "lose", lambda air, *_: {"H2": 1.0, "CO2": 1.0e3, "NH3": 1.0e2, "N2": 8.1e4})
    setting = SimpleNamespace(giant=False, far_uv=1.0, age_gy=1.0, gravity_m_s2=9.8, star_temperature_k=5800.0)
    blended = surface_air._photolysed(partials, setting, surface_air.PHOTOCHEMICAL_MIN_K)
    assert list(blended) == ["N2", "NH3", "CO2", "H2"]
