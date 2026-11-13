# tests/test_life_remakes_water.py

"""GEN.1590: when life changes the air, the water is made again at the temperature and pressure the air left, and a thaw remakes the class after the ground."""

import types

from planetgen.generation import life
from planetgen.physics import planets


def _body(temperature, pressure):
    return types.SimpleNamespace(surface_temperature=temperature, atmospheric_pressure=pressure)


def test_the_water_is_made_again_when_the_air_moved_the_temperature(monkeypatch):
    calls = []
    monkeypatch.setattr(planets, "refresh_at_final_age", lambda body: calls.append(body))
    body = _body(432.0, 63000.0)
    life._remake_water_if_moved(body, (293.0, 63000.0))
    assert calls == [body]


def test_the_water_is_made_again_when_the_air_moved_the_pressure(monkeypatch):
    calls = []
    monkeypatch.setattr(planets, "refresh_at_final_age", lambda body: calls.append(body))
    life._remake_water_if_moved(_body(293.0, 10000.0), (293.0, 63000.0))
    assert len(calls) == 1


def test_a_body_whose_temperature_and_pressure_did_not_move_keeps_its_water(monkeypatch):
    calls = []
    monkeypatch.setattr(planets, "refresh_at_final_age", lambda body: calls.append(body))
    life._remake_water_if_moved(_body(293.0, 63000.0), (293.0, 63000.0))
    assert calls == []
