# tests/test_store_air_haze.py

"""GEN.1589: `biogenic_oxygen.store_air` takes the haze of the air as the air solver does (a hazy air stored is not warmer by the share of light the haze keeps off the ground), and
`rebalance_co2` finds the CO2 its sink holds with the temperature that CO2 gives, instead of stopping after eight halvings."""

from types import SimpleNamespace

import pytest

from planetgen.physics import aerosols, atmosphere, atmosphere_state, biogenic_oxygen as bo, constants
from planetgen.generation.build_record import BuildRecord

AIR = {"N2": 70000.0, "CO2": 30000.0}


def _body(haze_draws=None, internal=0.0):
    kept = None if haze_draws is None else (haze_draws, 5800.0, None, internal)
    return SimpleNamespace(gravity=1.0, radius=6371.0, surface_temperature=288.0, greenhouse_factor=0.5, mantle_delta_iw=None, _air_haze=kept)


def _no_side_effects(monkeypatch):
    monkeypatch.setattr(atmosphere, "store_air", lambda *args, **kwargs: None)
    monkeypatch.setattr(bo.envelope_share, "fill", lambda body: None)


def test_a_stored_air_takes_the_haze_share_the_air_solver_takes(monkeypatch):
    _no_side_effects(monkeypatch)
    monkeypatch.setattr(aerosols, "haze_of", lambda *args, **kwargs: aerosols.Haze(1.0, 0.3, 500.0))
    body = _body(haze_draws=object())
    gravity = constants.EARTH_GRAVITY
    radiated = constants.STEFAN_BOLTZMANN_CONSTANT * 288.0 ** 4 / 1.5
    hazy = atmosphere_state.greenhouse_warming(AIR, gravity, radiated, 0.3, 500.0)
    clear = atmosphere_state.greenhouse_warming(AIR, gravity, radiated)
    assert abs(hazy.surface_temperature_k - clear.surface_temperature_k) > 1.0         # the haze matters for this air, or the test proves nothing
    bo.store_air(body, AIR)
    assert body.surface_temperature == pytest.approx(hazy.surface_temperature_k, abs=0.05)
    assert body.surface_temperature != pytest.approx(clear.surface_temperature_k, abs=0.5)


def test_a_body_without_haze_draws_or_an_air_step_gets_the_clear_warming(monkeypatch):
    _no_side_effects(monkeypatch)
    gravity = constants.EARTH_GRAVITY
    radiated = constants.STEFAN_BOLTZMANN_CONSTANT * 288.0 ** 4 / 1.5
    clear = atmosphere_state.greenhouse_warming(AIR, gravity, radiated)
    for body in (_body(haze_draws=None), SimpleNamespace(gravity=1.0, radius=6371.0, surface_temperature=288.0, greenhouse_factor=0.5, mantle_delta_iw=None)):
        bo.store_air(body, AIR)
        assert body.surface_temperature == pytest.approx(clear.surface_temperature_k, abs=1e-9)


def test_the_haze_takes_a_share_of_the_sunlight_not_of_the_heat_a_giant_radiates_from_inside(monkeypatch):
    seen = []
    monkeypatch.setattr(aerosols, "haze_of", lambda *args, **kwargs: aerosols.Haze(1.0, 0.4, 500.0))
    monkeypatch.setattr(atmosphere_state, "greenhouse_warming",
                        lambda air, gravity, radiated, aloft=0.0, pressure=None: seen.append((radiated, aloft)) or SimpleNamespace(surface_temperature_k=288.0, factor=0.5))
    body = _body(haze_draws=object(), internal=100.0)
    radiated = bo._warming(body, AIR, constants.EARTH_GRAVITY, 400.0)
    assert seen[-1][1] == pytest.approx(0.4 * (400.0 - 100.0) / 400.0)
    assert radiated.surface_temperature_k == 288.0


def _sink():
    class Sink:
        supply_bar = 1.0
        solved_bar = 0.001
        solved_temperature_k = 300.0

        @staticmethod
        def kept_at(temperature_k, burial=0.0):
            return 0.001 * 2.0 ** ((300.0 - temperature_k) / 10.0)
    return Sink()


def _fake_air(monkeypatch, body, greenhouse_k_per_kpa, base_k=280.0):
    monkeypatch.setattr(atmosphere, "gas_partials_pa", lambda b: dict(b.air))

    def store(b, air):
        b.air = dict(air)
        b.surface_temperature = base_k + greenhouse_k_per_kpa * air.get("CO2", 0.0) / 1000.0
    monkeypatch.setattr(bo, "store_air", store)


def test_the_co2_settles_where_the_sink_holds_it_even_when_the_greenhouse_is_strong(monkeypatch):
    """The 8 passes of halving the gap ended far from the root when the greenhouse is steep (562 K falling to 432 K, still falling): the root is found whatever the slope."""
    for slope in (5.0, 20.0, 60.0, 200.0):
        body = SimpleNamespace(air={"N2": 80000.0, "O2": 20000.0, "CO2": 100000.0}, surface_temperature=300.0, build=BuildRecord(co2_sink=_sink()))
        _fake_air(monkeypatch, body, slope)
        bo.rebalance_co2(body)
        held = max(100000.0 + (_sink().kept_at(body.surface_temperature) - 0.001) * 1.0e5, 0.0)
        assert body.air.get("CO2", 0.0) == pytest.approx(held, rel=0.02, abs=2.0), slope


def test_an_air_that_is_already_settled_is_left_as_it_was_stored(monkeypatch):
    body = SimpleNamespace(air={"N2": 80000.0, "O2": 20000.0, "CO2": 100.0}, surface_temperature=300.0, greenhouse_factor=0.7, build=BuildRecord(co2_sink=_sink()))
    _fake_air(monkeypatch, body, 0.0, base_k=300.0)
    assert bo.rebalance_co2(body) is False
    assert body.air["CO2"] == 100.0 and body.surface_temperature == 300.0 and body.greenhouse_factor == 0.7


def test_a_world_above_the_boiling_point_at_its_pressure_has_no_ocean_for_the_sink_and_is_left_as_it_is(monkeypatch):
    """GEN.1821: an Earth copy of 442 K under 615 kPa (its water boils above 431 K at that pressure) has no liquid surface, so the weathering of its sink is not asked and no trial air is stored."""
    class Cold:
        solved_bar = 0.001

        @staticmethod
        def kept_at(temperature_k, burial=0.0):
            return 1.0e-12
    body = SimpleNamespace(air={"N2": 180000.0, "H2O": 435000.0, "CO2": 25.0}, surface_temperature=442.0, atmospheric_pressure=615000.0, build=BuildRecord(co2_sink=Cold()))
    monkeypatch.setattr(atmosphere, "gas_partials_pa", lambda b: dict(b.air))
    monkeypatch.setattr(bo, "store_air", lambda *a, **k: pytest.fail("a steaming world is not asked about its sink"))
    assert bo.rebalance_co2(body, 0.02) is False and body.surface_temperature == 442.0


def test_a_trial_air_past_the_critical_point_of_water_is_not_a_root_either_way(monkeypatch):
    """GEN.1821: where a trial air with less CO2 than the world began with runs away past water's critical point (647 K) the root is above it; more CO2 than it began with that does is too much for the cycle (the root is below)."""
    class Cold:
        solved_bar = 0.001

        @staticmethod
        def kept_at(temperature_k, burial=0.0):
            return 1.0e-12
    body = SimpleNamespace(air={"N2": 100000.0, "CO2": 25.0}, surface_temperature=300.0, atmospheric_pressure=100000.0, build=BuildRecord(co2_sink=Cold()))
    monkeypatch.setattr(atmosphere, "gas_partials_pa", lambda b: dict(b.air))

    def store(b, air):
        b.air = dict(air)
        co2 = air.get("CO2", 0.0)
        b.surface_temperature = 5000.0 if co2 < 3.0 or co2 > 200.0 else 300.0 + 2.0 * (co2 - 25.0) if co2 > 25.0 else 300.0
    monkeypatch.setattr(bo, "store_air", store)
    bo.rebalance_co2(body, 0.02)
    assert body.surface_temperature < 647.0 and 3.0 <= body.air["CO2"] <= 200.0
