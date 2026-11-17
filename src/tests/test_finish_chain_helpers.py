# tests/test_finish_chain_helpers.py

"""
GEN.1652 slice 1: `planets.end_stages_after_air` and `planets.settle_surroundings` are the blocks that `remake_end_stages`, `_remake_what_reads_the_air`, `refinish_rocky` and
`validation.reconcile_moved_planet` each repeated. They make the same calls in the same order, and read the birth spin from the body's build record, which a body restored from the
database or stored for display answers with the record's default (`record_of`, GEN.1651).
"""

import types

from planetgen.generation.build_record import BuildRecord, HasBuildRecord
from planetgen.physics import daylight, ocean_tides, planets


def test_the_end_stages_after_the_air_run_in_the_first_draw_order_on_the_birth_spin(monkeypatch):
    calls = []
    monkeypatch.setattr(planets, "prepare_class", lambda body: calls.append(("prepare_class", body)))
    monkeypatch.setattr(planets, "generate_spin", lambda body, mass, birth=None: calls.append(("generate_spin", body, mass, birth)))
    monkeypatch.setattr(planets, "classify", lambda body: calls.append(("classify", body)))
    body = types.SimpleNamespace(build=BuildRecord())
    body.build.birth_spin = (1.0e-4, 30.0, None)
    planets.end_stages_after_air(body, 5.0e24)
    assert [c[0] for c in calls] == ["prepare_class", "generate_spin", "classify"]
    assert calls[1][2:] == (5.0e24, (1.0e-4, 30.0, None))


def test_a_body_restored_with_no_record_ends_its_stages_on_the_default_birth_spin(monkeypatch):
    class Restored(HasBuildRecord):
        pass
    seen = []
    monkeypatch.setattr(planets, "prepare_class", lambda body: None)
    monkeypatch.setattr(planets, "generate_spin", lambda body, mass, birth=None: seen.append(birth))
    monkeypatch.setattr(planets, "classify", lambda body: None)
    planets.end_stages_after_air(Restored(), 1.0e24)
    assert seen == [BuildRecord().birth_spin]


def test_the_surroundings_are_settled_daylight_first_then_every_tide_and_the_planets_own_daylight_only_when_asked(monkeypatch):
    calls = []
    monkeypatch.setattr(daylight, "generate", lambda body, planet=None: calls.append(("daylight", body.name, getattr(planet, "name", None))))
    monkeypatch.setattr(ocean_tides, "generate", lambda body, planet=None: calls.append(("tides", body.name, getattr(planet, "name", None))))
    moon = types.SimpleNamespace(name="moon")
    planet = types.SimpleNamespace(name="planet", moons=[moon])
    planets.settle_surroundings(planet, own_daylight=False)
    assert calls == [("daylight", "moon", "planet"), ("tides", "planet", None), ("tides", "moon", "planet")]
    calls.clear()
    planets.settle_surroundings(planet, own_daylight=True)
    assert calls == [("daylight", "planet", None), ("daylight", "moon", "planet"), ("tides", "planet", None), ("tides", "moon", "planet")]
