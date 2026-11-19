# tests/test_life_oxygen_after_methane.py

"""GEN.1745 slice 3: the oxygen is asked again at the temperature the life's methane left (`life._build_air`)."""

from types import SimpleNamespace

from planetgen.generation import life
from planetgen.physics import atmosphere, life_air, planet_type


def _run(monkeypatch, oxygen_answers, methane_changes):
    calls = []
    answers = iter(oxygen_answers)

    def oxygen(planet, living):
        calls.append(("oxygen", planet.surface_temperature))
        return next(answers)

    def gases(planet, living, rose):
        calls.append(("gases", rose))
        if methane_changes and not rose:
            planet.surface_temperature += 15.0         # the methane warms it
        return methane_changes if not rose else False

    monkeypatch.setattr(life, "_build_oxygen", oxygen)
    monkeypatch.setattr(life_air, "apply_gases", gases)
    monkeypatch.setattr(life, "_living_of", lambda cands: ["living"])
    monkeypatch.setattr(life, "_rebuild", lambda planet, remake_water=False: calls.append(("rebuild", remake_water)))
    monkeypatch.setattr(planet_type, "has_ground", lambda planet: True)
    monkeypatch.setattr(atmosphere, "partial_kpa", lambda planet, gas: 1.0)
    planet = SimpleNamespace(surface_temperature=255.0, atmospheric_pressure=1.0e5)
    life._build_air(planet, [SimpleNamespace(row=SimpleNamespace(id="1a"))])
    return calls


def test_a_world_the_methane_warmed_is_asked_for_its_oxygen_again_at_the_warmer_temperature(monkeypatch):
    calls = _run(monkeypatch, [False, True], methane_changes=True)
    assert calls[:4] == [("oxygen", 255.0), ("gases", False), ("oxygen", 270.0), ("gases", True)]
    assert calls[-1] == ("rebuild", True)


def test_a_world_the_methane_did_not_change_is_asked_once(monkeypatch):
    calls = _run(monkeypatch, [False], methane_changes=False)
    assert [name for name, _ in calls].count("oxygen") == 1


def test_a_world_whose_oxygen_rose_on_the_first_ask_is_not_asked_again(monkeypatch):
    calls = _run(monkeypatch, [True], methane_changes=True)
    assert [name for name, _ in calls].count("oxygen") == 1
    assert ("gases", True) in calls


def test_a_world_the_warmer_air_still_gives_no_oxygen_keeps_its_methane_air(monkeypatch):
    calls = _run(monkeypatch, [False, False], methane_changes=True)
    assert [name for name, _ in calls].count("oxygen") == 2
    assert ("gases", True) not in calls
