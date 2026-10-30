# tests/test_habitable_request.py

"""GEN.911, GEN.916: a request for a habitable world asks for conditions (temperate, air, a surface ocean), never for a class letter, and
aims its orbit by the temperature the draw solved."""

import types

import pytest

from planetgen import tuning
from planetgen.generation import habitable_request as hr
from planetgen.generation import system as systemData
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw


def _body(temperature=290.0, pressure=101325.0, hydrosphere="surface ocean", planet_class="C", moons=(), body_type="t"):
    return types.SimpleNamespace(surface_temperature=temperature, atmospheric_pressure=pressure, hydrosphere=hydrosphere,
                                 planet_class=planet_class, moons=list(moons), body_type=body_type)


def test_a_temperate_world_with_air_and_a_sea_meets_the_request_whatever_its_class_letter():
    for letter in ("A", "C", "M", "Z"):
        assert hr.meets(_body(planet_class=letter))


@pytest.mark.parametrize("change", [
    {"temperature": 250.0}, {"temperature": 340.0}, {"temperature": None},
    {"pressure": 5000.0}, {"pressure": None}, {"hydrosphere": "ice-covered ocean"}, {"hydrosphere": "hycean"},
    {"hydrosphere": None}, {"body_type": "g"},
])
def test_a_world_missing_any_condition_does_not_meet_it_even_in_a_class_that_can_carry_life(change):
    assert tuning.HABITABLE_PLANET_CLASSES[0] and not hr.meets(_body(planet_class="M", **change))


def test_the_edges_of_the_range_and_the_pressure_floor_count():
    low, high = tuning.REQUEST_TEMPERATURE_K
    assert hr.meets(_body(temperature=low)) and hr.meets(_body(temperature=high))
    assert hr.meets(_body(pressure=tuning.REQUEST_MIN_PRESSURE_PA))
    assert not hr.meets(_body(pressure=tuning.REQUEST_MIN_PRESSURE_PA - 1.0))


def test_a_moon_counts_and_a_belt_is_skipped():
    moon = _body()
    planet = _body(temperature=200.0, moons=[moon])
    belt = types.SimpleNamespace(body_type="a", moons=[])
    assert hr.meets_with_moons(planet) and not hr.meets(planet)
    assert hr.count([planet, belt, _body()]) == 2


def test_the_class_letter_table_does_not_decide_the_count(monkeypatch):
    monkeypatch.setattr(tuning, "HABITABLE_PLANET_CLASSES", [])
    assert hr.count([_body(planet_class="M")]) == 1
    monkeypatch.setattr(tuning, "HABITABLE_PLANET_CLASSES", ["M"])
    assert hr.count([_body(temperature=200.0, planet_class="M")]) == 0


def test_a_cold_draw_aims_closer_a_hot_one_farther_and_the_zone_bounds_it():
    target = tuning.REQUEST_TARGET_K
    assert hr.aimed_distance(2.0, target, 0.5, 5.0) == pytest.approx(2.0)
    assert hr.aimed_distance(2.0, target / 2.0, 0.5, 5.0) == pytest.approx(0.5)      # a quarter of the distance, kept to the zone's inner edge
    assert hr.aimed_distance(2.0, target * 0.9, 0.5, 5.0) == pytest.approx(2.0 * 0.81)
    assert hr.aimed_distance(2.0, target * 1.2, 0.5, 5.0) == pytest.approx(2.0 * 1.44)
    assert hr.aimed_distance(2.0, target * 3.0, 0.5, 5.0) == 5.0
    assert hr.aimed_distance(2.0, None, 0.5, 5.0) == 2.0 and hr.aimed_distance(2.0, 0.0, 0.5, 5.0) == 2.0


@pytest.fixture(scope="module")
def system():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.BINARY_SYSTEM = False
    draw.set_run_seed(5)
    return StarSystem(system_config=cfg)


def test_a_habitable_request_leaves_its_first_draws_where_they_were_and_then_aims_by_the_solved_temperature(system, monkeypatch):
    seen = []

    class Recorder:
        planet_class = "M"                      # a class that can carry life is no success without the conditions
        moons = []
        body_type = "t"
        surface_temperature = 0.75 * tuning.REQUEST_TARGET_K
        atmospheric_pressure = 0.0
        hydrosphere = None

        def __init__(self, config, star, zone, distance, **kwargs):
            seen.append(distance)

        def add_moons(self, count):
            pass

    monkeypatch.setattr(systemData, "Planet", Recorder)
    inner, outer = system.star.habitable_zone
    start = (inner + outer) / 2.0
    system.build_planet(system.star, system.star.habitable_zone, start, life=True)
    after = tuning.HABITABLE_REQUEST_AIM_AFTER_ATTEMPTS
    assert seen[:after] == [start] * after
    assert seen[after] == pytest.approx(max(start * 0.5625, max(inner, system._orbit_floor_au(system.star))))
    assert all(inner * 0.999 <= d <= outer for d in seen[after:]) or system._orbit_floor_au(system.star) > inner
    assert system._request_unreachable


def test_a_request_for_a_letter_alone_does_not_move_the_orbit(system, monkeypatch):
    seen = []

    class Recorder:
        planet_class = "C"
        moons = []
        body_type = "t"
        surface_temperature = 100.0

        def __init__(self, config, star, zone, distance, **kwargs):
            seen.append(distance)

        def add_moons(self, count):
            pass

    monkeypatch.setattr(systemData, "Planet", Recorder)
    system.build_planet(system.star, system.star.habitable_zone, 1.0, planet_class="D")
    assert set(seen) == {1.0}


def test_a_requested_habitable_world_is_judged_by_its_conditions_in_the_whole_system():
    found = 0
    for seed in range(6):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "K5V"
        cfg.BINARY_SYSTEM = False
        cfg.HABITABLE_WORLD = True
        draw.set_run_seed(seed)
        built = StarSystem(system_config=cfg)
        met = hr.count(built.planets) > 0
        assert met == ("a habitable world" not in built.unmet_requirements)
        found += met
    assert found >= 1
