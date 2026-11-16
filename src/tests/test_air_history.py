# tests/test_air_history.py

"""GEN.1745 slice 1: the step schedule and the clocks of a world's air history (`physics/air_history.py`)."""

import pytest

from planetgen.generation import evolution
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import air_history, runaway_history
from tests.fuzz_support import deterministic_entropy


def _worlds(count, **overrides):
    for seed in range(count):
        config = SystemConfig()
        for key, value in overrides.items():
            setattr(config, key, value)
        with deterministic_entropy(7100 + seed):
            system = StarSystem(system_config=config)
        for planet in system.planets + system.secondary_planets:
            if planet.body_type != "a":
                yield planet


def _a_rocky_world_with_life():
    for body in _worlds(60, HABITABLE_WORLD=True):
        if body.body_type == "t" and body.life_class and body.star.age:
            return body
    pytest.fail("no generated rocky world with a life class")


def test_the_schedule_is_logarithmic_then_even_cut_where_asked_and_ends_at_the_star_s_age():
    ages = air_history.schedule(4.6, cuts=(0.7, 2.3, None, 9.0, -1.0))
    assert ages == sorted(set(ages)) and ages[-1] == 4.6
    assert ages[0] == air_history.FIRST_STEP_GYR
    assert 0.7 in ages and 2.3 in ages and 9.0 not in ages
    logarithmic = [age for age in ages if age < air_history.LOG_UNTIL_GYR]
    assert all(later / earlier == pytest.approx(10.0 ** (1.0 / air_history.LOG_STEPS_PER_DECADE), rel=1e-9) for earlier, later in zip(logarithmic, logarithmic[1:]) if later not in (0.7,) and earlier not in (0.7,))
    after = [age for age in ages if age >= air_history.LOG_UNTIL_GYR]
    assert max(later - earlier for earlier, later in zip(after, after[1:])) <= air_history.MAX_STEP_GYR + 1e-12
    assert air_history.schedule(0.5)[-1] == 0.5 and air_history.schedule(0.5)[0] == air_history.FIRST_STEP_GYR
    assert air_history.schedule(None) == [] and air_history.schedule(0.0) == []


def test_a_cut_at_the_star_s_age_or_next_to_a_step_does_not_make_a_second_step():
    ages = air_history.schedule(3.0, cuts=(3.0, 1.25 + 1e-12, 1.0))
    assert ages.count(3.0) == 1 and len(ages) == len(set(ages))
    assert sum(1 for age in ages if abs(age - 1.25) < 1e-6) == 1


def test_the_clocks_of_a_world_run_from_a_dim_hot_young_star_to_its_age_now():
    body = _a_rocky_world_with_life()
    steps = air_history.run(body, cuts=tuple(getattr(body, field) for field in evolution.SCHEDULE_FIELDS))
    assert steps and steps[-1].age_gyr == body.star.age
    last = steps[-1]
    assert last.reductant_ratio == pytest.approx(1.0) and (last.heat_ratio is None or last.heat_ratio == pytest.approx(1.0))
    assert last.luminosity_ratio is None or last.luminosity_ratio == pytest.approx(1.0)
    ages = [step.age_gyr for step in steps]
    assert ages == sorted(ages)
    heat = [step.heat_ratio for step in steps if step.heat_ratio is not None]
    assert heat and heat[0] >= heat[-1], "the interior gives off more heat when it is young"
    lost = [step.water_lost_kg for step in steps]
    assert lost == sorted(lost) and lost[0] >= 0.0
    assert all(step.oxygen_law in (None, *air_history.OXYGEN_LAWS) for step in steps)


def test_the_water_lost_by_the_last_step_is_the_water_the_runaway_history_gives():
    body = _a_rocky_world_with_life()
    steps = air_history.run(body)
    distance = body.build.star_distance_au or body.distance
    assert steps[-1].water_lost_kg == pytest.approx(runaway_history.loss_kg(body.star, distance, body.mass, body.radius), rel=1e-9)


def test_the_trace_of_the_runaway_history_changes_no_answer_and_the_light_of_the_last_age_is_the_light_now():
    body = _a_rocky_world_with_life()
    distance = body.build.star_distance_au or body.distance
    trace = []
    with_trace = runaway_history.loss_kg(body.star, distance, body.mass, body.radius, trace=trace)
    assert with_trace == runaway_history.loss_kg(body.star, distance, body.mass, body.radius)
    assert (not trace and with_trace == 0.0) or trace[-1][1] == pytest.approx(with_trace)
    light = runaway_history.light_at(body.star, body.star.age)
    assert light is None or light == (body.star.luminosity, body.star.temperature)
    early = runaway_history.light_at(body.star, 0.01)
    assert early is None or early[0] > 0.0


def test_a_body_with_no_star_age_has_no_history():
    class Bare:
        star = None
    assert air_history.run(Bare()) == []
