# tests/test_air_history.py

"""GEN.1745 slice 1: the step schedule and the clocks of a world's air history (`physics/air_history.py`)."""

import copy
import functools

import pytest

from planetgen.generation import evolution
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import air_history, atmosphere, runaway_history
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


@functools.lru_cache(maxsize=None)
def _a_world_with_its_air_input():
    """`(body, (setting, start_k, entry_state))` of a generated rocky world with life: the input of its final air solve, caught as it was made."""
    from planetgen.physics import surface_air
    from planetgen.util import draw
    caught = []
    real = surface_air.solve

    def spy(setting, start):
        state = draw.getstate()
        air = real(setting, start)
        caught.append((setting, start, state, air))
        return air

    surface_air.solve = spy
    try:
        for body in _worlds(60, HABITABLE_WORLD=True):
            if body.body_type != "t" or not body.life_class or not body.star.age:
                continue
            for setting, start, state, air in reversed(caught):
                if abs(setting.radius_m - body.radius * 1000.0) < 1e-6 * setting.radius_m and air_history.aged_setting(setting, body, body.star.age) is not None:
                    return body, (setting, start, state)
    finally:
        surface_air.solve = real
    pytest.fail("no generated rocky world with a life class and an air solve")


def test_the_aged_setting_of_the_star_s_age_is_the_setting_of_the_final_air():
    body, (setting, _, _) = _a_world_with_its_air_input()
    now = air_history.aged_setting(setting, body, body.star.age)
    assert now.age_gy == body.star.age
    assert now.starlight_w_m2 == pytest.approx(setting.starlight_w_m2, rel=1e-9) and now.absorbed_flux_w_m2 == pytest.approx(setting.absorbed_flux_w_m2, rel=1e-9)
    assert now.heat_flux_w_m2 == pytest.approx(setting.heat_flux_w_m2, rel=1e-6)
    young = air_history.aged_setting(setting, body, 0.2)
    assert young.heat_flux_w_m2 > now.heat_flux_w_m2 and young.released_j < now.released_j
    assert young.water_kg == setting.water_kg and young.mantle_kg == setting.mantle_kg


def test_the_anchors_are_the_first_and_the_last_step_and_at_most_the_cap_and_move_where_light_or_heat_move():
    scale = [1.0 + 0.0005 * k for k in range(40)]
    steps = list(range(40))
    assert air_history._anchor_indices(steps, scale, [1.0] * 40) == [0, 39]
    moving = [1.0 * 1.2 ** k for k in range(40)]
    chosen = air_history._anchor_indices(steps, moving, [1.0] * 40)
    assert chosen[0] == 0 and chosen[-1] == 39 and chosen == sorted(set(chosen)) and len(chosen) <= air_history.MAX_ANCHORS
    heat = [1.0 * 1.5 ** k for k in range(40)]
    assert len(air_history._anchor_indices(steps, [1.0] * 40, heat)) > 2


def test_a_temperature_that_follows_the_light_is_read_exactly_between_anchors_and_solved_only_at_them(monkeypatch):
    body, air_input = _a_world_with_its_air_input()
    setting = air_input[0]
    solved = []

    def fake(aged, start_k, entry_state):
        solved.append(aged.age_gy)
        return 200.0 + 0.5 * aged.absorbed_flux_w_m2

    monkeypatch.setattr(air_history, "_solved_temperature", fake)
    steps = air_history.run(body, air_input=air_input)
    assert steps and all(step.temperature_k is not None for step in steps)
    assert solved and len(solved) == sum(1 for step in steps if step.anchor) <= air_history.MAX_ANCHORS
    assert steps[0].anchor and steps[-1].anchor and steps[-1].age_gyr == body.star.age
    for step in steps:
        aged = air_history.aged_setting(setting, body, step.age_gyr)
        if step.anchor:
            assert step.temperature_k == pytest.approx(200.0 + 0.5 * aged.absorbed_flux_w_m2)
        else:
            assert step.temperature_k == pytest.approx(200.0 + 0.5 * aged.absorbed_flux_w_m2, rel=0.02)


def test_a_history_asked_for_no_air_has_no_temperatures_and_an_anchor_with_no_air_leaves_none(monkeypatch):
    body, air_input = _a_world_with_its_air_input()
    assert all(step.temperature_k is None and not step.anchor for step in air_history.run(body))
    monkeypatch.setattr(air_history, "_solved_temperature", lambda *args: None)
    assert all(step.temperature_k is None for step in air_history.run(body, air_input=air_input))


def test_an_anchor_is_the_air_the_body_has_now_solved_under_the_light_and_heat_of_its_age():
    body, air_input = _a_world_with_its_air_input()
    setting, start_k, state = air_input
    from planetgen.physics import surface_air
    from planetgen.util import draw
    stream = draw.Stream(0)
    stream.setstate(state)
    with draw.bound(stream):
        final = surface_air.solve(setting, start_k)
    now = air_history._solved_temperature(air_history.aged_setting(setting, body, body.star.age), start_k, state)
    assert final is not None and now == pytest.approx(final.temperature_k, abs=0.5)
    assert air_history._solved_temperature(air_history.aged_setting(setting, body, 0.1), start_k, state) != now


@functools.lru_cache(maxsize=None)
def _a_world_the_oxygen_rises_on():
    """`(body, living, rise_gyr)` of a generated rocky world on which class 1a producers of 0.3 Gyr raise the oxygen before the star's age."""
    from planetgen.physics import biogenic_oxygen, life_air, life_classes
    living = [life_air.Living(life_classes.BY_ID["1a"], 1.0, 0.3)]
    for body in _worlds(60, HABITABLE_WORLD=True):
        if body.body_type != "t" or not body.star.age or not getattr(body, "radius", None):
            continue
        rise, _ = life_air.oxygen_state(body, living)
        if rise is None:
            continue
        have = biogenic_oxygen.air_moles(body).get("O2", 0.0)
        ages = air_history.schedule(body.star.age, air_history.oxygen_cuts(body, living))
        if air_history.oxygen_trajectory(body, living, ages)[-1][1] > have:
            return body, living, rise
        # a rise date with nothing to raise: the air already holds more oxygen than the balance gives (the oxygen of its lost water, 6.6e19 mol on the only world of 60 systems that has a rise
        # since GEN.1732), and the trajectory rightly stays at what the air has, since it never takes oxygen out; the same world with no abiotic oxygen is the world the rise acts on
        bare = copy.copy(body)
        setattr(bare, atmosphere.GAS_FIELD, {formula: pascal for formula, pascal in atmosphere.gas_partials_pa(body).items() if formula != "O2"})
        return bare, living, rise
    pytest.fail("no generated rocky world on which the oxygen of class 1a rises")


def test_the_oxygen_trajectory_is_the_air_as_found_until_the_rise_and_never_falls_after_it():
    """GEN.1745 slice 2a: before the rise the oxygen the air has (the abiotic oxygen, none for most), from the rise the plateau then the balance as the reductant sink falls with the interior's heat."""
    from planetgen.physics import biogenic_oxygen
    body, living, rise = _a_world_the_oxygen_rises_on()
    cuts = air_history.oxygen_cuts(body, living)
    assert rise in cuts and 0.3 in cuts
    ages = air_history.schedule(body.star.age, cuts)
    trajectory = air_history.oxygen_trajectory(body, living, ages)
    have = biogenic_oxygen.air_moles(body).get("O2", 0.0)
    assert [age for age, _ in trajectory] == ages
    assert all(value == have for age, value in trajectory if age < rise)
    after = [value for age, value in trajectory if age >= rise]
    assert after and all(b >= a - 1e-9 * max(a, 1.0) for a, b in zip(after, after[1:]))        # the reductant only falls, so the oxygen never does
    assert after[-1] > have                                                                       # the rise put oxygen into the air


def test_the_oxygen_trajectory_ends_at_the_oxygen_the_air_is_given_at_the_star_s_age():
    """The last value is what `stage_target_mol` gives, held to the fire limit and not under what the burnt air holds (what `biogenic_oxygen.apply` puts in); a world with no rise stays as found."""
    from planetgen.physics import biogenic_oxygen, life_air
    body, living, rise = _a_world_the_oxygen_rises_on()
    ages = air_history.schedule(body.star.age, air_history.oxygen_cuts(body, living))
    last = air_history.oxygen_trajectory(body, living, ages)[-1][1]
    rise_again, target = life_air.oxygen_state(body, living)
    burnt, _ = biogenic_oxygen.burn(biogenic_oxygen.air_moles(body))
    assert rise_again == rise
    assert last == pytest.approx(max(burnt.get("O2", 0.0), min(target, biogenic_oxygen._fire_limit_mol(body, burnt))), rel=1e-9)
    none = [life_air.Living(living[0].life_class, 1.0, 1.0e3)]            # producers older than the star can be: no rise
    assert air_history.oxygen_cuts(body, none) == []
    flat = air_history.oxygen_trajectory(body, none, ages)
    assert len({value for _, value in flat}) == 1
