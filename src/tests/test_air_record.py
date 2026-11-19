# tests/test_air_record.py

"""GEN.1745 slice 4a: the record a body keeps of its air history (`physics/air_record.py`)."""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import air_history, air_record, atmosphere
from tests.fuzz_support import deterministic_entropy


def _a_world_with_air():
    for seed in range(60):
        config = SystemConfig()
        config.HABITABLE_WORLD = True
        with deterministic_entropy(7400 + seed):
            system = StarSystem(system_config=config)
        for planet in system.primary.orbiters:
            if air_record.eligible(planet) and planet.star.age > 2.0 and planet.build.air_input is not None:
                return planet
    pytest.fail("no generated rocky world with an air of 500 Pa or more")


def test_a_world_with_an_air_gets_a_record_whose_last_step_is_the_world_s_own_air_and_age():
    body = _a_world_with_air()
    record = air_record.record_of(body)
    assert record is not None and 2 <= len(record.steps) <= air_record.MAX_STEPS
    ages = [step.age_gyr for step in record.steps]
    assert ages == sorted(ages) and ages[-1] == pytest.approx(body.star.age)
    last = record.steps[-1]
    stored = atmosphere.gas_partials_pa(body)
    assert set(last.partials_pa) == {formula for formula, pascal in stored.items() if pascal >= air_record.PARTIAL_FLOOR_PA}
    assert last.partials_pa == {formula: stored[formula] for formula in last.partials_pa}      # no life: the stored air at every age
    assert all(step.partials_pa == last.partials_pa for step in record.steps)
    assert all(step.temperature_k is None or step.temperature_k > 0.0 for step in record.steps)


def test_the_last_temperature_is_the_final_solve_s_own_to_a_kelvin_when_life_did_not_move_it():
    body = _a_world_with_air()
    record = air_record.record_of(body)
    assert record.steps[-1].temperature_k == pytest.approx(body.surface_temperature, abs=1.0)


def test_a_record_packs_to_a_short_text_and_back_to_the_digits_it_kept():
    body = _a_world_with_air()
    record = air_record.record_of(body)
    text = air_record.pack(record)
    assert len(text) < 8000
    again = air_record.unpack(text)
    assert len(again.steps) == len(record.steps)
    for one, other in zip(record.steps, again.steps):
        assert other.age_gyr == pytest.approx(one.age_gyr, rel=1e-4) and other.event == one.event and other.oxygen_law == one.oxygen_law
        assert set(other.partials_pa) == set(one.partials_pa)
        if one.temperature_k is not None:
            assert other.temperature_k == pytest.approx(one.temperature_k, rel=1e-4)


def test_a_thin_air_a_giant_or_a_world_with_no_star_age_gets_no_record():
    body = _a_world_with_air()
    assert air_record.record_of(body) is not None
    saved = body.atmospheric_pressure
    body.atmospheric_pressure = air_record.MIN_AIR_PA - 1.0
    try:
        assert air_record.record_of(body) is None
    finally:
        body.atmospheric_pressure = saved
    body.body_type = "g"
    try:
        assert air_record.record_of(body) is None
    finally:
        body.body_type = "t"


def test_the_steps_are_thinned_to_the_limit_keeping_the_first_the_last_and_the_special_ones():
    steps = [air_history.Step(age_gyr=0.01 * (index + 1), luminosity_ratio=1.0, heat_w_m2=None, heat_ratio=None, reductant_ratio=1.0, water_lost_kg=0.0, oxygen_law=None) for index in range(150)]
    special = {steps[40].age_gyr, steps[90].age_gyr}
    kept = air_record._thinned(steps, special)
    assert len(kept) <= air_record.MAX_STEPS + 1
    assert kept[0] is steps[0] and kept[-1] is steps[-1] and steps[40] in kept and steps[90] in kept
    assert [step.age_gyr for step in kept] == sorted(step.age_gyr for step in kept)
