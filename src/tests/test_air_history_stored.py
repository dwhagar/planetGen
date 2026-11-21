# tests/test_air_history_stored.py

"""GEN.1745 slice 4b: the record of a body's air history is made at the end of the life pass and kept on the body (`generation/life._history_of`, the `air_history` column)."""

import pytest

from planetgen.generation import life
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import air_record
from tests.fuzz_support import deterministic_entropy


def _bodies():
    for seed in range(40):
        config = SystemConfig()
        config.HABITABLE_WORLD = True
        with deterministic_entropy(7500 + seed):
            system = StarSystem(system_config=config)
        for planet in system.primary.orbiters:
            yield planet
            yield from getattr(planet, "moons", ())


def _first(condition):
    for body in _bodies():
        if condition(body):
            return body
    pytest.fail("no generated body of the kind the test needs")


def test_a_world_with_life_keeps_the_record_its_life_pass_made():
    body = _first(lambda b: air_record.eligible(b) and b.life_class is not None and b.build.air_input is not None)
    assert body.air_history is not None and body.air_history.startswith('{"s":')
    record = air_record.unpack(body.air_history)
    assert record.ages and record.steps[-1].age_gyr == pytest.approx(body.star.age, rel=1e-4)
    assert body.air_history == air_record.pack(air_record.record_of(body, body.build.air_living, body.build.air_found, body.build.air_risen, body.build.air_rise_gyr,
                                                                    min(one.started(body) for one in body.build.air_living)))
    assert body.build.air_found and all(amount >= 0.0 for amount in body.build.air_found.values())


def test_every_rocky_body_of_500_pa_or_more_has_a_record_and_a_giant_or_a_thin_air_has_none():
    seen = {"record": False, "none": False}
    for body in _bodies():
        if air_record.eligible(body):
            assert body.air_history is not None
            seen["record"] = True
        else:
            assert body.air_history is None
            seen["none"] = True
        if all(seen.values()):
            return
    pytest.fail(f"the bodies generated did not hold both kinds: {seen}")


def test_an_edit_keeps_the_record_the_body_has():
    body = _first(lambda b: air_record.eligible(b) and b.life_class is not None)
    before = body.air_history
    life.apply_life_data(body, keep_stored=True)
    assert body.air_history == before
