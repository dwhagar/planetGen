# tests/test_life_schedule.py

"""
DB.100: a world keeps the whole schedule it drew (GEN.750: the age it reaches each of the five milestones), has a pace only with a
life class, and an edit that re-works its life keeps the class and the pace it has while the climate still allows them.
"""

import pytest

from planetgen import tuning
from planetgen.generation import evolution
from planetgen.generation import life as planet_life
from planetgen.generation import validation
from planetgen.physics import atmosphere, biogenic_oxygen, life_classes, life_timeline
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy


def _bodies(count, **overrides):
    for seed in range(count):
        config = SystemConfig()
        for key, value in overrides.items():
            setattr(config, key, value)
        with deterministic_entropy(7100 + seed):
            system = StarSystem(system_config=config)
        for planet in system.planets + system.secondary_planets:
            if planet.body_type == "a":
                continue
            yield planet
            yield from getattr(planet, "moons", ())


def _a_world_with_life(**overrides):
    for body in _bodies(60, HABITABLE_WORLD=True, **overrides):
        if body.life_class:
            return body
    pytest.fail("no generated world with a life class")


def _no_draw(monkeypatch):
    def refuse(*_args, **_kwargs):
        raise AssertionError("a draw was made")
    for name in ("random", "uniform", "randrange", "randint", "choice", "choices", "shuffle", "sample", "gauss", "gamma",
                 "dirichlet", "getrandbits"):
        monkeypatch.setattr(draw, name, refuse)


def _ages(body):
    return tuple(getattr(body, field) for field in evolution.SCHEDULE_FIELDS)


def _in_order(ages):
    """Each age that is set is later than the one before it; once one is `None` (the oxygen never came) the rest are."""
    reached = [age for age in ages if age is not None]
    assert ages[:len(reached)] == tuple(reached), "an age is missing before a later one"
    return reached == sorted(reached) and (not reached or reached[0] >= tuning.LIFE_START_AGE_GYR)


def test_the_schedule_fields_are_one_age_for_each_milestone():
    assert evolution.SCHEDULE_FIELDS == tuple(f"life_age_{key}_gyr" for key in tuning.LIFE_MILESTONES)
    for factor in (0.3, 1.0, 2.5):
        for oxygen_delay in (0.0, 2.0, None):
            ages = life_timeline.schedule(life_timeline.unit_waits(), factor, oxygen_delay)
            assert len(ages) == len(evolution.SCHEDULE_FIELDS)
            assert _in_order(ages), "a milestone is never reached before an easier one"


def _short_of_oxygen(body):
    """Whether the world's class breathes oxygen without making any and its air holds too little for complex cells."""
    row = life_classes.describe(body.life_class_id)
    return row.oxidant == "o2" and not row.releases_oxygen and atmosphere.partial_kpa(body, "o2") * 1000.0 < biogenic_oxygen.EUKARYOTE_MINIMUM_PA


def test_a_body_has_a_pace_and_a_schedule_exactly_when_it_has_a_class():
    with_life = without = 0
    for body in _bodies(40, HABITABLE_WORLD=True, MOONS=True):
        stored = _ages(body)
        if body.life_class:
            assert body.evolutionary_pace in tuning.EVOLUTIONARY_TIMELINES, body.name
            assert stored[0] is not None and stored[1] is not None and _in_order(stored), body.name
            if _short_of_oxygen(body):
                assert stored[2:] == (None, None, None), f"{body.name}: complex cells without oxygen"
            with_life += 1
        else:
            assert body.evolutionary_pace is None and stored == (None,) * len(stored), body.name
            without += 1
    assert with_life > 3 and without > 3


def test_a_body_no_class_supports_loses_the_pace_it_had(monkeypatch):
    body = _a_world_with_life()
    assert body.evolutionary_pace is not None
    monkeypatch.setattr(planet_life.life_classes, "score_classes", lambda *args, **kwargs: {})
    planet_life.apply_life_data(body)
    assert body.life_class is None and body.evolutionary_pace is None and body.life_stage is None
    assert all(getattr(body, field) is None for field in evolution.SCHEDULE_FIELDS)


def test_an_edit_keeps_the_class_and_the_pace_while_they_still_fit(monkeypatch):
    body = _a_world_with_life()
    before = (body.life_class_id, body.evolutionary_pace, body.life_stage, body.life_milestone_age_gyr,
              body.life_system_age_gyr, *(getattr(body, field) for field in evolution.SCHEDULE_FIELDS))
    _no_draw(monkeypatch)
    validation.reapply_life(body)
    after = (body.life_class_id, body.evolutionary_pace, body.life_stage, body.life_milestone_age_gyr,
             body.life_system_age_gyr, *(getattr(body, field) for field in evolution.SCHEDULE_FIELDS))
    assert after == before


def test_an_edit_draws_a_class_again_only_when_the_old_one_no_longer_fits(monkeypatch):
    body = _a_world_with_life()
    old = body.life_class_id
    real = life_classes.score_classes

    def without_old(*args, **kwargs):
        found = real(*args, **kwargs)
        found.pop(old, None)
        return found

    monkeypatch.setattr(planet_life.life_classes, "score_classes", without_old)
    validation.reapply_life(body)
    assert body.life_class_id != old
    if body.life_class_id is not None:
        assert body.evolutionary_pace in tuning.EVOLUTIONARY_TIMELINES
        assert _in_order(_ages(body))
        assert body.reflection_spectrum_visible == (life_classes.describe(body.life_class_id).life["reflection_spectrum_visible"] or None)


def test_a_class_that_breathes_oxygen_but_makes_none_stops_at_photosynthesis_in_an_air_without_it():
    for body in _bodies(60, HABITABLE_WORLD=True):
        if body.life_class and _short_of_oxygen(body):
            assert _ages(body)[2:] == (None, None, None), body.name
            assert body.life_stage in (None, "abiogenesis", "photosynthesis"), body.name
            return
    pytest.skip("no generated world with an oxygen-breathing class that makes no oxygen")


def test_the_first_light_comes_earlier_for_a_simpler_way_of_harvesting_it():
    """The first-step scale keeps anoxygenic photosynthesis (Class 2) before the oxygenic kind (Class 1)."""
    assert life_classes.describe("2").life["first_step_scale"] < life_classes.describe("1a").life["first_step_scale"]
