# tests/test_life_schedule.py

"""
DB.100: a world keeps the whole schedule of its pace (the age it reaches each of the five milestones), has a pace only with a
life chemistry, and an edit that re-works its life keeps the chemistry and the pace it has while the climate still allows them.
"""

import pytest

from planetgen import tuning
from planetgen.generation import evolution
from planetgen.generation import life as planet_life
from planetgen.generation import validation
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
        if body.life_chemical:
            return body
    pytest.fail("no generated world with a life chemistry")


def _no_draw(monkeypatch):
    def refuse(*_args, **_kwargs):
        raise AssertionError("a draw was made")
    for name in ("random", "uniform", "randrange", "randint", "choice", "choices", "shuffle", "sample", "gauss", "gamma",
                 "dirichlet", "getrandbits"):
        monkeypatch.setattr(draw, name, refuse)


def test_the_schedule_fields_are_one_age_for_each_milestone():
    assert evolution.SCHEDULE_FIELDS == tuple(f"life_age_{key}_gyr" for key in tuning.LIFE_MILESTONES)
    for scale in tuning.EVOLUTIONARY_TIMELINES:
        ages = evolution.schedule(scale)
        assert len(ages) == len(evolution.SCHEDULE_FIELDS)
        assert list(ages) == sorted(ages), "a milestone is never reached before an easier one"


def test_a_body_has_a_pace_and_a_schedule_exactly_when_it_has_a_chemistry():
    with_life = without = 0
    for body in _bodies(40, HABITABLE_WORLD=True, MOONS=True):
        stored = tuple(getattr(body, field) for field in evolution.SCHEDULE_FIELDS)
        if body.life_chemical:
            assert body.evolutionary_pace is not None, body.name
            assert stored == evolution.schedule(body.evolutionary_pace), body.name
            with_life += 1
        else:
            assert body.evolutionary_pace is None and stored == (None,) * len(stored), body.name
            without += 1
    assert with_life > 3 and without > 3


def test_a_body_the_index_gives_no_microbial_life_loses_the_pace_it_had():
    body = _a_world_with_life()
    assert body.evolutionary_pace is not None
    body.phi_bio = 0.0
    planet_life.apply_life_data(body)
    assert body.life_chemical is None and body.evolutionary_pace is None and body.life_stage is None
    assert all(getattr(body, field) is None for field in evolution.SCHEDULE_FIELDS)


def test_an_edit_keeps_the_chemistry_and_the_pace_while_they_still_fit(monkeypatch):
    body = _a_world_with_life()
    before = (body.life_chemical, body.evolutionary_pace, body.life_stage, body.life_milestone_age_gyr,
              body.life_system_age_gyr, *(getattr(body, field) for field in evolution.SCHEDULE_FIELDS))
    _no_draw(monkeypatch)
    validation.reapply_life(body)
    after = (body.life_chemical, body.evolutionary_pace, body.life_stage, body.life_milestone_age_gyr,
             body.life_system_age_gyr, *(getattr(body, field) for field in evolution.SCHEDULE_FIELDS))
    assert after == before


def test_an_edit_draws_a_chemistry_again_only_when_the_old_one_no_longer_fits(monkeypatch):
    body = _a_world_with_life()
    old = body.life_chemical
    other = next(name for name in tuning.LIFE_CHEMICALS if name != old)
    monkeypatch.setattr(planet_life, "get_viable_life_chemicals", lambda _planet, _spectral=None: {other: 1.0})
    validation.reapply_life(body)
    assert body.life_chemical == other
    assert body.reflection_spectrum_visible == tuning.LIFE_CHEMICALS[other]["reflection_spectrum_visible"]
    assert body.evolutionary_pace in evolution.candidate_scales(body.star, planet_life.get_chemical_scales(other))[0]
    assert tuple(getattr(body, field) for field in evolution.SCHEDULE_FIELDS) == evolution.schedule(body.evolutionary_pace)


def test_an_edit_draws_a_pace_again_only_when_the_old_one_no_longer_fits():
    body = _a_world_with_life()
    allowed = evolution.candidate_scales(body.star, planet_life.get_chemical_scales(body.life_chemical))[0]
    unfit = next((scale for scale in tuning.EVOLUTIONARY_TIMELINES if scale not in allowed), None)
    if unfit is None:
        pytest.skip("this star and chemistry allow every pace")
    body.evolutionary_pace = unfit
    validation.reapply_life(body)
    assert body.evolutionary_pace in allowed
    assert tuple(getattr(body, field) for field in evolution.SCHEDULE_FIELDS) == evolution.schedule(body.evolutionary_pace)


def test_the_pace_drawn_is_one_of_the_candidates():
    seen = 0
    for body in _bodies(40, HABITABLE_WORLD=True):
        if body.life_chemical:
            assert body.evolutionary_pace in evolution.candidate_scales(
                body.star, planet_life.get_chemical_scales(body.life_chemical))[0]
            seen += 1
    assert seen > 3
