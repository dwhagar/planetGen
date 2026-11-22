"""GEN.1746 slice 4: what a mass extinction does to the animals and the civilization of a world (`physics.extinction_effects`, `generation.life._reset_by_extinctions`)."""

from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import life
from planetgen.physics import air_events, extinction_effects, life_ages

BACK = air_events.RECOVERY_GYR + extinction_effects.REGROWTH_GYR


def _ages(animals, civilization=None):
    return (0.3, 1.0, 2.5, animals, civilization)


def test_no_extinction_leaves_the_schedule_as_it_was():
    ages = _ages(3.0, 3.5)
    assert extinction_effects.regrown(ages, []) == ages


def test_an_extinction_before_the_animals_came_changes_nothing():
    ages = _ages(3.0, 3.5)
    assert extinction_effects.regrown(ages, [(2.0, 50.0)]) == ages


def test_the_animals_return_one_lag_and_one_regrowth_after_the_extinction_and_the_civilization_after_the_same_wait():
    ages = _ages(2.0, 2.4)          # the animals at 2.0, a civilization 0.4 Gyr later
    out = extinction_effects.regrown(ages, [(3.0, 50.0)])
    assert out[:3] == ages[:3]      # the microbes and the complex cells are never touched
    assert out[extinction_effects.ANIMALS] == pytest.approx(3.0 + BACK)
    assert out[extinction_effects.CIVILIZATION] == pytest.approx(3.0 + BACK + 0.4)


def test_a_civilization_that_had_not_come_yet_starts_again_the_same_wait_after_the_animals_are_back():
    ages = _ages(2.0, 5.0)          # the civilization comes 3.0 Gyr after the animals, so after the extinction at 3.0
    out = extinction_effects.regrown(ages, [(3.0, 50.0)])
    assert out[extinction_effects.CIVILIZATION] == pytest.approx(3.0 + BACK + 3.0)


def test_a_world_that_never_has_a_civilization_has_none_after_either():
    out = extinction_effects.regrown(_ages(2.0, None), [(3.0, 50.0)])
    assert out[extinction_effects.CIVILIZATION] is None and out[extinction_effects.ANIMALS] == pytest.approx(3.0 + BACK)


def test_a_world_that_never_has_animals_is_left_alone():
    ages = _ages(None, None)
    assert extinction_effects.regrown(ages, [(3.0, 50.0)]) == ages


def test_an_extinction_while_the_animals_are_still_returning_does_not_delay_them_again():
    out = extinction_effects.regrown(_ages(2.0, 2.4), [(3.0, 50.0), (3.0 + BACK / 2, 60.0)])
    assert out[extinction_effects.ANIMALS] == pytest.approx(3.0 + BACK)


def test_two_extinctions_after_the_animals_are_back_reset_them_twice():
    out = extinction_effects.regrown(_ages(2.0, 2.4), [(3.0, 50.0), (3.5, 60.0)])
    assert out[extinction_effects.ANIMALS] == pytest.approx(3.5 + BACK)
    assert out[extinction_effects.CIVILIZATION] == pytest.approx(3.5 + BACK + 0.4)


def test_the_size_that_resets_is_the_one_the_earth_has_once_in_a_gigayear_and_is_bigger_than_any_extinction_of_its_animal_record():
    km = extinction_effects.reset_impact_km()
    assert 30.0 < km < 45.0                      # about 37 km on the size law; the largest impactor among the Big Five is about 10 km
    assert km > life_ages.smallest_impact_km(0.65)
    per_gyr = air_events.IMPACT_FIVE_KM_PER_GYR * (km / 5.0) ** -air_events.IMPACT_SIZE_EXPONENT
    assert per_gyr == pytest.approx(1.0)
    assert tuning.LIFE_MILESTONES[extinction_effects.ANIMALS] == "multicellularity"
    assert tuning.LIFE_MILESTONES[extinction_effects.CIVILIZATION] == "technological_civilization"


def test_only_the_impacts_of_the_reset_size_count_and_flood_basalts_never_do(monkeypatch):
    km = extinction_effects.reset_impact_km()
    events = [air_events.Event(1.0, air_events.IMPACT, km - 1.0, 0.65), air_events.Event(2.0, air_events.IMPACT, km, 0.65),
              air_events.Event(3.0, air_events.FLOOD_BASALT, 2.5e6, 0.43), air_events.Event(3.5, air_events.IMPACT, 90.0, 0.65)]
    seen = []
    monkeypatch.setattr(air_events, "events_of", lambda body, begin, end, min_impact_km: seen.append((begin, end, min_impact_km)) or events)
    body = SimpleNamespace(star=SimpleNamespace(age=4.0))
    assert extinction_effects.extinctions_of(body, 1.5) == [(2.0, km), (3.5, 90.0)]
    assert seen == [(1.5, 4.0, life_ages.smallest_impact_km())]       # the same events the life ages are cut at
    assert extinction_effects.extinctions_of(body, None) == [] and extinction_effects.extinctions_of(SimpleNamespace(), 1.0) == []


def _world(monkeypatch, extinctions, star_age=4.0):
    stored = []
    monkeypatch.setattr(extinction_effects, "extinctions_of", lambda body, begin, end_gyr=None: extinctions)
    monkeypatch.setattr(life, "_store_schedule", lambda planet, cands: stored.append([c.ages for c in cands]))
    return SimpleNamespace(star=SimpleNamespace(age=star_age)), stored


def test_the_life_pass_resets_the_schedules_it_drew_and_stores_the_dominant_one_again(monkeypatch):
    planet, stored = _world(monkeypatch, [(3.0, 50.0)])
    drawn = SimpleNamespace(waits=(1, 1, 1, 1, 1), ages=_ages(2.0, 2.4))
    minor = SimpleNamespace(waits=(1, 1, 1, 1, 1), ages=(0.4, 1.2, 2.6, 2.5, None))
    life._reset_by_extinctions(planet, [drawn, minor])
    assert drawn.ages[extinction_effects.ANIMALS] == pytest.approx(3.0 + BACK)
    assert minor.ages[extinction_effects.ANIMALS] == pytest.approx(3.0 + BACK)
    assert len(stored) == 1


def test_a_schedule_kept_from_a_stored_body_is_not_reset_again(monkeypatch):
    planet, stored = _world(monkeypatch, [(3.0, 50.0)])
    kept = SimpleNamespace(waits=None, ages=_ages(2.0, 2.4))
    life._reset_by_extinctions(planet, [kept])
    assert kept.ages == _ages(2.0, 2.4) and stored == []


def test_a_world_with_no_animals_by_the_stars_age_asks_for_no_extinctions(monkeypatch):
    planet, stored = _world(monkeypatch, [(3.0, 50.0)])
    monkeypatch.setattr(extinction_effects, "extinctions_of", lambda *a, **k: pytest.fail("the extinctions are drawn only for a world that has animals"))
    cand = SimpleNamespace(waits=(1, 1, 1, 1, 1), ages=_ages(5.0, None))
    life._reset_by_extinctions(planet, [cand])
    assert cand.ages == _ages(5.0, None) and stored == []


def test_the_extinctions_are_read_from_the_age_the_first_producers_began(monkeypatch):
    planet, _ = _world(monkeypatch, [])
    seen = []
    monkeypatch.setattr(extinction_effects, "extinctions_of", lambda body, begin, end_gyr=None: seen.append(begin) or [])
    first = SimpleNamespace(waits=(1,) * 5, ages=(0.3, 1.0, 2.5, 2.0, None))
    other = SimpleNamespace(waits=(1,) * 5, ages=(0.4, 0.8, 2.6, 2.5, None))
    life._reset_by_extinctions(planet, [first, other])
    assert seen == [0.8]


def test_a_world_with_no_extinction_of_that_size_is_left_alone_and_stores_nothing(monkeypatch):
    planet, stored = _world(monkeypatch, [])
    cand = SimpleNamespace(waits=(1,) * 5, ages=_ages(2.0, 2.4))
    life._reset_by_extinctions(planet, [cand])
    assert cand.ages == _ages(2.0, 2.4) and stored == []
