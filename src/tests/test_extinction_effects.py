"""GEN.1746 slice 4: what a mass extinction does to the animals and the civilization of a world (`physics.extinction_effects`, `generation.life._reset_by_extinctions`)."""

from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import life
from planetgen.physics import air_events, extinction_effects, life_ages

BACK = air_events.RECOVERY_GYR + extinction_effects.REGROWTH_GYR


def _ages(animals, civilization=None):
    return (0.3, 1.0, 2.5, animals, civilization)


def _entry(end, share, number=1):
    return {"number": number, "start_gyr": 0.0, "end_gyr": end, "cause": "impact", "kill_share": share}


def test_only_the_extinctions_that_reach_the_reset_share_count_and_the_open_age_is_not_one():
    entries = [_entry(2.0, 0.40), _entry(1.0, 0.65), _entry(3.0, 0.90), _entry(None, None), _entry(2.5, None)]
    assert extinction_effects.reset_extinctions(entries) == [(1.0, 0.65), (3.0, 0.90)]
    assert extinction_effects.reset_extinctions(None) == []


def test_no_extinction_leaves_the_schedule_as_it_was():
    ages = _ages(3.0, 3.5)
    assert extinction_effects.regrown(ages, []) == ages


def test_an_extinction_before_the_animals_came_changes_nothing():
    ages = _ages(3.0, 3.5)
    assert extinction_effects.regrown(ages, [(2.0, 0.7)]) == ages


def test_the_animals_return_one_lag_and_one_regrowth_after_the_extinction_and_the_civilization_after_the_same_wait():
    ages = _ages(2.0, 2.4)          # the animals at 2.0, a civilization 0.4 Gyr later
    out = extinction_effects.regrown(ages, [(3.0, 0.7)])
    assert out[:3] == ages[:3]      # the microbes and the complex cells are never touched
    assert out[extinction_effects.ANIMALS] == pytest.approx(3.0 + BACK)
    assert out[extinction_effects.CIVILIZATION] == pytest.approx(3.0 + BACK + 0.4)


def test_a_civilization_that_had_not_come_yet_starts_again_the_same_wait_after_the_animals_are_back():
    ages = _ages(2.0, 5.0)          # the civilization comes 3.0 Gyr after the animals, so after the extinction at 3.0
    out = extinction_effects.regrown(ages, [(3.0, 0.7)])
    assert out[extinction_effects.CIVILIZATION] == pytest.approx(3.0 + BACK + 3.0)


def test_a_world_that_never_has_a_civilization_has_none_after_either():
    out = extinction_effects.regrown(_ages(2.0, None), [(3.0, 0.7)])
    assert out[extinction_effects.CIVILIZATION] is None and out[extinction_effects.ANIMALS] == pytest.approx(3.0 + BACK)


def test_a_world_that_never_has_animals_is_left_alone():
    ages = _ages(None, None)
    assert extinction_effects.regrown(ages, [(3.0, 0.9)]) == ages


def test_an_extinction_while_the_animals_are_still_returning_does_not_delay_them_again():
    out = extinction_effects.regrown(_ages(2.0, 2.4), [(3.0, 0.7), (3.0 + BACK / 2, 0.8)])
    assert out[extinction_effects.ANIMALS] == pytest.approx(3.0 + BACK)


def test_two_extinctions_after_the_animals_are_back_reset_them_twice():
    out = extinction_effects.regrown(_ages(2.0, 2.4), [(3.0, 0.7), (3.5, 0.9)])
    assert out[extinction_effects.ANIMALS] == pytest.approx(3.5 + BACK)
    assert out[extinction_effects.CIVILIZATION] == pytest.approx(3.5 + BACK + 0.4)


def test_the_reset_share_is_the_top_of_the_kill_curve():
    assert extinction_effects.RESET_KILL_SHARE == air_events.KILL_CURVE[-1][1]
    assert tuning.LIFE_MILESTONES[extinction_effects.ANIMALS] == "multicellularity"
    assert tuning.LIFE_MILESTONES[extinction_effects.CIVILIZATION] == "technological_civilization"


def _world(monkeypatch, entries, star_age=4.0):
    stored = []
    monkeypatch.setattr(life_ages, "ages_of", lambda body, begin, end_gyr=None, threshold=life_ages.AGE_KILL_THRESHOLD: ("ages", begin))
    monkeypatch.setattr(life_ages, "entries_of", lambda ages: entries)
    monkeypatch.setattr(life, "_store_schedule", lambda planet, cands: stored.append([c.ages for c in cands]))
    return SimpleNamespace(star=SimpleNamespace(age=star_age)), stored


def test_the_life_pass_resets_the_schedules_it_drew_and_stores_the_dominant_one_again(monkeypatch):
    planet, stored = _world(monkeypatch, [_entry(3.0, 0.7)])
    drawn = SimpleNamespace(waits=(1, 1, 1, 1, 1), ages=_ages(2.0, 2.4))
    minor = SimpleNamespace(waits=(1, 1, 1, 1, 1), ages=(0.4, 1.2, 2.6, 2.5, None))
    life._reset_by_extinctions(planet, [drawn, minor])
    assert drawn.ages[extinction_effects.ANIMALS] == pytest.approx(3.0 + BACK)
    assert minor.ages[extinction_effects.ANIMALS] == pytest.approx(3.0 + BACK)
    assert len(stored) == 1


def test_a_schedule_kept_from_a_stored_body_is_not_reset_again(monkeypatch):
    planet, stored = _world(monkeypatch, [_entry(3.0, 0.7)])
    kept = SimpleNamespace(waits=None, ages=_ages(2.0, 2.4))
    life._reset_by_extinctions(planet, [kept])
    assert kept.ages == _ages(2.0, 2.4) and stored == []


def test_a_world_with_no_animals_by_the_stars_age_has_nothing_to_reset_and_asks_for_no_extinctions(monkeypatch):
    planet, stored = _world(monkeypatch, [_entry(3.0, 0.7)])
    monkeypatch.setattr(life_ages, "ages_of", lambda *a, **k: pytest.fail("the extinctions are drawn only for a world that has animals"))
    cand = SimpleNamespace(waits=(1, 1, 1, 1, 1), ages=_ages(5.0, None))
    life._reset_by_extinctions(planet, [cand])
    assert cand.ages == _ages(5.0, None) and stored == []


def test_the_extinctions_are_read_from_the_age_the_first_producers_began(monkeypatch):
    planet, _ = _world(monkeypatch, [_entry(3.0, 0.7)])
    seen = []
    monkeypatch.setattr(life_ages, "ages_of", lambda body, begin, **kw: seen.append(begin) or ())
    first = SimpleNamespace(waits=(1,) * 5, ages=(0.3, 1.0, 2.5, 2.0, None))
    other = SimpleNamespace(waits=(1,) * 5, ages=(0.4, 0.8, 2.6, 2.5, None))
    life._reset_by_extinctions(planet, [first, other])
    assert seen == [0.8]
