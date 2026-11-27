# tests/test_life_ages.py

"""GEN.1746 slice 2: the life ages of a world, the stretches between its mass extinctions (`physics/life_ages.py`)."""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import air_events, life_ages
from tests.fuzz_support import deterministic_entropy


def _event(age, kind=air_events.IMPACT, share=0.4):
    return air_events.Event(age, kind, 10.0, share)


def _a_rocky_world():
    for seed in range(40):
        config = SystemConfig()
        config.HABITABLE_WORLD = True
        with deterministic_entropy(7300 + seed):
            system = StarSystem(system_config=config)
        for planet in system.primary.orbiters:
            if planet.body_type == "t" and planet.star.age and planet.star.age > 2.0:
                return planet
    pytest.fail("no generated rocky world")


def test_the_impact_that_kills_the_threshold_share_is_the_one_the_earth_has_at_the_curves_interval():
    diameter = life_ages.smallest_impact_km(0.30)
    assert air_events.impact_kill_share(diameter) == pytest.approx(0.30)
    assert air_events.impact_rate_per_gyr(diameter, air_events.constants.EARTH_RADIUS_M) == pytest.approx(100.0)          # once per 10 Myr
    assert life_ages.smallest_impact_km(0.65) > life_ages.smallest_impact_km(0.30) > life_ages.smallest_impact_km(0.15)
    assert life_ages.AGE_KILL_THRESHOLD == 0.60          # the Boss's pick of 11:18 PM EST 11-26 (GEN.1852), the Big Five's rate: once per 72 Myr, an impactor of about 10 km
    assert 10.0 < life_ages.smallest_impact_km() < 11.0


def test_only_the_events_that_kill_enough_end_a_life_age_and_the_ones_inside_a_recovery_lag_are_one_extinction():
    events = [_event(1.0, share=0.2), _event(1.5, air_events.FLOOD_BASALT, 0.4), _event(1.5 + 0.5 * air_events.RECOVERY_GYR, share=0.5), _event(2.0, share=0.35), _event(9.0, share=0.9)]
    found = life_ages.extinctions_of(events, 0.5, 4.0, threshold=0.30)
    assert found == [(1.5, air_events.IMPACT, 0.5), (2.0, air_events.IMPACT, 0.35)]
    assert life_ages.extinctions_of(events, 0.5, 10.0) == [(9.0, air_events.IMPACT, 0.9)]          # the default is 0.60 (GEN.1852): only the end-Permian-size kill is left
    assert life_ages.extinctions_of(events, 0.5, 4.0, threshold=0.45) == [(1.5 + 0.5 * air_events.RECOVERY_GYR, air_events.IMPACT, 0.5)]
    assert life_ages.extinctions_of(events, 1.6, 1.9, threshold=0.30) == []


def test_the_ages_run_from_the_start_of_life_one_recovery_lag_after_each_extinction_to_the_star_s_age(monkeypatch):
    body = _a_rocky_world()
    events = [_event(1.0), _event(2.0, air_events.FLOOD_BASALT, 0.43)]
    monkeypatch.setattr(air_events, "events_of", lambda body_, start, end, min_impact_km=5.0: events)
    ages = life_ages.ages_of(body, 0.7, threshold=0.30)
    assert [age.number for age in ages] == [1, 2, 3]
    assert ages[0] == life_ages.LifeAge(1, 0.7, 1.0, air_events.IMPACT, 0.4)
    assert ages[1].start_gyr == pytest.approx(1.0 + air_events.RECOVERY_GYR) and ages[1].end_gyr == 2.0 and ages[1].cause == air_events.FLOOD_BASALT
    assert ages[2].start_gyr == pytest.approx(2.0 + air_events.RECOVERY_GYR) and ages[2].end_gyr is None and ages[2].cause is None and ages[2].kill_share is None


def test_a_body_with_no_life_or_no_age_has_no_ages_and_a_quiet_one_has_one_open_age(monkeypatch):
    body = _a_rocky_world()
    assert life_ages.ages_of(body, None) == [] and life_ages.ages_of(body, body.star.age + 1.0) == []
    monkeypatch.setattr(air_events, "events_of", lambda *args, **kwargs: [])
    only = life_ages.ages_of(body, 1.0)
    assert only == [life_ages.LifeAge(1, 1.0, None, None, None)]


def test_the_ages_of_a_generated_world_tile_its_life_in_order_and_repeat():
    body = _a_rocky_world()
    begin = 1.0
    ages = life_ages.ages_of(body, begin)
    assert ages == life_ages.ages_of(body, begin)
    assert ages and ages[0].start_gyr == begin and ages[-1].end_gyr is None
    assert all(age.kill_share >= life_ages.AGE_KILL_THRESHOLD for age in ages[:-1])
    for earlier, later in zip(ages, ages[1:]):
        assert later.start_gyr == pytest.approx(earlier.end_gyr + air_events.RECOVERY_GYR) and later.number == earlier.number + 1
