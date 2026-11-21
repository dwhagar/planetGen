# tests/test_air_events.py

"""GEN.1746 slice 1: the shocks of a world's history and the share of life each kills (`physics/air_events.py`)."""

import statistics

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import air_events, constants, surface_pressure
from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy


def _a_rocky_world():
    for seed in range(40):
        config = SystemConfig()
        config.HABITABLE_WORLD = True
        with deterministic_entropy(7300 + seed):
            system = StarSystem(system_config=config)
        for planet in system.planets:
            if planet.body_type == "t" and planet.star.age and planet.star.age > 2.0:
                return planet
    pytest.fail("no generated rocky world")


def test_the_kill_curve_passes_through_its_three_points_straight_in_the_log_and_is_held_at_the_ends():
    assert air_events.kill_share(1.0) == pytest.approx(0.10)
    assert air_events.kill_share(10.0) == pytest.approx(0.30)
    assert air_events.kill_share(100.0) == pytest.approx(0.65)
    assert air_events.kill_share(0.01) == 0.10 and air_events.kill_share(1.0e5) == 0.65
    assert 0.10 < air_events.kill_share(3.0) < 0.30 < air_events.kill_share(30.0) < 0.65
    assert air_events.kill_share(10.0 ** 0.5) == pytest.approx(0.20)


def test_a_larger_event_recurs_less_often_and_kills_more():
    kills = [air_events.impact_kill_share(diameter) for diameter in (1.0, 2.0, 5.0, 10.0, 30.0)]
    assert kills == sorted(kills) and kills[0] == pytest.approx(0.10) and kills[-1] > kills[2]
    rates = [air_events.impact_rate_per_gyr(diameter, constants.EARTH_RADIUS_M) for diameter in (1.0, 5.0, 10.0)]
    assert rates[0] == pytest.approx(2000.0) and rates[1] == pytest.approx(1000.0 / 15.0) and rates[2] < rates[1]
    assert air_events.flood_basalt_kill_share() == pytest.approx(air_events.kill_share(25.0))


def test_the_impact_flux_is_the_earth_s_now_from_the_earth_s_age_on_and_hundreds_of_times_that_at_3_9_ga_before_it():
    now = tuning.STELLAR_SOLAR_SYSTEM_AGE_GY
    assert air_events.impact_flux_ratio(now) == pytest.approx(1.0) and air_events.impact_flux_ratio(now + 3.0) == pytest.approx(1.0)
    assert 200.0 < air_events.impact_flux_ratio(now - 3.9) < 300.0
    ratios = [air_events.impact_flux_ratio(age) for age in (4.0, 3.0, 2.0, 1.0, 0.5)]
    assert ratios == sorted(ratios)


def test_flood_basalts_follow_the_heat_of_the_interior_and_stop_without_it():
    assert air_events.flood_basalt_rate_per_gyr(surface_pressure.HEAT_FLUX_EARTH_W_M2) == pytest.approx(air_events.FLOOD_BASALT_PER_GYR)
    assert air_events.flood_basalt_rate_per_gyr(2.0 * surface_pressure.HEAT_FLUX_EARTH_W_M2) == pytest.approx(2.0 * air_events.FLOOD_BASALT_PER_GYR)
    assert air_events.flood_basalt_rate_per_gyr(0.0) == 0.0 and air_events.flood_basalt_rate_per_gyr(None) == 0.0


def test_a_body_melted_by_tides_erupts_all_the_time_and_not_billions_of_times_faster():
    """A tidally heated body (25 MW/m2 in the 9101 sweep) asked for 6.5e10 flood basalts and the draw never ended (12 GB): one takes a million years, so the rate stops at one at a time."""
    assert air_events.FLOOD_BASALT_MAX_PER_GYR == pytest.approx(1000.0)
    assert air_events.flood_basalt_rate_per_gyr(2.5e7) == air_events.FLOOD_BASALT_MAX_PER_GYR
    assert len(air_events._waits(draw.Stream("tidal-melt"), [air_events.flood_basalt_rate_per_gyr(2.5e7)], [4.0])) < 6000


def test_a_poisson_process_has_its_expected_count_over_a_constant_and_a_stepped_rate():
    counts = []
    for seed in range(300):
        stream = draw.Stream(f"air-events-test:{seed}")
        counts.append(len(air_events._waits(stream, [10.0, 0.0, 30.0], [1.0, 1.0, 0.5])))
    assert statistics.mean(counts) == pytest.approx(25.0, rel=0.08)
    stream = draw.Stream("air-events-test-ordered")
    ages = air_events._waits(stream, [10.0, 0.0, 30.0], [1.0, 1.0, 0.5])
    assert ages == sorted(ages) and all(0.0 <= age < 2.5 for age in ages) and not any(1.0 <= age < 2.0 for age in ages)
    assert air_events._waits(draw.Stream("none"), [0.0, 0.0], [1.0, 1.0]) == []


def test_the_spans_cover_the_range_without_a_gap_and_stay_short():
    pairs = air_events.spans(0.3, 4.5)
    assert sum(span for _, span in pairs) == pytest.approx(4.2)
    assert max(span for _, span in pairs) <= air_events.EVENT_STEP_GYR + 1e-12
    assert pairs[0][0] > 0.3 and pairs[-1][0] < 4.5 and air_events.spans(2.0, 1.0) == []


def test_the_events_of_a_world_are_ordered_sized_and_drawn_on_the_body_s_own_stream():
    body = _a_rocky_world()
    with deterministic_entropy(11):
        first = draw.random()
    with deterministic_entropy(11):
        events = air_events.events_of(body, 0.5)
        after = draw.random()
    assert first == after                                           # the events used none of the run's draws
    assert events == air_events.events_of(body, 0.5)                # and the same body has the same events again
    assert [event.age_gyr for event in events] == sorted(event.age_gyr for event in events)
    assert all(0.5 <= event.age_gyr <= body.star.age for event in events)
    flood = [event for event in events if event.kind == air_events.FLOOD_BASALT]
    impacts = [event for event in events if event.kind == air_events.IMPACT]
    assert all(event.size == air_events.FLOOD_BASALT_VOLUME_KM3 for event in flood)
    assert all(5.0 <= event.size <= air_events.IMPACT_LARGEST_KM for event in impacts)
    assert all(0.10 <= event.kill_share <= 0.65 for event in events)
    assert len(flood) > 0 or body.star.age < 1.0
    assert len(air_events.events_of(body, 0.5, min_impact_km=20.0)) <= len(events)
