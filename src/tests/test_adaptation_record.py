"""GEN.1808: what the life of a finished body had to become to live where it does (`planetgen.physics.adaptation_record`)."""

from types import SimpleNamespace

import pytest

from planetgen.physics import adaptation_record as ar, life_adaptations as la, surface_refuge

from tests.test_life_classes import world
from tests.test_surface_refuge import _dry, _frozen_sea


def _needs(values, class_id, age=3.0):
    return ar.needs(surface_refuge.refuge_of(values), values, class_id, age)


def test_an_open_sea_and_a_body_with_no_refuge_need_no_adaptation():
    assert _needs(world(), "1a") == []
    assert _needs(world(hydrosphere="dry"), "1a") == []


def test_a_class_without_a_record_gets_none():
    assert _needs(_dry(0.0, temperature_k=260.0), "4a") == [] and _needs(_dry(0.0, temperature_k=260.0), "3a") == []


def test_a_seep_in_the_ordinary_range_needs_nothing():
    """Wet all year at 300 K: inside the cold and hot limits of an ordinary life, with the water of the whole year."""
    assert _needs(_dry(0.0), "1a") == []


def test_a_cold_seep_needs_the_lowest_grade_that_reaches_its_temperature_and_the_rows_that_give_it():
    need, = _needs(_dry(0.0, temperature_k=260.0), "1a")
    assert (need.axis, need.habitat, need.grade) == ("t_cold", "seep", 3)
    assert need.demand == pytest.approx(260.0, abs=0.5) and need.rows == ("C8",)          # the lichen's 253 K, at the extreme grade
    assert la.limit("t_cold", 2, class_id="1a", life_age_gyr=3.0) > need.demand >= la.limit("t_cold", 3, class_id="1a", life_age_gyr=3.0)
    warmer, = _needs(_dry(0.0, temperature_k=265.0), "1a")
    assert warmer.grade == 2 and "C8" in warmer.rows


def test_a_life_too_young_for_the_adaptation_is_past_the_table_and_an_age_that_has_it_is_not():
    values = _dry(0.0, temperature_k=260.0)
    young, = _needs(values, "1a")
    assert young.grade == 3
    early, = ar.needs(surface_refuge.refuge_of(values), values, "1a", 0.1)          # the lichen needs land mats, 0.6 Gyr after the first life
    assert early.grade is None and early.rows == ()


def test_a_class_the_table_has_no_rows_for_is_past_the_table_not_unneeded():
    need, = _needs(_dry(0.0, temperature_k=260.0), "2")
    assert need.axis == "t_cold" and need.grade is None and need.rows == ()


def test_the_sea_under_the_ice_needs_the_light_floor_and_the_habitat_rows():
    values = _frozen_sea(270.0)
    refuge = surface_refuge.refuge_of(values)
    found = ar.needs(refuge, values, "1b", 3.0)
    light = next(need for need in found if need.axis == "light_floor")
    assert light.habitat == "water_under_ice" and light.demand == pytest.approx(refuge.under_ice_w_m2 * ar.UMOL_PER_W)
    assert light.grade == 3 and light.rows == ("C1",)          # 0.08 umol: the tolerant grade's 0.2 does not reach it, the extreme grade's 0.04 does
    assert {need.rows[0] for need in found if need.axis is None} == {"C2", "C3", "C4"}
    assert all(need.demand is None and need.grade is None for need in found if need.axis is None)


def test_the_refuge_carries_the_conditions_the_record_reads():
    assert surface_refuge.refuge_of(world()).habitats == ("open_water",)
    seep = surface_refuge.refuge_of(_dry(0.0, temperature_k=288.0))
    assert seep.habitats == ("seep",) and seep.coldest_k == pytest.approx(288.0, abs=0.5) and seep.hottest_k == pytest.approx(288.0, abs=0.5)
    under = surface_refuge.refuge_of(_frozen_sea(270.0))
    assert under.habitats == ("water_under_ice",) and under.under_ice_w_m2 >= surface_refuge.LIGHT_FLOOR_W_M2
    assert surface_refuge.refuge_of(world(hydrosphere="dry")) == surface_refuge.NONE


def test_the_temperature_range_holds_most_of_the_wet_time_and_not_the_one_hottest_step():
    import numpy as np

    temperature = np.linspace(270.0, 330.0, 1001)
    cold, hot = surface_refuge._range_of([(temperature, np.ones_like(temperature))])
    assert cold == pytest.approx(271.5, abs=0.2) and hot == pytest.approx(328.5, abs=0.2)
    assert surface_refuge._range_of([]) == (None, None)
    heavy = surface_refuge._range_of([(np.array([280.0, 340.0]), np.array([99.0, 1.0]))])      # a hot place of 1 percent of the weight is left out
    assert heavy == (280.0, 280.0) or heavy[1] < 340.0


def test_the_record_is_a_list_of_plain_dicts_and_none_for_none():
    assert ar.entries_of([], "1a") is None
    entries = ar.entries_of(_needs(_dry(0.0, temperature_k=260.0), "1a"), "1a")
    assert entries == [{"life_class": "1a", "axis": "t_cold", "habitat": "seep", "demand": pytest.approx(260.0, abs=0.5), "grade": 3, "rows": ["C8"]}]


def test_the_light_of_the_table_and_of_the_gate_agree_at_the_gates_floor():
    assert surface_refuge.LIGHT_FLOOR_W_M2 * ar.UMOL_PER_W == pytest.approx(0.01)


def test_the_record_of_a_body_reads_its_refuge_once_and_has_none_for_a_body_with_no_surface_or_a_class_with_no_record():
    body = SimpleNamespace(surface_temperature=0.0)
    assert ar.record_of(body, ["1a"]) is None and ar.record_of(body, ["4a"]) is None
