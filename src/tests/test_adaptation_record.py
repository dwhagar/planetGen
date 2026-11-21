"""GEN.1808: what the life of a finished body had to become to live where it does (`planetgen.physics.adaptation_record`)."""

from types import SimpleNamespace

import pytest

from planetgen.physics import acid_medium, adaptation_record as ar, life_adaptations as la, life_classes, surface_refuge

from tests.test_life_classes import world
from tests.test_surface_refuge import _dry, _frozen_sea


def _needs(values, class_id, age=3.0):
    return ar.needs(surface_refuge.refuge_of(values), values, class_id, age)


def test_an_open_sea_and_a_body_with_no_refuge_need_no_adaptation():
    assert _needs(world(), "1a") == []
    assert _needs(world(hydrosphere="dry"), "1a") == []


def test_a_class_id_that_is_no_class_gets_none():
    assert _needs(_dry(0.0, temperature_k=260.0), "x") == []


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


def test_the_record_of_a_body_with_no_surface_is_none():
    body = SimpleNamespace(surface_temperature=0.0)
    assert ar.record_of(body, ["1a"]) is None and ar.record_of(body, ["4a"]) is None


# --- the chemical classes (GEN.1811) ---

def _chemical(class_id, age=3.0, **over):
    values = world(**over)
    return ar.needs(None, values, class_id, age)


def test_every_class_of_the_rule_file_has_a_record_and_the_table_knows_every_class():
    assert set(ar.CLASSES) == set(life_classes.IDS)
    assert set(la.KNOWN_CLASSES) == set(life_classes.IDS)
    assert set(la.CHEMICAL_CLASSES) == {c.id for c in life_classes.CLASSES if c.energy["source"] != "light"}


@pytest.mark.parametrize("class_id", la.WATER_CHEMICAL_CLASSES)
def test_a_chemical_life_in_ordinary_water_needs_nothing(class_id):
    assert _chemical(class_id, surface_temperature_k=300.0) == []


def test_a_chemical_life_in_cold_water_needs_the_cold_rows_and_the_one_row_below_the_growth_edge():
    brine = _chemical("3a", surface_temperature_k=265.0)
    need, = brine
    assert (need.axis, need.habitat, need.grade, need.rows) == ("t_cold", "brine", 2, ("C16",))          # the measured growth edge, -10 C, with a cold-active enzyme set
    deep, = _chemical("4b", surface_temperature_k=255.0)
    assert deep.grade == 3 and deep.rows == ("C17",) and deep.demand == pytest.approx(255.0)
    past, = _chemical("3b", surface_temperature_k=240.0)
    assert past.grade is None and past.rows == ()          # past the metabolism of the coldest brine microbe


def test_a_chemical_life_in_hot_water_needs_hyperthermophily_up_to_the_hottest_growth():
    hot, = _chemical("3a", surface_temperature_k=380.0)
    assert (hot.axis, hot.habitat, hot.grade, hot.rows) == ("t_hot", "hot_pool", 3, ("H3",))
    warm, = _chemical("4a", surface_temperature_k=350.0)
    assert warm.grade == 2 and warm.rows == ("H3",)
    over, = _chemical("3a", surface_temperature_k=400.0)
    assert over.grade is None          # 122 C is the record; the gate's 395 K is at it, a hotter habitat is not carried


def test_the_cold_and_hot_rows_belong_to_the_classes_of_water_only():
    assert la.limit("t_cold", 3, class_id="3a") == pytest.approx(248.15) and la.limit("t_hot", 3, class_id="3a") == pytest.approx(395.15)
    assert la.limit("t_hot", 3, class_id="5b") < 330.0 and la.limit("t_cold", 3, class_id="7") > 270.0


def _acid(monkeypatch, weight_percent, temperature_k):
    monkeypatch.setattr(acid_medium, "medium", lambda values, cloud_only=False: acid_medium.AcidMedium(
        level=acid_medium.CLOUD, layer="cloud", temperature_k=temperature_k, pressure_pa=1.0e5, weight_percent=weight_percent))


def test_the_acid_classes_are_past_the_table_on_strength_and_water_and_the_record_says_so(monkeypatch):
    _acid(monkeypatch, 85.0, 330.0)
    found = _chemical("5a")
    by_axis = {need.axis: need for need in found}
    assert set(by_axis) == {"t_hot", "acid_strength", "water_activity"}
    assert all(need.habitat == ar.CHEMICAL_HABITAT and need.grade is None and need.rows == () for need in found)
    assert by_axis["acid_strength"].demand == pytest.approx(85.0) and by_axis["water_activity"].demand < 0.02          # the lowest water activity any cell has multiplied at is 0.6


def test_a_mild_acid_is_inside_the_table_for_an_acid_class(monkeypatch):
    _acid(monkeypatch, 3.0, 295.0)
    need, = _chemical("5b")
    assert (need.axis, need.grade, need.rows) == ("acid_strength", 2, ("A1",))          # pH 0 growth, 4.7 percent
    _acid(monkeypatch, 10.0, 295.0)
    heavy = {need.axis: need for need in _chemical("5b")}["acid_strength"]
    assert heavy.grade == 3 and heavy.rows == ("A1",)


def test_a_body_with_no_acid_gives_an_acid_class_no_record(monkeypatch):
    monkeypatch.setattr(acid_medium, "medium", lambda values, cloud_only=False: acid_medium.AcidMedium())
    assert _chemical("5a") == [] and _chemical("9") == []


def test_dense_co2_and_a_hydrocarbon_sea_have_no_water_and_no_row_reaches_that():
    for class_id, temperature in (("6", 300.0), ("7", 100.0)):
        found = {need.axis: need for need in _chemical(class_id, surface_temperature_k=temperature)}
        assert found["water_activity"].demand == 0.0 and found["water_activity"].grade is None and found["water_activity"].habitat == ar.CHEMICAL_HABITAT
    assert "t_cold" in {need.axis for need in _chemical("7", surface_temperature_k=100.0)}


def test_an_ammonia_water_life_needs_an_ammonia_tolerance_the_measured_life_does_not_have():
    found = {need.axis: need for need in _chemical("8", surface_temperature_k=250.0, ammonia_water_mass_fraction=0.1)}
    assert found["ammonia_fraction"].grade is None and found["ammonia_fraction"].demand == pytest.approx(0.1)          # the gate asks 10 percent, the measured limit is 0.43
    assert found["t_cold"].grade is None and found["t_cold"].habitat == ar.CHEMICAL_HABITAT
    mild = {need.axis: need for need in _chemical("8", surface_temperature_k=290.0, ammonia_water_mass_fraction=0.003)}
    assert (mild["ammonia_fraction"].grade, mild["ammonia_fraction"].rows) == (3, ("N1",))
    assert _chemical("8", surface_temperature_k=290.0, ammonia_water_mass_fraction=0.0005) == []
