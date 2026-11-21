# tests/test_stored_air.py

"""API.51: the life-air answer made from a stored air history (`physics/stored_air.py`, the record of `physics/air_record.py`), and the Life and air chart that draws every gas of it."""

import re

from planetgen.physics import air_events, air_record, life_ages, stored_air
from planetgen.web.lib import life_air_chart


def _record():
    steps = (
        air_record.RecordStep(0.01, 1800.0, 9.0e6, {"H2O": 6.0e6, "CO2": 3.0e6}, None, None),
        air_record.RecordStep(0.5, 288.0, 1.1e5, {"N2": 7.0e4, "CO2": 3.5e4, "CH4": 2.0}, "reduced", None),
        air_record.RecordStep(2.0, 290.0, 1.0e5, {"N2": 7.9e4, "O2": 2.1e4, "CO2": 40.0}, "plateau", None),
        air_record.RecordStep(3.0, 289.0, 9.9e4, {"N2": 7.9e4, "O2": 2.0e4, "CO2": 30.0, "Ar": 0.5}, "plateau", (air_events.FLOOD_BASALT, 0.4)),
    )
    ages = (life_ages.LifeAge(1, 0.5, 3.0, air_events.FLOOD_BASALT, 0.4), life_ages.LifeAge(2, 3.2, None, None, None))
    return air_record.Record(steps, ages)


def _history(rebuilt=None):
    text = air_record.pack(_record())
    return stored_air.from_text(text, 4.0, {"abiogenesis": 0.4, "photosynthesis": 0.5, "multicellularity": 9.0}, rebuilt)


def test_every_gas_of_the_record_is_in_the_answer_the_biggest_first_with_zero_where_it_is_under_the_floor():
    history = _history()
    assert [gas["formula"] for gas in history["gases"]] == ["H2O", "CO2", "N2", "O2", "CH4", "Ar"]
    by = {gas["formula"]: gas for gas in history["gases"]}
    assert by["O2"]["steps"] == [[0.01, 0.0], [0.5, 0.0], [2.0, 21.0], [3.0, 20.0]]
    assert by["H2O"]["steps"][0] == [0.01, 6000.0] and by["H2O"]["steps"][-1] == [3.0, 0.0]
    assert by["Ar"]["name"] == "Argon" and by["CO2"]["name"] == "Carbon dioxide"
    assert history["source"] == "stored" and history["age_gyr"] == 4.0


def test_the_temperature_the_pressure_the_extinctions_and_the_life_ages_come_with_it():
    history = _history()
    assert history["temperature_k"] == [[0.01, 1800.0], [0.5, 288.0], [2.0, 290.0], [3.0, 289.0]]
    assert history["pressure_kpa"][-1] == [3.0, 99.0]
    assert history["events"] == [{"age_gyr": 3.0, "cause": air_events.FLOOD_BASALT, "kill_share": 0.4}]
    assert history["life_ages"] == [{"number": 1, "start_gyr": 0.5, "end_gyr": 3.0, "cause": air_events.FLOOD_BASALT, "kill_share": 0.4},
                                    {"number": 2, "start_gyr": 3.2, "end_gyr": None, "cause": None, "kill_share": None}]


def test_the_stages_the_star_is_old_enough_for_are_marked_and_the_rebuilt_answer_supplies_the_classes_and_what_each_made():
    rebuilt = {"producers_gyr": 0.5, "rose_gyr": 1.8, "classes": [{"life_class": 1, "variant": "a", "share": 1.0, "producers_gyr": 0.5}], "speculative": {"tier": "T7", "tweak": "L4"},
               "gases": [{"formula": "CH4", "name": "Methane", "steps": [[0.0, 0.0]], "by_class": [{"life_class": 1, "variant": "a", "share": 1.0}]}]}
    history = _history(rebuilt)
    assert [stage["key"] for stage in history["stages"]] == ["abiogenesis", "photosynthesis"], "a stage the star is too young for is left out"
    assert history["producers_gyr"] == 0.5 and history["rose_gyr"] == 1.8 and history["classes"] == rebuilt["classes"] and history["speculative"] == rebuilt["speculative"]
    by = {gas["formula"]: gas for gas in history["gases"]}
    assert by["CH4"]["by_class"] == rebuilt["gases"][0]["by_class"] and by["O2"]["by_class"] == []
    plain = _history()
    assert plain["classes"] == [] and plain["producers_gyr"] is None and plain["rose_gyr"] is None and plain["speculative"] is None


def test_no_text_is_no_history_and_an_unnamed_gas_keeps_its_formula(monkeypatch):
    assert stored_air.from_text(None, 4.0, {}) is None and stored_air.from_text("", 4.0, {}) is None
    monkeypatch.setattr(stored_air.gases, "has", lambda formula: False)
    assert stored_air._name("XeF2") == "XeF2"


def test_the_chart_draws_a_line_for_every_gas_each_other_gas_in_a_style_of_its_own_and_marks_the_extinctions():
    view = life_air_chart.chart(_history(), "1a")
    svg = view["svg"]
    classes = re.findall(r'<polyline class="life-air-line ([a-z0-9-]+)"', svg)
    assert len(classes) == 6
    others = [name for name in classes if name.startswith("life-air-g")]
    assert others == ["life-air-g0", "life-air-g1", "life-air-g2"], "H2O, N2 and Ar take the other styles in the order of the history"
    assert {"life-air-o2", "life-air-co2", "life-air-ch4"} <= set(classes)
    assert svg.count('class="life-air-event"') == 1
    assert view["stored"] is True
    assert view["extinctions"] == [{"age": "3.00", "cause": "Flood basalt eruptions", "killed": "40 percent"}]
    assert view["life_ages"][0] == {"number": 1, "span": "0.50 to 3.00", "ended": "Flood basalt eruptions, 40 percent of species"}
    assert view["life_ages"][1] == {"number": 2, "span": "3.20 to now", "ended": ""}
    assert view["table"]["header"][0] == "Age (Gyr)" and len(view["table"]["header"]) == 7


def test_a_rebuilt_history_draws_as_it_did_with_no_extinctions_and_the_note_it_had():
    plain = {"age_gyr": 4.0, "producers_gyr": 0.5, "rose_gyr": 1.0, "stages": [], "speculative": None,
             "gases": [{"formula": "O2", "name": "Oxygen", "steps": [[0.0, 0.0], [1.0, 21.0]]}]}
    view = life_air_chart.chart(plain, "1a")
    assert view["stored"] is False and view["extinctions"] == [] and view["life_ages"] == [] and 'class="life-air-event"' not in view["svg"]
