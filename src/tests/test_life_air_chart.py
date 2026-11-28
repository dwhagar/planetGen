# tests/test_life_air_chart.py

"""UX.237: `web.lib.life_air_chart.chart` draws a Life and air history as an SVG with the stages marked, and as a table."""

import re

from planetgen.web.lib import life_air_chart

HISTORY = {
    "age_gyr": 4.0, "producers_gyr": 0.5, "rose_gyr": 1.0,
    "gases": [
        {"formula": "O2", "name": "Oxygen", "steps": [[0.0, 0.0], [0.5, 0.0], [1.0, 21.0]]},
        {"formula": "CO2", "name": "Carbon dioxide", "steps": [[0.0, 30.0], [0.5, 29.0], [1.0, 0.04]]},
        {"formula": "CH4", "name": "Methane", "steps": [[0.0, 0.0], [0.5, 1.0], [1.0, 0.0]]},
    ],
    "stages": [{"key": "abiogenesis", "age_gyr": 0.2}, {"key": "photosynthesis", "age_gyr": 0.5}, {"key": "multicellularity", "age_gyr": 3.0}],
    "speculative": None,
}


def test_there_is_nothing_to_draw_without_a_history_or_gases():
    assert life_air_chart.chart(None, "1a") is None
    assert life_air_chart.chart({**HISTORY, "gases": []}, "1a") is None


def test_the_chart_has_one_line_per_gas_a_mark_per_stage_and_says_what_it_shows():
    view = life_air_chart.chart(HISTORY, "1a", 2)
    svg = view["svg"]
    assert svg.startswith("<svg") and svg.endswith("</svg>") and 'role="img"' in svg and "<title" in svg and "<desc" in svg
    assert len(re.findall(r'<polyline class="life-air-line life-air-(o2|co2|ch4)"', svg)) == 3
    assert svg.count('class="life-air-stage"') == 3 and [stage["number"] for stage in view["stages"]] == [1, 2, 3]
    assert [stage["age"] for stage in view["stages"]] == ["0.20", "0.50", "3.00"]
    assert all(stage["name"] for stage in view["stages"])
    assert view["legend"].count("<li>") == 3 and view["note"] == ""


def test_the_chart_names_its_axes_with_units_and_says_in_plain_words_what_the_lines_are():
    """UX.322: the axes say what they hold and in what unit, the key names the vertical lines, and three plain sentences say how to read the chart."""
    view = life_air_chart.chart({**HISTORY, "events": [{"age_gyr": 2.5, "cause": "impact", "kill_share": 0.7}]}, "1a")
    assert "Time since the star formed (Gyr, billions of years)" in view["svg"] and "Pressure of the gas (kPa)" in view["svg"]
    assert len(view["explain"]) == 3 and all(sentence.endswith(".") for sentence in view["explain"])
    text = " ".join(view["explain"])
    assert "kilopascals (kPa)" in text and "billions of years (Gyr)" in text and "ten times" in text and "mass extinctions" in text and "stages" in text
    assert view["marks"].count("<li>") == 2 and "A stage the life reached" in view["marks"] and "A mass extinction" in view["marks"]
    plain = life_air_chart.chart({**HISTORY, "stages": []}, "1a")
    assert plain["marks"] == ""


def test_every_point_of_every_line_lies_inside_the_plot():
    view = life_air_chart.chart(HISTORY, "1a")
    for points in re.findall(r'points="([^"]+)"', view["svg"]):
        for pair in points.split():
            x, y = (float(part) for part in pair.split(","))
            assert life_air_chart.LEFT - 0.1 <= x <= life_air_chart.WIDTH - life_air_chart.RIGHT + 0.1, pair
            assert life_air_chart.TOP - 0.1 <= y <= life_air_chart.HEIGHT - life_air_chart.BOTTOM + 0.1, pair


def test_the_table_gives_each_span_of_ages_with_its_pressures_and_none_for_a_gas_not_there():
    view = life_air_chart.chart(HISTORY, "1a")
    assert view["table"]["header"] == ["Age (Gyr)", "Oxygen", "Carbon dioxide", "Methane"]
    assert view["table"]["rows"] == [["0.00 to 0.50", "none", "30 kPa", "none"], ["0.50 to 1.00", "none", "29 kPa", "1 kPa"],
                                     ["1.00 to 4.00", "21 kPa", "0.04 kPa", "none"]]


def test_a_less_sure_class_links_to_its_references_section_and_a_small_pressure_is_given_in_pascals():
    history = {**HISTORY, "speculative": {"tier": "T7", "tweak": "L4"}}
    history["gases"] = [{"formula": "CH4", "name": "Methane", "steps": [[0.0, 0.0], [0.5, 0.0035]]}]
    view = life_air_chart.chart(history, "6")
    assert "less sure" in view["note"] and 'href="/references/speculative#class-6"' in view["note"] and "T7" not in view["note"] and "L4" not in view["note"]
    assert view["table"]["rows"][-1][1] == "3.5 Pa"


def test_the_age_axis_has_a_round_step_whatever_the_age():
    for top in (0.8, 3.0, 4.6, 9.5, 13.8):
        ticks = life_air_chart._axis_ticks(top)
        assert ticks[0] == 0 and 3 <= len(ticks) <= 8 and ticks[-1] <= top


def test_a_world_with_several_classes_says_what_each_made_of_each_gas_and_a_world_with_one_says_nothing():
    several = {**HISTORY, "classes": [{"life_class": 1, "variant": "a", "share": 0.7, "producers_gyr": 0.5}, {"life_class": 3, "variant": None, "share": 0.3, "producers_gyr": 0.8}],
               "gases": [{**HISTORY["gases"][0], "by_class": [{"life_class": 1, "variant": "a", "share": 1.0}]},
                         {**HISTORY["gases"][2], "by_class": [{"life_class": 3, "variant": None, "share": 0.25}, {"life_class": 1, "variant": "a", "share": 0.75}]}]}
    assert life_air_chart.chart(several, "1a")["makers"] == [{"gas": "Methane", "text": "Class 1a 75 percent, Class 3 25 percent"}]
    assert life_air_chart.chart({**several, "classes": several["classes"][:1]}, "1a")["makers"] == []
    assert life_air_chart.chart(HISTORY, "1a")["makers"] == []


def test_the_page_lists_only_the_mass_extinctions_and_the_ages_name_the_one_that_ended_them():
    """UX.321: a shock that killed less than the share that ends an age of life (60 percent) is no mass extinction: the list and the chart's dashed lines hold the ones that killed that much or more."""
    events = [{"age_gyr": 1.0, "cause": "impact", "kill_share": 0.12}, {"age_gyr": 2.5, "cause": "impact", "kill_share": 0.7}, {"age_gyr": 2.8, "cause": "flood_basalt", "kill_share": 0.45}, {"age_gyr": 3.0, "cause": "flood_basalt", "kill_share": 0.8}]
    ages = [{"number": 1, "start_gyr": 0.2, "end_gyr": 2.5, "cause": "impact", "kill_share": 0.7}, {"number": 2, "start_gyr": 2.5, "end_gyr": None, "cause": None, "kill_share": 0.0}]
    view = life_air_chart.chart({**HISTORY, "events": events, "life_ages": ages}, "1a")
    assert [event["age"] for event in view["extinctions"]] == ["2.50", "3.00"] and view["svg"].count('class="life-air-event"') == 2
    assert view["extinction_floor"] == "60 percent"
    assert view["life_ages"][0]["ended"] == "the mass extinction at 2.50 Gyr" and view["life_ages"][1]["ended"] == ""
