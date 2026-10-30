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


def test_a_tuned_class_is_marked_speculative_and_a_small_pressure_is_given_in_pascals():
    history = {**HISTORY, "speculative": {"tier": "T7", "tweak": "L4"}}
    history["gases"] = [{"formula": "CH4", "name": "Methane", "steps": [[0.0, 0.0], [0.5, 0.0035]]}]
    view = life_air_chart.chart(history, "6")
    assert "Speculative chemistry" in view["note"] and "T7" in view["note"] and "L4" in view["note"]
    assert view["table"]["rows"][-1][1] == "3.5 Pa"


def test_the_age_axis_has_a_round_step_whatever_the_age():
    for top in (0.8, 3.0, 4.6, 9.5, 13.8):
        ticks = life_air_chart._axis_ticks(top)
        assert ticks[0] == 0 and 3 <= len(ticks) <= 8 and ticks[-1] <= top
