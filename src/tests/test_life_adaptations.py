# tests/test_life_adaptations.py

"""GEN.1725: the table of adaptations of life (`physics/data/life_adaptations.json`, `physics/life_adaptations.py`) and the traits cold and heat that read it (`physics/adaptation_traits.py`, `adaptation_draw.py`)."""

import copy
import json
import pathlib

import pytest

from planetgen import tuning
from planetgen.physics import adaptation_draw, adaptation_traits as traits, habitability as hab, life_adaptations as la, life_classes
from planetgen.util import draw

DATA = json.loads((pathlib.Path(la.__file__).parent / "data" / "life_adaptations.json").read_text(encoding="utf-8"))


def test_the_file_loads_and_every_row_has_what_a_row_needs():
    assert len(la.ADAPTATIONS) >= 30
    assert len({row.id for row in la.ADAPTATIONS}) == len(la.ADAPTATIONS)
    for row in la.ADAPTATIONS:
        assert row.tier in la.TIERS and row.tweak in la.TWEAKS and row.source_rows and row.name


def test_the_ordinary_limits_are_the_documented_ones():
    assert la.ordinary("t_cold") == pytest.approx(273.15)
    assert la.ordinary("t_hot") == pytest.approx(313.15)
    assert la.ordinary("light_floor") == pytest.approx(1.0)
    assert la.ordinary("water_hours") == 8760.0
    assert la.ordinary("ground_water_depth") == 1.0


def test_the_extreme_grade_reaches_the_measured_edges_of_oxygen_making_life():
    assert la.limit("t_cold", 3, class_id="1a") == pytest.approx(253.15)       # lichen, field, -20 C (GEN.1727)
    assert la.limit("t_hot", 3, class_id="1a") == pytest.approx(346.15)        # cyanobacteria, 73 C
    assert la.limit("light_floor", 3, class_id="1a") == pytest.approx(0.04)    # under the Arctic pack (GEN.1722 RF1)
    assert la.limit("water_hours", 3, class_id="1a") == pytest.approx(75.0)    # the hypolith dry limit (GEN.1722 RF3)
    assert la.limit("water_activity", 2, class_id="1a") == pytest.approx(0.75)
    assert la.limit("ground_water_depth", 3) == pytest.approx(68.0)


def test_a_grade_never_narrows_a_limit_and_a_wider_grade_never_falls_short_of_the_one_below():
    for axis, row in la.AXES.items():
        sign = -1.0 if row["wider"] == "lower" else 1.0
        values = [la.limit(axis, grade) for grade in (1, 2, 3)]
        assert values[0] == row["ordinary"]
        assert sign * (values[1] - values[0]) >= 0.0 and sign * (values[2] - values[1]) >= 0.0, axis


def test_a_young_life_does_not_have_what_a_late_stage_brings():
    assert la.limit("ground_water_depth", 3, life_age_gyr=1.0) == la.ordinary("ground_water_depth")      # rooted plants come 3.4 Gyr after the first life
    assert la.limit("ground_water_depth", 3, life_age_gyr=3.5) == pytest.approx(68.0)
    assert la.limit("t_cold", 3, class_id="1a", life_age_gyr=0.1) > la.limit("t_cold", 3, class_id="1a", life_age_gyr=1.0)   # the lichen needs land mats (0.6 Gyr)


def test_an_animal_stage_row_waits_for_animals():
    row = la.BY_ID["H8"]
    assert la.stage_delay_gyr(row.stage_needed) is None
    assert not la.available(row, 10.0) and la.available(row, 10.0, animals=True)


def test_a_dormant_row_gives_survival_not_activity():
    assert la.BY_ID["C15"].activity == "dormant"
    assert la.limit("t_cold", 3, life_types=["animal"], class_id=None, animals=True) > la.limit("t_cold", 3, life_types=["animal"], animals=True, dormant=True)


def test_a_class_that_makes_no_oxygen_borrows_nothing_from_a_row_of_the_oxygen_classes():
    assert la.limit("t_hot", 3, class_id="3") == la.ordinary("t_hot") or la.limit("t_hot", 3, class_id="3") < la.limit("t_hot", 3, class_id="1a")


def test_the_rows_of_the_chemical_classes_reach_the_measured_edges_and_no_further():
    """GEN.1810 and GEN.1811: the cold edge of a brine microbe, the hot record, the most acid a cell has been adapted to and the ammonia a cell has grown in."""
    assert la.limit("t_cold", 3, class_id="3a") == pytest.approx(248.15) and la.limit("t_hot", 3, class_id="4b") == pytest.approx(395.15)
    assert la.limit("acid_strength", 3, class_id="5a") == pytest.approx(11.0) and la.limit("acid_strength", 3, class_id="3a") == la.ordinary("acid_strength")
    assert la.limit("ammonia_fraction", 3, class_id="8") == pytest.approx(0.0043) and la.limit("ammonia_fraction", 3, class_id="5a") == la.ordinary("ammonia_fraction")
    assert la.limit("acid_strength", 3, class_id="5a") < 81.0 and la.limit("ammonia_fraction", 3, class_id="8") < 0.1 / 20          # the gates sit past them
    assert la.limit("water_activity", 3, class_id="5a") == pytest.approx(0.585)          # the lowest a cell has divided at; the acid's is 0.004 or less


def test_every_class_a_row_names_exists():
    for row in la.ADAPTATIONS:
        for class_id in row.classes:
            assert class_id == "all" or class_id in life_classes.BY_ID


def test_every_habitat_a_row_opens_is_listed_and_every_habitat_is_opened_or_open_water():
    opened = {row.opens_habitat for row in la.ADAPTATIONS if row.opens_habitat}
    assert opened <= set(la.HABITATS)
    assert set(la.HABITATS) - opened <= {"open_water"}
    for habitat in la.HABITATS.values():
        assert all(row_id in la.BY_ID for row_id in habitat["opened_by"])


def test_the_loader_refuses_a_bad_file(tmp_path):
    def broken(change):
        data = copy.deepcopy(DATA)
        change(data)
        path = tmp_path / "bad.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        return path

    with pytest.raises(la.AdaptationError):
        la._load(broken(lambda d: d["adaptations"][0].update(axis="t_sideways")))
    with pytest.raises(la.AdaptationError):
        la._load(broken(lambda d: d["adaptations"].append(dict(d["adaptations"][0]))))
    with pytest.raises(la.AdaptationError):
        la._load(broken(lambda d: d["adaptations"][0].update(stage_needed="the future")))
    with pytest.raises(la.AdaptationError):
        la._load(broken(lambda d: d["adaptations"][0].update(tier="T0")))
    with pytest.raises(la.AdaptationError):
        la._load(broken(lambda d: next(r for r in d["adaptations"] if r["id"] == "C17").update(classes=["10"])))      # not a life class
    cold = next(r for r in DATA["adaptations"] if r["id"] == "C8")
    with pytest.raises(la.AdaptationError):
        la._load(broken(lambda d: next(r for r in d["adaptations"] if r["id"] == "C8").update(grade2=290.0)))      # warmer than the ordinary limit is not wider
    assert cold["grade3"] < cold["grade2"]


# --- the traits ------------------------------------------------------------------------------------------------

def test_cold_and_heat_come_after_the_three_earlier_stresses():
    assert traits.STRESSES[:3] == ("radiation", "sulfide", "carbon_dioxide") and traits.STRESSES[3:] == ("cold", "heat")


def test_the_text_of_a_cold_and_a_heat_trait_round_trips():
    text = traits.text_of({"cold": 3, "heat": 2})
    assert text == "cold:3,heat:2"
    assert traits.parse(text) == {"cold": 3, "heat": 2}
    assert traits.sentence(text) == "cold (extreme), heat (tolerant)"


def test_the_trait_limit_reads_the_table():
    assert traits.limit(None, "cold") == la.ordinary("t_cold")
    assert traits.limit({"cold": 3}, "cold", class_id="1a") == pytest.approx(253.15)
    assert traits.limit({"heat": 3}, "heat", class_id="1a") == pytest.approx(346.15)
    with pytest.raises(ValueError):
        traits.limit({"radiation": 3}, "radiation")


def test_a_world_in_the_ordinary_range_draws_no_cold_or_heat_trait():
    for seed in range(300):
        found = adaptation_draw.draw({"surface_temperature_k": 288.0}, draw.Stream(seed))
        assert "cold" not in found and "heat" not in found


def test_a_cold_world_draws_cold_more_often_and_never_heat():
    cold = {"surface_temperature_k": 255.0}
    found = [adaptation_draw.draw(cold, draw.Stream(seed)) for seed in range(600)]
    share = sum(1 for f in found if "cold" in f) / len(found)
    assert 0.2 < share <= traits.P_CAP + 0.06
    assert all("heat" not in f for f in found)
    mild = [adaptation_draw.draw({"surface_temperature_k": 268.0}, draw.Stream(seed)) for seed in range(600)]
    assert sum(1 for f in mild if "cold" in f) < sum(1 for f in found if "cold" in f)


def test_a_hot_world_draws_heat_and_never_cold():
    hot = {"surface_temperature_k": 340.0}
    found = [adaptation_draw.draw(hot, draw.Stream(seed)) for seed in range(600)]
    assert sum(1 for f in found if "heat" in f) > 0 and all("cold" not in f for f in found)


def test_the_three_earlier_traits_of_a_body_are_what_they_were():
    values = {"surface_dose_msv_yr": 1.0e9, "gas_partial_pa": {"co2": 9.0e4, "h2s": 5.0e3}, "atmospheric_pressure_pa": 1.0e5, "surface_temperature_k": 255.0}
    plain = {k: v for k, v in values.items() if k != "surface_temperature_k"}
    for seed in range(40):
        with_cold = adaptation_draw.draw(values, draw.Stream(seed))
        without = adaptation_draw.draw(plain, draw.Stream(seed))
        assert {k: v for k, v in with_cold.items() if k in ("radiation", "sulfide", "carbon_dioxide")} == without


def test_the_score_of_a_shortfall_is_a_straight_line_from_zero_to_the_table_width():
    assert adaptation_draw.stress_score("cold", 0.0) == 0.0
    assert adaptation_draw.stress_score("cold", adaptation_draw.FULL["cold"] / 2) == pytest.approx(0.5)
    assert adaptation_draw.stress_score("heat", 1000.0) == 1.0
    assert adaptation_draw.FULL["cold"] == pytest.approx(20.0) and adaptation_draw.FULL["heat"] == pytest.approx(33.0)


def test_a_cold_or_hot_trait_never_changes_a_human_score():
    world = hab.World(pressure_kpa=101.3, temperature_c=-20, gases={"O2": 0.21, "N2": 0.78})
    tough = hab.World(pressure_kpa=101.3, temperature_c=-20, gases={"O2": 0.21, "N2": 0.78}, traits={"cold": 3, "heat": 3})
    assert hab.phi4(world) == hab.phi4(tough)
    assert hab.phi_bio(world) == hab.phi_bio(tough)      # nothing reads the new traits yet (GEN.1723 will)


def test_the_longest_trait_text_a_body_can_draw_fits_the_column():
    """`life_traits` is VARCHAR(48). Cold needs a mean temperature under the ordinary cold limit and heat one over the hot limit, so a body never has both."""
    assert la.ordinary("t_cold") < la.ordinary("t_hot")
    longest = traits.text_of({"radiation": 3, "sulfide": 3, "carbon_dioxide": 3, "heat": 3})
    assert len(longest) <= 48
    for temperature in (200.0, 250.0, 273.0, 288.0, 313.0, 340.0, 500.0):
        amounts = adaptation_draw.stresses({"surface_temperature_k": temperature})
        assert amounts["cold"] == 0.0 or amounts["heat"] == 0.0
