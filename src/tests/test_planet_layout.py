# tests/test_planet_layout.py

"""
UX.217 to UX.222: the pieces of the redesigned planet page: the zone's name and range, the written summary, the sorted
description blocks, the air by mole fraction, the body by mass and the civilization paragraph.
"""

import pytest

from planetgen.describe import model, sheet
from planetgen.web.lib import planet_layout, planet_summary


def test_a_zone_letter_reads_as_its_name_and_the_range_it_covers():
    specs = sheet.specs_for("planet")
    stored = {"habitable_zone_inner_km": 1.4e8, "habitable_zone_outer_km": 2.0e8}
    ecosphere = sheet._zone_value("E", stored, specs, "planets")
    assert ecosphere[0] == "Ecosphere / Habitable Zone, " and ecosphere[2] == " to " and ecosphere[-1] == " from its star"
    assert ecosphere[1]["qty"]["value"] == 1.4e8 and ecosphere[3]["qty"]["value"] == 2.0e8
    hot, cold = (sheet._zone_value(letter, stored, specs, "planets") for letter in ("h", "c"))
    assert hot[0] == "Hot Zone, closer than " and hot[1]["qty"]["value"] == 1.4e8
    assert cold[0] == "Cold Zone, farther than " and cold[1]["qty"]["value"] == 2.0e8
    assert sheet._zone_value("e", {}, specs, "planets") == ["Ecosphere / Habitable Zone"]
    assert sheet._zone_value("x", stored, specs, "planets") is None


ROW = {"planet_class": "M", "body_type": "t", "zone": "e", "atmospheric_pressure_pa": 101325.0, "phi4": 0.8,
       "tier_suitability": 1, "life_microbial": 3, "life_flora": 0, "life_fauna": 0, "life_sentient": 0,
       "equipment": {"required": ["a gas-mix mask", "a cold suit"]}}
VALUES = {"planet_class": '<a href="/classes/planet/M">M</a>', "surface_temperature_k": "288 K",
          "atmospheric_pressure_pa": "101 kPa", "gravity_g": "1 g"}


def test_the_summary_is_written_from_the_data():
    text = " ".join(planet_summary.summary_sentences("Kepler b", ROW, VALUES))
    assert text.startswith('Kepler b is a Class <a href="/classes/planet/M">M</a> world in the Ecosphere / Habitable Zone.')
    assert "Its surface is 288 K, the air at the surface is 101 kPa, and gravity is 1 g." in text
    assert "PHI-6 habitability rating is Green, so people can live here, and its Suitability score is Green" in text
    assert "needs a gas-mix mask and a cold suit." in text
    assert "life" not in text.lower().replace("people can live", "")       # the life of the world is its own paragraph (UX.251)


def test_a_lifeless_airless_world_says_so():
    row = {**ROW, "atmospheric_pressure_pa": 0, "life_microbial": 0, "phi4": 0.1, "tier_suitability": 3, "equipment": {"required": []}}
    text = " ".join(planet_summary.summary_sentences("Rock", row, {**VALUES, "atmospheric_pressure_pa": "0 Pa"}))
    assert "it has no atmosphere" in text and "No life" not in text
    assert "needs no safety equipment" in text and "so people can live here" not in text


def test_a_body_with_no_class_or_scores_still_gets_a_sentence():
    assert planet_summary.summary_sentences("Bare", {}, {}) == ["Bare is a world."]


def test_the_descriptions_blocks_are_sorted_into_their_sections():
    blocks = [model.facts("Planet", []), model.prose("description", "It is temperate."), model.composition("Atmosphere", [], "air"),
              model.ratings({"name": "PHI-6", "score": 0.8, "tier": "Green"}, [], []), model.facts("Radiation", []),
              model.facts("Survival equipment", []), model.prose("life-class", "Class 1a life lives here."), model.prose("life", "Life looks green."),
              model.prose("timeline", "Old."), model.prose("milestone", "Latest."), model.prose("flavour", "Sensors show.")]
    parts = planet_layout.split_blocks(blocks)
    assert [b["type"] for b in parts["phi6"]] == ["ratings", "facts", "facts"]
    assert [b["role"] for b in parts["life"]] == ["life"] and [b["role"] for b in parts["rest"]] == ["flavour"]


def _row(key, label, raw, value=None):
    return {"key": key, "label": label, "value": value or [str(raw)], "text": str(raw), "raw": raw, "stored": True}


def test_the_air_is_a_table_of_mole_fractions_with_partial_pressures():
    groups = [{"key": "atmosphere", "label": "Atmosphere", "rows": [
        _row("atmospheric_pressure_pa", "Surface pressure", 101325.0), _row("p_N2", "Nitrogen partial pressure", 79000.0),
        _row("p_O2", "Oxygen partial pressure", 21000.0), _row("p_Ar", "Argon partial pressure", 12.0)]}]
    table = planet_layout.gas_table({"groups": groups})
    assert [gas["name"] for gas in table] == ["Nitrogen", "Oxygen", "Argon"]
    assert table[0]["fraction"] == "79 percent" and table[1]["fraction"] == "21 percent"
    assert table[2]["fraction"] == "0.012 percent"
    assert planet_layout.gas_table({"groups": []}) == []


def test_the_body_is_a_table_of_its_mass_fractions():
    groups = [{"key": "composition", "label": "Composition", "rows": [
        _row("metal_mass_fraction", "Metal", 0.3), _row("rock_mass_fraction", "Rock", 0.69), _row("ice_mass_fraction", "Ice", 0.0),
        _row("envelope_mass_fraction", "Envelope", None), _row("porosity", "Porosity", 0.1)]}]
    table = planet_layout.mass_table({"groups": groups})
    assert [part["name"] for part in table] == ["Metal (iron and nickel core)", "Rock (silicate mantle and crust)"]
    assert table[0]["share"] == "30 percent of its mass"
    assert planet_layout.MASS_SHARED_KEYS >= {"metal_mass_fraction", "envelope_mass_fraction"}


def test_the_civilization_paragraph_exists_only_when_there_is_something_to_say():
    assert planet_layout.civilization_sentences({"inhabited": False}, None) == []
    assert planet_layout.civilization_sentences(None, None) == []
    both = planet_layout.civilization_sentences({"inhabited": True}, {"name": "Velarix <b>", "url": "/species/3"})
    assert len(both) == 2 and 'href="/species/3"' in both[0] and "Velarix &lt;b&gt;" in both[0]


EARTH_DAY = {**VALUES, "day_share_mean": "50.1 %", "day_share_equator": "50.2 %", "solar_day_hours": "24 h", "daylight_hours": "12 h"}


def test_the_summary_gives_the_time_the_star_is_up_in_one_sentence():
    """UX.205: the equator's share, the hours of daylight in the solar day and the share over the whole surface."""
    text = planet_summary.summary_sentences("Kepler b", ROW, EARTH_DAY)
    assert ("Its star is up 50.2 % of the time at the equator, 12 h of daylight in each solar day of 24 h, and 50.1 % "
            "averaged over the whole surface.") in text


def test_a_tidally_locked_world_says_it_has_a_day_side_a_terminator_and_a_night_side():
    locked = {**VALUES, "day_share_mean": "50 %"}      # a mean share but no solar day: it turns once an orbit
    text = " ".join(planet_summary.summary_sentences("Kepler b", ROW, locked))
    assert "keeps one face to its star" in text and "day side" in text and "terminator" in text and "night side" in text
    assert "solar day" not in text


def test_a_moon_in_its_planets_shadow_gets_the_extra_night_and_the_longest_eclipse():
    values = {**EARTH_DAY, "parent_shadow_share": "4.6 %", "longest_eclipse_hours": "2.3 h"}
    text = " ".join(planet_summary.summary_sentences("Io", ROW, values, {"parent_shadow_share": 0.046}))
    assert "Its planet's shadow adds a night: it is in the shadow 4.6 % of the time, the longest eclipse lasting 2.3 h." in text
    assert "shadow" not in " ".join(planet_summary.summary_sentences("Io", ROW, values, {"parent_shadow_share": 0.0}))


def test_a_body_with_no_day_and_night_says_nothing_of_them():
    """A star, or a free body with no rotation stored, has no shares: nothing is said."""
    text = " ".join(planet_summary.summary_sentences("Bare", ROW, VALUES))
    assert "star is up" not in text and "one face" not in text and "shadow" not in text


# --- the life class lines (UX.237) ----------------------------------------------------------------------------------

LIVING = {"life_class": 1, "life_class_variant": "a", "life_class_id": "1a", "life_class_best": "1a", "life_class_best_score": 7,
          "life_class_scores": {"1a": 7, "3a": 2, "2": 4}, "life_stage": "multicellularity", "life_animal_grade": 2}


def test_a_world_with_life_says_what_it_supports_what_else_and_the_stage_it_has_reached():
    view = planet_layout.life_class_view(LIVING)
    assert view["supports"] == "Best suited to Class 1a life (oxygenic light-eaters), suitability 7 of 9."
    assert [text.split(" (")[0] for text in view["also"]] == ["Class 2", "Class 3a"]      # best score first
    assert view["also"][0].endswith(", suitability 4 of 9)") and view["also"][1].endswith(", suitability 2 of 9)")
    assert view["has"] == "Life here is Class 1a, at the stage of small mobile animals (animal grade A2)."
    assert view["words"]["class"] == "1a" and view["words"]["lives"]
    assert [label for label, _text in view["words"]["rows"]] == ["Drinks", "Eats", "Breathes", "Gives off"]
    assert dict(view["words"]["rows"])["Drinks"] == "water"


def test_a_world_that_supports_a_class_but_has_none_still_says_what_that_class_lives_on():
    row = {**LIVING, "life_class": None, "life_class_id": None, "life_stage": None, "life_class_scores": {"1a": 7}}
    view = planet_layout.life_class_view(row)
    assert view["supports"].startswith("Best suited to Class 1a life") and view["also"] == [] and view["has"] is None
    assert view["words"]["class"] == "1a" and not view["words"]["lives"]


def test_a_stage_without_animals_names_no_grade_and_a_world_without_classes_has_no_lines():
    view = planet_layout.life_class_view({**LIVING, "life_stage": "photosynthesis", "life_animal_grade": 0})
    assert view["has"] == "Life here is Class 1a, at the stage of oxygenating mats."
    assert planet_layout.life_class_view({"life_class_scores": {}, "life_class_id": None}) is None
    assert planet_layout.life_class_view(None) is None


READINGS = [
    {"kind": "microbial", "label": "Microbes", "score": 8, "abundance": "Abundant", "trend": "Stable", "security": "Secure",
     "years_left_gyr": None, "text": "Microbes: abundant, stable, secure."},
    {"kind": "flora", "label": "Flora", "score": 3, "abundance": "Occasional", "trend": "Declining", "security": "Precarious",
     "years_left_gyr": 0.4, "text": "Flora: occasional, declining, precarious."},
    {"kind": "fauna", "label": "Fauna", "score": 0, "abundance": "Absent", "trend": None, "security": None,
     "years_left_gyr": None, "text": "Fauna: absent."},
    {"kind": "sentient", "label": "Sentient life", "score": 0, "abundance": "Absent", "trend": None, "security": None,
     "years_left_gyr": None, "text": "Sentient life: absent."},
]


def test_the_life_paragraph_leads_with_what_lives_then_the_classes_that_suit():
    view = planet_layout.life_class_view({**LIVING, "life_microbial": 8, "life_flora": 3}, READINGS)
    text = view["paragraph"]
    assert text.startswith("Life here is Class 1a, at the stage of small mobile animals (animal grade A2). Conditions here suit Class 1a life")
    assert text.endswith("best, suitability 7 of 9; they also suit " + ", ".join(view["also"]) + ".")
    assert view["also"][0].startswith("Class 2 (") and view["also"][1].startswith("Class 3a (")
    assert "By kind" not in text
    bare = planet_layout.life_class_view({**LIVING, "life_class_id": None, "life_stage": None})["paragraph"]
    assert bare.startswith("Nothing lives here yet.")


def test_at_a_glance_the_life_paragraph_gives_the_words_of_each_kind_that_is_there_then_the_classes_that_suit():
    glance = planet_layout.life_class_view(LIVING, READINGS)["glance"]
    assert glance == ("Microbes: abundant, stable, secure. Flora: occasional, declining, precarious. Conditions suit Class 1a (oxygenic "
                      "light-eaters) best, then Class 2 (" + planet_layout._short_name("2") + ") and Class 3a (" + planet_layout._short_name("3a") + ").")
    assert "of 9" not in glance and "Fauna" not in glance


def test_the_glance_names_two_other_classes_and_says_others_for_the_rest():
    scores = {"1a": 7, "2": 6, "3a": 5, "3b": 4, "4a": 3}
    glance = planet_layout.life_class_view({**LIVING, "life_class_scores": scores}, READINGS)["glance"]
    assert glance.endswith(", then Class 2 (" + planet_layout._short_name("2") + "), Class 3a (" + planet_layout._short_name("3a")
                           + ") and others.")


def test_a_world_that_suits_a_class_but_has_no_life_says_nothing_lives_there_yet():
    row = {**LIVING, "life_class": None, "life_class_id": None, "life_stage": None, "life_class_scores": {"1a": 7}}
    view = planet_layout.life_class_view(row)
    assert view["glance"] == "Nothing lives here yet. Conditions suit Class 1a (oxygenic light-eaters) best."
    assert view["readings"] == []


def test_a_world_with_no_life_and_no_class_that_suits_has_no_life_paragraph():
    row = {"life_class_scores": {}, "life_class_id": None, "life_class_best": None, "life_class_best_score": 0}
    assert planet_layout.life_class_view(row) is None


def test_the_life_section_has_a_row_for_each_kind_and_the_stored_facts():
    row = {**LIVING, "evolutionary_pace": "normal", "life_system_age_gyr": 4.6, "life_milestone_age_gyr": 3.6}
    view = planet_layout.life_class_view(row, READINGS)
    assert [r["kind"] for r in view["readings"]] == ["Microbes", "Flora", "Fauna", "Sentient life"]
    assert view["readings"][0] == {"kind": "Microbes", "abundance": "Abundant (8 of 9)", "trend": "Stable", "security": "Secure", "left": ""}
    assert view["readings"][1]["left"] == "about 0.4 billion years"
    assert view["readings"][2] == {"kind": "Fauna", "abundance": "Absent", "trend": "", "security": "", "left": ""}
    assert dict(view["facts"]) == {"Stage": "small mobile animals", "Animal grade": "A2", "Pace of evolution": "standard",
                                   "Age of the system": "4.60 billion years", "Stage reached at": "3.60 billion years"}
