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
    assert "needs a gas-mix mask and a cold suit." in text and "microbial life 3 of 9" in text


def test_a_lifeless_airless_world_says_so():
    row = {**ROW, "atmospheric_pressure_pa": 0, "life_microbial": 0, "phi4": 0.1, "tier_suitability": 3, "equipment": {"required": []}}
    text = " ".join(planet_summary.summary_sentences("Rock", row, {**VALUES, "atmospheric_pressure_pa": "0 Pa"}))
    assert "it has no atmosphere" in text and "No life has been found here." in text
    assert "needs no safety equipment" in text and "so people can live here" not in text


def test_a_body_with_no_class_or_scores_still_gets_a_sentence():
    assert planet_summary.summary_sentences("Bare", {}, {}) == ["Bare is a world."]


def test_the_descriptions_blocks_are_sorted_into_their_sections():
    blocks = [model.facts("Planet", []), model.prose("description", "It is temperate."), model.composition("Atmosphere", [], "air"),
              model.ratings({"name": "PHI-6", "score": 0.8, "tier": "Green"}, [], []), model.facts("Radiation", []),
              model.facts("Survival equipment", []), model.prose("life", "Life lives."), model.prose("timeline", "Old."),
              model.prose("flavour", "Sensors show.")]
    parts = planet_layout.split_blocks(blocks)
    assert [b["type"] for b in parts["phi6"]] == ["ratings", "facts", "facts"]
    assert [b["role"] for b in parts["life"]] == ["life", "timeline"] and [b["role"] for b in parts["rest"]] == ["flavour"]


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
