"""MAP.271: the look record `planetgen/web/maps/bodylook.py` builds for each body on the System Map."""

import json

from planetgen.web.maps import bodylook as bl


def _planet(**over):
    row = {"id": 7, "planet_class": "o", "body_type": "t", "atmosphere": "Nitrogen-Oxygen",
           "atmospheric_pressure_pa": 101_325.0, "surface_temperature_k": 288.0, "radius_km": 6371.0,
           "hydrosphere": "surface ocean", "ocean_fraction": 0.71, "land_fraction": 0.29,
           "ocean_class": "neutral", "ocean_depth_km": 3.7, "ice_shell_km": None}
    row.update(over)
    return row


def test_seed_is_stable_and_differs_by_kind_and_id():
    assert bl.seed_for("planet", 7) == bl.seed_for("planet", 7)
    assert bl.seed_for("planet", 7) != bl.seed_for("moon", 7)
    assert bl.seed_for("planet", 7) != bl.seed_for("planet", 8)
    assert 0 <= bl.seed_for("planet", 7) < 2 ** 32


def test_solid_planet_look_carries_the_stored_water_fractions():
    look = bl.planet_look(_planet())
    assert look["body"] == "solid" and look["cls"] == "O" and look["atm"] is True
    assert look["ocean"] == 0.71 and look["land"] == 0.29 and look["hydro"] == "surface ocean"
    assert look["oceanColor"] == bl.OCEAN_COLORS["neutral"]
    assert look["pkpa"] == 101.325 and look["tempK"] == 288.0
    assert look["quality"] == "estimate"
    assert "iceKm" not in look  # no stored value: left out, not zero


def test_airless_world_has_no_atmosphere_or_glow():
    look = bl.planet_look(_planet(atmosphere="None", atmospheric_pressure_pa=None))
    assert "atm" not in look and "glow" not in look and "pkpa" not in look


def test_gas_giant_look_has_no_surface_fields_and_a_temperature_colour():
    cold = bl.planet_look(_planet(body_type="g", surface_temperature_k=100.0))
    hot = bl.planet_look(_planet(body_type="g", surface_temperature_k=1600.0))
    assert cold["body"] == "gas" and "ocean" not in cold and "hydro" not in cold
    assert cold["color"] != hot["color"]


def test_land_colour_follows_temperature_bands():
    assert bl.land_color(100.0) == "#b8c4cc"
    assert bl.land_color(288.0) == "#a39070"
    assert bl.land_color(5000.0) == "#6a2c20"
    assert bl.land_color(None) == "#8a8f9c"


def test_star_look_uses_the_star_colour_and_solar_luminosity():
    star = {"id": 3, "star_type": "G2V", "temperature_k": 5778.0, "luminosity_w": 3.828e26, "radius_km": 696_000}
    look = bl.star_look(star)
    assert look["body"] == "star" and look["lumSolar"] == 1.0 and look["spec"] == "G2V"
    assert look["color"].startswith("#")


def test_look_json_is_compact_and_round_trips():
    look = bl.planet_look(_planet())
    text = bl.look_json(look)
    assert ", " not in text and ": " not in text and json.loads(text) == look
    assert bl.look_json(look) == text


def test_look_is_the_same_for_the_same_row():
    assert bl.look_json(bl.planet_look(_planet())) == bl.look_json(bl.planet_look(_planet()))


def _earth_air(**over):
    row = _planet(gravity_g=1.0, _star_temperature_k=5772.0, p_n2_kpa=79.1, p_o2_kpa=21.2, p_ar_kpa=0.95, p_co2_kpa=0.04)
    row.update(over)
    return row


def _channels(hex_colour):
    return [int(hex_colour[i:i + 2], 16) for i in (1, 3, 5)]


def test_the_glow_is_the_atmosphere_colour_of_the_stored_air_under_its_star():
    """MAP.244: the temperature gesture is gone; the limb colour follows the air and the star."""
    sun = bl.planet_look(_earth_air())["glow"]
    dwarf = bl.planet_look(_earth_air(_star_temperature_k=3300.0))["glow"]
    r, g, b = _channels(sun)
    assert b == 255 and b > g > r and sun != "#bfe3ff"
    assert _channels(dwarf)[2] < b                          # a red dwarf's air is not blue
    assert bl.planet_look(_earth_air(surface_temperature_k=400.0))["glow"] == sun        # temperature no longer decides it


def test_thick_carbon_dioxide_glows_orange_and_a_row_without_a_star_has_no_glow():
    venus = bl.planet_look(_earth_air(atmosphere="Carbon Dioxide", atmospheric_pressure_pa=9.2e6, gravity_g=0.9, p_co2_kpa=8880.0,
                                      p_n2_kpa=320.0, p_o2_kpa=0.0, p_ar_kpa=0.0))
    r, g, b = _channels(venus["glow"])
    assert r > g > b
    unknown = bl.planet_look(_earth_air(_star_temperature_k=None))
    assert unknown["atm"] is True and "glow" not in unknown


def test_a_giant_glows_with_the_colour_of_the_air_above_its_clouds():
    giant = bl.planet_look(_planet(body_type="g", planet_class="j", gravity_g=2.5, atmosphere="Hydrogen-Helium", atmospheric_pressure_pa=1e5,
                                   surface_temperature_k=130.0, _star_temperature_k=5772.0, p_h2_kpa=86.0))
    assert giant["body"] == "gas" and len(giant["glow"]) == 7
