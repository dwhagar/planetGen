"""MAP.271: the look record `planetgen/web/maps/bodylook.py` builds for each body on the System Map."""

import json

from planetgen.web.maps import bodylook as bl

from tests.air_support import fold_air


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
    return fold_air(row)


def _channels(hex_colour):
    return [int(hex_colour[i:i + 2], 16) for i in (1, 3, 5)]


def test_the_glow_is_the_atmosphere_colour_of_the_stored_air_under_its_star():
    """MAP.244: the temperature gesture is gone; the limb colour follows the air and the star."""
    sun = bl.planet_look(_earth_air())["glow"]
    dwarf = bl.planet_look(_earth_air(_star_temperature_k=3300.0))["glow"]
    r, g, b = _channels(sun)
    assert b == 255 and b > g > r and sun != "#bfe3ff"
    assert _channels(dwarf)[2] - _channels(dwarf)[0] < 0.5 * (b - r)  # a red dwarf's air is whiter than the Sun's (GEN.971: its model atmosphere)
    assert bl.planet_look(_earth_air(surface_temperature_k=400.0))["glow"] == sun        # temperature no longer decides it


def test_thick_carbon_dioxide_glows_orange_and_a_row_without_a_star_has_no_glow():
    venus = bl.planet_look(_earth_air(atmosphere="Carbon Dioxide", atmospheric_pressure_pa=9.2e6, gravity_g=0.9, p_co2_kpa=8880.0,
                                      p_n2_kpa=320.0, p_o2_kpa=0.0, p_ar_kpa=0.0))
    r, g, b = _channels(venus["glow"])
    assert r > g > b
    unknown = bl.planet_look(_earth_air(_star_temperature_k=None))
    assert unknown["atm"] is True and "glow" not in unknown


def test_a_giant_glows_with_the_colour_of_the_air_above_its_clouds():
    giant = bl.planet_look(fold_air(_planet(body_type="g", planet_class="j", gravity_g=2.5, atmosphere="Hydrogen-Helium", atmospheric_pressure_pa=1e5,
                                            surface_temperature_k=130.0, _star_temperature_k=5772.0, p_h2_kpa=86.0)))
    assert giant["body"] == "gas" and len(giant["glow"]) == 7


def test_star_colour_is_the_model_atmosphere_colour_of_its_temperature_and_gravity():
    from planetgen.physics import stellar_library
    star = {"id": 3, "star_type": "M3V", "yerkes_class": "V", "temperature_k": 3400.0, "luminosity_w": 1e24}
    look = bl.star_look(star)
    assert look["color"] == stellar_library.star_colour(3400.0).hex
    dense = bl.star_look({**star, "mass_kg": 4.0e29, "radius_km": 3.0e5})
    assert dense["color"] == stellar_library.star_colour(3400.0, stellar_library.surface_logg(4.0e29, 3.0e5)).hex
    assert look["yerkes"] == "V"
    hot = bl.star_look({**star, "temperature_k": 20000.0})
    assert hot["color"] != look["color"]


def test_star_without_a_temperature_falls_back_to_the_spectral_colour():
    star = {"id": 4, "star_type": "K2V", "temperature_k": None, "luminosity_w": 1e26}
    assert bl.star_look(star)["color"].startswith("#")


# --- Phenomena (MAP.278) -----------------------------------------------------

def test_asteroid_grain_colour_follows_the_composition_family():
    assert bl.asteroid_grain_color("B3") == bl.ASTEROID_GRAIN_COLORS["carbonaceous"]
    assert bl.asteroid_grain_color("E2") == bl.ASTEROID_GRAIN_COLORS["stony"]
    assert bl.asteroid_grain_color("H1") == bl.ASTEROID_GRAIN_COLORS["metallic"]
    assert bl.asteroid_grain_color("L4") == bl.ASTEROID_GRAIN_COLORS["icy"]
    assert bl.asteroid_grain_color(None) == "#b89a6e"
    assert bl.asteroid_grain_color("Z9") == "#b89a6e"


def test_rogue_planet_look_is_icy_when_cold_dry_when_warm_and_has_no_sea():
    cold = bl.rogue_planet_look({"id": 1, "planet_type": "t", "surface_temperature_k": 60.0, "radius_km": 5000})
    warm = bl.rogue_planet_look({"id": 1, "planet_type": "t", "surface_temperature_k": 300.0, "radius_km": 5000})
    assert cold["hydro"] == "ice" and warm["hydro"] == "dry"
    assert "ocean" not in cold and "glow" not in cold
    assert cold["kind"] == "rogue_planet" and cold["seed"] == bl.seed_for("rogue_planet", 1)
    buried = bl.rogue_planet_look({"id": 2, "planet_type": "t", "surface_temperature_k": 80.0,
                                   "has_liquid_water": 1, "ice_shell_thickness_km": 30.0})
    assert buried["hydro"] == "ice-covered ocean"


def test_rogue_giant_uses_the_giant_colour_and_falls_back_to_effective_temperature():
    giant = bl.rogue_planet_look({"id": 3, "planet_type": "g", "surface_temperature_k": None,
                                  "effective_temperature_k": 100.0})
    assert giant["body"] == "gas" and giant["color"] == bl.giant_color(100.0)


def test_comet_nucleus_is_pale_when_active_and_dark_when_not():
    assert bl.comet_nucleus_color({"is_active": 1}) != bl.comet_nucleus_color({"is_active": 0})


def test_blackbody_hex_rejects_unusable_temperatures():
    assert bl.blackbody_hex(None) is None and bl.blackbody_hex(0) is None
    assert bl.blackbody_hex(1.5e7).startswith("#")


def test_phenomenon_view_carries_the_painters_data():
    from planetgen.web.maps import phenomenonrender as pr
    params, _ = pr._params("rogue_planet", {"id": 9, "planet_type": "t", "surface_temperature_k": 50.0})
    assert params["look"]["body"] == "solid" and params["seed"] == params["look"]["seed"]
    params, _ = pr._params("interstellar_comet", {"id": 4, "is_active": 1})
    assert params["nucleus"] == bl.comet_nucleus_color({"is_active": 1}) and "seed" in params
    params, _ = pr._params("black_hole", {"id": 5, "has_accretion_disk": 1, "temperature_k": 1e7, "mass_class": "stellar"})
    assert params["diskColor"].startswith("#")
    params, _ = pr._params("black_hole", {"id": 5, "has_accretion_disk": 0, "temperature_k": 1e7})
    assert params["diskColor"] is None
    params, _ = pr._params("neutron_star", {"id": 6, "surface_temperature_k": 1e6, "spin_period_ms": 33})
    assert params["surface"].startswith("#")


# --- MAP.328: the spin a body is drawn turning with ---------------------------------

def test_the_look_carries_the_stored_spin_for_the_page_to_turn_the_body_by():
    look = bl.planet_look(_planet(rotation_period_hours=23.9345, spin_axis_x=0.0, spin_axis_y=3.0, spin_axis_z=4.0,
                                  axial_tilt_deg=23.4392))
    assert look["spinH"] == 23.93
    assert look["axis"] == [0.0, 0.6, 0.8], "a unit vector"
    assert look["tiltDeg"] == 23.44


def test_a_retrograde_body_keeps_its_axis_against_the_normal_and_its_tilt_past_ninety():
    look = bl.planet_look(_planet(rotation_period_hours=5832.0, spin_axis_x=0.0, spin_axis_y=0.0, spin_axis_z=-1.0,
                                  axial_tilt_deg=177.36))
    assert look["axis"] == [0.0, 0.0, -1.0] and look["tiltDeg"] > 90


def test_a_body_with_no_stored_spin_is_left_without_one_not_given_a_made_up_period():
    for row in (_planet(), _planet(rotation_period_hours=None), _planet(rotation_period_hours=0.0),
                _planet(rotation_period_hours=float("nan"), spin_axis_x=None, spin_axis_y=None, spin_axis_z=None)):
        look = bl.planet_look(row)
        assert not {"spinH", "axis", "tiltDeg"} & set(look), look


def test_a_star_look_carries_its_spin_too():
    star = {"id": 3, "star_type": "G2V", "temperature_k": 5772.0, "luminosity_w": 3.828e26, "radius_km": 695_700.0,
            "rotation_period_hours": 609.0, "spin_axis_x": 0.0, "spin_axis_y": 0.0, "spin_axis_z": 1.0,
            "axial_tilt_deg": 7.25}
    look = bl.star_look(star)
    assert look["spinH"] == 609.0 and look["axis"] == [0.0, 0.0, 1.0] and look["tiltDeg"] == 7.25
    assert "spinH" not in bl.star_look(dict(star, rotation_period_hours=None))


def test_a_body_whose_day_is_its_orbit_is_marked_locked_and_the_caption_says_so():
    """UX.244: `locked` is the generator's own lock (`regional_temperature.is_locked`); a free rotator is not marked."""
    locked = bl.planet_look(_planet(rotation_period_hours=24.0 * 365.25 * 0.1, period_years=0.1))
    free = bl.planet_look(_planet(rotation_period_hours=17.0, period_years=0.1))
    unknown = bl.planet_look(_planet(rotation_period_hours=17.0))
    assert locked["locked"] is True and "locked" not in free and "locked" not in unknown
    assert "Tidally locked" in bl.look_note(locked) and "its star" in bl.look_note(locked)
    assert "its planet" in bl.look_note(dict(locked, kind="moon"))
    assert "turning much faster" in bl.look_note(free)
    assert bl.look_note({"quality": "x"}) is None


# --- Land colour from the stored ground types (MAP.274) ---------------------------------------------------------------------

def _ground(**shares):
    return {f"ground_{name}_fraction": share for name, share in shares.items()}


def _channels(hexcolour):
    return tuple(int(hexcolour[i:i + 2], 16) for i in (1, 3, 5))


def test_a_world_with_stored_ground_takes_the_lab_colour_of_its_land_not_the_temperature_band():
    look = bl.planet_look(dict(_planet(), **_ground(ocean=0.71, basaltic=0.29)))
    assert look["quality"] == "composition" and look["landKnown"] == 1.0
    assert look["color"] != bl.land_color(288.0)


def test_sulfur_land_is_yellow_and_the_lighter_rocks_are_lighter():
    colours = {name: _channels(bl.planet_look(dict(_planet(), **_ground(**{name: 1.0})))["color"])
               for name in ("sulfur", "basaltic", "anorthositic", "carbon")}
    red, green, blue = colours["sulfur"]
    assert red > blue and blue < 0.8 * green  # the 450 nm band of sulfur is far below the other two
    assert sum(colours["anorthositic"]) > sum(colours["basaltic"]) > sum(colours["carbon"])


def test_ocean_and_ice_are_not_land_and_do_not_change_the_land_colour():
    dry = bl.planet_look(dict(_planet(), **_ground(felsic=0.4)))
    wet = bl.planet_look(dict(_planet(), **_ground(ocean=0.5, ice=0.1, felsic=0.4)))
    assert dry["color"] == wet["color"]


def test_land_with_no_measured_colour_is_left_out_and_too_little_known_falls_back_to_the_estimate():
    half = bl.planet_look(dict(_planet(), **_ground(basaltic=0.5, metal=0.5)))
    assert half["quality"] == "composition" and half["landKnown"] == 0.5
    assert half["color"] == bl.planet_look(dict(_planet(), **_ground(basaltic=1.0)))["color"]
    mostly_metal = bl.planet_look(dict(_planet(), **_ground(basaltic=0.2, metal=0.8)))
    assert mostly_metal["quality"] == "estimate" and mostly_metal["color"] == bl.land_color(288.0)


def test_a_world_with_no_stored_ground_or_no_land_keeps_the_estimate_and_a_giant_is_untouched():
    assert bl.planet_look(_planet())["quality"] == "estimate"
    assert bl.planet_look(dict(_planet(), **_ground(ocean=1.0)))["quality"] == "estimate"
    giant = bl.planet_look(dict(_planet(), body_type="g", **_ground(basaltic=1.0)))
    assert giant["quality"] == "estimate" and giant["color"] == bl.giant_color(288.0)


def test_the_note_says_where_the_land_colour_came_from_and_what_was_left_out():
    whole = bl.planet_look(dict(_planet(), **_ground(basaltic=1.0)))
    assert "stored ground types" in bl.look_note(whole)
    partial = bl.planet_look(dict(_planet(), **_ground(basaltic=0.6, metal=0.4)))
    assert "40 percent of the land has no measured colour" in bl.look_note(partial)


def test_the_colour_account_names_the_inputs_the_picture_used_and_calls_the_estimate_partial():
    """UX.158: the info panel line says which stored values gave the colour and that it is partial while the composition is not used."""
    earth = bl.look_note(bl.planet_look(_planet()))
    assert earth.startswith("Partial estimate: Ground colour from its mean surface temperature (288 K) alone; "
                            "sea colour from its ocean class (neutral); limb glow") or "no limb glow" in earth
    assert "sea colour from its ocean class (neutral)" in earth
    assert "Not used yet: surface composition, clouds and haze (not stored)" in earth
    bare = bl.colour_account(bl.planet_look(_planet(surface_temperature_k=None, ocean_class=None, atmosphere="None")))
    assert "Ground colour is a default: no temperature is stored" in bare
    assert "sea colour is a default: no ocean class is stored" in bare and "limb" not in bare
    giant = bl.colour_account(bl.planet_look(_planet(body_type="g", surface_temperature_k=165.0)))
    assert giant.startswith("Partial estimate: Cloud-top colour from its mean cloud-top temperature (165 K) alone")
    assert "chemistry bands and haze" in giant and "sea colour" not in giant
    rogue = bl.colour_account(bl.planet_look(_planet(), "rogue_planet"))
    assert "limb" not in rogue


def test_the_colour_account_of_a_star_is_full_with_a_temperature_and_partial_without():
    star = {"id": 3, "star_type": "G2V", "temperature_k": 5772.0, "luminosity_w": 3.828e26, "radius_km": 695_700.0}
    assert bl.colour_account(bl.star_look(star)) == "Full: colour of a model atmosphere at its effective temperature (5,772 K) and gravity"
    assert bl.colour_account(bl.star_look(dict(star, temperature_k=None))).startswith("Partial: colour from its spectral class alone")
    assert bl.colour_account({"body": "nebula"}) is None
    # MAP.274: a solid world coloured from its stored ground types no longer lacks the surface composition, but the cloud and haze are still not stored
    composed = bl.colour_account(dict(bl.planet_look(_planet()), quality="composition"))
    assert composed.startswith("Partial estimate: Ground colour from the stored ground types (lab reflectance)")
    assert "Not used yet: clouds and haze (not stored)" in composed and "surface composition" not in composed


def _air_look(pressure_kpa, gravity_g=1.0, **over):
    return bl.planet_look(_earth_air(atmospheric_pressure_pa=pressure_kpa * 1000.0, gravity_g=gravity_g, **over))


def test_the_glow_grows_with_the_column_of_the_air_from_none_to_thick():
    """MAP.277: the Earth's air glows clearly, Mars's (0.6 kPa at 0.38 g) hardly at all, Venus's and Titan's most; no air has no glow fields."""
    earth = _air_look(101.325)
    mars = _air_look(0.6, 0.38)
    venus = _air_look(9200.0, 0.9)
    titan = _air_look(146.7, 0.138)
    assert 0.5 < earth["glowStrength"] < bl.GLOW_STRENGTH_MAX
    assert bl.GLOW_STRENGTH_MIN <= mars["glowStrength"] < 0.05
    assert venus["glowStrength"] > earth["glowStrength"] and titan["glowStrength"] > earth["glowStrength"]
    assert venus["glowStrength"] <= bl.GLOW_STRENGTH_MAX
    airless = bl.planet_look(_planet(atmosphere="None", atmospheric_pressure_pa=0.0))
    assert not {"glow", "glowStrength", "glowWidth", "twilight"} & set(airless)


def test_an_air_too_thin_to_show_has_no_glow_at_all():
    look = _air_look(0.01, 0.38)
    assert look["atm"] is True and not {"glow", "glowStrength", "glowWidth", "twilight"} & set(look)


def test_the_twilight_grows_with_the_glow_and_a_row_with_no_star_keeps_it_and_its_strength():
    assert _air_look(0.6, 0.38)["twilight"] < _air_look(101.325)["twilight"] < _air_look(9200.0, 0.9)["twilight"] <= bl.TWILIGHT_MAX
    unknown = bl.planet_look(_earth_air(_star_temperature_k=None))
    assert "glow" not in unknown and unknown["twilight"] > 0.0 and unknown["glowStrength"] > 0.0      # the air still glows and wraps the light; no star colour to glow in


def test_the_glow_reaches_farther_past_the_limb_for_a_light_hot_air_on_a_small_body():
    earth = _air_look(101.325)["glowWidth"]
    titan = _air_look(146.7, 0.138, radius_km=2575.0, surface_temperature_k=94.0)["glowWidth"]
    assert bl.GLOW_WIDTH_RANGE[0] <= earth <= bl.GLOW_WIDTH_RANGE[1]
    assert titan > earth                                       # a scale height that is a larger part of the radius
    assert _air_look(101.325, surface_temperature_k=400.0)["glowWidth"] > earth      # hotter air is thicker
    # with no gases stored the strength comes from the column, and there is no width to work out
    assert "glowWidth" not in bl.planet_look(_planet(gravity_g=1.0, _star_temperature_k=5772.0))


def test_a_giants_glow_is_the_fixed_strong_one():
    giant = bl.planet_look(fold_air(_planet(body_type="g", planet_class="j", gravity_g=2.5, atmosphere="Hydrogen-Helium", atmospheric_pressure_pa=1e5,
                                            surface_temperature_k=130.0, _star_temperature_k=5772.0, p_h2_kpa=86.0)))
    assert giant["glowStrength"] == bl.GIANT_GLOW_STRENGTH


def test_a_rogue_planet_has_no_glow_or_twilight():
    look = bl.rogue_planet_look({"id": 1, "planet_type": "t", "surface_temperature_k": 300.0, "radius_km": 6000.0, "atmospheric_pressure_pa": 1e5,
                                 "gravity_g": 1.0, "atmosphere": "Nitrogen-Oxygen"})
    assert not {"glow", "glowStrength", "glowWidth", "twilight"} & set(look)


def test_a_solid_body_carries_its_stored_roundness_and_a_giant_does_not():
    assert bl.planet_look(_planet(roundness_ratio=0.31234567))["roundness"] == 0.3123
    assert "roundness" not in bl.planet_look(_planet())
    assert "roundness" not in bl.planet_look(_planet(body_type="g", roundness_ratio=None))
