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


def test_star_colour_is_the_blackbody_colour_of_its_temperature():
    from planetgen.physics import photometry
    star = {"id": 3, "star_type": "M3V", "yerkes_class": "V", "temperature_k": 3400.0, "luminosity_w": 1e24}
    look = bl.star_look(star)
    assert look["color"] == photometry.star_colour(3400.0).hex
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
