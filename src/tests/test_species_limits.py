# tests/test_species_limits.py

"""
POP.12: the environment limits of a species and the lists of what could live where. The first half needs no
database (the human table with rows replaced, scored through the same six domains); the second runs the lists and
the API against a real throwaway database.
"""

import dataclasses

import pytest

from planetgen import tuning
from planetgen.cli import species_check
from planetgen.db import species_fit, store
from planetgen.physics import atmosphere, habitability, species_limits
from planetgen.physics.protection_limits import LIMIT_ROWS, limit
from tests.publicids import pid
from tests.test_admin_edits import _saved_system
from tests.test_facility_types import anonymous  # noqa: F401 (fixture)
from tests.test_population import galaxy as civilized_galaxy  # noqa: F401 (fixture)


def _world(**changes):
    """An Earth-like `habitability.World`, with `changes`."""
    world = habitability.World(
        pressure_kpa=101.3, temperature_c=15.0, gases={"O2": 0.21, "CO2": 0.0004}, gravity_ms2=9.8, flux_variation=0.0,
        relative_humidity=0.5, solvent="water", water_activity=1.0, chaotropicity_kj_kg=0.0, ph=8.0,
        inventories={"C": 1.0, "N": 1.0, "H": 1.0}, phosphate_umol_l=2.3, energy_flux_w_m2=80.0,
        surface_dose_sv_yr=0.002, habitat_dose_sv_yr=0.002, water_source="liquid", life_stage="biosphere",
        life_on_surface=True)
    return dataclasses.replace(world, **changes)


# --- The human table is unchanged --------------------------------------------------------------------------------

def test_the_human_bands_are_the_ones_the_module_always_had():
    human = habitability.HUMAN
    assert (human.pressure, human.temperature, human.dose) == (
        habitability.PRESSURE_BAND, habitability.TEMPERATURE_BAND, habitability.DOSE_BAND)
    assert human.o2_high == habitability.O2_HIGH_BAND
    assert (human.armstrong, human.o2_min, human.wet_bulb) == (
        habitability.ARMSTRONG_KPA, habitability.INSPIRED_O2_MIN_KPA, habitability.WET_BULB_CRIT_C)


@pytest.mark.parametrize("world", [
    _world(), _world(pressure_kpa=3.0, temperature_c=-60.0), _world(pressure_kpa=9000.0, temperature_c=80.0),
    _world(surface_dose_sv_yr=3.0), _world(gases={"O2": 0.4, "CO2": 0.02}),
])
def test_an_empty_override_scores_exactly_as_the_human_index_does(world):
    assert species_limits.assess(world, {})["domains"] == habitability.phi4(world)


def test_a_replaced_row_moves_only_the_domain_it_belongs_to():
    world = _world(temperature_c=48.0, relative_humidity=0.05)
    human = habitability.phi4(world)
    heat_adapted = habitability.phi4(world, habitability.bands({"comfort_high_c": 55.0, "heat_garment_c": 65.0}))
    assert heat_adapted["Temperature"][1] == "Blue" and human["Temperature"][1] != "Blue"
    for domain in ("Pressure", "Air", "Water", "Radiation"):
        assert heat_adapted[domain] == human[domain]


# --- Validation ------------------------------------------------------------------------------------------------

def test_validate_accepts_the_human_rows_and_numbers_come_back_as_floats():
    assert species_limits.validate({row.name: row.value for row in LIMIT_ROWS}) == {
        row.name: float(row.value) for row in LIMIT_ROWS}
    assert species_limits.validate(None) == {}


@pytest.mark.parametrize("overrides, text", [
    ({"nope_kpa": 1.0}, "not an environment limit"),
    ({"comfort_high_c": "hot"}, "must be a number"),
    ({"comfort_high_c": True}, "must be a number"),
    ({"pressure_no_gear_kpa": 0.0}, "above zero"),
    ({"pressure_no_gear_kpa": 400.0}, "out of order"),
    ({"comfort_low_c": 40.0}, "out of order"),
    ({"dose_green_sv_yr": 5.0}, "out of order"),
])
def test_validate_says_what_is_wrong(overrides, text):
    with pytest.raises(ValueError, match=text):
        species_limits.validate(overrides)


def test_a_temperature_limit_may_be_below_zero():
    assert species_limits.validate({"cold_outer_c": -200.0, "cold_suit_c": -150.0})["cold_outer_c"] == -200.0


# --- Generated from the homeworld -------------------------------------------------------------------------------

def test_an_earth_like_homeworld_gives_the_human_limits():
    assert species_limits.generated(_world()) == {}


HOMEWORLDS = [
    _world(pressure_kpa=1500.0, temperature_c=70.0, surface_dose_sv_yr=0.5, gases={"O2": 0.05}),
    _world(pressure_kpa=12.0, temperature_c=-40.0, gases={"O2": 0.6}),
    _world(pressure_kpa=30.0, temperature_c=-5.0, surface_dose_sv_yr=0.001, gases={"O2": 0.3}),
    _world(pressure_kpa=8000.0, temperature_c=110.0, surface_dose_sv_yr=40.0, gases={"O2": 0.02}),
    _world(pressure_kpa=101.3, temperature_c=25.0, surface_dose_sv_yr=7.0),
    _world(pressure_kpa=101.3, temperature_c=30.0, relative_humidity=0.95, gases={"O2": 0.21}),
    _world(pressure_kpa=60.0, temperature_c=20.0, flux_variation=0.8),
]


@pytest.mark.parametrize("world", HOMEWORLDS)
def test_a_species_passes_the_five_core_domains_of_its_own_homeworld(world):
    overrides = species_limits.generated(world)
    species_limits.validate(overrides)
    result = species_limits.assess(world, overrides)
    for domain in ("Pressure", "Temperature", "Air", "Radiation"):
        assert result["domains"][domain][1] != "Red", (domain, overrides)
    assert result["unaided"], result["failing"]


@pytest.mark.parametrize("world", HOMEWORLDS[:5])
def test_a_widened_side_keeps_the_width_of_its_green_and_yellow_bands(world):
    overrides = species_limits.generated(world)
    changed = [name for name in overrides if name.startswith("pressure_") and name.endswith("_kpa")]
    if changed:
        edges = habitability.bands(overrides).pressure
        human = habitability.HUMAN.pressure
        high = [edge / old for edge, old in zip(edges[4:], human[4:])]
        low = [edge / old for edge, old in zip(edges[:4], human[:4])]
        assert max(high) - min(high) < 1e-6 or max(low) - min(low) < 1e-6


def test_the_oxygen_limits_follow_the_air_of_the_homeworld():
    rich = _world(gases={"O2": 0.6}, pressure_kpa=101.3)
    overrides = species_limits.generated(rich)
    assert overrides["inspired_o2_chronic_kpa"] > limit("inspired_o2_chronic_kpa")
    assert overrides["inspired_o2_acute_kpa"] / overrides["inspired_o2_chronic_kpa"] == pytest.approx(
        limit("inspired_o2_acute_kpa") / limit("inspired_o2_chronic_kpa"))
    assert species_limits.assess(rich, {})["domains"]["Air"][1] != "Blue"
    assert species_limits.assess(rich, overrides)["domains"]["Air"][1] == "Blue"


# --- Could it live there? ---------------------------------------------------------------------------------------

def test_a_human_like_species_could_live_on_earth_unaided_and_not_on_a_hot_high_dose_world():
    assert species_limits.assess(_world(), {}, 5.0)["unaided"]
    harsh = _world(temperature_c=150.0, surface_dose_sv_yr=50.0)
    result = species_limits.assess(harsh, {}, 5.0)
    assert not result["unaided"] and {"Temperature", "Radiation"} <= set(result["failing"])


def test_with_equipment_needs_the_tech_level_for_the_kit():
    harsh = _world(temperature_c=150.0, surface_dose_sv_yr=50.0)
    kit_level = tuning.SPECIES_EQUIPMENT_MIN_TECH[4]
    assert species_limits.assess(harsh, {}, kit_level)["equipped"]
    assert not species_limits.assess(harsh, {}, kit_level - 0.1)["equipped"]
    assert not species_limits.assess(harsh, {}, None)["equipped"]


def test_no_kit_makes_a_world_without_drinkable_water_habitable():
    dry = _world(solvent=None, water_source=None)
    result = species_limits.assess(dry, {}, 7.0)
    assert result["domains"]["Water"][1] == "Red" and not result["equipped"] and not result["unaided"]


# --- The database half ------------------------------------------------------------------------------------------

def _earthlike(conn, planet_id, **changes):
    """Makes planet `planet_id` a scored, Earth-like rocky world with a biosphere, then applies `changes`."""
    values = {
        **{f"p_{gas}_kpa": 0.0 for gas in atmosphere.SPECIES}, "ocean_class": None,
        "body_type": "t", "atmospheric_pressure_pa": 101325.0, "surface_temperature_k": 288.15, "gravity_g": 1.0,
        "hydrosphere": "surface ocean", "ocean_ph": 8.0, "water_activity": 1.0, "phosphorus": "limited",
        "water_mass_fraction": 0.0002, "surface_dose_msv_yr": 2.0, "dose_ground_msv_yr": 2.0, "mass_kg": 5.97e24,
        "radius_km": 6371.0, "energy_flux_w_m2": 80.0, "flux_variation": 0.0, "life_chemical": "chlorophyll",
        "p_o2_kpa": 21.2, "p_n2_kpa": 79.0, "p_co2_kpa": 0.04, "p_h2o_kpa": 1.0, "hab_note": None, **changes,
    }
    conn.execute(f"UPDATE planets SET {', '.join(f'{name} = ?' for name in values)} WHERE id = ?",
                 (*values.values(), planet_id))
    conn.execute("DELETE FROM planet_evolutionary_paragraphs WHERE planet_id = ?", (planet_id,))
    conn.execute("INSERT INTO planet_evolutionary_paragraphs (planet_id, position, paragraph) VALUES (?, 0, ?)",
                 (planet_id, tuning.EVOLUTIONARY_TEXT["technological_civilization"]))


def _species(conn, planet_id, system_id, name, tech):
    conn.execute(
        "INSERT INTO species (name, homeworld_planet_id, star_system_id, life_chemical, life_stage, build, climate, "
        "size, civilization_age_years, era, spacefaring, tech_level) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (name, planet_id, system_id, "chlorophyll", "technological_civilization", "medium", "temperate", "medium",
         60000.0, "established", 1, tech))
    return conn.execute("SELECT id FROM species WHERE name = ?", (name,)).fetchone()["id"]


@pytest.fixture
def galaxy(mysql_config):
    """One saved system with three rocky planets: an Earth-like one that is a species' homeworld, a hot one that is
    another species' homeworld, and a harsh one. Returns the ids."""
    _sector_id, system_id = _saved_system(mysql_config, want=lambda s: sum(p.body_type != "a" for p in s.planets) >= 3)
    conn = store.get_connection(mysql_config)
    try:
        planets = [r["id"] for r in conn.execute(
            "SELECT id FROM planets WHERE star_system_id = ? AND body_type <> 'a' ORDER BY id", (system_id,)).fetchall()][:3]
        earth, hot, harsh = planets
        conn.execute("DELETE FROM moons WHERE star_system_id = ?", (system_id,))
        conn.execute(f"DELETE FROM planets WHERE star_system_id = ? AND id NOT IN ({', '.join('?' for _ in planets)})",
                     (system_id, *planets))
        _earthlike(conn, earth)
        _earthlike(conn, hot, surface_temperature_k=423.15, atmospheric_pressure_pa=3_000_000.0,
                   surface_dose_msv_yr=40_000.0, p_o2_kpa=150.0, p_n2_kpa=2850.0, p_h2o_kpa=20.0)
        _earthlike(conn, harsh, surface_temperature_k=900.0, atmospheric_pressure_pa=101325.0,
                   surface_dose_msv_yr=5_000_000.0)
        humans = _species(conn, earth, system_id, "Terrans", 5.0)
        hot_ones = _species(conn, hot, system_id, "Scalders", 5.0)
        conn.commit()
    finally:
        conn.close()
    return {"system_id": system_id, "earth": earth, "hot": hot, "harsh": harsh, "humans": humans, "hot_ones": hot_ones}


def test_each_species_gets_limits_that_fit_its_own_homeworld(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        limits = species_fit.limits_for_species(conn)
        assert limits[galaxy["humans"]] == {} and limits[galaxy["hot_ones"]]
        for species_id, planet_id in ((galaxy["humans"], galaxy["earth"]), (galaxy["hot_ones"], galaxy["hot"])):
            assert planet_id in [b["id"] for b in species_fit.bodies_for_species(
                conn, species_id, system_id=galaxy["system_id"])["unaided"] if b["kind"] == "planet"]
    finally:
        conn.close()


def test_a_body_lists_the_species_that_could_live_there(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        earth = species_fit.species_for_body(conn, "planet", galaxy["earth"])
        names = [s["name"] for s in earth["unaided"]]
        assert earth["scored"] and "Terrans" in names
        hot = species_fit.species_for_body(conn, "planet", galaxy["hot"])
        assert [s["name"] for s in hot["unaided"]] == ["Scalders"]
        assert "Terrans" in [s["name"] for s in hot["equipped"]]  # tech 5 makes full life support
        assert all(s["unaided"] is False and s["equipped"] for s in hot["equipped"])
        harsh = species_fit.species_for_body(conn, "planet", galaxy["harsh"])
        assert harsh["unaided"] == [] and len(harsh["equipped"]) == 2
        assert all("temperature" in s["failing"] for s in harsh["equipped"])
        assert {"temperature", "radiation"} <= set(next(s for s in harsh["equipped"] if s["name"] == "Terrans")["failing"])
        assert species_fit.species_for_body(conn, "planet", 10 ** 12) is None
    finally:
        conn.close()


def test_a_species_without_a_tech_level_has_no_equipped_list(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE species SET tech_level = NULL")
        conn.commit()
        assert species_fit.species_for_body(conn, "planet", galaxy["harsh"])["equipped"] == []
    finally:
        conn.close()


def test_a_gas_giant_has_no_species(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE planets SET body_type = 'g' WHERE id = ?", (galaxy["harsh"],))
        conn.commit()
        found = species_fit.species_for_body(conn, "planet", galaxy["harsh"])
        assert found["scored"] is False and found["unaided"] == [] and found["equipped"] == []
    finally:
        conn.close()


def test_the_other_way_round_lists_the_worlds_of_a_system_or_a_sector(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        sector_id = conn.execute("SELECT sector_id FROM star_systems WHERE id = ?",
                                 (galaxy["system_id"],)).fetchone()["sector_id"]
        for scope in ({"system_id": galaxy["system_id"]}, {"sector_id": sector_id}):
            found = species_fit.bodies_for_species(conn, galaxy["humans"], **scope)
            assert found["scored"] >= 3 and not found["truncated"]
            assert [b["id"] for b in found["unaided"] if b["kind"] == "planet"] == [galaxy["earth"]]
            assert galaxy["hot"] in [b["id"] for b in found["equipped"]]
        assert species_fit.bodies_for_species(conn, 10 ** 12, system_id=galaxy["system_id"]) is None
    finally:
        conn.close()


def test_the_limits_table_marks_the_rows_that_differ_from_the_human_ones():
    rows = species_fit.limits_table({"comfort_high_c": 55.0})
    assert [row["name"] for row in rows] == [row.name for row in LIMIT_ROWS]
    changed = [row for row in rows if row["changed"]]
    assert [(row["name"], row["human"], row["value"]) for row in changed] == [("comfort_high_c", 31.0, 55.0)]


# --- The API ----------------------------------------------------------------------------------------------------

def test_the_api_lists_species_for_a_planet_and_worlds_for_a_species(mysql_config, galaxy, anonymous):
    planet = pid("planet", galaxy["hot"])
    body = anonymous.get(f"/api/planets/{planet}/habitable-species").get_json()
    assert body["id"] == planet and body["kind"] == "planet"
    assert [s["name"] for s in body["unaided"]] == ["Scalders"]
    assert body["unaided"][0]["homeworld_planet_id"] == planet
    assert set(body["unaided"][0]["domains"]) == {"pressure", "temperature", "air", "water", "radiation", "suitability"}
    system = pid("system", galaxy["system_id"])
    worlds = anonymous.get(f"/api/species/{galaxy['humans']}/habitable-worlds/systems/{system}").get_json()
    assert worlds["unaided"][0]["id"] == pid("planet", galaxy["earth"]) and worlds["unaided"][0]["kind"] == "planet"
    assert worlds["unaided"][0]["star_system_id"] == system
    assert anonymous.get(f"/api/species/{galaxy['humans']}/habitable-worlds/systems/{system}").status_code == 200


def test_the_api_404s_for_unknown_things(mysql_config, galaxy, anonymous):
    assert anonymous.get(f"/api/species/999999/habitable-worlds/systems/{pid('system', galaxy['system_id'])}").status_code == 404
    assert anonymous.get("/api/moons/FE81000A2B-0000005-000/habitable-species").status_code == 404


def test_species_detail_carries_its_limits(mysql_config, galaxy, anonymous):
    body = anonymous.get(f"/api/species/{galaxy['hot_ones']}").get_json()
    rows = {row["name"]: row for row in body["limits"]}
    assert set(rows) == {row.name for row in LIMIT_ROWS}
    assert rows["comfort_high_c"]["changed"] and rows["comfort_high_c"]["value"] > rows["comfort_high_c"]["human"]
    human = anonymous.get(f"/api/species/{galaxy['humans']}").get_json()
    assert not any(row["changed"] for row in human["limits"])


# --- POP.14: a species is only placed where it can live ------------------------------------------------------------

def test_a_species_on_its_own_world_is_not_judged_on_suitability():
    abiotic = _world(life_stage="none")
    assert species_limits.assess(abiotic, {})["failing"] == ["Suitability"]
    assert species_limits.assess(abiotic, {}, native=True)["failing"] == []
    assert species_limits.assess(abiotic, {}, native=True)["unaided"]


def test_the_reasons_are_plain_words_with_no_numbers():
    cases = [
        (_world(pressure_kpa=1.0), "the air is too thin"), (_world(pressure_kpa=90000.0), "crushingly dense"),
        (_world(temperature_c=-150.0), "too cold"), (_world(temperature_c=180.0), "too hot"),
        (_world(surface_dose_sv_yr=5000.0), "radiation"), (_world(solvent=None, water_source=None), "no water"),
    ]
    for world, text in cases:
        failing = species_limits.assess(world, {})["failing"]
        reasons = " ".join(species_limits.explain(world, {}, failing))
        assert text in reasons and not any(ch.isdigit() for ch in reasons), (text, reasons)


def test_a_placement_is_refused_with_the_failing_conditions(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        assert species_fit.check_placement(conn, galaxy["humans"], "planet", galaxy["earth"])["unaided"]
        with pytest.raises(species_fit.PlacementError) as refused:
            species_fit.check_placement(conn, galaxy["humans"], "planet", galaxy["harsh"])
        assert str(refused.value).startswith("The Terrans could not live on ")
        assert "too hot" in str(refused.value) and set(refused.value.failing) >= {"temperature", "radiation"}
    finally:
        conn.close()


def test_a_habitat_may_be_placed_where_only_equipment_lets_the_species_live(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        with pytest.raises(species_fit.PlacementError):
            species_fit.check_placement(conn, galaxy["humans"], "planet", galaxy["hot"])
        assert species_fit.check_placement(conn, galaxy["humans"], "planet", galaxy["hot"], habitat=True)["equipped"]
        conn.execute("UPDATE species SET tech_level = 1.0")  # no kit for a hot world at tech 1
        conn.commit()
        with pytest.raises(species_fit.PlacementError):
            species_fit.check_placement(conn, galaxy["humans"], "planet", galaxy["hot"], habitat=True)
    finally:
        conn.close()


def test_a_gas_giant_a_missing_species_and_a_missing_body_are_refused_an_unscored_body_passes(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE planets SET body_type = 'g' WHERE id = ?", (galaxy["harsh"],))
        conn.execute("UPDATE planets SET energy_flux_w_m2 = NULL WHERE id = ?", (galaxy["hot"],))
        conn.commit()
        with pytest.raises(species_fit.PlacementError, match="gas giant"):
            species_fit.check_placement(conn, galaxy["humans"], "planet", galaxy["harsh"])
        with pytest.raises(species_fit.PlacementError, match="no species"):
            species_fit.check_placement(conn, 999999, "planet", galaxy["earth"])
        with pytest.raises(species_fit.PlacementError, match="no planet"):
            species_fit.check_placement(conn, galaxy["humans"], "planet", 999999)
        assert species_fit.check_placement(conn, galaxy["humans"], "planet", galaxy["hot"]) == {}
    finally:
        conn.close()


def test_the_report_lists_a_species_whose_homeworld_no_longer_fits(mysql_config, galaxy, capsys, monkeypatch):
    # Limits made from the homeworld follow it; the report matters once they are stored (DB.30, POP.13).
    conn = store.get_connection(mysql_config)
    try:
        stored = species_fit.limits_for_species(conn)
        monkeypatch.setattr(species_fit, "stored_overrides", lambda conn, ids: {i: stored[i] for i in ids})
        assert species_fit.inconsistencies(conn) == []
        # Its homeworld's radiation is refreshed far past anything its limits were made for.
        conn.execute("UPDATE planets SET surface_dose_msv_yr = 9e9 WHERE id = ?", (galaxy["earth"],))
        conn.commit()
        found = species_fit.inconsistencies(conn)
        assert [(f["name"], f["failing"]) for f in found] == [("Terrans", ["radiation"])]
        assert found[0]["id"] == galaxy["earth"] and "radiation" in found[0]["reasons"]
    finally:
        conn.close()
    argv = ["--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port), "--mysql-user",
            mysql_config.user, "--mysql-password", mysql_config.password, "--mysql-database", mysql_config.database]
    assert species_check.main(argv) == 1
    out = capsys.readouterr().out
    assert "Terrans" in out and "1 of 2 species is placed outside its limits." in out


def test_the_report_is_clean_for_a_consistent_galaxy(mysql_config, galaxy, capsys):
    argv = ["--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port), "--mysql-user",
            mysql_config.user, "--mysql-password", mysql_config.password, "--mysql-database", mysql_config.database]
    assert species_check.main(argv) == 0
    assert "All 2 species are placed inside their limits." in capsys.readouterr().out


def test_a_species_whose_own_homeworld_is_an_abiotic_surface_is_still_consistent(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM planet_evolutionary_paragraphs WHERE planet_id = ?", (galaxy["earth"],))
        conn.execute("UPDATE planets SET life_chemical = NULL WHERE id = ?", (galaxy["earth"],))
        conn.commit()
        assert species_fit.inconsistencies(conn) == []
        # ...while another species placed there is judged on Suitability as well.
        with pytest.raises(species_fit.PlacementError, match="food"):
            species_fit.check_placement(conn, galaxy["hot_ones"], "planet", galaxy["earth"])
    finally:
        conn.close()


def test_creating_a_species_checks_its_homeworld(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        species_fit.check_homeworld(conn, galaxy["earth"])
        conn.execute("UPDATE planets SET hydrosphere = NULL, water_mass_fraction = NULL, p_h2o_kpa = 0 WHERE id = ?",
                     (galaxy["harsh"],))
        conn.commit()
        with pytest.raises(species_fit.PlacementError, match="no water"):
            species_fit.check_homeworld(conn, galaxy["harsh"])
        assert species_fit.check_homeworld(conn, 999999) is None
    finally:
        conn.close()


def test_a_generated_galaxys_species_are_all_inside_their_limits(mysql_config, civilized_galaxy):
    """Boss's POP.14: the consistency check over a galaxy the population pass has just made."""
    from planetgen.population import model as population
    conn = store.get_connection(mysql_config)
    try:
        assert population.run_pass(conn)["new_species"] >= 2
        assert species_fit.inconsistencies(conn) == []
        for row in conn.execute("SELECT id, homeworld_planet_id FROM species").fetchall():
            species_fit.check_placement(conn, row["id"], "planet", row["homeworld_planet_id"], native=True)
    finally:
        conn.close()
