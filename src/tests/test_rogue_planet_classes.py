"""
Rogue planet classes (GEN.8): every `PLANET_CLASSES` class carries an
`"r"` zone flag, a rogue planet draws its class from the flagged classes
that fit its type, radius and mass, and the class is stored (schema v47)
and shown on the phenomenon page, the sector's Contents and the Sector
Map.
"""

import pytest

from planetgen.web.app import create_app
from planetgen.api.config import Config

from planetgen.db import store
from planetgen.physics import constants
from planetgen import tuning
from planetgen.generation import phenomena_plausibility as pp
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import (RoguePlanet, default_rogue_planet_class, rogue_planet_class_candidates,
                                            rogue_planet_classes)

from planetgen.web.maps import starmap  # noqa: E402
from planetgen.web import sector_page, system_pages  # noqa: E402


class _Config(Config):
    TESTING = True
    SECRET_KEY = "test"


def test_every_class_has_a_rogue_flag_and_no_rogue_class_has_life():
    for code, data in tuning.PLANET_CLASSES.items():
        assert isinstance(data.get("r"), bool), code
        if data["r"]:
            assert not data.get("life_chemical"), code
            assert code not in tuning.HABITABLE_PLANET_CLASSES
    assert rogue_planet_classes("t") == ["C", "D", "S"]
    assert rogue_planet_classes("g") == ["I", "J", "T"]


def test_generated_rogues_get_a_fitting_rogue_class():
    for _ in range(300):
        planet = RoguePlanet(SystemConfig())
        assert planet.planet_class in rogue_planet_classes(planet.planet_type)
        assert planet.planet_class in rogue_planet_class_candidates(planet.planet_type, planet.radius_km,
                                                                    planet.mass_kg)
        assert f"Class {planet.planet_class}" in planet.to_paragraph_list()[0]


def test_brown_dwarfs_have_no_planet_class():
    assert RoguePlanet(SystemConfig(), mass_bin="brown-dwarf").planet_class is None


def test_a_rogue_too_big_for_every_class_gets_the_nearest():
    """A rocky rogue past every rocky class's ceiling gets the nearest:
    class S, the rocky super-Earth (GEN.38), not class C."""
    mass = 15 * constants.EARTH_MASS_TO_KG
    assert rogue_planet_class_candidates("t", 17000.0, mass) == ["S"]
    assert default_rogue_planet_class("t", 17000.0, mass) == "S"
    assert default_rogue_planet_class("g", 70000.0, constants.JUPITER_MASS_TO_KG) == "J"


def test_rocky_rogues_over_10000_km_are_class_s():
    """GEN.38: a rocky rogue bigger than Class C's 10,000 km ceiling is a
    rocky super-Earth (S), and one that fits S's radius and mass ranges
    may only be S or nothing smaller."""
    mass = 5 * constants.EARTH_MASS_TO_KG
    assert rogue_planet_class_candidates("t", 10500.0, mass) == ["S"]
    big = []
    for _ in range(600):
        planet = RoguePlanet(SystemConfig(), mass_bin="sub-neptune")
        if planet.planet_type == "t" and planet.radius_km > tuning.PLANET_CLASSES["C"]["radius_range"][1]:
            big.append(planet)
            assert planet.planet_class == "S"
    assert big


def test_from_dict_without_a_class_gives_the_default():
    planet = RoguePlanet(SystemConfig(), mass_bin="terrestrial")
    data = planet.to_dict()
    assert data["planet_class"] == planet.planet_class
    data.pop("planet_class")
    assert RoguePlanet.from_dict(data, SystemConfig()).planet_class == "C"


def test_plausibility_flags_a_non_rogue_class():
    record = {"phenomenon_type": "rogue-planet", "mass_kg": 1e27, "radius_km": 70000.0, "planet_type": "g",
              "planet_class": "M", "galactic_orbital_speed_kms": 205.7, "galactic_orbital_period_gy": 0.236,
              "galactic_orbital_phase_deg": 180.0, "galactic_min_update_interval_years": 1e-9}
    assert any("planet_class" in issue for issue in pp.check_hard_invariants(record))
    record["planet_class"] = "J"
    assert not any("planet_class" in issue for issue in pp.check_hard_invariants(record))


def test_phenomenon_page_contents_and_sector_map_show_the_class():
    with create_app(_Config).test_request_context("/"):
        fields = {label: (text, url) for label, text, url in
                  system_pages.phenomenon_fields("rogue_planet", {"planet_class": "C", "planet_type": "t"})}
        assert fields["Planet Class"][0] == "C"
        assert fields["Planet Class"][1]

        rogue = {"type": "rogue_planet", "id": 1, "name": "Wanderer", "descriptor": "terrestrial", "class": "C",
                 "radius_ly": 0, "distance_ly": 1.0}
        group = sector_page._rogue_group_row([rogue, {**rogue, "id": 2, "class": None, "descriptor": "brown dwarf"}], lambda system_id: f"/s/{system_id}")
        assert [m["details"] for m in group["members"]] == ["Class C, Terrestrial", "Brown dwarf"]

    data = starmap._cloud_data(lambda *a, **k: "/x", rogue, 0, 0, 0, 1)
    assert data["typeLabel"] == "Rogue Planet (Class C, terrestrial)"
    data = starmap._cloud_data(lambda *a, **k: "/x", {**rogue, "class": None}, 0, 0, 0, 1)
    assert data["typeLabel"] == "Rogue Planet (Terrestrial)"


def test_class_is_stored_and_backfilled_by_the_v47_migration(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            planet_id = store.insert_rogue_planet(conn, RoguePlanet(SystemConfig(), mass_bin="jupiter"))
            dwarf_id = store.insert_rogue_planet(conn, RoguePlanet(SystemConfig(), mass_bin="brown-dwarf"))
            rocky_id = store.insert_rogue_planet(conn, RoguePlanet(SystemConfig(), mass_bin="terrestrial"))
        row = conn.execute("SELECT planet_class FROM rogue_planets WHERE id = ?", (planet_id,)).fetchone()
        assert row["planet_class"] == "J"
        conn.execute("ALTER TABLE rogue_planets DROP COLUMN planet_class")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (46)")
        conn.commit()
    finally:
        conn.close()

    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION

    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        classes = {row["id"]: row["planet_class"]
                   for row in conn.execute("SELECT id, planet_class FROM rogue_planets").fetchall()}
    finally:
        conn.close()
    assert classes == {planet_id: "J", dwarf_id: None, rocky_id: "C"}
