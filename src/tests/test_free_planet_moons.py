# tests/test_free_planet_moons.py

"""
DB.94: a moon of a free planet (DB.73, `moons.rogue_planet_id`) is read where a moon of a system is: the object index
files it under its rogue planet and the planet's sector, `resolve_object` names the rogue planet as its parent, and the
Search page's moons panel lists it with the rogue planet where the planet's name goes.
"""

import pytest

from planetgen.db import objectindex, query, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.system import StarSystem
from planetgen.galaxy.sector import SpaceSector
from planetgen.util import draw

pytestmark = pytest.mark.fresh_database

PLACEMENT = {"center_x_pc": 5.0, "center_y_pc": 1.0, "center_z_pc": 2.0, "galactic_radius_pc": 5.5}


def _moons_of_some_planet():
    for seed in range(200):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        for planet in system.planets:
            if len(getattr(planet, "moons", ())) >= 2:
                return planet.moons
    pytest.fail("no planet with two moons in 200 seeds")


@pytest.fixture
def world(mysql_config):
    """A saved sector holding a rogue planet with two moons, and the sector's id."""
    sector_id = store.save_sector(SpaceSector("Free Sector", edge_ly=10.0), config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            rogue_id = store.insert_rogue_planet(conn, RoguePlanet(SystemConfig()), sector_id=sector_id, placement=PLACEMENT)
            store.assign_uids(conn, phenomenon=("rogue_planets", rogue_id))
            moon_ids = [store.insert_moon_of_rogue_planet(conn, moon, rogue_id, position)
                        for position, moon in enumerate(_moons_of_some_planet()[:2])]
        return mysql_config, sector_id, rogue_id, moon_ids
    finally:
        conn.close()


def test_the_object_index_files_a_free_moon_under_its_rogue_planet(world):
    config, sector_id, rogue_id, moon_ids = world
    conn = store.get_connection(config)
    try:
        objectindex.build(conn, ["moon", "rogue_planet"])
        rogue_uid = conn.execute("SELECT uid FROM rogue_planets WHERE id = %s", (rogue_id,)).fetchone()["uid"]
        rows = conn.execute("SELECT * FROM object_index WHERE kind = 'moon' ORDER BY source_id").fetchall()
        assert [row["source_id"] for row in rows] == moon_ids
        assert {row["parent_uid"] for row in rows} == {rogue_uid}
        assert {row["sector_id"] for row in rows} == {sector_id} and {row["system_id"] for row in rows} == {None}
        assert {row["placed"] for row in rows} == {1}
        report = objectindex.check(conn)
        assert (report["moon"]["missing"], report["moon"]["extra"], report["moon"]["different"]) == (0, 0, 0)
    finally:
        conn.close()


def test_a_refresh_sees_an_edit_to_a_free_moon_and_to_its_rogue_planet(world):
    config, _sector_id, rogue_id, moon_ids = world
    conn = store.get_connection(config)
    try:
        objectindex.build(conn, ["moon"], config.database)
        conn.execute("UPDATE moons SET name = 'Renamed Moon' WHERE id = %s", (moon_ids[0],))
        conn.commit()
        assert objectindex.refresh(conn, ["moon"], config.database)["moon"]["changed"] >= 1
        assert conn.execute("SELECT name FROM object_index WHERE kind = 'moon' AND source_id = %s",
                            (moon_ids[0],)).fetchone()["name"] == "Renamed Moon"
        conn.execute("UPDATE rogue_planets SET center_x_pc = NULL, center_y_pc = NULL, center_z_pc = NULL,"
                     " galactic_radius_pc = NULL WHERE id = %s", (rogue_id,))
        conn.commit()
        objectindex.refresh(conn, ["moon"], config.database)
        assert {row["placed"] for row in conn.execute("SELECT placed FROM object_index WHERE kind = 'moon'").fetchall()} == {0}
    finally:
        conn.close()


def test_resolve_object_names_the_rogue_planet_as_the_parent(world):
    config, _sector_id, rogue_id, moon_ids = world
    conn = store.get_connection(config)
    try:
        found = query.resolve_object(conn, "moon", moon_ids[0])
        assert [parent["kind"] for parent in found["parents"]] == ["galaxy", "rogue_planet"]
        assert found["parents"][1]["ref"] == query.object_ref.format("rogue_planet", rogue_id)
        assert found["positions"]["galaxy_pc"] == [5.0, 1.0, 2.0] and found["positions"]["system_km"] is None
        assert found["siblings"] == [query.object_ref.format("moon", moon_ids[1])]
    finally:
        conn.close()


def test_the_search_page_lists_a_free_moon_with_its_rogue_planet(world):
    config, sector_id, rogue_id, moon_ids = world
    conn = store.get_connection(config)
    try:
        tags = {facet: set() for facet in query.SEARCH_TAG_FACETS}
        panel = query.search(conn, {"moon_q": ""}, tags, only=["moons"])["results"]["moons"]
        rows = {row["id"]: row for row in panel["rows"]}
        assert set(rows) == set(moon_ids)
        row = rows[moon_ids[0]]
        assert row["system_name"] is None and row["star_system_id"] is None and row["sector_id"] == sector_id
        assert row["planet_name"] == conn.execute("SELECT name FROM rogue_planets WHERE id = %s",
                                                  (rogue_id,)).fetchone()["name"]
        assert [parent["kind"] for parent in row["parents"]] == ["sector", "rogue_planet"]
    finally:
        conn.close()
