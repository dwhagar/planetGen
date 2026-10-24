# tests/test_parent_link.py

"""
DB.73: every kind of object says what it orbits in `parent_type` and `parent_id`; a moon can orbit a free (rogue) planet
and is stored, numbered, loaded, fingerprinted and deleted with it; a planet that leaves its star keeps its row and ID.
"""

import pytest

from planetgen.db import fingerprint, store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.system import StarSystem
from planetgen.util import draw

pytestmark = pytest.mark.fresh_database


def _system_with_everything():
    for seed in range(200):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if system.comets and any(getattr(p, "moons", None) for p in system.planets) and any(p.body_type == "a" for p in system.planets):
            return system, cfg
    pytest.fail("no system with a belt, a comet and a moon in 200 seeds")


def _moon():
    """The moons of some generated planet, at least two of them."""
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


def _rogue_planet(conn, moons=0):
    """A saved rogue planet with its own ID and `moons` moons of a generated system under it."""
    rogue_id = store.insert_rogue_planet(conn, RoguePlanet(SystemConfig()))
    store.assign_uids(conn, phenomenon=("rogue_planets", rogue_id))
    ids = [store.insert_moon_of_rogue_planet(conn, moon, rogue_id, position)
           for position, moon in enumerate(_moon()[:moons])]
    return rogue_id, ids


def test_every_kind_a_system_holds_says_what_it_orbits(mysql_config):
    system, cfg = _system_with_everything()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        star = conn.execute("SELECT id, parent_type, parent_id FROM stars WHERE star_system_id = ?",
                            (system_id,)).fetchone()
        assert (star["parent_type"], star["parent_id"]) == ("system", system_id)
        for table in ("planets", "asteroid_belts", "comets"):
            rows = conn.execute(f"SELECT star_id, parent_type, parent_id FROM {table} WHERE star_system_id = ?",
                                (system_id,)).fetchall()
            assert rows, table
            for row in rows:
                assert (row["parent_type"], row["parent_id"]) == ("star", row["star_id"]), table
        moons = conn.execute("SELECT planet_id, parent_type, parent_id, rogue_planet_id FROM moons"
                             " WHERE star_system_id = ?", (system_id,)).fetchall()
        assert moons and all((m["parent_type"], m["parent_id"], m["rogue_planet_id"]) == ("planet", m["planet_id"], None)
                             for m in moons)
    finally:
        conn.rollback()
        conn.close()


def test_a_moon_of_a_rogue_planet_stores_loads_and_is_numbered_under_it(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            rogue_id, moon_ids = _rogue_planet(conn, moons=2)
            store.assign_uids(conn, phenomenon=("rogue_planets", rogue_id))
        assert len(moon_ids) == 2
        rows = conn.execute("SELECT * FROM moons WHERE rogue_planet_id = ? ORDER BY orbital_index",
                            (rogue_id,)).fetchall()
        assert [(r["parent_type"], r["parent_id"], r["planet_id"], r["star_system_id"], r["star_id"]) for r in rows] == [
            ("rogue_planet", rogue_id, None, None, None)] * 2
        rogue_uid = conn.execute("SELECT uid FROM rogue_planets WHERE id = ?", (rogue_id,)).fetchone()["uid"]
        from planetgen.galaxy import object_uid
        parent = object_uid.from_bytes(rogue_uid)
        bodies = sorted(object_uid.body_of(object_uid.from_bytes(r["uid"])) for r in rows)
        assert bodies == [1, 2]
        assert all(object_uid.with_body(parent, object_uid.body_of(object_uid.from_bytes(r["uid"])))
                   == object_uid.from_bytes(r["uid"]) for r in rows)
        loaded = store.load_moons_of_rogue_planet(conn, rogue_id)
        assert [m.db_id for m in loaded] == moon_ids and all(m.is_moon for m in loaded)
        assert [m.name for m in loaded] == [r["name"] for r in rows]
    finally:
        conn.rollback()
        conn.close()


def test_a_rogue_planets_moons_go_with_it_and_with_its_sector(mysql_config):
    sector = SpaceSector("Wanderer Sector", edge_ly=10.0)
    sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            rogue_id = conn.execute("SELECT id FROM rogue_planets WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
            for position, moon in enumerate(_moon()[:2]):
                store.insert_moon_of_rogue_planet(conn, moon, rogue_id, position)
            store.assign_uids(conn, sector_id=sector_id)
        assert conn.execute("SELECT COUNT(*) AS n FROM moons WHERE uid IS NOT NULL AND rogue_planet_id = ?",
                            (rogue_id,)).fetchone()["n"] == 2
        with conn:
            conn.execute("DELETE FROM sectors WHERE id = ?", (sector_id,))
        assert conn.execute("SELECT COUNT(*) AS n FROM moons").fetchone()["n"] == 0
        assert conn.execute("SELECT COUNT(*) AS n FROM moon_reflection_spectrum").fetchone()["n"] == 0
    finally:
        conn.rollback()
        conn.close()


def test_a_sector_digest_reads_the_moons_of_its_rogue_planets(mysql_config):
    sector = SpaceSector("Digest Sector", edge_ly=10.0)
    sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        before = fingerprint.sector_digest(conn, sector_id=sector_id)
        with conn:
            rogue_id = conn.execute("SELECT id FROM rogue_planets WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
            store.insert_moon_of_rogue_planet(conn, _moon()[0], rogue_id, 0)
            store.assign_uids(conn, sector_id=sector_id)
        after = fingerprint.sector_digest(conn, sector_id=sector_id)
        assert after != before
        assert fingerprint.sector_digest(conn, sector_id=sector_id) == after
    finally:
        conn.close()


def test_a_planet_that_leaves_its_star_keeps_its_row_and_id(mysql_config):
    system, cfg = _system_with_everything()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
            store.assign_uids(conn, system_ids=[system_id])
        before = conn.execute("SELECT * FROM planets WHERE star_system_id = ? AND body_type <> 'a' ORDER BY id LIMIT 1",
                              (system_id,)).fetchone()
        with conn:
            store.set_parent(conn, "planets", before["id"], "free")
        after = conn.execute("SELECT * FROM planets WHERE id = ?", (before["id"],)).fetchone()
        assert (after["parent_type"], after["parent_id"]) == ("free", None)
        assert {k: v for k, v in after.items() if k not in ("parent_type", "parent_id")} == \
            {k: v for k, v in before.items() if k not in ("parent_type", "parent_id")}
        assert after["uid"] == before["uid"] and after["uid"] is not None
        with conn:
            store.set_parent(conn, "planets", before["id"], "star", before["star_id"])
        again = conn.execute("SELECT parent_type, parent_id FROM planets WHERE id = ?", (before["id"],)).fetchone()
        assert (again["parent_type"], again["parent_id"]) == ("star", before["star_id"])
    finally:
        conn.rollback()
        conn.close()


def test_a_parent_must_be_there_and_say_what_it_is(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            rogue_id, _moons = _rogue_planet(conn)
        for arguments, message in (
                (("planets", 1, "galaxy"), "no parent type"),
                (("planets", 1, "free", 5), "no parent id"),
                (("planets", 1, "star"), "every other parent has one"),
                (("belts", 1, "free"), "no parent link"),
                (("planets", 999999, "free"), "no planets row"),
                (("rogue_planets", rogue_id, "rogue_planet", rogue_id), "cannot orbit itself")):
            with pytest.raises(ValueError, match=message):
                store.set_parent(conn, *arguments)
    finally:
        conn.rollback()
        conn.close()
