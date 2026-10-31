# tests/test_asteroid_type.py

"""
GEN.1160: a body under a roundness ratio of 1 is an Asteroid, with no class letter. It is stored with a NULL class, its system, object and map pages draw, and it validates.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation import validation
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import planet_classes
from tests.publicids import pid
from tests.test_admin_edits import client  # noqa: F401


def _has_asteroid_moon(system):
    return any(m.planet_class is None for p in system.planets if getattr(p, "moons", None) for m in p.moons)


def _saved_with_asteroid(mysql_config):
    for _ in range(300):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MOONS = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(cfg)
        if _has_asteroid_moon(system):
            assert not validation.check_star_system(system)
            sector = SpaceSector("Asteroid Sector", edge_ly=10.0)
            sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
            sector_id = store.save_sector(sector, config=mysql_config)
            conn = store.get_connection(mysql_config)
            try:
                system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
            finally:
                conn.close()
            return system_id
    pytest.skip("no system with an asteroid moon generated")


def test_an_asteroid_moon_is_stored_without_a_class_and_every_page_of_its_system_draws(client, mysql_config):  # noqa: F811
    system_id = _saved_with_asteroid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT m.id, m.planet_class, m.roundness_ratio, m.description FROM moons m JOIN planets p ON p.id = m.planet_id WHERE p.star_system_id = ?", (system_id,)).fetchall()
    finally:
        conn.close()
    asteroids = [row for row in rows if row["planet_class"] is None]
    assert asteroids and all(row["roundness_ratio"] < 1.0 for row in asteroids)
    assert all(row["description"] == planet_classes.ASTEROID.sentence for row in asteroids)
    assert all(row["roundness_ratio"] is None or row["roundness_ratio"] >= 1.0 for row in rows if row["planet_class"] is not None)    # a giant moon has no ratio
    assert client.get(f"/system/{pid('system', system_id)}").status_code == 200
    assert client.get(f"/system/{pid('system', system_id)}/scene").status_code == 200
    for row in asteroids:
        assert client.get(f"/object/moon:{pid('moon', row['id'])}", follow_redirects=True).status_code == 200
        assert client.get(f"/object/moon:{pid('moon', row['id'])}/tree", follow_redirects=True).status_code == 200


def test_the_summary_sentence_names_an_asteroid():
    from planetgen.web.lib import planet_summary
    assert planet_summary._class_sentence("Rock", {"body_type": "t", "roundness_ratio": 0.4, "planet_class": None}, {}) == "Rock is an asteroid."
    assert planet_summary._class_sentence("Mimas", {"body_type": "t", "roundness_ratio": 2.0, "planet_class": "D"}, {}) == "Mimas is a Class D world."
