# tests/test_db_persistence.py

"""
End-to-end save-path tests for `stellarObjects._db.insert_star_system`,
specifically that a freshly generated system's planets and moons land in
their respective schema-v2 tables (`planets` vs. `moons`, see
`schema.sql`'s "v2" header note) rather than sharing one table the way
schema v1 did. There's otherwise no test coverage of `_db.py`'s save path
at all (`db/README.md`'s own "no read path yet" status), so this is
deliberately a real generate-then-save-then-query round trip rather than
a unit test against a hand-built object graph -- it's the same shape of
bug (an `INSERT`'s column list and value tuple silently drifting out of
count/order) that unit tests calling `insert_moon` directly with
hand-picked arguments would be less likely to catch.
"""

import sqlite3

import pytest

from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.systemData import StarSystem


def _make_system_with_moons_and_belt():
    """Retries generation (bounded) until a system with at least one moon
    and one asteroid belt comes out -- MOONS/ASTEROID_BELT=True bias
    generation heavily toward both, but neither is unconditionally
    guaranteed for every planet/system."""
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.MOONS = True
        cfg.MAX_PLANETS = True
        cfg.ASTEROID_BELT = True
        system = StarSystem(system_config=cfg)
        planets = [obj for obj in system.planets if obj.body_type != "a"]
        belts = [obj for obj in system.planets if obj.body_type == "a"]
        if belts and any(p.moons for p in planets):
            return system, cfg
    pytest.fail("could not generate a system with both a moon and an asteroid belt")


def test_insert_star_system_splits_planets_and_moons_into_their_own_tables(tmp_path):
    system, cfg = _make_system_with_moons_and_belt()
    db_path = str(tmp_path / "test.db")

    conn = _db.get_connection(db_path)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
    finally:
        conn.close()

    expected_planets = [obj for obj in system.planets if obj.body_type != "a"]
    expected_belts = [obj for obj in system.planets if obj.body_type == "a"]
    expected_moon_names = sorted(moon.name for planet in expected_planets for moon in planet.moons)
    assert expected_moon_names, "test fixture must actually contain moons"

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        db_planets = conn.execute("SELECT * FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()
        assert len(db_planets) == len(expected_planets)
        assert "is_moon" not in db_planets[0].keys()
        assert "parent_planet_id" not in db_planets[0].keys()

        db_moons = conn.execute("SELECT * FROM moons WHERE star_system_id = ?", (system_id,)).fetchall()
        assert sorted(row["name"] for row in db_moons) == expected_moon_names

        planet_ids_by_name = {row["name"]: row["id"] for row in db_planets}
        for planet in expected_planets:
            moon_rows = [row for row in db_moons if row["planet_id"] == planet_ids_by_name[planet.name]]
            assert sorted(row["name"] for row in moon_rows) == sorted(m.name for m in planet.moons)

        belts = conn.execute("SELECT * FROM asteroid_belts WHERE star_system_id = ?", (system_id,)).fetchall()
        assert len(belts) == len(expected_belts)

        assert conn.execute("PRAGMA user_version").fetchone()[0] == _db.SCHEMA_VERSION
    finally:
        conn.close()
