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
from stellarObjects.doubleStar import BinaryStarProxy
from stellarObjects.spaceSector import SpaceSector
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


# ---------------------------------------------------------------------------
# Read path (_db.load_star_system / load_sector / load_system_config) --
# previously entirely untested (db/README.md's "no read path yet" status,
# now built). Covers the gaps this module's own docstring called out as
# still missing: binary systems, lifespan_gy round-tripping, and
# PRAGMA foreign_key_check, plus save_sector/insert_sector and
# insert_system_config's SLOTS child rows.
# ---------------------------------------------------------------------------

def test_load_star_system_round_trips_single_star_system_exactly(tmp_path):
    system, cfg = _make_system_with_moons_and_belt()
    db_path = str(tmp_path / "test.db")

    conn = _db.get_connection(db_path)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert str(reloaded) == str(system)
    assert reloaded.planet_count == system.planet_count
    assert reloaded.belt_count == system.belt_count
    assert reloaded.moon_count == system.moon_count
    assert reloaded.hab_count == system.hab_count
    assert reloaded.m_count == system.m_count
    # Shared back-references: the SAME star/system_config instance everywhere.
    for obj in reloaded.planets:
        assert obj.system_config is reloaded.system_config
        if obj.body_type != "a":
            assert obj.star is reloaded.star


def test_load_star_system_round_trips_binary_system_exactly(tmp_path):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.BINARY_SYSTEM = True
    cfg.PLANETS = False
    system = StarSystem(system_config=cfg)
    db_path = str(tmp_path / "test.db")

    conn = _db.get_connection(db_path)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert isinstance(reloaded.star, BinaryStarProxy)
    assert len(reloaded.stars) == 2
    assert str(reloaded) == str(system)
    assert reloaded.star.mass == system.star.mass
    assert reloaded.star.luminosity == system.star.luminosity
    # The DB schema has exactly one system_config_id per star_systems row,
    # so the secondary star's generation-time-only deep-copied config
    # (see TODO.md's resolved "secondary-star system_config asymmetry"
    # question) is naturally collapsed to the one shared config on reload.
    assert reloaded.secondary_star.system_config is reloaded.system_config
    assert reloaded.primary_star.system_config is reloaded.system_config


def test_lifespan_gy_null_round_trips_to_infinite_lifespan(tmp_path):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "M2VII"  # white dwarf: lifespan == float('inf')
    cfg.PLANETS = False
    system = StarSystem(system_config=cfg)
    assert system.star.lifespan == float('inf')
    db_path = str(tmp_path / "test.db")

    conn = _db.get_connection(db_path)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        row = conn.execute(
            "SELECT lifespan_gy FROM stars WHERE star_system_id = ?", (system_id,)
        ).fetchone()
        assert row["lifespan_gy"] is None  # NULL in storage, never the JSON "Infinity" token

        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert reloaded.star.lifespan == float('inf')
    assert str(reloaded) == str(system)


def test_foreign_key_check_is_clean_after_insert(tmp_path):
    system, cfg = _make_system_with_moons_and_belt()
    db_path = str(tmp_path / "test.db")

    conn = _db.get_connection(db_path)
    try:
        with conn:
            _db.insert_star_system(conn, system, cfg)
        violations = conn.execute("PRAGMA foreign_key_check").fetchall()
    finally:
        conn.close()

    assert violations == []


def test_save_sector_and_load_sector_round_trip(tmp_path):
    sector = SpaceSector("Persisted Sector", edge_ly=20.0)

    system_a, cfg_a = _make_system_with_moons_and_belt()
    sector.add_system(system_a, position=(1.0, 2.0, 3.0), system_config=cfg_a)

    cfg_b = SystemConfig()
    cfg_b.STAR_TYPE = "M5V"
    system_b = StarSystem(system_config=cfg_b)
    sector.add_system(system_b, position=(-4.0, 0.0, 5.5), system_config=cfg_b)

    db_path = str(tmp_path / "test.db")
    sector_id = _db.save_sector(sector, db_path=db_path)

    conn = _db.get_connection(db_path)
    try:
        reloaded_sector = _db.load_sector(conn, sector_id)
    finally:
        conn.close()

    assert reloaded_sector.name == "Persisted Sector"
    assert reloaded_sector.edge_ly == pytest.approx(20.0)
    assert len(reloaded_sector) == 2

    reloaded_by_name = {entry.star_system.star.name: entry for entry in reloaded_sector.entries}
    for original_system, original_position in ((system_a, (1.0, 2.0, 3.0)), (system_b, (-4.0, 0.0, 5.5))):
        entry = reloaded_by_name[original_system.star.name]
        assert str(entry.star_system) == str(original_system)
        # Round-tripped through ly -> milliparsecs -> ly, so exact equality
        # isn't expected -- only floating-point-close.
        assert entry.position == pytest.approx(original_position)


def test_insert_system_config_round_trips_slots_child_rows(tmp_path):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.SLOTS = [
        None,
        {"type": "planet", "planet_class": "M", "moons": 2},
        {"type": "asteroid_belt", "planet_class": None, "moons": None},
    ]
    db_path = str(tmp_path / "test.db")

    conn = _db.get_connection(db_path)
    try:
        with conn:
            config_id = _db.insert_system_config(conn, cfg)
        reloaded = _db.load_system_config(conn, config_id)
    finally:
        conn.close()

    assert reloaded.SLOTS == cfg.SLOTS
    assert reloaded.STAR_TYPE == cfg.STAR_TYPE
