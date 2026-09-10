# tests/test_db_persistence.py

"""
End-to-end save-path tests for `stellarObjects._db.insert_star_system`,
specifically that a freshly generated system's planets and moons land in
their respective schema-v2 tables (`planets` vs. `moons`, see
`schema.sql`'s "v2" header note) rather than sharing one table the way
schema v1 did. There's otherwise no test coverage of `_db.py`'s save path
at all (`docs/database-schema.md`'s own "no read path yet" status,
historical, since resolved), so this is deliberately a real
generate-then-save-then-query round trip rather than a unit test against
a hand-built object graph -- it's the same shape of bug (an `INSERT`'s
column list and value tuple silently drifting out of count/order) that
unit tests calling `insert_moon` directly with hand-picked arguments
would be less likely to catch.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
they're skipped, not failed, when no MySQL test server is configured/
reachable.
"""

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


def test_insert_star_system_splits_planets_and_moons_into_their_own_tables(mysql_config):
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
    finally:
        conn.close()

    expected_planets = [obj for obj in system.planets if obj.body_type != "a"]
    expected_belts = [obj for obj in system.planets if obj.body_type == "a"]
    expected_moon_names = sorted(moon.name for planet in expected_planets for moon in planet.moons)
    assert expected_moon_names, "test fixture must actually contain moons"

    conn = _db.get_connection(mysql_config)
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

        version_row = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()
        assert version_row["version"] == _db.SCHEMA_VERSION
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Read path (_db.load_star_system / load_sector / load_system_config).
# Covers the gaps this module's own docstring called out as still missing:
# binary systems, lifespan_gy round-tripping, and foreign-key integrity,
# plus save_sector/insert_sector and insert_system_config's SLOTS child rows.
# ---------------------------------------------------------------------------

def test_load_star_system_round_trips_single_star_system_exactly(mysql_config):
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
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


def test_load_star_system_round_trips_binary_system_exactly(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.BINARY_SYSTEM = True
    cfg.PLANETS = False
    system = StarSystem(system_config=cfg)

    conn = _db.get_connection(mysql_config)
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


def test_lifespan_gy_null_round_trips_to_infinite_lifespan(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "M2VII"  # white dwarf: lifespan == float('inf')
    cfg.PLANETS = False
    system = StarSystem(system_config=cfg)
    assert system.star.lifespan == float('inf')

    conn = _db.get_connection(mysql_config)
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


def test_insert_star_system_respects_foreign_keys(mysql_config):
    # MySQL/InnoDB enforces every foreign key eagerly, at INSERT time --
    # unlike SQLite (which needs a separate `PRAGMA foreign_key_check`
    # pass after the fact to catch a constraint left unenforced mid-
    # transaction), an insert violating a foreign key here would already
    # have raised `pymysql.err.IntegrityError` above rather than
    # completing silently. This test's real assertion is simply that the
    # insert-then-reload round trip above completes without one.
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_star_system(conn, system, cfg)
    finally:
        conn.close()


def test_save_sector_and_load_sector_round_trip(mysql_config):
    sector = SpaceSector("Persisted Sector", edge_ly=20.0)

    system_a, cfg_a = _make_system_with_moons_and_belt()
    sector.add_system(system_a, position=(1.0, 2.0, 3.0), system_config=cfg_a)

    cfg_b = SystemConfig()
    cfg_b.STAR_TYPE = "M5V"
    system_b = StarSystem(system_config=cfg_b)
    sector.add_system(system_b, position=(-4.0, 0.0, 5.5), system_config=cfg_b)

    sector_id = _db.save_sector(sector, config=mysql_config)

    conn = _db.get_connection(mysql_config)
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


def test_orbital_motion_fields_round_trip_exactly(mysql_config):
    """
    Regression test for the exact bug class this module's own docstring
    warns about (an `INSERT`'s column list and value tuple silently
    drifting out of count): `insert_moon`'s `VALUES` clause was originally
    missing one placeholder relative to its column list when the v9
    orbital-motion columns were added, which pymysql surfaced as a generic
    "not all arguments converted during string formatting" `TypeError` --
    a symptom easy to mistake for something else. Confirms every
    orbital-motion field survives a save/load round trip exactly, for both
    planets and moons.
    """
    system, cfg = _make_system_with_moons_and_belt()
    planets = [obj for obj in system.planets if obj.body_type != "a"]
    assert any(p.moons for p in planets), "test fixture must actually contain moons"

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    reloaded_planets_by_name = {p.name: p for p in reloaded.planets if p.body_type != "a"}
    for planet in planets:
        reloaded_planet = reloaded_planets_by_name[planet.name]
        assert reloaded_planet.orbital_inclination_deg == pytest.approx(planet.orbital_inclination_deg)
        assert reloaded_planet.orbital_ascending_node_deg == pytest.approx(planet.orbital_ascending_node_deg)
        assert reloaded_planet.orbital_phase_deg == pytest.approx(planet.orbital_phase_deg)
        assert reloaded_planet.rotation_period_hours == pytest.approx(planet.rotation_period_hours)

        reloaded_moons_by_name = {m.name: m for m in reloaded_planet.moons}
        for moon in planet.moons:
            reloaded_moon = reloaded_moons_by_name[moon.name]
            assert reloaded_moon.orbital_inclination_deg == pytest.approx(moon.orbital_inclination_deg)
            assert reloaded_moon.orbital_ascending_node_deg == pytest.approx(moon.orbital_ascending_node_deg)
            assert reloaded_moon.orbital_phase_deg == pytest.approx(moon.orbital_phase_deg)
            assert reloaded_moon.rotation_period_hours == pytest.approx(moon.rotation_period_hours)


# ---------------------------------------------------------------------------
# Orbital motion updates (stellarObjects._db.advance_orbital_phases /
# get_orbit_update_elapsed_years) -- see updateOrbits.py.
# ---------------------------------------------------------------------------

def test_get_orbit_update_elapsed_years_is_none_before_first_update(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        assert _db.get_orbit_update_elapsed_years(conn) is None
    finally:
        conn.close()


def test_advance_orbital_phases_applies_the_correct_delta_and_leaves_other_fields_untouched(mysql_config):
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_star_system(conn, system, cfg)

        before = conn.execute(
            "SELECT id, orbital_phase_deg, orbital_inclination_deg, orbital_ascending_node_deg, "
            "rotation_period_hours, period_years FROM planets"
        ).fetchall()

        planets_updated, moons_updated = _db.advance_orbital_phases(conn, elapsed_years=0.5)
        assert planets_updated > 0
        assert moons_updated > 0

        after = conn.execute(
            "SELECT id, orbital_phase_deg, orbital_inclination_deg, orbital_ascending_node_deg, "
            "rotation_period_hours, period_years FROM planets"
        ).fetchall()
    finally:
        conn.close()

    after_by_id = {row["id"]: row for row in after}
    for row in before:
        updated = after_by_id[row["id"]]
        # Orientation/rotation are fixed at generation time -- only phase moves.
        assert updated["orbital_inclination_deg"] == pytest.approx(row["orbital_inclination_deg"])
        assert updated["orbital_ascending_node_deg"] == pytest.approx(row["orbital_ascending_node_deg"])
        assert updated["rotation_period_hours"] == pytest.approx(row["rotation_period_hours"])

        expected_phase = (row["orbital_phase_deg"] + (0.5 / row["period_years"]) * 360) % 360
        assert updated["orbital_phase_deg"] == pytest.approx(expected_phase, abs=1e-6)

    # get_orbit_update_elapsed_years should now report ~0 elapsed time
    # (the call above just set last_updated_at to NOW()), not None.
    conn = _db.get_connection(mysql_config)
    try:
        elapsed = _db.get_orbit_update_elapsed_years(conn)
    finally:
        conn.close()
    assert elapsed is not None
    assert 0 <= elapsed < 0.01


def test_advance_orbital_phases_rejects_negative_elapsed_years(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError):
            _db.advance_orbital_phases(conn, elapsed_years=-1.0)
    finally:
        conn.close()


def test_migrate_v8_to_v9_adds_orbital_motion_columns(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v9) -- so this simulates an existing v8
    database by tearing the v9 additions back out (the new columns, the
    new `orbit_simulation_state` table, and the `schema_migrations` v9
    row) before calling `migrate_database`, and confirms it puts them
    back and reports the database as current again.
    """
    conn = _db.get_connection(mysql_config)
    try:
        for table in ("planets", "moons"):
            conn.execute(
                f"ALTER TABLE {table} "
                f"DROP COLUMN orbital_inclination_deg, DROP COLUMN orbital_ascending_node_deg, "
                f"DROP COLUMN orbital_phase_deg, DROP COLUMN rotation_period_hours"
            )
        conn.execute("DROP TABLE orbit_simulation_state")
        conn.execute("DELETE FROM schema_migrations WHERE version = 9")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (8)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 8
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION == 9

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM planets").fetchall()}
        assert {
            "orbital_inclination_deg", "orbital_ascending_node_deg",
            "orbital_phase_deg", "rotation_period_hours",
        } <= columns
        assert _db.get_orbit_update_elapsed_years(conn) is None  # table exists, no row yet
    finally:
        conn.close()

    # Idempotent: running it again against an already-v9 database is a no-op.
    assert _db.migrate_database(mysql_config) == 9


def test_insert_system_config_round_trips_slots_child_rows(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.SLOTS = [
        None,
        {"type": "planet", "planet_class": "M", "moons": 2},
        {"type": "asteroid_belt", "planet_class": None, "moons": None},
    ]

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            config_id = _db.insert_system_config(conn, cfg)
        reloaded = _db.load_system_config(conn, config_id)
    finally:
        conn.close()

    assert reloaded.SLOTS == cfg.SLOTS
    assert reloaded.STAR_TYPE == cfg.STAR_TYPE
