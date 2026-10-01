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

import datetime
import math
import time

import pytest

from stellarObjects import _db, program_constants
from stellarObjects import physical_constants as pc
from stellarObjects import program_constants
from stellarObjects.cometData import Comet
from stellarObjects.compactRemnant import BlackHole, NeutronStar
from stellarObjects.config import SystemConfig
from stellarObjects.doubleStar import BinaryStarProxy
from stellarObjects.planetPhysics import calculate_orbital_period_years
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem
from stellarObjects.utils import circular_orbital_speed_kms, minimum_update_interval_years, orbital_position_au


def _drop_v17_phenomenon_columns(conn):
    """
    Drops the v17 `galactic_orbital_*` columns from the six pre-existing
    exotic-phenomenon tables (`black_holes`, `neutron_stars`, `nebulae`,
    `supernova_remnants`, `rogue_planets`, `interstellar_comets`) -- shared
    by every `test_migrate_vN_to_vN+1_*` test below that simulates a
    database older than v16 (the version these tables were introduced in).
    `mysql_config` always bootstraps a fresh database at the CURRENT
    schema (today, v17), so those tables and their v17 columns already
    exist even though a real pre-v16 database would have neither -- left
    undropped, `migrate_database`'s replayed `_migrate_v16_to_v17` would
    try to `ADD COLUMN` ones that already exist and fail with a duplicate-
    column error. See `schema.sql`'s "v17" header note.
    """
    for table in ("black_holes", "neutron_stars", "nebulae", "supernova_remnants",
                  "rogue_planets", "interstellar_comets"):
        conn.execute(
            f"ALTER TABLE {table} "
            "DROP COLUMN galactic_orbital_speed_kms, DROP COLUMN galactic_orbital_period_gy, "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years"
        )


def _drop_v18_phenomenon_columns(conn):
    """
    Drops the v18 galaxy-frame placement columns from `nebulae`/
    `asteroid_fields` -- same "already exists on a freshly-bootstrapped
    test database, so migrate_database's replayed _migrate_v17_to_v18
    would hit a duplicate-column error" reasoning
    `_drop_v17_phenomenon_columns` gives for its own six tables. See
    `schema.sql`'s "v18" header note.

    The named `chk_{table}_placement` CHECK constraint must be dropped
    first (`DROP CONSTRAINT`, portable across MySQL 8.0.19+/MariaDB --
    MySQL's own `DROP CHECK` alias is not): MySQL refuses `DROP COLUMN` on
    a column a CHECK still references (error 3959), and a fresh test
    database already has this constraint (`schema.sql`'s own
    `CREATE TABLE` body), unlike a real pre-v18 database, which never had
    it either.
    """
    for table in ("nebulae", "asteroid_fields"):
        conn.execute(
            f"ALTER TABLE {table} "
            f"DROP CONSTRAINT chk_{table}_placement, "
            f"DROP INDEX idx_{table}_galactic_radius_pc, "
            "DROP COLUMN center_x_pc, DROP COLUMN center_y_pc, "
            "DROP COLUMN center_z_pc, DROP COLUMN galactic_radius_pc"
        )


def _drop_v20_trajectory_columns(conn):
    """
    Drops the v20 proper-two-body/barycentric-trajectory columns from the
    three pre-existing tables they were added to (`stars`, `planets`,
    `star_systems`) -- the same "shared by every `test_migrate_vN_to_vN+1_*`
    test below that simulates a database older than v20" reasoning
    `_drop_v17_phenomenon_columns` gives for its own six tables. See
    `schema.sql`'s "v20" header note.
    """
    for table in ("stars", "planets"):
        conn.execute(
            f"ALTER TABLE {table} "
            "DROP COLUMN reflex_offset_x_km, DROP COLUMN reflex_offset_y_km, "
            "DROP COLUMN reflex_offset_z_km"
        )
    conn.execute(
        "ALTER TABLE star_systems "
        "DROP COLUMN binary_primary_position_x_km, DROP COLUMN binary_primary_position_y_km, "
        "DROP COLUMN binary_primary_position_z_km, "
        "DROP COLUMN binary_secondary_position_x_km, DROP COLUMN binary_secondary_position_y_km, "
        "DROP COLUMN binary_secondary_position_z_km, "
        "DROP COLUMN binary_secondary_mass_fraction, "
        "DROP COLUMN binary_planetary_wobble_x_km, DROP COLUMN binary_planetary_wobble_y_km, "
        "DROP COLUMN binary_planetary_wobble_z_km"
    )


def _drop_v21_phenomenon_columns(conn):
    """
    Drops the v21 sector-placement columns (`sector_id`, `center_x/y/z_pc`,
    `galactic_radius_pc`) from `black_holes`/`neutron_stars` -- the same
    "already exists on a freshly-bootstrapped test database, so
    migrate_database's replayed `_migrate_v20_to_v21` would hit a
    duplicate-column error" reasoning `_drop_v17_phenomenon_columns`/
    `_drop_v18_phenomenon_columns` give for their own tables. See
    `schema.sql`'s "v21" header note.

    The named `chk_{table}_placement` CHECK constraint (and the `sector_id`
    FK) must be dropped first -- same ordering reasoning
    `_drop_v18_phenomenon_columns` gives for its own identical constraint.
    """
    for table in ("black_holes", "neutron_stars"):
        conn.execute(
            f"ALTER TABLE {table} "
            f"DROP CONSTRAINT chk_{table}_placement, "
            f"DROP FOREIGN KEY fk_{table}_sector, "
            f"DROP INDEX idx_{table}_sector_id, "
            f"DROP INDEX idx_{table}_galactic_radius_pc, "
            "DROP COLUMN sector_id, DROP COLUMN center_x_pc, DROP COLUMN center_y_pc, "
            "DROP COLUMN center_z_pc, DROP COLUMN galactic_radius_pc"
        )


def _drop_v22_search_indexes(conn):
    """
    Drops the v22 search-facing indexes `_migrate_v21_to_v22` adds -- same
    "already exists on a freshly-bootstrapped test database, so replaying
    the migration would hit a duplicate-key error" reasoning every other
    `_drop_vNN_...` helper here gives for its own table. See `schema.sql`'s
    "v22" header note.
    """
    conn.execute("ALTER TABLE sectors DROP INDEX idx_sectors_name")
    conn.execute("ALTER TABLE star_systems DROP INDEX idx_star_systems_name")
    conn.execute(
        "ALTER TABLE stars "
        "DROP INDEX idx_stars_name, DROP INDEX idx_stars_yerkes_class"
    )
    conn.execute(
        "ALTER TABLE planets "
        "DROP INDEX idx_planets_name, DROP INDEX idx_planets_planet_class, "
        "DROP INDEX idx_planets_body_type, DROP INDEX idx_planets_life_chemical"
    )
    conn.execute(
        "ALTER TABLE moons "
        "DROP INDEX idx_moons_name, DROP INDEX idx_moons_planet_class, "
        "DROP INDEX idx_moons_body_type, DROP INDEX idx_moons_life_chemical"
    )
    conn.execute("ALTER TABLE asteroid_belts DROP INDEX idx_asteroid_belts_density")


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
        # Pinned False: left at its own default, `BINARY_SYSTEM` now rolls
        # real chance (StarSystem._should_generate_binary), and this
        # helper's only caller counts expected rows from `system.planets`
        # alone -- an S-type (wide) pair's own `secondary_planets` would
        # add real DB rows that count wouldn't account for.
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        planets = [obj for obj in system.planets if obj.body_type != "a"]
        belts = [obj for obj in system.planets if obj.body_type == "a"]
        if belts and any(p.moons for p in planets):
            return system, cfg
    pytest.fail("could not generate a system with both a moon and an asteroid belt")


def _make_system_with_comets():
    """Retries generation (bounded) until a system with at least one
    elliptical AND at least one parabolic comet comes out -- COMETS=True
    biases generation heavily toward having some, but not toward any
    particular orbit_type mix (see program_constants.COMET_PARABOLIC_CHANCE)."""
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.COMETS = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        orbit_types = {c.orbit_type for c in system.comets}
        if {"elliptical", "parabolic"} <= orbit_types:
            return system, cfg
    pytest.fail("could not generate a system with both an elliptical and a parabolic comet")


def test_insert_star_system_persists_comets_in_their_own_table(mysql_config):
    system, cfg = _make_system_with_comets()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)

        db_comets = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchall()
        assert len(db_comets) == system.comet_count

        star_row = conn.execute(
            "SELECT id FROM stars WHERE star_system_id = ? AND role = 'single'", (system_id,)
        ).fetchone()

        db_comets_by_name = {row["name"]: row for row in db_comets}
        for comet in system.comets:
            row = db_comets_by_name[comet.name]
            assert row["orbit_type"] == comet.orbit_type
            # A single star DOES get a real stars.id (unlike a 'close'
            # binary's merged proxy, which has none) -- see insert_comet's
            # docstring and insert_star_system's own primary_star_id handling.
            assert row["star_id"] == star_row["id"]
            assert math.isclose(row["perihelion_distance_km"], comet.perihelion_distance_au * pc.AU_TO_KM, rel_tol=1e-9)
            assert math.isclose(row["eccentricity"], comet.eccentricity, rel_tol=1e-9)
            assert row["composition_summary"] == comet.get_composition_summary()

            comp_rows = conn.execute(
                "SELECT component FROM comet_composition WHERE comet_id = ? ORDER BY position", (row["id"],)
            ).fetchall()
            assert [r["component"] for r in comp_rows] == list(comet.composition)

            if comet.orbit_type == "elliptical":
                assert row["orbital_period_years"] is not None
                assert row["mean_anomaly_deg"] is not None
                assert row["parabolic_mean_anomaly"] is None
            else:
                assert row["orbital_period_years"] is None
                assert row["mean_anomaly_deg"] is None
                assert row["parabolic_mean_anomaly"] is not None
    finally:
        conn.close()


def test_load_star_system_round_trips_comets_with_composition(mysql_config):
    system, cfg = _make_system_with_comets()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert reloaded.comet_count == system.comet_count
    original_by_name = {c.name: c for c in system.comets}
    for comet in reloaded.comets:
        original = original_by_name[comet.name]
        assert comet.orbit_type == original.orbit_type
        assert comet.period_class == original.period_class
        assert math.isclose(comet.perihelion_distance_au, original.perihelion_distance_au, rel_tol=1e-9)
        assert math.isclose(comet.eccentricity, original.eccentricity, rel_tol=1e-9)
        assert math.isclose(comet.distance_au, original.distance_au, rel_tol=1e-9)
        assert math.isclose(comet.orbital_speed_kms, original.orbital_speed_kms, rel_tol=1e-9)
        assert comet.composition == original.composition
        assert comet.system_config is reloaded.system_config
    assert str(reloaded) == str(system)


def _make_wide_binary_system_with_comets():
    """Retries generation (bounded) until a 'wide' (S-type) binary comes
    out with at least one comet on EACH star -- needed to exercise
    insert_star_system's/load_star_system's star_id-based grouping for
    comets (see insert_comet's and load_star_system's own docstrings),
    which _make_system_with_comets' single-star fixture above can't."""
    for _ in range(50):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.COMETS = True
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = True
        system = StarSystem(system_config=cfg)
        if system.binary_type == "wide" and system.comets and system.secondary_comets:
            return system, cfg
    pytest.fail("could not generate a wide binary with comets on both stars")


def test_insert_and_load_star_system_round_trips_comets_for_a_wide_binary(mysql_config):
    """
    A 'wide' binary's two stars each get their own, independently-rolled
    comet population (StarSystem._generate_comets), disambiguated in the
    database by each row's own star_id (see insert_comet's docstring) --
    unlike the single-star case test_insert_star_system_persists_comets_in_their_own_table
    and test_load_star_system_round_trips_comets_with_composition cover
    above, this is the only place that star_id-based split is actually
    exercised through a real insert/load round trip, rather than just the
    in-memory to_dict/from_dict path StarSystem's own serialization tests
    (test_serialization.py) cover.
    """
    system, cfg = _make_wide_binary_system_with_comets()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)

        primary_row = conn.execute(
            "SELECT id FROM stars WHERE star_system_id = ? AND role = 'primary'", (system_id,)
        ).fetchone()
        secondary_row = conn.execute(
            "SELECT id FROM stars WHERE star_system_id = ? AND role = 'secondary'", (system_id,)
        ).fetchone()

        db_comets = conn.execute("SELECT name, star_id FROM comets WHERE star_system_id = ?", (system_id,)).fetchall()
        primary_names = {row["name"] for row in db_comets if row["star_id"] == primary_row["id"]}
        secondary_names = {row["name"] for row in db_comets if row["star_id"] == secondary_row["id"]}
        assert primary_names == {c.name for c in system.comets}
        assert secondary_names == {c.name for c in system.secondary_comets}
        # Every row's star_id must land in exactly one of the two buckets
        # above -- confirms there's no third value (e.g. a stray NULL,
        # only ever correct for a 'close' binary's merged proxy) hiding.
        assert primary_names | secondary_names == {row["name"] for row in db_comets}

        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert {c.name for c in reloaded.comets} == {c.name for c in system.comets}
    assert {c.name for c in reloaded.secondary_comets} == {c.name for c in system.secondary_comets}
    assert reloaded.comet_count == system.comet_count


def _make_close_binary_system_with_comets():
    """Retries generation (bounded) until a 'close' (P-type) binary comes
    out with at least one comet -- a close pair's comets are bound to the
    merged BinaryStarProxy (StarSystem._generate_comets), not either
    individually-stored star row, so insert_comet stores star_id as NULL
    for them (see that function's own docstring); untested by any other
    comet fixture above, all of which use a real star_id."""
    for _ in range(50):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.COMETS = True
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = False
        system = StarSystem(system_config=cfg)
        if system.binary_type == "close" and system.comets:
            return system, cfg
    pytest.fail("could not generate a close binary with comets")


def test_insert_and_load_star_system_round_trips_comets_for_a_close_binary(mysql_config):
    system, cfg = _make_close_binary_system_with_comets()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)

        db_comets = conn.execute("SELECT name, star_id FROM comets WHERE star_system_id = ?", (system_id,)).fetchall()
        assert len(db_comets) == system.comet_count
        assert all(row["star_id"] is None for row in db_comets)

        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert {c.name for c in reloaded.comets} == {c.name for c in system.comets}
    assert reloaded.secondary_comets == []
    assert reloaded.comet_count == system.comet_count


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
    cfg.WIDE_BINARY = False  # this test specifically asserts BinaryStarProxy-only behavior
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


def test_galactic_orbit_fields_round_trip_exactly(mysql_config):
    """
    Regression test for the v10 `galactic_orbital_speed_kms`/
    `galactic_orbital_period_gy` columns, for both a single star and a
    binary pair -- the same "INSERT column list vs. value tuple drift"
    bug class `test_orbital_motion_fields_round_trip_exactly` guards
    against for the v9 columns.
    """
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert reloaded.star.galactic_orbital_speed_kms == pytest.approx(system.star.galactic_orbital_speed_kms)
    assert reloaded.star.galactic_orbital_period_gy == pytest.approx(system.star.galactic_orbital_period_gy)

    binary_cfg = SystemConfig()
    binary_cfg.STAR_TYPE = "G2V"
    binary_cfg.BINARY_SYSTEM = True
    binary_cfg.WIDE_BINARY = False  # pin to close/P-type -- these assertions are BinaryStarProxy-specific
    binary_cfg.PLANETS = False
    binary_system = StarSystem(system_config=binary_cfg)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            binary_system_id = _db.insert_star_system(conn, binary_system, binary_cfg)
        reloaded_binary = _db.load_star_system(conn, binary_system_id)
    finally:
        conn.close()

    assert reloaded_binary.star.galactic_orbital_speed_kms == pytest.approx(binary_system.star.galactic_orbital_speed_kms)
    assert reloaded_binary.star.galactic_orbital_period_gy == pytest.approx(binary_system.star.galactic_orbital_period_gy)


def test_star_motion_fields_round_trip_exactly(mysql_config):
    """
    Regression test for the v13 star-motion columns and the v14 binary
    mutual-orbit position columns, for both a single star and a binary
    pair -- the same "INSERT column list vs. value tuple drift" bug class
    `test_orbital_motion_fields_round_trip_exactly`/
    `test_galactic_orbit_fields_round_trip_exactly` guard against for the
    v9/v10 columns.
    """
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert reloaded.star.galactic_orbital_phase_deg == pytest.approx(system.star.galactic_orbital_phase_deg)
    assert reloaded.star.galactic_min_update_interval_years == pytest.approx(
        system.star.galactic_min_update_interval_years
    )

    binary_cfg = SystemConfig()
    binary_cfg.STAR_TYPE = "G2V"
    binary_cfg.BINARY_SYSTEM = True
    binary_cfg.WIDE_BINARY = False  # pin to close/P-type -- these assertions are BinaryStarProxy-specific
    binary_cfg.PLANETS = False
    binary_system = StarSystem(system_config=binary_cfg)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            binary_system_id = _db.insert_star_system(conn, binary_system, binary_cfg)
        reloaded_binary = _db.load_star_system(conn, binary_system_id)
    finally:
        conn.close()

    proxy = binary_system.star
    reloaded_proxy = reloaded_binary.star
    assert reloaded_proxy.galactic_orbital_phase_deg == pytest.approx(proxy.galactic_orbital_phase_deg)
    assert reloaded_proxy.galactic_min_update_interval_years == pytest.approx(
        proxy.galactic_min_update_interval_years
    )
    # Both stars of a binary pair, and the proxy, always share one galactic phase.
    assert reloaded_binary.primary_star.galactic_orbital_phase_deg == pytest.approx(proxy.galactic_orbital_phase_deg)
    assert reloaded_binary.secondary_star.galactic_orbital_phase_deg == pytest.approx(
        proxy.galactic_orbital_phase_deg
    )

    assert reloaded_proxy.binary_mutual_orbital_period_years == pytest.approx(
        proxy.binary_mutual_orbital_period_years
    )
    assert reloaded_proxy.binary_mutual_orbital_speed_kms == pytest.approx(proxy.binary_mutual_orbital_speed_kms)
    assert reloaded_proxy.binary_mutual_orbital_inclination_deg == pytest.approx(
        proxy.binary_mutual_orbital_inclination_deg
    )
    assert reloaded_proxy.binary_mutual_orbital_ascending_node_deg == pytest.approx(
        proxy.binary_mutual_orbital_ascending_node_deg
    )
    assert reloaded_proxy.binary_mutual_orbital_phase_deg == pytest.approx(proxy.binary_mutual_orbital_phase_deg)
    assert reloaded_proxy.binary_mutual_min_update_interval_years == pytest.approx(
        proxy.binary_mutual_min_update_interval_years
    )
    assert reloaded_proxy.binary_mutual_position_x == pytest.approx(proxy.binary_mutual_position_x)
    assert reloaded_proxy.binary_mutual_position_y == pytest.approx(proxy.binary_mutual_position_y)
    assert reloaded_proxy.binary_mutual_position_z == pytest.approx(proxy.binary_mutual_position_z)

    # The stored position must actually match orbital_position_au at the
    # stored separation/inclination/ascending_node/phase, not just survive
    # a round trip unchanged.
    expected = orbital_position_au(
        proxy.binary_separation_au, proxy.binary_mutual_orbital_inclination_deg,
        proxy.binary_mutual_orbital_ascending_node_deg, proxy.binary_mutual_orbital_phase_deg,
    )
    assert (proxy.binary_mutual_position_x, proxy.binary_mutual_position_y, proxy.binary_mutual_position_z) == \
        pytest.approx(expected)


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

    reloaded_by_name = {entry.star_system.name: entry for entry in reloaded_sector.entries}
    for original_system, original_position in ((system_a, (1.0, 2.0, 3.0)), (system_b, (-4.0, 0.0, 5.5))):
        entry = reloaded_by_name[original_system.name]
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

    # Matched positionally, not by name: `insert_moon`/`insert_planet` persist
    # each body's list position via `orbital_index`, and `load_star_system`
    # loads both planets and moons `ORDER BY orbital_index`, so `reloaded`
    # preserves the original ordering exactly. Name-keyed dicts would break
    # here -- moon (and planet) names come from `generate_phoneme_salad_name`,
    # which picks names independently per body with no check against its own
    # siblings, so two moons on the same planet can occasionally collide on
    # the same generated name and silently overwrite each other in a dict.
    reloaded_planets = [p for p in reloaded.planets if p.body_type != "a"]
    assert len(reloaded_planets) == len(planets)
    for planet, reloaded_planet in zip(planets, reloaded_planets):
        assert reloaded_planet.orbital_inclination_deg == pytest.approx(planet.orbital_inclination_deg)
        assert reloaded_planet.orbital_ascending_node_deg == pytest.approx(planet.orbital_ascending_node_deg)
        assert reloaded_planet.orbital_phase_deg == pytest.approx(planet.orbital_phase_deg)
        assert reloaded_planet.position_x == pytest.approx(planet.position_x)
        assert reloaded_planet.position_y == pytest.approx(planet.position_y)
        assert reloaded_planet.position_z == pytest.approx(planet.position_z)
        assert reloaded_planet.orbital_speed_kms == pytest.approx(planet.orbital_speed_kms)
        assert reloaded_planet.min_update_interval_years == pytest.approx(planet.min_update_interval_years)
        assert reloaded_planet.rotation_period_hours == pytest.approx(planet.rotation_period_hours)

        assert len(reloaded_planet.moons) == len(planet.moons)
        for moon, reloaded_moon in zip(planet.moons, reloaded_planet.moons):
            assert reloaded_moon.orbital_inclination_deg == pytest.approx(moon.orbital_inclination_deg)
            assert reloaded_moon.orbital_ascending_node_deg == pytest.approx(moon.orbital_ascending_node_deg)
            assert reloaded_moon.orbital_phase_deg == pytest.approx(moon.orbital_phase_deg)
            assert reloaded_moon.position_x == pytest.approx(moon.position_x)
            assert reloaded_moon.position_y == pytest.approx(moon.position_y)
            assert reloaded_moon.position_z == pytest.approx(moon.position_z)
            assert reloaded_moon.orbital_speed_kms == pytest.approx(moon.orbital_speed_kms)
            assert reloaded_moon.min_update_interval_years == pytest.approx(moon.min_update_interval_years)
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
            "rotation_period_hours, period_years, distance_km, "
            "position_x_km, position_y_km, position_z_km, orbital_speed_kms FROM planets"
        ).fetchall()

        counts = _db.advance_orbital_phases(conn, elapsed_years=0.5)
        assert counts["planets"] > 0
        assert counts["moons"] > 0
        assert counts["stars"] > 0

        after = conn.execute(
            "SELECT id, orbital_phase_deg, orbital_inclination_deg, orbital_ascending_node_deg, "
            "rotation_period_hours, period_years, distance_km, "
            "position_x_km, position_y_km, position_z_km, orbital_speed_kms FROM planets"
        ).fetchall()
    finally:
        conn.close()

    after_by_id = {row["id"]: row for row in after}
    for row in before:
        updated = after_by_id[row["id"]]
        # Orientation/rotation are fixed at generation time -- only phase
        # (and position, in lockstep with it) moves. Speed is constant
        # around a circular orbit.
        assert updated["orbital_inclination_deg"] == pytest.approx(row["orbital_inclination_deg"])
        assert updated["orbital_ascending_node_deg"] == pytest.approx(row["orbital_ascending_node_deg"])
        assert updated["rotation_period_hours"] == pytest.approx(row["rotation_period_hours"])
        assert updated["orbital_speed_kms"] == pytest.approx(row["orbital_speed_kms"])

        expected_phase = (row["orbital_phase_deg"] + (0.5 / row["period_years"]) * 360) % 360
        assert updated["orbital_phase_deg"] == pytest.approx(expected_phase, abs=1e-6)

        expected_x_au, expected_y_au, expected_z_au = orbital_position_au(
            row["distance_km"] / pc.AU_TO_KM, row["orbital_inclination_deg"],
            row["orbital_ascending_node_deg"], expected_phase,
        )
        assert updated["position_x_km"] == pytest.approx(expected_x_au * pc.AU_TO_KM, abs=1e-3)
        assert updated["position_y_km"] == pytest.approx(expected_y_au * pc.AU_TO_KM, abs=1e-3)
        assert updated["position_z_km"] == pytest.approx(expected_z_au * pc.AU_TO_KM, abs=1e-3)

    # get_orbit_update_elapsed_years should now report ~0 elapsed time
    # (the call above just set last_updated_at to NOW()), not None.
    conn = _db.get_connection(mysql_config)
    try:
        elapsed = _db.get_orbit_update_elapsed_years(conn)
    finally:
        conn.close()
    assert elapsed is not None
    assert 0 <= elapsed < 0.01


def test_advance_orbital_phases_advances_binary_mutual_orbit_position_in_lockstep(mysql_config):
    """
    Regression test for the v14 addition: `advance_orbital_phases` must
    move `binary_mutual_orbital_phase_deg` *and* recompute
    `binary_mutual_position_x/y/z_km` from that new phase in the same
    pass -- the same "position has no independent update of its own"
    treatment `test_advance_orbital_phases_applies_the_correct_delta_and_leaves_other_fields_untouched`
    already verifies for planets/moons.
    """
    binary_cfg = SystemConfig()
    binary_cfg.STAR_TYPE = "G2V"
    binary_cfg.BINARY_SYSTEM = True
    binary_cfg.WIDE_BINARY = False  # pin to close/P-type -- these assertions are BinaryStarProxy-specific
    binary_cfg.PLANETS = False
    binary_system = StarSystem(system_config=binary_cfg)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            binary_id = _db.insert_star_system(conn, binary_system, binary_cfg)

        before = conn.execute(
            "SELECT binary_separation_km, binary_mutual_orbital_period_years, "
            "binary_mutual_orbital_inclination_deg, binary_mutual_orbital_ascending_node_deg, "
            "binary_mutual_orbital_phase_deg "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()

        # Big enough elapsed time to clearly move the (fast) mutual orbit,
        # regardless of its randomly generated period.
        elapsed_years = before["binary_mutual_orbital_period_years"] * 137.25
        _db.advance_orbital_phases(conn, elapsed_years=elapsed_years)

        after = conn.execute(
            "SELECT binary_mutual_orbital_phase_deg, binary_mutual_position_x_km, "
            "binary_mutual_position_y_km, binary_mutual_position_z_km "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()
    finally:
        conn.close()

    expected_phase = (
        before["binary_mutual_orbital_phase_deg"]
        + (elapsed_years / before["binary_mutual_orbital_period_years"]) * 360
    ) % 360
    assert after["binary_mutual_orbital_phase_deg"] == pytest.approx(expected_phase, abs=1e-6)
    assert after["binary_mutual_orbital_phase_deg"] != pytest.approx(before["binary_mutual_orbital_phase_deg"])

    expected_x_au, expected_y_au, expected_z_au = orbital_position_au(
        before["binary_separation_km"] / pc.AU_TO_KM, before["binary_mutual_orbital_inclination_deg"],
        before["binary_mutual_orbital_ascending_node_deg"], expected_phase,
    )
    assert after["binary_mutual_position_x_km"] == pytest.approx(expected_x_au * pc.AU_TO_KM, abs=1e-3)
    assert after["binary_mutual_position_y_km"] == pytest.approx(expected_y_au * pc.AU_TO_KM, abs=1e-3)
    assert after["binary_mutual_position_z_km"] == pytest.approx(expected_z_au * pc.AU_TO_KM, abs=1e-3)


def test_advance_orbital_phases_advances_both_binary_members_barycenter_offsets(mysql_config):
    """
    Regression test for the v18 addition: both `binary_primary_position_*`
    and `binary_secondary_position_*` must move together with the freshly-
    advanced `binary_mutual_position_*` -- not just one star sitting fixed
    while the other orbits it -- and `secondary_position == primary_position
    + mutual_position` must keep holding after every advance, not just at
    generation time.
    """
    binary_cfg = SystemConfig()
    binary_cfg.STAR_TYPE = "G2V"
    binary_cfg.BINARY_SYSTEM = True
    binary_cfg.WIDE_BINARY = False
    binary_cfg.PLANETS = False
    binary_system = StarSystem(system_config=binary_cfg)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            binary_id = _db.insert_star_system(conn, binary_system, binary_cfg)

        before = conn.execute(
            "SELECT binary_mutual_orbital_period_years, binary_primary_position_x_km, "
            "binary_secondary_position_x_km, binary_secondary_mass_fraction "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()

        elapsed_years = before["binary_mutual_orbital_period_years"] * 137.25
        counts = _db.advance_orbital_phases(conn, elapsed_years=elapsed_years)
        assert counts["binary_mutual_orbits"] > 0

        after = conn.execute(
            "SELECT binary_mutual_position_x_km, binary_mutual_position_y_km, binary_mutual_position_z_km, "
            "binary_primary_position_x_km, binary_primary_position_y_km, binary_primary_position_z_km, "
            "binary_secondary_position_x_km, binary_secondary_position_y_km, binary_secondary_position_z_km "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()
    finally:
        conn.close()

    # Both members actually moved -- neither is the trivial "stayed put" case.
    assert after["binary_primary_position_x_km"] != pytest.approx(before["binary_primary_position_x_km"])
    assert after["binary_secondary_position_x_km"] != pytest.approx(before["binary_secondary_position_x_km"])

    mu = before["binary_secondary_mass_fraction"]
    assert after["binary_primary_position_x_km"] == pytest.approx(-mu * after["binary_mutual_position_x_km"], abs=1e-3)
    assert after["binary_primary_position_y_km"] == pytest.approx(-mu * after["binary_mutual_position_y_km"], abs=1e-3)
    assert after["binary_primary_position_z_km"] == pytest.approx(-mu * after["binary_mutual_position_z_km"], abs=1e-3)
    assert after["binary_secondary_position_x_km"] == pytest.approx(
        after["binary_primary_position_x_km"] + after["binary_mutual_position_x_km"], abs=1e-3
    )
    assert after["binary_secondary_position_y_km"] == pytest.approx(
        after["binary_primary_position_y_km"] + after["binary_mutual_position_y_km"], abs=1e-3
    )
    assert after["binary_secondary_position_z_km"] == pytest.approx(
        after["binary_primary_position_z_km"] + after["binary_mutual_position_z_km"], abs=1e-3
    )


def test_advance_orbital_phases_recomputes_star_and_planet_reflex_offsets(mysql_config):
    """
    Regression test for the v18 addition: `stars.reflex_offset_*_km`/
    `planets.reflex_offset_*_km` must be recomputed from each hosted
    body's freshly-advanced position on every call (they have no
    independent guard interval of their own -- see
    `advance_orbital_phases`' own docstring), not just set once at
    generation time.
    """
    system = None
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.MOONS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        candidate = StarSystem(system_config=cfg)
        real_planets = [p for p in candidate.planets if p.body_type != "a"]
        if real_planets and any(p.moons for p in real_planets):
            system = candidate
            break
    if system is None:
        pytest.fail("could not generate a single-star system with planets and at least one moon")

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, system.system_config)

        star_id = conn.execute(
            "SELECT id FROM stars WHERE star_system_id = ? AND role = 'single'", (system_id,)
        ).fetchone()["id"]
        before_offset = conn.execute(
            "SELECT reflex_offset_x_km FROM stars WHERE id = ?", (star_id,)
        ).fetchone()["reflex_offset_x_km"]

        # Big enough elapsed time to clearly move every planet, regardless
        # of its randomly generated period.
        max_period = conn.execute(
            "SELECT MAX(period_years) AS p FROM planets WHERE star_id = ?", (star_id,)
        ).fetchone()["p"]
        counts = _db.advance_orbital_phases(conn, elapsed_years=max_period * 137.25)
        assert counts["star_reflex_offsets"] > 0
        assert counts["planet_reflex_offsets"] > 0

        star_row = conn.execute(
            "SELECT mass_kg, reflex_offset_x_km, reflex_offset_y_km, reflex_offset_z_km "
            "FROM stars WHERE id = ?", (star_id,)
        ).fetchone()
        planet_rows = conn.execute(
            "SELECT mass_kg, position_x_km, position_y_km, position_z_km FROM planets WHERE star_id = ?",
            (star_id,),
        ).fetchall()
    finally:
        conn.close()

    # The offset actually changed -- not accidentally still the
    # generation-time value.
    assert star_row["reflex_offset_x_km"] != pytest.approx(before_offset)

    expected_x = -sum(
        (p["mass_kg"] / (star_row["mass_kg"] + p["mass_kg"])) * p["position_x_km"] for p in planet_rows
    )
    assert star_row["reflex_offset_x_km"] == pytest.approx(expected_x, abs=1e-3)


def test_advance_orbital_phases_rejects_negative_elapsed_years(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError):
            _db.advance_orbital_phases(conn, elapsed_years=-1.0)
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Comet orbital motion updates (stellarObjects._db.advance_comet_orbits) --
# a separate Python-loop call from advance_orbital_phases above, since a
# comet's position isn't a linear function of elapsed time the way a
# circular planet/moon orbit's is -- see that function's own docstring.
# ---------------------------------------------------------------------------

def test_advance_comet_orbits_advances_elliptical_mean_anomaly_and_recomputes_position(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.COMETS = False
    system = StarSystem(system_config=cfg)
    comet = Comet(cfg, primary_mass_solar=system.star.mass / pc.SOLAR_MASS_TO_KG, orbit_type="elliptical")
    comet.mean_anomaly_deg = 10.0
    comet.update_orbital_state()
    system.comets = [comet]

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)

        before = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
        step_years = before["orbital_period_years"] * 0.05  # 18 degrees of mean anomaly

        updated = _db.advance_comet_orbits(conn, step_years)
        assert updated == 1

        after = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
    finally:
        conn.close()

    assert math.isclose(after["mean_anomaly_deg"], 28.0, rel_tol=1e-9)  # 10 + 18
    assert after["parabolic_mean_anomaly"] is None
    # Moved further from perihelion (mean anomaly 10 -> 28, still well
    # short of aphelion) -- distance must have increased.
    assert after["distance_km"] > before["distance_km"]
    # Position is a pure function of distance/orientation/anomaly, so it
    # must have moved along with mean_anomaly_deg/distance_km.
    assert (after["position_x_km"], after["position_y_km"], after["position_z_km"]) != (
        before["position_x_km"], before["position_y_km"], before["position_z_km"]
    )


def test_advance_comet_orbits_advances_parabolic_mean_anomaly_without_wrapping(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.COMETS = False
    system = StarSystem(system_config=cfg)
    comet = Comet(cfg, primary_mass_solar=system.star.mass / pc.SOLAR_MASS_TO_KG, orbit_type="parabolic")
    comet.parabolic_mean_anomaly = -1.0  # still approaching perihelion
    comet.update_orbital_state()
    system.comets = [comet]

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)

        before = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
        assert before["mean_anomaly_deg"] is None
        assert before["min_update_interval_years"] is None

        # A large elapsed time -- since a parabolic anomaly doesn't wrap
        # (unlike mean_anomaly_deg's MOD 360), this should NOT be clamped
        # or wrapped, just added linearly.
        updated = _db.advance_comet_orbits(conn, elapsed_years=50.0)
        assert updated == 1

        after = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
    finally:
        conn.close()

    assert after["mean_anomaly_deg"] is None
    assert after["parabolic_mean_anomaly"] > before["parabolic_mean_anomaly"]
    assert after["distance_km"] > before["distance_km"]


def test_advance_comet_orbits_skips_elliptical_rows_below_min_update_interval(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.COMETS = False
    system = StarSystem(system_config=cfg)
    comet = Comet(cfg, primary_mass_solar=system.star.mass / pc.SOLAR_MASS_TO_KG, orbit_type="elliptical")
    system.comets = [comet]

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)

        before = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
        # Comfortably below this comet's own min_update_interval_years floor.
        tiny_elapsed = before["min_update_interval_years"] / 2
        updated = _db.advance_comet_orbits(conn, tiny_elapsed)
        assert updated == 0

        after = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
    finally:
        conn.close()

    assert after["mean_anomaly_deg"] == before["mean_anomaly_deg"]
    assert after["distance_km"] == before["distance_km"]


def test_advance_comet_orbits_rejects_negative_elapsed_years(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError):
            _db.advance_comet_orbits(conn, elapsed_years=-1.0)
    finally:
        conn.close()


def test_migrate_v8_to_v9_adds_orbital_motion_columns(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing v8
    database by tearing the v9 through v15 additions back out (a
    real v8 database would have none of them), then confirms
    `migrate_database` puts everything back and reports the database as
    fully current again (not just v9 -- `migrate_database` applies every
    step up to `SCHEMA_VERSION` in one call, so a v8 database now lands on
    v15 directly).
    """
    conn = _db.get_connection(mysql_config)
    try:
        for table in ("planets", "moons"):
            conn.execute(
                f"ALTER TABLE {table} "
                f"DROP COLUMN orbital_inclination_deg, DROP COLUMN orbital_ascending_node_deg, "
                f"DROP COLUMN orbital_phase_deg, "
                f"DROP COLUMN position_x_km, DROP COLUMN position_y_km, DROP COLUMN position_z_km, "
                f"DROP COLUMN orbital_speed_kms, DROP COLUMN min_update_interval_years, "
                f"DROP COLUMN rotation_period_hours"
            )
        conn.execute("DROP TABLE orbit_simulation_state")
        conn.execute(
            "ALTER TABLE stars "
            "DROP COLUMN galactic_orbital_speed_kms, DROP COLUMN galactic_orbital_period_gy, "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years, "
            "DROP COLUMN wide_binary_a_crit_km"
        )
        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_galactic_orbital_speed_kms, DROP COLUMN binary_galactic_orbital_period_gy, "
            "DROP COLUMN binary_galactic_orbital_phase_deg, DROP COLUMN binary_galactic_min_update_interval_years, "
            "DROP COLUMN binary_mutual_orbital_period_years, DROP COLUMN binary_mutual_orbital_speed_kms, "
            "DROP COLUMN binary_mutual_orbital_inclination_deg, DROP COLUMN binary_mutual_orbital_ascending_node_deg, "
            "DROP COLUMN binary_mutual_orbital_phase_deg, DROP COLUMN binary_mutual_min_update_interval_years, "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        _drop_v20_trajectory_columns(conn)
        _drop_v21_phenomenon_columns(conn)
        _drop_v22_search_indexes(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (8)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 8
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        planet_columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM planets").fetchall()}
        assert {
            "orbital_inclination_deg", "orbital_ascending_node_deg",
            "orbital_phase_deg", "rotation_period_hours",
            "position_x_km", "position_y_km", "position_z_km", "orbital_speed_kms",
            "min_update_interval_years",
        } <= planet_columns
        assert _db.get_orbit_update_elapsed_years(conn) is None  # table exists, no row yet

        star_columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM stars").fetchall()}
        assert {
            "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
            "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
        } <= star_columns
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION


def test_migrate_v9_to_v10_adds_galactic_orbit_columns(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing v9
    database by tearing the v10 through v15 additions back out (the new
    `stars`/`star_systems`/`planets`/`moons`/`asteroid_belts` columns and the
    `schema_migrations` v10/v11/v12/v13/v14/v15 rows) before calling
    `migrate_database`, and confirms it puts everything back and reports
    the database as current again, without disturbing the v9
    orbital-motion columns already in place.
    """
    conn = _db.get_connection(mysql_config)
    try:
        conn.execute(
            "ALTER TABLE stars "
            "DROP COLUMN galactic_orbital_speed_kms, DROP COLUMN galactic_orbital_period_gy, "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years, "
            "DROP COLUMN wide_binary_a_crit_km"
        )
        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_galactic_orbital_speed_kms, DROP COLUMN binary_galactic_orbital_period_gy, "
            "DROP COLUMN binary_galactic_orbital_phase_deg, DROP COLUMN binary_galactic_min_update_interval_years, "
            "DROP COLUMN binary_mutual_orbital_period_years, DROP COLUMN binary_mutual_orbital_speed_kms, "
            "DROP COLUMN binary_mutual_orbital_inclination_deg, DROP COLUMN binary_mutual_orbital_ascending_node_deg, "
            "DROP COLUMN binary_mutual_orbital_phase_deg, DROP COLUMN binary_mutual_min_update_interval_years, "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        for table in ("planets", "moons"):
            conn.execute(
                f"ALTER TABLE {table} "
                f"DROP COLUMN position_x_km, DROP COLUMN position_y_km, DROP COLUMN position_z_km, "
                f"DROP COLUMN orbital_speed_kms, DROP COLUMN min_update_interval_years"
            )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        _drop_v20_trajectory_columns(conn)
        _drop_v21_phenomenon_columns(conn)
        _drop_v22_search_indexes(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (9)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 9
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        star_columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM stars").fetchall()}
        assert {
            "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
            "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
        } <= star_columns

        star_system_columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM star_systems").fetchall()}
        assert {
            "binary_galactic_orbital_speed_kms", "binary_galactic_orbital_period_gy",
            "binary_galactic_orbital_phase_deg", "binary_galactic_min_update_interval_years",
            "binary_mutual_orbital_period_years", "binary_mutual_orbital_speed_kms",
            "binary_mutual_orbital_inclination_deg", "binary_mutual_orbital_ascending_node_deg",
            "binary_mutual_orbital_phase_deg", "binary_mutual_min_update_interval_years",
            "binary_mutual_position_x_km", "binary_mutual_position_y_km", "binary_mutual_position_z_km",
        } <= star_system_columns

        planet_columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM planets").fetchall()}
        assert {
            "position_x_km", "position_y_km", "position_z_km", "orbital_speed_kms",
            "min_update_interval_years",
        } <= planet_columns
        # v9's orbital-motion columns are untouched by these migration steps.
        assert "orbital_phase_deg" in planet_columns
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION


def test_migrate_v10_to_v11_adds_and_backfills_position_columns(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing
    v10 database by tearing the v11 through v15 additions back out (the new
    `planets`/`moons`/`stars`/`star_systems`/`asteroid_belts` columns and the
    `schema_migrations` v11/v12/v13/v14/v15 rows) before calling
    `migrate_database`, and confirms it not only adds the columns back but
    backfills them with real derived values (not an arbitrary placeholder)
    from each row's own already-stored `distance_km`/`period_years`/
    orbital-motion columns -- see `_migrate_v10_to_v11`'s docstring.
    """
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)

        for table in ("planets", "moons"):
            conn.execute(
                f"ALTER TABLE {table} "
                f"DROP COLUMN position_x_km, DROP COLUMN position_y_km, DROP COLUMN position_z_km, "
                f"DROP COLUMN orbital_speed_kms, DROP COLUMN min_update_interval_years"
            )
        # A "v10" database also predates v13/v14/v15's stars/star_systems/asteroid_belts columns.
        conn.execute(
            "ALTER TABLE stars "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years, "
            "DROP COLUMN wide_binary_a_crit_km"
        )
        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_galactic_orbital_phase_deg, DROP COLUMN binary_galactic_min_update_interval_years, "
            "DROP COLUMN binary_mutual_orbital_period_years, DROP COLUMN binary_mutual_orbital_speed_kms, "
            "DROP COLUMN binary_mutual_orbital_inclination_deg, DROP COLUMN binary_mutual_orbital_ascending_node_deg, "
            "DROP COLUMN binary_mutual_orbital_phase_deg, DROP COLUMN binary_mutual_min_update_interval_years, "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        _drop_v20_trajectory_columns(conn)
        _drop_v21_phenomenon_columns(conn)
        _drop_v22_search_indexes(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (10)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 10
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = conn.execute(
            "SELECT distance_km, orbital_inclination_deg, orbital_ascending_node_deg, orbital_phase_deg, "
            "period_years, position_x_km, position_y_km, position_z_km, orbital_speed_kms "
            "FROM planets WHERE star_system_id = ?", (system_id,)
        ).fetchall()
        assert rows, "test fixture must actually contain planets"
        for row in rows:
            expected_x_au, expected_y_au, expected_z_au = orbital_position_au(
                row["distance_km"] / pc.AU_TO_KM, row["orbital_inclination_deg"],
                row["orbital_ascending_node_deg"], row["orbital_phase_deg"],
            )
            assert row["position_x_km"] == pytest.approx(expected_x_au * pc.AU_TO_KM, abs=1e-3)
            assert row["position_y_km"] == pytest.approx(expected_y_au * pc.AU_TO_KM, abs=1e-3)
            assert row["position_z_km"] == pytest.approx(expected_z_au * pc.AU_TO_KM, abs=1e-3)

            expected_speed = (2 * math.pi * row["distance_km"]) / (row["period_years"] * pc.SECONDS_PER_YEAR)
            assert row["orbital_speed_kms"] == pytest.approx(expected_speed, rel=1e-9)
    finally:
        conn.close()


def test_migrate_v11_to_v12_adds_and_backfills_min_update_interval_years(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing
    v11 database by tearing the v12 through v15 additions back out (the
    new `planets`/`moons.min_update_interval_years` column, `stars`/
    `star_systems`/`asteroid_belts`' v13/v14/v15 columns, and the
    `schema_migrations` v12/v13/v14/v15 rows) before calling
    `migrate_database`, and confirms it not only adds the column back but
    backfills it with a real derived value (not an arbitrary placeholder)
    from each row's own already-stored `period_years` -- see
    `_migrate_v11_to_v12`'s docstring. Scoped to `planets`/`moons` only for
    v12: `stars`/`star_systems` didn't gain a column there (see
    `_migrate_v11_to_v12`'s docstring for why), only later at v13.
    """
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)

        for table in ("planets", "moons"):
            conn.execute(f"ALTER TABLE {table} DROP COLUMN min_update_interval_years")
        conn.execute(
            "ALTER TABLE stars "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years, "
            "DROP COLUMN wide_binary_a_crit_km"
        )
        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_galactic_orbital_phase_deg, DROP COLUMN binary_galactic_min_update_interval_years, "
            "DROP COLUMN binary_mutual_orbital_period_years, DROP COLUMN binary_mutual_orbital_speed_kms, "
            "DROP COLUMN binary_mutual_orbital_inclination_deg, DROP COLUMN binary_mutual_orbital_ascending_node_deg, "
            "DROP COLUMN binary_mutual_orbital_phase_deg, DROP COLUMN binary_mutual_min_update_interval_years, "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        _drop_v20_trajectory_columns(conn)
        _drop_v21_phenomenon_columns(conn)
        _drop_v22_search_indexes(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (11)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 11
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        planet_rows = conn.execute(
            "SELECT period_years, min_update_interval_years "
            "FROM planets WHERE star_system_id = ?", (system_id,)
        ).fetchall()
        assert planet_rows, "test fixture must actually contain planets"
        for row in planet_rows:
            expected = minimum_update_interval_years(row["period_years"])
            assert row["min_update_interval_years"] == pytest.approx(expected, rel=1e-9)

        moon_rows = conn.execute(
            "SELECT period_years, min_update_interval_years "
            "FROM moons WHERE star_system_id = ?", (system_id,)
        ).fetchall()
        assert moon_rows, "test fixture must actually contain moons"
        for row in moon_rows:
            expected = minimum_update_interval_years(row["period_years"])
            assert row["min_update_interval_years"] == pytest.approx(expected, rel=1e-9)
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION


def test_migrate_v12_to_v13_adds_and_backfills_star_motion_columns(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing
    v12 database by tearing the v13 through v15 additions back out (the new
    `stars`/`star_systems`/`asteroid_belts` columns and the
    `schema_migrations` v13/v14/v15 rows) before calling `migrate_database`,
    and confirms it not only adds the
    columns back but backfills every real-derived one (the two
    `*_min_update_interval_years` guards, plus the binary pair's
    `binary_mutual_orbital_period_years`/`_speed_kms`) with real values
    from each row's own already-stored `galactic_orbital_period_gy`/
    `binary_separation_km`/`binary_effective_mass_kg` -- see
    `_migrate_v12_to_v13`'s docstring. The phase/orientation columns with
    no derivable "correct" value (`galactic_orbital_phase_deg` and the
    three `binary_mutual_orbital_{inclination,ascending_node,phase}_deg`
    columns) get the arbitrary `0` placeholder instead.
    """
    system, cfg = _make_system_with_moons_and_belt()
    binary_config = SystemConfig()
    binary_config.BINARY_SYSTEM = True
    binary_config.WIDE_BINARY = False  # pin to close/P-type -- these assertions are BinaryStarProxy-specific
    binary_system = StarSystem(binary_config)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
            binary_id = _db.insert_star_system(conn, binary_system, binary_config)

        conn.execute(
            "ALTER TABLE stars "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years, "
            "DROP COLUMN wide_binary_a_crit_km"
        )
        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_galactic_orbital_phase_deg, DROP COLUMN binary_galactic_min_update_interval_years, "
            "DROP COLUMN binary_mutual_orbital_period_years, DROP COLUMN binary_mutual_orbital_speed_kms, "
            "DROP COLUMN binary_mutual_orbital_inclination_deg, DROP COLUMN binary_mutual_orbital_ascending_node_deg, "
            "DROP COLUMN binary_mutual_orbital_phase_deg, DROP COLUMN binary_mutual_min_update_interval_years, "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        _drop_v20_trajectory_columns(conn)
        _drop_v21_phenomenon_columns(conn)
        _drop_v22_search_indexes(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (12)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 12
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        star_rows = conn.execute(
            "SELECT galactic_orbital_phase_deg, galactic_orbital_period_gy, galactic_min_update_interval_years "
            "FROM stars WHERE star_system_id = ?", (system_id,)
        ).fetchall()
        assert star_rows, "test fixture must actually contain a star"
        for row in star_rows:
            assert row["galactic_orbital_phase_deg"] == 0
            expected = minimum_update_interval_years(row["galactic_orbital_period_gy"] * 1e9)
            assert row["galactic_min_update_interval_years"] == pytest.approx(expected, rel=1e-9)

        binary_row = conn.execute(
            "SELECT binary_galactic_orbital_phase_deg, binary_galactic_orbital_period_gy, "
            "binary_galactic_min_update_interval_years, binary_separation_km, binary_effective_mass_kg, "
            "binary_mutual_orbital_period_years, binary_mutual_orbital_speed_kms, "
            "binary_mutual_orbital_inclination_deg, binary_mutual_orbital_ascending_node_deg, "
            "binary_mutual_orbital_phase_deg, binary_mutual_min_update_interval_years "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()
        assert binary_row["binary_galactic_orbital_phase_deg"] == 0
        expected_galactic_guard = minimum_update_interval_years(
            binary_row["binary_galactic_orbital_period_gy"] * 1e9
        )
        assert binary_row["binary_galactic_min_update_interval_years"] == pytest.approx(
            expected_galactic_guard, rel=1e-9
        )

        expected_mutual_period = calculate_orbital_period_years(
            binary_row["binary_separation_km"] / pc.AU_TO_KM, binary_row["binary_effective_mass_kg"]
        )
        assert binary_row["binary_mutual_orbital_period_years"] == pytest.approx(expected_mutual_period, rel=1e-9)

        expected_mutual_speed = circular_orbital_speed_kms(
            binary_row["binary_separation_km"] / pc.AU_TO_KM, expected_mutual_period
        )
        assert binary_row["binary_mutual_orbital_speed_kms"] == pytest.approx(expected_mutual_speed, rel=1e-9)

        assert binary_row["binary_mutual_orbital_inclination_deg"] == 0
        assert binary_row["binary_mutual_orbital_ascending_node_deg"] == 0
        assert binary_row["binary_mutual_orbital_phase_deg"] == 0

        expected_mutual_guard = minimum_update_interval_years(expected_mutual_period)
        assert binary_row["binary_mutual_min_update_interval_years"] == pytest.approx(
            expected_mutual_guard, rel=1e-9
        )
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION


def test_migrate_v13_to_v14_adds_and_backfills_binary_mutual_position(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing
    v13 database by tearing the v14 *and* v15 additions back out (the new
    `star_systems.binary_mutual_position_x_km`/`_y_km`/`_z_km` columns, the
    v15 wide-binary columns, and the `schema_migrations` v14/v15 rows)
    before calling `migrate_database`, and confirms it not only adds the
    columns back but backfills them with a real derived value (not an
    arbitrary placeholder) from each row's own already-stored
    `binary_separation_km` and v13
    `binary_mutual_orbital_{inclination,ascending_node,phase}_deg` columns
    -- see `_migrate_v13_to_v14`'s docstring.
    """
    binary_config = SystemConfig()
    binary_config.BINARY_SYSTEM = True
    binary_config.WIDE_BINARY = False  # pin to close/P-type -- these assertions are BinaryStarProxy-specific
    binary_system = StarSystem(binary_config)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            binary_id = _db.insert_star_system(conn, binary_system, binary_config)

        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute("ALTER TABLE stars DROP COLUMN wide_binary_a_crit_km")
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        _drop_v20_trajectory_columns(conn)
        _drop_v21_phenomenon_columns(conn)
        _drop_v22_search_indexes(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (13)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 13
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute(
            "SELECT binary_separation_km, binary_mutual_orbital_inclination_deg, "
            "binary_mutual_orbital_ascending_node_deg, binary_mutual_orbital_phase_deg, "
            "binary_mutual_position_x_km, binary_mutual_position_y_km, binary_mutual_position_z_km "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()
        expected_x_au, expected_y_au, expected_z_au = orbital_position_au(
            row["binary_separation_km"] / pc.AU_TO_KM, row["binary_mutual_orbital_inclination_deg"],
            row["binary_mutual_orbital_ascending_node_deg"], row["binary_mutual_orbital_phase_deg"],
        )
        assert row["binary_mutual_position_x_km"] == pytest.approx(expected_x_au * pc.AU_TO_KM, abs=1e-3)
        assert row["binary_mutual_position_y_km"] == pytest.approx(expected_y_au * pc.AU_TO_KM, abs=1e-3)
        assert row["binary_mutual_position_z_km"] == pytest.approx(expected_z_au * pc.AU_TO_KM, abs=1e-3)
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION


def test_migrate_v19_to_v20_backfills_star_and_planet_reflex_offsets(mysql_config):
    """
    Simulates an existing v19 database (tearing v20's `stars`/`planets`/
    `star_systems` trajectory columns back out) for a single-star system
    with both star-hosted planets and a moon-having planet, and confirms
    `migrate_database` backfills `stars.reflex_offset_*_km`/
    `planets.reflex_offset_*_km` with real derived values -- the exact
    same `utils.calculate_reflex_offset` formula generation time uses --
    rather than leaving them `NULL`. See `_migrate_v19_to_v20`'s docstring.
    """
    system = None
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.MOONS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        candidate = StarSystem(system_config=cfg)
        real_planets = [p for p in candidate.planets if p.body_type != "a"]
        if real_planets and any(p.moons for p in real_planets):
            system = candidate
            break
    if system is None:
        pytest.fail("could not generate a single-star system with planets and at least one moon")

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, system.system_config)
        _drop_v20_trajectory_columns(conn)
        _drop_v21_phenomenon_columns(conn)
        _drop_v22_search_indexes(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (19)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 19
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        star_row = conn.execute(
            "SELECT id, mass_kg, reflex_offset_x_km, reflex_offset_y_km, reflex_offset_z_km "
            "FROM stars WHERE star_system_id = ? AND role = 'single'", (system_id,)
        ).fetchone()
        planet_rows = conn.execute(
            "SELECT mass_kg, position_x_km, position_y_km, position_z_km "
            "FROM planets WHERE star_id = ?", (star_row["id"],)
        ).fetchall()
        expected_x, expected_y, expected_z = -sum(
            (p["mass_kg"] / (star_row["mass_kg"] + p["mass_kg"])) * p["position_x_km"] for p in planet_rows
        ), -sum(
            (p["mass_kg"] / (star_row["mass_kg"] + p["mass_kg"])) * p["position_y_km"] for p in planet_rows
        ), -sum(
            (p["mass_kg"] / (star_row["mass_kg"] + p["mass_kg"])) * p["position_z_km"] for p in planet_rows
        )
        assert star_row["reflex_offset_x_km"] == pytest.approx(expected_x, abs=1e-3)
        assert star_row["reflex_offset_y_km"] == pytest.approx(expected_y, abs=1e-3)
        assert star_row["reflex_offset_z_km"] == pytest.approx(expected_z, abs=1e-3)

        moon_planet_row = conn.execute(
            "SELECT p.id, p.mass_kg, p.reflex_offset_x_km, p.reflex_offset_y_km, p.reflex_offset_z_km "
            "FROM planets p WHERE p.star_id = ? "
            "AND EXISTS (SELECT 1 FROM moons m WHERE m.planet_id = p.id)",
            (star_row["id"],),
        ).fetchone()
        moon_rows = conn.execute(
            "SELECT mass_kg, position_x_km, position_y_km, position_z_km FROM moons WHERE planet_id = ?",
            (moon_planet_row["id"],),
        ).fetchall()
        expected_px = -sum(
            (m["mass_kg"] / (moon_planet_row["mass_kg"] + m["mass_kg"])) * m["position_x_km"] for m in moon_rows
        )
        expected_py = -sum(
            (m["mass_kg"] / (moon_planet_row["mass_kg"] + m["mass_kg"])) * m["position_y_km"] for m in moon_rows
        )
        expected_pz = -sum(
            (m["mass_kg"] / (moon_planet_row["mass_kg"] + m["mass_kg"])) * m["position_z_km"] for m in moon_rows
        )
        assert moon_planet_row["reflex_offset_x_km"] == pytest.approx(expected_px, abs=1e-3)
        assert moon_planet_row["reflex_offset_y_km"] == pytest.approx(expected_py, abs=1e-3)
        assert moon_planet_row["reflex_offset_z_km"] == pytest.approx(expected_pz, abs=1e-3)
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION


def test_migrate_v19_to_v20_backfills_binary_trajectory_columns(mysql_config):
    """
    Same simulated-old-database shape as the reflex-offset migration test
    above, for a 'close' (P-type) binary with circumbinary planets instead
    -- confirms `binary_secondary_mass_fraction`, `binary_primary/
    secondary_position_*_km`, and `binary_planetary_wobble_*_km` are all
    backfilled with real derived values from data the row already had
    (both stars' own `mass_kg`, the already-stored
    `binary_mutual_position_*_km`, and the circumbinary planets'
    `mass_kg`/`position_*_km`), not left `NULL`.
    """
    system = None
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = False
        cfg.MAX_PLANETS = True
        candidate = StarSystem(system_config=cfg)
        if [p for p in candidate.planets if p.body_type != "a"]:
            system = candidate
            break
    if system is None:
        pytest.fail("could not generate a close-binary system with circumbinary planets")

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, system.system_config)
        _drop_v20_trajectory_columns(conn)
        _drop_v21_phenomenon_columns(conn)
        _drop_v22_search_indexes(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (19)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 19
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        ss_row = conn.execute(
            "SELECT binary_effective_mass_kg, binary_mutual_position_x_km, binary_mutual_position_y_km, "
            "binary_mutual_position_z_km, binary_secondary_mass_fraction, "
            "binary_primary_position_x_km, binary_primary_position_y_km, binary_primary_position_z_km, "
            "binary_secondary_position_x_km, binary_secondary_position_y_km, binary_secondary_position_z_km, "
            "binary_planetary_wobble_x_km, binary_planetary_wobble_y_km, binary_planetary_wobble_z_km "
            "FROM star_systems WHERE id = ?", (system_id,)
        ).fetchone()
        primary_mass = conn.execute(
            "SELECT mass_kg FROM stars WHERE star_system_id = ? AND role = 'primary'", (system_id,)
        ).fetchone()["mass_kg"]
        secondary_mass = conn.execute(
            "SELECT mass_kg FROM stars WHERE star_system_id = ? AND role = 'secondary'", (system_id,)
        ).fetchone()["mass_kg"]

        expected_mu = secondary_mass / (primary_mass + secondary_mass)
        assert ss_row["binary_secondary_mass_fraction"] == pytest.approx(expected_mu, rel=1e-9)

        assert ss_row["binary_primary_position_x_km"] == pytest.approx(
            -expected_mu * ss_row["binary_mutual_position_x_km"], abs=1e-3
        )
        assert ss_row["binary_secondary_position_x_km"] == pytest.approx(
            (1 - expected_mu) * ss_row["binary_mutual_position_x_km"], abs=1e-3
        )
        # secondary_position == primary_position + mutual_position always.
        assert ss_row["binary_secondary_position_x_km"] == pytest.approx(
            ss_row["binary_primary_position_x_km"] + ss_row["binary_mutual_position_x_km"], abs=1e-3
        )

        circumbinary_rows = conn.execute(
            "SELECT mass_kg, position_x_km, position_y_km, position_z_km "
            "FROM planets WHERE star_system_id = ? AND star_id IS NULL", (system_id,)
        ).fetchall()
        assert circumbinary_rows  # MAX_PLANETS=True on a close binary should guarantee at least one
        total_mass = ss_row["binary_effective_mass_kg"]
        expected_wobble_x = -sum(
            (p["mass_kg"] / (total_mass + p["mass_kg"])) * p["position_x_km"] for p in circumbinary_rows
        )
        assert ss_row["binary_planetary_wobble_x_km"] == pytest.approx(expected_wobble_x, abs=1e-3)
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION


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


def test_migrate_v20_to_v21_adds_sector_placement_columns(mysql_config):
    """
    Simulates a database created under schema v20 (before black_holes/
    neutron_stars gained sector_id/center_x/y/z_pc/galactic_radius_pc),
    with a pre-existing standalone black hole and neutron star already in
    it, then checks that migrate_database brings it up to v21: the new
    columns exist, are NULL on the pre-existing rows (no backfill is
    possible -- a pre-v21 row was always fully standalone, with no
    placement to recover), and are fully usable for a phenomenon inserted
    after the migration runs.
    """
    cfg = SystemConfig()
    bh = BlackHole(cfg)
    ns = NeutronStar(cfg)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            pre_existing_bh_id = _db.insert_black_hole(conn, bh)
            pre_existing_ns_id = _db.insert_neutron_star(conn, ns)

        _drop_v21_phenomenon_columns(conn)
        _drop_v22_search_indexes(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (20)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 20

        bh_columns_before = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM black_holes").fetchall()}
        assert "sector_id" not in bh_columns_before
        assert "center_x_pc" not in bh_columns_before
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        bh_columns_after = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM black_holes").fetchall()}
        ns_columns_after = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM neutron_stars").fetchall()}
        for columns in (bh_columns_after, ns_columns_after):
            assert {"sector_id", "center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc"} <= columns

        # No backfill possible -- the pre-existing rows simply gain NULL columns.
        pre_bh_row = conn.execute(
            "SELECT sector_id, center_x_pc FROM black_holes WHERE id = ?", (pre_existing_bh_id,)
        ).fetchone()
        pre_ns_row = conn.execute(
            "SELECT sector_id, center_x_pc FROM neutron_stars WHERE id = ?", (pre_existing_ns_id,)
        ).fetchone()
        assert pre_bh_row["sector_id"] is None and pre_bh_row["center_x_pc"] is None
        assert pre_ns_row["sector_id"] is None and pre_ns_row["center_x_pc"] is None

        # The migrated schema is fully usable going forward: a sector-linked,
        # placed black hole round-trips through the new columns correctly.
        sector = SpaceSector("Migration Test Sector", edge_ly=11.5)
        galaxy_position = {
            "center_x_pc": 5.0, "center_y_pc": -3.0, "center_z_pc": 1.0,
            "galactic_radius_pc": math.sqrt(5.0 ** 2 + 3.0 ** 2 + 1.0 ** 2),
        }
        sector_id = _db.save_sector(sector, config=mysql_config, galaxy_position=galaxy_position)

        new_bh = BlackHole(SystemConfig())
        placement = {
            "center_x_pc": 5.5, "center_y_pc": -3.0, "center_z_pc": 1.0,
            "galactic_radius_pc": math.sqrt(5.5 ** 2 + 3.0 ** 2 + 1.0 ** 2),
        }
        with conn:
            new_bh_id = _db.insert_black_hole(conn, new_bh, sector_id=sector_id, placement=placement)
        new_bh_row = conn.execute("SELECT * FROM black_holes WHERE id = ?", (new_bh_id,)).fetchone()
        assert new_bh_row["sector_id"] == sector_id
        assert new_bh_row["center_x_pc"] == pytest.approx(5.5)
        assert new_bh_row["galactic_radius_pc"] == pytest.approx(placement["galactic_radius_pc"])
    finally:
        conn.close()


def _drop_v27_timestamp_columns(conn):
    """
    Drops v27's `created_at`/`modified_at` columns and `modified_at`
    indexes from every `_db.TIMESTAMPED_TABLES` table (keeping
    `star_systems.created_at`, which predates v27) -- the same "already
    exists on a freshly-bootstrapped test database" reasoning
    `_drop_v17_phenomenon_columns` gives for its own tables. See
    `schema.sql`'s "v27" header note.
    """
    for table in _db.TIMESTAMPED_TABLES:
        drops = [f"DROP INDEX idx_{table}_modified_at", "DROP COLUMN modified_at"]
        if table != "star_systems":
            drops.append("DROP COLUMN created_at")
        conn.execute(f"ALTER TABLE {table} " + ", ".join(drops))


def test_migrate_v26_to_v27_adds_and_backfills_row_timestamps(mysql_config, monkeypatch):
    """
    Simulates a database created under schema v26 (no row timestamps
    beyond `star_systems.created_at`), with a sector holding two systems,
    an empty sector and a standalone black hole already in it, then checks
    that migrate_database brings it up to v27: every top-level table has
    both columns and the `modified_at` index; each system's `modified_at`
    is backfilled from its own `created_at`; the populated sector takes
    its oldest system's `created_at` for both; and rows with nothing to
    recover from (the empty sector, the black hole) get the migration's
    own time. A batch size of 1 makes the backfill cross batch
    boundaries.
    """
    monkeypatch.setattr(_db, "_V27_BACKFILL_BATCH_SIZE", 1)
    first_created = datetime.datetime(2020, 1, 1, 0, 0, 0)
    second_created = datetime.datetime(2021, 6, 15, 12, 30, 0)

    conn = _db.get_connection(mysql_config)
    try:
        sector_id = _db.save_sector(SpaceSector("Timestamp Migration Sector", edge_ly=11.5), config=mysql_config)
        empty_sector_id = _db.save_sector(SpaceSector("Empty Timestamp Sector", edge_ly=11.5), config=mysql_config)
        with conn:
            system_ids = []
            for _ in range(2):
                system, cfg = _make_system_with_moons_and_belt()
                system_ids.append(_db.insert_star_system(conn, system, cfg))
            bh_id = _db.insert_black_hole(conn, BlackHole(SystemConfig()))

        _drop_v27_timestamp_columns(conn)
        for system_id, created in zip(system_ids, (second_created, first_created)):
            conn.execute(
                "UPDATE star_systems SET sector_id = ?, created_at = ? WHERE id = ?",
                (sector_id, created, system_id),
            )
        conn.execute("DELETE FROM schema_migrations WHERE version IN (27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (26)")
        conn.commit()

        sector_columns_before = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM sectors").fetchall()}
        assert "modified_at" not in sector_columns_before
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in _db.TIMESTAMPED_TABLES:
            columns = {row["Field"] for row in conn.execute(f"SHOW COLUMNS FROM {table}").fetchall()}
            assert {"created_at", "modified_at"} <= columns, table
            indexes = {row["Key_name"] for row in conn.execute(f"SHOW INDEX FROM {table}").fetchall()}
            assert f"idx_{table}_modified_at" in indexes, table

        def timestamps(table, row_id):
            return conn.execute(f"SELECT created_at, modified_at FROM {table} WHERE id = ?", (row_id,)).fetchone()

        for system_id, created in zip(system_ids, (second_created, first_created)):
            row = timestamps("star_systems", system_id)
            assert row["created_at"] == created
            assert row["modified_at"] == created

        row = timestamps("sectors", sector_id)
        assert row["created_at"] == first_created
        assert row["modified_at"] == first_created

        for table, row_id in (("sectors", empty_sector_id), ("black_holes", bh_id)):
            row = timestamps(table, row_id)
            assert row["created_at"] > second_created, table
            assert row["modified_at"] > second_created, table
    finally:
        conn.close()

    # Running it again on an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION



def _drop_v28_placement_columns(conn):
    """
    Drops v28's placement columns, indexes and CHECKs from every
    `_db.V28_PLACED_TABLES` table -- see `schema.sql`'s "v28" header note
    and `_drop_v27_timestamp_columns`'s reasoning.
    """
    for table in _db.V28_PLACED_TABLES:
        conn.execute(f"ALTER TABLE {table} DROP CONSTRAINT chk_{table}_placement")
        conn.execute(
            f"ALTER TABLE {table} DROP INDEX idx_{table}_center, DROP INDEX idx_{table}_galactic_radius_pc, "
            "DROP COLUMN center_x_pc, DROP COLUMN center_y_pc, DROP COLUMN center_z_pc, "
            "DROP COLUMN galactic_radius_pc"
        )


def test_migrate_v27_to_v28_adds_and_backfills_phenomenon_placement(mysql_config, monkeypatch):
    """
    Simulates a v27 database holding a supernova remnant (with an embedded
    black hole), a rogue planet and an interstellar comet linked to a
    galaxy-placed sector, plus a rogue planet in a never-placed sector and
    a comet with no sector at all. After migrating: the three tables have
    v28's columns, indexes and CHECK; each row in the placed sector sits
    inside that sector's own cube; the remnant's black hole shares
    its remnant's sector and center; and the rows with no placed sector
    stay unplaced. A batch size of 1 makes the backfill cross batch
    boundaries.
    """
    from stellarObjects.roguePlanetData import InterstellarComet, RoguePlanet
    from stellarObjects.supernovaRemnantData import SupernovaRemnant
    from stellarObjects.utils import ly_to_pc

    monkeypatch.setattr(_db, "_V27_BACKFILL_BATCH_SIZE", 1)
    center_pc = (4000.0, 3000.0, 20.0)
    galaxy_position = {
        "center_x_pc": center_pc[0], "center_y_pc": center_pc[1], "center_z_pc": center_pc[2],
        "galactic_radius_pc": math.sqrt(sum(c * c for c in center_pc)),
    }
    placed_sector_id = _db.save_sector(
        SpaceSector("Placement Migration Sector", edge_ly=11.5), config=mysql_config,
        galaxy_position=galaxy_position,
    )
    unplaced_sector_id = _db.save_sector(SpaceSector("Unplaced Migration Sector", edge_ly=11.5), config=mysql_config)

    remnant = SupernovaRemnant(SystemConfig())
    remnant.compact_remnant = BlackHole(SystemConfig())
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            snr_id = _db.insert_supernova_remnant(conn, remnant, sector_id=placed_sector_id)
            rogue_id = _db.insert_rogue_planet(conn, RoguePlanet(SystemConfig()), sector_id=placed_sector_id)
            comet_id = _db.insert_interstellar_comet(
                conn, InterstellarComet(SystemConfig()), sector_id=placed_sector_id,
            )
            stray_rogue_id = _db.insert_rogue_planet(
                conn, RoguePlanet(SystemConfig()), sector_id=unplaced_sector_id,
            )
            loose_comet_id = _db.insert_interstellar_comet(conn, InterstellarComet(SystemConfig()))
        bh_id = conn.execute(
            "SELECT compact_remnant_black_hole_id AS id FROM supernova_remnants WHERE id = ?", (snr_id,),
        ).fetchone()["id"]

        _drop_v28_placement_columns(conn)
        conn.execute("UPDATE black_holes SET sector_id = NULL WHERE id = ?", (bh_id,))
        conn.execute("DELETE FROM schema_migrations WHERE version IN (28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (27)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    half_edge_pc = ly_to_pc(11.5) / 2
    # No grid address, so placement falls back to an axis-aligned cube.
    axes = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in _db.V28_PLACED_TABLES:
            indexes = {row["Key_name"] for row in conn.execute(f"SHOW INDEX FROM {table}").fetchall()}
            assert {f"idx_{table}_center", f"idx_{table}_galactic_radius_pc"} <= indexes, table
            assert _db._has_constraint(conn, table, f"chk_{table}_placement"), table

        def placement(table, row_id):
            return conn.execute(
                f"SELECT sector_id, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc"
                f" FROM {table} WHERE id = ?",
                (row_id,),
            ).fetchone()

        for table, row_id in (
            ("supernova_remnants", snr_id), ("rogue_planets", rogue_id), ("interstellar_comets", comet_id),
        ):
            row = placement(table, row_id)
            assert row["center_x_pc"] is not None, table
            offset = (row["center_x_pc"] - center_pc[0], row["center_y_pc"] - center_pc[1],
                      row["center_z_pc"] - center_pc[2])
            for axis in axes:
                assert abs(sum(o * a for o, a in zip(offset, axis))) <= half_edge_pc + 1e-6, table
            assert row["galactic_radius_pc"] == pytest.approx(
                math.sqrt(row["center_x_pc"] ** 2 + row["center_y_pc"] ** 2 + row["center_z_pc"] ** 2)
            )

        snr = placement("supernova_remnants", snr_id)
        bh = placement("black_holes", bh_id)
        assert bh["sector_id"] == placed_sector_id
        assert (bh["center_x_pc"], bh["center_y_pc"], bh["center_z_pc"]) == (
            snr["center_x_pc"], snr["center_y_pc"], snr["center_z_pc"],
        )

        assert placement("rogue_planets", stray_rogue_id)["center_x_pc"] is None
        assert placement("interstellar_comets", loose_comet_id)["center_x_pc"] is None
    finally:
        conn.close()

    # Running it again on an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

def test_modified_at_tracks_edits_but_not_orbit_ticks(mysql_config):
    """
    v27: `modified_at` moves when a row is edited (MySQL's `ON UPDATE`)
    or one of a system's child rows changes (`touch_star_system`), but
    NOT when `advance_orbital_phases` ticks the simulation clock forward
    -- see `schema.sql`'s "v27" header note.
    """
    binary_cfg = SystemConfig()
    binary_cfg.STAR_TYPE = "G2V"
    binary_cfg.BINARY_SYSTEM = True
    binary_cfg.WIDE_BINARY = False
    binary_cfg.PLANETS = False
    binary_system = StarSystem(system_config=binary_cfg)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, binary_system, binary_cfg)
            bh_id = _db.insert_black_hole(conn, BlackHole(SystemConfig()))

        def modified(table, row_id):
            return conn.execute(f"SELECT modified_at FROM {table} WHERE id = ?", (row_id,)).fetchone()["modified_at"]

        system_before = modified("star_systems", system_id)
        bh_before = modified("black_holes", bh_id)
        time.sleep(0.05)

        counts = _db.advance_orbital_phases(conn, elapsed_years=1e5)
        assert counts["binary_mutual_orbits"] > 0
        assert counts["black_holes"] > 0
        assert modified("star_systems", system_id) == system_before
        assert modified("black_holes", bh_id) == bh_before

        with conn:
            conn.execute("UPDATE star_systems SET name = ? WHERE id = ?", ("Renamed For Test", system_id))
        system_renamed = modified("star_systems", system_id)
        assert system_renamed > system_before

        time.sleep(0.05)
        with conn:
            _db.touch_star_system(conn, system_id)
        assert modified("star_systems", system_id) > system_renamed
    finally:
        conn.close()


def test_migrate_v29_to_v30_cleans_up_surface_conditions(mysql_config):
    """
    v30 NULLs a stale atmosphere left on airless bodies and raises surface
    temperatures below the cosmic background to it, on planets and moons.
    """
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    cfg.PLANETS = True
    cfg.MOONS = True
    _db.save_system(StarSystem(system_config=cfg), cfg, config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        for table in ("planets", "moons"):
            row = conn.execute(f"SELECT id FROM {table} ORDER BY id LIMIT 1").fetchone()
            if row is None:
                continue
            conn.execute(
                f"UPDATE {table} SET atmosphere = 'None', atm_density = 1.5, atm_molar_density = 0.03,"
                " scale_height_km = 8.0, surface_temperature_k = 0.4 WHERE id = ?",
                (row["id"],),
            )
        conn.execute("DELETE FROM schema_migrations WHERE version IN (30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (29)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            stale = conn.execute(
                f"SELECT COUNT(*) AS n FROM {table} WHERE atmosphere = 'None' AND"
                " (atm_density IS NOT NULL OR atm_molar_density IS NOT NULL OR scale_height_km IS NOT NULL)"
            ).fetchone()["n"]
            assert stale == 0, table
            coldest = conn.execute(f"SELECT MIN(surface_temperature_k) AS t FROM {table}").fetchone()["t"]
            assert coldest is None or coldest >= pc.COSMIC_BACKGROUND_TEMPERATURE_K, table
    finally:
        conn.close()


def _sector_with_one_system(name):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    sector = SpaceSector(name, edge_ly=11.5)
    sector.add_system(StarSystem(system_config=cfg), position=(0.0, 0.0, 0.0), system_config=cfg)
    return sector


def test_migrate_v31_to_v32_regenerates_the_galaxy_on_the_cylindrical_grid(mysql_config):
    """
    v32 replaces shell addressing with the ring/layer/slot grid: every
    galaxy-placed sector is deleted with its systems and phenomena, a
    never-placed sector keeps everything, the old tables/columns go, and
    the skeleton is rebuilt from the stored shape (by v33, which runs
    right after).
    """
    from stellarObjects.galaxyDensity import build_galaxy_shape
    from stellarObjects.roguePlanetData import RoguePlanet

    shape = build_galaxy_shape(
        disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
        bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
    )
    _db.save_galaxy_shape(shape, edge_pc=3.526, outer_ring_index=7,
                          expected_system_count_at_density_1=20.0, config=mysql_config)
    placed_id = _db.save_sector(_sector_with_one_system("Placed"), config=mysql_config, galaxy_position={
        "center_x_pc": 1.0, "center_y_pc": 2.0, "center_z_pc": 0.5, "galactic_radius_pc": 2.29,
    })
    unplaced_id = _db.save_sector(_sector_with_one_system("Unplaced"), config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_rogue_planet(conn, RoguePlanet(SystemConfig()), sector_id=placed_id)
            kept_rogue_id = _db.insert_rogue_planet(conn, RoguePlanet(SystemConfig()), sector_id=unplaced_id)
        # Put the database back in its v31 shape.
        conn.execute("ALTER TABLE sectors DROP INDEX uq_sectors_address, DROP COLUMN ring_index, "
                     "DROP COLUMN layer_index, DROP COLUMN ring_slot_index, ADD COLUMN shell_index INT, "
                     "ADD COLUMN shell_slot_index INT, ADD UNIQUE KEY shell_index (shell_index, shell_slot_index)")
        conn.execute("UPDATE sectors SET shell_index = 0, shell_slot_index = 1 WHERE id = ?", (placed_id,))
        conn.execute("CREATE TABLE sector_vertices (sector_id BIGINT UNSIGNED NOT NULL, "
                     "FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE)")
        conn.execute("INSERT INTO sector_vertices (sector_id) VALUES (?)", (placed_id,))
        conn.execute("CREATE TABLE galaxy_shell_band (shell_index INT)")
        conn.execute("ALTER TABLE galaxy_shape CHANGE COLUMN outer_ring_index outer_shell_index INT NOT NULL")
        conn.execute("DELETE FROM galaxy_layer")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (31)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        sector_ids = {row["id"] for row in conn.execute("SELECT id FROM sectors").fetchall()}
        assert sector_ids == {unplaced_id}
        system_sectors = [row["sector_id"] for row in conn.execute("SELECT sector_id FROM star_systems").fetchall()]
        assert system_sectors == [unplaced_id]
        rogue_ids = [row["id"] for row in conn.execute("SELECT id FROM rogue_planets").fetchall()]
        assert rogue_ids == [kept_rogue_id]

        assert not _db._has_column(conn, "sectors", "shell_index")
        assert _db._has_column(conn, "sectors", "ring_index")
        assert _db._has_index(conn, "sectors", "uq_sectors_address")
        for table in ("sector_vertices", "galaxy_shell_band"):
            assert conn.execute(
                "SELECT 1 FROM information_schema.tables WHERE table_schema = DATABASE() AND table_name = ?",
                (table,),
            ).fetchone() is None, table

        # v33 then rebuilds the skeleton per layer at the standard edge.
        skeleton = _db.get_galaxy_shape(conn)
        assert skeleton.edge_pc == program_constants.DEFAULT_SECTOR_EDGE_PC
        assert _db.get_galaxy_layers(conn)
        assert _db.get_galaxy_layer_outer_ring(conn, 0) == skeleton.outer_ring_index
    finally:
        conn.close()

    # Running it again on an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION


def test_migrate_v32_to_v33_moves_to_the_sector_standard(mysql_config):
    """
    v33 moves the galaxy to 4 pc sectors, round(2*pi*(i + 1/2)) slots per
    ring and a per-layer skeleton: every galaxy-placed sector is deleted
    with its systems and phenomena, a never-placed sector keeps
    everything, `galaxy_ring_band` goes, and `galaxy_layer` is rebuilt
    from the stored shape at the standard edge.
    """
    from stellarObjects.galaxyDensity import build_galaxy_shape
    from stellarObjects.galaxySkeleton import build_layer_extents, column_extents
    from stellarObjects.roguePlanetData import RoguePlanet

    shape = build_galaxy_shape(
        disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
        bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
    )
    _db.save_galaxy_shape(shape, edge_pc=3.526, outer_ring_index=7,
                          expected_system_count_at_density_1=20.0, config=mysql_config)
    placed_id = _db.save_sector(_sector_with_one_system("Placed"), config=mysql_config, galaxy_position={
        "center_x_pc": 1.0, "center_y_pc": 2.0, "center_z_pc": 0.5, "galactic_radius_pc": 2.29,
        "ring_index": 0, "layer_index": 0, "ring_slot_index": 1,
    })
    unplaced_id = _db.save_sector(_sector_with_one_system("Unplaced"), config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_rogue_planet(conn, RoguePlanet(SystemConfig()), sector_id=placed_id)
            kept_rogue_id = _db.insert_rogue_planet(conn, RoguePlanet(SystemConfig()), sector_id=unplaced_id)
        # Put the database back in its v32 shape.
        conn.execute("CREATE TABLE galaxy_ring_band (ring_index INT NOT NULL PRIMARY KEY, "
                     "layer_index_min INT NOT NULL, layer_index_max INT NOT NULL)")
        conn.execute("INSERT INTO galaxy_ring_band VALUES (0, -1, 1)")
        conn.execute("DELETE FROM galaxy_layer")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (32)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        sector_ids = {row["id"] for row in conn.execute("SELECT id FROM sectors").fetchall()}
        assert sector_ids == {unplaced_id}
        system_sectors = [row["sector_id"] for row in conn.execute("SELECT sector_id FROM star_systems").fetchall()]
        assert system_sectors == [unplaced_id]
        rogue_ids = [row["id"] for row in conn.execute("SELECT id FROM rogue_planets").fetchall()]
        assert rogue_ids == [kept_rogue_id]
        assert conn.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_schema = DATABASE() AND table_name = ?",
            ("galaxy_ring_band",),
        ).fetchone() is None

        skeleton = _db.get_galaxy_shape(conn)
        assert skeleton.edge_pc == program_constants.DEFAULT_SECTOR_EDGE_PC
        assert skeleton.expected_system_count_at_density_1 == pytest.approx(
            SpaceSector("x", edge_ly=program_constants.DEFAULT_SECTOR_EDGE_LY).expected_system_count()
        )
        expected, outer, _confirmed = build_layer_extents(
            shape, program_constants.DEFAULT_SECTOR_EDGE_PC, 1.0 / skeleton.expected_system_count_at_density_1,
        )
        assert _db.get_galaxy_layers(conn) == expected
        assert skeleton.outer_ring_index == outer
        columns = [
            (row["ring_index"], row["layer_index_min"], row["layer_index_max"])
            for row in conn.execute("SELECT * FROM galaxy_column ORDER BY ring_index").fetchall()
        ]
        assert columns == column_extents(expected)
        assert _db.get_galaxy_column(conn, 0) == columns[0][1:]
        bounds = _db.get_galaxy_bounds(conn)
        assert bounds.outer_ring == dict(expected)
        assert bounds.contains(outer, 0) and not bounds.contains(outer + 1, 0)
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION


def test_migrate_v33_to_v34_drops_the_planet_and_moon_name_registry(mysql_config):
    """v34 names planets and moons from their system (`bodyNames.py`), so
    `body_name_registry` goes; the rest of the database is untouched."""
    system_id = _db.save_system(StarSystem(SystemConfig()), SystemConfig(), config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        conn.execute(
            "CREATE TABLE body_name_registry (id BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY, "
            "base_name VARCHAR(255) NOT NULL, occurrence_count INT NOT NULL, "
            "first_body_kind VARCHAR(8) NOT NULL, first_body_id BIGINT UNSIGNED NOT NULL, suffix_index INT)"
        )
        conn.execute("DELETE FROM schema_migrations WHERE version IN (34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (33)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_schema = DATABASE() AND table_name = ?",
            ("body_name_registry",),
        ).fetchone() is None
        assert conn.execute("SELECT id FROM star_systems").fetchone()["id"] == system_id
    finally:
        conn.close()


def test_migrate_v34_to_v35_deletes_sectors_in_changed_rings_and_keeps_the_skeleton(mysql_config):
    """v35's master-wedge slot rule changes what a stored slot index
    means in every ring whose count changed, so those sectors go with
    their systems and phenomena. A sector in an unchanged ring (ring 1
    holds 9 slots under both rules) and a never-placed sector keep
    everything, and the skeleton stays."""
    from stellarObjects.galaxyDensity import build_galaxy_shape
    from stellarObjects.roguePlanetData import RoguePlanet

    shape = build_galaxy_shape(
        disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
        bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
    )
    _db.save_galaxy_shape(shape, edge_pc=4.0, outer_ring_index=7,
                          expected_system_count_at_density_1=20.0, config=mysql_config)

    def placed(name, ring):
        return _db.save_sector(_sector_with_one_system(name), config=mysql_config, galaxy_position={
            "center_x_pc": 4.0 * ring + 2.0, "center_y_pc": 0.1, "center_z_pc": 0.0,
            "galactic_radius_pc": 4.0 * ring + 2.0, "ring_index": ring, "layer_index": 0, "ring_slot_index": 0,
        })

    changed_id = placed("Changed", 2)  # 16 slots before, 15 after
    unchanged_id = placed("Unchanged", 1)
    unplaced_id = _db.save_sector(_sector_with_one_system("Unplaced"), config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_rogue_planet(conn, RoguePlanet(SystemConfig()), sector_id=changed_id)
            kept_rogue_id = _db.insert_rogue_planet(conn, RoguePlanet(SystemConfig()), sector_id=unplaced_id)
        layers_before = _db.get_galaxy_layers(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (34)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        kept = {unchanged_id, unplaced_id}
        assert {row["id"] for row in conn.execute("SELECT id FROM sectors").fetchall()} == kept
        assert {row["sector_id"] for row in conn.execute("SELECT sector_id FROM star_systems").fetchall()} == kept
        assert [row["id"] for row in conn.execute("SELECT id FROM rogue_planets").fetchall()] == [kept_rogue_id]
        assert _db.get_galaxy_layers(conn) == layers_before
        assert _db.get_galaxy_shape(conn).outer_ring_index == 7
    finally:
        conn.close()


def test_migrate_v35_to_v36_adds_black_hole_mass_classes(mysql_config):
    """v36 adds `black_holes.mass_class`, filled from each row's mass."""
    from stellarObjects.compactRemnant import BlackHole

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            ids = [_db.insert_black_hole(conn, BlackHole(SystemConfig())) for _ in range(3)]
        for row_id, mass in zip(ids, (10.0, 5e3, 4.3e6)):
            conn.execute("UPDATE black_holes SET mass_solar = ? WHERE id = ?", (mass, row_id))
        conn.execute("ALTER TABLE black_holes DROP CONSTRAINT chk_black_holes_mass_class")
        conn.execute("ALTER TABLE black_holes DROP COLUMN mass_class")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (35)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        classes = [conn.execute("SELECT mass_class FROM black_holes WHERE id = ?", (i,)).fetchone()["mass_class"]
                   for i in ids]
        assert classes == ["stellar", "intermediate", "supermassive"]
        assert _db._has_constraint(conn, "black_holes", "chk_black_holes_mass_class")
    finally:
        conn.close()


def test_migrate_v36_to_v37_adds_rogue_mass_bins_and_runaway_columns(mysql_config):
    """v37 adds `rogue_planets.mass_bin` (filled from mass) and
    `star_systems.runaway_class`/`runaway_speed_kms`."""
    from stellarObjects.roguePlanetData import RoguePlanet

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            ids = [_db.insert_rogue_planet(conn, RoguePlanet(SystemConfig())) for _ in range(5)]
        earth = pc.EARTH_MASS_TO_KG
        masses = (1.0 * earth, 10.0 * earth, 100.0 * earth, 1000.0 * earth, 30 * pc.JUPITER_MASS_TO_KG)
        for row_id, mass in zip(ids, masses):
            conn.execute("UPDATE rogue_planets SET mass_kg = ? WHERE id = ?", (mass, row_id))
        conn.execute("ALTER TABLE rogue_planets DROP COLUMN mass_bin")
        conn.execute("ALTER TABLE star_systems DROP COLUMN runaway_class, DROP COLUMN runaway_speed_kms")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (36)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        bins = [conn.execute("SELECT mass_bin FROM rogue_planets WHERE id = ?", (i,)).fetchone()["mass_bin"]
                for i in ids]
        assert bins == ["terrestrial", "sub-neptune", "saturn", "jupiter", "brown-dwarf"]
        assert _db._has_column(conn, "star_systems", "runaway_class")
        assert _db._has_column(conn, "star_systems", "runaway_speed_kms")
    finally:
        conn.close()


def test_runaway_flags_round_trip(mysql_config):
    system = StarSystem(SystemConfig())
    system.runaway_class, system.runaway_speed_kms = "runaway", 55.0
    system_id = _db.save_system(system, SystemConfig(), config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        loaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()
    assert (loaded.runaway_class, loaded.runaway_speed_kms) == ("runaway", 55.0)
    assert StarSystem.from_dict(loaded.to_dict()).runaway_class == "runaway"


def test_migrate_v37_to_v38_classes_existing_nebulae_remnants_and_fields(mysql_config):
    """v38 adds letter classes and contents, inferring them for old rows,
    and widens `nebula_type` to the `diffuse` family."""
    from stellarObjects.asteroidFieldData import AsteroidField
    from stellarObjects.nebulaData import Nebula
    from stellarObjects.supernovaRemnantData import SupernovaRemnant

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            nebula_id = _db.insert_nebula(conn, Nebula(SystemConfig(), nebula_class="F"))
            remnant = SupernovaRemnant(SystemConfig())
            remnant_id = _db.insert_supernova_remnant(conn, remnant)
            field_id = _db.insert_asteroid_field(conn, AsteroidField(SystemConfig()))
        conn.execute("UPDATE nebulae SET radius_ly = 5.0 WHERE id = ?", (nebula_id,))
        conn.execute("UPDATE supernova_remnants SET morphology = 'plerion', progenitor_type = 'core-collapse' "
                     "WHERE id = ?", (remnant_id,))
        conn.execute("UPDATE asteroid_fields SET density = 'dense', radius_ly = 0.5 WHERE id = ?", (field_id,))
        conn.execute("ALTER TABLE nebulae DROP CONSTRAINT chk_nebulae_type, DROP COLUMN nebula_class, "
                     "DROP COLUMN dominant_species, DROP COLUMN density_cm3, DROP COLUMN temperature_k, "
                     "DROP COLUMN extinction_av")
        conn.execute("ALTER TABLE nebulae ADD CONSTRAINT nebulae_chk_old CHECK "
                     "(nebula_type IN ('emission', 'reflection', 'planetary', 'dark'))")
        conn.execute("ALTER TABLE supernova_remnants DROP COLUMN remnant_class, DROP COLUMN dominant_species, "
                     "DROP COLUMN density_cm3, DROP COLUMN temperature_k, DROP COLUMN extinction_av")
        conn.execute("ALTER TABLE asteroid_fields DROP COLUMN field_class, DROP COLUMN composition_family")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (37)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        nebula = conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
        assert nebula["nebula_class"] == "F"
        assert 100 <= nebula["density_cm3"] <= 1000
        remnant_row = conn.execute("SELECT * FROM supernova_remnants WHERE id = ?", (remnant_id,)).fetchone()
        assert remnant_row["remnant_class"] == "T"
        assert remnant_row["dominant_species"]
        field = conn.execute("SELECT * FROM asteroid_fields WHERE id = ?", (field_id,)).fetchone()
        assert (field["composition_family"], field["field_class"]) == ("mixed", "S4")
        with conn:
            _db.insert_nebula(conn, Nebula(SystemConfig(), nebula_class="A"))
    finally:
        conn.close()


def test_classed_phenomena_round_trip_their_contents(mysql_config):
    from stellarObjects.asteroidFieldData import AsteroidField
    from stellarObjects.nebulaData import Nebula

    nebula = Nebula(SystemConfig(), nebula_class="Q")
    field = AsteroidField(SystemConfig())
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            nebula_id = _db.insert_nebula(conn, nebula)
            field_id = _db.insert_asteroid_field(conn, field)
        row = conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
        assert (row["nebula_class"], row["nebula_type"]) == ("Q", "dark")
        assert row["temperature_k"] == pytest.approx(nebula.temperature_k)
        field_row = conn.execute("SELECT * FROM asteroid_fields WHERE id = ?", (field_id,)).fetchone()
        assert field_row["field_class"] == field.field_class
    finally:
        conn.close()


def test_innermost_container_picks_the_smallest_cloud_holding_the_point():
    big = {"column": "inside_nebula_id", "id": 1, "center": (0.0, 0.0, 0.0), "radius_pc": 10.0}
    small = {"column": "inside_remnant_id", "id": 2, "center": (1.0, 0.0, 0.0), "radius_pc": 2.0}
    assert _db.innermost_container((1.5, 0.0, 0.0), [big, small]) is small
    assert _db.innermost_container((5.0, 0.0, 0.0), [big, small]) is big
    assert _db.innermost_container((20.0, 0.0, 0.0), [big, small]) is None
    # A nebula only nests in a larger cloud, and never in itself.
    assert _db.innermost_container((1.0, 0.0, 0.0), [big, small], own_radius_pc=3.0) is big
    assert _db.innermost_container((0.0, 0.0, 0.0), [big], own_radius_pc=0.0,
                                   own=("inside_nebula_id", 1)) is None


def _placed_nebula(conn, sector_id, center_pc, radius_ly, nebula_class="D"):
    from stellarObjects.nebulaData import Nebula
    nebula = Nebula(SystemConfig(), nebula_class=nebula_class)
    nebula.radius_ly = radius_ly
    x, y, z = center_pc
    return _db.insert_nebula(conn, nebula, sector_id=sector_id, placement={
        "center_x_pc": x, "center_y_pc": y, "center_z_pc": z, "galactic_radius_pc": math.hypot(x, y, z),
    })


def test_systems_inside_a_nebula_point_at_it_and_the_innermost_wins(mysql_config):
    sector_id = _db.save_sector(_sector_with_one_system("Cloudy"), config=mysql_config, galaxy_position={
        "center_x_pc": 100.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 100.0,
    })
    far_id = _db.save_sector(_sector_with_one_system("Clear"), config=mysql_config, galaxy_position={
        "center_x_pc": 400.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 400.0,
    })
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            big = _placed_nebula(conn, sector_id, (110.0, 0.0, 0.0), radius_ly=100.0, nebula_class="E")
        row = conn.execute("SELECT inside_nebula_id, inside_remnant_id FROM star_systems WHERE sector_id = ?",
                           (sector_id,)).fetchone()
        assert (row["inside_nebula_id"], row["inside_remnant_id"]) == (big, None)
        far = conn.execute("SELECT inside_nebula_id FROM star_systems WHERE sector_id = ?", (far_id,)).fetchone()
        assert far["inside_nebula_id"] is None

        with conn:
            small = _placed_nebula(conn, sector_id, (100.5, 0.0, 0.0), radius_ly=3.0, nebula_class="C")
        row = conn.execute("SELECT inside_nebula_id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()
        assert row["inside_nebula_id"] == small
        nested = conn.execute("SELECT inside_nebula_id FROM nebulae WHERE id = ?", (small,)).fetchone()
        assert nested["inside_nebula_id"] == big

        with conn:
            conn.execute("DELETE FROM nebulae WHERE id = ?", (small,))
            _db.refresh_containment(conn, [sector_id])
        row = conn.execute("SELECT inside_nebula_id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()
        assert row["inside_nebula_id"] == big
        import queryDb
        inside = queryDb.sector_detail(conn, sector_id)["systems"][0]["inside"]
        assert (inside["type"], inside["id"], inside["class"]) == ("nebula", big, "E")
    finally:
        conn.close()


def test_migrate_v38_to_v39_fills_containment(mysql_config):
    sector_id = _db.save_sector(_sector_with_one_system("Old Cloud"), config=mysql_config, galaxy_position={
        "center_x_pc": 50.0, "center_y_pc": 50.0, "center_z_pc": 0.0, "galactic_radius_pc": 70.7,
    })
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            nebula_id = _placed_nebula(conn, sector_id, (50.0, 50.0, 1.0), radius_ly=30.0)
        for table in _db.CONTAINABLE_TABLES:
            conn.execute(f"ALTER TABLE {table} DROP FOREIGN KEY fk_{table}_inside_nebula, "
                         f"DROP FOREIGN KEY fk_{table}_inside_remnant")
            conn.execute(f"ALTER TABLE {table} DROP COLUMN inside_nebula_id, DROP COLUMN inside_remnant_id")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (38)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT inside_nebula_id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()
        assert row["inside_nebula_id"] == nebula_id
    finally:
        conn.close()


def _named_system(name):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    system.name = name
    return system, cfg


def _remnant_with_core(name):
    from stellarObjects.supernovaRemnantData import SupernovaRemnant
    for _ in range(200):
        remnant = SupernovaRemnant(SystemConfig(), name=name)
        if remnant.compact_remnant is not None:
            return remnant
    pytest.fail("no supernova remnant with a detectable core")


def test_phenomena_share_the_system_name_registry(mysql_config):
    """GEN.13 (v40): a nebula whose name clashes with a system's is
    decorated like a second system would be, and the first holder is
    renamed whichever kind it is."""
    from stellarObjects.nebulaData import Nebula

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system, cfg = _named_system("Kelvaro")
            system_id = _db.insert_star_system(conn, system, cfg)
            nebula = Nebula(SystemConfig(), name="Kelvaro")
            nebula_id = _db.insert_nebula(conn, nebula)
        assert conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"] == "Alpha Kelvaro"
        assert conn.execute("SELECT name FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()["name"] == "Beta Kelvaro"
        assert nebula.name == "Beta Kelvaro"

        with conn:
            first = _db.insert_rogue_planet(conn, __import__(
                "stellarObjects.roguePlanetData", fromlist=["RoguePlanet"]).RoguePlanet(SystemConfig(), name="Ossandre"))
            system, cfg = _named_system("Ossandre")
            _db.insert_star_system(conn, system, cfg)
        assert conn.execute("SELECT name FROM rogue_planets WHERE id = ?", (first,)).fetchone()["name"] == "Alpha Ossandre"
        assert system.name == "Beta Ossandre"
        assert _db.name_in_use(conn, "Alpha Ossandre") == "rogue_planets"
    finally:
        conn.close()


def test_a_renamed_remnant_takes_its_core_along(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            remnant_id = _db.insert_supernova_remnant(conn, _remnant_with_core("Thessavel"))
            _db.insert_nebula(conn, __import__(
                "stellarObjects.nebulaData", fromlist=["Nebula"]).Nebula(SystemConfig(), name="Thessavel"))
        row = conn.execute("SELECT * FROM supernova_remnants WHERE id = ?", (remnant_id,)).fetchone()
        assert row["name"] == "Alpha Thessavel"
        core_table = "black_holes" if row["compact_remnant_black_hole_id"] else "neutron_stars"
        core_id = row["compact_remnant_black_hole_id"] or row["compact_remnant_neutron_star_id"]
        core = conn.execute(f"SELECT name FROM {core_table} WHERE id = ?", (core_id,)).fetchone()
        assert core["name"] == "Alpha Thessavel Core"
        registry = conn.execute(
            "SELECT occurrence_count, first_object_table FROM system_name_registry WHERE base_name = 'Thessavel'"
        ).fetchone()
        # The core isn't registered on its own.
        assert (registry["occurrence_count"], registry["first_object_table"]) == (2, "supernova_remnants")
    finally:
        conn.close()


def test_comets_carry_designations_that_follow_their_host(mysql_config):
    from stellarObjects.cometData import PERIODIC_COMET_MAX_PERIOD_YEARS

    system, cfg = _make_system_with_comets()
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        rows = conn.execute(
            "SELECT name, orbit_type, orbital_period_years FROM comets WHERE star_system_id = ? ORDER BY id",
            (system_id,),
        ).fetchall()
        for n, row in enumerate(rows, start=1):
            periodic = row["orbit_type"] == "elliptical" and row["orbital_period_years"] < PERIODIC_COMET_MAX_PERIOD_YEARS
            assert row["name"] == f"{'P' if periodic else 'C'}/{system.name}-{n}"
        with conn:
            _db.rename_star_system(conn, system_id, "Neraloth")
        renamed = conn.execute("SELECT name FROM comets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()
        assert [r["name"][2:] for r in renamed] == [f"Neraloth-{n}" for n in range(1, len(rows) + 1)]
    finally:
        conn.close()


def test_interstellar_comets_and_asteroid_fields_are_designated_by_sector(mysql_config):
    from stellarObjects.asteroidFieldData import AsteroidField
    from stellarObjects.roguePlanetData import InterstellarComet

    sector_id = _db.save_sector(_sector_with_one_system("Designated"), config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_interstellar_comet(conn, InterstellarComet(SystemConfig()), sector_id=sector_id)
            _db.insert_interstellar_comet(conn, InterstellarComet(SystemConfig()), sector_id=sector_id)
            field = AsteroidField(SystemConfig())
            _db.insert_asteroid_field(conn, field, sector_id=sector_id)
        names = [r["name"] for r in conn.execute(
            "SELECT name FROM interstellar_comets WHERE sector_id = ? ORDER BY id", (sector_id,)).fetchall()]
        sector_name = conn.execute("SELECT name FROM sectors WHERE id = ?", (sector_id,)).fetchone()["name"]
        assert names == [f"I/{sector_name}-1", f"I/{sector_name}-2"]
        assert field.name == f"AF {field.field_class}-{sector_name}-01"
    finally:
        conn.close()


def test_migrate_v39_to_v40_registers_and_designates_existing_rows(mysql_config):
    system, cfg = _make_system_with_comets()
    system.name = "Morrowen"
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
            conn.execute("UPDATE comets SET name = 'Old Comet' WHERE star_system_id = ?", (system_id,))
            from stellarObjects.nebulaData import Nebula
            nebula_id = _db.insert_nebula(conn, Nebula(SystemConfig(), name="Unregistered"))
            # Before v40 a phenomenon's name was never checked.
            conn.execute("UPDATE nebulae SET name = 'Morrowen' WHERE id = ?", (nebula_id,))
            for table in _db.NAMED_PHENOMENON_TABLES:
                conn.execute(f"ALTER TABLE {table} DROP KEY idx_{table}_name")
            conn.execute("DELETE FROM system_name_registry WHERE first_star_system_id IS NULL")
            conn.execute("ALTER TABLE system_name_registry DROP COLUMN first_object_table, DROP COLUMN first_object_id, "
                         "MODIFY first_star_system_id BIGINT UNSIGNED NOT NULL")
            conn.execute("DELETE FROM schema_migrations WHERE version IN (40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
            conn.execute("INSERT INTO schema_migrations (version) VALUES (39)")
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"] == "Alpha Morrowen"
        assert conn.execute("SELECT name FROM nebulae").fetchone()["name"] == "Beta Morrowen"
        comet_names = [r["name"] for r in conn.execute(
            "SELECT name FROM comets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
        assert all(name[:2] in ("P/", "C/") and name[2:].startswith("Alpha Morrowen-") for name in comet_names)
    finally:
        conn.close()


def test_system_grid_nearest_matches_brute_force():
    import random as _random
    rng = _random.Random(7)
    systems = [(i, (rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(-5, 5))) for i in range(300)]
    grid = _db._SystemGrid(systems)
    for _ in range(50):
        point = (rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(-5, 5))
        brute = sorted((math.dist(point, p), i) for i, p in systems if math.dist(point, p) <= 4.0)[:3]
        assert grid.nearest(point) == brute


def test_galaxy_to_local_undoes_local_to_galaxy():
    from stellarObjects.galaxyGeometry import galaxy_to_local_pc, local_to_galaxy_pc
    center = (30.0, -40.0, 2.0)
    offset = (1.2, -0.7, 0.3)
    assert galaxy_to_local_pc(center, local_to_galaxy_pc(center, offset)) == pytest.approx(offset)


def _sector_with_systems(name, positions_ly):
    sector = SpaceSector(name, edge_ly=13.0)
    for position in positions_ly:
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = False
        cfg.BINARY_SYSTEM = False
        sector.add_system(StarSystem(system_config=cfg), position=position, system_config=cfg)
    return sector


def test_nearest_systems_cross_sector_boundaries(mysql_config):
    """UX.18 (v41): a system at a sector's edge lists a system just
    across the boundary once that sector is generated."""
    import queryDb
    from stellarObjects.utils import pc_to_ly

    edge_ly = pc_to_ly(4.0)
    first = _db.save_sector(_sector_with_systems("Westmark", [(edge_ly / 2 - 0.5, 0.0, 0.0), (-5.0, 0.0, 0.0)]),
                            config=mysql_config, galaxy_position={
        "center_x_pc": 100.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 100.0,
    })
    conn = _db.get_connection(mysql_config)
    try:
        border_id, far_id = [r["id"] for r in conn.execute(
            "SELECT id FROM star_systems WHERE sector_id = ? ORDER BY position_x_mpc DESC", (first,)).fetchall()]
        before = queryDb.nearest_systems(conn, "star_systems", [border_id])[border_id]
        assert [n["id"] for n in before] == [far_id]
    finally:
        conn.close()

    second = _db.save_sector(_sector_with_systems("Eastmark", [(-edge_ly / 2 + 0.5, 0.0, 0.0)]),
                             config=mysql_config, galaxy_position={
        "center_x_pc": 104.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 104.0,
    })
    conn = _db.get_connection(mysql_config)
    try:
        across = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (second,)).fetchone()["id"]
        after = queryDb.nearest_systems(conn, "star_systems", [border_id])[border_id]
        assert [n["id"] for n in after] == [across, far_id]
        assert after[0]["distance_ly"] == pytest.approx(1.0, abs=0.01)
        own = queryDb.nearest_systems(conn, "star_systems", [across])[across]
        assert own[0]["id"] == border_id
        detail = queryDb.sector_detail(conn, first)
        by_id = {system["id"]: system for system in detail["systems"]}
        assert by_id[border_id]["nearest"][0]["id"] == across
    finally:
        conn.close()


def test_phenomena_store_their_octant_and_nearest_systems(mysql_config):
    import queryDb
    sector_id = _db.save_sector(_sector_with_systems("Octmark", [(1.0, 1.0, 1.0)]), config=mysql_config,
                                galaxy_position={
        "center_x_pc": 0.0, "center_y_pc": 200.0, "center_z_pc": 0.0, "galactic_radius_pc": 200.0,
    })
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            # Local +X points away from the galactic axis (+y here), local +Y along -x.
            nebula_id = _placed_nebula(conn, sector_id, (0.5, 200.5, -0.5), radius_ly=0.5)
            _db.refresh_nearest_systems(conn, [sector_id])
        quadrant = conn.execute("SELECT quadrant FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()["quadrant"]
        from stellarObjects.spaceSector import classify_octant
        assert quadrant == classify_octant((0.5, -0.5, -0.5))[0]
        listed = [p for p in queryDb.phenomena_near_sector(conn, sector_id) if p["type"] == "nebula"][0]
        assert listed["octant"] == quadrant
        assert len(listed["nearest"]) == 1
    finally:
        conn.close()


def test_migrate_v40_to_v41_fills_nearest_systems(mysql_config):
    sector_id = _db.save_sector(_sector_with_systems("Oldmark", [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)]),
                                config=mysql_config, galaxy_position={
        "center_x_pc": 60.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 60.0,
    })
    conn = _db.get_connection(mysql_config)
    try:
        conn.execute("DROP TABLE nearest_systems")
        for table in _db.PLACED_PHENOMENON_TABLES:
            # MariaDB drops a column-level CHECK with its column.
            conn.execute(f"ALTER TABLE {table} DROP COLUMN quadrant")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (41, 42, 43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (40)")
        conn.commit()
    finally:
        conn.close()

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = conn.execute("SELECT COUNT(*) AS n FROM nearest_systems WHERE sector_id = ?", (sector_id,)).fetchone()
        assert rows["n"] == 2
    finally:
        conn.close()


def _bright_row(ring=3, layer=0, slot=1, luminosity_sol=800.0):
    return (ring, layer, slot, 12000, 3000, 0, "young", "B2V", "V", 1.4e31, 3.0e6, 22000.0,
            luminosity_sol * 3.828e26,
            0.02, 0.03, 7.0, 0.03, 42)


def test_bright_stars_store_and_clear(mysql_config):
    """Bright-star pre-placement storage (v43): bulk insert, per-sector
    lookup brightest first, fill link, and a plan re-run clearing it."""
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc,"
                         " bulge_scale_radius_pc, bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude,"
                         " spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
                         " expected_system_count_at_density_1, outer_ring_index)"
                         " VALUES (1, 1, 1, 1, 1, 2, 0.2, 0.3, 1, 0, 1, 4, 10, 5)")
            assert _db.bright_star_scatter_settings(conn) is None
            written = _db.insert_bright_stars(conn, [_bright_row(luminosity_sol=600.0), _bright_row(),
                                                     _bright_row(slot=2)], batch_size=2)
            _db.record_bright_star_scatter(conn, 500.0, 7)
        assert written == 3
        assert _db.bright_star_scatter_settings(conn) == (500.0, 7)
        found = _db.bright_stars_for_sector(conn, 3, 0, 1)
        assert len(found) == 2 and found[0]["luminosity_w"] > found[1]["luminosity_w"]

        system, cfg = _named_system("Beaconholm")
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
            _db.mark_bright_star_filled(conn, found[0]["id"], system_id)
        assert len(_db.bright_stars_for_sector(conn, 3, 0, 1)) == 1

        _db.clear_bright_stars(conn)
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 0
        assert _db.bright_star_scatter_settings(conn) is None
    finally:
        conn.close()


def test_bright_star_web_queries(mysql_config):
    """The Stats page's placed/filled counts, the sector page's per-cell and
    per-box listings (unfilled only by default) and the Generate page's
    scatter status all agree with what was stored."""
    import adminStats
    import queryDb

    conn = _db.get_connection(mysql_config)
    try:
        assert adminStats.bright_star_counts(conn) == {"placed": 0, "filled": 0, "unfilled": 0}
        status = queryDb.bright_star_scatter_status(conn)
        assert status["scattered"] is False and status["min_luminosity_sol"] is None and status["seed"] is None
        assert status["default_min_luminosity_sol"] == program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL == 1000.0
        assert queryDb.bright_stars_in_sector(conn, 3, 0, 1) == []
        with conn:
            conn.execute("INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc,"
                         " bulge_scale_radius_pc, bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude,"
                         " spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
                         " expected_system_count_at_density_1, outer_ring_index)"
                         " VALUES (1, 1, 1, 1, 1, 2, 0.2, 0.3, 1, 0, 1, 4, 10, 5)")
            _db.insert_bright_stars(conn, [_bright_row(luminosity_sol=600.0), _bright_row(),
                                           _bright_row(slot=2)], batch_size=2)
            _db.record_bright_star_scatter(conn, 100.0, 9)
        listed = queryDb.bright_stars_in_sector(conn, 3, 0, 1)
        assert [s["luminosity_sol"] for s in listed] == pytest.approx([800.0, 600.0], rel=1e-2)
        assert listed[0]["x"] == pytest.approx(12.0) and listed[0]["yerkes_class"] == "V"
        assert listed[0]["system_id"] is None

        system, cfg = _named_system("Lanternfall")
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
            _db.mark_bright_star_filled(conn, listed[0]["id"], system_id)
        assert adminStats.bright_star_counts(conn) == {"placed": 3, "filled": 1, "unfilled": 2}
        assert [s["id"] for s in queryDb.bright_stars_in_sector(conn, 3, 0, 1)] == [listed[1]["id"]]
        every = queryDb.bright_stars_in_sector(conn, 3, 0, 1, unfilled_only=False)
        assert [(s["id"], s["system_id"]) for s in every] == [(listed[0]["id"], system_id), (listed[1]["id"], None)]

        lo, hi = (11.0, 2.0, -1.0), (13.0, 4.0, 1.0)
        boxed = queryDb.galaxy_bright_stars_in_box(conn, lo, hi, 4.0)
        assert system_id in [s["system_id"] for s in boxed]
        assert queryDb.galaxy_bright_stars_in_box(conn, lo, hi, 4.0, unfilled_only=True) \
            == [s for s in boxed if s["system_id"] is None]

        status = queryDb.bright_star_scatter_status(conn)
        assert status == {"scattered": True, "min_luminosity_sol": 100.0, "seed": 9,
                          "default_min_luminosity_sol": 1000.0}

        # A re-scatter empties the table and restarts the ids, so the id
        # span still counts exactly.
        _db.clear_bright_stars(conn)
        with conn:
            _db.insert_bright_stars(conn, [_bright_row(slot=s) for s in range(5)])
        assert adminStats.bright_star_counts(conn) == {"placed": 5, "filled": 0, "unfilled": 5}
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 5
    finally:
        conn.close()


def test_migrate_v42_to_v43_adds_bright_star_storage(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        conn.execute("DROP TABLE bright_stars")
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN bright_star_min_luminosity_sol, DROP COLUMN bright_star_seed")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (43, 44, 45, 46, 47, 48, 49, 50)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (42)")
        conn.commit()
    finally:
        conn.close()
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION
    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 0
        assert _db._has_column(conn, "galaxy_shape", "bright_star_seed")
    finally:
        conn.close()


def test_migrate_v48_to_v49_adds_bright_star_blocks(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        conn.execute("DROP TABLE bright_star_blocks")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (49, 50)")
        conn.commit()
    finally:
        conn.close()
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION
    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_star_blocks").fetchone()["n"] == 0
    finally:
        conn.close()
