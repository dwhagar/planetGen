# tests/test_db_migration.py

"""
Tests for `stellarObjects._db.migrate_database` -- the schema v1 -> v2
converter (moons split out of the shared `planets` table into their own
`moons` table, see `schema.sql`'s "v2" header note). This is the one
piece of `_db.py` that runs against a real user's already-deployed
database rather than freshly generated data, so it's worth verifying
directly rather than only smoke-testing it by hand.

`V1_SCHEMA_SQL` is the literal schema.sql DDL from immediately before the
v2 change -- kept here rather than derived from the current schema.sql --
so these tests build a synthetic "old" database exactly the shape a real
v1 deployment's database would be, independent of any future edits to the
current schema.
"""

import os
import sqlite3

import pytest

from stellarObjects import _db

V1_SCHEMA_SQL = """
PRAGMA foreign_keys = ON;
PRAGMA user_version = 1;

CREATE TABLE IF NOT EXISTS sectors (
    id        INTEGER PRIMARY KEY,
    name      TEXT NOT NULL,
    edge_mpc  REAL NOT NULL
);

CREATE TABLE IF NOT EXISTS system_configs (
    id                INTEGER PRIMARY KEY,
    markdown          INTEGER NOT NULL DEFAULT 0 CHECK (markdown IN (0, 1)),
    habitable_world   INTEGER CHECK (habitable_world IN (0, 1)),
    asteroid_belt     INTEGER CHECK (asteroid_belt IN (0, 1)),
    large_star        INTEGER CHECK (large_star IN (0, 1)),
    moons             INTEGER CHECK (moons IN (0, 1)),
    max_planets       INTEGER CHECK (max_planets IN (0, 1)),
    planets           INTEGER CHECK (planets IN (0, 1)),
    star_type         TEXT,
    name              TEXT,
    age               TEXT CHECK (age IN ('young', 'old')),
    intelligent_life  INTEGER CHECK (intelligent_life IN (0, 1)),
    binary_system     INTEGER CHECK (binary_system IN (0, 1)),
    num_orbits        INTEGER
);

CREATE TABLE IF NOT EXISTS system_config_slots (
    id            INTEGER PRIMARY KEY,
    config_id     INTEGER NOT NULL REFERENCES system_configs(id) ON DELETE CASCADE,
    orbit_index   INTEGER NOT NULL,
    type          TEXT CHECK (type IN ('planet', 'asteroid_belt')),
    planet_class  TEXT,
    moons         INTEGER
);

CREATE TABLE IF NOT EXISTS star_systems (
    id                     INTEGER PRIMARY KEY,
    sector_id              INTEGER REFERENCES sectors(id) ON DELETE SET NULL,
    system_config_id       INTEGER NOT NULL REFERENCES system_configs(id),
    name                   TEXT NOT NULL,
    position_x_mpc          REAL,
    position_y_mpc          REAL,
    position_z_mpc          REAL,
    quadrant                TEXT CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')),
    is_binary              INTEGER NOT NULL DEFAULT 0 CHECK (is_binary IN (0, 1)),
    binary_separation_km            REAL,
    binary_type                     TEXT,
    binary_temperature_k            REAL,
    binary_radius_km                REAL,
    binary_effective_mass_kg        REAL,
    binary_effective_luminosity_w   REAL,
    binary_age_gy                   REAL,
    binary_lifespan_gy              REAL,
    binary_habitable_zone_inner_km  REAL,
    binary_habitable_zone_outer_km  REAL,
    binary_system_perimeter_km      REAL,
    binary_heliosphere_radius_km    REAL,
    binary_table_type        TEXT,
    binary_table_mass        TEXT,
    binary_table_lum         TEXT,
    binary_table_hab         TEXT,
    binary_table_separation  TEXT,
    binary_table_loc         TEXT,
    system_flavor_text   TEXT,
    schema_version       INTEGER NOT NULL DEFAULT 1,
    wikitext_content     TEXT,
    markdown_content     TEXT,
    mediawiki_url        TEXT,
    wikijs_url           TEXT,
    created_at           TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS stars (
    id                        INTEGER PRIMARY KEY,
    star_system_id            INTEGER NOT NULL REFERENCES star_systems(id) ON DELETE CASCADE,
    role                      TEXT NOT NULL CHECK (role IN ('primary', 'secondary', 'single')),
    name                      TEXT NOT NULL,
    star_type                 TEXT NOT NULL,
    yerkes_class              TEXT NOT NULL,
    mass_kg                   REAL NOT NULL,
    radius_km                 REAL NOT NULL,
    temperature_k             REAL NOT NULL,
    luminosity_w              REAL NOT NULL,
    age_gy                    REAL NOT NULL,
    lifespan_gy               REAL,
    habitable_zone_inner_km   REAL NOT NULL,
    habitable_zone_outer_km   REAL NOT NULL,
    system_perimeter_km       REAL NOT NULL,
    heliosphere_radius_km     REAL NOT NULL,
    table_type    TEXT NOT NULL,
    table_radius  TEXT NOT NULL,
    table_mass    TEXT NOT NULL,
    table_temp    TEXT NOT NULL,
    table_lum     TEXT NOT NULL,
    table_hab     TEXT NOT NULL,
    table_loc     TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS planets (
    id                        INTEGER PRIMARY KEY,
    star_system_id            INTEGER NOT NULL REFERENCES star_systems(id) ON DELETE CASCADE,
    star_id                   INTEGER REFERENCES stars(id) ON DELETE SET NULL,
    parent_planet_id          INTEGER REFERENCES planets(id) ON DELETE CASCADE,
    orbital_index             INTEGER NOT NULL,
    is_moon                   INTEGER NOT NULL DEFAULT 0 CHECK (is_moon IN (0, 1)),
    body_type                 TEXT NOT NULL CHECK (body_type IN ('t', 'g')),
    name                      TEXT NOT NULL,
    planet_class              TEXT,
    distance_km               REAL NOT NULL,
    radius_km                 REAL NOT NULL,
    mass_kg                   REAL NOT NULL,
    volume_km3                REAL NOT NULL,
    period_years              REAL NOT NULL,
    zone                      TEXT CHECK (zone IN ('h', 'e', 'c')),
    description               TEXT,
    gravity_g                 REAL,
    surface_temperature_k     REAL,
    density_g_cm3             REAL,
    atmosphere                TEXT,
    atm_density               REAL,
    atm_molar_density         REAL,
    atmospheric_pressure_pa   REAL,
    composition               TEXT,
    scale_height_km           REAL,
    hill_radius_km            REAL,
    min_orbit_distance_km     REAL,
    habitable_zone_inner_km   REAL NOT NULL,
    habitable_zone_outer_km   REAL NOT NULL,
    life_chemical             TEXT,
    evolutionary_speed        TEXT,
    flavor_text               TEXT,
    flavor_text_count         INTEGER NOT NULL DEFAULT 0,
    table_class     TEXT,
    table_distance  TEXT NOT NULL,
    table_period    TEXT NOT NULL,
    table_radius    TEXT NOT NULL,
    table_gravity   TEXT
);

CREATE TABLE IF NOT EXISTS planet_evolutionary_paragraphs (
    id          INTEGER PRIMARY KEY,
    planet_id   INTEGER NOT NULL REFERENCES planets(id) ON DELETE CASCADE,
    position    INTEGER NOT NULL,
    paragraph   TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS planet_reflection_spectrum (
    id             INTEGER PRIMARY KEY,
    planet_id      INTEGER NOT NULL REFERENCES planets(id) ON DELETE CASCADE,
    spectrum_type  TEXT NOT NULL CHECK (spectrum_type IN ('visible', 'non_visible')),
    position       INTEGER NOT NULL,
    value          TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS asteroid_belts (
    id                   INTEGER PRIMARY KEY,
    star_system_id       INTEGER NOT NULL REFERENCES star_systems(id) ON DELETE CASCADE,
    orbital_index        INTEGER NOT NULL,
    distance_km          REAL NOT NULL,
    lower_limit_km       REAL NOT NULL,
    upper_limit_km       REAL NOT NULL,
    density              TEXT NOT NULL CHECK (density IN ('dense', 'sparse', 'typical')),
    composition_summary  TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS asteroid_belt_composition (
    id             INTEGER PRIMARY KEY,
    belt_id        INTEGER NOT NULL REFERENCES asteroid_belts(id) ON DELETE CASCADE,
    position       INTEGER NOT NULL,
    component      TEXT NOT NULL,
    concentration  TEXT NOT NULL CHECK (concentration IN ('high', 'moderate', 'small', 'trace'))
);
"""


def _build_v1_database(path):
    """Creates a synthetic v1 database at `path`: one sector, one system
    config, one star system, one star, a top-level planet with its own
    evolutionary paragraph, and a moon of that planet with its own
    evolutionary paragraph *and* reflection-spectrum row (specifically
    exercising the "moon-owned child row" split, not just the
    planet-owned one) -- plus one asteroid belt with a composition row."""
    conn = sqlite3.connect(path)
    conn.executescript(V1_SCHEMA_SQL)

    conn.execute("INSERT INTO sectors (id, name, edge_mpc) VALUES (1, 'Test Sector', 1000.0)")
    conn.execute("INSERT INTO system_configs (id, markdown) VALUES (1, 0)")
    conn.execute(
        "INSERT INTO star_systems (id, sector_id, system_config_id, name) VALUES (1, 1, 1, 'Test System')"
    )
    conn.execute(
        """
        INSERT INTO stars (
            id, star_system_id, role, name, star_type, yerkes_class, mass_kg, radius_km,
            temperature_k, luminosity_w, age_gy, habitable_zone_inner_km, habitable_zone_outer_km,
            system_perimeter_km, heliosphere_radius_km, table_type, table_radius, table_mass,
            table_temp, table_lum, table_hab, table_loc
        ) VALUES (1, 1, 'single', 'Test Star', 'G2V Yellow Main Sequence Star', 'V', 1.0, 1.0,
                  5778, 1.0, 4.6, 1.0, 1.5, 2.0, 3.0, 'G2V', '1 R', '1 M', '5778 K', '1 L',
                  '1-1.5 AU', 'Test Star')
        """
    )
    # Top-level planet, id=10.
    conn.execute(
        """
        INSERT INTO planets (
            id, star_system_id, star_id, parent_planet_id, orbital_index, is_moon, body_type,
            name, planet_class, distance_km, radius_km, mass_kg, volume_km3, period_years,
            habitable_zone_inner_km, habitable_zone_outer_km, table_distance, table_period, table_radius
        ) VALUES (10, 1, 1, NULL, 0, 0, 't', 'Test Planet', 'M', 1.5e8, 6371, 5.97e24, 1.08e12, 1.0,
                  1.0, 1.5, '1 AU', '1 yr', '1 R')
        """
    )
    # Moon of the planet above, id=20.
    conn.execute(
        """
        INSERT INTO planets (
            id, star_system_id, star_id, parent_planet_id, orbital_index, is_moon, body_type,
            name, planet_class, distance_km, radius_km, mass_kg, volume_km3, period_years,
            habitable_zone_inner_km, habitable_zone_outer_km, table_distance, table_period, table_radius
        ) VALUES (20, 1, 1, 10, 0, 1, 't', 'Test Moon', 'D', 3.8e5, 1737, 7.3e22, 2.2e10, 0.08,
                  1.0, 1.5, '0.0025 AU', '0.08 yr', '0.27 R')
        """
    )
    conn.execute(
        "INSERT INTO planet_evolutionary_paragraphs (planet_id, position, paragraph) VALUES (10, 0, 'Planet paragraph.')"
    )
    conn.execute(
        "INSERT INTO planet_evolutionary_paragraphs (planet_id, position, paragraph) VALUES (20, 0, 'Moon paragraph.')"
    )
    conn.execute(
        "INSERT INTO planet_reflection_spectrum (planet_id, spectrum_type, position, value) "
        "VALUES (20, 'visible', 0, 'green')"
    )
    conn.execute(
        """
        INSERT INTO asteroid_belts (
            id, star_system_id, orbital_index, distance_km, lower_limit_km, upper_limit_km,
            density, composition_summary
        ) VALUES (1, 1, 1, 4.0e8, 3.5e8, 4.5e8, 'typical', 'trace amounts of iron')
        """
    )
    conn.execute(
        "INSERT INTO asteroid_belt_composition (belt_id, position, component, concentration) VALUES (1, 0, 'iron', 'trace')"
    )
    conn.commit()
    conn.close()


def test_migrate_v1_to_v2_splits_moons_into_their_own_table(tmp_path):
    db_path = str(tmp_path / "sector.db")
    _build_v1_database(db_path)

    backup_path = _db.migrate_database(db_path)

    assert backup_path is not None
    assert os.path.exists(backup_path)

    # The backup is untouched v1 data, not the migrated result.
    backup_conn = sqlite3.connect(backup_path)
    assert backup_conn.execute("PRAGMA user_version").fetchone()[0] == 1
    assert backup_conn.execute("SELECT COUNT(*) FROM planets").fetchone()[0] == 2
    backup_conn.close()

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        assert conn.execute("PRAGMA user_version").fetchone()[0] == _db.SCHEMA_VERSION

        planets = conn.execute("SELECT * FROM planets").fetchall()
        assert len(planets) == 1
        assert planets[0]["id"] == 10
        assert planets[0]["name"] == "Test Planet"

        moons = conn.execute("SELECT * FROM moons").fetchall()
        assert len(moons) == 1
        assert moons[0]["id"] == 20
        assert moons[0]["planet_id"] == 10
        assert moons[0]["name"] == "Test Moon"
        assert moons[0]["star_system_id"] == 1
        assert moons[0]["star_id"] == 1

        planet_paragraphs = conn.execute("SELECT * FROM planet_evolutionary_paragraphs").fetchall()
        assert len(planet_paragraphs) == 1
        assert planet_paragraphs[0]["planet_id"] == 10
        assert planet_paragraphs[0]["paragraph"] == "Planet paragraph."

        moon_paragraphs = conn.execute("SELECT * FROM moon_evolutionary_paragraphs").fetchall()
        assert len(moon_paragraphs) == 1
        assert moon_paragraphs[0]["moon_id"] == 20
        assert moon_paragraphs[0]["paragraph"] == "Moon paragraph."

        assert conn.execute("SELECT COUNT(*) FROM planet_reflection_spectrum").fetchone()[0] == 0
        moon_spectrum = conn.execute("SELECT * FROM moon_reflection_spectrum").fetchall()
        assert len(moon_spectrum) == 1
        assert moon_spectrum[0]["moon_id"] == 20
        assert moon_spectrum[0]["value"] == "green"

        # Unrelated tables copied through untouched.
        assert conn.execute("SELECT COUNT(*) FROM sectors").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM asteroid_belts").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM asteroid_belt_composition").fetchone()[0] == 1

        object_types = {row[0] for row in conn.execute("SELECT DISTINCT object_type FROM sector_objects")}
        assert object_types == {"star", "planet", "moon", "asteroid_belt"}
    finally:
        conn.close()


def test_migrate_database_is_noop_when_already_current(tmp_path):
    db_path = str(tmp_path / "current.db")
    conn = _db.get_connection(db_path)
    conn.close()

    result = _db.migrate_database(db_path)

    assert result is None
    # No backup file should have been created for a no-op.
    assert sorted(os.listdir(tmp_path)) == ["current.db"]


def test_migrate_database_rejects_unknown_schema_version(tmp_path):
    db_path = str(tmp_path / "future.db")
    conn = sqlite3.connect(db_path)
    conn.execute(f"PRAGMA user_version = {_db.SCHEMA_VERSION + 1}")
    conn.close()

    with pytest.raises(_db.UnsupportedSchemaVersionError):
        _db.migrate_database(db_path)
