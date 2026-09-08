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

import gzip
import os
import shutil
import sqlite3
import sys

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

    # The backup is gzip-compressed, not a plain SQLite file copy -- and
    # its name carries `BACKUP_MARKER` and doesn't end in bare `.db`, so
    # it can never be mistaken for a live database by a naive `*.db` glob
    # (see `html/lib/dbutil.py` and `migrateDb.py`).
    assert backup_path.endswith(".db.gz")
    assert _db.BACKUP_MARKER in os.path.basename(backup_path)
    assert not backup_path.endswith(".db")
    with open(backup_path, "rb") as f:
        assert f.read(2) == b"\x1f\x8b"  # gzip magic bytes

    # The backup is untouched v1 data, not the migrated result.
    decompressed_path = str(tmp_path / "decompressed.db")
    with gzip.open(backup_path, "rb") as f_in, open(decompressed_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    backup_conn = sqlite3.connect(decompressed_path)
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


# --- html/lib/dbutil.py: the database picker must never surface a backup ---
#
# `src/html/lib` isn't part of the installed `stellarObjects` package
# (CGI-only plumbing, see `dbutil.py`'s own module docstring), so it's
# added to `sys.path` here the same way the CGI scripts themselves do.
# This file lives at src/tests/, two levels under the repo root (src
# layout), not one -- two dirname() calls reach `src/`, then down into
# html/lib.
_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_HTML_LIB_DIR = os.path.join(_SRC_DIR, "html", "lib")
if _HTML_LIB_DIR not in sys.path:
    sys.path.insert(0, _HTML_LIB_DIR)

import dbutil  # noqa: E402


def test_list_databases_excludes_migration_backups(tmp_path, monkeypatch):
    monkeypatch.setattr(dbutil, "DEFAULT_DB_DIR", str(tmp_path))
    monkeypatch.delenv(dbutil.DB_DIR_ENV_VAR, raising=False)

    real_db_path = tmp_path / "planetgen.db"
    _db.get_connection(str(real_db_path)).close()

    old_db_path = tmp_path / "old.db"
    old_db_path.write_bytes(b"not a real db, just needs to exist")
    # Build a stand-in backup the same way `migrate_database` names one,
    # without needing a real old-schema database for this test's purposes.
    backup_path = str(tmp_path / f"old.db.v2{_db.BACKUP_MARKER}20260101T000000.db.gz")
    with open(backup_path, "wb") as f:
        f.write(b"\x1f\x8b" + b"0" * 10)  # gzip magic + filler, contents unused here

    listed_names = {entry["name"] for entry in dbutil.list_databases()}

    assert listed_names == {"planetgen.db", "old.db"}
    assert os.path.basename(backup_path) not in listed_names


def test_resolve_db_path_refuses_a_backup_even_by_exact_name(tmp_path, monkeypatch):
    monkeypatch.setattr(dbutil, "DEFAULT_DB_DIR", str(tmp_path))
    monkeypatch.delenv(dbutil.DB_DIR_ENV_VAR, raising=False)

    backup_name = f"planetgen.db.v2{_db.BACKUP_MARKER}20260101T000000.db.gz"
    (tmp_path / backup_name).write_bytes(b"\x1f\x8b" + b"0" * 10)

    with pytest.raises(dbutil.NotFoundError):
        dbutil.resolve_db_path(backup_name)


# --- migrateDb.py: re-running the CLI must not re-migrate its own backup ---

import migrateDb  # noqa: E402


def test_running_migrate_cli_twice_does_not_re_migrate_the_backup(tmp_path):
    db_path = str(tmp_path / "sector.db")
    _build_v1_database(db_path)

    old_argv = sys.argv
    try:
        sys.argv = ["migrateDb.py", str(tmp_path)]
        migrateDb.main()
        first_run_files = sorted(os.listdir(tmp_path))
        backups_after_first_run = [f for f in first_run_files if _db.BACKUP_MARKER in f]
        assert len(backups_after_first_run) == 1

        backup_path = tmp_path / backups_after_first_run[0]
        mtime_after_first_run = backup_path.stat().st_mtime
        size_after_first_run = backup_path.stat().st_size

        # Second run: the live db is already current (a no-op for it), and
        # the backup from the first run must be skipped entirely rather
        # than handed back into `migrate_database` (which would treat its
        # still-v1 `PRAGMA user_version` as needing migration all over
        # again, producing a nested second backup).
        migrateDb.main()
    finally:
        sys.argv = old_argv

    second_run_files = sorted(os.listdir(tmp_path))
    backups_after_second_run = [f for f in second_run_files if _db.BACKUP_MARKER in f]

    assert len(backups_after_second_run) == 1
    assert backups_after_second_run == backups_after_first_run
    assert backup_path.stat().st_mtime == mtime_after_first_run
    assert backup_path.stat().st_size == size_after_first_run


# ---------------------------------------------------------------------------
# v3 -> v4 (Track C: galaxy coordinate system) -- sectors gains six nullable
# galaxy-placement columns (center_x/y/z_pc, galactic_radius_pc,
# shell_index, shell_slot_index). V3_SCHEMA_SQL is the schema exactly as it
# stood immediately before that change: identical to the current
# schema.sql in every other table (v2's moons split and v3's
# star_systems.location are both already present), with `sectors` still
# the plain 3-column table.
# ---------------------------------------------------------------------------

V3_SCHEMA_SQL = """
PRAGMA foreign_keys = ON;
PRAGMA user_version = 3;

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
    location                TEXT,
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
    orbital_index             INTEGER NOT NULL,
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

CREATE TABLE IF NOT EXISTS moons (
    id                        INTEGER PRIMARY KEY,
    planet_id                 INTEGER NOT NULL REFERENCES planets(id) ON DELETE CASCADE,
    star_system_id            INTEGER NOT NULL REFERENCES star_systems(id) ON DELETE CASCADE,
    star_id                   INTEGER REFERENCES stars(id) ON DELETE SET NULL,
    orbital_index             INTEGER NOT NULL,
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

CREATE TABLE IF NOT EXISTS moon_evolutionary_paragraphs (
    id          INTEGER PRIMARY KEY,
    moon_id     INTEGER NOT NULL REFERENCES moons(id) ON DELETE CASCADE,
    position    INTEGER NOT NULL,
    paragraph   TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS moon_reflection_spectrum (
    id             INTEGER PRIMARY KEY,
    moon_id        INTEGER NOT NULL REFERENCES moons(id) ON DELETE CASCADE,
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


def _build_v3_database(path):
    """Creates a synthetic v3 database at `path`: one sector (with a
    system placed in it, so `location`/`quadrant`/position are exercised),
    one star, one top-level planet with a moon, and an asteroid belt --
    the plain 3-column `sectors` table this migration must extend."""
    conn = sqlite3.connect(path)
    conn.executescript(V3_SCHEMA_SQL)

    conn.execute("INSERT INTO sectors (id, name, edge_mpc) VALUES (1, 'Legacy Sector', 3526.0)")
    conn.execute("INSERT INTO system_configs (id, markdown) VALUES (1, 0)")
    conn.execute(
        """
        INSERT INTO star_systems (
            id, sector_id, system_config_id, name,
            position_x_mpc, position_y_mpc, position_z_mpc, quadrant, location
        ) VALUES (1, 1, 1, 'Legacy System', 100.0, 200.0, 300.0, 'I', 'Legacy Sector')
        """
    )
    conn.execute(
        """
        INSERT INTO stars (
            id, star_system_id, role, name, star_type, yerkes_class, mass_kg, radius_km,
            temperature_k, luminosity_w, age_gy, habitable_zone_inner_km, habitable_zone_outer_km,
            system_perimeter_km, heliosphere_radius_km, table_type, table_radius, table_mass,
            table_temp, table_lum, table_hab, table_loc
        ) VALUES (1, 1, 'single', 'Legacy Star', 'G2V Yellow Main Sequence Star', 'V', 1.0, 1.0,
                  5778, 1.0, 4.6, 1.0, 1.5, 2.0, 3.0, 'G2V', '1 R', '1 M', '5778 K', '1 L',
                  '1-1.5 AU', 'Legacy Star')
        """
    )
    conn.execute(
        """
        INSERT INTO planets (
            id, star_system_id, star_id, orbital_index, body_type,
            name, planet_class, distance_km, radius_km, mass_kg, volume_km3, period_years,
            habitable_zone_inner_km, habitable_zone_outer_km, table_distance, table_period, table_radius
        ) VALUES (10, 1, 1, 0, 't', 'Legacy Planet', 'M', 1.5e8, 6371, 5.97e24, 1.08e12, 1.0,
                  1.0, 1.5, '1 AU', '1 yr', '1 R')
        """
    )
    conn.execute(
        """
        INSERT INTO moons (
            id, planet_id, star_system_id, star_id, orbital_index, body_type,
            name, planet_class, distance_km, radius_km, mass_kg, volume_km3, period_years,
            habitable_zone_inner_km, habitable_zone_outer_km, table_distance, table_period, table_radius
        ) VALUES (20, 10, 1, 1, 0, 't', 'Legacy Moon', 'D', 3.8e5, 1737, 7.3e22, 2.2e10, 0.08,
                  1.0, 1.5, '0.0025 AU', '0.08 yr', '0.27 R')
        """
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


def test_migrate_v3_to_v4_adds_null_galaxy_columns_to_existing_sectors(tmp_path):
    db_path = str(tmp_path / "sector.db")
    _build_v3_database(db_path)

    backup_path = _db.migrate_database(db_path)

    assert backup_path is not None
    assert os.path.exists(backup_path)
    assert backup_path.endswith(".db.gz")
    assert _db.BACKUP_MARKER in os.path.basename(backup_path)

    # Same gzip-compressed backup convention as the v1 -> v2 test above
    # (`test_migrate_v1_to_v2_splits_moons_into_their_own_table`) --
    # `migrate_database` always writes a gzip-compressed backup regardless
    # of which source version it's migrating from, so this can't
    # `sqlite3.connect` the backup path directly.
    decompressed_backup_path = str(tmp_path / "decompressed_v3_backup.db")
    with gzip.open(backup_path, "rb") as f_in, open(decompressed_backup_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    backup_conn = sqlite3.connect(decompressed_backup_path)
    assert backup_conn.execute("PRAGMA user_version").fetchone()[0] == 3
    backup_conn.close()

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        assert conn.execute("PRAGMA user_version").fetchone()[0] == _db.SCHEMA_VERSION

        sectors = conn.execute("SELECT * FROM sectors").fetchall()
        assert len(sectors) == 1
        sector = sectors[0]
        assert sector["id"] == 1
        assert sector["name"] == "Legacy Sector"
        assert sector["edge_mpc"] == pytest.approx(3526.0)
        # The whole point of this migration: the new galaxy-placement
        # columns exist and are NULL for a pre-v4 sector -- there is no
        # way to recover an intended galaxy position after the fact.
        assert sector["center_x_pc"] is None
        assert sector["center_y_pc"] is None
        assert sector["center_z_pc"] is None
        assert sector["galactic_radius_pc"] is None
        assert sector["shell_index"] is None
        assert sector["shell_slot_index"] is None

        # Every other table copied through untouched.
        system_row = conn.execute("SELECT * FROM star_systems WHERE id = 1").fetchone()
        assert system_row["name"] == "Legacy System"
        assert system_row["location"] == "Legacy Sector"
        assert system_row["quadrant"] == "I"

        assert conn.execute("SELECT COUNT(*) FROM stars").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM planets").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM moons").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM asteroid_belts").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM asteroid_belt_composition").fetchone()[0] == 1

        violations = conn.execute("PRAGMA foreign_key_check").fetchall()
        assert violations == []
    finally:
        conn.close()


def test_migrate_v3_to_v4_new_sector_saved_afterward_can_have_a_galaxy_position(tmp_path):
    """After migrating an old database, saving a brand-new, galaxy-placed
    sector into it must work -- confirms the migrated `sectors` table
    really does have the new columns (and their CHECK constraint), not
    just that the old row's columns read back as NULL."""
    db_path = str(tmp_path / "sector.db")
    _build_v3_database(db_path)
    _db.migrate_database(db_path)

    from stellarObjects.spaceSector import SpaceSector

    sector = SpaceSector("New Placed Sector", edge_ly=11.5)
    galaxy_position = {
        "center_x_pc": 1.0, "center_y_pc": 2.0, "center_z_pc": 3.0,
        "galactic_radius_pc": (1.0 ** 2 + 2.0 ** 2 + 3.0 ** 2) ** 0.5,
        "shell_index": 0, "shell_slot_index": 0,
    }
    sector_id = _db.save_sector(sector, db_path=db_path, galaxy_position=galaxy_position)

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        row = conn.execute("SELECT * FROM sectors WHERE id = ?", (sector_id,)).fetchone()
        assert row["center_x_pc"] == pytest.approx(1.0)
        assert row["center_y_pc"] == pytest.approx(2.0)
        assert row["center_z_pc"] == pytest.approx(3.0)
        assert row["galactic_radius_pc"] == pytest.approx((14.0) ** 0.5)
        assert row["shell_index"] == 0
        assert row["shell_slot_index"] == 0
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# v4 -> v5 (pre-rendered display text removed) -- every `table_*`/
# `binary_table_*` TEXT column dropped from `star_systems`/`stars`/
# `planets`/`moons` (see schema.sql's "v5" header note): these held a
# formatted display string (e.g. wikitext "{{Exp|4.20|5}}" or HTML
# "4.20 &times; 10<sup>5</sup>") baked in at generation time, not data --
# every number it was derived from already has its own proper column
# (mass_kg, radius_km, ...), so the migration just drops them.
# V4_SCHEMA_SQL is `V3_SCHEMA_SQL` plus v4's six nullable galaxy-placement
# `sectors` columns (identical to today's `schema.sql` in every table).
# ---------------------------------------------------------------------------

V4_SCHEMA_SQL = V3_SCHEMA_SQL.replace(
    "PRAGMA user_version = 3;",
    "PRAGMA user_version = 4;",
).replace(
    """CREATE TABLE IF NOT EXISTS sectors (
    id        INTEGER PRIMARY KEY,
    name      TEXT NOT NULL,
    edge_mpc  REAL NOT NULL
);""",
    """CREATE TABLE IF NOT EXISTS sectors (
    id                  INTEGER PRIMARY KEY,
    name                TEXT NOT NULL,
    edge_mpc            REAL NOT NULL,
    center_x_pc         REAL,
    center_y_pc         REAL,
    center_z_pc         REAL,
    galactic_radius_pc  REAL,
    shell_index         INTEGER,
    shell_slot_index    INTEGER,
    CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    )
);""",
)
assert "PRAGMA user_version = 4;" in V4_SCHEMA_SQL
assert "galactic_radius_pc" in V4_SCHEMA_SQL


def _build_v4_database(path):
    """Creates a synthetic v4 database at `path`: one sector, one binary
    star system (so every `binary_table_*` column is exercised, not just
    the per-star `table_*` ones) with two stars, a top-level planet with a
    moon, and an asteroid belt -- every `table_*`/`binary_table_*` column
    populated with the wikitext-template display strings a real pre-v5
    deployment would have stored."""
    conn = sqlite3.connect(path)
    conn.executescript(V4_SCHEMA_SQL)

    conn.execute("INSERT INTO sectors (id, name, edge_mpc) VALUES (1, 'Test Sector', 1000.0)")
    conn.execute("INSERT INTO system_configs (id, markdown) VALUES (1, 0)")
    conn.execute(
        """
        INSERT INTO star_systems (
            id, sector_id, system_config_id, name, is_binary,
            binary_separation_km, binary_type, binary_temperature_k, binary_radius_km,
            binary_effective_mass_kg, binary_effective_luminosity_w, binary_age_gy,
            binary_habitable_zone_inner_km, binary_habitable_zone_outer_km,
            binary_system_perimeter_km, binary_heliosphere_radius_km,
            binary_table_type, binary_table_mass, binary_table_lum, binary_table_hab,
            binary_table_separation, binary_table_loc
        ) VALUES (
            1, 1, 1, 'Test Binary System', 1,
            1.5e8, 'G2V + K5V Binary', 5000, 1.2,
            2.5e30, 1.5e26, 4.6,
            1.0, 1.5, 3.0, 4.0,
            'G2V + K5V', '{{Exp|2.5|30}} kg (100.00% of Sol)', '{{Exp|1.5|26}} W (100.00% of Sol)',
            'Between 1.0 and 1.5 AU', '{{Exp|1.5|8}} km (1.00 AU)', 'Test Binary System'
        )
        """
    )
    conn.execute(
        """
        INSERT INTO stars (
            id, star_system_id, role, name, star_type, yerkes_class, mass_kg, radius_km,
            temperature_k, luminosity_w, age_gy, habitable_zone_inner_km, habitable_zone_outer_km,
            system_perimeter_km, heliosphere_radius_km, table_type, table_radius, table_mass,
            table_temp, table_lum, table_hab, table_loc
        ) VALUES (1, 1, 'primary', 'Test Star A', 'G2V Yellow Main Sequence Star', 'V', 2.0e30, 700000,
                  5778, 3.8e26, 4.6, 1.0, 1.5, 2.0, 3.0, 'G2V', '{{Exp|7.0|5}} km', '{{Exp|2.0|30}} kg (100.00% of Sol)',
                  '5778 K', '{{Exp|3.8|26}} W (100.00% of Sol)', '1-1.5 AU', 'Test Star A')
        """
    )
    conn.execute(
        """
        INSERT INTO stars (
            id, star_system_id, role, name, star_type, yerkes_class, mass_kg, radius_km,
            temperature_k, luminosity_w, age_gy, habitable_zone_inner_km, habitable_zone_outer_km,
            system_perimeter_km, heliosphere_radius_km, table_type, table_radius, table_mass,
            table_temp, table_lum, table_hab, table_loc
        ) VALUES (2, 1, 'secondary', 'Test Star B', 'K5V Orange Main Sequence Star', 'V', 5.0e29, 500000,
                  4500, 1.2e26, 4.6, 1.0, 1.5, 2.0, 3.0, 'K5V', '{{Exp|5.0|5}} km', '{{Exp|5.0|29}} kg (25.00% of Sol)',
                  '4500 K', '{{Exp|1.2|26}} W (31.58% of Sol)', '1-1.5 AU', 'Test Star B')
        """
    )
    conn.execute(
        """
        INSERT INTO planets (
            id, star_system_id, star_id, orbital_index, body_type,
            name, planet_class, distance_km, radius_km, mass_kg, volume_km3, period_years, gravity_g,
            habitable_zone_inner_km, habitable_zone_outer_km,
            table_class, table_distance, table_period, table_radius, table_gravity
        ) VALUES (10, 1, NULL, 0, 't', 'Test Planet', 'M', 1.5e8, 6371, 5.97e24, 1.08e12, 1.0, 1.0,
                  1.0, 1.5, 'M', '1.000 AU', '1 year', '6371.00 km', '1.0 g')
        """
    )
    conn.execute(
        """
        INSERT INTO moons (
            id, planet_id, star_system_id, star_id, orbital_index, body_type,
            name, planet_class, distance_km, radius_km, mass_kg, volume_km3, period_years, gravity_g,
            habitable_zone_inner_km, habitable_zone_outer_km,
            table_class, table_distance, table_period, table_radius, table_gravity
        ) VALUES (20, 10, 1, NULL, 0, 't', 'Test Moon', 'D', 3.8e5, 1737, 7.3e22, 2.2e10, 0.08, 0.17,
                  1.0, 1.5, 'D', '{{Exp|3.8|5}} km', '29 days', '1737.00 km', '0.17 g')
        """
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


def test_migrate_v4_to_v5_drops_table_columns_and_preserves_raw_numbers(tmp_path):
    db_path = str(tmp_path / "sector.db")
    _build_v4_database(db_path)

    backup_path = _db.migrate_database(db_path)

    assert backup_path is not None
    assert os.path.exists(backup_path)

    decompressed_backup_path = str(tmp_path / "decompressed_v4_backup.db")
    with gzip.open(backup_path, "rb") as f_in, open(decompressed_backup_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    backup_conn = sqlite3.connect(decompressed_backup_path)
    assert backup_conn.execute("PRAGMA user_version").fetchone()[0] == 4
    # The backup is untouched v4 data, still carrying the old display strings.
    assert backup_conn.execute("SELECT table_mass FROM stars WHERE id = 1").fetchone()[0] == "{{Exp|2.0|30}} kg (100.00% of Sol)"
    backup_conn.close()

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        assert conn.execute("PRAGMA user_version").fetchone()[0] == _db.SCHEMA_VERSION == 5

        # The columns are genuinely gone, not just left NULL/unused.
        for table, column in (
            ("stars", "table_mass"), ("stars", "table_type"),
            ("planets", "table_class"), ("moons", "table_distance"),
            ("star_systems", "binary_table_mass"),
        ):
            with pytest.raises(sqlite3.OperationalError):
                conn.execute(f"SELECT {column} FROM {table}").fetchall()

        # The raw numbers these display strings were derived from survive
        # the migration untouched.
        star_a = conn.execute("SELECT * FROM stars WHERE id = 1").fetchone()
        assert star_a["mass_kg"] == pytest.approx(2.0e30)
        assert star_a["radius_km"] == pytest.approx(700000)
        assert star_a["temperature_k"] == pytest.approx(5778)

        system_row = conn.execute("SELECT * FROM star_systems WHERE id = 1").fetchone()
        assert system_row["binary_effective_mass_kg"] == pytest.approx(2.5e30)
        assert system_row["binary_separation_km"] == pytest.approx(1.5e8)

        planet = conn.execute("SELECT * FROM planets WHERE id = 10").fetchone()
        assert planet["distance_km"] == pytest.approx(1.5e8)
        assert planet["gravity_g"] == pytest.approx(1.0)

        moon = conn.execute("SELECT * FROM moons WHERE id = 20").fetchone()
        assert moon["distance_km"] == pytest.approx(3.8e5)

        assert conn.execute("SELECT COUNT(*) FROM asteroid_belts").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM asteroid_belt_composition").fetchone()[0] == 1

        violations = conn.execute("PRAGMA foreign_key_check").fetchall()
        assert violations == []
    finally:
        conn.close()
