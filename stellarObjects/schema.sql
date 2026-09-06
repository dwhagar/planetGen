-- stellarObjects/schema.sql
--
-- SQLite schema for planetGen persistence (TODO.md Phase 2). Plain SQL DDL,
-- no ORM. Column lists are verified against the actual generator source
-- (stellarObjects/config.py, starData.py, doubleStar.py, planetData.py,
-- asteroidData.py, spaceSector.py) as of schema version 1 -- not against
-- TODO.md's earlier field-count estimates.
--
-- Searchable-field principle: the primary fields a caller will search on
-- are exactly the ones each object already exposes in its own generated
-- text -- the `*_properties` dict built in `to_paragraph_list()` for
-- stars/binaries/planets (`table_*` columns below), and the equivalent
-- prose facts for asteroid belts, which have no properties dict at all
-- (density, the distance range, and composition -- see the asteroid_belts
-- section below). Raw generative scalar columns (mass_kg, radius_km, etc.)
-- are kept alongside for full object-graph fidelity (Phase 1 merge
-- decision), but the `table_*`/prose-derived columns are the ones meant
-- for search, since they're exactly what's published.
--
-- Two independent version numbers, per TODO.md's "Schema versioning"
-- decision:
--   - PRAGMA user_version below is the DDL-level structure version.
--   - star_systems.schema_version (per row) is the version of the
--     serialized object-graph shape (Phase 1's to_dict()) that produced it.
-- Both start at 1.
PRAGMA foreign_keys = ON;
PRAGMA user_version = 1;

-- ---------------------------------------------------------------------
-- sectors
-- ---------------------------------------------------------------------
CREATE TABLE sectors (
    id       INTEGER PRIMARY KEY,
    name     TEXT NOT NULL,
    edge_ly  REAL NOT NULL
);

-- ---------------------------------------------------------------------
-- system_configs -- one row per SystemConfig "recipe"
-- (stellarObjects/config.py:29-33 SERIALIZABLE_FIELDS)
-- ---------------------------------------------------------------------
CREATE TABLE system_configs (
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

-- Child table for the variable-length SLOTS recipe list
-- (config.py:135-151).
CREATE TABLE system_config_slots (
    id            INTEGER PRIMARY KEY,
    config_id     INTEGER NOT NULL REFERENCES system_configs(id) ON DELETE CASCADE,
    orbit_index   INTEGER NOT NULL,
    type          TEXT CHECK (type IN ('planet', 'asteroid_belt')),
    planet_class  TEXT,
    moons         INTEGER
);
CREATE INDEX idx_system_config_slots_config_id ON system_config_slots(config_id);

-- ---------------------------------------------------------------------
-- star_systems
-- ---------------------------------------------------------------------
CREATE TABLE star_systems (
    id                     INTEGER PRIMARY KEY,
    sector_id              INTEGER REFERENCES sectors(id) ON DELETE SET NULL,
    system_config_id       INTEGER NOT NULL REFERENCES system_configs(id),
    name                   TEXT NOT NULL,

    -- Position within its sector (spaceSector.py SectorSystemEntry.position),
    -- light-years, relative to the sector's cubic center. NULL iff this
    -- system was never placed in a sector.
    position_x_ly          REAL,
    position_y_ly          REAL,
    position_z_ly          REAL,

    is_binary              INTEGER NOT NULL DEFAULT 0 CHECK (is_binary IN (0, 1)),

    -- BinaryStarProxy-derived fields (doubleStar.py) -- all NULL for a
    -- single-star system, stored rather than re-derived since
    -- _effective_mass/_effective_luminosity are computed once at
    -- generation time.
    binary_separation_au           REAL,
    binary_type                    TEXT,
    binary_temperature_k           REAL,
    binary_radius_km               REAL,
    binary_effective_mass_kg       REAL,
    binary_effective_luminosity_w  REAL,
    binary_age_gy                  REAL,
    binary_lifespan_gy             REAL,  -- NULL = float('inf')
    binary_habitable_zone_inner_au REAL,
    binary_habitable_zone_outer_au REAL,
    binary_system_perimeter_au     REAL,
    binary_heliosphere_radius_au   REAL,

    -- "Binary System Data" table (doubleStar.py:158-170), one column per
    -- key -- the one properties table with no owning row elsewhere, since
    -- BinaryStarProxy is never itself stored as a `stars` row. NULL unless
    -- is_binary.
    binary_table_type        TEXT,
    binary_table_mass        TEXT,
    binary_table_lum         TEXT,
    binary_table_hab         TEXT,
    binary_table_separation  TEXT,
    binary_table_loc         TEXT,

    system_flavor_text   TEXT,
    schema_version       INTEGER NOT NULL DEFAULT 1,

    -- Fully rendered page text, ready to paste/upload. Both must be
    -- rendered from the SAME already-generated StarSystem object
    -- (toggle system_config.MARKDOWN, render, toggle back, render again)
    -- -- generation mixes unseedable `secrets` with seedable `random`
    -- (see spaceSector.py's module docstring), so the other format can
    -- never be faithfully regenerated later.
    wikitext_content     TEXT,
    markdown_content     TEXT,

    -- One system = one wiki page (stars/planets/moons are sections within
    -- it, per StarSystem.__str__), so exactly one URL per wiki target.
    mediawiki_url        TEXT,
    wikijs_url           TEXT,

    created_at           TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX idx_star_systems_sector_id ON star_systems(sector_id);
CREATE INDEX idx_star_systems_system_config_id ON star_systems(system_config_id);

-- ---------------------------------------------------------------------
-- stars -- one row per single star, two rows (primary/secondary) per
-- binary. Never a row for the BinaryStarProxy itself (its snapshot lives
-- on star_systems above). Fields verified against starData.py:399-445 /
-- Star.__init__.
-- ---------------------------------------------------------------------
CREATE TABLE stars (
    id                        INTEGER PRIMARY KEY,
    star_system_id            INTEGER NOT NULL REFERENCES star_systems(id) ON DELETE CASCADE,
    role                      TEXT NOT NULL CHECK (role IN ('primary', 'secondary', 'single')),
    name                      TEXT NOT NULL,
    star_type                 TEXT NOT NULL,   -- Star.type, e.g. "G2V Yellow Main Sequence Star"
    yerkes_class              TEXT NOT NULL,
    mass_kg                   REAL NOT NULL,
    radius_km                 REAL NOT NULL,
    temperature_k             REAL NOT NULL,
    luminosity_w              REAL NOT NULL,
    age_gy                    REAL NOT NULL,
    lifespan_gy               REAL,            -- NULL = float('inf'), white dwarfs
    habitable_zone_inner_au   REAL NOT NULL,
    habitable_zone_outer_au   REAL NOT NULL,
    system_perimeter_au       REAL NOT NULL,
    heliosphere_radius_au     REAL NOT NULL,

    -- "Star Data" table (starData.py:488-509), one column per key --
    -- always present, every constituent Star renders its own individual
    -- table even inside a binary.
    table_type    TEXT NOT NULL,
    table_radius  TEXT NOT NULL,
    table_mass    TEXT NOT NULL,
    table_temp    TEXT NOT NULL,
    table_lum     TEXT NOT NULL,
    table_hab     TEXT NOT NULL,
    table_loc     TEXT NOT NULL
);
CREATE INDEX idx_stars_star_system_id ON stars(star_system_id);

-- ---------------------------------------------------------------------
-- planets -- covers moons via self-reference. Fields verified against
-- planetData.py:117-200 / Planet.__init__ (27 scalar fields; excludes
-- system_config/star back-refs, the `moons` list [this self-reference],
-- `evolutionary_data` [-> planet_evolutionary_paragraphs], and
-- reflection_spectrum_visible/non_visible [-> planet_reflection_spectrum]).
-- ---------------------------------------------------------------------
CREATE TABLE planets (
    id                        INTEGER PRIMARY KEY,
    star_system_id            INTEGER NOT NULL REFERENCES star_systems(id) ON DELETE CASCADE,
    -- The specific star this planet orbits, when that's a real stored
    -- `stars` row -- true for every single-star system. NULL for a
    -- binary's planets: systemData.py always generates planets against
    -- `self.star`, which for binaries is the `BinaryStarProxy` (never one
    -- individual constituent star -- there is no S-type/circumbinary
    -- choice in the current generator), and the proxy is deliberately not
    -- stored as its own `stars` row (see star_systems.binary_table_*
    -- above). `star_system_id` above is always the reliable owning link
    -- regardless of star_id.
    star_id                   INTEGER REFERENCES stars(id) ON DELETE SET NULL,
    parent_planet_id          INTEGER REFERENCES planets(id) ON DELETE CASCADE,  -- set for moons
    orbital_index             INTEGER NOT NULL,
    is_moon                   INTEGER NOT NULL DEFAULT 0 CHECK (is_moon IN (0, 1)),
    body_type                 TEXT NOT NULL CHECK (body_type IN ('t', 'g')),
    name                      TEXT NOT NULL,
    planet_class              TEXT,
    distance_au               REAL NOT NULL,   -- from star (top-level) or parent planet (moon)
    radius_km                 REAL NOT NULL,
    mass_kg                   REAL NOT NULL,
    volume_km3                REAL NOT NULL,   -- derived, stored not recomputed
    period_years              REAL NOT NULL,   -- derived, stored not recomputed
    zone                      TEXT CHECK (zone IN ('h', 'e', 'c')),
    description               TEXT,
    gravity_g                 REAL,
    surface_temperature_k     REAL,
    density_g_cm3             REAL,
    atmosphere                TEXT,
    atm_density               REAL,
    atm_molar_density         REAL,
    atmospheric_pressure_pa   REAL,
    composition               TEXT,            -- descriptive string (contrast asteroid_belt_composition)
    scale_height_km           REAL,
    hill_radius_km            REAL,
    min_orbit_distance_au     REAL,
    habitable_zone_inner_au   REAL NOT NULL,   -- copied from host star at generation time
    habitable_zone_outer_au   REAL NOT NULL,
    life_chemical             TEXT,
    evolutionary_speed        TEXT,
    flavor_text               TEXT,
    flavor_text_count         INTEGER NOT NULL DEFAULT 0,

    -- "Planet Data" / "Class Data" table (planetData.py:291-302), one
    -- column per key -- same dict shape for planets and moons.
    table_class     TEXT,
    table_distance  TEXT NOT NULL,
    table_period    TEXT NOT NULL,
    table_radius    TEXT NOT NULL,
    table_gravity   TEXT
);
CREATE INDEX idx_planets_star_system_id ON planets(star_system_id);
CREATE INDEX idx_planets_star_id ON planets(star_id);
CREATE INDEX idx_planets_parent_planet_id ON planets(parent_planet_id);

-- Child table for the variable-length evolutionary narrative list
-- (Planet.evolutionary_data).
CREATE TABLE planet_evolutionary_paragraphs (
    id          INTEGER PRIMARY KEY,
    planet_id   INTEGER NOT NULL REFERENCES planets(id) ON DELETE CASCADE,
    position    INTEGER NOT NULL,
    paragraph   TEXT NOT NULL
);
CREATE INDEX idx_planet_evolutionary_paragraphs_planet_id ON planet_evolutionary_paragraphs(planet_id);

-- Child table replacing reflection_spectrum_visible/non_visible -- no
-- JSON column, same position+value shape as the paragraphs/composition
-- child tables (this schema is a pure relational design, not a
-- JSON-blob-in-a-column hybrid).
CREATE TABLE planet_reflection_spectrum (
    id             INTEGER PRIMARY KEY,
    planet_id      INTEGER NOT NULL REFERENCES planets(id) ON DELETE CASCADE,
    spectrum_type  TEXT NOT NULL CHECK (spectrum_type IN ('visible', 'non_visible')),
    position       INTEGER NOT NULL,
    value          TEXT NOT NULL
);
CREATE INDEX idx_planet_reflection_spectrum_planet_id ON planet_reflection_spectrum(planet_id);

-- ---------------------------------------------------------------------
-- asteroid_belts -- belts have no properties-dict data table
-- (asteroidData.py:93-141 is prose only), so their searchable columns are
-- the facts that prose always states instead: density, the distance range
-- from the star, and composition.
-- ---------------------------------------------------------------------
CREATE TABLE asteroid_belts (
    id                   INTEGER PRIMARY KEY,
    star_system_id       INTEGER NOT NULL REFERENCES star_systems(id) ON DELETE CASCADE,
    orbital_index        INTEGER NOT NULL,
    distance_au          REAL NOT NULL,
    lower_limit_au       REAL NOT NULL,   -- distance, low (asteroidData.py's "distance_text" range)
    upper_limit_au       REAL NOT NULL,   -- distance, high
    density              TEXT NOT NULL CHECK (density IN ('dense', 'sparse', 'typical')),

    -- Searchable composition summary, built the same way as the prose
    -- sentence in AsteroidBelt.to_paragraph_list() (asteroidData.py:119-137),
    -- e.g. "high concentrations of iron, moderate concentrations of nickel,
    -- and trace amounts of platinum" -- a single column so a belt's
    -- composition can be searched without a join, alongside the
    -- structured per-component breakdown in asteroid_belt_composition
    -- below for queries that need one specific component.
    composition_summary  TEXT NOT NULL
);
CREATE INDEX idx_asteroid_belts_star_system_id ON asteroid_belts(star_system_id);

-- Structured per-component detail behind composition_summary above --
-- the belt's (component, concentration) list.
CREATE TABLE asteroid_belt_composition (
    id             INTEGER PRIMARY KEY,
    belt_id        INTEGER NOT NULL REFERENCES asteroid_belts(id) ON DELETE CASCADE,
    position       INTEGER NOT NULL,
    component      TEXT NOT NULL,
    concentration  TEXT NOT NULL CHECK (concentration IN ('high', 'moderate', 'small', 'trace'))
);
CREATE INDEX idx_asteroid_belt_composition_belt_id ON asteroid_belt_composition(belt_id);

-- ---------------------------------------------------------------------
-- sector_objects -- unified "every stellar object in a sector" search.
--
-- Standard relational choice here: a VIEW that UNIONs the three typed
-- tables above, not a fourth physical table duplicating their rows. A
-- real table would need to be kept in sync on every insert/update/delete
-- to the tables it mirrors (or drift out of sync); a view has no storage
-- and no sync problem -- SQLite resolves it against current data on every
-- query, and each row still traces back to its real table via
-- (object_type, object_id). Query with e.g.
-- `SELECT * FROM sector_objects WHERE sector_id = ?`.
-- ---------------------------------------------------------------------
CREATE VIEW sector_objects AS
    SELECT
        'star'                          AS object_type,
        s.id                            AS object_id,
        s.star_system_id                AS star_system_id,
        ss.sector_id                    AS sector_id,
        s.name                          AS name,
        s.table_type                    AS summary,
        NULL                            AS orbital_index
    FROM stars s
    JOIN star_systems ss ON ss.id = s.star_system_id

    UNION ALL

    SELECT
        'planet'                        AS object_type,
        p.id                            AS object_id,
        p.star_system_id                AS star_system_id,
        ss.sector_id                    AS sector_id,
        p.name                          AS name,
        COALESCE(p.table_class, p.body_type) AS summary,
        p.orbital_index                 AS orbital_index
    FROM planets p
    JOIN star_systems ss ON ss.id = p.star_system_id

    UNION ALL

    SELECT
        'asteroid_belt'                 AS object_type,
        ab.id                           AS object_id,
        ab.star_system_id               AS star_system_id,
        ss.sector_id                    AS sector_id,
        'Asteroid Belt'                 AS name,
        ab.density                      AS summary,
        ab.orbital_index                AS orbital_index
    FROM asteroid_belts ab
    JOIN star_systems ss ON ss.id = ab.star_system_id;
