-- stellarObjects/schema.sql
--
-- MySQL (InnoDB) schema for planetGen persistence (TODO.md Phase 2, ported
-- off SQLite for Phase 5's MySQL migration -- see docs/TODO.md). Plain SQL
-- DDL, no ORM. Column lists are verified against the actual generator
-- source (stellarObjects/config.py, starData.py, doubleStar.py,
-- planetData.py, asteroidData.py, spaceSector.py) as of schema version 2 --
-- not against TODO.md's earlier field-count estimates.
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
-- Distance-unit convention -- TWO standard units, by scale:
--   - SECTOR-SCALE PLACEMENT (`_mpc` suffix): star_systems.position_x/y/z
--     and sectors.edge -- where a system sits within its sector. Kept in
--     milliparsecs rather than km specifically because this is the one
--     place km's numbers get unwieldy (an 11.5-ly sector edge is
--     ~1.09e14 km vs. ~3526 mpc) and because sector geometry has no other
--     unit competing for consistency the way orbital distances do (every
--     orbital-scale quantity already shares km with radius/volume, so
--     there's no benefit standardizing sector geometry to km too).
--   - EVERYTHING ELSE (`_km`/`_km3` suffix): every other distance/length
--     column -- both "distance between things" (orbital distances,
--     habitable zones, system perimeter, heliosphere radius, hill radius,
--     scale height, binary separation) AND "size of an object" (radius,
--     volume) -- is kilometers, one standard unit instead of the
--     generator's native mix of km/AU, so search/comparison across
--     system-scale quantities never requires unit-aware query logic.
-- Every `table_*` column (the rendered wiki-table snapshot strings) is
-- unaffected by either convention -- those are copies of already-formatted
-- display text (still AU/ly/km as the wiki shows today), independent of
-- the raw column's storage unit. Conversion helpers live in
-- stellarObjects/utils.py: `ly_to_milliparsecs`/`milliparsecs_to_ly` for
-- the sector-scale columns; AU-to-km needs no helper for the km columns,
-- it's a single multiply by the existing `physical_constants.AU_TO_KM`.
-- Generation/physics code is untouched either way and keeps using its own
-- native km/AU/ly throughout; only the persistence layer converts, at the
-- moment of writing to (or reading from) these columns.
--
-- star_systems.quadrant stores the sector octant label (Roman numerals
-- I-VIII, one of the 8 sign(x)/sign(y)/sign(z) combinations -- see
-- spaceSector.py's `classify_octant`/`program_constants.SECTOR_OCTANT_LABELS`)
-- alongside the raw x/y/z -- purely derived from position, but worth
-- storing (not just recomputing on read) so it's directly
-- queryable/indexable without a UDF or generated-column expression. The
-- column name is kept as `quadrant` for schema stability, but every
-- human-facing label now displays it as "Octant" (html/sector.py,
-- html/system.py, html/static/sectormap.js) so it doesn't collide with
-- the unrelated, galaxy-scale "Quadrant" concept `sectors.center_x/y/z_pc`
-- and html/lib/galaxymap.py introduced later (4 azimuthal regions
-- spanning many sectors, not this column's 8 sign-combination regions
-- within one sector's own cube).
--
-- v3: star_systems.location stores a human-readable "sector name + nearest
-- neighbors" summary, e.g. "Voranthis Kelmoor -- nearest: Alpha Prime
-- (4.2 ly), Beta Cerise (7.8 ly), Gamma Ost (9.1 ly)" -- the same
-- "derive once at write time, persist for queryability" treatment
-- `quadrant` gets above, so a system's neighborhood is directly searchable
-- without recomputing nearest-neighbor distances on every read. NULL under
-- exactly the same condition as `quadrant`/`position_x/y/z_mpc`: a system
-- never placed in a sector. Computed in `stellarObjects/_db.py` from
-- `SpaceSector.nearest_neighbors` (up to 3, nearest first) at
-- `insert_sector` time; `migrate_database`'s `_migrate_v1_to_v2`/
-- `_migrate_v2_to_v3` backfill it for a database migrated from an older
-- version, applying the same nearest-neighbor logic directly to the
-- already-stored rows instead of a live `SpaceSector` object.
--
-- Two independent version numbers, per TODO.md's "Schema versioning"
-- decision:
--   - `schema_migrations` (below) is the DDL-level structure version --
--     replaces SQLite's `PRAGMA user_version` (no MySQL equivalent) with a
--     real tracking table, per TODO.md's Phase 5 MySQL-migration note.
--   - star_systems.schema_version (per row) is the version of the
--     serialized object-graph shape (Phase 1's to_dict()) that produced it.
-- Both started at 1; the DDL version below is CURRENT_SCHEMA_VERSION in
-- `_db.py`.
--
-- v4: added six nullable galaxy-frame placement columns to `sectors` --
-- center_x/y/z_pc (Cartesian center, parsecs, galactic-origin-relative),
-- galactic_radius_pc (sqrt(x^2+y^2+z^2), persisted for the same
-- can't-index-an-expression reason star_systems.quadrant is persisted
-- rather than recomputed), and shell_index/shell_slot_index (the
-- sector's stable address within the shell/Fibonacci-sphere tiling
-- scheme -- see docs/design/galaxy-coordinate-system.md). All six are
-- NULL together for a sector never placed in a galaxy -- every sector
-- generated by sectorGen.py's own standalone (non-galaxy) tooling, and
-- any migrated from an older database -- enforced by the CHECK below on
-- the first four (shell_index/shell_slot_index are deliberately excluded
-- from that CHECK: they're bookkeeping for the placement *algorithm*, not
-- physically implied by a center point, so a sector could in principle
-- have a hand-authored galaxy position without this scheme's own shell
-- addressing).
--
-- v2: moons split out of `planets` into their own `moons` table (with a
-- `planet_id` FK to the planet they orbit), instead of self-referencing
-- via `planets.parent_planet_id`/`is_moon` -- the same "give it its own
-- table" treatment `star_systems`/`planets` already have, so a query or a
-- UI facet can tell a Class M planet from a Class M moon without an
-- `is_moon` filter. `planet_evolutionary_paragraphs`/
-- `planet_reflection_spectrum` gained moon-owned counterparts
-- (`moon_evolutionary_paragraphs`/`moon_reflection_spectrum`) for the same
-- reason. Moons never generate their own moons (`Planet.__init__` only
-- calls `generate_moons` `if not self.is_moon`), so `moons` has no
-- self-reference of its own. `stellarObjects/_db.py`'s `migrate_database`
-- converts an existing v1 database in place (backing up the original
-- first); `migrateDb.py` runs it over every database in a directory, and
-- `install.sh` (and so `update.sh`, which calls it) does this on every
-- deploy.
-- v5: dropped every pre-rendered `table_*`/`binary_table_*` TEXT column
-- (stars.table_type/table_radius/table_mass/table_temp/table_lum/table_hab/
-- table_loc; planets/moons.table_class/table_distance/table_period/
-- table_radius/table_gravity; star_systems.binary_table_*). These held a
-- display string (e.g. "{{Exp|4.20|5}} kg (1.00% of Sol)") rendered once at
-- generation time from `Star`/`Planet`/`BinaryStarProxy.get_table_properties()`,
-- baked in whichever of wikitext-template or HTML form `SystemConfig.MARKDOWN`
-- happened to be at insert time -- wrong for any consumer wanting the other
-- form (the interactive HTML viewer in particular, which was displaying raw
-- unrendered `{{Exp|...}}` wikitext template syntax). No data is lost: every
-- number these strings were derived from already has its own proper
-- numeric/text column on the same row (mass_kg, radius_km, temperature_k,
-- luminosity_w, habitable_zone_inner/outer_km, distance_km, period_years,
-- gravity_g, star_type, name, planet_class, binary_effective_mass_kg,
-- binary_effective_luminosity_w, binary_separation_km, binary_type, ...) --
-- see this file's own header comment above. Consumers now compute display
-- formatting on demand from those columns instead of reading a frozen
-- pre-rendered copy (`html/lib/tabledisplay.py` for the HTML viewer;
-- `get_table_properties()` is still used, unchanged, to build the actual
-- wiki-page text in `wikitext_content`/`markdown_content`).
-- v6/v7: give every galaxy-placed sector exact vertices -- built from a
-- local spherical Voronoi tessellation among its own same-shell
-- neighbors, extruded radially between the shell's inner and outer
-- bounding spheres (`stellarObjects/sectorGeometry.prism_vertices`) --
-- genuinely gap-free against same-shell neighbors (not merely reduced;
-- see that module's docstring for why exact circumcenters, not nudged
-- approximations, make this possible) and area-matched (not vertex-
-- matched) against the shells in front of and behind it. Vertex/face
-- count varies per sector (typically 5-7, not a fixed number) because it
-- has to: a cube tiling of a sphere cannot be gap-free in general (the
-- same reason a soccer ball needs pentagons mixed with hexagons).
-- Supersedes an earlier, never-released "relax a fixed 8-vertex cube
-- toward its neighbors" approach, which only approximately closed gaps.
--   - v6 stored this as `sectors.vertices_pc`, a JSON `{"inner": [...],
--     "outer": [...]}` blob -- reconsidered almost immediately in favor
--     of v7's plain relational table below, so no released version ever
--     depended on the JSON shape.
--   - v7 replaces `sectors.vertices_pc` with the `sector_vertices` table
--     (see that table's own comment) -- one row per vertex, ordinary
--     columns throughout, no serialized blob anywhere in the schema.
--
-- v8: the galaxy-wide density "skeleton" (`galaxyPlan.py`), stored as
-- compact structural facts rather than one row per sector. Position,
-- density, and vertices are pure deterministic functions of `(shell_index,
-- shell_slot_index)` and a handful of galaxy-wide shape parameters
-- (`galaxyDensity.GalaxyShape`, `sectorGeometry.prism_vertices`) -- none
-- of that needs to be stored per candidate sector, since it's cheaper to
-- recompute on demand than to look up. What genuinely needs precomputing
-- is *where the galaxy has any content at all*, which is why this adds:
--   - `galaxy_shape`: one singleton row holding the whole galaxy's shape
--     parameters, its calibration constant, its sector edge length, and
--     its outer edge (the last shell with any qualifying content) -- the
--     handful of numbers `stellarObjects.galaxyDensity`/`galaxySkeleton`
--     need to recompute any sector's exact position/density on demand.
--   - `galaxy_shell_band`: one row per contiguous *candidate* slot-index
--     band per shell (`stellarObjects.galaxySkeleton.find_shell_bands`) --
--     the safe (possibly slightly wider than exact) range of slots that
--     might hold qualifying content in that shell, found via the exact
--     closed-form phi<->slot-index inverse
--     (`galaxyGeometry.slot_index_bounds_for_phi_range`) already used
--     elsewhere in this codebase. Deliberately NOT one row per qualifying
--     sector, and deliberately NOT storing that sector's position/density/
--     vertices -- an individual slot's exact qualification is checked live
--     (cheap: a single `relative_density` evaluation) only when that slot
--     is actually visited (see `galaxyGen.ensure_sector_generated`), which
--     is also the only time its full content, vertices, and everything
--     else about it get generated and stored in `sectors`/`sector_vertices`
--     /the star-system tables -- lazily, once, never recomputed again.
--     Almost every shell has exactly one band (this density model is
--     symmetric about and peaks at the galactic plane for any realistic
--     parameter choice -- verified directly across a full real-Milky-Way-
--     scale build), but the schema allows more than one per shell rather
--     than assuming it, since nothing about the model rules it out for
--     every possible parameter choice.
--   - `sectors` also gains a `UNIQUE (shell_index, shell_slot_index)`
--     constraint (see that table's own comment) -- now that a sector can
--     be lazily generated the moment it's visited
--     (`galaxyGen.ensure_sector_generated`), rather than only ever
--     through a single-process batch script, two concurrent visits to the
--     same never-before-generated address are possible; this constraint
--     turns that race into a clear `IntegrityError` the visit path
--     recovers from (re-fetch and return the sector the other visit just
--     created) instead of a silent duplicate row at the same address.
--
-- v9: orbital motion. `planets`/`moons` each gain three fixed-at-
--   generation-time orbital-orientation columns (`orbital_inclination_deg`,
--   `orbital_ascending_node_deg`) plus one that changes over time
--   (`orbital_phase_deg`, this body's current position angle around its
--   otherwise-circular orbit) and one static descriptive stat
--   (`rotation_period_hours`, axial "day length" -- this generator doesn't
--   track rotational phase, only period). `updateOrbits.py` is a new,
--   separately-run script that advances every body's `orbital_phase_deg`
--   in place based on `period_years` and real elapsed time, using the new
--   `orbit_simulation_state` singleton row (one per database, same
--   pattern as `galaxy_shape` above) to track when it last ran --
--   `stellarObjects.galaxyGeometry`'s "derive, don't store" principle
--   doesn't apply here: the whole point of this feature is a value that
--   changes with real-world time rather than being a pure function of the
--   body's other (fixed) generative facts, so it has to be stored and
--   periodically updated, not recomputed on read.
--
-- MySQL port -- type mapping and idempotency notes (TODO.md Phase 5):
--   - SQLite's `INTEGER PRIMARY KEY` (a 64-bit rowid alias) becomes
--     `BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY` throughout, with every
--     referencing FK column matching that same `BIGINT UNSIGNED` type --
--     InnoDB requires matching column types/signedness for a foreign key,
--     and BIGINT headroom means this never needs revisiting even at
--     Phase 4 galaxy scale (tens of millions of stars/planets/moons).
--   - SQLite's `REAL` becomes `DOUBLE` (MySQL's `REAL` is a synonym for
--     `DOUBLE` under default SQL modes, but this spells it out rather than
--     relying on a server-side mode setting).
--   - SQLite's untyped `TEXT` splits by actual content length: short
--     identifier/enum-like/name fields become `VARCHAR(n)` (indexable,
--     and MySQL's utf8mb4 index prefix limits make an unbounded TEXT index
--     impractical anyway); genuinely long generated prose (descriptions,
--     flavor text, composition summaries, evolutionary paragraphs) stays
--     `TEXT`; the two full-page rendered columns (`wikitext_content`/
--     `markdown_content`) that can exceed `TEXT`'s 64KiB ceiling for a
--     large system become `LONGTEXT`.
--   - `created_at` becomes a native `TIMESTAMP NOT NULL DEFAULT
--     CURRENT_TIMESTAMP` instead of `TEXT DEFAULT CURRENT_TIMESTAMP`.
--   - Every `CHECK` constraint carries over unchanged (MySQL enforces
--     `CHECK` from 8.0.16 -- this schema targets MySQL 8.0.16+ /
--     compatible MariaDB); multi-column CHECKs (`sectors`) are equally
--     supported.
--   - SQLite's column-level `REFERENCES` shorthand is parsed but NOT
--     enforced by MySQL/InnoDB -- every foreign key here is instead a
--     named table-level `CONSTRAINT ... FOREIGN KEY (...) REFERENCES
--     ...` clause, which InnoDB does enforce.
--   - SQLite's separate `CREATE INDEX IF NOT EXISTS` statements don't
--     translate directly: MySQL's `CREATE INDEX` has no `IF NOT EXISTS`
--     form, which would make re-running this script against an existing
--     database fail with "Duplicate key name" on the second run. Every
--     index is instead declared inline as a `KEY` clause inside its
--     table's `CREATE TABLE IF NOT EXISTS` -- since the whole table
--     (indexes included) is only ever created once, this script stays
--     safe to run repeatedly against an existing database, matching the
--     SQLite schema's original idempotency guarantee.
--   - MySQL has no `CREATE VIEW IF NOT EXISTS`; `sector_objects` below
--     uses `CREATE OR REPLACE VIEW` instead -- equally idempotent (and
--     "replace" is the more correct intent for a view, whose definition
--     might legitimately need updating across a migration in a way a
--     physical table's columns wouldn't).
--   - Every table is explicit `ENGINE=InnoDB` (foreign keys require it)
--     and `utf8mb4`/`utf8mb4_unicode_ci` (full Unicode, including
--     generated names/flavor text outside the Basic Multilingual Plane,
--     rather than MySQL's legacy 3-byte `utf8`).

-- ---------------------------------------------------------------------
-- schema_migrations -- DDL-level structure version tracking. Replaces
-- SQLite's `PRAGMA user_version` (see the MySQL port note above);
-- `_db.py`'s `migrate_database` reads `MAX(version)` from this table
-- and inserts a row per migration step it applies.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS schema_migrations (
    version      INT NOT NULL PRIMARY KEY,
    applied_at   TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- sectors
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS sectors (
    id                  BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    name                VARCHAR(255) NOT NULL,
    edge_mpc            DOUBLE NOT NULL,

    -- Galaxy-frame placement (v4) -- see the header comment's "v4" note.
    -- NULL together: a sector never placed in a galaxy (today's
    -- sectorGen.py standalone tooling, or migrated from a pre-v4
    -- database).
    center_x_pc         DOUBLE,
    center_y_pc         DOUBLE,
    center_z_pc         DOUBLE,
    galactic_radius_pc  DOUBLE,

    -- The sector's stable (shell, slot) address within the shell/
    -- Fibonacci-sphere tiling scheme -- deliberately not part of the
    -- CHECK below (see the header comment's "v4" note).
    shell_index         INT,
    shell_slot_index    INT,

    CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    ),

    -- v8: at most one sector per (shell_index, shell_slot_index) address --
    -- see the header comment's "v8" note (galaxyGen.ensure_sector_generated
    -- relies on this to make a lazy-generation race produce a clear
    -- IntegrityError rather than a silent duplicate sector at the same
    -- address). NULL-together sectors (never placed in a galaxy) don't
    -- collide with each other or with a placed sector -- ordinary SQL NULL
    -- semantics for UNIQUE.
    UNIQUE (shell_index, shell_slot_index),

    KEY idx_sectors_galactic_radius_pc (galactic_radius_pc),
    KEY idx_sectors_shell_index (shell_index)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- This sector's own exact prism vertices (v7 -- see the header comment's
-- "v6/v7" note), one row per vertex rather than a serialized blob:
-- `sectorGeometry.prism_vertices` returns a variable number of vertices
-- per sector (typically 5-7, not fixed), split into an "inner" ring (on
-- the shell's inner bounding sphere) and an "outer" ring (on its outer
-- bounding sphere), both in the same cyclic order. `vertex_index` is that
-- cyclic position (0-based) within its own ring, not a global ordering --
-- pairing `(sector_id, vertex_index)` across the two rings gives the
-- lateral edge each pair of inner/outer vertices spans. No row exists for
-- a sector never placed in a galaxy (mirrors, at the application level
-- rather than a cross-table CHECK -- SQLite can't express "rows exist in
-- another table" as a CHECK constraint -- the same NULL-together
-- condition `sectors`'s own galaxy-placement columns enforce directly).
CREATE TABLE IF NOT EXISTS sector_vertices (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id     BIGINT UNSIGNED NOT NULL,
    ring          VARCHAR(8) NOT NULL CHECK (ring IN ('inner', 'outer')),
    vertex_index  INT NOT NULL,
    x_pc          DOUBLE NOT NULL,
    y_pc          DOUBLE NOT NULL,
    z_pc          DOUBLE NOT NULL,

    CONSTRAINT fk_sector_vertices_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    UNIQUE (sector_id, ring, vertex_index),
    KEY idx_sector_vertices_sector_id (sector_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- galaxy_shape / galaxy_shell_band -- the galaxy-wide density "skeleton"
-- (v8 -- see the header comment's "v8" note). `galaxy_shape` is a
-- singleton (`id` pinned to 1, enforced by the CHECK) -- there is exactly
-- one galaxy. Building or rebuilding the skeleton (`galaxyPlan.py`)
-- replaces this row and every `galaxy_shell_band` row wholesale; neither
-- table is ever partially updated.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS galaxy_shape (
    id                          BIGINT UNSIGNED PRIMARY KEY CHECK (id = 1),

    -- galaxyDensity.GalaxyShape fields, verbatim -- see that module for
    -- what each one means and how `k_norm` is calibrated.
    disk_scale_length_pc        DOUBLE NOT NULL,
    disk_scale_height_pc        DOUBLE NOT NULL,
    bulge_scale_radius_pc       DOUBLE NOT NULL,
    bulge_amplitude             DOUBLE NOT NULL,
    arm_count                   INT NOT NULL,
    pitch_angle_rad             DOUBLE NOT NULL,
    arm_amplitude               DOUBLE NOT NULL,
    spiral_reference_radius_pc  DOUBLE NOT NULL,
    spiral_reference_angle_rad  DOUBLE NOT NULL,
    k_norm                      DOUBLE NOT NULL,

    -- The sector edge length this skeleton was built at, in parsecs
    -- (galaxyGeometry's own scope -- every shell/sector-address formula
    -- takes this as a parameter rather than assuming a fixed constant).
    edge_pc                     DOUBLE NOT NULL,

    -- SpaceSector(edge_ly=...).expected_system_count() at relative_density
    -- = 1 -- cached because every qualification check
    -- (predicted_star_count >= 1.0) needs it, and it's a fixed galaxy-wide
    -- fact, not worth re-deriving from edge_pc on every call.
    expected_system_count_at_density_1  DOUBLE NOT NULL,

    -- The last shell index with any qualifying content -- this galaxy's
    -- real edge, discovered by `galaxyPlan.py` (a run of consecutive empty
    -- shells beyond it), not picked as an arbitrary radius.
    outer_shell_index           INT NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- One contiguous candidate slot-index band per shell (almost always
-- exactly one -- see the header comment's "v8" note). `band_index` orders
-- multiple bands within the same shell (0-based); a shell with no
-- qualifying content at all has no rows here.
CREATE TABLE IF NOT EXISTS galaxy_shell_band (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    shell_index      INT NOT NULL,
    band_index       INT NOT NULL,
    slot_index_min   INT NOT NULL,
    slot_index_max   INT NOT NULL,

    UNIQUE (shell_index, band_index),
    KEY idx_galaxy_shell_band_shell_index (shell_index)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Singleton row (same pattern as galaxy_shape above) tracking when
-- updateOrbits.py last advanced every planet's/moon's orbital_phase_deg in
-- this database, so the next run knows how much real time has actually
-- elapsed since then. Absent entirely until updateOrbits.py's first run
-- against a given database (it creates this row itself).
CREATE TABLE IF NOT EXISTS orbit_simulation_state (
    id               BIGINT UNSIGNED PRIMARY KEY CHECK (id = 1),
    last_updated_at  TIMESTAMP NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- system_configs -- one row per SystemConfig "recipe"
-- (stellarObjects/config.py:29-33 SERIALIZABLE_FIELDS)
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS system_configs (
    id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    markdown          TINYINT(1) NOT NULL DEFAULT 0 CHECK (markdown IN (0, 1)),
    habitable_world   TINYINT(1) CHECK (habitable_world IN (0, 1)),
    asteroid_belt     TINYINT(1) CHECK (asteroid_belt IN (0, 1)),
    large_star        TINYINT(1) CHECK (large_star IN (0, 1)),
    moons             TINYINT(1) CHECK (moons IN (0, 1)),
    max_planets       TINYINT(1) CHECK (max_planets IN (0, 1)),
    planets           TINYINT(1) CHECK (planets IN (0, 1)),
    star_type         VARCHAR(64),
    name              VARCHAR(255),
    age               VARCHAR(16) CHECK (age IN ('young', 'old')),
    intelligent_life  TINYINT(1) CHECK (intelligent_life IN (0, 1)),
    binary_system     TINYINT(1) CHECK (binary_system IN (0, 1)),
    num_orbits        INT
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Child table for the variable-length SLOTS recipe list
-- (config.py:135-151).
CREATE TABLE IF NOT EXISTS system_config_slots (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    config_id     BIGINT UNSIGNED NOT NULL,
    orbit_index   INT NOT NULL,
    type          VARCHAR(16) CHECK (type IN ('planet', 'asteroid_belt')),
    planet_class  VARCHAR(16),
    moons         INT,

    CONSTRAINT fk_system_config_slots_config
        FOREIGN KEY (config_id) REFERENCES system_configs(id) ON DELETE CASCADE,
    KEY idx_system_config_slots_config_id (config_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- star_systems
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS star_systems (
    id                     BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id              BIGINT UNSIGNED,
    system_config_id       BIGINT UNSIGNED NOT NULL,
    name                   VARCHAR(255) NOT NULL,

    -- Position within its sector (spaceSector.py SectorSystemEntry.position),
    -- milliparsecs, relative to the sector's cubic center. NULL iff this
    -- system was never placed in a sector.
    position_x_mpc          DOUBLE,
    position_y_mpc          DOUBLE,
    position_z_mpc          DOUBLE,

    -- Sector octant label derived from the position above (Roman numeral
    -- I-VIII, see the header comment) -- NULL iff position is NULL.
    quadrant                VARCHAR(4) CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')),

    -- Human-readable "sector name + nearest neighbors" summary, see the
    -- header comment's "v3" note -- NULL iff position is NULL.
    location                TEXT,

    is_binary              TINYINT(1) NOT NULL DEFAULT 0 CHECK (is_binary IN (0, 1)),

    -- BinaryStarProxy-derived fields (doubleStar.py) -- all NULL for a
    -- single-star system, stored rather than re-derived since
    -- _effective_mass/_effective_luminosity are computed once at
    -- generation time.
    binary_separation_km            DOUBLE,
    binary_type                     VARCHAR(64),
    binary_temperature_k            DOUBLE,
    binary_radius_km                DOUBLE,
    binary_effective_mass_kg        DOUBLE,
    binary_effective_luminosity_w   DOUBLE,
    binary_age_gy                   DOUBLE,
    binary_lifespan_gy              DOUBLE,  -- NULL = float('inf')
    binary_habitable_zone_inner_km  DOUBLE,
    binary_habitable_zone_outer_km  DOUBLE,
    binary_system_perimeter_km      DOUBLE,
    binary_heliosphere_radius_km    DOUBLE,

    system_flavor_text   TEXT,
    schema_version       INT NOT NULL DEFAULT 1,

    -- Fully rendered page text, ready to paste/upload. Both must be
    -- rendered from the SAME already-generated StarSystem object
    -- (toggle system_config.MARKDOWN, render, toggle back, render again)
    -- -- generation mixes unseedable `secrets` with seedable `random`
    -- (see spaceSector.py's module docstring), so the other format can
    -- never be faithfully regenerated later. LONGTEXT, not TEXT -- see
    -- the MySQL port note above.
    wikitext_content     LONGTEXT,
    markdown_content     LONGTEXT,

    -- One system = one wiki page (stars/planets/moons are sections within
    -- it, per StarSystem.__str__), so exactly one URL per wiki target.
    mediawiki_url        VARCHAR(2048),
    wikijs_url           VARCHAR(2048),

    created_at           TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT fk_star_systems_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE SET NULL,
    CONSTRAINT fk_star_systems_config
        FOREIGN KEY (system_config_id) REFERENCES system_configs(id),
    KEY idx_star_systems_sector_id (sector_id),
    KEY idx_star_systems_system_config_id (system_config_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- stars -- one row per single star, two rows (primary/secondary) per
-- binary. Never a row for the BinaryStarProxy itself (its snapshot lives
-- on star_systems above). Fields verified against starData.py:399-445 /
-- Star.__init__.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS stars (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    star_system_id            BIGINT UNSIGNED NOT NULL,
    role                      VARCHAR(16) NOT NULL CHECK (role IN ('primary', 'secondary', 'single')),
    name                      VARCHAR(255) NOT NULL,
    star_type                 VARCHAR(64) NOT NULL,   -- Star.type, e.g. "G2V Yellow Main Sequence Star"
    yerkes_class              VARCHAR(16) NOT NULL,
    mass_kg                   DOUBLE NOT NULL,
    radius_km                 DOUBLE NOT NULL,
    temperature_k             DOUBLE NOT NULL,
    luminosity_w              DOUBLE NOT NULL,
    age_gy                    DOUBLE NOT NULL,
    lifespan_gy               DOUBLE,            -- NULL = float('inf'), white dwarfs
    habitable_zone_inner_km   DOUBLE NOT NULL,
    habitable_zone_outer_km   DOUBLE NOT NULL,
    system_perimeter_km       DOUBLE NOT NULL,
    heliosphere_radius_km     DOUBLE NOT NULL,

    CONSTRAINT fk_stars_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    KEY idx_stars_star_system_id (star_system_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- planets -- top-level planets only as of schema v2 (see the header
-- comment's "v2" note); moons live in the `moons` table below instead of
-- self-referencing here. Fields verified against planetData.py:117-200 /
-- Planet.__init__ (27 scalar fields; excludes system_config/star
-- back-refs, the `moons` list [-> the moons table], `evolutionary_data`
-- [-> planet_evolutionary_paragraphs], and reflection_spectrum_visible/
-- non_visible [-> planet_reflection_spectrum]).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS planets (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    star_system_id            BIGINT UNSIGNED NOT NULL,
    -- The specific star this planet orbits, when that's a real stored
    -- `stars` row -- true for every single-star system. NULL for a
    -- binary's planets: systemData.py always generates planets against
    -- `self.star`, which for binaries is the `BinaryStarProxy` (never one
    -- individual constituent star -- there is no S-type/circumbinary
    -- choice in the current generator), and the proxy is deliberately not
    -- stored as its own `stars` row (see star_systems.binary_* above).
    -- `star_system_id` above is always the reliable owning link
    -- regardless of star_id.
    star_id                   BIGINT UNSIGNED,
    orbital_index             INT NOT NULL,
    body_type                 VARCHAR(4) NOT NULL CHECK (body_type IN ('t', 'g')),
    name                      VARCHAR(255) NOT NULL,
    planet_class              VARCHAR(16),
    distance_km               DOUBLE NOT NULL,   -- from the star
    radius_km                 DOUBLE NOT NULL,
    mass_kg                   DOUBLE NOT NULL,
    volume_km3                DOUBLE NOT NULL,   -- derived, stored not recomputed
    period_years              DOUBLE NOT NULL,   -- derived, stored not recomputed
    zone                      VARCHAR(4) CHECK (zone IN ('h', 'e', 'c')),
    description               TEXT,
    gravity_g                 DOUBLE,
    surface_temperature_k     DOUBLE,
    density_g_cm3             DOUBLE,
    atmosphere                TEXT,
    atm_density               DOUBLE,
    atm_molar_density         DOUBLE,
    atmospheric_pressure_pa   DOUBLE,
    composition               TEXT,            -- descriptive string (contrast asteroid_belt_composition)
    scale_height_km           DOUBLE,
    hill_radius_km            DOUBLE,
    min_orbit_distance_km     DOUBLE,
    habitable_zone_inner_km   DOUBLE NOT NULL,   -- copied from host star at generation time
    habitable_zone_outer_km   DOUBLE NOT NULL,
    life_chemical             VARCHAR(64),
    evolutionary_speed        VARCHAR(64),
    flavor_text               TEXT,
    flavor_text_count         INT NOT NULL DEFAULT 0,
    -- Orbital motion (v9, see header comment) -- inclination/ascending
    -- node are fixed at generation time; phase changes over time, advanced
    -- in place by updateOrbits.py. rotation_period_hours is a separate,
    -- static "day length" stat.
    orbital_inclination_deg     DOUBLE NOT NULL,
    orbital_ascending_node_deg  DOUBLE NOT NULL,
    orbital_phase_deg           DOUBLE NOT NULL,
    rotation_period_hours       DOUBLE NOT NULL,

    CONSTRAINT fk_planets_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_planets_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL,
    KEY idx_planets_star_system_id (star_system_id),
    KEY idx_planets_star_id (star_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Child table for the variable-length evolutionary narrative list
-- (Planet.evolutionary_data).
CREATE TABLE IF NOT EXISTS planet_evolutionary_paragraphs (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    planet_id   BIGINT UNSIGNED NOT NULL,
    position    INT NOT NULL,
    paragraph   TEXT NOT NULL,

    CONSTRAINT fk_planet_evolutionary_paragraphs_planet
        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
    KEY idx_planet_evolutionary_paragraphs_planet_id (planet_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Child table replacing reflection_spectrum_visible/non_visible -- no
-- JSON column, same position+value shape as the paragraphs/composition
-- child tables (this schema is a pure relational design, not a
-- JSON-blob-in-a-column hybrid).
CREATE TABLE IF NOT EXISTS planet_reflection_spectrum (
    id             BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    planet_id      BIGINT UNSIGNED NOT NULL,
    spectrum_type  VARCHAR(16) NOT NULL CHECK (spectrum_type IN ('visible', 'non_visible')),
    position       INT NOT NULL,
    value          VARCHAR(255) NOT NULL,

    CONSTRAINT fk_planet_reflection_spectrum_planet
        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
    KEY idx_planet_reflection_spectrum_planet_id (planet_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- moons -- one row per moon, as of schema v2 (see the header comment's
-- "v2" note). Exactly the same column shape as `planets` (moons are
-- `Planet` instances too, with `is_moon=True`), except `planet_id`
-- replaces `star_system_id`/`star_id`'s role as "what this body directly
-- orbits" -- `star_system_id`/`star_id` are still carried here too
-- (redundant with the owning planet's own columns) purely so a moon can
-- be queried/joined to its system or star without an extra hop through
-- `planets`. No self-reference: moons never generate their own moons
-- (`Planet.__init__` only calls `generate_moons` `if not self.is_moon`).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS moons (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    planet_id                 BIGINT UNSIGNED NOT NULL,
    star_system_id            BIGINT UNSIGNED NOT NULL,
    star_id                   BIGINT UNSIGNED,
    orbital_index             INT NOT NULL,   -- position in the parent planet's moons list
    body_type                 VARCHAR(4) NOT NULL CHECK (body_type IN ('t', 'g')),
    name                      VARCHAR(255) NOT NULL,
    planet_class              VARCHAR(16),
    distance_km               DOUBLE NOT NULL,   -- from the parent planet
    radius_km                 DOUBLE NOT NULL,
    mass_kg                   DOUBLE NOT NULL,
    volume_km3                DOUBLE NOT NULL,
    period_years              DOUBLE NOT NULL,
    zone                      VARCHAR(4) CHECK (zone IN ('h', 'e', 'c')),
    description               TEXT,
    gravity_g                 DOUBLE,
    surface_temperature_k     DOUBLE,
    density_g_cm3             DOUBLE,
    atmosphere                TEXT,
    atm_density               DOUBLE,
    atm_molar_density         DOUBLE,
    atmospheric_pressure_pa   DOUBLE,
    composition               TEXT,
    scale_height_km           DOUBLE,
    hill_radius_km            DOUBLE,
    min_orbit_distance_km     DOUBLE,
    habitable_zone_inner_km   DOUBLE NOT NULL,
    habitable_zone_outer_km   DOUBLE NOT NULL,
    life_chemical             VARCHAR(64),
    evolutionary_speed        VARCHAR(64),
    flavor_text               TEXT,
    flavor_text_count         INT NOT NULL DEFAULT 0,
    -- Orbital motion (v9, see header comment) -- inclination/ascending
    -- node are fixed at generation time; phase changes over time, advanced
    -- in place by updateOrbits.py. rotation_period_hours is a separate,
    -- static "day length" stat.
    orbital_inclination_deg     DOUBLE NOT NULL,
    orbital_ascending_node_deg  DOUBLE NOT NULL,
    orbital_phase_deg           DOUBLE NOT NULL,
    rotation_period_hours       DOUBLE NOT NULL,

    CONSTRAINT fk_moons_planet
        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
    CONSTRAINT fk_moons_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_moons_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL,
    KEY idx_moons_planet_id (planet_id),
    KEY idx_moons_star_system_id (star_system_id),
    KEY idx_moons_star_id (star_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Moon counterpart of planet_evolutionary_paragraphs.
CREATE TABLE IF NOT EXISTS moon_evolutionary_paragraphs (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    moon_id     BIGINT UNSIGNED NOT NULL,
    position    INT NOT NULL,
    paragraph   TEXT NOT NULL,

    CONSTRAINT fk_moon_evolutionary_paragraphs_moon
        FOREIGN KEY (moon_id) REFERENCES moons(id) ON DELETE CASCADE,
    KEY idx_moon_evolutionary_paragraphs_moon_id (moon_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Moon counterpart of planet_reflection_spectrum.
CREATE TABLE IF NOT EXISTS moon_reflection_spectrum (
    id             BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    moon_id        BIGINT UNSIGNED NOT NULL,
    spectrum_type  VARCHAR(16) NOT NULL CHECK (spectrum_type IN ('visible', 'non_visible')),
    position       INT NOT NULL,
    value          VARCHAR(255) NOT NULL,

    CONSTRAINT fk_moon_reflection_spectrum_moon
        FOREIGN KEY (moon_id) REFERENCES moons(id) ON DELETE CASCADE,
    KEY idx_moon_reflection_spectrum_moon_id (moon_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- asteroid_belts -- belts have no properties-dict data table
-- (asteroidData.py:93-141 is prose only), so their searchable columns are
-- the facts that prose always states instead: density, the distance range
-- from the star, and composition.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS asteroid_belts (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    star_system_id       BIGINT UNSIGNED NOT NULL,
    orbital_index        INT NOT NULL,
    distance_km          DOUBLE NOT NULL,
    lower_limit_km       DOUBLE NOT NULL,   -- distance, low (asteroidData.py's "distance_text" range)
    upper_limit_km       DOUBLE NOT NULL,   -- distance, high
    density              VARCHAR(16) NOT NULL CHECK (density IN ('dense', 'sparse', 'typical')),

    -- Searchable composition summary, built the same way as the prose
    -- sentence in AsteroidBelt.to_paragraph_list() (asteroidData.py:119-137),
    -- e.g. "high concentrations of iron, moderate concentrations of nickel,
    -- and trace amounts of platinum" -- a single column so a belt's
    -- composition can be searched without a join, alongside the
    -- structured per-component breakdown in asteroid_belt_composition
    -- below for queries that need one specific component.
    composition_summary  TEXT NOT NULL,

    CONSTRAINT fk_asteroid_belts_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    KEY idx_asteroid_belts_star_system_id (star_system_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Structured per-component detail behind composition_summary above --
-- the belt's (component, concentration) list.
CREATE TABLE IF NOT EXISTS asteroid_belt_composition (
    id             BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    belt_id        BIGINT UNSIGNED NOT NULL,
    position       INT NOT NULL,
    component      VARCHAR(64) NOT NULL,
    concentration  VARCHAR(16) NOT NULL CHECK (concentration IN ('high', 'moderate', 'small', 'trace')),

    CONSTRAINT fk_asteroid_belt_composition_belt
        FOREIGN KEY (belt_id) REFERENCES asteroid_belts(id) ON DELETE CASCADE,
    KEY idx_asteroid_belt_composition_belt_id (belt_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- sector_objects -- unified "every stellar object in a sector" search.
--
-- Standard relational choice here: a VIEW that UNIONs the four typed
-- tables above, not a fifth physical table duplicating their rows. A
-- real table would need to be kept in sync on every insert/update/delete
-- to the tables it mirrors (or drift out of sync); a view has no storage
-- and no sync problem -- MySQL resolves it against current data on every
-- query, and each row still traces back to its real table via
-- (object_type, object_id). Query with e.g.
-- `SELECT * FROM sector_objects WHERE sector_id = %s`.
-- ---------------------------------------------------------------------
CREATE OR REPLACE VIEW sector_objects AS
    SELECT
        'star'                          AS object_type,
        s.id                            AS object_id,
        s.star_system_id                AS star_system_id,
        ss.sector_id                    AS sector_id,
        s.name                          AS name,
        s.star_type                     AS summary,
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
        COALESCE(p.planet_class, p.body_type) AS summary,
        p.orbital_index                 AS orbital_index
    FROM planets p
    JOIN star_systems ss ON ss.id = p.star_system_id

    UNION ALL

    SELECT
        'moon'                          AS object_type,
        m.id                            AS object_id,
        m.star_system_id                AS star_system_id,
        ss.sector_id                    AS sector_id,
        m.name                          AS name,
        COALESCE(m.planet_class, m.body_type) AS summary,
        m.orbital_index                 AS orbital_index
    FROM moons m
    JOIN star_systems ss ON ss.id = m.star_system_id

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
