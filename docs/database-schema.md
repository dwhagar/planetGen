# planetGen Database Format

This document describes the MySQL database schema defined in
[`src/stellarObjects/schema.sql`](../src/stellarObjects/schema.sql). It's the
reference for anyone reading, querying, or extending the database — table by
table, every column's meaning and unit, and the conventions that hold the
schema together.

**Status**: the schema is implemented, with both a write and a read path.
[`src/stellarObjects/_db.py`](../src/stellarObjects/_db.py) (private — leading
underscore, not part of the package's public generation API) writes
already-generated `StarSystem`/`SpaceSector` objects straight into these
tables (`generate.py sector` calls it automatically on every run), and
reconstructs them back into live objects from rows (`load_star_system`/
`load_sector`/`load_system_config`), inverting every unit conversion the
write path applies. `queryDb.py` is the "list what's stored" CLI (`TODO.md`
Phase 3). The database itself lives on a MySQL server (TODO.md Phase 5 —
this schema previously targeted SQLite; see "MySQL port" in `schema.sql`'s
own header comment for the type-mapping/idempotency notes that move brought),
reachable via the `$PLANETGEN_MYSQL_HOST`/`$PLANETGEN_MYSQL_PORT`/
`$PLANETGEN_MYSQL_USER`/`$PLANETGEN_MYSQL_PASSWORD`/`$PLANETGEN_MYSQL_DATABASE`
environment variables (or the equivalent `--mysql-*` CLI flags every entry
point accepts — see `stellarObjects._db.MySQLConfig`), with its tables created
automatically on first connection. See `TODO.md` for the full roadmap.

## Persistence layer

`src/stellarObjects/_db.py` owns every unit conversion at the point of writing
a value into the database; nothing else in the package imports it, and it
never mutates the generation/physics code's own native units. Its shape:

- `get_connection(config=None)` / `save_sector(sector, config=None)` —
  the two entry points most callers need. `save_sector` connects (pooled —
  see `MySQLConfig`/`_get_pool`), applies the schema if needed (idempotent —
  `CREATE TABLE IF NOT EXISTS`/`CREATE OR REPLACE VIEW` throughout
  `schema.sql`), and writes the whole sector in one transaction.
  `config` (a `MySQLConfig`) defaults to `DEFAULT_MYSQL_CONFIG`, itself
  built from the `$PLANETGEN_MYSQL_*` environment variables above.
- `insert_sector` / `insert_star_system` / `insert_star` / `insert_planet`
  (calls `insert_moon` for each of a planet's moons) / `insert_moon` /
  `insert_asteroid_belt` / `insert_system_config` — the per-table building
  blocks, callable individually for a standalone system with no sector
  (`insert_star_system(conn, star_system, system_config, sector_id=None,
  position=None)`).
- `migrate_database(config=None)` — brings a database's `schema_migrations`
  bookkeeping up to `SCHEMA_VERSION`, applying any migration step in
  between. A brand-new database is created at the current schema, so this
  is a no-op there; a database created by an older release gets one
  `_migrate_vN_to_vN+1` step per version it is behind (see "Versioning"
  below). `src/migrateDb.py` is its CLI. There is no import from a
  pre-MySQL-port SQLite database (see "Versioning" below).
- `load_star_system(conn, star_system_id)` / `load_sector(conn, sector_id)`
  / `load_system_config(conn, config_id)` — the read-path counterparts,
  reconstructing a live `StarSystem`/`SpaceSector`/`SystemConfig` from
  rows. Each row is mapped to the flat dict shape the corresponding
  class's own `from_dict` already accepts (`Star`/`BinaryStarProxy`/
  `Planet`/`AsteroidBelt`/`SystemConfig`, from `TODO.md`'s Phase 1
  serialization), rather than a second, independent object-construction
  path — `StarSystem` itself is assembled directly (mirroring
  `StarSystem.from_dict`'s own wiring) since the data is normalized across
  many flat tables rather than naturally nested the way a JSON export is.
- Every `composition_summary` column is read from a
  `get_composition_summary()` method (`AsteroidBelt`, asteroid fields,
  comets), the same text the rendered page uses. The old `table_*`
  display-string columns were dropped in v5 and the stored page text in
  v29: display formatting is computed on demand from the data columns
  (`html/lib/tabledisplay.py`, `stellarObjects/systemRender.py`).

## How to read this document

- Every table is listed with its columns, each column's type, nullability,
  and unit (where relevant), and a short note on where the value comes from
  in the generator.
- `FK -> table.column` marks a foreign key.
- Types are given loosely in the older tables below (`INTEGER`, `TEXT`).
  In MySQL every `id` and foreign key is `BIGINT UNSIGNED` (ids
  `AUTO_INCREMENT`), most short strings (names, types, classes) are
  `VARCHAR(n)`, booleans are `TINYINT(1)`, and URL columns are
  `VARCHAR(2048)`. `schema.sql` is authoritative for exact types.
- Conventions that apply across many tables (units, the searchable-field
  principle, versioning) are explained once, up front, rather than repeated
  per table.

## Conventions

### Two distance units, by scale

Every distance/length-shaped column uses one of two units, chosen by scale
rather than forced into a single unit everywhere:

- **Milliparsecs (`_mpc` suffix)** — sector-scale placement only:
  `star_systems.position_x/y/z_mpc` and `sectors.edge_mpc`. This is
  specifically where kilometers get unwieldy (an 11.5-light-year sector
  edge is ~1.09×10^14 km, but only ~3526 mpc), and sector geometry doesn't
  share a unit with anything else the way orbital-scale quantities share
  kilometers with radius.
- **Kilometers (`_km`/`_km3` suffix)** — everything else: orbital
  distances, habitable zones, system perimeter, heliosphere radius, hill
  radius, scale height, binary separation, asteroid belt distances, *and*
  the size of an object (radius, volume). One unit for all of it means
  comparing or filtering across these columns never needs unit-aware query
  logic.

Neither convention touches the generator itself — every attribute in
`src/stellarObjects/starData.py`, `doubleStar.py`, `planetData.py`,
`asteroidData.py`, and `spaceSector.py` keeps its own native unit (km, AU,
or ly) exactly as today. Conversion only happens at the persistence
boundary, once it's built: `src/stellarObjects/utils.py` provides
`ly_to_milliparsecs`/`milliparsecs_to_ly` for the sector-scale columns;
AU-to-km needs no helper, since it's a single multiply by the existing
`physical_constants.AU_TO_KM`.

The `table_*` columns that once held copies of already-formatted display
text (e.g. `"1.2 R☉"`) were the one exception to both conventions; v5
dropped them (see "The searchable-field principle" below).

### The searchable-field principle

*Historical: v5 dropped every `table_*`/`binary_table_*` column described
in the next paragraph, and display strings are now computed on demand
from the raw columns (`html/lib/tabledisplay.py`). The belt columns in
the second paragraph remain.*

Every `*Data`/`*_properties` dict that `to_paragraph_list()` builds in the
generator — the exact values that appear in each object's rendered wiki
table — gets its own set of `table_*` columns (one column per dict key),
rather than being reconstructed from the raw scalar columns at query time.
This matters because the raw scalar columns don't always match the display
text one-to-one (units, rounding, phrasing can differ or change across
versions); the `table_*` columns are the "as-published" snapshot, kept
alongside the raw generative fields for full object-graph fidelity.

Asteroid belts have no such dict (`AsteroidBelt.to_paragraph_list()` is
prose only) — their searchable columns instead capture the same facts the
prose always states: `density`, the distance range (`lower_limit_km`/
`upper_limit_km`), and `composition_summary`.

### Versioning

Two independent version numbers:

- `schema_migrations` (one row per applied DDL migration step) — the DDL
  structure version, `MAX(version)` in that table (this schema is version
  `51`, `_db.SCHEMA_VERSION`). Replaces SQLite's `PRAGMA user_version`, which has no MySQL
  equivalent — see `schema.sql`'s "MySQL port" header note.
- `star_systems.schema_version` (per row) — the version of the serialized
  object-graph shape (Phase 1's `to_dict()`) that produced that row.
  Independent of the DDL version because a JSON export/import could bring
  an older object-graph shape into a newer database.

The schema evolved through several versions while still SQLite-backed;
each version's structural change is recorded in `schema.sql`'s own header
comment ("v2" through "v50" notes) rather than duplicated here, since that
file is the one place both the current column list and the historical
rationale for it live together. In brief: v1→v2 split moons out of the
shared `planets` table into their own `moons` table; v2→v3 added
`star_systems.location`; v3→v4 added `sectors`' galaxy-frame placement
columns; v4→v5 dropped every pre-rendered `table_*`/`binary_table_*`
display-string column (superseded by computing display formatting on
demand from the underlying data columns, e.g. `html/lib/tabledisplay.py`);
v6/v7 gave every galaxy-placed sector exact vertices (`sector_vertices`,
built from an exact local spherical Voronoi tessellation among its
same-shell neighbors — see `stellarObjects/galaxyGeometry.py`); v8 added
the galaxy-wide density "skeleton" (`galaxy_shape`/`galaxy_layer`,
built by `generate.py plan` — see `stellarObjects/galaxyDensity.py`/
`galaxySkeleton.py`) plus a `UNIQUE (shell_index, shell_slot_index)`
constraint on `sectors`, turning a concurrent lazy-generation race
(`generate.ensure_sector_generated`) into a recoverable `IntegrityError`
instead of a silent duplicate row; v9 added orbital motion —
`orbital_inclination_deg`/`orbital_ascending_node_deg`/
`orbital_phase_deg`/`rotation_period_hours` on `planets`/`moons`, plus the
`orbit_simulation_state` singleton row `updateOrbits.py` uses to track
elapsed time between runs (see that table's own section above); v10 added
`galactic_orbital_speed_kms`/`galactic_orbital_period_gy` to `stars` (and
their `binary_galactic_orbital_*` counterparts on `star_systems`) — a star
system's circular orbital speed/period around the galactic center, from a
simple rotation-curve model (see `stellarObjects/physical_constants.py`'s
`GALACTIC_ROTATION_FLAT_VELOCITY_KMS` comment); v11 added
`position_x_km`/`_y_km`/`_z_km`/`orbital_speed_kms` to `planets`/`moons` —
each body's Cartesian position relative to its orbital anchor (the star,
or a binary's combined center, for a planet; the parent planet for a
moon), derived from `distance_km` and the v9 orbital-motion columns (see
`stellarObjects/utils.py`'s `orbital_position_au`); v12 added
`planets`/`moons.min_update_interval_years` — a floating-point update
guard, not a narrative stat: the shortest `elapsed_years` worth calling
`_db.advance_orbital_phases` for, below which the phase delta added is
smaller than `orbital_phase_deg`'s own IEEE 754 double-precision
resolution and so is guaranteed to round back to the exact value already
stored (see `stellarObjects/utils.py`'s `minimum_update_interval_years`).
Scoped to `planets`/`moons` only at v12 — `stars`' galactic-orbit values
were, at that point, fixed forever at generation time, with no periodic
update mechanism to guard; v13 changed that (see below). v13 added star
motion: `stars.galactic_orbital_phase_deg`/`galactic_min_update_interval_years`
(the same phase/guard pair planets/moons have, now advanced by
`_db.advance_orbital_phases` too, based on `galactic_orbital_period_gy`),
`star_systems`' matching `binary_galactic_orbital_phase_deg`/
`binary_galactic_min_update_interval_years` (always identical to both
constituent stars' own values — a binary pair's negligible AU-scale
separation next to its light-year-scale galactic orbit means the pair
moves around the galaxy together, not independently — see
`StarSystem.__init__`), plus the binary pair's own *mutual* orbit around
each other (entirely separate from, and vastly faster than, the galactic
orbit above): `binary_mutual_orbital_period_years`/`_speed_kms`
(`planetPhysics.calculate_orbital_period_years`/
`utils.circular_orbital_speed_kms`, the same Kepler/circular-orbit formulas
a planet's orbit around its star already uses, applied to the pair's
`binary_separation_km`/`binary_effective_mass_kg`), `_inclination_deg`/
`_ascending_node_deg`/`_phase_deg` (the same `utils.orbital_position_au`
orbital-element convention planets/moons use, drawn from the full
`[0, 180)`/`[0, 360)` range with no small-tilt bias — a binary's mutual
orbital plane has no protoplanetary-disk reason to prefer any alignment,
unlike a planet's), and `binary_mutual_min_update_interval_years`. v14
added `binary_mutual_position_x_km`/`_y_km`/`_z_km` — the secondary's
Cartesian position relative to the primary, derived from
`binary_separation_km` and the v13 `binary_mutual_orbital_*` columns via
`utils.orbital_position_au`, the same "position relative to whatever this
orbit is around" convention `planets`/`moons.position_x/y/z_km` already
use (see v11 above) — recomputed by `_db.advance_orbital_phases` in
lockstep every time `binary_mutual_orbital_phase_deg` advances.

v15 added S-type (wide) binary support, planetGen's second real binary
configuration alongside the P-type/close pair every `binary_*` column
above already covers (see `stellarObjects/wideBinary.py`'s module
docstring for the physics). `star_systems` gains `binary_configuration`
(`'close'` | `'wide'` | NULL for a single star — the authoritative
discriminator going forward; `is_binary` is kept and now means "this
system has two stars", true for either configuration) plus
`binary_eccentricity`/`binary_periapsis_km`/`binary_apoapsis_km` (0/
`binary_separation_km` for a `'close'` pair, real values for a `'wide'`
one — see `utils.holman_wiegert_critical_semimajor_axis`). The existing
`binary_separation_km`/`binary_mutual_orbital_*`/
`binary_mutual_position_*_km` columns are reused UNCHANGED for a `'wide'`
pair's own mutual orbit. NULL for a `'wide'` pair, unlike a `'close'` one
(no merged effective star exists to describe): the pre-existing
`binary_type` column (note this is a different thing from the new
`binary_configuration` — `binary_type` holds the close pair's merged
spectral-summary string, e.g. `"Binary (G/K)"`, kept under its original
name rather than repurposed), `binary_temperature_k`, `binary_radius_km`,
`binary_effective_mass_kg`, `binary_effective_luminosity_w`,
`binary_age_gy`, `binary_lifespan_gy`,
`binary_habitable_zone_inner_km`/`_outer_km`, `binary_system_perimeter_km`,
`binary_heliosphere_radius_km`, and the four `binary_galactic_orbital_*`
columns — a `'wide'` pair's two stars already carry their own galactic
orbit and habitable-zone columns individually on their own `stars` rows.
`stars` gains `wide_binary_a_crit_km` — this star's own Holman & Wiegert
(1999) critical semi-major axis, the maximum orbit distance that stays
long-term stable given its companion's perturbation; NULL except for a
constituent of a `'wide'` pair. `asteroid_belts` gains `star_id`, the same
"which specific star this orbits" column `planets`/`moons` already have
(see that column's own row below) — a `'wide'` pair's two stars can now
each have their own asteroid belts, not just their own planets.

v16 added six tables for `generate.py phenomenon`'s separate, rarer exotic-
phenomenon generation mode (`stellarObjects/compactRemnant.py`/
`nebulaData.py`/`supernovaRemnantData.py`/`roguePlanetData.py`) — never
produced by `generate.py system`/`generate.py sector`'s normal generation odds, so no
pre-existing table's shape changes. `black_holes`/`neutron_stars` are
satellite tables extending a `stars` row (nullable `star_id`) when a
compact remnant anchors a full `StarSystem` (`generate.py
phenomenon --anchor-system` — the owning `stars.yerkes_class` is then the literal
marker `'BH'`/`'NS'` rather than a real Yerkes class); `star_id` is NULL
for a remnant generated standalone. `nebulae`/`supernova_remnants`/
`rogue_planets`/`interstellar_comets` are always standalone, with a
nullable `sector_id` reserved for a future sector-context encounter
(unused by `generate.py phenomenon` today). `supernova_remnants` references at
most one of `black_holes`/`neutron_stars` (`compact_remnant_kind`
discriminator) for a core-collapse progenitor whose collapsed core is
still detectable — always both NULL for a Type Ia progenitor, which
leaves nothing behind. `interstellar_comets` has its own
`interstellar_comet_composition` child table, mirroring
`asteroid_belt_composition`'s per-component breakdown minus a
concentration level.

v17 added galactic-orbital motion for every standalone exotic phenomenon
— a black hole/neutron star with no owning system, a nebula, a supernova
remnant, a rogue planet, an interstellar comet, or a standalone asteroid
field is still gravitationally part of the galaxy despite being bound to
no star, so each now carries the same `galactic_orbital_speed_kms`/
`_period_gy`/`_phase_deg`/`_min_update_interval_years` quartet a lone star
has (new shared `utils.generate_galactic_orbit_fields`/
`format_galactic_orbit` helpers). `black_holes`/`neutron_stars` gain these
four columns nullable (populated only when `star_id IS NULL` — an
anchored remnant's motion already lives on its own `stars` row);
`nebulae`/`supernova_remnants`/`rogue_planets`/`interstellar_comets` gain
them `NOT NULL` (always standalone, always populated). Also added a
seventh phenomenon type, standalone asteroid fields (new `asteroid_fields`
table plus child `asteroid_field_composition`, mirroring
`asteroid_belts`/`asteroid_belt_composition`'s shape) — physically the
same object as `asteroid_belts` (density + mineral composition, generated
via the same shared `asteroidData.generate_asteroid_composition`/
`format_composition_summary` helpers) but standalone, drifting in open
space rather than orbiting a star.

v18 gave `nebulae`/`asteroid_fields` an actual galaxy placement — the
first real use of `sector_id`, reserved-but-unused since v16/v17. A
nebula/asteroid field is frequently far larger than a single sector's
cube (an emission nebula can span up to 200 ly; the default sector edge
is 4 pc, ~13 ly), so unlike `star_systems.position_x/y/z_mpc` (relative to one
owning sector's own center) these are placed directly in the same
galaxy-frame Cartesian space `sectors.center_x/y/z_pc` already uses — a
sphere (new `center_x/y/z_pc`/`galactic_radius_pc` columns, NULL together)
that may overlap zero, one, or several sectors' cubes, not a single
sector-relative offset. `sector_id` is now actually populated
(`generate.py phenomenon --sector-id`) with the *nearest* already-generated
sector to the phenomenon's own center — a convenience "home" link for
browsing (`ON DELETE SET NULL` now, not `CASCADE`: deleting that sector
shouldn't delete a phenomenon merely linked to it), not the authoritative
geometry (`queryDb.phenomena_near_sector` finds every phenomenon whose
sphere overlaps a given sector's cube by real distance, regardless of
which sector it's linked to). `stellarObjects._db.compute_phenomenon_placement`
picks the center: the given sector's own stored galaxy position plus a
uniform random jitter within that sector's own cube half-extent. The
null-together CHECK on each table is named explicitly
(`chk_nebulae_placement`/`chk_asteroid_fields_placement`), unlike
`sectors`' own identical v4 CHECK — MySQL auto-names an anonymous CHECK
opaquely and refuses `DROP COLUMN` on a column it still references, so an
explicit name is what lets `_migrate_v17_to_v18` add (and, for a test
simulating an older database, drop) the exact same constraint by name.

v20 added a proper two-body (barycentric) trajectory treatment for every
orbital pair where the orbited body isn't overwhelmingly more massive than
what orbits it: binary stars (a secondary is sampled at 0.1-0.8x the
primary's mass), star↔planet, and planet↔moon (a moon can reach 1/10 its
parent planet's mass). Every existing "relative position" column
(`planets`/`moons.position_x/y/z_km`,
`star_systems.binary_mutual_position_x/y/z_km`) is unchanged — still the
true separation a large amount of existing physics (insolation, Hill
sphere, tidal locking) depends on. New columns instead add the ORBITED
body's own small "reflex offset"/"wobble" away from its nominal fixed
point (see `utils.calculate_reflex_offset`): `stars`/`planets` each gain
`reflex_offset_x/y/z_km` (a star's from the planets orbiting it directly,
a planet's from its own moons — NULL/0 with none); `star_systems` gains
`binary_primary_position_x/y/z_km`/`binary_secondary_position_x/y/z_km`
(each binary member's own offset from the pair's barycenter —
`secondary_position = primary_position + binary_mutual_position` always
holds, but both are stored explicitly), the constant
`binary_secondary_mass_fraction`, and `binary_planetary_wobble_x/y/z_km`
(a 'close' pair's combined pull from its own circumbinary planets, which
orbit the merged proxy rather than either individual star — modeled as
one shared wobble rather than split between primary/secondary, since that
would need a real 3+-body solve). These are new columns on already-
existing tables (`star_systems`, `stars`, `planets`), so
`_migrate_v19_to_v20` needs real `ALTER TABLE` steps — and, unlike v17's
own migration, every value is fully derivable from data already on
existing rows, so it backfills real values rather than leaving them NULL.

v21 gave every generated sector its own realistically sparse population of
exotic phenomena (`generate.generate_sector_phenomena`, sampled from
`program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM`, replaced in v37 by
per-star research densities) instead of leaving
all seven types reachable only through `generate.py phenomenon`'s separate,
on-demand CLI. `nebulae`/`supernova_remnants`/`rogue_planets`/
`interstellar_comets`/`asteroid_fields` needed no schema change at all —
`insert_sector` simply started populating their pre-existing `sector_id`
column (and, for the two galaxy-placeable types, `center_x/y/z_pc`) for
real. `black_holes`/`neutron_stars` gained the same `sector_id`/
`center_x/y/z_pc`/`galactic_radius_pc` shape `nebulae`/`asteroid_fields`
got in v18, since a standalone compact remnant is a real, stellar-mass
gravitating body and needed to be placeable too — but unlike every other
phenomenon type, its own galaxy-frame position is usually the *exact* spot
`SpaceSector.add_phenomenon`'s Hill-sphere-aware placement computed at
generation time (never within a neighboring star system's or another
compact remnant's own Hill sphere), converted from a sector-relative
offset rather than independently re-randomized.

The SQLite-specific machinery that once converted an existing database
between these versions in place (gzip-compressed file backups, a
`_migrate_vN_to_vN+1` function per version) was removed during the MySQL
port (TODO.md Phase 5) on the assumption that every MySQL database this
project creates starts at the current schema directly, with no "upgrade
an older MySQL database" case to handle — true until v9, whose
`_migrate_v8_to_v9` (`stellarObjects/_db.py`) is the first real migration
function of the MySQL era, reviving the same per-version-step pattern
(minus the file backups, which made no sense for a live database anyway)
for a database created under the v8 schema; `_migrate_v9_to_v10` follows
the same pattern for the v10 galactic-orbit columns, `_migrate_v10_to_v11`
for the v11 planet/moon position columns, `_migrate_v11_to_v12` for
the v12 floating-point update-guard column, `_migrate_v12_to_v13` for the
v13 star-motion/binary-mutual-orbit columns, `_migrate_v13_to_v14` for
the v14 binary-mutual-orbit-position columns, `_migrate_v14_to_v15`
for the v15 S-type/wide-binary columns, `_migrate_v15_to_v16` for
v16's exotic-phenomenon tables (needing no `ALTER TABLE` at all — six
brand-new tables are already created by `_ensure_schema`'s
`CREATE TABLE IF NOT EXISTS` regardless of the database's recorded
version; this step exists purely to keep the `schema_migrations`
bookkeeping counter itself accurate), `_migrate_v16_to_v17` for v17's
galactic-orbital-motion columns (real `ALTER TABLE` steps this time, since
v16's tables already existed with a fixed shape), `_migrate_v17_to_v18`
for v18's galaxy-frame placement columns on `nebulae`/`asteroid_fields`
(also real `ALTER TABLE` steps, left NULL on every pre-existing row --
there's no way to recover a legacy phenomenon's intended galaxy placement
after the fact, the same situation `sectors.center_x/y/z_pc` is in for a
`_migrate_v3_to_v4`-migrated sector), `_migrate_v18_to_v19` for v19's
new `comets`/`comet_composition` tables (needing no `ALTER TABLE` either,
for the same "brand-new tables, bookkeeping-only step" reason as
`_migrate_v15_to_v16`), `_migrate_v19_to_v20` for v20's binary/star/
planet trajectory columns (also real `ALTER TABLE` steps, on
`star_systems`/`stars`/`planets`), `_migrate_v20_to_v21` for v21's
sector-placement columns on `black_holes`/`neutron_stars` (also real
`ALTER TABLE` steps, left NULL on every pre-existing row — same
"nothing to recover" situation v18's migration is in), and
`_migrate_v21_to_v22` for v22's search-facing indexes on `sectors`/
`star_systems`/`stars`/`planets`/`moons`.`name` and the facet/filter
columns `GET /api/search` groups/filters by (`ALTER TABLE ... ADD KEY`
steps only — no new columns, nothing to backfill), and so on, one step per
version, through `_migrate_v50_to_v51`. `migrate_database`
applies whatever steps are needed to reach `SCHEMA_VERSION`, one call
`migrateDb.py` wraps as a CLI (also run automatically by
`install.sh`/`update.sh` on every deploy). A SQLite database from before
the MySQL port can't be brought in: the one-time import script that used
to be here (`src/migrateSqliteToMysql.py`) was retired (TEST.61), because
it only accepted a file already at the current `SCHEMA_VERSION`, which no
SQLite database ever was (SQLite stopped at v5, and the MySQL migration
steps start at v8). Generate the galaxy again instead.

**v19 to v26, in brief.** v19 added star-bound comets (`comets`,
`comet_composition`; see that table below). v22 added the search-facing
indexes above. v23 wired up wiki publishing (`sectors.wiki_url`, and
`star_systems.wikijs_url`/`mediawiki_url` for existing databases). v24
added the name-uniqueness registries (`sector_name_registry`,
`system_name_registry` and a planet/moon `body_name_registry` that v34
later dropped), so no two sectors or systems share a display name. v25
and v26 added composite `(center_x_pc, center_y_pc, center_z_pc)` indexes
on `sectors` and the placeable phenomenon tables, so bounding-box queries
range-scan instead of reading every row.

**Row timestamps (v27).** `sectors`, `star_systems` and the seven
exotic-phenomenon tables each carry `created_at` and `modified_at`
(millisecond precision), with an index on `modified_at`. MySQL bumps
`modified_at` itself on any `UPDATE` that changes the row. Child rows
(stars, planets, moons, belts, comets) have no timestamps of their own;
changing one bumps its system's `modified_at` (`_db.touch_star_system`),
and deleting a system bumps its sector's (`_db.touch_sector`). Orbit
ticks from `updateOrbits.py` deliberately don't count as a change; see
`orbit_simulation_state.last_updated_at` for those. `_migrate_v26_to_v27`
adds the columns with `ALGORITHM=INSTANT` where the server supports it
and builds the indexes online, so it's safe on a large, live database.
It then backfills in small batches: each system's `modified_at` starts at
its own `created_at`, and each sector takes its oldest system's
`created_at` for both columns. Phenomena and empty sectors have no record
of when they were made, so they keep the migration's own time. The full
reasoning is in `schema.sql`'s "v27" header note.

**Phenomenon placement (v28).** `supernova_remnants`, `rogue_planets`
and `interstellar_comets` gain the same `center_x_pc`/`center_y_pc`/
`center_z_pc`/`galactic_radius_pc` columns the other four phenomenon
tables have, so every phenomenon type can appear on the Sector Map and be
a NAV endpoint. A supernova remnant's embedded black hole or neutron star
takes the remnant's own sector and center. `_migrate_v27_to_v28` adds the
columns and indexes online like v27, then gives each existing row that is
linked to a galaxy-placed sector a random point inside that sector's
cube (its real generated position was never saved), seeded by the row's
id. Rows with no placed sector stay unplaced. The full reasoning is in
`schema.sql`'s "v28" header note.

**Page text rendered on demand (v29).** `_migrate_v28_to_v29` drops
`star_systems.wikitext_content`/`markdown_content`; see "Rendered wiki
text and URLs" below.

**Data cleanup (v30).** No shape change. `_migrate_v29_to_v30` clears the
atmosphere-only columns on planets and moons reclassified into an
airless class (and on `atmosphere = 'None'` rows), and raises any surface
temperature below the cosmic microwave background (2.725 K) to it.

**Quasars (v31).** A new `quasars` table holds a galaxy's active nucleus
(`quasarData.Quasar`). A quasar is only ever generated at the galactic
center, by the core sector at ring 0, layer 0, slot 0 (`generate.add_galactic_nucleus`,
rolled against `program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE`), so a
placed row always sits at (0, 0, 0) and a galaxy has at most one. It has
the usual placement columns and row timestamps, but no galactic-orbit
columns: it is the point everything else orbits. A brand-new table needs
no `ALTER TABLE`, so `_migrate_v30_to_v31` only records the version.

**Cylindrical sector grid (v32).** Galaxy-placed sectors moved from
spherical shells to rings, layers and slots: `sectors.shell_index`/
`shell_slot_index` became `ring_index`/`layer_index`/`ring_slot_index`,
`sector_vertices` was dropped, `galaxy_shell_band` became
`galaxy_ring_band`, and `galaxy_shape.outer_shell_index` became
`outer_ring_index`. Old addresses have no matching cell, so
`_migrate_v31_to_v32` deletes every galaxy-placed sector with its systems
and phenomena (never-placed sectors stay). See `schema.sql`'s "v32"
header note.

**One sector standard (v33).** The sector edge became a whole number of
parsecs, 4 pc (~13.05 ly) instead of 11.5 ly; ring `i` holds
`round(2*pi*(i + 1/2))` slots on every layer instead of a multiple of 4;
and `galaxy_ring_band` became `galaxy_layer`, one row per layer with the
last ring it reaches, plus `galaxy_column`, one row per ring with the
highest and lowest layer it reaches. Almost every address moves, so `_migrate_v32_to_v33`
again deletes every galaxy-placed sector with its systems and phenomena,
and rebuilds the skeleton from the stored shape at 4 pc (rewriting
`galaxy_shape.edge_pc`, `expected_system_count_at_density_1` and
`outer_ring_index`). See `schema.sql`'s "v33" header note and
`docs/design/galaxy-coordinate-system.md`.

**Planets and moons named from their system (v34).** Only the system
draws a generated name. Planets are `<star> I`, `<star> II` in orbit
order, moons add a letter (`<star> IIa`), and a binary's stars are
`<system> <word>`. The system name is already unique, so these are too,
and `body_name_registry` (v24's planet/moon registry, which gave
colliding names a companion suffix like `"Kin"`) is dropped by
`_migrate_v33_to_v34`. Existing rows keep their names until regenerated.
See `stellarObjects/bodyNames.py`.

**Hybrid master-wedge slots (v35).** Ring `i` now holds the multiple of
its master wedge count (3 at the center, doubling once each master wedge
would hold 8 slots) nearest `2*pi*(i + 1/2)`: 3, 9, 15, 21, 27, 36, ...
(`galaxyGeometry.ring_sector_count`/`ring_master_count`). Slot boundaries
then line up on the master lines from the center out, which the Galaxy
Map's mega-blocks cut on. No column changes, but a stored
`ring_slot_index` changes meaning in all but 15 of 3,856 rings, so
`_migrate_v34_to_v35` deletes every sector in a changed ring with its
systems and phenomena; the skeleton (`galaxy_layer`, `galaxy_column`) is
keyed by ring and layer and stays. See `schema.sql`'s "v35" header note
and `docs/design/galaxy-coordinate-system.md`.

**A black hole at every galaxy's center (v36).** When the nucleus roll
finds no quasar, `generate.add_galactic_nucleus` places a quiescent
supermassive black hole there instead, so every galaxy has one.
`black_holes.mass_class` records whether a row is stellar-mass,
intermediate-mass or supermassive; `_migrate_v35_to_v36` fills it from each
row's mass.

**Id blocks (v45, PERF.13).** `id_blocks` holds one row per table in
`_db.ID_BLOCK_TABLES` (sectors, star systems, their configs, stars,
planets, moons, belts, comets and every phenomenon table): `table_name`
and the `next_id` no writer has reserved yet. A new row in one of those
tables gets its `id` from this process's current block
(`_db._allocate_id`), not from AUTO_INCREMENT, so a sector's rows can be
written many at a time (`Connection.batched`) with each child already
knowing its parent's id. A block is reserved on its own autocommitted
connection with `UPDATE ... SET next_id = LAST_INSERT_ID(GREATEST(next_id,
MAX(id) + 1) + n)`, 64 ids at first and doubling up to 4,096, so writers
never wait on each other's transactions; unused ids and rolled-back
sectors only leave gaps. `_migrate_v44_to_v45` creates it empty, and a
database without it yet simply falls back to AUTO_INCREMENT.

**Full-text name search (v46, PERF.16).** `sectors`, `star_systems`,
`stars`, `planets` and `moons` each carry a `FULLTEXT KEY
ft_<table>_name (name)` beside their plain `idx_<table>_name`. Search
(`queryDb._name_match`) matches whole words in boolean mode, so a
search never scans every row with `LIKE '%term%'`; words below the
server's `innodb_ft_min_token_size` or on its stopword list fall back
to a whole-word `REGEXP` on the rows the other words narrowed.
`_migrate_v45_to_v46` adds the indexes in place (`ALGORITHM=INPLACE,
LOCK=SHARED`): reads keep working, but writes to those tables wait
until each index is built, which on a large galaxy can take minutes.

**Rogue planet classes (v47, GEN.8).** `rogue_planets.planet_class` is a
rogue's `PLANET_CLASSES` letter, drawn from the classes whose `"r"` flag
says they can be rogue (C, D, I, J and T) and that fit its type, radius
and mass. NULL for a brown dwarf. `_migrate_v46_to_v47` gives every
stored rogue its most probable fitting class.

**Rogue planet surface conditions (v48).** `rogue_planets` gains
`age_gy`, `internal_heat_flux_w_m2`, `effective_temperature_k`,
`surface_regime`, `surface_temperature_k`, `surface_pressure_pa`,
`ice_shell_thickness_km`, `ocean_depth_km` and `has_liquid_water`, from
`rogueSurface.rogue_surface_conditions` (see
`docs/design/rogue-planet-surface.md`). `_migrate_v47_to_v48` fills them
for every stored rogue, seeded by its name, and resets
`has_internal_heat` to the computed answer.

**Population and politics (v44).** Filled by `generate.py population`
(or `--population` on a `sector`/`galaxy` run, or the optional question
the install and update scripts ask) from what is already stored; see
docs/design/population-and-politics.md. `species` holds one dominant
species per life world (a planet whose evolutionary milestone is
multicellularity or a technological civilization): a galaxy-unique
`name`, its homeworld (`homeworld_planet_id`, `star_system_id`, both
cascading), `life_chemical`, `life_stage`, `build`/`climate`/`size`, and
for a civilization its `civilization_age_years`, `era` and `spacefaring`.
`polities` holds one government per spacefaring species (`species_id`
and `capital_system_id` cascade) with its `government`, map `color` and
territory `reach_ly`. `system_owners` maps each owned system to its
polity and `distance_ly` from the capital, rebuilt from scratch by every
pass. `population_state` keeps the highest planet id scanned.
`_migrate_v43_to_v44` creates all four, empty.

**Bright-star pre-placement (v43).** `bright_stars` holds every star at
least `galaxy_shape.bright_star_min_luminosity_sol` bright, generated and
placed galaxy-wide by `generate.py plan` before any sector is filled: its
cell (`ring_index`, `layer_index`, `ring_slot_index`), galaxy-frame
position in milliparsecs, population (young, intermediate, old or
bulge), the finished star's properties (including `initial_mass_sol`,
which a companion's mass is drawn from, and `phase_end_age_gy`) and a
seed for its later system. Rows carry no name until their sector is
filled, which builds a system around each and sets `star_system_id`
(`ON DELETE SET NULL`). `galaxy_shape.bright_star_seed` and the threshold
record the scatter; a fill reads them, not the constant. A plan re-run
truncates the table (`_db.clear_bright_stars`). `_migrate_v42_to_v43`
adds both, empty.

**Bright-star backfill per sector block (v49, GEN.23).** Generating any
galaxy sector first backfills the stars around it
(`generate.backfill_bright_stars`): every sector block (a level-3 block
of `galaxyDrill`, 3 rings by 3 layers by its wedge's slots) with a sector
within 100 ly gets every star from 100 L_sun up to the level it already
holds, in its sectors that aren't filled yet, as new `bright_stars` rows.
`bright_star_blocks` keeps, per block (`block_ring`, `block_wedge`,
`block_slab`), the dimmest luminosity it now goes down to
(`min_luminosity_sol`), so a block already at 100 L_sun is skipped and no
star is drawn twice. A block with no row is at the galaxy scatter's
threshold. A sector's fill caps its own dim stars at its block's level
(`_db.bright_star_fill_level`). A backfill takes the block's row lock
(`INSERT ... ON DUPLICATE KEY UPDATE`, then `SELECT ... FOR UPDATE`,
level NULL until it commits; parallel workers queue on it), so two generators never draw the same block. A plan re-run
truncates it with `bright_stars`. `_migrate_v48_to_v49` creates it empty.

**The galaxy seed (v51).** `_migrate_v50_to_v51` adds
`galaxy_shape.galaxy_seed` (`BINARY(16)`: half the bytes of a
32-character hex string and no letter case to normalize), NULL on an
existing galaxy, which has no seed and can't be given one: GEN.39 starts
from a wiped galaxy.

**Matching a new database (v50).** `_migrate_v49_to_v50` adds
`system_configs.comets` and `.wide_binary`, and brings a database
migrated from an old version to exactly the shape `schema.sql` gives a
new one: it drops the placeholder DEFAULTs older steps left on NOT NULL
columns and makes `nebulae`/`asteroid_fields`' `sector_id` keys ON
DELETE SET NULL. Every migration step can be re-run after a crash
(`_MigrationConnection` skips the clauses already applied).
`tests/test_db_old_schemas.py` migrates every released schema and
compares it with a new database.

**Facilities (v42).** `facilities` holds starbases, colonies and outposts,
each on one host named by `host_type`: a star, planet, moon, asteroid belt,
asteroid field, or open space in a sector (`sector_id` plus a galaxy-frame
center). `star_system_id` is set for every in-system host. An orbital
facility stores its circular orbit (`orbit_distance_km`,
`orbit_period_years`, `orbital_speed_kms`, `orbit_phase_deg`). Every host
key cascades. `program_constants.FACILITY_RULES` decides which kinds go
where, checked by `_db.add_facility` (MySQL refuses a CHECK on a cascading
column). `_migrate_v41_to_v42` creates the table.

**Octants and nearest systems (v41).** Every placeable phenomenon
table gains `quadrant`, the sector octant its center sits in (the labels
`star_systems.quadrant` uses, in the `sector_orientation` frame).
`nearest_systems` holds up to 3 rows per placed system or phenomenon:
`object_table`/`object_id`, `neighbor_rank`, `neighbor_system_id` and
`distance_pc`, searched across sector boundaries out to
`_db.NEAREST_SYSTEMS_SEARCH_PC` (4 pc). `sector_id` (the object's
sector) and `star_system_id` (set for a system) cascade deletes.
`insert_sector` fills a new sector's rows and merges its systems into
its neighbors' lists; `refresh_nearest_systems` recomputes whole
sectors. `_migrate_v40_to_v41` adds both and fills them.

**Names (v40).** Black holes, neutron stars, nebulae, supernova
remnants, rogue planets and quasars draw their names through
`system_name_registry`, the same collision rules star systems use, so
none shares a name with a system, a sector or each other.
`system_name_registry.first_star_system_id` is nullable and
`first_object_table`/`first_object_id` name a phenomenon holder. A
remnant's core is named `"<remnant> Core"` and follows its remnant; an
anchored black hole or neutron star shares its system's name. Comets and
asteroid fields get designations instead: `P/<host>-<n>` (periodic,
under 200 years) or `C/<host>-<n>` for a star-bound comet, which follow
a rename of their star, `I/<sector>-<n>` for an interstellar comet and
`AF <field_class>-<sector>-<nn>` for an asteroid field, where
`<sector>` is the sector's grid designation (or its name off the grid).
Each named phenomenon table indexes `name`, and `_db.name_in_use`
searches them. `_migrate_v39_to_v40` renames and registers existing rows.

**Containment (v39).** `star_systems`, `rogue_planets`,
`interstellar_comets`, `black_holes`, `neutron_stars`, `asteroid_fields`
and `nebulae` get `inside_nebula_id` and `inside_remnant_id`, foreign keys
(`ON DELETE SET NULL`) to the innermost nebula or supernova remnant whose
sphere holds the object; at most one is set. `_db.refresh_containment`
sets them by a 3D distance test when a sector is generated and, for
every sector it reaches, when a nebula or remnant is placed.
`_migrate_v38_to_v39` adds the columns and fills them. A nebula only
nests inside a larger cloud. `queryDb.containing_cloud` names the cloud,
and `sector_detail` returns it per system as `inside`. `schema.sql` turns
foreign key checks off while it runs, since these tables are created
before `nebulae`.

**Letter classes (v38).** Nebulae (A-Q) and supernova remnants (R-W)
get a class from `program_constants.NEBULA_CLASSES` with their contents
(`dominant_species`, `density_cm3`, `temperature_k`, `extinction_av`),
and asteroid fields a `field_class` such as `C3` (letter from composition
and density, digit from size). `_migrate_v37_to_v38` infers classes for
existing rows (old asteroid fields become the `mixed` family) and fills
each class's typical contents. See
`docs/design/nebula-and-asteroid-field-classes.md`.

**Research-based interstellar rates (v37).** Phenomena are now drawn per
star from Boss's research densities (`program_constants.
PHENOMENON_DENSITY_PC3`, `docs/design/interstellar-object-rates.md`).
`rogue_planets.mass_bin` records a rogue's mass bin (`terrestrial`,
`sub-neptune`, `saturn`, `jupiter`) or `brown-dwarf` (free-floating brown
dwarfs share the table), and `star_systems.runaway_class`/
`runaway_speed_kms` flag runaway and hypervelocity stars.
`_migrate_v36_to_v37` fills `mass_bin` from each existing row's mass.

**This versioning is independent of the control schema's own.** Admin
logins/sessions/API keys/the write-action audit log live in a separate
MySQL schema entirely (`stellarObjects/control_schema.sql`,
`control_schema_migrations`, currently version 7) — see "The control
schema" below. `SCHEMA_VERSION`/`schema_migrations` above only ever
describe the per-galaxy content schema this whole document is otherwise
about.

## The control schema

A second, deployment-global MySQL schema (`PLANETGEN_CONTROL_DATABASE`,
default `planetgen_control`) holds everything about *who can administer
this deployment*, separate from every per-galaxy content schema this
document otherwise describes — see `stellarObjects/control_schema.sql`'s
header comment for the full rationale (in short: a deployment can host
several galaxy databases sharing one MySQL server, and admin identities
describe the deployment, not any one galaxy, so they aren't duplicated
into each content schema's `schema.sql`).

Thirteen tables, versioned independently via `control_schema_migrations`
(currently version 7, mirroring `schema_migrations`'s own shape; v2 added
`login_throttle`, v3 `admin_devices`, v4 `admin_totp` and
`admin_recovery_codes`, v5 the work queue's three tables, v6
`generation_stats` and `generation_size`, and v7 the work queue's job
tree and pause columns; v2 to v6 are new tables, which `CREATE TABLE IF
NOT EXISTS` adds to an older schema on the next `migrateDb.py` run, and
v7 adds columns, see below):

- **`admin_users`** — one row per admin (`username`, `password_hash`,
  `must_change_credentials`). No roles/permissions column — every admin
  has the same full access (see `docs/TODO.md`; this project deliberately
  has no general user-accounts system, just a handful of admins).
- **`admin_sessions`** — web-UI login sessions (`token_hash`, a SHA-256
  digest of the actual cookie value — never the raw token itself;
  `expires_at`, a fixed lifetime set at creation, no sliding renewal).
- **`admin_api_keys`** — API keys for programmatic callers (`key_hash`,
  same "hash only, never the raw key" treatment; `revoked_at` rather than
  a hard delete, so a revoked key's history stays visible).
- **`admin_audit_log`** — one row per write/admin action (`admin_user_id`
  + a denormalized `admin_username` snapshot, `action`, `target`,
  `detail`, `created_at`) — written by `html/api/routes.py`'s write
  routes (`html/api/authz.audit`) after each one actually succeeds,
  plus one per refused sign-in (`login.failed`, `login.locked`,
  `password.failed`, with `target` `ip:<address>`; those are deleted
  after 90 days).
- **`login_throttle`** (v2, SEC.1/SEC.21) — failed-login counts and
  lockouts: `scope` (`ip`, or `user`) and `subject` (the address, an IPv6
  one by its /64, or the case-folded username) as the key, then
  `failures`, `level` (lockouts so far, which double the next), and
  `locked_until`/`last_failure_at`/`last_lockout_at` in Unix seconds.
  Read and written by `stellarObjects/loginThrottle.py`; idle rows are
  deleted after a week.
- **`admin_devices`** (v3, SEC.22) — trusted browsers: `admin_user_id`,
  `token_hash` (SHA-256 of the `pg_admin_device` cookie, never the raw
  value), `expires_at` (90 days after the login that made it),
  `last_used_at`. A login from a browser holding a valid one for that
  username skips the per-username lock. Deleted on a credentials change,
  with `src/loginLockouts.py --forget-devices <user>`, and (expired ones)
  when the admin's next device is made.
- **`admin_totp`** (v4, SEC.26) — one row per admin who has set up an
  authenticator app: `secret` (the base32 key itself, as sensitive as a
  password, since a code can only be checked against the key),
  `enabled_at` (NULL until the first code confirms the setup), and
  `last_step` (the newest 30-second step a code was accepted for, so no
  code works twice).
- **`admin_recovery_codes`** (v4) — ten single-use codes per admin with
  two-factor sign-in on: `code_hash` (SHA-256), `used_at`.
- **`work_jobs`**, **`work_tasks`**, **`work_lease`** (v5, PERF.8) — the
  generation work queue (`stellarObjects/workQueue.py`, see
  [`cli.md`](cli.md#parallel-generation)). `work_jobs` is one row per
  run that used worker processes (`title`, `holder` host:pid:token,
  `state` waiting/running/done/failed/cancelled, `workers`, task
  counts, timings, `heartbeat_at`); `work_tasks` one per sector it
  handed out (`kind`, `task_key` such as "ring,layer,slot", `state`,
  the worker's `seconds`, a short JSON `result` or the `error`);
  `work_lease` the single row naming the run whose workers are using
  the machine, refreshed every 5 s. A lease or job not refreshed for
  30 s belongs to a run that died: the next run takes the lease and
  marks that run's unfinished tasks cancelled. Jobs older than a week
  are deleted with their tasks.
- **Job tree columns** (v7, ADM.12) — every job is a tree. `work_jobs`
  gained `parent_id` (the node above, deleting a node deletes its
  subtree), `root_id` (the top of its tree), `kind` (`web-job` for a
  Generate page job, `step` for one of its steps, a `generate.py`
  command such as `plan` or `galaxy` for a run, a phase such as
  `skeleton`, `bright-stars` or `population`, and `queue` for a work
  queue, whose leaves are its `work_tasks`), `seconds` (its own
  duration), `tasks_total` (tasks the queue said it would queue, for
  the ETA), `web_job_id`, `database_name`, `argv` (the run's command
  line as JSON, without the `--mysql-*` and `--debug` options) and
  `control` (`pause` or `cancel`, asked from the admin queue page).
  Every node, not only those with workers, is recorded and timed; a
  parent's totals are added up from its children when the tree is read
  (`workQueue.load_tree`). A node whose run stopped refreshing it for
  30 s reads as interrupted. Finished trees older than a week are
  deleted whole. `work_lease` gained `paused`, `paused_by` and
  `paused_at` for "Pause the queue" (ADM.10). The columns are added to
  an older control schema by `update.sh` (`_db._add_control_columns`).
- **`generation_stats`**, **`generation_size`** (v6, PERF.3, PERF.10) —
  how fast this server generates and how big a galaxy gets
  (`stellarObjects/generationStats.py`, see
  [`cli.md`](cli.md#size-and-time-estimates)). `generation_stats` has
  one row per `kind` (`sector` fill, or a `scatter` layer of bright
  stars) and density `bucket` (`floor(2 * log10(density / 0.01))`: two
  per decade from 0.01, with no top edge, `density_low`/`density_high`
  its edges), each column a decaying average over every task that ever
  finished in it (each new task weighs 5%): `seconds_per_task` (a
  worker's wall time), `seconds_per_system`, `systems_per_task`,
  `stars_per_system`, plus `samples` and the densest task seen
  (`max_density`). `generation_size` has one row per galaxy database:
  `bytes_per_system` (its tables' data and index bytes over its star
  systems, measured after each run), `systems`, `total_bytes`. A galaxy
  reset keeps both. These and the work queue's tables are the only
  control tables the generator writes; it reaches them with its own
  MySQL account and runs without them when it can't.

`stellarObjects/adminAuth.py` is the only code that reads/writes the
admin tables directly — `bootstrap_control_schema` creates the schema and,
when `admin_users` is empty, seeds an `admin` row with a random first
password that `migrateDb.py` prints once (`migrateDb.py` calls this
automatically, alongside its usual content-schema migration), and every
other function there implements one piece of the login/session/API-key/
audit lifecycle `html/api/auth.py`'s routes expose.

### Booleans and tri-state flags

MySQL has no dedicated boolean type either. Plain booleans are `TINYINT(1)`
`0`/`1` with a `CHECK` constraint. `SystemConfig`'s tri-state flags (`True`/
`False`/`None` in Python — force-on / force-off / random) are nullable
`INTEGER` `0`/`1`/`NULL`.

### Rendered wiki text and URLs

No rendered page text is stored (v29 dropped `star_systems.wikitext_content`/
`markdown_content`). Every value the page shows has its own column, so
`_db.load_star_system` rebuilds the generation object graph and
`stellarObjects/systemRender.py` renders either format on demand — the
same `StarSystem.__str__` generation used to run once before saving, so
the text is identical, except that it now follows later changes (renames,
names made unique, orbit ticks' current wobble and comet positions).
`src/checkRenderParity.py` compares the two on a database still at v28.
`mediawiki_url`/`wikijs_url` (v23 — see `schema.sql`'s header comment)
record where that page lives on each wiki, once `POST
/api/systems/<id>/wiki` (`src/wikiClient/`, `web/system_pages.py`'s "Upload to
Wiki" form) has actually uploaded it there — one system is one wiki page;
individual stars/planets/moons are sections within that one page, not
separate pages. Existence on a given wiki is never a separate stored
flag — it's exactly "the matching URL column is not NULL"; `web/system_pages.py`
links to that page (opening in a new tab) whenever either is set. `sectors.wiki_url` is the
per-sector equivalent (see that table's own column doc above) — a single
column, since a sector's page is generated fresh at upload time rather
than persisted the way a system's is.

## Planned changes (not built)

Schema changes the open TODO items plan, as of 2026-10-02 (7.132.433,
galaxy schema v50, control schema v7). None of these exists yet; each
lands as a numbered migration when its item is built, and this section
moves into the table descriptions below. Writers of the galaxy schema go
one at a time in this order: GEN.44, PERF.11 with MAP.86, NAV.10,
API.11; of the control schema, API.9 first, then USR.2, USR.4, USR.7 and
NAV.19 ([plan/notes.md](plan/notes.md)). Why and how the reproducibility
pieces fit together is in
[design/reproducible-galaxies.md](design/reproducible-galaxies.md).

Galaxy schema:

| Item | Phase | Change |
|---|---|---|
| DB.6 | 0 | `galaxy_shape` gains, next to the galaxy seed (v51), the version that made the galaxy: the 22-hex-digit version key, the full MAJOR.REVISION.BUILD string, the Python version and the platform. A new `generation_runs` table holds one row per run that changes the galaxy: command and options, version, start and end times, outcome. No migration of old galaxies: GEN.39 starts from a wiped galaxy. |
| GEN.44, PERF.11 | 1 | One per-sector stats table keyed by sector address (an unfilled sector has no `sectors` row): the backfill level (-1 never backfilled, the dimmest L_sun reached, 0 fully generated), expected and actual density, and the galaxy-wide expected-against-actual decaying average. MAP.86's per-sector colour, saturation and lightness go in the same table. |
| DB.7 | 1 | Each `sectors` row records the version key and full version string that generated it. |
| DB.9 | 1 | Each sector gets a content checksum (the hash of GEN.58's fingerprint); the Reed-Solomon parity itself lives in a file outside the database. |
| NAV.10 | 1 | Indexes on positions, for the route corridor query. |
| API.10, API.11 | 2 | Run reservations (claimed sectors and id ranges per run) and staging tables for uploaded batches. |

Control schema:

| Item | Phase | Change |
|---|---|---|
| API.15 | 0 | Every API call logged: time, route, account (the key's owner, the signed-in admin, or "god" for the console), how it came in (API key, web session or console) and the HTTP response code. Where the rows are kept is settled when it is built. |
| OPS.13 | 1 | A version-key history table: one row per galaxy per update with the galaxy seed, the version key, SHA-256 hashes of nltk's `words` corpus files, `offensive_words.txt`, any name lists and `requirements.lock`, and the date; only the last 10 rows per galaxy are kept. OPS.15 (phase 2) adds the fingerprint of a small fixed region to each row. |
| GEN.59 | 1 | A pending-delta table: admin edits, deletes and regenerate seeds, by stable address path, as they happen, with the positional-update epoch. The daily merge (GEN.61, phase 2) folds them into a new settings file and clears them only after the file is written and read back. |
| API.9 | 1 | `admin_api_keys` gains a scope column (read, admin, upload). |
| USR.2, USR.4, USR.7, NAV.19 | 3+ | Accounts with roles, invite links, per-account bookmarks and `user_courses`. |

## Tables

### `schema_migrations`

One row per applied DDL migration step — see "Versioning" above. Not tied
to any generated object; exists purely to track `MAX(version)` against
`SCHEMA_VERSION`.

| Column | Type | Null | Notes |
|---|---|---|---|
| `version` | INT | PK | |
| `applied_at` | TIMESTAMP | NOT NULL, default `CURRENT_TIMESTAMP` | |

### `sectors`

One row per generated sector.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `name` | TEXT | NOT NULL | e.g. `"Voranthis Sector"` |
| `edge_mpc` | DOUBLE | NOT NULL | Cube edge length, milliparsecs. Native generator value is `SpaceSector.edge_ly` (light-years). |
| `center_x_pc`, `center_y_pc`, `center_z_pc` | DOUBLE | nullable | The sector's center, in a galaxy-frame Cartesian coordinate system whose origin is the galactic center (parsecs — see `docs/design/galaxy-coordinate-system.md`). NULL together iff this sector has never been placed in a galaxy (`generate.py sector`'s own standalone CLI, or a sector migrated from a pre-v4 database). |
| `galactic_radius_pc` | DOUBLE | nullable | `sqrt(x^2+y^2+z^2)`, persisted (not just derivable) so "sectors within radius R of the core" is a plain indexed range scan — same treatment `star_systems.quadrant` gets. NULL iff the center columns are NULL. |
| `ring_index`, `layer_index`, `ring_slot_index` | INT | nullable | This sector's address on the cylindrical grid (v32, `stellarObjects/galaxyGeometry.py`; see `docs/design/galaxy-coordinate-system.md`, "Cylindrical sector grid"): the radial ring (one sector edge, 4 pc, wide), the height layer (layer 0 centered on the plane) and the angular slot within the ring. Independently nullable from the center/radius columns above (not part of the same CHECK) — a sector could in principle have a hand-authored galaxy position without this particular placement algorithm's own addressing. |
| `wiki_url` | VARCHAR(2048) | nullable | Where this sector's summary page lives on a wiki (v23) — set either by `POST /api/sectors/<id>/wiki` (uploading `html/api/routes.py`'s `_sector_wiki_content`) or directly via `PATCH /api/sectors/<id>` (the `/admin` page's manual-link form, `web/admin_pages.py`). A single column, not one per backend the way `star_systems.mediawiki_url`/`wikijs_url` are — a sector has no persisted rendered page of its own to independently re-upload to a second backend, so only one link is ever tracked at a time. NULL means no page yet. |
| `created_at` | TIMESTAMP | NOT NULL | Added in v27. See "Row timestamps (v27)" above. |
| `modified_at` | TIMESTAMP(3) | NOT NULL, `ON UPDATE CURRENT_TIMESTAMP(3)` | Added in v27. Indexed. |

A `CHECK` constraint enforces `center_x_pc`/`center_y_pc`/`center_z_pc`/
`galactic_radius_pc` being NULL together (see "v3 → v4" in "Schema
history" above for why that addition needed a wholesale table rewrite
rather than an incremental `ALTER TABLE`). A cell's corners are
closed-form (`galaxyGeometry.sector_cell_vertices_pc`), so none are
stored (v32 dropped the old `sector_vertices` table). A
`UNIQUE (ring_index, layer_index, ring_slot_index)` constraint
(`uq_sectors_address`) guarantees at most one sector per galaxy address —
NULL-together rows (never placed in a galaxy) don't collide with each
other or with a placed sector, ordinary SQL `NULL` semantics for `UNIQUE`.
This is what turns a lazy-generation race (`generate.ensure_sector_generated`)
into a clear `IntegrityError` its caller recovers from, instead of a silent
duplicate row at the same address.

### `galaxy_shape`

The galaxy-wide density "skeleton" (v8) — a singleton row (`id` pinned to
`1`) holding everything needed to recompute any sector's exact position
and density on demand, built by `generate.py plan`. Deliberately **not** one
row per sector: a sector's position (`stellarObjects.galaxyGeometry.
sector_position_pc`), density (`stellarObjects.galaxyDensity.
relative_density`), and corners (`galaxyGeometry.sector_cell_vertices_pc`)
are all pure deterministic functions of its `(ring_index, layer_index,
ring_slot_index)`
address plus this handful of galaxy-wide numbers — cheaper to recompute
on demand (sub-millisecond per sector) than to look up, so none of it is
stored per address. What genuinely needs precomputing — where the galaxy
has any content at all — lives in `galaxy_layer` below instead.
Building or rebuilding the skeleton replaces this row (and every
`galaxy_layer` row) wholesale; there is no partial update, since a
full build takes a few milliseconds even at real Milky-Way scale (a
closed-form walk over ~3,900 rings and ~320 layers, not a per-sector scan over
billions of candidates).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK, `CHECK (id = 1)` | Pinned to `1` — there is exactly one galaxy. |
| `disk_scale_length_pc`, `disk_scale_height_pc`, `bulge_scale_radius_pc`, `bulge_amplitude`, `pitch_angle_rad`, `arm_amplitude`, `spiral_reference_radius_pc`, `spiral_reference_angle_rad`, `k_norm` | DOUBLE | NOT NULL | `stellarObjects.galaxyDensity.GalaxyShape`'s own fields, verbatim — see that module for what each means and how `k_norm` is calibrated. |
| `arm_count` | INT | NOT NULL | Same source. |
| `edge_pc` | DOUBLE | NOT NULL | The sector edge length this skeleton was built at, parsecs: always the standard `program_constants.DEFAULT_SECTOR_EDGE_PC` (4) since v33. |
| `expected_system_count_at_density_1` | DOUBLE | NOT NULL | `SpaceSector(edge_ly=...).expected_system_count()` at `relative_density = 1` — cached since every qualification check needs it. |
| `outer_ring_index` | INT | NOT NULL | The last ring with any qualifying content — this galaxy's real edge, found by `generate.py plan` (the plane's layer reaches farthest), not an arbitrary radius. Named `outer_shell_index` before v32. |
| `bright_star_min_luminosity_sol` | DOUBLE | nullable | Added in v43. The luminosity threshold the bright-star scatter used (see `bright_stars` below). A sector fill reads this, not the current constant. NULL until a scatter has run. |
| `bright_star_seed` | BIGINT UNSIGNED | nullable | Added in v43. The scatter's seed: since v51 the top 63 bits of the galaxy seed's unit seed `bright-stars:scatter` (GEN.39). |
| `galaxy_seed` | BINARY(16) | nullable | Added in v51 (GEN.39). The galaxy's 128-bit seed, shown and typed as 32 hex digits (`generate.py plan --seed`), written by the first plan and kept by every later one. Every sector, bright-star scatter, band and backfill block draws from SHA-256 of it and the unit's `kind:address` (`galaxySeed.unit_seed`). NULL only on a galaxy planned before v51. |

### `galaxy_layer`

The galaxy's outline, one row per layer that can hold content (v33,
replacing v32's `galaxy_ring_band`), from the highest layer to the lowest.
Each layer is a circular slice holding rings 0 through `outer_ring_index`,
the last ring whose sector centers could clear the qualification
threshold (`stellarObjects.galaxySkeleton.build_layer_extents`). It comes
from an exact upper bound over every spiral-arm angle, so a sector inside
a layer's extent isn't guaranteed to qualify, but one outside it is
guaranteed not to. The exact per-sector answer (one `relative_density`
evaluation) is made when the address is visited
(`generate.ensure_sector_generated`). The bound falls steadily with radius
and height, so the layers are symmetric about the plane and shrink away
from it. A layer with no qualifying content has no row.

| Column | Type | Null | Notes |
|---|---|---|---|
| `layer_index` | INT | PK | The layer (signed; 0 is centered on the plane). |
| `outer_ring_index` | INT | NOT NULL | The last ring this layer reaches; it holds rings 0 through this one. |

### `galaxy_column`

The same outline seen from the side (v33): one row per ring, holding the
highest and lowest layer that ring's column of sectors reaches (its stack
bound). It is derived from `galaxy_layer`
(`galaxySkeleton.column_extents`) and rewritten with it whenever
`generate.py plan` runs.

| Column | Type | Null | Notes |
|---|---|---|---|
| `ring_index` | INT | PK | The ring. |
| `layer_index_min` | INT | NOT NULL | The lowest layer this ring reaches. |
| `layer_index_max` | INT | NOT NULL | The highest layer this ring reaches. |

Together these two tables are the galaxy's bounds.
`_db.get_galaxy_bounds` loads them as a `galaxySkeleton.GalaxyBounds`,
and every generation path checks an address against it before generating
anything there, so a sector is never placed outside the galaxy.

### `orbit_simulation_state`

Singleton row (same pattern as `galaxy_shape` above) added in v9, tracking
when `updateOrbits.py` last advanced every planet's/moon's
`orbital_phase_deg` in this database — absent entirely until that
script's first run against a given database (it creates this row
itself). `stellarObjects._db.get_orbit_update_elapsed_years` reads it
(via `TIMESTAMPDIFF` against `NOW()`, server-side, rather than trusting
the calling process' own clock) to compute how much simulated time has
passed since the last update.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK, `CHECK (id = 1)` | Always `1` — singleton. |
| `last_updated_at` | TIMESTAMP | NOT NULL | When `updateOrbits.py` last ran against this database. |

Meant to run on a schedule, not on every deploy (`install.sh`/`update.sh`
don't call it) -- e.g. a monthly cron entry:

```
0 3 1 * * cd /var/lib/planetGen && python3 src/updateOrbits.py >> /var/log/planetgen-orbits.log 2>&1
```

or, on a systemd-based (Ubuntu/Debian) host, the equivalent systemd timer
under [`../examples/maintenance/`](../examples/maintenance/) -- journald
captures the run's output automatically, with no logfile/logrotate entry
to maintain, and (unless installed with `--skip-update-timer`) `sudo
./update.sh` itself is scheduled too, 30 minutes ahead of the orbit
update on the same monthly run:

```
sudo examples/maintenance/install-maintenance-timer.sh [database ...]
```

`updateOrbits.py` mutates rows, so it needs the same read-write database
account `generate.py` and the web app use, not a `SELECT`-only one you may
have made for `queryDb.py`.

### `system_configs`

One row per `SystemConfig` "recipe" — the generation parameters a
`StarSystem` was built from (`src/stellarObjects/config.py`,
`SERIALIZABLE_FIELDS`).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `markdown` | INTEGER (0/1) | NOT NULL, default 0 | Whether this recipe renders Markdown (1) or wikitext (0). Not tri-state — always a concrete choice. |
| `habitable_world` | INTEGER (0/1) | tri-state | Force/forbid/random a habitable world. |
| `asteroid_belt` | INTEGER (0/1) | tri-state | Force/forbid/random an asteroid belt. |
| `comets` | INTEGER (0/1) | tri-state | Force/forbid/random a star-bound comet (v50). |
| `large_star` | INTEGER (0/1) | tri-state | Force/forbid/random a larger star. |
| `moons` | INTEGER (0/1) | tri-state | Force/forbid/random moon generation. |
| `max_planets` | INTEGER (0/1) | tri-state | Force max vs. min planet count. |
| `planets` | INTEGER (0/1) | tri-state | Force/forbid at least one planet or belt. |
| `star_type` | TEXT | nullable | Explicit spectral type, e.g. `"G2V"`. |
| `name` | TEXT | nullable | Forced system name, if any. |
| `age` | TEXT | nullable, CHECK IN ('young','old') | |
| `intelligent_life` | INTEGER (0/1) | tri-state | |
| `binary_system` | INTEGER (0/1) | tri-state | |
| `wide_binary` | INTEGER (0/1) | tri-state | Force a wide (S-type) or close binary (v50). |
| `num_orbits` | INTEGER | nullable | Explicit orbital slot count override. |

### `system_config_slots`

Child table for a config's variable-length `SLOTS` list (per-orbit
overrides).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `config_id` | INTEGER | FK -> `system_configs.id`, `ON DELETE CASCADE` | |
| `orbit_index` | INTEGER | NOT NULL | Position in the `SLOTS` list. |
| `type` | TEXT | nullable, CHECK IN ('planet','asteroid_belt') | |
| `planet_class` | TEXT | nullable | Specific class letter, planet slots only. |
| `moons` | INTEGER | nullable | Exact moon count override for this slot. |

### `star_systems`

One row per generated system (single-star or binary).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE SET NULL`, nullable | NULL for a standalone system never placed in a sector. |
| `system_config_id` | INTEGER | FK -> `system_configs.id`, NOT NULL | The recipe this system was generated from. For binaries, this is the *shared* config the primary/proxy/planets all use — the secondary star is generated with `LARGE_STAR` temporarily forced off (`systemData.py`), a transient change this schema doesn't capture, so there's no second config row. |
| `name` | TEXT | NOT NULL | The system's display name, e.g. `"Voranthis"`. Unique across systems and sectors (v24). Its stars, planets and moons are named from it (v34, see below); renaming it (`PATCH /api/systems/<id>`, or a uniqueness decoration) renames every one still carrying it. |
| `position_x_mpc`, `position_y_mpc`, `position_z_mpc` | DOUBLE | nullable | Position relative to the sector's cubic center. NULL iff not placed in a sector. |
| `quadrant` | TEXT | nullable, CHECK IN ('I'..'VIII') | The sector octant label derived from the position above (see "Quadrant labeling" below). NULL iff position is NULL. |
| `location` | TEXT | nullable | Human-readable "sector name + nearest neighbors" summary, e.g. `"Voranthis Kelmoor — nearest: Alpha Vesta (4.2 ly), Beta Cerise (7.8 ly), Gamma Ost (9.1 ly)"` — up to 3 neighbors, nearest first, computed once at write time from `SpaceSector.nearest_neighbors` (see "v3" in "Schema history" above). NULL iff position is NULL. |
| `is_binary` | INTEGER (0/1) | NOT NULL, default 0 | True for either binary configuration as of v15 (see `binary_configuration`) — before v15, true only ever meant a `'close'` pair. |
| `binary_configuration` | TEXT | nullable, CHECK IN ('close','wide') | Added in v15. `'close'` (P-type/circumbinary, `doubleStar.BinaryStarProxy`) or `'wide'` (S-type, `wideBinary.WideBinaryPair`) — NULL for a single star. The authoritative discriminator going forward; see "Binary configuration" below for which other columns apply to which value. |
| `binary_separation_km` | DOUBLE | nullable | Orbital separation between the two stars. NULL for single-star systems. Reused unchanged for both binary configurations (see "Binary configuration" below). |
| `binary_eccentricity` | DOUBLE | nullable | Added in v15. The pair's own orbital eccentricity — `0` for a `'close'` pair (tidal circularization is a legitimate simplification at its 0.05-0.25 AU separations), a real, generally non-zero value (sampled from a "thermal" distribution) for a `'wide'` pair, which never circularizes. |
| `binary_periapsis_km`, `binary_apoapsis_km` | DOUBLE | nullable | Added in v15. Closest/farthest separation between the two stars — `binary_separation_km * (1 ∓ binary_eccentricity)`. Both equal `binary_separation_km` for a `'close'` pair (e=0). |
| `binary_type` | TEXT | nullable | e.g. `"Binary (G/K)"`. **Not** the same thing as `binary_configuration` above — this pre-existing column holds a `'close'` pair's merged spectral-summary string; NULL for a `'wide'` pair (no merged effective star exists to summarize). |
| `binary_temperature_k` | DOUBLE | nullable | Average of the two stars' temperatures. |
| `binary_radius_km` | DOUBLE | nullable | The larger constituent star's radius (used as an approximation). |
| `binary_effective_mass_kg` | DOUBLE | nullable | Sum of both stars' masses. |
| `binary_effective_luminosity_w` | DOUBLE | nullable | Sum of both stars' luminosities. |
| `binary_age_gy`, `binary_lifespan_gy` | DOUBLE | nullable | Max of the two stars' age/lifespan. `binary_lifespan_gy` NULL = infinite (a white-dwarf constituent). |
| `binary_habitable_zone_inner_km`, `_outer_km` | DOUBLE | nullable | Computed from the pair's combined luminosity. |
| `binary_system_perimeter_km` | DOUBLE | nullable | Hill sphere, combined mass. |
| `binary_heliosphere_radius_km` | DOUBLE | nullable | |
| `binary_galactic_orbital_speed_kms`, `binary_galactic_orbital_period_gy` | DOUBLE | nullable | Added in v10. Circular orbital speed/period around the galactic center (see `stars.galactic_orbital_speed_kms` below) — independent of mass, so identical to the primary/secondary stars' own values, just mirrored here for the combined-pair row. |
| `binary_galactic_orbital_phase_deg`, `binary_galactic_min_update_interval_years` | DOUBLE | nullable | Added in v13. Same pair as `stars.galactic_orbital_phase_deg`/`galactic_min_update_interval_years` below, mirrored here — always identical to both constituent stars' own values (see "Schema history" above for why). |
| `binary_mutual_orbital_period_years`, `_speed_kms` | DOUBLE | nullable | Added in v13. The pair's own mutual orbit around each other — entirely separate from, and vastly faster than, the galactic orbit above. Kepler's third law / circular-orbit speed (`planetPhysics.calculate_orbital_period_years`/`utils.circular_orbital_speed_kms`) applied to `binary_separation_km`/`binary_effective_mass_kg`. |
| `binary_mutual_orbital_inclination_deg`, `_ascending_node_deg`, `_phase_deg` | DOUBLE | nullable | Added in v13. Orients the mutual orbit in 3D and tracks the pair's current position within it — same `utils.orbital_position_au` convention as `planets.orbital_inclination_deg`/etc, but drawn from the full `[0, 180)`/`[0, 360)` range (no small-tilt bias — a binary's mutual orbital plane has no preferred alignment the way a planet's protoplanetary-disk-derived orbit does). `_phase_deg` is advanced by `_db.advance_orbital_phases`, guarded by the interval below. |
| `binary_mutual_min_update_interval_years` | DOUBLE | nullable | Added in v13. Floating-point update guard for `binary_mutual_orbital_phase_deg`, same formula as `planets.min_update_interval_years`. |
| `binary_mutual_position_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v14. The secondary's Cartesian position relative to the primary — same "position relative to whatever this orbit is around" convention as `planets.position_x/y/z_km` (`utils.orbital_position_au`, applied to `binary_separation_km` and the mutual-orbit orientation columns above). Recomputed by `_db.advance_orbital_phases` in lockstep every time `binary_mutual_orbital_phase_deg` advances. **Unchanged in meaning by v20** — still the true separation, not a barycenter-reduced value. |
| `binary_primary_position_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v20. The primary star's own offset from the pair's barycenter — `-binary_secondary_mass_fraction * binary_mutual_position_*_km`. Reused unchanged for both binary configurations, like `binary_mutual_position_*_km` above. |
| `binary_secondary_position_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v20. The secondary star's own offset from the barycenter — always equal to `binary_primary_position_* + binary_mutual_position_*`, but stored explicitly rather than derived on read (same convention `binary_mutual_position_*_km` itself already set). |
| `binary_secondary_mass_fraction` | DOUBLE | nullable | Added in v20. `secondary_mass / (primary_mass + secondary_mass)`, constant since masses don't change — stored so `_db.advance_orbital_phases` never needs to join back to `stars` for either mass. |
| `binary_planetary_wobble_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v20. NULL unless `binary_configuration = 'close'`. Additional pair-wide wobble from circumbinary planets (`planets.star_id IS NULL`) — modeled as one shared wobble applied to the whole pair rather than split unevenly between primary/secondary, which would need a real 3+-body solve. |
| `system_flavor_text` | TEXT | nullable | Decided once at generation time (Phase 0 fix). |
| `runaway_class` | VARCHAR(16) | nullable, `runaway` or `hypervelocity` | Added in v37. NULL for an ordinary star; set by `generate.flag_fast_stars`. |
| `runaway_speed_kms` | DOUBLE | nullable | Added in v37. The star's speed relative to its neighbors when `runaway_class` is set. |
| `schema_version` | INTEGER | NOT NULL, default 1 | See "Versioning" above. |
| `mediawiki_url` | TEXT | nullable | Where this system's page lives (or should live) on MediaWiki. |
| `wikijs_url` | TEXT | nullable | Where this system's page lives (or should live) on Wiki.js. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | FK -> `nebulae.id` / `supernova_remnants.id`, `ON DELETE SET NULL`, nullable | Added in v39. The innermost nebula or supernova remnant whose sphere holds this system; at most one is set. See "Containment (v39)" above. |
| `created_at` | TIMESTAMP | NOT NULL, default `CURRENT_TIMESTAMP` | |
| `modified_at` | TIMESTAMP(3) | NOT NULL, `ON UPDATE CURRENT_TIMESTAMP(3)` | Added in v27. Also bumped when a child row (star, planet, moon, belt, comet) changes. |

**Quadrant labeling.** `quadrant` reuses the generator's own octant scheme
(`src/stellarObjects/spaceSector.py`'s `classify_octant`, backed by
`program_constants.SECTOR_OCTANT_LABELS`): each axis's sign (`x >= 0`,
`y >= 0`, `z >= 0`) picks one of 8 Roman-numeral labels, `I` through
`VIII` — the same labels the generator's own `format_named_location`
produces (e.g. `"Quadrant III (2.10, 4.40, 1.05 ly from center)"`). The
code calls these "quadrants" even though three signed axes make it a 3D
*octant* scheme, not a 2D quadrant one — this schema keeps that same
terminology and label set for consistency with the generator's own output,
rather than introducing a different name for the same thing. It's stored
rather than only computed on read because it's what a query like "every
system in Quadrant III" filters on directly, without a UDF or generated
column.

**Binary configuration.** `binary_configuration` (added v15) decides which
other `binary_*` columns are populated: `binary_separation_km`,
`binary_eccentricity`, `binary_periapsis_km`, `binary_apoapsis_km`, and
every `binary_mutual_*`/`binary_mutual_position_*_km` column are reused
UNCHANGED for both `'close'` and `'wide'` pairs (a `'wide'` pair's own
mutual orbit fits the exact same shape a `'close'` pair's already
occupies). Every other `binary_*` column (`binary_type`,
`binary_temperature_k`, `binary_radius_km`, `binary_effective_mass_kg`,
`binary_effective_luminosity_w`, `binary_age_gy`, `binary_lifespan_gy`,
`binary_habitable_zone_inner_km`/`_outer_km`, `binary_system_perimeter_km`,
`binary_heliosphere_radius_km`, and the four `binary_galactic_orbital_*`
columns) describes the `'close'` pair's *merged* `BinaryStarProxy` and is
NULL for a `'wide'` pair — a `'wide'` pair's two stars never merge, so
there's no combined effective star for these to describe; each star's own
equivalent data lives on its own `stars` row instead (including its own
`wide_binary_a_crit_km`, see below).

### `stars`

One row per individual star: one row for a single-star system, two rows
(`primary`/`secondary`) for either binary configuration. There is
**never** a row for a `'close'` pair's `BinaryStarProxy` itself — its
combined-pair values live on `star_systems` above (the `binary_*`
columns), not here. A `'wide'` pair's two stars, by contrast, are two
ordinary rows just like a single star's, since neither one merges into
anything.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `star_system_id` | INTEGER | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | |
| `role` | TEXT | NOT NULL, CHECK IN ('primary','secondary','single') | |
| `name` | TEXT | NOT NULL | A single star shares the system's name (renaming one renames the other). A binary's stars put their own word after it, e.g. `"Voranthis Kelmoor"` and `"Voranthis Ostra"`, with no A/B letters (v34). |
| `star_type` | TEXT | NOT NULL | Full descriptive string, e.g. `"G2V Yellow Main Sequence Star"` — unrelated to `planets.body_type`'s single-character code. |
| `yerkes_class` | TEXT | NOT NULL | e.g. `"V"`, `"VII"` (white dwarf). |
| `mass_kg`, `radius_km`, `luminosity_w` | DOUBLE | NOT NULL | |
| `temperature_k` | DOUBLE | NOT NULL | |
| `age_gy` | DOUBLE | NOT NULL | |
| `lifespan_gy` | DOUBLE | nullable | **NULL = `float('inf')`** — white dwarfs (`yerkes_class` `'VII'`/`'D'`). Never the JSON `Infinity` token. |
| `habitable_zone_inner_km`, `_outer_km` | DOUBLE | NOT NULL | |
| `system_perimeter_km` | DOUBLE | NOT NULL | Hill sphere relative to the galaxy. |
| `heliosphere_radius_km` | DOUBLE | NOT NULL | |
| `galactic_orbital_speed_kms` | DOUBLE | NOT NULL | Added in v10. Circular orbital speed around the galactic center (`utils.calculate_galactic_orbit`), from this system's actual distance from the galactic center where known (a sector-placed system), or the fixed `physical_constants.GALACTIC_CENTER_DISTANCE_LY` fallback otherwise — same fallback convention as `system_perimeter_km`. |
| `galactic_orbital_period_gy` | DOUBLE | NOT NULL | Added in v10. Orbital period for the circular orbit above, in billions of years (Gy) — the same unit `age_gy`/`lifespan_gy` use. |
| `galactic_orbital_phase_deg` | DOUBLE | NOT NULL | Added in v13. This star's current angular position around its galactic orbit — the same role `planets.orbital_phase_deg` plays, advanced by `_db.advance_orbital_phases`. Both stars of a binary pair always carry the identical value (see "Schema history" above for why — `StarSystem.__init__` rolls it once and threads it to primary/secondary/proxy alike). |
| `galactic_min_update_interval_years` | DOUBLE | NOT NULL | Added in v13. Floating-point update guard for `galactic_orbital_phase_deg`, same formula as `planets.min_update_interval_years` (`utils.minimum_update_interval_years`), applied to `galactic_orbital_period_gy * 1e9` years. |
| `wide_binary_a_crit_km` | DOUBLE | nullable | Added in v15. This star's own Holman & Wiegert (1999) critical semi-major axis (`utils.holman_wiegert_critical_semimajor_axis`) — the maximum orbit distance that stays long-term stable given its companion's perturbation. NULL for a single star or either constituent of a `'close'` pair; populated for both stars of a `'wide'` pair. |
| `reflex_offset_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v20. This star's own displacement from its nominal fixed point, from the combined pull of every planet orbiting it directly (`planets.star_id`) — see `utils.calculate_reflex_offset`. NULL/0 with no planets. Class-blind: applies identically to an anchored `black_holes`/`neutron_stars` row. |

### `planets`

One row per **top-level** planet (as of schema v2 — moons live in their
own `moons` table below, not here; see "Schema history" above). Covers
both terrestrial and gas-giant bodies (`body_type`).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `star_system_id` | INTEGER | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | Always the reliable owning link, regardless of `star_id`. |
| `star_id` | INTEGER | FK -> `stars.id`, `ON DELETE SET NULL`, nullable | The specific star this planet orbits, when that's a real stored `stars` row — true for every single-star system, and, as of v15, every `'wide'` binary's planets too (each orbits one specific constituent star). Still **NULL for a `'close'` binary's planets**: the generator builds those against the merged `BinaryStarProxy`, never one individual constituent star, and the proxy has no `stars` row to point at. |
| `orbital_index` | INTEGER | NOT NULL | Position in the star's ordered `planets` list. |
| `body_type` | TEXT | NOT NULL, CHECK IN ('t','g') | Terrestrial or gas giant. Unrelated to `stars.star_type`. |
| `name` | TEXT | NOT NULL | Generated as the star it orbits plus a roman numeral in orbit order, e.g. `"Voranthis II"` (v34); asteroid belts take no number. A close pair's planets use the system name, a wide pair's their own star's. Can be renamed (`PATCH /api/planets/<id>`). |
| `planet_class` | TEXT | nullable | e.g. `"M"`. |
| `distance_km` | DOUBLE | NOT NULL | From the star. |
| `radius_km`, `mass_kg` | DOUBLE | NOT NULL | |
| `volume_km3` | DOUBLE | NOT NULL | Derived, stored not recomputed. |
| `period_years` | DOUBLE | NOT NULL | Derived, stored not recomputed. |
| `zone` | TEXT | nullable, CHECK IN ('h','e','c') | Hot / ecosphere / cold. |
| `description` | TEXT | nullable | |
| `gravity_g` | DOUBLE | nullable | Surface gravity, g's. |
| `surface_temperature_k` | DOUBLE | nullable | |
| `density_g_cm3` | DOUBLE | nullable | |
| `atmosphere` | TEXT | nullable | Composition description. |
| `atm_density`, `atm_molar_density`, `atmospheric_pressure_pa` | DOUBLE | nullable | |
| `composition` | TEXT | nullable | Descriptive string — contrast `asteroid_belt_composition`'s structured (component, concentration) rows. |
| `scale_height_km`, `hill_radius_km`, `min_orbit_distance_km` | DOUBLE | nullable | |
| `habitable_zone_inner_km`, `_outer_km` | DOUBLE | NOT NULL | Copied from the host star at generation time, never recomputed. |
| `life_chemical`, `evolutionary_speed` | TEXT | nullable | Set by life-data generation, if any. |
| `flavor_text` | TEXT | nullable | |
| `flavor_text_count` | INTEGER | NOT NULL, default 0 | |
| `orbital_inclination_deg`, `orbital_ascending_node_deg` | DOUBLE | NOT NULL | Added in v9. Fixed at generation time — together they orient this (circular) orbital plane in 3D. |
| `orbital_phase_deg` | DOUBLE | NOT NULL | Added in v9. This body's current position angle around its orbit — the one orbital-motion column that changes over time, advanced in place by `updateOrbits.py` (see `orbit_simulation_state` below). |
| `position_x_km`, `_y_km`, `_z_km` | DOUBLE | NOT NULL | Added in v11. This body's Cartesian position relative to its orbital anchor — the star (or a binary's combined center) for a planet — derived from `distance_km` and the three orbital-motion columns above (`utils.orbital_position_au`). Changes in lockstep with `orbital_phase_deg` as `updateOrbits.py` advances it. |
| `orbital_speed_kms` | DOUBLE | NOT NULL | Added in v11. Constant circular-orbit speed (`utils.circular_orbital_speed_kms`, `v = 2*pi*r/T`). Only changes if `distance_km`/`period_years` do (e.g. `StarSystem.validate_system` resolving an orbital overlap at generation time), never from phase advancing alone. |
| `min_update_interval_years` | DOUBLE | NOT NULL | Added in v12. Not a narrative stat -- a floating-point update guard for `_db.advance_orbital_phases`: the shortest `elapsed_years` worth calling it for, below which the phase delta added is smaller than `orbital_phase_deg`'s own double-precision resolution and so is guaranteed to be a no-op write (`utils.minimum_update_interval_years`, `period_years * math.ulp(360.0) / 360`). Like `orbital_speed_kms`, only changes if `distance_km`/`period_years` do. |
| `rotation_period_hours` | DOUBLE | NOT NULL | Added in v9. Axial rotation ("day length") — a static descriptive stat; no rotational phase is tracked. |
| `reflex_offset_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v20. This planet's own displacement from its nominal fixed point, from the combined pull of its own moons (`moons.planet_id`) — see `utils.calculate_reflex_offset`. NULL/0 with no moons. **Not present on `moons`** — a moon never hosts its own moons. |

### `planet_evolutionary_paragraphs`

Child table for a planet's variable-length evolutionary narrative
(`Planet.evolutionary_data`).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `planet_id` | INTEGER | FK -> `planets.id`, `ON DELETE CASCADE`, NOT NULL | |
| `position` | INTEGER | NOT NULL | List order. |
| `paragraph` | TEXT | NOT NULL | |

### `planet_reflection_spectrum`

Child table for a planet's reflection-spectrum descriptor lists
(`Planet.reflection_spectrum_visible`/`_non_visible`). Deliberately a
normalized table, not a JSON column — this schema has no JSON-blob columns
anywhere.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `planet_id` | INTEGER | FK -> `planets.id`, `ON DELETE CASCADE`, NOT NULL | |
| `spectrum_type` | TEXT | NOT NULL, CHECK IN ('visible','non_visible') | |
| `position` | INTEGER | NOT NULL | List order. |
| `value` | TEXT | NOT NULL | |

### `moons`

One row per moon — added in schema v2, in place of the v1 design where a
moon was a `planets` row self-referencing via `parent_planet_id`. Exactly
the same column shape as `planets` (a moon is a `Planet` instance too,
with `is_moon=True`), except `planet_id` names the specific planet it
orbits. No self-reference here: moons never generate their own moons
(`Planet.__init__` only calls `generate_moons` `if not self.is_moon`). A moon is generated as its planet's numeral plus a letter in orbit order,
e.g. `"Voranthis IIa"` (v34), and can be renamed (`PATCH /api/moons/<id>`).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `planet_id` | INTEGER | FK -> `planets.id`, `ON DELETE CASCADE`, NOT NULL | The planet this moon orbits. |
| `star_system_id` | INTEGER | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | Redundant with the owning planet's own `star_system_id` — kept here too so a moon can be queried/joined to its system without an extra hop through `planets`. |
| `star_id` | INTEGER | FK -> `stars.id`, `ON DELETE SET NULL`, nullable | Same value as the owning planet's `star_id` (see that column's note above — NULL for a binary system). |
| `orbital_index` | INTEGER | NOT NULL | Position in the parent planet's `moons` list. |
| `body_type`, `name`, `planet_class`, `distance_km` (from the parent planet), `radius_km`, `mass_kg`, `volume_km3`, `period_years`, `zone`, `description`, `gravity_g`, `surface_temperature_k`, `density_g_cm3`, `atmosphere`, `atm_density`, `atm_molar_density`, `atmospheric_pressure_pa`, `composition`, `scale_height_km`, `hill_radius_km`, `min_orbit_distance_km`, `habitable_zone_inner_km`, `_outer_km`, `life_chemical`, `evolutionary_speed`, `flavor_text`, `flavor_text_count`, `orbital_inclination_deg`, `orbital_ascending_node_deg`, `orbital_phase_deg`, `position_x_km`, `_y_km`, `_z_km`, `orbital_speed_kms`, `min_update_interval_years`, `rotation_period_hours` | — | — | Identical meaning/type/nullability to the same-named column on `planets` above, except `position_x/y/z_km` are relative to *this moon's* orbital anchor — its parent planet, not the star. |

### `moon_evolutionary_paragraphs`

Moon counterpart of `planet_evolutionary_paragraphs` — identical shape,
`moon_id` FK -> `moons.id` (`ON DELETE CASCADE`) instead of `planet_id`.

### `moon_reflection_spectrum`

Moon counterpart of `planet_reflection_spectrum` — identical shape,
`moon_id` FK -> `moons.id` (`ON DELETE CASCADE`) instead of `planet_id`.

### `asteroid_belts`

One row per asteroid belt. Belts have no properties-dict data table
(`AsteroidBelt.to_paragraph_list()` is prose only), so their searchable
columns capture the facts the prose always states instead.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `star_system_id` | INTEGER | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | |
| `star_id` | INTEGER | FK -> `stars.id`, `ON DELETE SET NULL`, nullable | Added in v15. Same "which specific star this orbits" semantics as `planets.star_id` above — NULL for a single star's or a `'close'` binary's belt, set to the owning star for a `'wide'` binary's. |
| `orbital_index` | INTEGER | NOT NULL | |
| `distance_km` | DOUBLE | NOT NULL | Average distance from the star. |
| `lower_limit_km`, `upper_limit_km` | DOUBLE | NOT NULL | The belt's inner/outer boundary. |
| `density` | TEXT | NOT NULL, CHECK IN ('dense','sparse','typical') | |
| `composition_summary` | TEXT | NOT NULL | Human-readable summary built the same way as the prose sentence (`asteroidData.py:119-137`), e.g. `"high concentrations of iron, moderate concentrations of nickel, and trace amounts of platinum"` — searchable without a join, alongside the structured breakdown below. |

### `asteroid_belt_composition`

Structured per-component detail behind `composition_summary` above.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `belt_id` | INTEGER | FK -> `asteroid_belts.id`, `ON DELETE CASCADE`, NOT NULL | |
| `position` | INTEGER | NOT NULL | List order. |
| `component` | TEXT | NOT NULL | e.g. `"iron"`. |
| `concentration` | TEXT | NOT NULL, CHECK IN ('high','moderate','small','trace') | |

### `comets`

Added in v19. One row per star-bound comet (`cometData.Comet`), propagated
via real two-body Kepler/Barker orbital mechanics (`keplerMotion.py`) — see
`docs/design/comet-orbital-realism.md`. Contrast `interstellar_comets`
below, an always-standalone, unbound object on a fixed hyperbolic
trajectory. Deliberately has no `orbital_index`: a comet's `distance_km`
is its current, continuously-varying position along its orbit, not a
fixed slot in the `planets`/`asteroid_belts` orbital-spacing sequence.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `star_system_id` | INTEGER | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | |
| `star_id` | INTEGER | FK -> `stars.id`, `ON DELETE SET NULL`, nullable | Same real semantics as `planets.star_id` above — a real `stars.id` for a single star or a `'wide'` binary's comet; NULL only for a `'close'` binary's, which orbits the merged pair rather than one individually-stored star row. |
| `name` | VARCHAR(255) | NOT NULL | Added in v40. A designation, `P/<host>-<n>` (periodic, under 200 years) or `C/<host>-<n>`, that follows a rename of its star. See "Names (v40)" above. |
| `orbit_type` | TEXT | NOT NULL, CHECK IN ('elliptical','parabolic') | |
| `period_class` | TEXT | nullable, CHECK IN ('jupiter_family','halley_type','long_period') | Only set for `orbit_type = 'elliptical'` — flavor/plausibility metadata only. |
| `nucleus_diameter_km` | DOUBLE | NOT NULL | |
| `composition_summary` | TEXT | NOT NULL | Searchable summary; structured breakdown lives in `comet_composition`. |
| `perihelion_distance_km` | DOUBLE | NOT NULL | |
| `eccentricity`, `inclination_deg`, `arg_periapsis_deg`, `ascending_node_deg` | DOUBLE | NOT NULL | Full 3D orbit orientation. |
| `orbital_period_years`, `mean_anomaly_deg` | DOUBLE | nullable | Only set for `orbit_type = 'elliptical'`. |
| `parabolic_mean_anomaly` | DOUBLE | nullable | Only set for `orbit_type = 'parabolic'`. |
| `min_update_interval_years` | DOUBLE | nullable | Only set for `orbit_type = 'elliptical'` — a parabolic comet's anomaly doesn't wrap, so it has no periodic floating-point-resolution floor to guard. |
| `primary_mass_solar` | DOUBLE | NOT NULL | |
| `is_active` | BOOLEAN | NOT NULL | Coma/tail activity. |
| `distance_km`, `position_x_km`, `position_y_km`, `position_z_km`, `orbital_speed_kms` | DOUBLE | NOT NULL | Current derived orbital state — recomputed by `_db.advance_comet_orbits` as `mean_anomaly_deg`/`parabolic_mean_anomaly` advance over time. |

### `comet_composition`

Structured per-component detail behind `composition_summary` above — a
plain component list (no concentration gradient), mirroring
`interstellar_comet_composition`'s shape rather than
`asteroid_belt_composition`'s.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `comet_id` | INTEGER | FK -> `comets.id`, `ON DELETE CASCADE`, NOT NULL | |
| `position` | INTEGER | NOT NULL | List order. |
| `component` | TEXT | NOT NULL | e.g. `"water ice"`. |

### `black_holes` / `neutron_stars`

Added in v16, for `generate.py phenomenon`'s separate, rarer exotic-phenomenon
generation mode (`stellarObjects/compactRemnant.py`). Satellite tables
extending a `stars` row (the same "extra detail alongside an existing row"
shape `asteroid_belt_composition` has to `asteroid_belts`) when a compact
remnant anchors a full `StarSystem` (`generate.py phenomenon --anchor-system`)
— the owning `stars.yerkes_class` is then the literal marker `'BH'`/`'NS'`
rather than a real Yerkes class. `star_id` is NULL for a remnant generated
standalone (no owning `StarSystem` at all).

**`black_holes`**

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `star_id` | INTEGER | FK -> `stars.id`, `ON DELETE CASCADE`, nullable | NULL for a standalone black hole. |
| `name` | TEXT | NOT NULL | |
| `mass_class` | VARCHAR(16) | NOT NULL, default `'stellar'`, CHECK | Added in v36. `'stellar'` (5-20 Msun), `'intermediate'` (1e2-1e5 Msun, `BLACK_HOLE_INTERMEDIATE_MASS_CHANCE` of rolls) or `'supermassive'` (a galaxy's quiescent central black hole, 1e6-1e8 Msun, at the galactic center with NULL galactic-orbit columns; see `generate.add_galactic_nucleus`). |
| `mass_solar` | DOUBLE | NOT NULL | |
| `event_horizon_radius_km` | DOUBLE | NOT NULL | Schwarzschild radius, `2GM/c^2` — also `BlackHole.radius`. |
| `spin` | DOUBLE | NOT NULL | Dimensionless spin parameter a* in [0, 1). |
| `has_accretion_disk` | BOOLEAN | NOT NULL | |
| `temperature_k`, `luminosity_w` | DOUBLE | NOT NULL | Both 0 unless `has_accretion_disk`. |
| `age_gy` | DOUBLE | NOT NULL | Time since the core-collapse supernova. |
| `galactic_orbital_speed_kms`, `_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | nullable | Added in v17. NULL when `star_id` is set (an anchored remnant's motion lives on its own `stars` row instead); populated only for a standalone black hole. See `stars`' identical columns below. |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE SET NULL`, nullable | Added in v21. The sector this standalone black hole was generated as part of (`generate.generate_sector_phenomena`) or placed near (`generate.py phenomenon --sector-id`), the same convention `nebulae.sector_id` uses. Always NULL when `star_id` is set. |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v21. This black hole's own galaxy-frame center, the same shape `nebulae`'s identical columns use (see below) -- but here it's usually the *exact* position `SpaceSector.add_phenomenon`'s Hill-sphere-aware in-sector placement computed at generation time (never within a neighboring star system's or another compact remnant's own Hill sphere), converted to galaxy-frame coordinates, rather than `compute_phenomenon_placement`'s independent random jitter (still used for `generate.py phenomenon`'s own standalone `--sector-id`, which has no specific in-sector position to convert). |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | nullable, `ON DELETE SET NULL` | Added in v39. The innermost cloud holding it, as on `star_systems`. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |

**`neutron_stars`**

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `star_id` | INTEGER | FK -> `stars.id`, `ON DELETE CASCADE`, nullable | NULL for a standalone neutron star. |
| `name` | TEXT | NOT NULL | |
| `mass_solar`, `radius_km` | DOUBLE | NOT NULL | |
| `spin_period_ms`, `magnetic_field_gauss` | DOUBLE | NOT NULL | |
| `pulsar_type` | TEXT | NOT NULL, CHECK IN ('young','millisecond','non-pulsing') | |
| `surface_temperature_k`, `luminosity_w` | DOUBLE | NOT NULL | `luminosity_w` is thermal blackbody emission (Stefan-Boltzmann), always nonzero. |
| `age_gy` | DOUBLE | NOT NULL | Time since the core-collapse supernova. |
| `galactic_orbital_speed_kms`, `_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | nullable | Added in v17. Same nullable-when-anchored convention as `black_holes` above. |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE SET NULL`, nullable | Added in v21. Same convention as `black_holes.sector_id` above. |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v21. Same convention as `black_holes`' identical columns above. |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | nullable, `ON DELETE SET NULL` | Added in v39. The innermost cloud holding it, as on `star_systems`. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |

### `nebulae`

Added in v16. Always standalone (nothing in this generator places a
nebula within a `StarSystem`); `sector_id` was reserved for a future
sector-context encounter through v17, unused (always NULL) by
`generate.py phenomenon` until v18 gave it one.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE SET NULL`, nullable | v18: the nearest already-generated sector to `center_x/y/z_pc` below -- a convenience "home" link, not this nebula's real geometry (its sphere may overlap several sectors, or none). NULL iff `center_x/y/z_pc` are NULL. `ON DELETE SET NULL` (not `CASCADE`, unlike v16/v17): deleting that sector doesn't delete a nebula that merely happens to be near it. |
| `name` | TEXT | NOT NULL | |
| `nebula_class` | CHAR(1) | NOT NULL | Added in v38: letter class A-Q (`program_constants.NEBULA_CLASSES`). |
| `nebula_type` | TEXT | NOT NULL, CHECK IN ('diffuse','emission','reflection','planetary','dark') | The class's family (`diffuse` added in v38). |
| `dominant_species`, `density_cm3`, `temperature_k`, `extinction_av` | VARCHAR(255) / DOUBLE | NOT NULL | Added in v38: what the cloud holds, its particle density nH (cm⁻³), gas temperature (K) and optical extinction (magnitudes). |
| `radius_ly` | DOUBLE | NOT NULL | |
| `composition`, `formation_cause` | TEXT | NOT NULL | Descriptive strings, one per `nebula_type` (`program_constants.NEBULA_TYPES`). |
| `galactic_orbital_speed_kms`, `_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | NOT NULL | Added in v17. Always populated (a nebula is always standalone). See `stars`' identical columns above. |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v18 (`generate.py phenomenon --sector-id`): this nebula's own galaxy-frame center, in the same Cartesian space `sectors.center_x/y/z_pc` uses -- a sphere (this + `radius_ly`), not a sector-relative offset, since a nebula (up to 200 ly across) is frequently far larger than one sector (default edge 4 pc, ~13 ly) and may overlap several. NULL together: never placed in the galaxy (still the default -- `--sector-id` is optional). See `queryDb.phenomena_near_sector`/`galaxy_placed_phenomena` for how this is read back, and `docs/design/galaxy-coordinate-system.md` for the coordinate system itself. |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | nullable, `ON DELETE SET NULL` | Added in v39. Set only when this nebula nests inside a larger cloud. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |

### `supernova_remnants`

Added in v16. Always standalone. At most one of
`compact_remnant_black_hole_id`/`compact_remnant_neutron_star_id` is ever
set (see `compact_remnant_kind`), and only for a core-collapse progenitor
whose collapsed core is still detectable — never for a Type Ia progenitor,
which leaves nothing behind (thermonuclear disruption of the progenitor
white dwarf).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE CASCADE`, nullable | The sector it was generated as part of (`generate.py sector`) or placed near (`generate.py phenomenon --sector-id`); NULL for one generated standalone. |
| `name` | TEXT | NOT NULL | |
| `remnant_class` | CHAR(1) | NOT NULL | Added in v38: letter class R-W (`program_constants.NEBULA_CLASSES`). |
| `morphology` | TEXT | NOT NULL, CHECK IN ('shell','plerion','composite') | Fixed by the class since v38. |
| `dominant_species`, `density_cm3`, `temperature_k`, `extinction_av` | VARCHAR(255) / DOUBLE | NOT NULL | Added in v38, as on `nebulae`. |
| `age_years` | DOUBLE | NOT NULL | |
| `radius_ly` | DOUBLE | NOT NULL | Derived from `age_years` via the Sedov-Taylor blast-wave relation (`radius ∝ age^(2/5)`). |
| `progenitor_type` | TEXT | NOT NULL, CHECK IN ('Type Ia','core-collapse') | |
| `compact_remnant_kind` | TEXT | nullable, CHECK IN ('black_hole','neutron_star') | NULL when no compact remnant is embedded. |
| `compact_remnant_black_hole_id` | INTEGER | FK -> `black_holes.id`, `ON DELETE SET NULL`, nullable | |
| `compact_remnant_neutron_star_id` | INTEGER | FK -> `neutron_stars.id`, `ON DELETE SET NULL`, nullable | |
| `galactic_orbital_speed_kms`, `_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | NOT NULL | Added in v17. Always populated (a supernova remnant is always standalone). |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v28: this remnant's own galaxy-frame center, the same shape `nebulae` uses. NULL together: never placed in the galaxy. |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |

### `rogue_planets`

Added in v16. Always standalone — a rogue planet is by definition unbound
from any star.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE CASCADE`, nullable | The sector it was generated as part of (`generate.py sector`) or placed near (`generate.py phenomenon --sector-id`); NULL for one generated standalone. |
| `name` | TEXT | NOT NULL | |
| `planet_type` | TEXT | NOT NULL, CHECK IN ('t','g') | Same letters as `planets.body_type`. |
| `planet_class` | VARCHAR(4) | nullable | Added in v47: its `PLANET_CLASSES` letter (GEN.8); NULL for a brown dwarf. |
| `mass_bin` | VARCHAR(16) | NOT NULL | Added in v37: `terrestrial`, `sub-neptune`, `saturn`, `jupiter` or `brown-dwarf`. |
| `mass_kg`, `radius_km` | DOUBLE | NOT NULL | |
| `composition` | TEXT | NOT NULL | Descriptive bulk-composition string. |
| `has_internal_heat`, `has_moons` | BOOLEAN | NOT NULL | `has_internal_heat`: still geologically active (computed since v48; a giant always is). |
| `age_gy`, `internal_heat_flux_w_m2`, `effective_temperature_k` | DOUBLE | nullable | Added in v48: its age and its own heat (W/m^2, K). |
| `surface_regime` | VARCHAR(24) | nullable | Added in v48: `bare-rock`, `frozen-atmosphere`, `ice-shell-ocean`, `ice-world`, `hydrogen-envelope`, `gas-giant` or `brown-dwarf`. |
| `surface_temperature_k`, `surface_pressure_pa` | DOUBLE | nullable | Added in v48: at the ground, or a giant's at 1 bar. |
| `ice_shell_thickness_km`, `ocean_depth_km`, `has_liquid_water` | DOUBLE, DOUBLE, BOOLEAN | nullable | Added in v48: its ice lid and buried (or envelope-warmed) ocean; NULL when there is none. |
| `galactic_orbital_speed_kms`, `_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | NOT NULL | Added in v17. Always populated (a rogue planet is always standalone). |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v28: this rogue planet's own galaxy-frame center, the same shape `nebulae` uses. NULL together: never placed in the galaxy. |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | nullable, `ON DELETE SET NULL` | Added in v39. The innermost cloud holding it, as on `star_systems`. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |

### `quasars`

Added in v31. A galaxy's active nucleus: its central supermassive black
hole, accreting near its Eddington limit. See `schema.sql`'s "v31" header
note.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE CASCADE`, nullable | The core sector that generated it (or `phenomenon --sector-id`, which only accepts a ring-0, layer-0 sector); NULL for one generated standalone. |
| `name` | TEXT | NOT NULL | |
| `black_hole_mass_solar` | DOUBLE | NOT NULL | Log-uniform 1e8-1e10. |
| `event_horizon_radius_km` | DOUBLE | NOT NULL | Schwarzschild radius. |
| `eddington_ratio` | DOUBLE | NOT NULL | Log-uniform 0.1-1. |
| `luminosity_w` | DOUBLE | NOT NULL | `eddington_ratio` x the Eddington limit for its mass. |
| `accretion_rate_solar_per_year` | DOUBLE | NOT NULL | `luminosity_w / (0.1 c^2)`. |
| `broad_line_region_light_days` | DOUBLE | NOT NULL | From the reverberation-mapped radius-luminosity relation. |
| `is_radio_loud` | BOOLEAN | NOT NULL | ~10% launch relativistic jets. |
| `jet_length_ly` | DOUBLE | nullable | Set exactly when `is_radio_loud`. |
| `active_age_years` | DOUBLE | NOT NULL | How long this episode of activity has run. |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Always the galactic center (all 0) when placed; NULL together when never placed. |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `created_at`, `modified_at` | TIMESTAMP | NOT NULL | Row timestamps, as v27 gave every other phenomenon table. |

### `interstellar_comets`

Added in v16. Always standalone — an interstellar comet is by definition
unbound from any star.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE CASCADE`, nullable | The sector it was generated as part of (`generate.py sector`) or placed near (`generate.py phenomenon --sector-id`); NULL for one generated standalone. |
| `name` | TEXT | NOT NULL | |
| `nucleus_diameter_km`, `velocity_kms` | DOUBLE | NOT NULL | `velocity_kms` is hyperbolic excess speed relative to any star it passes — a separate, non-advancing descriptive stat from the `galactic_orbital_*` columns below (its own bulk galactic motion). |
| `is_active` | BOOLEAN | NOT NULL | Whether it currently shows a coma/tail. |
| `composition_summary` | TEXT | NOT NULL | Human-readable summary, same role as `asteroid_belts.composition_summary`. |
| `galactic_orbital_speed_kms`, `_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | NOT NULL | Added in v17. Always populated (an interstellar comet is always standalone). |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v28: this comet's own galaxy-frame center, the same shape `nebulae` uses. NULL together: never placed in the galaxy. |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | nullable, `ON DELETE SET NULL` | Added in v39. The innermost cloud holding it, as on `star_systems`. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |

### `interstellar_comet_composition`

Added in v16. Structured per-component detail behind `composition_summary`
above, mirroring `asteroid_belt_composition` minus a concentration level
(`InterstellarComet`'s own composition list carries no such gradient).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `comet_id` | INTEGER | FK -> `interstellar_comets.id`, `ON DELETE CASCADE`, NOT NULL | |
| `position` | INTEGER | NOT NULL | List order. |
| `component` | TEXT | NOT NULL | e.g. `"water ice"`. |

### `asteroid_fields`

Added in v17, the seventh exotic phenomenon. Always standalone — a field
drifting in open space, as opposed to `asteroid_belts`, which always
orbits a star. `sector_id` was reserved for a future sector-context
encounter through v17, unused (always NULL) by `generate.py phenomenon` until
v18 gave it one.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE SET NULL`, nullable | v18: the nearest already-generated sector to `center_x/y/z_pc` below -- see `nebulae`'s identical "v18" column note above. |
| `name` | TEXT | NOT NULL | |
| `field_class` | VARCHAR(4) | NOT NULL | Added in v38: letter from composition and density plus a size digit, floor(log10(radius in AU)), e.g. `C3` (`program_constants.ASTEROID_FIELD_COMPOSITIONS`). |
| `composition_family` | VARCHAR(16) | NOT NULL | Added in v38: `carbonaceous`, `stony`, `metallic`, `icy`, `basaltic`, `mixed`, `dust` or `collisional`. |
| `density` | TEXT | NOT NULL, CHECK IN ('dense','sparse','typical') | Same three levels `asteroid_belts.density` uses. |
| `radius_ly` | DOUBLE | NOT NULL | |
| `composition_summary` | TEXT | NOT NULL | Human-readable summary, same role as `asteroid_belts.composition_summary` — generated via the same shared `asteroidData.generate_asteroid_composition`/`format_composition_summary` helpers. |
| `galactic_orbital_speed_kms`, `_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | NOT NULL | Always populated (an asteroid field is always standalone). |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v18 -- see `nebulae`'s identical "v18" column note above (an asteroid field's `radius_ly` tops out much smaller, `program_constants.ASTEROID_FIELD_RADIUS_RANGE_LY` = 0.001-1.0 ly, so it usually stays within a single sector, but the same galaxy-frame-sphere model is used for consistency). |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | nullable, `ON DELETE SET NULL` | Added in v39. The innermost cloud holding it, as on `star_systems`. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |

### `asteroid_field_composition`

Added in v17. Structured per-component detail behind `composition_summary`
above, mirroring `asteroid_belt_composition` exactly (same per-component/
concentration shape — both are generated via the same shared
`asteroidData.generate_asteroid_composition`).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `field_id` | INTEGER | FK -> `asteroid_fields.id`, `ON DELETE CASCADE`, NOT NULL | |
| `position` | INTEGER | NOT NULL | List order. |
| `component` | TEXT | NOT NULL | e.g. `"iron"`. |
| `concentration` | TEXT | NOT NULL, CHECK IN ('high','moderate','small','trace') | |

### `sector_name_registry` / `system_name_registry`

Name-uniqueness bookkeeping (v24, `stellarObjects/nameUniqueness.py`).
One row per base name (the name with every decoration this project adds
stripped off) that has collided at least once; a name only ever used once
has no row. `first_*` names the row that first used the base name, the
one renamed as collisions happen (bare, then Alpha, ...).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `base_name` | VARCHAR(255) | NOT NULL, UNIQUE | |
| `occurrence_count` | INT | NOT NULL | How many rows have used this base name. |
| `first_sector_id` | BIGINT UNSIGNED | FK -> `sectors.id`, `ON DELETE CASCADE`, NOT NULL | `sector_name_registry` only. |
| `first_star_system_id` | BIGINT UNSIGNED | FK -> `star_systems.id`, `ON DELETE CASCADE`, nullable | `system_name_registry` only. NULL when the first holder is a phenomenon. |
| `first_object_table`, `first_object_id` | VARCHAR(32), BIGINT UNSIGNED | nullable | `system_name_registry` only (v40): a uniquely named phenomenon holder, e.g. `('nebulae', 12)`. No foreign key. |
| `diminutive_index` | INT | nullable | `system_name_registry` only: the diminutive prefix already applied against a colliding sector name. |

### `nearest_systems`

Added in v41. Up to 3 rows per placed star system or phenomenon: the
nearest star systems, searched across sector boundaries out to
`_db.NEAREST_SYSTEMS_SEARCH_PC`. Filled when a sector is generated (which
also updates its neighbors' lists) and by `updateOrbits.py`. See
"Octants and nearest systems (v41)" above.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `sector_id` | BIGINT UNSIGNED | FK -> `sectors.id`, `ON DELETE CASCADE`, NOT NULL | The object's sector. |
| `object_table`, `object_id` | VARCHAR(32), BIGINT UNSIGNED | NOT NULL | The object: `'star_systems'` or a phenomenon table, and its id. UNIQUE with `neighbor_rank`. |
| `star_system_id` | BIGINT UNSIGNED | FK -> `star_systems.id`, `ON DELETE CASCADE`, nullable | Repeats `object_id` for a system, so deleting it removes its rows. |
| `neighbor_rank` | TINYINT | NOT NULL | 1 for the nearest. |
| `neighbor_system_id` | BIGINT UNSIGNED | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | |
| `distance_pc` | DOUBLE | NOT NULL | |

### `facilities`

Added in v42. Starbases, colonies and outposts, each on exactly one host
named by `host_type` (checked by `_db.add_facility`, since MySQL refuses
a CHECK on a cascading column). See "Facilities (v42)" above and the API's
"Facilities" section ([`api.md`](api.md#facilities)).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `name` | VARCHAR(255) | NOT NULL | Indexed. |
| `kind` | VARCHAR(16) | NOT NULL, CHECK IN ('colony','outpost','mining-colony','station','starbase') | |
| `placement` | VARCHAR(16) | NOT NULL, CHECK IN ('terrestrial','orbital','asteroid','standalone') | |
| `host_type` | VARCHAR(16) | NOT NULL | `star`, `planet`, `moon`, `asteroid_belt`, `asteroid_field` or `space`. |
| `star_system_id` | BIGINT UNSIGNED | FK -> `star_systems.id`, `ON DELETE CASCADE`, nullable | Set for every host inside a system. |
| `star_id`, `planet_id`, `moon_id`, `asteroid_belt_id`, `asteroid_field_id` | BIGINT UNSIGNED | FK, `ON DELETE CASCADE`, nullable | The host row. `star_id` is NULL for a close pair's shared orbit. |
| `sector_id` | BIGINT UNSIGNED | FK -> `sectors.id`, `ON DELETE CASCADE`, nullable | For a stand-alone facility in open space. |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable | A stand-alone facility's galaxy-frame position. |
| `orbit_distance_km`, `orbit_period_years`, `orbital_speed_kms`, `orbit_phase_deg` | DOUBLE | nullable | An orbital facility's circular orbit, from its host's mass; an asteroid facility in a belt gets one too, around its star from a random spot in the belt (no schema change). `updateOrbits.py` advances both. |
| `description` | TEXT | nullable | |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | |

### `bright_stars`

Added in v43. Every star at least `galaxy_shape.bright_star_min_luminosity_sol`
bright, generated and placed galaxy-wide by `generate.py plan`'s scatter
before any sector is filled. See "Bright-star pre-placement (v43)" above.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `ring_index`, `layer_index`, `ring_slot_index` | INT / SMALLINT / INT | NOT NULL | The sector cell it falls in. Indexed together. |
| `position_x_mpc`, `position_y_mpc`, `position_z_mpc` | BIGINT | NOT NULL | Galaxy-frame position, milliparsecs. |
| `population` | VARCHAR(12) | NOT NULL, CHECK IN ('young','intermediate','old','bulge') | |
| `star_type`, `yerkes_class` | VARCHAR | NOT NULL | As on `stars`. |
| `mass_kg`, `radius_km`, `temperature_k`, `luminosity_w`, `age_gy` | DOUBLE | NOT NULL | The finished star. `luminosity_w` is indexed. |
| `lifespan_gy`, `phase_end_age_gy` | DOUBLE | nullable | |
| `initial_mass_sol` | DOUBLE | NOT NULL | A companion's mass is drawn from it. |
| `seed` | BIGINT UNSIGNED | NOT NULL | Seed for the system built around it later. |
| `star_system_id` | BIGINT UNSIGNED | FK -> `star_systems.id`, `ON DELETE SET NULL`, nullable | Set when its sector is filled. |
| `created_at` | TIMESTAMP | NOT NULL | |

### `bright_star_blocks`

Added in v49 (GEN.23). How deep the bright-star backfill has gone in each
sector block. See "Bright-star backfill per sector block (v49, GEN.23)"
above.

| Column | Type | Null | Notes |
|---|---|---|---|
| `block_ring`, `block_wedge`, `block_slab` | INT / INT / SMALLINT | PK | The level-3 block (`galaxyDrill.DrillBlock(3, ring, wedge, slab)`). |
| `min_luminosity_sol` | DOUBLE | nullable | The dimmest luminosity (L_sun) the block's unfilled sectors hold every star down to. NULL only while a backfill holds the row. |
| `updated_at` | TIMESTAMP | NOT NULL | |

### `species` / `polities` / `system_owners` / `population_state`

Added in v44, filled only by the population pass
(`stellarObjects/population.py`); see "Population and politics (v44)"
above and `docs/design/population-and-politics.md`.

**`species`**: one dominant species per life world.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `name` | VARCHAR(255) | NOT NULL, UNIQUE | |
| `homeworld_planet_id` | BIGINT UNSIGNED | FK -> `planets.id`, `ON DELETE CASCADE`, NOT NULL, UNIQUE | |
| `star_system_id` | BIGINT UNSIGNED | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | |
| `life_chemical` | VARCHAR(64) | nullable | |
| `life_stage` | VARCHAR(32) | NOT NULL, CHECK IN ('multicellularity','technological_civilization') | |
| `build`, `climate`, `size` | VARCHAR(16) | NOT NULL | |
| `civilization_age_years`, `era` | DOUBLE, VARCHAR(16) | nullable | NULL without a civilization. |
| `spacefaring` | TINYINT(1) | NOT NULL, default 0 | Era Interstellar or later. Indexed. |
| `created_at` | TIMESTAMP | NOT NULL | |

**`polities`**: one government per spacefaring species.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `name` | VARCHAR(255) | NOT NULL, UNIQUE | |
| `species_id` | BIGINT UNSIGNED | FK -> `species.id`, `ON DELETE CASCADE`, NOT NULL, UNIQUE | |
| `capital_system_id` | BIGINT UNSIGNED | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | The homeworld's system. |
| `government` | VARCHAR(32) | NOT NULL | |
| `color` | CHAR(7) | NOT NULL | `#rrggbb`, for the map. |
| `reach_ly` | DOUBLE | NOT NULL | How far its territory extends. |
| `created_at` | TIMESTAMP | NOT NULL | |

**`system_owners`**: the polity with the strongest claim on each system
inside some polity's reach, rebuilt from scratch by every pass.

| Column | Type | Null | Notes |
|---|---|---|---|
| `star_system_id` | BIGINT UNSIGNED | PK, FK -> `star_systems.id`, `ON DELETE CASCADE` | |
| `polity_id` | BIGINT UNSIGNED | FK -> `polities.id`, `ON DELETE CASCADE`, NOT NULL | |
| `distance_ly` | DOUBLE | NOT NULL | From the capital. |

**`population_state`**: a singleton row (`id = 1`) holding
`scanned_planet_id`, the highest `planets.id` the pass has scanned.

### `sector_objects` (view, not a table)

`UNION ALL` across `stars`, `planets`, `moons`, `asteroid_belts`, and
`comets` (each joined to `star_systems` for `sector_id`), for "every
stellar object in this sector" queries without hand-writing the union
each time:

```sql
SELECT * FROM sector_objects WHERE sector_id = ?;
```

Columns: `object_type` (`'star'`/`'planet'`/`'moon'`/`'asteroid_belt'`/
`'comet'`), `object_id` (the row's real id in its own table),
`star_system_id`, `sector_id`, `name` (the literal `'Asteroid Belt'` for
a belt), `summary` (a short type-appropriate label — a star's
`star_type`, a planet's/moon's class or body type, a belt's density, a
comet's `orbit_type`), `orbital_index` (NULL for stars and comets).

This is a view rather than a sixth physical table specifically to avoid
write-side upkeep: a real table would need to be kept in sync on every
insert/update/delete to the tables it mirrors, or drift out of sync. A view
has no storage and resolves against current data on every query.

## Relationships at a glance

```
sectors 1───* star_systems *───1 system_configs 1───* system_config_slots
                  │
                  ├──1───* stars
                  │
                  ├──1───* planets
                  │             │
                  │             ├──1───* planet_evolutionary_paragraphs
                  │             ├──1───* planet_reflection_spectrum
                  │             └──1───* moons
                  │                       │
                  │                       ├──1───* moon_evolutionary_paragraphs
                  │                       └──1───* moon_reflection_spectrum
                  │
                  └──1───* asteroid_belts 1───* asteroid_belt_composition

planets.star_id ──────────> stars.id   (nullable; NULL for a 'close' binary's planets, set for a single star's or a 'wide' binary's)
moons.star_id   ──────────> stars.id   (nullable; same convention as planets.star_id)
asteroid_belts.star_id ───> stars.id   (nullable; same convention as planets.star_id, added in v15)
moons.star_system_id ─────> star_systems.id   (redundant with moons.planet_id's own owner)
```

`black_holes`/`neutron_stars` (v16) also stand apart from the tree above
when generated standalone (`star_id` NULL) — only when `star_id` is set do
they hang off a `stars` row the same way `asteroid_belt_composition` hangs
off `asteroid_belts`. `supernova_remnants`/`rogue_planets`/
`interstellar_comets` (v16/v17) stand apart entirely — always standalone,
with only a reserved, currently-unused `sector_id` FK. `nebulae`/
`asteroid_fields` (v16/v17) are the same shape but, as of v18, their
`sector_id` is a real (if non-authoritative) link to the nearest
already-generated sector, alongside their own galaxy-frame
`center_x/y/z_pc`/`galactic_radius_pc` — see this file's "v18" note above.

`facilities` hangs off whichever host it names (a star, planet, moon,
belt, asteroid field, or a sector for open space), and `species` off a
planet, with `polities` and `system_owners` hanging off `species` and
`star_systems`. `nearest_systems` links systems and phenomena to their
nearest systems; `bright_stars` points at the system built around each
star once its sector is filled.

`galaxy_shape`/`galaxy_layer`/`galaxy_column` stand apart from the tree above —
neither has a foreign key to `sectors` or anything else. They describe the
galaxy as a whole (its shape parameters, and which addresses could hold
content), not any individual sector; `sectors.ring_index` is the ring
both an actual generated sector and a `galaxy_layer` row are
independently expressed in, not an FK
relationship.
