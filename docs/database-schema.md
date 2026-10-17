# planetGen Database Format

This document describes the MySQL database schema defined in
[`src/planetgen/db/schema.sql`](../src/planetgen/db/schema.sql). It's the
reference for anyone reading, querying, or extending the database — table by
table, every column's meaning and unit, and the conventions that hold the
schema together.

**Status**: the schema is implemented, with both a write and a read path.
[`src/planetgen/db/store.py`](../src/planetgen/db/store.py) (private — leading
underscore, not part of the package's public generation API) writes
already-generated `StarSystem`/`SpaceSector` objects straight into these
tables (`planetgen sector` calls it automatically on every run), and
reconstructs them back into live objects from rows (`load_star_system`/
`load_sector`/`load_system_config`), inverting every unit conversion the
write path applies. `planetgen.db.query` is the "list what's stored" CLI (`TODO.md`
Phase 3). The database itself lives on a MySQL server (TODO.md Phase 5 —
this schema previously targeted SQLite; see "MySQL port" in `schema.sql`'s
own header comment for the type-mapping/idempotency notes that move brought),
reachable via the `$PLANETGEN_MYSQL_HOST`/`$PLANETGEN_MYSQL_PORT`/
`$PLANETGEN_MYSQL_USER`/`$PLANETGEN_MYSQL_PASSWORD`/`$PLANETGEN_MYSQL_DATABASE`
environment variables (or the equivalent `--mysql-*` CLI flags every entry
point accepts — see `planetgen.db.store.MySQLConfig`), with its tables created
automatically on first connection. See `TODO.md` for the full roadmap.

## Persistence layer

`src/planetgen/db/store.py` (imported as `_db` in the older code and docs below) owns every unit conversion at the point of writing
a value into the database; nothing else in the package imports it, and it
never mutates the generation/physics code's own native units. Its shape:

- `get_connection(config=None)` / `save_sector(sector, config=None)` —
  the two entry points most callers need. `save_sector` connects (pooled —
  see `MySQLConfig`/`_get_pool`), applies the schema if needed (idempotent —
  `CREATE TABLE IF NOT EXISTS`/`CREATE OR REPLACE VIEW` throughout
  `schema.sql`, which also adds any table an older database lacks), and writes the whole sector in one transaction.
  `config` (a `MySQLConfig`) defaults to `DEFAULT_MYSQL_CONFIG`, itself
  built from the `$PLANETGEN_MYSQL_*` environment variables above.
- `insert_sector` / `insert_star_system` / `insert_star` / `insert_planet`
  (calls `insert_moon` for each of a planet's moons) / `insert_moon` /
  `insert_asteroid_belt` / `insert_system_config` — the per-table building
  blocks, callable individually for a standalone system with no sector
  (`insert_star_system(conn, star_system, system_config, sector_id=None,
  position=None)`).
- `migrate_database(config=None)` — brings a database up to
  `SCHEMA_VERSION`, the number of the newest Alembic revision in
  `src/planetgen/db/migrations/versions/` (see "Versioning" below). A
  brand-new database is created at the current schema by `schema.sql`, so
  this is a no-op there; a database made by an older release gets one
  revision per version it is behind, and one older than the baseline (v61)
  is refused with `SchemaTooOldError`. `planetgen.cli.migrate` is its CLI
  (`install.sh` and `update.sh` run it on every deploy).
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
  (`planetgen/web/lib/tabledisplay.py`, `planetgen/db/render.py`).

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
`src/planetgen/generation/star.py`, `doubleStar.py`, `planetData.py`,
`asteroidData.py`, and `spaceSector.py` keeps its own native unit (km, AU,
or ly) exactly as today. Conversion only happens at the persistence
boundary, once it's built: `src/planetgen/physics/units.py` provides
`ly_to_milliparsecs`/`milliparsecs_to_ly` for the sector-scale columns;
AU-to-km needs no helper, since it's a single multiply by the existing
`physical_constants.AU_TO_KM`.

The `table_*` columns that once held copies of already-formatted display
text (e.g. `"1.2 R☉"`) were the one exception to both conventions; v5
dropped them (see "The searchable-field principle" below).

### The searchable-field principle

*Historical: v5 dropped every `table_*`/`binary_table_*` column described
in the next paragraph, and display strings are now computed on demand
from the raw columns (`planetgen/web/lib/tabledisplay.py`). The belt columns in
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

### The stored point of an object

Every object is one `SpatialPosition3D` (`planetgen/physics/position.py`): a position, a point mass and a movement vector at an epoch (GEN.294). `planetgen/db/point_columns.py` (GEN.335) is the one place that names the columns that hold those four parts, in a `PointLayout` per table: the position columns with their frame (galaxy, sector or system) and unit, the velocity columns (km/s), the mass column and its unit (or the table it is summed from, as a star system's is from its stars), the epoch and the due time. `planets` and `moons` (system frame, km, `mass_kg`), `bright_stars` (galaxy frame, mpc, `mass_kg`), `phenomenon_scatter` (galaxy frame, mpc, `mass_solar`) and `star_systems` (sector frame, mpc, with the galactic pc beside it, mass from `stars`) have all four. The SQL that reads and writes their point takes its column lists from the layout, and `spatial_from_row` / `values_from_spatial` convert a row to a `SpatialPosition3D` and back (a system-frame position is kept exactly as stored, not as a galactic sum).

The other tables that hold an object are on the layout's allow-list, each with the parts it lacks and the item that adds them (GEN.336 for stars, black holes, neutron stars, rogue planets, interstellar comets, asteroid fields, nebulae and remnants; GEN.337 for belts, the quasar, comets, facilities and scatter members). `test_point_columns.py` lists every table that holds an object, fails when one is on neither list, and fails when an allow-list entry names a part its table already has, so the list only shrinks and a new kind of object cannot ship without a point, a mass and a vector.

### Versioning

Two independent version numbers:

- `schema_migrations` (one row per applied DDL migration step) — the DDL
  structure version, `MAX(version)` in that table (this schema is version
  `91`, `store.SCHEMA_VERSION`, which is the number of the newest Alembic
  revision). `alembic_version` holds the same number as four digits
  (`0091`), kept in step by `alembic_runner.upgrade`. Replaces SQLite's
  `PRAGMA user_version`, which has no MySQL equivalent — see `schema.sql`'s
  "MySQL port" header note.
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
demand from the underlying data columns, e.g. `planetgen/web/lib/tabledisplay.py`);
v6/v7 gave every galaxy-placed sector exact vertices (`sector_vertices`,
built from an exact local spherical Voronoi tessellation among its
same-shell neighbors — the cylinder-era `galaxy/geometry.py`, removed in TEST.156); v8 added
the galaxy-wide density "skeleton" (`galaxy_shape`/`galaxy_layer`,
built by `planetgen plan` — see `planetgen/galaxy/density.py`/
`galaxySkeleton.py`) plus a `UNIQUE (shell_index, shell_slot_index)`
constraint on `sectors`, turning a concurrent lazy-generation race
(`generate.ensure_sector_generated`) into a recoverable `IntegrityError`
instead of a silent duplicate row; v9 added orbital motion —
`orbital_inclination_deg`/`orbital_ascending_node_deg`/
`orbital_phase_deg`/`rotation_period_hours` on `planets`/`moons`, plus the
`orbit_simulation_state` singleton row `planetgen.cli.orbits` uses to track
elapsed time between runs (see that table's own section above); v10 added
`galactic_orbital_speed_kms`/`galactic_orbital_period_gy` to `stars` (and
their `binary_galactic_orbital_*` counterparts on `star_systems`) — a star
system's circular orbital speed/period around the galactic center, from a
simple rotation-curve model (see `planetgen/physics/constants.py`'s
`GALACTIC_ROTATION_FLAT_VELOCITY_KMS` comment); v11 added
`position_x_km`/`_y_km`/`_z_km`/`orbital_speed_kms` to `planets`/`moons` —
each body's Cartesian position relative to its orbital anchor (the star,
or a binary's combined center, for a planet; the parent planet for a
moon), derived from `distance_km` and the v9 orbital-motion columns (see
`planetgen/physics/orbits.py`'s `orbital_position_au`); v12 added
`planets`/`moons.min_update_interval_years` — a floating-point update
guard, not a narrative stat: the shortest `elapsed_years` worth calling
`_db.advance_orbital_phases` for, below which the phase delta added is
smaller than `orbital_phase_deg`'s own IEEE 754 double-precision
resolution and so is guaranteed to round back to the exact value already
stored (see `planetgen/physics/orbits.py`'s `minimum_update_interval_years`).
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
`orbits.circular_orbital_speed_kms`, the same Kepler/circular-orbit formulas
a planet's orbit around its star already uses, applied to the pair's
`binary_separation_km`/`binary_effective_mass_kg`), `_inclination_deg`/
`_ascending_node_deg`/`_phase_deg` (the same `orbits.orbital_position_au`
orbital-element convention planets/moons use, drawn from the full
`[0, 180)`/`[0, 360)` range with no small-tilt bias — a binary's mutual
orbital plane has no protoplanetary-disk reason to prefer any alignment,
unlike a planet's), and `binary_mutual_min_update_interval_years`. v14
added `binary_mutual_position_x_km`/`_y_km`/`_z_km` — the secondary's
Cartesian position relative to the primary, derived from
`binary_separation_km` and the v13 `binary_mutual_orbital_*` columns via
`orbits.orbital_position_au`, the same "position relative to whatever this
orbit is around" convention `planets`/`moons.position_x/y/z_km` already
use (see v11 above) — recomputed by `_db.advance_orbital_phases` in
lockstep every time `binary_mutual_orbital_phase_deg` advances.

v15 added S-type (wide) binary support, planetGen's second real binary
configuration alongside the P-type/close pair every `binary_*` column
above already covers (see `planetgen/generation/wide_binary.py`'s module
docstring for the physics). `star_systems` gains `binary_configuration`
(`'close'` | `'wide'` | NULL for a single star — the authoritative
discriminator going forward; `is_binary` is kept and now means "this
system has two stars", true for either configuration) plus
`binary_eccentricity`/`binary_periapsis_km`/`binary_apoapsis_km` (0/
`binary_separation_km` for a `'close'` pair, real values for a `'wide'`
one — see `orbits.holman_wiegert_critical_semimajor_axis`). The existing
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

v16 added six tables for `planetgen phenomenon`'s separate, rarer exotic-
phenomenon generation mode (`planetgen/generation/phenomena/compact_remnant.py`/
`nebulaData.py`/`supernovaRemnantData.py`/`roguePlanetData.py`) — never
produced by `planetgen system`/`planetgen sector`'s normal generation odds, so no
pre-existing table's shape changes. `black_holes`/`neutron_stars` are
satellite tables extending a `stars` row (nullable `star_id`) when a
compact remnant anchors a full `StarSystem` (`planetgen
phenomenon --anchor-system` — the owning `stars.yerkes_class` is then the literal
marker `'BH'`/`'NS'` rather than a real Yerkes class); `star_id` is NULL
for a remnant generated standalone. `nebulae`/`supernova_remnants`/
`rogue_planets`/`interstellar_comets` are always standalone, with a
nullable `sector_id` reserved for a future sector-context encounter
(unused by `planetgen phenomenon` today). `supernova_remnants` references at
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
has (new shared `galactic_orbit.generate_galactic_orbit_fields`/
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
(`planetgen phenomenon --sector-id`) with the *nearest* already-generated
sector to the phenomenon's own center — a convenience "home" link for
browsing (`ON DELETE SET NULL` now, not `CASCADE`: deleting that sector
shouldn't delete a phenomenon merely linked to it), not the authoritative
geometry (`queryDb.phenomena_near_sector` finds every phenomenon whose
sphere overlaps a given sector's cube by real distance, regardless of
which sector it's linked to). `planetgen.db.store.compute_phenomenon_placement`
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
point (see `orbits.calculate_reflex_offset`): `stars`/`planets` each gain
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
all seven types reachable only through `planetgen phenomenon`'s separate,
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

**How a database is brought up to date (DB.11).** Every schema change is an
Alembic revision (`src/planetgen/db/migrations/versions/`, one file per
version: `0061_baseline.py` to `0091_spherical_sectors.py`), numbered by the
schema version it produces. `schema.sql` is kept at the current shape, so a
new database is created from it directly and recorded at `SCHEMA_VERSION`; an
older one gets its missing tables from `schema.sql` and then each revision
after its recorded version in order (`alembic_runner.upgrade`, run by
`store.migrate_database`), under a schema lock so a first connection
elsewhere waits instead of stamping the database current half way. The
per-version `_migrate_vN_to_vN+1` functions the entries below mention
(v8 to v60) are gone: v61 is the baseline, and a database older than that is
refused (`SchemaTooOldError`) and has to be generated again. A database newer
than the code is refused too (`SchemaTooNewError`). Every released schema from v61
is kept as a fixture (`src/tests/fixtures/old_schemas/`) and
`tests/test_db_old_schemas.py` migrates each to the current shape and checks
it against `db/models.py` (the SQLAlchemy models that mirror `schema.sql`).
`planetgen check-db` (`db/check.py`) compares the live tables against those
models and reports any difference. A SQLite database from before the MySQL
port can't be brought in (SQLite stopped at v5); generate the galaxy again.

The entries from "Schema history" down describe each version as it was when
it was built. In particular the cylindrical sector grid they mention (rings,
layers and slots; `galaxy_layer`, `galaxy_column`, the block ladder) was
replaced by shells, bands and wedges in v91, which the `sectors` section and
the v91 entry at the end of this document describe.

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
ticks from `planetgen.cli.orbits` deliberately don't count as a change; see
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
See `planetgen/names/bodies.py`.

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

**Population and politics (v44).** Filled by `planetgen population`
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
placed galaxy-wide by `planetgen plan` before any sector is filled: its
cell (`ring_index`, `layer_index`, `ring_slot_index`), galaxy-frame
position in milliparsecs, population (young, intermediate, old,
bulge, halo or cluster), the finished star's properties (including `initial_mass_sol`,
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
v53 replaced the table with per-sector levels in `sector_stats`.

**What made the galaxy, and every run (v52).** `_migrate_v51_to_v52`
adds `galaxy_shape`'s version columns (NULL on an existing galaxy) and
the empty `generation_runs` table (DB.6).

**Velocity (v60, GEN.121).** `planets`, `moons` and `comets` each gain
`velocity_x_kms`, `velocity_y_kms` and `velocity_z_kms`: the body's velocity,
km/s, relative to its orbital anchor and on the same axes as its
`position_x/y/z_km`. `advance_orbital_phases` (planets and moons, from the
new phase) and `advance_comet_orbits` (Kepler or Barker) move it with the
position. `_migrate_v59_to_v60` works it out for rows saved before.

**A system's galactic velocity (v61, GEN.121).** `star_systems` gains
`velocity_x_kms`, `velocity_y_kms` and `velocity_z_kms`, km/s on the galactic
axes: the rotation curve's tangent at the system's place
(`galactic_orbital_speed_kms`) plus, for a runaway or hypervelocity system,
`runaway_speed_kms` along a random direction (`StarSystem.runaway_direction`,
drawn by `flag_fast_stars`). `advance_galactic_positions` turns it with the
position and a loaded sector puts it on the entry, its stars and its bodies
(`system_position.place_system`). `_migrate_v60_to_v61` works it out for rows
saved before, drawing a runaway's direction from its id. A placed phenomenon
stores none: it only follows the rotation curve, which its
`galactic_orbital_speed_kms` gives.

**A nebula's shape (v56, GEN.75).** `nebulae.shape_*` hold the single
values of its shape (`galaxy/nebula_shape.py`: the ellipsoid's stretch
and turn, the noise warp, the isovalue and the scale that puts the
surface's farthest point at the radius) and `nebula_shape_balls` its 4-8
metaballs, one row each (`ON DELETE CASCADE`). `_migrate_v55_to_v56`
leaves them NULL for a nebula already saved; `queryDb.nebula_shape` then
draws the same shape from the nebula's own properties.

**Per-sector stats (v53, GEN.44, PERF.11).** `sector_stats` holds one
row per sector address a backfill reached or a fill generated (an
unfilled sector has no `sectors` row, so the key is the address). Its
`bright_level_sol` replaces v49's per-block levels: -1 untouched (the
sector follows the galaxy scatter's threshold), a positive L_sun for the
dimmest a backfill drew every star down to, 0 once the sector is
generated. The backfill (`generate.backfill_bright_stars`) works sector
by sector: each sector in range gets the stars from its own distance
tier's floor up to its level, a few hundred sectors per transaction under
their rows' locks (`_db.lock_sector_stats`), and its new level is written
with its stars, so a crash leaves either the old level or none. A sector
holding unbuilt stars with no level anywhere (no row level and no
finished scatter) is left over from a failed run: the backfill deletes
them and draws it whole. A band run (`plan --bright-stars-down-to`)
leaves sectors with their own level out of its layers and tops each up
from its own level. A sector's draws come in fixed luminosity bands
(`brightStars.canonical_bands`, eight to a decade), each from its own
stream, so down to 1000 L_sun and later to 500 gives exactly the stars of
one draw to 500. A fill (`_db.record_sector_stats`) sets level 0, keeps
the level it found in `level_before_fill_sol` (deleting the sector puts
it back), and records the sector's expected density and what it got.
`galaxy_shape.density_ratio_avg` is a decaying average of actual against
expected systems over every fill. `_migrate_v52_to_v53` creates the
table, gives every filled grid sector a row at level 0, moves each
`bright_star_blocks` level onto its block's sectors, and drops that
table. Density stats start empty for sectors filled before.

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

**Object IDs (GEN.64).** (A star system built around a bright-sweep star is named like any system once its sector is generated, and keeps its position ID as its unique ID: GEN.72.) Every placed interstellar object (a rogue
planet, standalone black hole or neutron star, nebula, supernova remnant
and its collapsed core, quasar, interstellar comet or asteroid field)
is named by its 76-bit position ID, 19 hex digits, instead of the rules below:
`planetgen/names/object_id.py` packs type, distance unit, distance,
bearing, mark and a 4-bit collision number from the galactic center,
and `_db._claim_object_ids` hands out collision numbers in generation
order. These names never go through `system_name_registry`. A name given
by hand and an unplaced object still follow the rules below, and rows
saved before this keep their names. See `docs/design/object-ids.md`.

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
every sector it reaches, when a nebula or remnant is placed. Since GEN.75 a
nebula holds a point only inside its stored shape (the sphere is the quick
first test, then `NebulaShape.contains`); a remnant is still a sphere.
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
MySQL schema entirely (`planetgen/db/control_schema.sql`,
`control_schema_migrations`, currently version 14) — see "The control
schema" below. `SCHEMA_VERSION`/`schema_migrations` above only ever
describe the per-galaxy content schema this whole document is otherwise
about.

## The control schema

A second, deployment-global MySQL schema (`PLANETGEN_CONTROL_DATABASE`,
default `planetgen_control`) holds everything about *who can administer
this deployment*, separate from every per-galaxy content schema this
document otherwise describes — see `planetgen/db/control_schema.sql`'s
header comment for the full rationale (in short: a deployment can host
several galaxy databases sharing one MySQL server, and admin identities
describe the deployment, not any one galaxy, so they aren't duplicated
into each content schema's `schema.sql`).

Eighteen tables, versioned independently via `control_schema_migrations`
(currently version 14, `store.CONTROL_SCHEMA_VERSION`, mirroring `schema_migrations`'s own shape; v2 added
`login_throttle`, which v8 dropped (SEC.30: the login lockouts live in Redis), v3 `admin_devices`, v4 `admin_totp` and
`admin_recovery_codes`, v5 the work queue's three tables, v6
`generation_stats` and `generation_size`, v7 the work queue's job
tree and pause columns, v9 `galaxy_naming`, v10 `work_job_args`, v11
`version_key_history`, v12 the `version_key` column of `generation_stats` and v14
`generation_stage_runs`; the new tables are added to an older schema by
`CREATE TABLE IF NOT EXISTS` on the next `planetgen.cli.migrate` run, and v7 and v12 add columns, see below):

- **`admin_users`** — one row per admin (`username`, `password_hash`,
  `must_change_credentials`). No roles/permissions column — every admin
  has the same full access (see `docs/TODO.md`; this project deliberately
  has no general user-accounts system, just a handful of admins).
- **`admin_sessions`** — web-UI login sessions (`token_hash`, a SHA-256
  digest of the actual cookie value — never the raw token itself;
  `expires_at`, a fixed lifetime set at creation, no sliding renewal).
- **`admin_api_keys`** — API keys for programmatic callers (`key_hash`,
  same "hash only, never the raw key" treatment; `revoked_at` rather than
  a hard delete, so a revoked key's history stays visible; `key_prefix`,
  its first eight characters, tells keys apart without the key;
  `expires_at` makes a key stop working, like a revoked one — control
  schema v13, API.9).
- **`admin_api_key_scopes`** — what each key may do: one row per
  `(key_id, scope)` with `scope` one of `read`, `generate`, `upload`,
  `admin` (`admin` implies all; `generate` and `upload` imply `read`).
  Rows go with their key (`ON DELETE CASCADE`); keys made before v13 were
  given `admin`. No CHECK on `scope`, so a new scope needs no migration.
- **`admin_audit_log`** — one row per write/admin action (`admin_user_id`
  + a denormalized `admin_username` snapshot, `action`, `target`,
  `detail`, `created_at`) — written by `planetgen/api/routes.py`'s write
  routes (`planetgen/api/authz.audit`) after each one actually succeeds,
  plus one per refused sign-in (`login.failed`, `login.locked`,
  `password.failed`, with `target` `ip:<address>`; those are deleted
  after 90 days).
- **`admin_devices`** (v3, SEC.22) — trusted browsers: `admin_user_id`,
  `token_hash` (SHA-256 of the `pg_admin_device` cookie, never the raw
  value), `expires_at` (90 days after the login that made it),
  `last_used_at`. A login from a browser holding a valid one for that
  username skips the per-username lock. Deleted on a credentials change,
  with `python3 -m planetgen.cli.lockouts --forget-devices <user>`, and (expired ones)
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
  generation work queue (`planetgen/queue/work.py`, see
  [`cli.md`](cli.md#parallel-generation)). `work_jobs` is one row per
  run that used worker processes (`title`, `holder` host:pid:token,
  `state` waiting/running/done/failed/cancelled, `workers`, task
  counts, timings, `heartbeat_at`); `work_tasks` one per sector it
  handed out (`kind`, `task_key` such as "shell,band,wedge", `state`,
  the worker's `seconds`, a short JSON `result` or the `error`);
  `work_lease` the single row naming the run whose workers are using
  the machine, refreshed every 5 s. A lease or job not refreshed for
  30 s belongs to a run that died: the next run takes the lease and
  marks that run's unfinished tasks cancelled. Jobs older than a week
  are deleted with their tasks.
- **Job tree columns** (v7, ADM.12) — every job is a tree. `work_jobs`
  gained `parent_id` (the node above, deleting a node deletes its
  subtree), `root_id` (the top of its tree), `kind` (`web-job` for a
  Generate page job, `step` for one of its steps, a `planetgen`
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
- **`galaxy_naming`** (v9, GEN.70) — one row per galaxy database:
  `naming_key` (8 uppercase hex digits), drawn from the galaxy seed when
  the galaxy is planned (`planetgen plan` draws a new one for a new seed
  and keeps an admin's over the same seed) and changeable by an admin
  (Stats page, `POST /api/admin/naming-key`); `codec_version` (the name
  codec's version when the key was set), `drawn_at`, and `changed_at` /
  `changed_by` for the last admin change. The names of the objects the
  phoneme codec names are computed from their ID and this key
  (`planetgen/names/naming_key.py`), so changing it rewrites no row. A
  galaxy reset keeps the row until the next plan draws a new key.
- **`work_job_args`** (v10, DB.13) — a job's command line, one row per
  argument in order (`job_id`, `position`, `value`), without the `--mysql-*`
  options. It replaced the JSON text column `work_jobs.argv`.
- **`version_key_history`** (v11, OPS.13) — the version key each update ran
  under, per galaxy database, newest last: `database_name`, `galaxy_seed`
  (32 hex digits), `version_key` (22), `planetgen_version`,
  `requirements_sha256` (NULL with no `requirements.lock`), `recorded_at`.
  Every update adds a row per galaxy and keeps the last 10 per database
  (`planetgen/galaxy/version_history.py`), so an admin can see which
  release and environment a seed was last known to make the same galaxy under.
- **`generation_stage_runs`** (v14, PERF.56) — how long each stage of a
  generation run took, with the settings it ran with
  (`planetgen/generation/stages.py` writes a row as a stage ends):
  `database_name`, `command` (`galaxy` or `plan`), `stage_key`, `stage_n` of
  `stage_total`, `label`, `skipped` and `skip_reason` (a skipped stage gets a row
  with no time), `started_at`, `finished_at`, `seconds`, `workers`, `settings`
  and `metrics` (JSON text: the arguments that shape the stage, and what it did)
  and `version_key`. Estimates use only rows of the running version key.
  Shown on the Stats page.
- **`generation_stats`**, **`generation_size`** (v6, PERF.3, PERF.10) —
  how fast this server generates and how big a galaxy gets
  (`planetgen/generation/stats.py`, see
  [`cli.md`](cli.md#size-and-time-estimates)). `generation_stats` has
  one row per `kind` (`sector` fill, or a `scatter` layer of bright
  stars) and density `bucket` (`floor(2 * log10(density / 0.01))`: two
  per decade from 0.01, with no top edge, `density_low`/`density_high`
  its edges), each column a decaying average over every task that ever
  finished in it with that many `workers` (each new task weighs 5%; `version_key`
  names the release, Python and platform that measured it, v12, PERF.32: rows
  of any other key are deleted by the next run, so a new release starts from
  nothing, and an admin can delete them all; a benchmark's rows use a `bench:`
  `kind` prefix and never feed a live estimate): `seconds_per_task` (a
  worker's wall time), `seconds_per_system`, `systems_per_task`,
  `stars_per_system`, plus `samples` and the densest task seen
  (`max_density`). `generation_size` has one row per galaxy database:
  `bytes_per_system` (its tables' data and index bytes over its star
  systems, measured after each run), `systems`, `total_bytes`. A galaxy
  reset keeps both. These and the work queue's tables are the only
  control tables the generator writes; it reaches them with its own
  MySQL account and runs without them when it can't.

`planetgen/admin/auth.py` is the only code that reads/writes the
admin tables directly — `bootstrap_control_schema` creates the schema and,
when `admin_users` is empty, seeds an `admin` row with a random first
password that `planetgen.cli.migrate` prints once (`planetgen.cli.migrate` calls this
automatically, alongside its usual content-schema migration), and every
other function there implements one piece of the login/session/API-key/
audit lifecycle `planetgen/api/auth.py`'s routes expose.

### Booleans and tri-state flags

MySQL has no dedicated boolean type either. Plain booleans are `TINYINT(1)`
`0`/`1` with a `CHECK` constraint. `SystemConfig`'s tri-state flags (`True`/
`False`/`None` in Python — force-on / force-off / random) are nullable
`INTEGER` `0`/`1`/`NULL`.

### Rendered wiki text and URLs

No rendered page text is stored (v29 dropped `star_systems.wikitext_content`/
`markdown_content`). Every value the page shows has its own column, so
`_db.load_star_system` rebuilds the generation object graph and
`planetgen/db/render.py` renders either format on demand — the
same `StarSystem.__str__` generation used to run once before saving, so
the text is identical, except that it now follows later changes (renames,
names made unique, orbit ticks' current wobble and comet positions).
`planetgen.cli.render_parity` compares the two on a database still at v28.
`mediawiki_url`/`wikijs_url` (v23 — see `schema.sql`'s header comment)
record where that page lives on each wiki, once `POST
/api/systems/<id>/wiki` (`src/planetgen/wiki/`, `web/system_pages.py`'s "Upload to
Wiki" form) has actually uploaded it there — one system is one wiki page;
individual stars/planets/moons are sections within that one page, not
separate pages. Existence on a given wiki is never a separate stored
flag — it's exactly "the matching URL column is not NULL"; `web/system_pages.py`
links to that page (opening in a new tab) whenever either is set. `sectors.wiki_url` is the
per-sector equivalent (see that table's own column doc above) — a single
column, since a sector's page is generated fresh at upload time rather
than persisted the way a system's is.

## Planned changes (not built)

Schema changes the open TODO items plan, as of 2026-10-15 (galaxy
schema v92, control schema v14). None of these exists yet; each lands as
a numbered Alembic revision when its item is built, and this section moves
into the table descriptions below. The galaxy schema's revisions are
written one lane at a time, and the second of two branches to reach main
renumbers its revision ([plan/notes.md](plan/notes.md)). Why and how the
reproducibility pieces fit together is in
[design/reproducible-galaxies.md](design/reproducible-galaxies.md).

Galaxy schema:

| Item | Change |
|---|---|
| DB.9 | Each sector gets a content checksum (the hash of GEN.58's fingerprint); the Reed-Solomon parity itself lives in a file outside the database. |
| DB.30 | One revision for the per-species environment limits. |
| DB.31 | Structure instead of prose for life and flavour text. |
| GEN.109 | A per-sector point-mass table for the orbital update. |
| PERF.29 | Which fill steps each sector has finished, so partial runs resume. |
| API.10, API.11 | Run reservations (claimed sectors and id ranges per run) and staging tables for uploaded batches. |

Control schema:

| Item | Change |
|---|---|
| PERF.100 | `work_jobs`, `work_job_args`, `work_tasks` and `work_lease` are dropped once the work queue lives in Redis only (PERF.86); the Redis keys are not part of this schema. |
| API.15 | Every API call logged: time, route, account (the key's owner, the signed-in admin, or "god" for the console), how it came in (API key, web session or console) and the HTTP response code. Where the rows are kept is decided with that item. |
| USR.2, USR.4, USR.7, USR.11, NAV.19, API.6 | Accounts with roles, invite links, per-account bookmarks and saved courses, a home field, and API keys owned by accounts. |

Everything else the earlier plan listed (the version key and unique IDs
on every object, the point-in-space object, nebula shapes, the
atmosphere, activity, radiation, hydrosphere and habitability columns, the
phenomenon scatter, spin and next-update columns, species tech levels,
facility types, the naming key, and the version-key history) is built:
see "Schema history" and the table descriptions.

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
| `center_x_pc`, `center_y_pc`, `center_z_pc` | DOUBLE | nullable | The sector's center, in a galaxy-frame Cartesian coordinate system whose origin is the galactic center (parsecs — see `docs/design/galaxy-coordinate-system.md`). NULL together iff this sector has never been placed in a galaxy (`planetgen sector`'s own standalone CLI, or a sector migrated from a pre-v4 database). |
| `galactic_radius_pc` | DOUBLE | nullable | `sqrt(x^2+y^2+z^2)`, persisted (not just derivable) so "sectors within radius R of the core" is a plain indexed range scan — same treatment `star_systems.quadrant` gets. NULL iff the center columns are NULL. |
| `shell_index`, `band_index`, `wedge_index` | INT | nullable | This sector's address on the spherical grid (v91, DB.32, `planetgen/galaxy/spherical.py`; see `docs/design/spherical-sector-coordinates.md`): the shell (distance from the galactic centre in 4 pc steps; shell 0 is the core ball), the band (latitude; band 0 holds the plane and a shell's end bands are its polar caps) and the wedge (longitude, counted from the +X meridian). Named `ring_index`, `layer_index` and `ring_slot_index` on the cylindrical grid until v90. Independently nullable from the center/radius columns above (not part of the same CHECK) — a sector could in principle have a hand-authored galaxy position without this particular placement algorithm's own addressing. |
| `wiki_url` | VARCHAR(2048) | nullable | Where this sector's summary page lives on a wiki (v23) — set either by `POST /api/sectors/<id>/wiki` (uploading `planetgen/api/routes.py`'s `_sector_wiki_content`) or directly via `PATCH /api/sectors/<id>` (the `/admin` page's manual-link form, `web/admin_pages.py`). A single column, not one per backend the way `star_systems.mediawiki_url`/`wikijs_url` are — a sector has no persisted rendered page of its own to independently re-upload to a second backend, so only one link is ever tracked at a time. NULL means no page yet. |
| `version_key`, `planetgen_version`, `python_version`, `platform` | CHAR(22), VARCHAR(32), VARCHAR(32), VARCHAR(64) | nullable, `version_key` indexed | Added in v67 (DB.7). The code that generated the sector: the packed version key (`planetgen/galaxy/version_key.py`) and the release, Python and platform strings it stands for, as `galaxy_shape` records them for the plan. NULL on a sector generated before v67. The Generate page and `planetgen galaxy` warn when the running version differs from a sector's (`galaxy/version_check.py`). |
| `created_at` | TIMESTAMP | NOT NULL | Added in v27. See "Row timestamps (v27)" above. |
| `modified_at` | TIMESTAMP(3) | NOT NULL, `ON UPDATE CURRENT_TIMESTAMP(3)` | Added in v27. Indexed. |

A `CHECK` constraint enforces `center_x_pc`/`center_y_pc`/`center_z_pc`/
`galactic_radius_pc` being NULL together (see "v3 → v4" in "Schema
history" above for why that addition needed a wholesale table rewrite
rather than an incremental `ALTER TABLE`). A cell's corners are
closed-form (`spherical.sector_cell_vertices_pc`), so none are
stored (v32 dropped the old `sector_vertices` table). A
`UNIQUE (shell_index, band_index, wedge_index)` constraint
(`uq_sectors_address`) guarantees at most one sector per galaxy address —
NULL-together rows (never placed in a galaxy) don't collide with each
other or with a placed sector, ordinary SQL `NULL` semantics for `UNIQUE`.
This is what turns a lazy-generation race (`generate.ensure_sector_generated`)
into a clear `IntegrityError` its caller recovers from, instead of a silent
duplicate row at the same address.

### `galaxy_shape`

The galaxy-wide density "skeleton" (v8) — a singleton row (`id` pinned to
`1`) holding everything needed to recompute any sector's exact position
and density on demand, built by `planetgen plan`. Deliberately **not** one
row per sector: a sector's position (`planetgen.galaxy.spherical.
sector_position_pc`), density (`planetgen.galaxy.density.
relative_density`), and corners (`spherical.sector_cell_vertices_pc`)
are all pure deterministic functions of its `(shell_index, band_index,
wedge_index)`
address plus this handful of galaxy-wide numbers — cheaper to recompute
on demand (sub-millisecond per sector) than to look up, so none of it is
stored per address. What genuinely needs precomputing — where the galaxy
has any content at all — lives in `galaxy_shell` below instead.
Building or rebuilding the skeleton replaces this row (and every
`galaxy_shell` row) wholesale; there is no partial update, since a
full build takes a few milliseconds even at real Milky-Way scale (a
closed-form walk over ~3,800 rings and ~1,020 layers, not a per-sector scan over
billions of candidates).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK, `CHECK (id = 1)` | Pinned to `1` — there is exactly one galaxy. |
| `disk_scale_length_pc`, `disk_scale_height_pc`, `bulge_scale_radius_pc`, `bulge_amplitude`, `pitch_angle_rad`, `arm_amplitude`, `arm_level`, `nuclear_amplitude`, `spiral_reference_radius_pc`, `spiral_reference_angle_rad`, `k_norm` | DOUBLE | NOT NULL | `planetgen.galaxy.density.GalaxyShape`'s own fields, verbatim — see that module for what each means and how `k_norm` is calibrated. |
| `arm_count` | INT | NOT NULL | Same source. |
| `edge_pc` | DOUBLE | NOT NULL | The sector edge length this skeleton was built at, parsecs: always the standard `program_constants.DEFAULT_SECTOR_EDGE_PC` (4) since v33. |
| `expected_system_count_at_density_1` | DOUBLE | NOT NULL | `SpaceSector(edge_ly=...).expected_system_count()` at `relative_density = 1` — cached since every qualification check needs it. |
| `outer_shell_index` | INT | NOT NULL | The last shell with any qualifying content — this galaxy's real edge, found by `planetgen plan` (the equatorial band reaches farthest), not an arbitrary radius. Named `outer_ring_index` from v32 to v90. |
| `bright_star_min_luminosity_sol` | DOUBLE | nullable | Added in v43. The luminosity threshold the bright-star scatter used (see `bright_stars` below). A sector fill reads this, not the current constant. NULL until a scatter has run. |
| `bright_star_seed` | BIGINT UNSIGNED | nullable | Added in v43. The scatter's seed: since v51 the top 63 bits of the galaxy seed's unit seed `bright-stars:scatter` (GEN.39). |
| `bright_star_mass_limit_sol` | DOUBLE | nullable | Added in v73 (GEN.185). The mass, in solar masses, the star scatter's mass pass placed every star from (the same limit as `phenomenon_min_mass_solar`): a sector's own draw and the luminosity pass are lighter than it. NULL on a scatter drawn before v73, which had no mass pass. |
| `phenomenon_scatter_seed` | BIGINT UNSIGNED | nullable | Added in v64 (GEN.100). The phenomenon scatter's seed (`phenomenon_scatter`); NULL until one has run, when a sector fill rolls its own black holes, neutron stars, planetary nebulae and supernova remnants. |
| `phenomenon_min_mass_solar` | DOUBLE | nullable | Added in v71 (GEN.167). The lowest mass of the neutron stars and black holes the phenomenon scatter placed (`planetgen plan --compact-min-mass`: 1, 2, 4 or 6 solar masses, or the stellar mass limit, GEN.195; before that `--phenomenon-min-mass`). A sector fill draws the ones below it from its own stream (GEN.168). NULL on a scatter drawn before v71, which placed every mass. |
| `galaxy_seed` | BINARY(16) | nullable | Added in v51 (GEN.39). The galaxy's 128-bit seed, shown and typed as 32 hex digits (`planetgen plan --seed`), written by the first plan and kept by every later one. Every sector, bright-star scatter, band and backfill block draws from SHA-256 of it and the unit's `kind:address` (`galaxySeed.unit_seed`). NULL only on a galaxy planned before v51. |
| `version_key` | CHAR(22) | nullable | Added in v52 (DB.6). The version key of the code that wrote `galaxy_seed` (`version_key.version_key`): MAJOR 4, REVISION 4, BUILD 6, Python major, minor and micro 2 each, OS 1 (0 Linux, 1 Windows, 2 other Unix, 3 macOS, F unknown) and architecture 1 (0 x86-64, 1 ARM64, 2 x86, 3 ARM32, 4 RISC-V 64, F unknown) hex digits, for example `0007007F000160030C0300`. Written with a new seed; a later plan keeping the seed keeps it. |
| `planetgen_version`, `python_version`, `platform` | VARCHAR(32) / VARCHAR(16) / VARCHAR(64) | nullable | Added in v52 (DB.6). The same code's full release (`7.127.352`), Python version (`3.12.3`) and `platform.system()` and `.machine()` (`Linux x86_64`). |
| `density_ratio_avg` | DOUBLE | nullable | Added in v53 (PERF.11). Systems a fill got against systems expected (`sector_stats`), as a mean over the first 1,000 fills, then a decaying average weighing each new fill 1 in 1,000 (`_db.DENSITY_RATIO_WINDOW`). NULL before any fill. |
| `density_ratio_samples` | BIGINT UNSIGNED | NOT NULL, default 0 | Added in v53 (PERF.11). How many fills went into `density_ratio_avg`. |

### `galaxy_shell`

The galaxy's outline on the spherical grid (v91, replacing `galaxy_layer` and
`galaxy_column`): one row per shell that can hold content, from the core out
to the last shell whose equatorial band could clear the qualification
threshold (`planetgen.galaxy.skeleton.build_shell_extents`). A shell holds
bands `-high_band_index .. high_band_index` (both polar caps when
`high_band_index` is the shell's top band), every wedge. Like the cylinder's
outline it comes from an exact upper bound, so a sector inside it isn't
guaranteed to qualify, but one outside it is guaranteed not to; the exact
per-sector answer is made when the address is visited
(`generate.ensure_sector_generated`). The outline is symmetric about the
plane. A shell with no row has no content at all.

| Column | Type | Null | Notes |
|---|---|---|---|
| `shell_index` | INT | PK | The shell (0 is the core ball). |
| `high_band_index` | INT | NOT NULL | The highest band the shell reaches; it holds bands `-high_band_index` to `high_band_index`. |
| `cell_count` | BIGINT UNSIGNED | NOT NULL | The cells that is (every wedge of every band), stored so a count of the galaxy's candidate sectors needs no recount. |

Together with `galaxy_shape` this is the galaxy's bounds.
`_db.get_galaxy_bounds` loads it as a `skeleton.GalaxyBounds`,
and every generation path checks an address against it before generating
anything there, so a sector is never placed outside the galaxy.

### `orbit_simulation_state`

Singleton row (same pattern as `galaxy_shape` above) added in v9, tracking
when `planetgen.cli.orbits` last advanced every planet's/moon's
`orbital_phase_deg` in this database — absent entirely until that
script's first run against a given database (it creates this row
itself). `planetgen.db.store.get_orbit_update_elapsed_years` reads it
(via `TIMESTAMPDIFF` against `NOW()`, server-side, rather than trusting
the calling process' own clock) to compute how much simulated time has
passed since the last update.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK, `CHECK (id = 1)` | Always `1` — singleton. |
| `last_updated_at` | TIMESTAMP | NOT NULL | When `planetgen.cli.orbits` last ran against this database. |

**Per-object clocks (v66, GEN.106).** Each object the update moves also
keeps its own clock: `epoch_unix`, when its stored position holds, and an
indexed `next_update_due`, when it will have moved far enough to store
again, both Unix seconds by the database server's clock. The columns are on
`planets`, `moons`, `comets`, `star_systems` (its galactic orbit, plus
`binary_epoch_unix`/`binary_next_update_due` for a pair's mutual orbit),
`facilities`, `black_holes`, `neutron_stars`, `nebulae`,
`supernova_remnants`, `rogue_planets`, `interstellar_comets` and
`asteroid_fields`. The thresholds are 0.01 mpc on a galactic orbit, 0.01 AU
in a system (planets, comets, a pair's mutual orbit, facilities round a
star or in a belt) and 100,000 km round a planet (moons, facilities round
a planet or moon); the due time is the epoch plus threshold / speed, capped
at a billion years (`planetgen.physics.position.update_interval_s`). A run
moves each due object from its own epoch to now, sets its epoch to now and
works out its next due time; an object that isn't due is neither moved nor
counted. NULL `epoch_unix` means the position holds at `last_updated_at`;
NULL `next_update_due` means not worked out yet (a new row, or an orbit
edited since), and the next run fills it in before choosing what is due.
The older `*min_update_interval_years` floating-point guards are no longer
read by the update.

Meant to run on a schedule, not on every deploy (`install.sh`/`update.sh`
don't call it) -- e.g. a monthly cron entry:

```
0 3 1 * * cd /var/lib/planetGen && python3 -m planetgen.cli.orbits >> /var/log/planetgen-orbits.log 2>&1
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

`planetgen.cli.orbits` mutates rows, so it needs the same read-write database
account `planetgen` and the web app use, not a `SELECT`-only one you may
have made for `planetgen.db.query`.

### `system_configs`

One row per `SystemConfig` "recipe" — the generation parameters a
`StarSystem` was built from (`src/planetgen/generation/config.py`,
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
| `prevalence_habitable_world`, `prevalence_asteroid_belt`, `prevalence_comets`, `prevalence_large_star`, `prevalence_moons`, `prevalence_max_planets`, `prevalence_intelligent_life`, `prevalence_binary_system`, `prevalence_wide_binary`, `prevalence_planets` | DOUBLE | NOT NULL, default 0 | Added in v54 (GEN.52). One per forcing option above: a sector or galaxy run's percentage deviation from that feature's normal chance, 0 for none. |

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
| `binary_galactic_orbital_period_gy` | DOUBLE | nullable | Added in v10. Circular orbital speed/period around the galactic center (see `stars.galactic_orbital_speed_kms` below) — independent of mass, so identical to the primary/secondary stars' own values, just mirrored here for the combined-pair row. Its speed column was dropped in v103 (GEN.336): the speed is the length of `velocity_x/y/z_kms` (a star and a close pair read their system's). |
| `binary_galactic_orbital_phase_deg`, `binary_galactic_min_update_interval_years` | DOUBLE | nullable | Added in v13. Same pair as `stars.galactic_orbital_phase_deg`/`galactic_min_update_interval_years` below, mirrored here — always identical to both constituent stars' own values (see "Schema history" above for why). |
| `binary_mutual_orbital_period_years`, `_speed_kms` | DOUBLE | nullable | Added in v13. The pair's own mutual orbit around each other — entirely separate from, and vastly faster than, the galactic orbit above. Kepler's third law / circular-orbit speed (`planetPhysics.calculate_orbital_period_years`/`orbits.circular_orbital_speed_kms`) applied to `binary_separation_km`/`binary_effective_mass_kg`. |
| `binary_mutual_orbital_inclination_deg`, `_ascending_node_deg`, `_phase_deg` | DOUBLE | nullable | Added in v13. Orients the mutual orbit in 3D and tracks the pair's current position within it — same `orbits.orbital_position_au` convention as `planets.orbital_inclination_deg`/etc, but drawn from the full `[0, 180)`/`[0, 360)` range (no small-tilt bias — a binary's mutual orbital plane has no preferred alignment the way a planet's protoplanetary-disk-derived orbit does). `_phase_deg` is advanced by `_db.advance_orbital_phases`, guarded by the interval below. |
| `binary_mutual_min_update_interval_years` | DOUBLE | nullable | Added in v13. Floating-point update guard for `binary_mutual_orbital_phase_deg`, same formula as `planets.min_update_interval_years`. |
| `binary_mutual_position_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v14. The secondary's Cartesian position relative to the primary — same "position relative to whatever this orbit is around" convention as `planets.position_x/y/z_km` (`orbits.orbital_position_au`, applied to `binary_separation_km` and the mutual-orbit orientation columns above). Recomputed by `_db.advance_orbital_phases` in lockstep every time `binary_mutual_orbital_phase_deg` advances. **Unchanged in meaning by v20** — still the true separation, not a barycenter-reduced value. |
| `binary_primary_position_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v20. The primary star's own offset from the pair's barycenter — `-binary_secondary_mass_fraction * binary_mutual_position_*_km`. Reused unchanged for both binary configurations, like `binary_mutual_position_*_km` above. |
| `binary_secondary_position_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v20. The secondary star's own offset from the barycenter — always equal to `binary_primary_position_* + binary_mutual_position_*`, but stored explicitly rather than derived on read (same convention `binary_mutual_position_*_km` itself already set). |
| `binary_secondary_mass_fraction` | DOUBLE | nullable | Added in v20. `secondary_mass / (primary_mass + secondary_mass)`, constant since masses don't change — stored so `_db.advance_orbital_phases` never needs to join back to `stars` for either mass. |
| `binary_planetary_wobble_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v20. NULL unless `binary_configuration = 'close'`. Additional pair-wide wobble from circumbinary planets (`planets.star_id IS NULL`) — modeled as one shared wobble applied to the whole pair rather than split unevenly between primary/secondary, which would need a real 3+-body solve. |
| `system_flavor_text` | TEXT | nullable | Decided once at generation time (Phase 0 fix). |
| `runaway_class` | VARCHAR(16) | nullable, `runaway` or `hypervelocity` | Added in v37. NULL for an ordinary star; set by `generate.flag_fast_stars`. |
| `runaway_speed_kms` | DOUBLE | nullable | Added in v37. The star's speed relative to its neighbors when `runaway_class` is set. |
| `velocity_x_kms`, `velocity_y_kms`, `velocity_z_kms` | DOUBLE | NOT NULL, default 0 | Added in v61 (GEN.121). The system's velocity, km/s on the galactic axes: the rotation curve's tangent at its place plus its runaway motion. `advance_galactic_positions` turns it with the position (a hypervelocity star instead flies a straight line, `p0 + v (t − t0)`, and keeps its velocity; GEN.137, `galaxy/straight_line.py`). 0 while the system has no place in the galaxy. |
| `epoch_unix`, `next_update_due` | DOUBLE | nullable, `next_update_due` indexed | Added in v66 (GEN.106). The galactic orbit's clock (see `orbit_simulation_state`): the system, its stars' `galactic_orbital_phase_deg` and a close pair's `binary_galactic_orbital_phase_deg` move together when it is due. `binary_epoch_unix`/`binary_next_update_due` are the same for a pair's mutual orbit. |
| `schema_version` | INTEGER | NOT NULL, default 1 | See "Versioning" above. |
| `mediawiki_url` | TEXT | nullable | Where this system's page lives (or should live) on MediaWiki. |
| `wikijs_url` | TEXT | nullable | Where this system's page lives (or should live) on Wiki.js. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | FK -> `nebulae.id` / `supernova_remnants.id`, `ON DELETE SET NULL`, nullable | Added in v39. The innermost nebula or supernova remnant whose sphere holds this system; at most one is set. See "Containment (v39)" above. |
| `created_at` | TIMESTAMP | NOT NULL, default `CURRENT_TIMESTAMP` | |
| `modified_at` | TIMESTAMP(3) | NOT NULL, `ON UPDATE CURRENT_TIMESTAMP(3)` | Added in v27. Also bumped when a child row (star, planet, moon, belt, comet) changes. |

**Quadrant labeling.** `quadrant` reuses the generator's own octant scheme
(`src/planetgen/galaxy/sector.py`'s `classify_octant`, backed by
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
| `name` | TEXT | NOT NULL | A single star shares the system's name (renaming one renames the other). A wide binary's stars are the system name's first word plus their own word, e.g. `"Voranthis Kelmoor"` and `"Voranthis Pikkita"`; a close pair's are `"Voranthis A"` and `"Voranthis B"` (GEN.62; names stored before keep theirs). |
| `star_type` | TEXT | NOT NULL | Full descriptive string, e.g. `"G2V Yellow Main Sequence Star"` — unrelated to `planets.body_type`'s single-character code. |
| `yerkes_class` | TEXT | NOT NULL | e.g. `"V"`, `"VII"` (white dwarf). |
| `mass_kg`, `radius_km`, `luminosity_w` | DOUBLE | NOT NULL | |
| `temperature_k` | DOUBLE | NOT NULL | |
| `age_gy` | DOUBLE | NOT NULL | |
| `lifespan_gy` | DOUBLE | nullable | **NULL = `float('inf')`** — white dwarfs (`yerkes_class` `'VII'`/`'D'`). Never the JSON `Infinity` token. |
| `habitable_zone_inner_km`, `_outer_km` | DOUBLE | NOT NULL | |
| `system_perimeter_km` | DOUBLE | NOT NULL | Hill sphere relative to the galaxy. |
| `heliosphere_radius_km` | DOUBLE | NOT NULL | |
| ~~`galactic_orbital_speed_kms`~~ | (dropped) | | Added in v10. Circular orbital speed around the galactic center (`galactic_orbit.calculate_galactic_orbit`), from this system's actual distance from the galactic center where known (a sector-placed system), or the fixed `physical_constants.GALACTIC_CENTER_DISTANCE_LY` fallback otherwise — same fallback convention as `system_perimeter_km`. Its speed column was dropped in v103 (GEN.336): the speed is the length of `velocity_x/y/z_kms` (a star and a close pair read their system's). |
| `galactic_orbital_period_gy` | DOUBLE | NOT NULL | Added in v10. Orbital period for the circular orbit above, in billions of years (Gy) — the same unit `age_gy`/`lifespan_gy` use. |
| `galactic_orbital_phase_deg` | DOUBLE | NOT NULL | Added in v13. This star's current angular position around its galactic orbit — the same role `planets.orbital_phase_deg` plays, advanced by `_db.advance_orbital_phases`. Both stars of a binary pair always carry the identical value (see "Schema history" above for why — `StarSystem.__init__` rolls it once and threads it to primary/secondary/proxy alike). |
| `galactic_min_update_interval_years` | DOUBLE | NOT NULL | Added in v13. Floating-point update guard for `galactic_orbital_phase_deg`, same formula as `planets.min_update_interval_years` (`orbits.minimum_update_interval_years`), applied to `galactic_orbital_period_gy * 1e9` years. |
| `wide_binary_a_crit_km` | DOUBLE | nullable | Added in v15. This star's own Holman & Wiegert (1999) critical semi-major axis (`orbits.holman_wiegert_critical_semimajor_axis`) — the maximum orbit distance that stays long-term stable given its companion's perturbation. NULL for a single star or either constituent of a `'close'` pair; populated for both stars of a `'wide'` pair. |
| `reflex_offset_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v20. This star's own displacement from its nominal fixed point, from the combined pull of every planet orbiting it directly (`planets.star_id`) — see `orbits.calculate_reflex_offset`. NULL/0 with no planets. Class-blind: applies identically to an anchored `black_holes`/`neutron_stars` row. |
| `rotation_period_hours` | DOUBLE | nullable | Added in v68 (GEN.104). Rotation period: gyrochronology for cool dwarfs, a log-normal speed under breakup for hot stars, slow for giants, about a day for white dwarfs; a black hole's horizon period, a neutron star's spin period. The spin vector is the axis times 2 pi over this period. |
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the galactic pole; see `planetgen/physics/spin.py` and docs/design/orbital-updates.md section 6. NULL on a row generated before v68. |
| `log_lx_lbol`, `l_xuv_w`, `xuv_saturated`, `flare_n33_per_yr`, `flare_alpha`, `xuv_fluence_j` | DOUBLE (`xuv_saturated` BOOLEAN) | nullable | Added in v72 (GEN.86). Stellar activity: log10 of the coronal X-ray share of the luminosity (saturated at 7.4e-4 for G, K and M stars, decaying after a mass-dependent saturation age; NULL for a white dwarf or compact remnant), the X-ray plus EUV luminosity in W (a hot photosphere's output above 13.6 eV included; a neutron star's thermal or spin-down X-rays; a black hole's disk), whether the corona is still saturated, flares above 1e33 erg a year and their energy slope, and the XUV given off over the star's age, J. See `planetgen/physics/activity.py` and docs/design/activity-magnetism-radiation-hydrosphere.md section 2. NULL on a row generated before v72. |
| `lethal_event_rate_per_gyr` | DOUBLE | nullable | Added in v75 (GEN.87). Supernovae within 8 pc per Gyr here: 1.5 at the Sun's 8 kpc, rising as exp((8 kpc - R) / 1.45 kpc) to at most 250 times that at the centre. NULL for a star with no galaxy placement and on a row generated before v75. |
| `metallicity_fe_h` | DOUBLE | nullable | Added in v95 (GEN.158). The star's [Fe/H] in dex against the Sun, drawn from its population and its distance from the galactic centre; a binary's secondary shares its primary's. NULL on a row generated before v95. |

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
| `orbital_phase_deg` | DOUBLE | NOT NULL | Added in v9. This body's current position angle around its orbit — the one orbital-motion column that changes over time, advanced in place by `planetgen.cli.orbits` (see `orbit_simulation_state` below). |
| `position_x_km`, `_y_km`, `_z_km` | DOUBLE | NOT NULL | Added in v11. This body's Cartesian position relative to its orbital anchor — the star (or a binary's combined center) for a planet — derived from `distance_km` and the three orbital-motion columns above (`orbits.orbital_position_au`). Changes in lockstep with `orbital_phase_deg` as `planetgen.cli.orbits` advances it. |
| `velocity_x_kms`, `_y_kms`, `_z_kms` | DOUBLE | NOT NULL, default 0 | Added in v60 (GEN.121). This body's velocity, km/s, relative to its orbital anchor and on the same axes as `position_x/y/z_km`: the tangent of its circular orbit at `orbital_phase_deg`, `orbital_speed_kms` long (`orbits.circular_orbital_velocity_au_per_year`). `advance_orbital_phases` moves it with the position. |
| `epoch_unix`, `next_update_due` | DOUBLE | nullable, `next_update_due` indexed | Added in v66 (GEN.106). When the stored position holds and when the body will have moved 0.01 AU (100,000 km for a moon) at its orbital speed; see `orbit_simulation_state`. An edit that changes the orbit sets `next_update_due` back to NULL. |
| `orbital_speed_kms` | DOUBLE | NOT NULL | Added in v11. Constant circular-orbit speed (`orbits.circular_orbital_speed_kms`, `v = 2*pi*r/T`). Only changes if `distance_km`/`period_years` do (e.g. `StarSystem.validate_system` resolving an orbital overlap at generation time), never from phase advancing alone. |
| `min_update_interval_years` | DOUBLE | NOT NULL | Added in v12; not read by the update since v66 (`next_update_due` decides). Not a narrative stat -- a floating-point update guard: the shortest `elapsed_years` worth calling it for, below which the phase delta added is smaller than `orbital_phase_deg`'s own double-precision resolution and so is guaranteed to be a no-op write (`orbits.minimum_update_interval_years`, `period_years * math.ulp(360.0) / 360`). Like `orbital_speed_kms`, only changes if `distance_km`/`period_years` do. |
| `rotation_period_hours` | DOUBLE | NOT NULL | Added in v9. Axial rotation ("day length") — a static descriptive stat; no rotational phase is tracked. |
| `reflex_offset_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v20. This planet's own displacement from its nominal fixed point, from the combined pull of its own moons (`moons.planet_id`) — see `orbits.calculate_reflex_offset`. NULL/0 with no moons. **Not present on `moons`** — a moon never hosts its own moons. |
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the orbit normal (the orbit frame: inclination and ascending node); see `planetgen/physics/spin.py` and docs/design/orbital-updates.md section 6. NULL on a row generated before v68. A tidally locked body (now planets as well as moons) has tilt 0 and `rotation_period_hours` equal to its orbit. |
| `mantle_redox` | VARCHAR(16) | nullable, CHECK reduced / intermediate / oxidized | Added in v70 (GEN.85). The upper mantle's redox state: reduced at or below the iron-wustite buffer, oxidized from IW+2.5. NULL for a gas giant. See `planetgen/physics/atmosphere.py`. |
| `mantle_delta_iw` | DOUBLE | nullable | Added in v70 (GEN.85). The upper mantle's oxygen fugacity, log units above IW, drawn from mass (Earth about +3.5, the Moon about -1). |
| `p_o2_kpa`, `p_co2_kpa`, `p_co_kpa`, `p_n2_kpa`, `p_ar_kpa`, `p_h2_kpa`, `p_h2o_kpa`, `p_ch4_kpa`, `p_h2s_kpa`, `p_so2_kpa` | DOUBLE | nullable | Added in v70 (GEN.85). Each gas's partial pressure, kPa, from the class's mix, the mantle redox and the cold side (a gas past its vapour pressure is ice). 0 for an airless body. `atmosphere` is written from them; a gas outside the ten (helium, ammonia, sodium) is named there but not stored, so on a giant they sum to less than `atmospheric_pressure_pa`. |
| `magnetic_moment_a_m2`, `dipole_class`, `magnetopause_rp` | DOUBLE, VARCHAR(16) (CHECK none / weak / earth-like / strong / multipolar), DOUBLE | nullable | Added in v72 (GEN.86). The magnetic dipole moment, A m2 (0 with no dynamo), its class (multipolar when slow rotation leaves under half of it in the dipole), and the magnetopause standoff in body radii against the star's wind (NULL with no field). See `planetgen/physics/magnetism.py` and docs/design/activity-magnetism-radiation-hydrosphere.md section 3. NULL on a row generated before v72. |
| `xuv_flux_earth`, `xuv_exposure_index`, `flare_irradiation_index` | DOUBLE | nullable | Added in v72 (GEN.86). The star's XUV flux here in units of Earth's today, its lifetime XUV here against a 1 solar mass star's at 1 AU over 5 Gyr (2 or less Earth-like, over 30 extreme), and its flares above 1e33 erg a year over the distance squared in AU (Earth 0.003). A moon takes its planet's distance. NULL on a row generated before v72. |
| `water_mass_fraction`, `hydrosphere`, `ocean_fraction`, `land_fraction` | DOUBLE, VARCHAR(20) (CHECK dry / vapour / ice / ice-covered ocean / surface ocean / hycean), DOUBLE, DOUBLE | nullable | Added in v74 (GEN.88). The body's water share of its mass (drawn by class), where its surface water is (gone, in the air, frozen, under ice, open, or under a hydrogen sky), and the shares of the surface the water (liquid or frozen) and the land cover. NULL for a gas giant. See `planetgen/physics/hydrosphere.py` and docs/design/activity-magnetism-radiation-hydrosphere.md section 5. NULL on a row generated before v74. |
| `ocean_depth_km`, `ice_shell_km`, `hp_ice_km` | DOUBLE | nullable | Added in v74 (GEN.88). The mean liquid depth where there is water (NULL without liquid), the ice lid, and the high-pressure ice (VI or VII) under the deepest liquid the liquidus allows (NULL when there is none). |
| `ocean_class`, `ocean_ph`, `water_activity`, `phosphorus` | VARCHAR(16) (CHECK ice-sealed / chloride brine / acid sulfate / soda / neutral), DOUBLE, DOUBLE, VARCHAR(8) (CHECK high / limited / starved) | nullable | Added in v74 (GEN.88). The ocean's chemistry class by rule (section 5.4), its pH and water activity drawn in the class's range, and its phosphorus supply. NULL without liquid water. |
| `surface_dose_msv_yr`, `dose_gcr_msv_yr`, `dose_sep_msv_yr`, `dose_ground_msv_yr`, `dose_helio_mult` | DOUBLE | nullable | Added in v75 (GEN.87). The yearly dose at the surface, mSv, and its parts: `surface_dose_msv_yr` = `dose_gcr_msv_yr` x `dose_helio_mult` + `dose_sep_msv_yr` + `dose_ground_msv_yr`. Cosmic rays through the air (a log-log curve through the Moon, Mars and a 100 g/cm2 point, then 164 g/cm2 attenuation) cut to 0.6 at most by a dipole; stellar particles (20 mSv/yr x sqrt(N33 / N33 of the Sun) / d_AU^2, stopped by the air and steered by the dipole); crust gamma and radon (0.48 mSv/yr x abundance x radiogenic heat at the star's age against 4.5 Gyr, plus 0.5 mSv/yr x abundance with air and land). `dose_helio_mult` is 1 to 2.5: 1 in open space, raised while a nebula or remnant presses on the star's heliosphere (`refresh_containment` and an edit update it and the total). NULL for a gas giant. |
| `uv_surface_index`, `ozone_loss_flag` | DOUBLE, BOOLEAN | nullable | Added in v75 (GEN.87). DNA-weighted surface UV against Earth's 1 (the star's 200 to 300 nm power against the Sun's, floored at 0.1, times its flux, times the ozone shielding (O3 / O3 of Earth)^-1.6 up to 1,000 when there is none), and whether ozone is lost: the layer exists but the star's `lethal_event_rate_per_gyr` is over 10 (one per 100 Myr) or particles give 100 mSv/yr at the ozone layer. NULL for a gas giant, and `uv_surface_index` NULL without a star temperature. |
| `phi4`, `phi4_pressure`, `phi4_temperature`, `phi4_air`, `phi4_radiation` | DOUBLE | nullable | Added in v76 (GEN.89); `phi4_chemistry` was renamed `phi4_air` in v87. The PHI-4 display score for a human visitor: the four domains each 0 to 1 (1 is Blue) and `phi4`, their geometric mean (0 when any domain is 0). NULL for a gas giant and on a row generated before v76. |
| `tier_pressure`, `tier_temperature`, `tier_air`, `tier_radiation` | TINYINT | nullable | Added in v76 (GEN.89); `tier_chemistry` was renamed `tier_air` in v87. Each domain's colour: 0 Blue, 1 Green, 2 Yellow, 3 Red. |
| `equipment_tier` | TINYINT | nullable | Added in v76 (GEN.89). The kit a human needs: 0 ideal (named shirtsleeve before 2026-10-10), 1 breathing mask, 2 mask with scrubber, 3 environment suit, 4 full life support with radiation hardening. A pulsar's, neutron star's or black hole's planet is always 4. |
| `phi_bio`, `phi_cpx`, `phi_tech`, `l_solv`, `l_chem`, `l_ener`, `l_rad` | DOUBLE | nullable | Added in v76 (GEN.89). The Xenobiology numbers: microbial life (the product of the four likelihoods `l_solv` solvent, `l_chem` chemistry, `l_ener` energy and `l_rad` radiation), complex life (`phi_bio` times the metazoa's gas tolerance) and human operability with equipment (`phi_tech`, the geometric mean of the pressure, temperature, water and radiation mitigation factors). |
| `hab_note` | VARCHAR(255) | nullable | Added in v76 (GEN.89). The domains short of Blue ("Pressure Red, Radiation Yellow"), "All four domains Blue", or "Planet of a pulsar, neutron star or black hole: lethal radiation, rated as 1,000 Sv/yr at the surface". |
| `energy_flux_w_m2` | DOUBLE | nullable | Added in v76 (GEN.89). The light or chemical power life can use: 6% of the starlight, or 1% of the body's own heat flow, whichever is more (the heat only, under an ice lid). The one input the score keeps, so the dose refresh can re-score from the row. |
| `flux_variation` | DOUBLE | nullable | Added in v87 (DB.26, GEN.148). The relative swing of the stellar flux over the orbit or pair cycle; 0 is steady. NULL until regenerated. |
| `induced_moment_a_m2` | DOUBLE | nullable | Added in v87 (DB.26, GEN.178), placed after `dipole_class`. The induced magnetic moment in A m2; NULL when not computed or none. |
| `phi4_water`, `phi4_suitability` | DOUBLE | nullable | Added in v87 (DB.26, GEN.207, PHI-6), beside `phi4_air`. The water score and the suitability (the lowest of water, food and shelter). NULL until regenerated. |
| `tier_water`, `tier_suitability` | TINYINT | nullable | Added in v87 (DB.26, GEN.207). The colour of each, 0 Blue to 3 Red like the other tiers. |

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
| `body_type`, `name`, `planet_class`, `distance_km` (from the parent planet), `radius_km`, `mass_kg`, `volume_km3`, `period_years`, `zone`, `description`, `gravity_g`, `surface_temperature_k`, `density_g_cm3`, `atmosphere`, `atm_density`, `atm_molar_density`, `atmospheric_pressure_pa`, `composition`, `scale_height_km`, `hill_radius_km`, `min_orbit_distance_km`, `habitable_zone_inner_km`, `_outer_km`, `life_chemical`, `evolutionary_speed`, `flavor_text`, `flavor_text_count`, `orbital_inclination_deg`, `orbital_ascending_node_deg`, `orbital_phase_deg`, `position_x_km`, `_y_km`, `_z_km`, `orbital_speed_kms`, `velocity_x_kms`, `_y_kms`, `_z_kms`, `epoch_unix`, `next_update_due`, `min_update_interval_years`, `rotation_period_hours`, `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg`, `mantle_redox`, `mantle_delta_iw`, `p_<gas>_kpa`, `magnetic_moment_a_m2`, `dipole_class`, `magnetopause_rp`, `xuv_flux_earth`, `xuv_exposure_index`, `flare_irradiation_index`, `water_mass_fraction`, `hydrosphere`, `ocean_fraction`, `land_fraction`, `ocean_depth_km`, `ice_shell_km`, `hp_ice_km`, `ocean_class`, `ocean_ph`, `water_activity`, `phosphorus`, `surface_dose_msv_yr`, `dose_gcr_msv_yr`, `dose_sep_msv_yr`, `dose_ground_msv_yr`, `dose_helio_mult`, `uv_surface_index`, `ozone_loss_flag`, `phi4`, `phi4_<domain>`, `tier_<domain>`, `phi_bio`, `phi_cpx`, `phi_tech`, `l_solv`, `l_chem`, `l_ener`, `l_rad`, `equipment_tier`, `hab_note`, `energy_flux_w_m2` | — | — | Identical meaning/type/nullability to the same-named column on `planets` above, except `position_x/y/z_km` are relative to *this moon's* orbital anchor — its parent planet, not the star. |

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
| `velocity_x_kms`, `_y_kms`, `_z_kms` | DOUBLE | NOT NULL, default 0 | Added in v60 (GEN.121). The velocity, km/s, relative to the star on the position's axes (`kepler.comet_orbital_state`), recomputed with the position by `advance_comet_orbits`. |
| `rotation_period_hours` | DOUBLE | nullable | Added in v68 (GEN.104). A log-normal small-body period, never under the 2.2-hour spin barrier. The spin vector is the axis times 2 pi over this period. |
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the orbit normal (YORP: near 10 or 170 degrees); see `planetgen/physics/spin.py` and docs/design/orbital-updates.md section 6. NULL on a row generated before v68. |

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

Added in v16, for `planetgen phenomenon`'s separate, rarer exotic-phenomenon
generation mode (`planetgen/generation/phenomena/compact_remnant.py`). Satellite tables
extending a `stars` row (the same "extra detail alongside an existing row"
shape `asteroid_belt_composition` has to `asteroid_belts`) when a compact
remnant anchors a full `StarSystem` (`planetgen phenomenon --anchor-system`)
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
| `galactic_orbital_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | nullable | Added in v17. NULL when `star_id` is set (an anchored remnant's motion lives on its own `stars` row instead); populated only for a standalone black hole. See `stars`' identical columns below. Its speed column was dropped in v103 (GEN.336): the speed is the length of `velocity_x/y/z_kms` (a star and a close pair read their system's). |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE SET NULL`, nullable | Added in v21. The sector this standalone black hole was generated as part of (`generate.generate_sector_phenomena`) or placed near (`planetgen phenomenon --sector-id`), the same convention `nebulae.sector_id` uses. Always NULL when `star_id` is set. |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v21. This black hole's own galaxy-frame center, the same shape `nebulae`'s identical columns use (see below) -- but here it's usually the *exact* position `SpaceSector.add_phenomenon`'s Hill-sphere-aware in-sector placement computed at generation time (never within a neighboring star system's or another compact remnant's own Hill sphere), converted to galaxy-frame coordinates, rather than `compute_phenomenon_placement`'s independent random jitter (still used for `planetgen phenomenon`'s own standalone `--sector-id`, which has no specific in-sector position to convert). |
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
| `pulsar_type` | TEXT | NOT NULL, CHECK IN ('young','millisecond','magnetar','non-pulsing') | `magnetar` since v96 (GEN.153). |
| `surface_temperature_k`, `luminosity_w` | DOUBLE | NOT NULL | `luminosity_w` is thermal blackbody emission (Stefan-Boltzmann), always nonzero. |
| `age_gy` | DOUBLE | NOT NULL | Time since the core-collapse supernova. |
| `galactic_orbital_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | nullable | Added in v17. Same nullable-when-anchored convention as `black_holes` above. Its speed column was dropped in v103 (GEN.336): the speed is the length of `velocity_x/y/z_kms` (a star and a close pair read their system's). |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE SET NULL`, nullable | Added in v21. Same convention as `black_holes.sector_id` above. |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v21. Same convention as `black_holes`' identical columns above. |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | nullable, `ON DELETE SET NULL` | Added in v39. The innermost cloud holding it, as on `star_systems`. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the galactic pole; see `planetgen/physics/spin.py` and docs/design/orbital-updates.md section 6. NULL on a row generated before v68. NULL for an anchored remnant: its spin is on its `stars` row. |

### `nebulae`

Added in v16. Always standalone (nothing in this generator places a
nebula within a `StarSystem`); `sector_id` was reserved for a future
sector-context encounter through v17, unused (always NULL) by
`planetgen phenomenon` until v18 gave it one.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE SET NULL`, nullable | v18: the nearest already-generated sector to `center_x/y/z_pc` below -- a convenience "home" link, not this nebula's real geometry (its sphere may overlap several sectors, or none). NULL iff `center_x/y/z_pc` are NULL. `ON DELETE SET NULL` (not `CASCADE`, unlike v16/v17): deleting that sector doesn't delete a nebula that merely happens to be near it. |
| `name` | TEXT | NOT NULL | |
| `nebula_class` | CHAR(1) | NOT NULL | Added in v38: letter class A-Q (`program_constants.NEBULA_CLASSES`). |
| `nebula_type` | TEXT | NOT NULL, CHECK IN ('diffuse','emission','reflection','planetary','dark') | The class's family (`diffuse` added in v38). |
| `dominant_species`, `density_cm3`, `temperature_k`, `extinction_av` | VARCHAR(255) / DOUBLE | NOT NULL | Added in v38: what the cloud holds, its particle density nH (cm⁻³), gas temperature (K) and optical extinction (magnitudes). |
| `age_myr`, `q_total`, `hosts` | DOUBLE, DOUBLE, VARCHAR(512) | nullable | Added in v87 (DB.26, GEN.150), after `extinction_av`. The nebula's age (Myr), the total ionizing photon rate (photons/s) and its host stars as counts by spectral type ("O5V x1, O9V x4"). NULL until regenerated and for classes with no hosts. |
| `radius_ly` | DOUBLE | NOT NULL | |
| `composition`, `formation_cause` | TEXT | NOT NULL | Descriptive strings, one per `nebula_type` (`program_constants.NEBULA_TYPES`). |
| `galactic_orbital_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | NOT NULL | Added in v17. Always populated (a nebula is always standalone). See `stars`' identical columns above. Its speed column was dropped in v103 (GEN.336): the speed is the length of `velocity_x/y/z_kms` (a star and a close pair read their system's). |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v18 (`planetgen phenomenon --sector-id`): this nebula's own galaxy-frame center, in the same Cartesian space `sectors.center_x/y/z_pc` uses -- a sphere (this + `radius_ly`), not a sector-relative offset, since a nebula (up to 200 ly across) is frequently far larger than one sector (default edge 4 pc, ~13 ly) and may overlap several. NULL together: never placed in the galaxy (still the default -- `--sector-id` is optional). See `queryDb.phenomena_near_sector`/`galaxy_placed_phenomena` for how this is read back, and `docs/design/galaxy-coordinate-system.md` for the coordinate system itself. |
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
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE CASCADE`, nullable | The sector it was generated as part of (`planetgen sector`) or placed near (`planetgen phenomenon --sector-id`); NULL for one generated standalone. |
| `name` | TEXT | NOT NULL | |
| `remnant_class` | CHAR(1) | NOT NULL | Added in v38: letter class R-W (`program_constants.NEBULA_CLASSES`). |
| `morphology` | TEXT | NOT NULL, CHECK IN ('shell','plerion','composite') | Fixed by the class since v38. |
| `dominant_species`, `density_cm3`, `temperature_k`, `extinction_av` | VARCHAR(255) / DOUBLE | NOT NULL | Added in v38, as on `nebulae`. |
| `age_years` | DOUBLE | NOT NULL | |
| `radius_ly` | DOUBLE | NOT NULL | Derived from `age_years` and `density_cm3` via the Sedov-Taylor blast-wave relation (`radius ∝ (1e51 erg / density)^(1/5) age^(2/5)`), held at the snowplow onset radius. |
| `progenitor_type` | TEXT | NOT NULL, CHECK IN ('Type Ia','core-collapse') | |
| `compact_remnant_kind` | TEXT | nullable, CHECK IN ('black_hole','neutron_star') | NULL when no compact remnant is embedded. |
| `compact_remnant_black_hole_id` | INTEGER | FK -> `black_holes.id`, `ON DELETE SET NULL`, nullable | |
| `compact_remnant_neutron_star_id` | INTEGER | FK -> `neutron_stars.id`, `ON DELETE SET NULL`, nullable | |
| `galactic_orbital_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | NOT NULL | Added in v17. Always populated (a supernova remnant is always standalone). Its speed column was dropped in v103 (GEN.336): the speed is the length of `velocity_x/y/z_kms` (a star and a close pair read their system's). |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v28: this remnant's own galaxy-frame center, the same shape `nebulae` uses. NULL together: never placed in the galaxy. |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |

### `rogue_planets`

Added in v16. Always standalone — a rogue planet is by definition unbound
from any star.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE CASCADE`, nullable | The sector it was generated as part of (`planetgen sector`) or placed near (`planetgen phenomenon --sector-id`); NULL for one generated standalone. |
| `name` | TEXT | NOT NULL | |
| `planet_type` | TEXT | NOT NULL, CHECK IN ('t','g') | Same letters as `planets.body_type`. |
| `planet_class` | VARCHAR(4) | nullable | Added in v47: its `PLANET_CLASSES` letter (GEN.8); NULL for a brown dwarf. |
| `mass_bin` | VARCHAR(16) | NOT NULL | Added in v37: `terrestrial`, `sub-neptune`, `saturn`, `jupiter` or `brown-dwarf`. |
| `mass_kg`, `radius_km` | DOUBLE | NOT NULL | |
| `composition` | TEXT | NOT NULL | Descriptive bulk-composition string. |
| `has_internal_heat`, `has_moons` | BOOLEAN | NOT NULL | `has_internal_heat`: still geologically active (computed since v48; a giant always is). |
| `age_gy`, `internal_heat_flux_w_m2`, `effective_temperature_k` | DOUBLE | nullable | Added in v48: its age and its own heat (W/m^2, K). |
| `surface_regime` | VARCHAR(24) | nullable | Added in v48: `bare-rock`, `frozen-atmosphere`, `ice-shell-ocean`, `ice-world`, `hydrogen-envelope`, `hycean` (v73: an ocean under the envelope with no high-pressure ice beneath), `gas-giant` or `brown-dwarf`. |
| `surface_temperature_k`, `surface_pressure_pa` | DOUBLE | nullable | Added in v48: at the ground, or a giant's at 1 bar. |
| `ice_shell_thickness_km`, `ocean_depth_km`, `has_liquid_water` | DOUBLE, DOUBLE, BOOLEAN | nullable | Added in v48: its ice lid and buried (or envelope-warmed) ocean; NULL when there is none. |
| `hp_ice_km` | DOUBLE | nullable | Added in v74 (GEN.88): high-pressure ice under the deepest liquid the liquidus allows; NULL when there is none. |
| `galactic_orbital_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | NOT NULL | Added in v17. Always populated (a rogue planet is always standalone). Its speed column was dropped in v103 (GEN.336): the speed is the length of `velocity_x/y/z_kms` (a star and a close pair read their system's). |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v28: this rogue planet's own galaxy-frame center, the same shape `nebulae` uses. NULL together: never placed in the galaxy. |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | nullable, `ON DELETE SET NULL` | Added in v39. The innermost cloud holding it, as on `star_systems`. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |
| `rotation_period_hours` | DOUBLE | nullable | Added in v68 (GEN.104). A planet's day for its type. The spin vector is the axis times 2 pi over this period. |
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the galactic pole; see `planetgen/physics/spin.py` and docs/design/orbital-updates.md section 6. NULL on a row generated before v68. |
| `log_lx_lbol`, `l_xuv_w`, `xuv_saturated`, `flare_n33_per_yr`, `flare_alpha`, `xuv_fluence_j` | DOUBLE (`xuv_saturated` BOOLEAN) | nullable | Added in v72 (GEN.86). Stellar activity: log10 of the coronal X-ray share of the luminosity (saturated at 7.4e-4 for G, K and M stars, decaying after a mass-dependent saturation age; NULL for a white dwarf or compact remnant), the X-ray plus EUV luminosity in W (a hot photosphere's output above 13.6 eV included; a neutron star's thermal or spin-down X-rays; a black hole's disk), whether the corona is still saturated, flares above 1e33 erg a year and their energy slope, and the XUV given off over the star's age, J. See `planetgen/physics/activity.py` and docs/design/activity-magnetism-radiation-hydrosphere.md section 2. NULL on a row generated before v72. |

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
| `sector_id` | INTEGER | FK -> `sectors.id`, `ON DELETE CASCADE`, nullable | The sector it was generated as part of (`planetgen sector`) or placed near (`planetgen phenomenon --sector-id`); NULL for one generated standalone. |
| `name` | TEXT | NOT NULL | |
| `nucleus_diameter_km`, `velocity_kms` | DOUBLE | NOT NULL | `velocity_kms` is hyperbolic excess speed relative to any star it passes — a separate, non-advancing descriptive stat from the `galactic_orbital_*` columns below (its own bulk galactic motion). |
| `is_active` | BOOLEAN | NOT NULL | Whether it currently shows a coma/tail. |
| `composition_summary` | TEXT | NOT NULL | Human-readable summary, same role as `asteroid_belts.composition_summary`. |
| `galactic_orbital_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | NOT NULL | Added in v17. Always populated (an interstellar comet is always standalone). Its speed column was dropped in v103 (GEN.336): the speed is the length of `velocity_x/y/z_kms` (a star and a close pair read their system's). |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v28: this comet's own galaxy-frame center, the same shape `nebulae` uses. NULL together: never placed in the galaxy. |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | nullable, `ON DELETE SET NULL` | Added in v39. The innermost cloud holding it, as on `star_systems`. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |
| `rotation_period_hours` | DOUBLE | nullable | Added in v68 (GEN.104). A log-normal small-body period, never under the 2.2-hour spin barrier. The spin vector is the axis times 2 pi over this period. |
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the galactic pole (YORP: near 10 or 170 degrees); see `planetgen/physics/spin.py` and docs/design/orbital-updates.md section 6. NULL on a row generated before v68. |

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
encounter through v17, unused (always NULL) by `planetgen phenomenon` until
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
| `galactic_orbital_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | NOT NULL | Always populated (an asteroid field is always standalone). Its speed column was dropped in v103 (GEN.336): the speed is the length of `velocity_x/y/z_kms` (a star and a close pair read their system's). |
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

Name-uniqueness bookkeeping (v24, `planetgen/names/uniqueness.py`).
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
also updates its neighbors' lists) and by `planetgen.cli.orbits`. See
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
| `velocity_x_kms`, `velocity_y_kms`, `velocity_z_kms` | DOUBLE | NOT NULL, default 0 | Added in v64 (GEN.125). A stand-alone facility's velocity, km/s on the galactic axes: the rotation curve's tangent at its place, set when it is added and turned with its position by `advance_galactic_positions`. 0 for a facility on a body, in orbit or in a belt, whose motion is the `orbit_*` columns below, relative to its host. |
| `orbit_distance_km`, `orbit_period_years`, `orbital_speed_kms`, `orbit_phase_deg` | DOUBLE | nullable | An orbital facility's circular orbit, from its host's mass; an asteroid facility in a belt gets one too, around its star from a random spot in the belt (no schema change). `planetgen.cli.orbits` advances both. |
| `description` | TEXT | nullable | |
| `type_id` | BIGINT UNSIGNED | FK -> `facility_types.id`, `ON DELETE RESTRICT`, nullable, indexed | Added in v89 (DB.28). The type that refines `kind`. |
| `affiliation_polity_id`, `affiliation_name` | BIGINT UNSIGNED, VARCHAR(255) | FK -> `polities.id`, `ON DELETE SET NULL`, nullable | Added in v89 (DB.28). The polity it belongs to, or a free-text affiliation. |
| `rating_crime`, `rating_housing`, `rating_resources`, `rating_maintenance`, `rating_health` | TINYINT UNSIGNED | nullable, CHECK BETWEEN 1 AND 3, each indexed | Added in v89 (DB.28). The Green (1), Yellow (2) or Red (3) rating; NULL until rated. |
| `ratings_source`, `ratings_note`, `ratings_set_at` | VARCHAR(16), TEXT, TIMESTAMP | nullable | Added in v89 (DB.28). How the ratings were set, a note, and when. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | |

### `facility_types`, `facility_type_fields`, `facility_field_values`

Added in v89 (DB.28, POP.10); design in `docs/design/population-and-politics.md`, "Data model" of the facility section. A type refines a kind.

| Table | Columns |
|---|---|
| `facility_types` | `id`, `name` (UNIQUE), `kind` (CHECK like `facilities.kind`), `icon`, `description`, `sort_order`, `active` (0 or 1) |
| `facility_type_fields` | `id`, `type_id` (FK `ON DELETE CASCADE`), `key`, `label`, `datatype` (CHECK IN 'int','decimal','text','bool','choice'), `unit`, `required`, `min_value`, `max_value`, `choices` (newline-separated), `default_text`, `sort_order`; UNIQUE (`type_id`, `key`) |
| `facility_field_values` | `facility_id` (FK `ON DELETE CASCADE`), `field_id` (FK `ON DELETE RESTRICT`: retire a field, do not delete it), `value_int`, `value_dec`, `value_text`, `value_bool`; PRIMARY KEY (`facility_id`, `field_id`) |

The key pattern `^[a-z][a-z0-9_]{0,31}$`, at most 24 fields per type and the value checks are made by the application, not a CHECK.

### `bright_stars`

Added in v43. Every star at least `galaxy_shape.bright_star_min_luminosity_sol`
bright, generated and placed galaxy-wide by `planetgen plan`'s scatter
before any sector is filled. See "Bright-star pre-placement (v43)" above.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `shell_index`, `band_index`, `wedge_index` | INT / INT / INT | NOT NULL | The sector cell it falls in. Indexed together. |
| `position_x_mpc`, `position_y_mpc`, `position_z_mpc` | BIGINT | NOT NULL | Galaxy-frame position, milliparsecs. |
| `population` | VARCHAR(12) | NOT NULL, CHECK IN ('young','intermediate','old','bulge','halo','cluster') | |
| `star_type`, `yerkes_class` | VARCHAR | NOT NULL | As on `stars`. |
| `mass_kg`, `radius_km`, `temperature_k`, `luminosity_w`, `age_gy` | DOUBLE | NOT NULL | The finished star. `luminosity_w` is indexed. |
| `lifespan_gy`, `phase_end_age_gy` | DOUBLE | nullable | |
| `initial_mass_sol` | DOUBLE | NOT NULL | A companion's mass is drawn from it. |
| `seed` | BIGINT UNSIGNED | NOT NULL | Seed for the system built around it later. |
| `star_system_id` | BIGINT UNSIGNED | FK -> `star_systems.id`, `ON DELETE SET NULL`, nullable | Set when its sector is filled. |
| `created_at` | TIMESTAMP | NOT NULL | |
| `metallicity_fe_h` | DOUBLE | nullable | Added in v95 (GEN.158). The star's [Fe/H] in dex, drawn at scatter time from its population and position on a stream of its own seed; the fill that builds its system keeps it. NULL on older rows and on a cluster star until GEN.159. |
| `off_plane` | TINYINT | virtual | Added in v57 (GEN.117): 1 when `|position_z_mpc|` is 250 pc or more. Computed; nothing writes it. |

Indexes: the address (`shell_index`, `band_index`, `wedge_index`),
`luminosity_w`, and `idx_bright_stars_population` (`population`,
`off_plane`, `luminosity_w`, v59), which the Galaxy Map's tile picks walk
once per population so the old giants of the bulge and the thick disk
are listed beside the young stars that outshine them.

### `sector_stats`

Added in v53 (GEN.44, PERF.11); replaces v49's `bright_star_blocks`. One
row per sector address a backfill reached or a fill generated. See
"Per-sector stats (v53, GEN.44, PERF.11)" above.

| Column | Type | Null | Notes |
|---|---|---|---|
| `shell_index`, `band_index`, `wedge_index` | INT / INT / INT | PK | The sector's grid address. |
| `bright_level_sol` | DOUBLE | NOT NULL, default -1 | How deep its bright stars go: -1 untouched (follows `galaxy_shape.bright_star_min_luminosity_sol`), a positive L_sun for a backfill's floor, 0 generated. Indexed. |
| `level_before_fill_sol` | DOUBLE | nullable | The level a fill found, put back when the sector is deleted. |
| `bright_mass_sol` | DOUBLE | nullable | Added in v79 (GEN.187). The initial mass, in solar masses, down to which a backfill placed every living star of the unfilled sector; NULL where none did, so the sector holds only the galaxy scatter's mass pass (stars born at or above `galaxy_shape.bright_star_mass_limit_sol`). A new plan scatter clears it. |
| `relative_density` | DOUBLE | nullable | The galaxy model's density at the sector's center (`density.relative_density`). |
| `expected_systems` | DOUBLE | nullable | `relative_density` times `galaxy_shape.expected_system_count_at_density_1`. |
| `actual_systems`, `actual_stars` | INT | nullable | What the fill gave it; NULL before. |
| `mean_temperature_k`, `mean_luminosity_sol` | DOUBLE | nullable | The mean temperature and luminosity (L_sun) of its systems' stars. |
| `mean_age_gy` | DOUBLE | nullable | The mean age (Gy) of its stars; NULL with no stars (v55, DB.14). The Galaxy Map's hue comes from it. |
| `total_luminosity_sol` | DOUBLE | nullable | The summed luminosity (L_sun) of its stars; NULL with no stars (v55, DB.14). The Galaxy Map's brightness comes from it. No color is stored: the page works it out (MAP.128). |
| `filled_at` | TIMESTAMP(3) | nullable | When the fill saved it. |
| `updated_at` | TIMESTAMP | NOT NULL | |

### `sector_paths`

Added in v62 (GEN.123). One row per star system, rogue planet or
interstellar comet: the path it takes from where it is now to where it
leaves its sector (`physics/sector_path.py`), recomputed by
`sector_paths.compute_sector_paths` for every sector holding one of them
at each orbit update, and for the sectors a `galaxy` run created and the
sectors around them as the last step of the run (`--no-settle` skips it).
Never sector by sector during the run: a path depends on the neighbours
that exist, which varies with the order workers fill sectors.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `sector_id` | BIGINT UNSIGNED | NOT NULL | `sectors.id`, cascade delete. |
| `object_table`, `object_id` | VARCHAR(24) / BIGINT UNSIGNED | NOT NULL | `star_systems`, `rogue_planets` or `interstellar_comets`, and the row's id. Unique together; no foreign key. |
| `exited` | TINYINT(1) | NOT NULL | 1 when the body leaves the sector, 0 when the path was cut off (at rest, held by a mass, step limit). |
| `duration_years` | DOUBLE | NOT NULL | How long the crossing takes. |
| `computed_at` | TIMESTAMP | NOT NULL | |

### `sector_path_knots`

Added in v62 (GEN.123). The cubic Hermite spline knots of a `sector_paths`
row, in order (`position` 0 is where the body is now; the last knot is where
it leaves, the next sector's entry).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `path_id` | BIGINT UNSIGNED | NOT NULL | `sector_paths.id`, cascade delete. Unique with `position`. |
| `position` | INT | NOT NULL | The knot's index. |
| `t_years` | DOUBLE | NOT NULL | Years after the path starts. |
| `x_pc`, `y_pc`, `z_pc` | DOUBLE | NOT NULL | Position, galactic parsecs. |
| `vx_kms`, `vy_kms`, `vz_kms` | DOUBLE | NOT NULL | Velocity, km/s, galactic axes. |

### `generation_runs`

Added in v52 (DB.6). One row per `planetgen` run that changes the
galaxy (every subcommand except `check-math`, `system --output` and an
`--estimate-only` run), written when it starts and finished when it ends.
A galaxy is built by a series of runs (plan, sectors, scatter,
backfill), not by its seed alone, so these rows are what a rebuild
replays. Writing them never stops a run.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK, AUTO_INCREMENT | |
| `command` | VARCHAR(32) | NOT NULL | The subcommand (`plan`, `galaxy`, ...). |
| `run_seed` | BINARY(16) | nullable | The run's own 128-bit seed: what draws outside a galaxy unit (a random start's address, a one-off system) come from. |
| `galaxy_seed` | BINARY(16) | nullable | The galaxy seed when the run started; NULL for the plan that drew it. |
| `version_key`, `planetgen_version`, `python_version`, `platform` | CHAR(22) / VARCHAR | NOT NULL | The running code, as in `galaxy_shape`. |
| `started_at`, `finished_at` | TIMESTAMP(3) | `finished_at` nullable | |
| `outcome` | VARCHAR(16) | nullable | `ok`, `failed` or `interrupted`; NULL while it runs, or when it died without saying. |

The command line is in `generation_run_arguments` (v63, DB.13), which replaced a JSON text column of `generation_runs`: `run_id` (FK -> `generation_runs.id`, `ON DELETE CASCADE`), `position` (the argument's place, from 0) and `value` (VARCHAR(1024)), primary key `(run_id, position)`. It leaves out the `--mysql-*` and `--debug` options (`generate._run_argv`).

### `species` / `polities` / `system_owners` / `population_state`

Added in v44, filled only by the population pass
(`planetgen/population/model.py`); see "Population and politics (v44)"
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
| `energy_index`, `materials_index`, `information_index`, `medical_index`, `propulsion_index`, `defense_index` | TINYINT UNSIGNED | nullable, CHECK 0 to 7 | Added in v88 (POP.9). The six tech indexes; NULL until drawn. |
| `tech_level` | DECIMAL(3,2) | nullable | Added in v88 (POP.9). The weighted sum of the indexes. Indexed (`idx_species_tech_level`) for range search. |
| `tech_source` | VARCHAR(16) | nullable, CHECK IN ('generated','admin') | Added in v88 (POP.9). 'admin' marks a hand-edited species so a pass never redraws it. |
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

### `phenomenon_scatter` / `phenomenon_scatter_classes`

Added in v64 (GEN.100) and v77. `planetgen plan` places every black hole, neutron star, planetary nebula and supernova remnant of the galaxy, its hypervelocity stars, its quasar and its nucleus galaxy-wide in `phenomenon_scatter`, just before the bright-star scatter and before any sector is filled. A row is a point in its sector's cell and the seed the object is built from; filling the sector builds the object at that point and sets `built_at`, so what was placed stays where it was put. A plan re-run or reset empties the table.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK, AUTO_INCREMENT | |
| `shell_index`, `band_index`, `wedge_index` | INT / INT / INT | NOT NULL | The sector cell the point is in. Indexed together (`idx_phenomenon_scatter_address`). |
| `kind` | VARCHAR(24) | NOT NULL, CHECK | `black-hole`, `neutron-star`, `planetary-nebula`, `supernova-remnant`, `hypervelocity-star` or `quasar`. |
| `subtype` | VARCHAR(16) | nullable | `supermassive` for the nucleus black hole, NULL for the rest (the object draws its own class from `seed`). |
| `position_x_mpc`, `position_y_mpc`, `position_z_mpc` | BIGINT | NOT NULL | The point, in milliparsecs, galaxy frame. |
| `velocity_x_kms`, `velocity_y_kms`, `velocity_z_kms` | DOUBLE | nullable | A hypervelocity star's galactic velocity (km/s); NULL for the rest. |
| `seed` | BIGINT UNSIGNED | NOT NULL | The 63-bit seed the object is built from. |
| `built_at` | TIMESTAMP | nullable | When a sector fill built the object; NULL while it is only a point. |
| `epoch_unix` | DOUBLE | nullable | Added in v69 (GEN.137). When the position holds (see the v69 entry). |
| `mass_solar` | DOUBLE | nullable | Added in v80 (MAP.165). The solar mass a scattered black hole or neutron star was drawn with. |
| `map_luminosity_sol` | FLOAT | nullable | Added in v84 (DB.23). The luminosity the Galaxy Map draws the point with (see the v84 entry). |

`idx_phenomenon_scatter_class (kind, subtype, mass_solar)` (v82, PERF.68) serves the Galaxy Map's coarse tiles.

`phenomenon_scatter_classes` (v77) counts the scatter per class so the Phenomena table counts and pages it without scanning its rows: `kind`, `subtype` (`''` where a class has none), `placed` (written when the scatter ends, `store.record_phenomenon_scatter_classes`) and `built` (goes up as `store.mark_phenomena_built` stamps rows); primary key `(kind, subtype)`. A plan re-run or reset empties it with the scatter.

### `sector_system_counts`

Added in v82 (PERF.74). How many star systems each sector holds, as of the last refresh (`query.refresh_sector_system_counts`), so the Sectors list and its sorts by systems and density read this table instead of counting every system per request. A sector with no systems has no row.

| Column | Type | Null | Notes |
|---|---|---|---|
| `sector_id` | BIGINT UNSIGNED | PK | FK-like to `sectors.id` (no constraint). |
| `system_count` | INT UNSIGNED | NOT NULL | The systems it holds. Indexed with `sector_id` (`idx_sector_system_counts_count`). |

### `object_timing_stats` / `object_timing_samples`

Added in v85 (PERF.92). How long each kind of generated object takes to build, exclusive of the objects inside it (`planetgen/generation/object_timing.py`), saved under the run (`run_id` -> `generation_runs.id`, `ON DELETE CASCADE`; `sector_id` -> `sectors.id`, `ON DELETE CASCADE`, NULL for work not tied to one sector). Written unless `PLANETGEN_GENERATION_STATS` is set to 0.

`object_timing_stats` has one row per object kind a run built in a sector:

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK, AUTO_INCREMENT | |
| `run_id`, `sector_id` | BIGINT UNSIGNED | `sector_id` nullable | As above. |
| `kind` | VARCHAR(24) | NOT NULL | The object kind (`object_timing.KINDS`). |
| `object_count` | INT UNSIGNED | NOT NULL | Objects of that kind built. |
| `total_us`, `min_us`, `max_us` | BIGINT UNSIGNED | NOT NULL | Their own time in microseconds: the total, the fastest and the slowest. |
| `histogram` | VARCHAR(255) | NOT NULL | 24 comma-separated counts in doubling-microsecond buckets (`object_timing.HISTOGRAM_BUCKETS`). |

A run's figures for a kind are the sum of its sector rows.

`object_timing_samples` keeps single objects with their full timing: `id`, `run_id`, `sector_id`, `kind`, `ordinal` (the object's number within its kind in that sector), `own_us` (its own time), `total_us` (including the objects inside it) and `sampled` (1 for the repeatable 1-in-1,000 sample, a hash of the unit, kind and `ordinal`; 0 for the extra rows a benchmark run asks for).

### `id_blocks` / `id_counters`

Bookkeeping that hands out ids. Neither describes a generated object, and `planetgen reset` keeps both.

`id_blocks` (v45, PERF.13): `table_name` (PK, a table in `store.ID_BLOCK_TABLES`) and `next_id` (BIGINT UNSIGNED), the next unreserved row id. `store._allocate_id` moves `next_id` up by a block at a time (never below the table's own `MAX(id) + 1`) on a separate autocommitted connection, so writers don't queue behind each other, and the writer names those ids in its INSERTs.

`id_counters` (v78, DB.20): `kind` (VARCHAR(8)), `scope` (VARBINARY(16)) and `next_value` (BIGINT UNSIGNED), primary key `(kind, scope)`. `kind = 'sector'`: `scope` is a birth sector's 5-byte address and the number is the next run-time serial there. `kind = 'system'`: `scope` is a system's 10-byte ID and the number the next body number in it. `kind = 'sdel'` (PERF.38), empty `scope`: the deletion epoch, one more with every sector deleted, so the Galaxy Map's stamp tells a deletion without counting the placed sectors. A row only moves up, on its own autocommitted connection (`store._reserve_counter`), and survives the deletion of the object or sector, so a number is never given twice.

### `nebula_shape_balls`

Added in v56 (GEN.75). The 4 to 8 metaballs of a nebula's shape, one row each, in order (`galaxy/nebula_shape.py`; the single shape values are the `nebulae.shape_*` columns): `nebula_id` (FK -> `nebulae.id`, `ON DELETE CASCADE`), `ball_index`, `center_x`, `center_y`, `center_z` and `radius` (DOUBLE, in the shape's own units), primary key `(nebula_id, ball_index)`. Absent for a nebula saved before v56; `queryDb.nebula_shape` then draws the same shape from the nebula's own properties.

### `alembic_version`

Alembic's own record of the revision a database is at (DB.11), declared in `schema.sql` as Alembic would create it so that a new and a migrated database have the same shape. One row, `version_num` (VARCHAR(32), PK): the schema version as four digits (`0091`). `alembic_runner.upgrade` keeps it in step with `schema_migrations`.

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

### `galaxy_components`

The discrete stellar components of a galaxy beyond the disk, one row each, of every kind in one table (v97, GEN.159): the globular clusters (GEN.159) and the satellites (GEN.237: the Magellanic Clouds, Sagittarius and eight dwarfs, with an `exponential` or `plummer` profile), and streams later. A component is a density, not a list of stars: a profile of its mass about a centre, from which a sector's stars are drawn when the sector is generated (`docs/design/globular-clusters.md`, `docs/design/intergalactic-and-halo-contents.md` section 5.5). The `king` profile is `planetgen/physics/king.py`: ten tables for central potentials 3 to 12, solved once and scaled by length; the `exponential` and `plummer` profiles are `planetgen/physics/component_profiles.py`, each cut at `tidal_radius_pc` (the column is set for every profile, `king_w0` for `king` alone). `planetgen/galaxy/components.py` has the `Component` record and `store.replace_galaxy_components` / `get_galaxy_components` write and read it. `planetgen plan` draws the clusters from the galaxy seed with `planetgen/generation/clusters.py` before it walks the outline (GEN.160), replacing the `cluster` rows, and places the Milky Way's eleven satellites with `planetgen/generation/satellites.py` (GEN.237), replacing the `satellite` rows; the read-back shape carries both as `GalaxyShape.components` (`planetgen/galaxy/component_field.py`), which adds their density to the sector gate.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `kind` | VARCHAR(12) | NOT NULL | `cluster`, `satellite` or `stream` (CHECK). |
| `name` | VARCHAR(255) | NOT NULL | Unique within its kind. |
| `catalogue_id` | VARCHAR(40) | nullable | The catalogue's own name ("NGC 104"), so a later import of a real catalogue can match rows. |
| `source` | VARCHAR(12) | NOT NULL | `catalogue` for a real object, `generated` for one drawn from the galaxy seed (CHECK). |
| `parent_component_id` | BIGINT UNSIGNED | nullable | Self-reference: a cluster inside a satellite, or a satellite of a larger galaxy. `ON DELETE SET NULL`. |
| `center_x_pc`, `center_y_pc`, `center_z_pc` | DOUBLE | NOT NULL | The centre in the galaxy frame (origin the galactic centre), parsecs. |
| `cell_key` | BIGINT UNSIGNED | generated, stored | The Morton key of the centre's 0.25 pc cube, like every other table with a galaxy-frame point (`planetgen/galaxy/cellkey.py`); NULL outside plus or minus 262,144 pc. Covering index `(cell_key, center_x_pc, center_y_pc, center_z_pc)`. |
| `velocity_x_kms`, `velocity_y_kms`, `velocity_z_kms` | DOUBLE | NOT NULL | The velocity in the galaxy frame, km/s, at `epoch_unix`. |
| `epoch_unix` | DOUBLE | NOT NULL | When the position and velocity hold. |
| `mass_solar` | DOUBLE | NOT NULL | Total mass. A cluster's systems number this over 0.4. |
| `half_mass_radius_pc` | DOUBLE | NOT NULL | The radius holding half the mass. |
| `profile` | VARCHAR(12) | NOT NULL | `king`, `exponential` or `plummer` (CHECK). |
| `king_w0` | DOUBLE | nullable | The King central potential, 3 to 12; set exactly when `profile` is `king` (`Component` checks this, not the table). |
| `tidal_radius_pc` | DOUBLE | nullable | The King tidal radius, never beyond the cluster's Jacobi radius; set exactly when `profile` is `king`. |
| `age_gy` | DOUBLE | NOT NULL | The stars' age. |
| `metallicity_fe_h` | DOUBLE | NOT NULL | The stars' [Fe/H] in dex, which a cluster star takes (GEN.158 left it NULL until this table). |
| `core_collapsed` | TINYINT(1) | NOT NULL, default 0 | Whether the core has collapsed (CHECK 0 or 1). |
| `seed` | BIGINT UNSIGNED | NOT NULL | The seed its stars are drawn from. |
| `created_at` | TIMESTAMP | NOT NULL | |

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

`generation_runs` owns its `generation_run_arguments`, `object_timing_stats`
and `object_timing_samples` rows (`ON DELETE CASCADE`); the two timing
tables also point at `sectors`. `nebula_shape_balls` hangs off `nebulae`.
`phenomenon_scatter` and `phenomenon_scatter_classes` carry a sector address
but no foreign key, like `sector_stats`.

`galaxy_shape`/`galaxy_shell` stand apart from the tree above —
neither has a foreign key to `sectors` or anything else. They describe the
galaxy as a whole (its shape parameters, and which addresses could hold
content), not any individual sector; `sectors.shell_index` is the shell
both an actual generated sector and a `galaxy_shell` row are
independently expressed in, not an FK
relationship.


### v69 (GEN.137): `phenomenon_scatter.epoch_unix`

`phenomenon_scatter.epoch_unix` (DOUBLE, nullable): the orbit epoch a scattered hypervelocity star's position holds at, stamped when the plan draws it (NULL when no orbit update has run yet, meaning the database's orbit epoch). The built star system inherits it as its own `epoch_unix`, so the next orbit update flies it on from there.

### v108 (DB.72): band numbers are INT in the scatter tables

Migration 0108 widens `band_index` from SMALLINT to INT in `bright_stars`, `sector_stats` and `phenomenon_scatter`. A shell's bands run from `-(h + 1)` to `h + 1`, which passes SMALLINT's 32,767 at shell 20,861 and reaches 102,942 at shell 65,535, the outer shell of the 262 kpc reach, so the whole-sphere scatter stopped at 97 percent with "Out of range value for column 'band_index'". `sectors`, `galaxy_shell` and `galaxy_block_stats` were INT already. Each table is rebuilt in place (its address index with it) and no value is touched; a column that is already INT is skipped, so a stopped run can be started again. The version is detected from `bright_stars.band_index` being INT.

### v107 (GEN.526): a satellite's young and intermediate stars

Migration 0107 lets `bright_stars.population` take `satellite_young` and `satellite_intermediate`, the stars a satellite galaxy's star formation history puts in the last 0.1 Gy and in the 3 Gy before that (`galaxy/satellite_history.py`); `satellite` is now the rest, 3 Gy and older. The column widens to VARCHAR(24) for the longer name (an in-place change that copies no row) and no row is touched. The version is detected from the check constraint `chk_bright_stars_population` listing `satellite_young`.

### v106 (DB.65): the cell key holds the last cube at the box edge

Migration 0106 recomputes the generated `cell_key` of the 14 tables that hold one (`star_systems`, the phenomenon tables, `bright_stars`, `phenomenon_scatter`, `facilities`, `nebula_field`, `galaxy_components`, `asteroid_fields`) with `LEAST(cube number, 2097151)` on each axis. A coordinate just under the 262,144 pc reach whose rounded cube number came to 2^21 took the key of the opposite edge (the 21-bit field dropped the top bit); it now takes the last cube's. The version is detected from `LEAST` in `star_systems.cell_key`'s generation expression. Altering a stored generated column rebuilds the table and its covering index, so the migration runs once per table with a progress line.

### v105 (GEN.324): satellite stars

Migration 0105 lets `bright_stars.population` take `satellite`, the stars of a satellite galaxy (the Clouds, Sagittarius, the dwarfs) that the component pass places (`generation/component_pass.py`); a globular cluster's stars were already `cluster` (v94). Nothing else changes and no row is touched. The version is detected from the check constraint `chk_bright_stars_population` listing `satellite`.

### v104 (DB.53): the orbit tree of a built system

Migration 0104 gives a built multi-star system the storage the scatter side already had, so the items of the multi-star train (GEN.242, GEN.264, GEN.269, GEN.273, GEN.274, GEN.275) need no revision of their own. `orbit_nodes` is keyed by `star_system_id` and `node_no` and has the columns of `scatter_nodes` by name and unit (`parent_node_no`, `slot`, `mass_solar`, `separation_km`, `eccentricity`, `inclination_deg`, `ascending_node_deg`, `arg_periapsis_deg`, `mean_anomaly_deg`, `epoch_unix`, `period_years`, `kind` of `close`, `wide`, `compact` or `group`, `interaction`, `flags`), so a fill copies a node column for column, and adds the disk descriptor (`disk_kind`, `disk_mass_solar`, `disk_inner_au`, `disk_outer_au`, all nullable; a disk holds no bodies of its own). A node is deleted with its system (foreign key, ON DELETE CASCADE). `stars` gets `parent_node_id` (the node a star sits in), the member columns of `scatter_members` (`state`, `flags`, `offset_x/y/z_au`, `velocity_x/y/z_kms`, `decay_age_yr`, `epoch_unix`) and the same disk descriptor; `planets` gets `host_node_id` (the node a circumbinary planet orbits, NULL for a planet of one star). A star is now a layout of `planetgen/db/point_columns.py` (system frame, offset in au, velocity in km/s, mass in kg, epoch, no clock of its own). A `close` or `wide` pair converts to one node (number 0, no parent) from its stored mutual orbit and the mass of its two stars, its stars take `parent_node_id` 0, and every star's offset is the stored binary position (zero for a single star); velocity and epoch stay NULL until the code that moves members fills them. Every step is skipped when done and the version is detected from `planets.host_node_id`. GEN.242 fills the disk descriptor of `stars` for a protostar or pre-main-sequence star (`disk_kind` of `class-0`, `class-i` or `class-ii`, the mass in solar masses, the inner and outer edge in AU; NULL for every other star); `yerkes_class` is `protostar` or `PMS` for those stars.

### v103 (GEN.336): the movement vector replaces the scalar galactic speed

Migration 0103 drops `galactic_orbital_speed_kms` from `stars`, `black_holes`, `neutron_stars`, `nebulae`, `supernova_remnants`, `rogue_planets`, `interstellar_comets` and `asteroid_fields`, and `binary_galactic_orbital_speed_kms` from `star_systems`. The vector of v101 (`velocity_x/y/z_kms`, galactic axes) is the one record of how an object moves; a speed is its length (`store._speed_kms`). A free object is written with the tangent of the galaxy's rotation at its place at the speed it was generated with (`store._free_vector`), taken at `store._DEFAULT_PLACE_PC` (the Sun's distance from the centre on the +x axis, where the generators work the speed out) while it has no place; the orbit update turns the vector with the position; the mass columns of v101 are filled when the object is made (a comet nucleus from its size, a nebula or remnant from its bounding sphere and density, see `physics/bulk_mass.py`). A star and a close pair have no vector of their own: their speed is the length of their system's. The version is detected from the speed column's absence on a database that has `sectors.max_mass_id`. Each column is dropped only when present.

### v102 (DB.55): the point of every kind of object

Migration 0102 adds, nullable and backfilled, the columns the layout in `planetgen/db/point_columns.py` says an object lacks, so the code items that fill them (GEN.336, GEN.337, GEN.289) need no revision of their own. A movement vector (`velocity_x/y/z_kms`, galactic axes) goes on black holes, neutron stars, rogue planets, interstellar comets, asteroid fields, nebulae, supernova remnants, asteroid belts and the quasar: a free object takes the tangent of the galaxy's rotation at its place at its old `galactic_orbital_speed_kms`, and the quasar and a belt take a stated zero. A mass (`mass_kg`) goes on interstellar comets, asteroid fields, nebulae, remnants, belts, comets and facilities: a comet nucleus takes its sphere's mass at the density of 67P (`tuning.COMET_NUCLEUS_DENSITY_KG_M3`), a nebula or remnant the mass of its bounding sphere at its particle density (an upper bound for a cloud that fills part of its sphere, `physics/bulk_mass.py`), a field or belt a stated zero that the gravity sum skips, and a facility the mass of what it sits on. Belts also get `position_x/y/z_km` (the system frame, at the centre), belts and the quasar get `epoch_unix` and `next_update_due`, and scatter members get `epoch_unix`. Asteroid fields get `formed_at` and `origin_event_id` (GEN.110; no foreign key until the event table, DB.56). Planets and moons get `parent_type` (`star`, `planet` or `system`, with checks), `parent_id` and `dynamical_state` (`bound` or `unbound`), backfilled from `star_id` and `planet_id`. `sectors` gets `max_mass_kg` and `max_mass_id` (the heaviest object and its uid, for GEN.109), NULL until worked out. A star has no columns of its own: it moves with its system, whose point, vector and epoch it reads. The scalar `galactic_orbital_speed_kms` stays until GEN.336's drop revision. Every step is skipped when done, and the version is detected from `sectors.max_mass_id`.

### v101 (DB.54): the scatter limits by place

Migration 0101 adds three columns to `galaxy_shape`: `core_min_luminosity_sol` and `core_mass_limit_sol` (`DOUBLE`, NULL until set, and then the core follows the baseline's limits times 1.1) and `halo_min_luminosity_sol` (`DOUBLE NOT NULL`, default 500 solar luminosities), with the CHECK `chk_galaxy_shape_scatter_limits` (the halo floor at least 50, the core's two more than 0 when set), by which the version is detected. The halo has no mass limit of its own. A galaxy planned before it takes the defaults, which change nothing until GEN.333 reads them. Each step is skipped when it is already done.

### v100 (GEN.281): the cell star ceiling defaults to 5,000

Migration 0100 sets the default of `galaxy_shape.cell_star_ceiling` to 5,000 systems (it was 50,000, which took 14 to 17 minutes a cell at 20 ms a system), moves a galaxy still at 50,000 to 5,000 (one planned with another value keeps it), and adds the CHECK `chk_galaxy_shape_cell_star_ceiling` (at least 1), by which the version is detected. No other row is touched. Each step is skipped when it is already done.

### v99 (GEN.280): the nuclear structures

Migration 0099 renames `galaxy_shape.core_amplitude` to `nuclear_amplitude` (DOUBLE, NOT NULL, default 0): how much of the observed nuclear star cluster (a Hernquist profile, 2.5e7 solar masses) and nuclear stellar disc (an exponential disc 90 pc long, 1.5e9) the galactic centre holds, scaled to the bulge, 1 being the observed pair and 0 none (`planetgen plan --nuclear-amplitude`, default 1). It replaces ADM.49's Gaussian core, which is gone. A galaxy planned before it had no such structures, so every row is set to 0 and its density is unchanged. `bright_stars.population` becomes VARCHAR(16) (an in-place change that copies no row) and takes `nuclear_cluster` and `nuclear_disc`. The version is detected from `galaxy_shape.nuclear_amplitude`.

### v98 (GEN.281): the cell star ceiling and the fill luminosity cap

Migration 0098 adds two columns to `galaxy_shape`, both `DOUBLE NOT NULL` and set for the whole galaxy: `cell_star_ceiling` (default 5,000 systems) and `fill_luminosity_cap_sol` (default 100 solar luminosities). A cell that expects more systems than the ceiling is filled only with the stars at or above the cap; the rest is unresolved, and its count and mass are worked out from the address by `planetgen/galaxy/resolution.py`, never stored. Nothing else changes and no row is touched. It is skipped when the columns are already there. The version is detected from `galaxy_shape.cell_star_ceiling`.

### v97 (GEN.159): `galaxy_components`

Migration 0097 creates the empty `galaxy_components` table (see its section above): the components of a galaxy beyond the disk, as densities with a centre, velocity, mass, age and [Fe/H]. It changes no other table. It is skipped when the table is already there, so a stopped run can be started again. The version is detected from the table.

### v96 (GEN.153): the magnetar pulsar type

Migration 0096 replaces the CHECK on `neutron_stars.pulsar_type` (`neutron_stars_chk_1`, the name MySQL gave it) with one that also allows `magnetar`: a neutron star younger than a magnetar's active life may be born one (`tuning.MAGNETAR_*`). The quadrant CHECK beside it stays `neutron_stars_chk_2` and no row is touched. The version is detected from the text of `neutron_stars_chk_1`.

### v95 (GEN.158): the stored [Fe/H]

Migration 0095 adds `metallicity_fe_h` (DOUBLE, nullable) to `stars` and to `bright_stars`: the star's metallicity in dex against the Sun, drawn from its population and its distance from the galactic centre (`planetgen/physics/metallicity.py`). A scattered star's value is drawn at scatter time, on a stream of its own seed, and the fill that builds its system keeps it; a binary's secondary shares its primary's. It is NULL on rows made before v95 and on a cluster star until the cluster table (GEN.159) gives it its cluster's value; a star built from a NULL draws its own as before. The columns go on the end of each table, so MySQL adds them without copying the rows. The version is detected from `bright_stars.metallicity_fe_h`. The sector page's map shows it in a star's info panel.

### v94 (GEN.234): the halo and cluster populations

Migration 0094 widens the CHECK on `bright_stars.population` from the four disk and bulge populations to also allow `halo` and `cluster`: the stellar halo replaces the flat density floor, so a pre-placed star can belong to it. Nothing is renumbered. The version is detected from that constraint's text (`store._check_marker`).

### v93 (GEN.230, GEN.233, GEN.223, GEN.260, DB.34, DB.35): 96-bit IDs, reach, scatter points, cell keys

One migration (0093) for the changes that renumber or reshape the same rows. It can be stopped and started again; the last step (the reach column) marks it done.

- **Object IDs are 96 bits** (GEN.230): the birth sector's shell, biased band and wedge are 19, 18 and 19 bits (the serial is 28 bits, the body 12), the layout is fixed for the whole galaxy (`object_uid.LAYOUT`; the outline no longer fits one), and every `uid` is `BINARY(12)`. `sectors.uid` (the designation) is the 56-bit sector field, and an unplaced sector counts up from `1 << 56`. `id_counters` scopes follow (a sector's is 7 bytes, a system's 12). The migration builds the new IDs in a `uid_new` column in pieces of 20,000 rows, then swaps the column in one rebuild.
- **Reach** (GEN.233): `galaxy_shape.reach_shell_index` (INT, default 65,535 = 262 kpc). The grid can be named out to that shell and a little past it (`spherical.MAX_ADDRESSABLE_SHELL`); `GalaxyBounds.qualifies` is the disk outline that limits what `generate` builds, `contains` the reach.
- **Scatter points** (GEN.223): `bright_stars` and `phenomenon_scatter` keep their point like every other object: `position_x/y/z_mpc` (DOUBLE, galaxy frame), `velocity_x/y/z_kms`, `epoch_unix`, an indexed `next_update_due`. The nightly orbit update moves every unbuilt point (`store.advance_scatter_points`) and refiles it when it crosses a cell; the quasar stays at the center and a built row is left to the object built from it.
- **Companions of a scattered system** (GEN.260): `member_count` and `system_mass_solar` on both scatters (the row's point is the barycentre), and `scatter_members` / `scatter_nodes` (docs/design/multi-star-systems.md sections 5, 13, 14), written by the scatter passes (GEN.262, `planetgen/generation/scatter_companions.py`: a star row always draws its tree, a black hole, neutron star or planetary nebula row at the table's odds, and deleting a scatter row or clearing a scatter deletes its members and nodes): a member's `kind`, `state`, `flags`, masses, seed, and for a young group its offset, velocity and `decay_age_yr`; a node's orbit, `kind` (`close`, `wide`, `compact`, `group`), `interaction` and `flags`; `slot` is SMALLINT on both. No foreign keys (the parent is one of two tables). The nightly update (GEN.265) moves only the row's point, the barycentre; it never writes a member or node row, and a member's place is computed on read from the node phase and epoch (`planetgen/physics/system_tree.py`, `planetgen/db/scatter_systems.py`).
- **Cell key** (DB.34, DB.35; docs/design/spherical-region-queries.md): `cell_key` BIGINT UNSIGNED, a stored generated column holding the Morton key of the object's 0.25 pc cube over plus or minus 262,144 pc (NULL outside), and the covering index `idx_<table>_cell (cell_key, x, y, z)`, on `star_systems`, the eight phenomenon tables, `facilities`, `nebula_field`, `bright_stars` and `phenomenon_scatter`. MySQL computes it from the position columns (`galaxy/cellkey.py` has the same function, tested against it), so it cannot disagree with them. `star_systems` keeps its position relative to its sector, so it also has `galaxy_x/y/z_pc`, its point in the galaxy frame, written by every save and by the orbit update.

### v92 (MAP.205): `galaxy_block_stats` on shells

Migration 0092 makes `galaxy_block_stats` again, for the Galaxy Map's block ladder on shells (`galaxy/drill.py`, `docs/design/spherical-block-ladder.md`). It holds one row per block of 243, 81, 27, 9 or 3 shells that holds a generated sector: `m` (SMALLINT UNSIGNED), `shell_index`, `band_index` and `wedge_index` (INT, the block's own address; together the primary key), `sector_count`, `systems` (BIGINT UNSIGNED), `expected_systems` (DOUBLE), `stars` (BIGINT UNSIGNED), `age_sum` (DOUBLE, each sector's mean age times its stars, so a mean is `age_sum / stars`) and `luminosity_sol` (DOUBLE). There is no `ladder_step` column any more: the ladder is 3-ary and a block is named by its shell, band and wedge. The figures come from `sectors` and `sector_stats`; a block with no generated sector has no row. `query.galaxy_block_stats(conn, m, cells)` reads a window's cells in one statement. `db/block_stats.py` keeps the table: a sector saved on its own (a visit) or deleted updates the blocks above it (`drill_chain_of`), and a `galaxy` run updates the blocks above all the sectors it created once at its end. The migration fills the table from the sectors already there, a shell at a time, and skips itself when the table exists, so a stopped run can be started again.

### v91 (DB.32): sectors on spherical shells

The sector grid moved from the cylinder's ring, layer and slot to the sphere's shell, band and wedge (`docs/design/spherical-sector-coordinates.md`). Migration 0091 renames `ring_index`, `layer_index` and `ring_slot_index` to `shell_index`, `band_index` and `wedge_index` on `sectors`, `bright_stars`, `sector_stats` and `phenomenon_scatter`, renames `galaxy_shape.outer_ring_index` to `outer_shell_index`, replaces `galaxy_layer` and `galaxy_column` with `galaxy_shell` and drops `galaxy_block_stats` (the block ladder is rebuilt on shells with MAP.205 and MAP.206).

It converts the data; nothing is deleted. Every star system, phenomenon and facility keeps its place in the galaxy (a system's stored offset is worked out again from its place and its new sector's centre), its velocity and its epoch. The sector rows are rebuilt: one for every cell that holds an object and one for every cell whose centre lay inside a generated sector; a new sector takes the name of the old sector holding its centre (second and later ones get " 2", " 3" ...), and a system's location text takes the new sector's name. Object IDs are renumbered for the new birth sector (the serial is the object's rank there, a body keeps its body number); `id_counters` is rebuilt to match. The bright stars and phenomenon scatter get the address of the cell holding their position, `sector_stats` is counted again for the generated sectors (the untouched and partly backfilled rows are dropped, as a new scatter does), `sector_paths` of the old sectors go with them and come back with the next orbit update, and `galaxy_shell` is built again from the stored shape. The address indexes of `bright_stars` and `phenomenon_scatter` are dropped for the bulk update and built again after: about 35 microseconds a row, so a hundred million scatter rows take about an hour.

### v90 (MAP.151b): `galaxy_block_stats.ladder_step` and the 3-ary ladder

`galaxy_block_stats.ladder_step` (TINYINT UNSIGNED, NOT NULL, default 3) is the step between the block ladder's levels (`galaxy.drill.DRILL_STEP`) the row's totals were computed for, and `query.galaxy_block_stats` reads only rows of the code's own step. The ladder went from 243, 27 and 3 sectors a side to 243, 81, 27, 9 and 3, and the same `m` names another block on each, so the migration empties the table and fills it again for the new ladder (the pass migration 0086 makes). The Galaxy Map reads a window of blocks bigger than a sector from this table in one query (`query.galaxy_window`).

### v89 (DB.28): facility types

`facility_types`, `facility_type_fields` and `facility_field_values` are new, and `facilities` gains `type_id`, `affiliation_polity_id`, `affiliation_name`, the five Green/Yellow/Red ratings (`rating_crime`, `rating_housing`, `rating_resources`, `rating_maintenance`, `rating_health`) and `ratings_source`, `ratings_note`, `ratings_set_at`. All NULL on facilities saved earlier. The migration skips each step already done, so a stopped run can be started again.

### v88 (POP.9): `species` tech level

`species` gains the six tech indexes (`energy_index` to `defense_index`, each 0 to 7), `tech_level` (DECIMAL(3,2), indexed) and `tech_source` ('generated' or 'admin'). All NULL for species saved earlier. The migration skips each step already done, so a stopped run can be started again.

### v87 (DB.26): habitability columns

`planets` and `moons` gain `induced_moment_a_m2` (after `dipole_class`), `flux_variation` (after `energy_flux_w_m2`), `phi4_water`, `tier_water`, `phi4_suitability` and `tier_suitability`, and `phi4_chemistry` / `tier_chemistry` are renamed `phi4_air` / `tier_air` (GEN.178, GEN.148, GEN.207). `nebulae` gains `age_myr`, `q_total` (ionizing photons per second) and `hosts` (host stars as counts by spectral type) after `extinction_av` (GEN.150). All NULL on rows generated earlier; no compatibility name is kept for the old columns. The migration skips each step already done, so a stopped run can be started again.

### v86 (DB.25): `galaxy_block_stats` (cylinder form, dropped in v91)

The table was dropped in v91 (DB.32) and made again in v92 on shells (see "v92" above). What it held, as it was:

`galaxy_block_stats` holds what a drill-down stage shows for each block of the Galaxy Map's block ladder above the sector (`galaxy.drill`, 243, 81, 27, 9 and 3 sectors a side): one row per block `(m, ring_index, wedge, slab)` that holds a generated sector, with `sector_count`, `systems`, `expected_systems`, `stars`, `age_sum` (each sector's mean age times its stars, so a mean is `age_sum / stars`) and `luminosity_sol` and `ladder_step` (v89). The figures come from `sectors` and `sector_stats`; a block with no generated sector has no row. `query.galaxy_block_stats(conn, m, cells)` reads a window's cells in one statement, and a test holds every row equal to what `query.galaxy_stage` computes for the same block. `db/block_stats.py` keeps it: a sector saved on its own (a visit) or deleted updates the blocks above it, a `galaxy` run updates the blocks above all the sectors it created once at its end (the stage "Update the block totals"), and `planetgen plan --fill-block-stats` rebuilds every row. The migration fills the table from the sectors already there. The ladder is 3-ary since MAP.151b; v90 emptied the table and filled it again for it. No foreign keys: a block is not a row of its own.

### v85 (PERF.92): `object_timing_stats` and `object_timing_samples`

`object_timing_stats` has one row per object kind a run built in a sector: `object_count`, `total_us` (the objects' own time in microseconds, the time spent building objects inside them taken out), `min_us`, `max_us` and `histogram` (24 comma-separated counts in doubling-microsecond buckets, the first for under a microsecond). A run's figure for a kind is the sum of its sector rows (`object_timing.run_totals`); the mean over every run is `object_timing.type_means`. `object_timing_samples` keeps single objects with their full timing: `ordinal` is the object's number within its kind in that sector, `sampled` is 1 for the repeatable 1-in-1,000 sample (a hash of the sector, kind and ordinal, so the same galaxy and seed sample the same objects) and 0 for the extra rows `PLANETGEN_OBJECT_TIMING_FULL=1` asks for in a benchmark run. Both tables cascade from `generation_runs` and `sectors`, and a fingerprint leaves them out. The migration only creates them empty.

### v84 (DB.23, GEN.201): `phenomenon_scatter.map_luminosity_sol`

`phenomenon_scatter.map_luminosity_sol` (FLOAT, nullable) is the luminosity, in solar luminosities, of a main-sequence star of the row's `mass_solar` by the Galaxy Map's relation (`stellar_evolution.map_luminosity_sol`, GEN.200), so a faint but massive neutron star or black hole is shown where a star of its mass would be. It is NULL for a row with no mass (planetary nebulae, supernova remnants, hypervelocity stars, the quasar). A scatter writes it with each row; `planetgen plan --fill-map-luminosity` fills it, in id-range batches of 200,000 that readers never wait for, in a galaxy scattered before v84 (resumable: it only fills rows still empty). The migration adds the column empty. No index. FLOAT keeps about seven digits and the server prints six, so a value read back agrees with the function to 1e-5 relative.

### v83 (DB.24): `nebula_field`

`nebula_field` holds the galaxy's field clouds of classes M and N (giant molecular clouds and dark clouds; P and Q stay implicit in the seed), one row each, written by `planetgen plan` from the galaxy seed over the 50 pc field cells (`galaxy/nebula_field.py`). Primary key `(cell_i, cell_j, cell_k, n)`: the field cell and the cloud's place in that cell's draw, so a row is found again from the seed alone. `nebula_class` (CHAR(1)), `radius_ly` and `center_x_pc`, `center_y_pc`, `center_z_pc` (DOUBLE, the exact numbers the seed gives) describe the cloud. `uid` (BINARY(12), nullable, indexed) and `is_stored` (TINYINT, default 0) stay empty until a sector's save stores the cloud as a `nebulae` row, which gives it its object ID: working an ID out for every cloud would cost a nebula shape each (about 0.1 s). The table is seed-fixed and rebuilt by every plan, by `plan --redo-scatters phenomena` and by `plan --redo-nebula-field apply`. The migration only creates the empty table.

### v82 (PERF.68, PERF.70, PERF.74): indexes and `sector_system_counts` for the pages that timed out

`phenomenon_scatter` gains `idx_phenomenon_scatter_class (kind, subtype, mass_solar)`, which the Galaxy Map's coarse tiles use for the few biggest classes instead of reading an address range (3 million scattered rows: 12.7 s to 0.25 s). `star_systems.idx_star_systems_quadrant` becomes `(quadrant, name, id)` and `idx_star_systems_binary_name (is_binary, name, id)` is new: the Systems list sorted by Octant, Binary or Sector counts group after group and reads the page from the group it lands in (`query._systems_by_group`) instead of sorting every system. `sector_system_counts (sector_id PK, system_count)` stores each sector's system count, with an index on the count; the migration fills it, and `query.refresh_sector_system_counts` (one `GROUP BY` on the sector index, run in the background through `countcache`, not more often than every 30 s) keeps it. A sector with no systems has no row. The Sectors list and its sorts by systems and density read it.

### v81 (ADM.49): `galaxy_shape.arm_level`, `galaxy_shape.core_amplitude`

`galaxy_shape.arm_level` (DOUBLE, NOT NULL, default 1) is the thin disk's level between its arm crest and its inter-arm trough: the arm factor is `arm_level * (1 + arm_amplitude * cos)`, so the crest density is `arm_level * (1 + arm_amplitude)` and the trough `arm_level * (1 - arm_amplitude)`. `galaxy_shape.core_amplitude` (DOUBLE, NOT NULL, default 0) was the galactic core's centre density, a Gaussian added to the bulge; v99 (GEN.280) replaced it with `nuclear_amplitude`. The defaults are the model before v81, so an upgraded galaxy keeps its density. The command line and the Generate page take the arm density, inter-arm density and (since v99) the nuclear amplitude instead (`planetgen plan --arm-density 1.4 --interarm-density 0.6 --nuclear-amplitude 1`), and `galaxy.density.arm_terms` turns the first two into `arm_amplitude` and `arm_level`.

### v80 (MAP.165): `phenomenon_scatter.mass_solar`

`phenomenon_scatter.mass_solar` (DOUBLE, nullable): the solar mass a scattered black hole (the central one too) or neutron star was drawn with, inside its class range and above the scatter's cut. The sector build gives the object exactly this mass, and the Galaxy Map sizes the point from it (a straight line in log mass from 0.2 at one solar mass to 1.0 at 1e8), so a heavier object is never drawn smaller. NULL for quasars, hypervelocity stars, nebulae and remnants, and on rows scattered before v80 (the map then sizes those by mass class).

### v71 (GEN.167): `galaxy_shape.phenomenon_min_mass_solar`

`galaxy_shape.phenomenon_min_mass_solar` (DOUBLE, nullable): the phenomenon scatter's mass cut. The scatter places only the neutron stars and black holes at or above it (with 20, only the intermediate-mass black holes: about 2.7e5 rows at a 0.1% share, about 800 at the 3e-6 of GEN.259, instead of 1.17e9); a sector fill draws the share below it on the stream `{seed}:phenomena-fill:{ring}:{layer}:{slot}`, so the two add up to the whole population. `planetgen plan --phenomena-only --phenomenon-min-mass M` re-scatters at a new cut; filled sectors already hold every mass and are left out. See docs/design/phenomenon-scatter-mass-cut.md.

### v72 (GEN.86): stellar activity and magnetic fields

`stars` gains `log_lx_lbol`, `l_xuv_w`, `xuv_saturated`, `flare_n33_per_yr`, `flare_alpha` and `xuv_fluence_j`; `planets` and `moons` gain `magnetic_moment_a_m2`, `dipole_class`, `magnetopause_rp`, `xuv_flux_earth`, `xuv_exposure_index` and `flare_irradiation_index`. They are the inputs of the surface radiation dose (GEN.87) and the habitability index (GEN.83, GEN.89).

### v74 (GEN.88): hydrosphere and ocean chemistry

`planets` and `moons` gain `water_mass_fraction`, `hydrosphere`, `ocean_fraction`, `land_fraction`, `ocean_depth_km`, `ice_shell_km`, `hp_ice_km`, `ocean_class`, `ocean_ph`, `water_activity` and `phosphorus`; `rogue_planets` gains `hp_ice_km` and the `hycean` regime. Rogue oceans are now capped by the liquidus, and their ice lid counts as 1.09 km of ice per km of water with a pressure-lowered base. See `planetgen/physics/hydrosphere.py`.

### v75 (GEN.87): surface radiation dose, UV and the galactic hazard

`planets` and `moons` gain `surface_dose_msv_yr`, `dose_gcr_msv_yr`, `dose_sep_msv_yr`, `dose_ground_msv_yr`, `dose_helio_mult`, `uv_surface_index` and `ozone_loss_flag`; `stars` gains `lethal_event_rate_per_gyr`. Generated with the planet's spin (it needs the air, water and magnetic field) and updated by `refresh_containment` when a system enters or leaves a nebula or remnant. See `planetgen/physics/radiation.py` and docs/design/activity-magnetism-radiation-hydrosphere.md section 4.

### v76 (GEN.89): the habitability score

`planets` and `moons` gain `phi4`, `phi4_pressure`, `phi4_temperature`, `phi4_chemistry`, `phi4_radiation`, `tier_pressure`, `tier_temperature`, `tier_chemistry`, `tier_radiation`, `phi_bio`, `phi_cpx`, `phi_tech`, `l_solv`, `l_chem`, `l_ener`, `l_rad`, `equipment_tier`, `hab_note` and `energy_flux_w_m2`, all NULL for a gas giant and on a row generated before v76. Computed last in a body's generation (after its air, water, dose and light) and again by `refresh_containment` when the dose changes. See `planetgen/physics/habitability_world.py` and docs/design/habitability-index.md section 8.

### v77: `phenomenon_scatter_classes`

`phenomenon_scatter_classes` (`kind`, `subtype`, `placed`, `built`; `subtype` is `''` where a class has none) holds how many scattered phenomena of each class the plan placed and how many a sector's save has built since. The Phenomena table counts and pages the unbuilt ones from it instead of scanning `phenomenon_scatter`, which holds over a hundred million rows in a real galaxy. Written when a phenomena scatter ends, kept by `mark_phenomena_built`, emptied with the scatter.

### v78 (DB.20): object IDs

`uid` is `BINARY(12)` (`BINARY(10)` before v93) and UNIQUE on its own on `star_systems`, `stars`, `planets`, `moons`, `asteroid_belts`, `comets`, the eight phenomenon tables (`black_holes`, `neutron_stars`, `nebulae`, `supernova_remnants`, `quasars`, `rogue_planets`, `interstellar_comets`, `asteroid_fields`) and `facilities`; the old `UNIQUE (star_system_id, uid)` and the 96- and 64-bit hashes are gone. An ID is 96 bits since v93 (80 before): the birth sector's address (shell 19, biased band 18, wedge 19 bits; ring 12, layer 12, slot 16 in the 80-bit form), a serial in that sector (28 bits; the top two say generated, run-time or field-drawn) and a body number (12 bits; 0 for the top-level object, 1 and up for the stars, planets, moons, belts and comets born in a system, one counter across the five kinds). Printed `FE81000A2B-0000005-000`; the layout is `planetgen/galaxy/object_uid.py` (GEN.170). The sector fill gives a sector's systems, then its phenomena, generated serials 0, 1, 2 ... in the order they are inserted (`store._UidIssuer`); a row saved later (an admin-added body, a system saved on its own) takes a run-time serial from the sector's `id_counters` row, or for a body the system's next body number (`store.assign_uids`). `sectors.uid` is unchanged (the designation).

`id_counters (kind, scope, next_value)`, primary key `(kind, scope)`: `kind = 'sector'` with the 5-byte birth sector as `scope` holds the next run-time serial there; `kind = 'system'` with the system's 10-byte ID holds the next body number in it (a system with no row starts past the bodies it holds). A row only moves up, on its own autocommitted connection (`store._reserve_counter`, like `id_blocks`), and `planetgen reset` keeps it, so a number is never given twice. The migration numbers existing rows by row order (a sector's systems, then its phenomena; a system's stars, planets, moons, belts, comets) and gives rows with no sector address run-time serials at the no-sector address (ring 0, layer 0, slot 0). `facilities.uid` is added NULL; GEN.172 fills it.
