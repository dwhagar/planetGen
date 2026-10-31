# planetGen Database Format

This document describes the MySQL database schema defined in
[`src/planetgen/db/schema.sql`](../../src/planetgen/db/schema.sql). It's the
reference for anyone reading, querying, or extending the database — table by
table, every column's meaning and unit, and the conventions that hold the
schema together.

**Status**: the schema is implemented, with both a write and a read path.
[`src/planetgen/db/store.py`](../../src/planetgen/db/store.py) (private — leading
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
  revision per version it is behind, and one older than the baseline (v130)
  is refused with `SchemaTooOldError`, and a database that has the galaxy
  tables but no recorded version (`schema_migrations` emptied or missing)
  is refused with `SchemaVersionMissingError` instead of guessed at from its
  tables. `planetgen.cli.migrate` is its CLI
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

### Units of a column

Every numeric column declares the unit it is stored in, in one list: `db/column_units.py` (DB.76). `TABLES` gives, for each
table, the columns of each unit; `COLUMNS` is the same list by `(table, column)`. A unit is a key of the unit registry
(`physics/unit_registry.py`, GEN.568), which holds the kind of each unit and its factor to the base unit of its kind, or of
`LOCAL_UNITS` for the few the registry does not hold yet. The registry is the only table of units (see
`docs/design/units-and-number-formatting.md`); this page states no factor.

Rules (Boss's cards Q3 and Q4 of `docs/design/unified-unit-quantity.md`):

- **The stored unit is chosen for precision, not forced to SI.** A position keeps its frame's unit (galaxy parsec, sector
  milliparsec, system kilometres or AU) because a double's spacing grows with the number.
- **The suffix stays in the name and equals the declaration.** `SUFFIXES` maps each suffix to its unit, and
  `tests/test_column_units.py` fails when a column's suffix and declaration disagree, when a numeric column has no
  declaration, or when a declaration names no column. A new table cannot ship without declaring its numeric columns in the
  same change.
- **A column with no unit of its own is declared by what it is:** `ratio`, `count`, `index`, `id`, `key`, `seed`, `flag`,
  `unit-vector`, `row-defined` or `unspecified` (a quantity whose code states no unit; it waits for a source, DB.86).
- **Values are not converted on write.** The object classes hold each field in the unit its `FIELD_UNITS` declares
  (`generation/mass_object.py`, for the families declared so far); a value reaches a column in that unit, and a `Quantity` (`physics/quantity.py`) is how it is
  read in any other unit at the edges (the describer, the API, the page).
- The point columns of an object (`point_columns.PointLayout`) read their position and mass units from this list.

### The searchable-field principle

*Historical: v5 dropped every `table_*`/`binary_table_*` column described
in the next paragraph, and display strings are now computed on demand
from the raw columns (`planetgen/web/lib/tabledisplay.py`). The belt columns in
the second paragraph remain.*

Every `*Data`/`*_properties` dict that the generator's old text builder (`to_paragraph_list()`, gone with UX.109) built in the
generator — the exact values that appear in each object's rendered
table — gets its own set of `table_*` columns (one column per dict key),
rather than being reconstructed from the raw scalar columns at query time.
This matters because the raw scalar columns don't always match the display
text one-to-one (units, rounding, phrasing can differ or change across
versions); the `table_*` columns are the "as-published" snapshot, kept
alongside the raw generative fields for full object-graph fidelity.

Asteroid belts have no such dict (their text was prose only) — their searchable columns instead capture the same facts the
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
constraint on `sectors`, turning a concurrent generation race
(two runs reaching one address) into a recoverable `IntegrityError`
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
version: `0130_baseline.py` onward), numbered by the
schema version it produces. `schema.sql` is kept at the current shape, so a
new database is created from it directly and recorded at `SCHEMA_VERSION`; an
older one gets its missing tables from `schema.sql` and then each revision
after its recorded version in order (`alembic_runner.upgrade`, run by
`store.migrate_database`), under a schema lock so a first connection
elsewhere waits instead of stamping the database current half way. The
per-version `_migrate_vN_to_vN+1` functions the entries below mention
(v8 to v60) are gone, and so are the revisions up to v130 (DB.98): v130 is
the baseline, and a database older than that is refused (`SchemaTooOldError`)
and has to be generated again. A database with tables but no recorded version is refused too
(`SchemaVersionMissingError`, DB.97). A database newer
than the code is refused too (`SchemaTooNewError`). Every released schema from v130
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
`star_systems.wikijs_url`/`mediawiki_url`; all three were dropped again in v149, OPS.69). v24
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
`star_systems.wikitext_content`/`markdown_content`; see "Rendered page text" below.

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
`docs/design/galaxy/galaxy-coordinate-system.md`.

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
and `docs/design/galaxy/galaxy-coordinate-system.md`.

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
`docs/design/systems/rogue-planet-surface.md`). `_migrate_v47_to_v48` fills them
for every stored rogue, seeded by its name, and resets
`has_internal_heat` to the computed answer.

**Population and politics (v44).** Filled by `planetgen population`
(or a `galaxy` run, or the optional question
the install and update scripts ask) from what is already stored; see
docs/design/habitability/population-and-politics.md. `species` holds one dominant
species per life world (a planet whose evolutionary milestone is
multicellularity or a technological civilization): a galaxy-unique
`name`, its homeworld (`homeworld_planet_id`, `star_system_id`, both
cascading), `life_class` (the id of its life class, v175), `life_stage`, `build`/`climate`/`size`, and
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
`position_x/y/z_km`. `advance_orbital_phases` (planets and moons, by `kepler.universal_step` on the state) and `advance_comet_orbits` (Kepler or Barker) move it with the
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
them and draws it whole. A band run (`scatter --band --down-to`)
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
saved before this keep their names. See `docs/design/architecture/object-ids.md`.

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
`docs/design/systems/nebula-and-asteroid-field-classes.md`.

**Research-based interstellar rates (v37).** Phenomena are now drawn per
star from Boss's research densities (`program_constants.
PHENOMENON_DENSITY_PC3`, `docs/design/systems/interstellar-object-rates.md`).
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

Twenty-three tables, versioned independently via `control_schema_migrations`
(currently version 21, `store.CONTROL_SCHEMA_VERSION`, mirroring `schema_migrations`'s own shape; v2 added
`login_throttle`, which v8 dropped (SEC.30: the login lockouts live in Redis), v3 `admin_devices`, v4 `admin_totp` and
`admin_recovery_codes`, v5 the work queue's three tables, v6
`generation_stats` and `generation_size`, v7 the work queue's job
tree and pause columns, v9 `galaxy_naming`, v10 `work_job_args`, v11
`version_key_history`, v12 the `version_key` column of `generation_stats`, v14
`generation_stage_runs`, v15 the account columns of `admin_users`, `account_tokens` and `invite_uses`, v16 the five upload tables v17 `generate_recipes` and `generate_recipe_fields` and v18 the epoch, battery and environment columns of `version_key_history`, and v21 the drop of the work queue's four tables (PERF.100); the new tables are added to an older schema by
`CREATE TABLE IF NOT EXISTS` on the next `planetgen.cli.migrate` run, and v7 and v12 add columns, see below):

- **`admin_users`** — one row per account (`username`, `password_hash`,
  `must_change_credentials`). v15 (DB.58) adds `email` (VARCHAR(254), NULL,
  unique), `role` (`user`, `moderator`, `admin` or `owner`, default `admin`,
  so every earlier row is an admin; a moderator is accepted and has no
  permissions until 10.0) and `status` (`active` or `disabled`, default
  `active`). There is exactly one owner: `owner_flag`, a virtual column that
  is 1 for an owner row and NULL otherwise, has a UNIQUE key, so MySQL refuses
  a second owner. The row the install seeds is the owner, and when v15 is
  applied to a control database that has accounts the lowest id becomes the
  owner (`store._promote_control_owner`). The account flows of USR.2 to USR.6
  read these columns; nothing else does yet.
- **`account_tokens`** — v15 (DB.58), the one-time tokens of the account
  flows (USR.28): `account_id` (removed with the account), `purpose`
  (`invite`, `first_password`, `reset`, `email_change`, `owner_transfer_old`
  or `owner_transfer_new`), `token_hash` (the SHA-256 hex digest, unique;
  the token itself is never stored), `created_at`, `expires_at`, `used_at`
  (set once when the token is redeemed) and `meta` (JSON for what a purpose
  needs, such as an invite's use count). Nothing writes it yet.
- **`invite_uses`** — v15 (DB.58), who registered through an invite link:
  `token_id` (the invite's `account_tokens` row) and `account_id`, one row
  per account, and `used_at`. Nothing writes it yet.
- **`upload_runs`, `upload_batches`, `upload_units`, `upload_claims`,
  `upload_id_ranges`** — v16 (DB.62), the bookkeeping for sectors uploaded
  from a remote run (API.11, [api-design-standards.md](../design/operations/api-design-standards.md)
  section 7.3); nothing writes or reads them yet. `upload_runs` is one remote run:
  `owner_id` (an account), `key_id` (the API key, NULL once it is removed),
  `galaxy_database`, `state` (set by API.11), `created_at` and `last_contact_at`.
  `upload_batches` is one row per `(run_id, seq)` with the `sha256` and `bytes`
  of the received batch. `upload_units` is one row per `(run_id, unit_key)` with
  its `state` (`open`, `sealed`, `verified`, `finalized` or `rejected`), the
  `manifest_hash` and the `errors` (JSON). `upload_claims` holds a sector address
  (`galaxy_database`, `shell_index`, `band_index`, `wedge_index`) for one run, so
  a sector can be claimed once. `upload_id_ranges` is the ids a run may use, one
  `first_id` to `last_id` range per table and lease. Every table goes with its run.
- **`admin_sessions`** — web-UI login sessions (`token_hash`, a SHA-256
  digest of the actual cookie value — never the raw token itself;
  `expires_at`, a fixed lifetime set at creation, no sliding renewal).
- **`admin_api_keys`** — API keys for programmatic callers (`key_hash`,
  same "hash only, never the raw key" treatment; `revoked_at` rather than
  a hard delete, so a revoked key's history stays visible; `key_prefix`,
  its first eight characters, tells keys apart without the key;
  `expires_at` makes a key stop working, like a revoked one — control
  schema v13, API.9).
- **`generate_recipes`** — an admin's named, saved Generate forms (`admin_id`,
  `name`; `UNIQUE (admin_id, name)`; deleted with the admin; control schema
  v17, ADM.125, `admin/recipes.py`). **`generate_recipe_fields`** holds a
  recipe's form, one row per `[name, value]` pair in the order the form sent
  them (`recipe_id`, `position`, `field_name`, `field_value`; the first pair
  is the `action`); they go with their recipe.
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
  with `python3 -m planetgen.cli.forget_devices <user>`, and (expired ones)
  when the admin's next device is made.
- **`admin_totp`** (v4, SEC.26) — one row per admin who has set up an
  authenticator app: `secret` (the base32 key itself, as sensitive as a
  password, since a code can only be checked against the key),
  `enabled_at` (NULL until the first code confirms the setup), and
  `last_step` (the newest 30-second step a code was accepted for, so no
  code works twice).
- **`admin_recovery_codes`** (v4) — ten single-use codes per admin with
  two-factor sign-in on: `code_hash` (SHA-256), `used_at`.
- **`work_jobs`**, **`work_tasks`**, **`work_lease`** and **`work_job_args`** (v5 to v10, PERF.8, ADM.12, DB.13) —
  gone in v21 (PERF.100). They were the MySQL job tree of a run, its tasks and the lease on the machine. A run is
  a record in Redis now (`planetgen/queue/runrecord.py`: stages, units, state, log), the pause of the queues is a Redis
  flag (`planetgen/queue/pause.py`), and the time each step of a Generate page job took is kept in
  `step-times.json` in the jobs folder. `update.sh` drops the four tables when it moves an older control schema to v21.

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
- **`version_key_history`** (v11, OPS.13) — the version key each update ran
  under, per galaxy database, newest last: `database_name`, `galaxy_seed`
  (32 hex digits), `version_key` (22), `planetgen_version`,
  `requirements_sha256` (NULL with no `requirements.lock`), `recorded_at`, and
  from v18 (OPS.58) `generator_epoch` (SMALLINT UNSIGNED), `battery_digest`
  (CHAR(64), the generator battery's overall digest, OPS.57), `battery_cases`
  (JSON, `{case: full digest}`), `environment` (JSON: Python, system, machine,
  libc and the numpy, astropy, scipy and scikit-image versions) and
  `last_seen_at` (DATETIME(6)); all five are NULL on a row recorded before v18
  (`last_seen_at` then reads as `recorded_at`). An update adds a row only when
  the seed, the key, the epoch, the battery digest, the lock hash or the
  environment differs from the newest row, and otherwise sets that row's
  `last_seen_at`; the last 10 rows per database are kept
  (`planetgen/galaxy/version_history.py`), so an admin can see which
  release and environment a seed was last known to make the same galaxy
  under and what changed between two rows.
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
when `admin_users` is empty, seeds an `admin` row (the owner) with a random first
password that `planetgen.cli.migrate` prints once (`planetgen.cli.migrate` calls this
automatically, alongside its usual content-schema migration), and every
other function there implements one piece of the login/session/API-key/
audit lifecycle `planetgen/api/auth.py`'s routes expose.

### Booleans and tri-state flags

MySQL has no dedicated boolean type either. Plain booleans are `TINYINT(1)`
`0`/`1` with a `CHECK` constraint. `SystemConfig`'s tri-state flags (`True`/
`False`/`None` in Python — force-on / force-off / random) are nullable
`INTEGER` `0`/`1`/`NULL`.

### Rendered page text

No rendered page text is stored (v29 dropped `star_systems.wikitext_content`/
`markdown_content`). Every value the page shows has its own column, so
`_db.load_star_system` rebuilds the generation object graph and
`planetgen/db/render.py` renders Markdown or plain text on demand — the
same `StarSystem.__str__` generation used to run once before saving, so
the text is identical, except that it now follows later changes (renames,
names made unique, orbit ticks' current wobble and comet positions).
(`planetgen.cli.render_parity`, since removed, compared the two on a database still at v28.)

The wiki links (`sectors.wiki_url`, `star_systems.wikijs_url`/`mediawiki_url`) went in v149 with the wiki upload (OPS.69).

## Planned changes (not built)

Schema changes the open TODO items plan, as of 2026-10-15 (galaxy
schema v92, control schema v14). None of these exists yet; each lands as
a numbered Alembic revision when its item is built, and this section moves
into the table descriptions below. The galaxy schema's revisions are
written one lane at a time, and the second of two branches to reach main
renumbers its revision ([plan/notes.md](../plan/notes.md)). Why and how the
reproducibility pieces fit together is in
[design/generation/reproducible-galaxies.md](../design/generation/reproducible-galaxies.md).

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
| PERF.100 | Done in control schema v21: `work_jobs`, `work_job_args`, `work_tasks` and `work_lease` are dropped now the work queue lives in Redis only (PERF.86); the Redis keys are not part of this schema. |
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
| `center_x_pc`, `center_y_pc`, `center_z_pc` | DOUBLE | nullable | The sector's center, in a galaxy-frame Cartesian coordinate system whose origin is the galactic center (parsecs — see `docs/design/galaxy/galaxy-coordinate-system.md`). NULL together iff this sector has never been placed in a galaxy (`planetgen sector`'s own standalone CLI, or a sector migrated from a pre-v4 database). |
| `galactic_radius_pc` | DOUBLE | nullable | `sqrt(x^2+y^2+z^2)`, persisted (not just derivable) so "sectors within radius R of the core" is a plain indexed range scan — same treatment `star_systems.quadrant` gets. NULL iff the center columns are NULL. |
| `shell_index`, `band_index`, `wedge_index` | INT | nullable | This sector's address on the spherical grid (v91, DB.32, `planetgen/galaxy/spherical.py`; see `docs/design/galaxy/spherical-sector-coordinates.md`): the shell (distance from the galactic centre in 4 pc steps; shell 0 is the core ball), the band (latitude; band 0 holds the plane and a shell's end bands are its polar caps) and the wedge (longitude, counted from the +X meridian). Named `ring_index`, `layer_index` and `ring_slot_index` on the cylindrical grid until v90. Independently nullable from the center/radius columns above (not part of the same CHECK) — a sector could in principle have a hand-authored galaxy position without this particular placement algorithm's own addressing. |
| `generator_epoch`, `run_id` | SMALLINT UNSIGNED, BIGINT UNSIGNED | nullable, `run_id` indexed, foreign key to `generation_runs.id` | Added in v117 (DB.16), replacing the four text columns of v67 (DB.7: `version_key`, `planetgen_version`, `python_version`, `platform`, dropped). `generator_epoch` is the epoch of `util/generator_epoch.py` the code ran under (NULL on a sector saved before v117, when it was not recorded); `run_id` is the `generation_runs` row of the run that saved the sector, which holds the packed version key (`planetgen/galaxy/version_key.py`) and the release, Python and platform strings. A sector saved outside a run points at one shared run row (command `save`) per version of the code. NULL `run_id` on a sector generated before v67. The Generate page and `planetgen galaxy` warn when the running version differs from a sector's (`galaxy/version_check.py`, `store.sector_versions`). |
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
This is what turns a generation race (two runs reaching one address)
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
| `generator_epoch`, `fp_spec` | SMALLINT UNSIGNED | nullable | Added in v131 (OPS.58). The generator epoch (`util/generator_epoch.py`) and the fingerprint format (`util/fingerprint_encoding.FP_SPEC`) of the code that wrote `galaxy_seed`. NULL on a galaxy planned before v131, which `version_check` says nothing about. |
| `battery_digest`, `battery_cases` | CHAR(64) / JSON | nullable | Added in v131 (OPS.58). The generator battery's overall digest (`generation/generator_battery.py`, OPS.57) when `galaxy_seed` was written, and `{case: full digest}` for each case, so a check names the first case that differs. NULL before v131. |
| `pinned_stellar_mass_solar` | DOUBLE | nullable | Added in v113 (DB.81). The stellar mass, solar masses, the plan was told to build the galaxy at (a pinned mass, or the one a pinned size and density hold); NULL when the seed drew it. A loaded shape draws its mass chain and structure again at this mass. |
| `density_ratio_avg` | DOUBLE | nullable | Added in v53 (PERF.11). Systems a fill got against systems expected (`sector_stats`), as a mean over the first 1,000 fills, then a decaying average weighing each new fill 1 in 1,000 (`_db.DENSITY_RATIO_WINDOW`). NULL before any fill. |
| `density_ratio_samples` | BIGINT UNSIGNED | NOT NULL, default 0 | Added in v53 (PERF.11). How many fills went into `density_ratio_avg`. |

### `galaxy_shape_pins`

Added in v113 (DB.81). The shape levers the stored plan was told (`planetgen.galaxy.dials`), one row per lever the plan pinned; a lever with no row was left to the seed (Auto). A plan replaces every row. The galaxy page and `GET /api/galaxy/shape` (`shape.pins`) say "mass pinned at X, the rest drawn" from it.

| Column | Type | Null | Notes |
|---|---|---|---|
| `name` | VARCHAR(40) | PK | The lever's name in the table (`mass`, `compactness`, `disk_scale_length`, ...). |
| `value` | DOUBLE | NOT NULL | The value given, in the table's units (lengths in parsecs). |

### `galaxy_shell`

The galaxy's outline on the spherical grid (v91, replacing `galaxy_layer` and
`galaxy_column`): one row per shell that can hold content, from the core out
to the last shell whose equatorial band could clear the qualification
threshold (`planetgen.galaxy.skeleton.build_shell_extents`). A shell holds
bands `-high_band_index .. high_band_index` (both polar caps when
`high_band_index` is the shell's top band), every wedge. Like the cylinder's
outline it comes from an exact upper bound, so a sector inside it isn't
guaranteed to qualify, but one outside it is guaranteed not to; the exact
per-sector answer is made from the address (`regenerate_sector_cell`, a run); a visit generates nothing (PERF.149). The outline is symmetric about the
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
state vector in this database — absent entirely until that
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
under [`../examples/maintenance/`](../../examples/maintenance/) -- journald
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
| `stars` | TINYINT UNSIGNED | nullable | The forced number of stars (`--stars N`), NULL when left to chance (v134, GEN.267). |
| `max_stars` | TINYINT UNSIGNED | NOT NULL, default 12 | The cap on the star count (`--max-stars`) (v134, GEN.267). |
| `prevalence_habitable_world`, `prevalence_asteroid_belt`, `prevalence_comets`, `prevalence_large_star`, `prevalence_moons`, `prevalence_max_planets`, `prevalence_intelligent_life`, `prevalence_binary_system`, `prevalence_wide_binary`, `prevalence_planets`, `prevalence_multi_star` | DOUBLE | NOT NULL, default 0 | Added in v54 (GEN.52). One per forcing option above: a sector or galaxy run's percentage deviation from that feature's normal chance, 0 for none. |

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
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the galactic pole; see `planetgen/physics/spin.py` and docs/design/orbits/orbital-updates.md section 6. NULL on a row generated before v68. |
| `log_lx_lbol`, `l_xuv_w`, `xuv_saturated`, `flare_n33_per_yr`, `flare_alpha`, `xuv_fluence_j` | DOUBLE (`xuv_saturated` BOOLEAN) | nullable | Added in v72 (GEN.86). Stellar activity: log10 of the coronal X-ray share of the luminosity (saturated at 7.4e-4 for G, K and M stars, decaying after a mass-dependent saturation age; NULL for a white dwarf or compact remnant), the X-ray plus EUV luminosity in W (a hot photosphere's output above 13.6 eV included; a neutron star's thermal or spin-down X-rays; a black hole's disk), whether the corona is still saturated, flares above 1e33 erg a year and their energy slope, and the XUV given off over the star's age, J. See `planetgen/physics/activity.py` and docs/design/habitability/activity-magnetism-radiation-hydrosphere.md section 2. NULL on a row generated before v72. |
| `lethal_event_rate_per_gyr` | DOUBLE | nullable | Added in v75 (GEN.87). Supernovae within 8 pc per Gyr here: 1.5 at the Sun's 8 kpc, rising as exp((8 kpc - R) / 1.45 kpc) to at most 250 times that at the centre. NULL for a star with no galaxy placement and on a row generated before v75. |
| `uv_c_share` | DOUBLE | nullable | Added in v145 (GEN.443). The share of the star's blackbody UV-C (100 to 280 nm) it gives off, log-uniform 0.1 to 1 (the Sun's measured 0.256 is one value of it), drawn once from the star's own stream (`physics/activity.py`); the radiation model's UV-B and UV-A shares follow from it (`physics/stellar_spectrum.py`). NULL for a remnant and for a row saved before v145, which the model reads as the blackbody (1). |
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
| `life_class`, `life_class_variant` | TINYINT UNSIGNED, CHAR(1) | nullable | v175 (GEN.889). The life class that lives on the body (`physics/life_classes.py`): its number and its letter (`1`, `a` is Class 1a; Class 6 has no letter). Set by life-data generation, NULL for none; replaces `life_chemical`. The pace life evolves at is `evolutionary_pace` (below, since v135). |
| `life_class_best`, `life_class_best_score` | VARCHAR(2), TINYINT UNSIGNED | nullable | v175. The class the body supports best and its score, 0 to 9; the classes it supports at all are the rows of `planet_life_class_scores` / `moon_life_class_scores` (one row for each class scoring 1 or more). |
| `life_animal_grade`, `life_complexity` | TINYINT UNSIGNED, DOUBLE | nullable | v175 (GEN.890). The animal grade (0 none, 1 sessile, 2 small and mobile, 3 large and active) and the complexity index of the class on the body (`physics/life_complexity.py`; 1 is an Earth-like biosphere). |
| `flavor_text_count` | INTEGER | NOT NULL, default 0 | |
| `orbital_inclination_deg`, `orbital_ascending_node_deg` | DOUBLE | NOT NULL | Added in v9. Fixed at generation time — together they orient this (circular) orbital plane in 3D. |
| `position_x_km`, `_y_km`, `_z_km` | DOUBLE | NOT NULL | Added in v11. This body's Cartesian position relative to its orbital anchor — the star (or a binary's combined center) for a planet — derived from `distance_km` and the three orbital-motion columns above (`orbits.orbital_position_au`). `planetgen.cli.orbits` moves it by `kepler.universal_step` (v109: there is no stored phase; the phase is derived from this position). |
| `velocity_x_kms`, `_y_kms`, `_z_kms` | DOUBLE | NOT NULL, default 0 | Added in v60 (GEN.121). This body's velocity, km/s, relative to its orbital anchor and on the same axes as `position_x/y/z_km`: the tangent of its circular orbit where the position is, `orbital_speed_kms` long (`orbits.circular_orbital_velocity_au_per_year`). `advance_orbital_phases` moves it with the position. |
| `epoch_unix`, `next_update_due` | DOUBLE | nullable, `next_update_due` indexed | Added in v66 (GEN.106). When the stored position holds and when the body will have moved 0.01 AU (100,000 km for a moon) at its orbital speed; see `orbit_simulation_state`. An edit that changes the orbit sets `next_update_due` back to NULL. |
| `orbital_speed_kms` | DOUBLE | NOT NULL | Added in v11. Constant circular-orbit speed (`orbits.circular_orbital_speed_kms`, `v = 2*pi*r/T`). Only changes if `distance_km`/`period_years` do (e.g. `StarSystem.validate_system` resolving an orbital overlap at generation time), never from phase advancing alone. |
| `min_update_interval_years` | DOUBLE | NOT NULL | Added in v12; not read by the update since v66 (`next_update_due` decides). Not a narrative stat -- a floating-point update guard: the shortest `elapsed_years` worth calling it for, below which the phase delta added is smaller than `orbital_phase_deg`'s own double-precision resolution and so is guaranteed to be a no-op write (`orbits.minimum_update_interval_years`, `period_years * math.ulp(360.0) / 360`). Like `orbital_speed_kms`, only changes if `distance_km`/`period_years` do. |
| `rotation_period_hours` | DOUBLE | NOT NULL | Added in v9. Axial rotation ("day length") — a static descriptive stat; no rotational phase is tracked. |
| `reflex_offset_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v20. This planet's own displacement from its nominal fixed point, from the combined pull of its own moons (`moons.planet_id`) — see `orbits.calculate_reflex_offset`. NULL/0 with no moons. **Not present on `moons`** — a moon never hosts its own moons. |
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the orbit normal (the orbit frame: inclination and ascending node); see `planetgen/physics/spin.py` and docs/design/orbits/orbital-updates.md section 6. NULL on a row generated before v68. A tidally locked body (now planets as well as moons) has tilt 0 and `rotation_period_hours` equal to its orbit. |
| `mantle_redox` | VARCHAR(16) | nullable, CHECK reduced / intermediate / oxidized | Added in v70 (GEN.85). The upper mantle's redox state: reduced at or below the iron-wustite buffer, oxidized from IW+2.5. NULL for a gas giant. See `planetgen/physics/atmosphere.py`. |
| `mantle_delta_iw` | DOUBLE | nullable | Added in v70 (GEN.85). The upper mantle's oxygen fugacity, log units above IW, drawn from mass (Earth about +3.5, the Moon about -1). |
| (the air) | | | Since v155 (DB.71) the partial pressures are not columns of `planets`/`moons` any more: see `planet_gases` and `moon_gases`. `atmosphere` is written from them. |
| `magnetic_moment_a_m2`, `dipole_class`, `magnetopause_rp` | DOUBLE, VARCHAR(16) (CHECK none / weak / earth-like / strong / multipolar), DOUBLE | nullable | Added in v72 (GEN.86). The magnetic dipole moment, A m2 (0 with no dynamo), its class (multipolar when slow rotation leaves under half of it in the dipole), and the magnetopause standoff in body radii against the star's wind (NULL with no field). See `planetgen/physics/magnetism.py` and docs/design/habitability/activity-magnetism-radiation-hydrosphere.md section 3. NULL on a row generated before v72. |
| `xuv_flux_earth`, `xuv_exposure_index`, `flare_irradiation_index` | DOUBLE | nullable | Added in v72 (GEN.86). The star's XUV flux here in units of Earth's today, its lifetime XUV here against a 1 solar mass star's at 1 AU over 5 Gyr (2 or less Earth-like, over 30 extreme), and its flares above 1e33 erg a year over the distance squared in AU (Earth 0.003). A moon takes its planet's distance. NULL on a row generated before v72. |
| `water_mass_fraction`, `hydrosphere`, `ocean_fraction`, `land_fraction` | DOUBLE, VARCHAR(20) (CHECK dry / vapour / ice / ice-covered ocean / surface ocean / brine seep / hycean), DOUBLE, DOUBLE | nullable | Added in v74 (GEN.88). The body's water share of its mass (drawn by class), where its surface water is (gone, in the air, frozen, under ice, open, or under a hydrogen sky), and the shares of the surface the water (liquid or frozen) and the land cover. NULL for a gas giant. See `planetgen/physics/hydrosphere.py` and docs/design/habitability/activity-magnetism-radiation-hydrosphere.md section 5. NULL on a row generated before v74. |
| `ocean_depth_km`, `ice_shell_km`, `hp_ice_km` | DOUBLE | nullable | Added in v74 (GEN.88). The mean liquid depth where there is water (NULL without liquid), the ice lid, and the high-pressure ice (VI or VII) under the deepest liquid the liquidus allows (NULL when there is none). |
| `water_base_pressure_pa` | DOUBLE | nullable | Added in v147 (GEN.461). The pressure at the bottom of the water layer over the covered area, Pa (liquid, ice or both), from the compressed column of `physics/water_layer.py`. NULL for a gas giant, a body with no water layer and a row made before v147, whose `ocean_depth_km`, `ice_shell_km` and `hp_ice_km` are water-equivalent km. From v147 those three are the real, compressed kilometres. |
| `mantle_eta0_pa_s`, `mantle_activation_kj_mol` | DOUBLE | nullable | Added in v153 (GEN.699). The two draws of a rocky body's mantle rheology, `eta = eta0 exp(Ea / R T)`: the viscosity prefactor in Pa s and the activation energy in kJ/mol (uniform over the ranges of Guimond et al. 2022). NULL for a gas giant and a row made before v152. |
| `internal_heat_flux_w_m2` | DOUBLE | nullable | Added in v153 (GEN.699). The heat flowing out of the interior, radiogenic plus primordial (the cooling of the mantle, `physics/interior_heat.py`), over the body's surface, W/m2, at the star's age when the body was last generated. NULL for a gas giant and a row made before v152. |
| `lid_regime` | VARCHAR(10) | nullable | Added in v154 (GEN.721). The tectonic regime of a rocky body: `mobile` (plate tectonics) or `stagnant` (one rigid plate), drawn once after the hydrosphere. NULL for a gas giant and for a row made before v154. |
| `ocean_class`, `ocean_ph`, `water_activity`, `phosphorus` | VARCHAR(16) (CHECK ice-sealed / chloride brine / acid sulfate / soda / neutral), DOUBLE, DOUBLE, VARCHAR(8) (CHECK high / limited / starved) | nullable | Added in v74 (GEN.88). The ocean's chemistry class by rule (section 5.4), its pH and water activity drawn in the class's range, and its phosphorus supply. NULL without liquid water. |
| `surface_dose_msv_yr`, `dose_gcr_msv_yr`, `dose_sep_msv_yr`, `dose_ground_msv_yr`, `dose_helio_mult` | DOUBLE | nullable | Added in v75 (GEN.87). The yearly dose at the surface, mSv, and its parts: `surface_dose_msv_yr` = `dose_gcr_msv_yr` x `dose_helio_mult` + `dose_sep_msv_yr` + `dose_ground_msv_yr`. Cosmic rays through the air (a log-log curve through the Moon, Mars and a 100 g/cm2 point, then 164 g/cm2 attenuation) cut to 0.6 at most by a dipole; stellar particles (20 mSv/yr x sqrt(N33 / N33 of the Sun) / d_AU^2, stopped by the air and steered by the dipole); crust gamma and radon (0.48 mSv/yr x abundance x radiogenic heat at the star's age against 4.5 Gyr, plus 0.5 mSv/yr x abundance with air and land). `dose_helio_mult` is 1 to 2.5: 1 in open space, raised while a nebula or remnant presses on the star's heliosphere (`refresh_containment` and an edit update it and the total). NULL for a gas giant. |
| `dose_belt_msv_yr` | DOUBLE | nullable | Added in v187 (GEN.983). The dose, mSv/yr, a moon takes at its orbit from the trapped-particle belts of its parent planet's field, part of `surface_dose_msv_yr`: 0 for a planet or a moon outside the belts, NULL for a giant and on a row made before v187. |
| `uv_surface_index`, `ozone_loss_flag` | DOUBLE, BOOLEAN | nullable | Added in v75 (GEN.87). DNA-weighted surface UV against Earth's 1 (the star's 200 to 300 nm power against the Sun's, floored at 0.1, times its flux, times the ozone shielding (`ozone_column_du` / 300 DU)^-1.6 up to 1,000 when there is none), and whether ozone is lost: the layer exists but the star's `lethal_event_rate_per_gyr` is over 10 (one per 100 Myr) or particles give 100 mSv/yr at the ozone layer. NULL for a gas giant, and `uv_surface_index` NULL without a star temperature. |
| `tau_haze`, `tau_dust`, `tau_cloud`, `cloud_cover`, `aerosol_draw` | DOUBLE | nullable | Added in v136 (GEN.414, DB.68). The haze, dust and in-cloud optical depths at 550 nm, the share of the sky under cloud (0 to 1) and the number in [0, 1) that picks the aerosol optics inside their brackets (`physics/aerosols.py`); drawn from the body's own stream after its water, NULL for a body with no air and for a row saved before v136. |
| `c_k0_ppm`, `c_th0_ppm`, `c_u0_ppm`, `al26_ratio`, `urey_ratio` | DOUBLE | nullable | Added in v141 (GEN.459). The radiogenic content of a rocky body: potassium, thorium and uranium in the bulk silicate at formation (ppm by mass of the element), the 26Al/27Al atom ratio at formation, and the Urey ratio today (an output of the thermal history since v146, not a draw). Drawn from the host star's abundances (Th and U against Mg, potassium by a depletion against the CI K/U); the interior heat flux and the ground dose both read them. NULL for a gas giant and on a row generated before v136. |
| `phi4`, `phi4_pressure`, `phi4_temperature`, `phi4_air`, `phi4_radiation` | DOUBLE | nullable | Added in v76 (GEN.89); `phi4_chemistry` was renamed `phi4_air` in v87. The PHI-4 display score for a human visitor: the four domains each 0 to 1 (1 is Blue) and `phi4`, their geometric mean (0 when any domain is 0). NULL for a gas giant and on a row generated before v76. |
| `tier_pressure`, `tier_temperature`, `tier_air`, `tier_radiation` | TINYINT | nullable | Added in v76 (GEN.89); `tier_chemistry` was renamed `tier_air` in v87. Each domain's colour: 0 Blue, 1 Green, 2 Yellow, 3 Red. |
| `eq_body`, `eq_body_kind`, `eq_body_rec`, `eq_body_rec_kind`, `eq_rad_suit`, `eq_rad_suit_rec`, `eq_mobility`, `eq_dust`, `eq_gas_filter`, `eq_co2_scrub`, `eq_supply`, `eq_supply_rec`, `eq_shelter`, `eq_shelter_rec` | TINYINT | nullable | Added in v136 (GEN.655, DB.93), replacing `equipment_tier` (v76). The survival equipment a human needs, as levels: `eq_body` is the required Body Protection (0 none, 1 clothing, 2 thermal suit, 3 sealed suit, 4 pressure suit, 5 armoured suit, 6 vehicle only) and `eq_body_kind` what set it (0 none, 1 cold, 2 heat, 3 pressure, 4 chemical, 5 diving, 6 UV, 7 gravity); `eq_body_rec` and `eq_body_rec_kind` the recommended level, above the required one only (else NULL); `eq_rad_suit` and `eq_rad_suit_rec` the radiation suit step (1 shielded garment, 2 storm vest, 3 hardened suit); `eq_mobility` the powered frame (1 recommended, 2 required); `eq_dust` the dust filter or seals; `eq_gas_filter` the gas filter grade 0 to 3; `eq_co2_scrub` the CO2 scrubber; `eq_supply` and `eq_supply_rec` the Air Supply 0 to 4; `eq_shelter` and `eq_shelter_rec` the radiation shelter (1 shielded shelter, 2 hardened habitat). Where an item has a recommended and a required level, the recommended is stored in the `_rec` column (a 1 in a field with only one encoding is recommended, 2 required). NULL for a gas giant, a body saved without a score and a row not yet re-scored (`planetgen rescore-habitability`). The name, the short code and the lists are derived (`physics/equipment.py`), never stored. |
| `shelter_column_g_cm2`, `outside_hours` | DOUBLE | nullable | Added in v136 (GEN.655). The shielding in g/cm2 a shelter needs to bring the year's dose to the limit, and the hours outside a year the listed suit allows (NULL when unlimited or no shelter is needed). |
| `band1_mean_k` to `band7_mean_k` | DOUBLE | nullable | Added in v196 (GEN.1108), replacing `temp_polar_k`, `temp_north_k`, `temp_equator_k`, `temp_south_k`, `temp_day_side_k`, `temp_terminator_k` and `temp_night_side_k` (v137). The annual mean temperature in kelvin of the seven bands of latitude of `physics/latitude_bands.py`, `band1` the north pole to `band7` the south pole, from the annual insolation of the stored axial tilt and the air column's heat transport (`physics/regional_temperature.py`). For a world that turns once an orbit the seven are rings from the point under the star (`band1`) to the point away from it (`band7`). NULL for a gas giant, an unscored body, a body with no rotation period and a row made before v195. |
| `band1_swing_k` to `band7_swing_k` | DOUBLE | nullable | Added in v196 (GEN.1108), replacing `temp_swing_k` and the swings of v182 by place. The day and night swing of the band, peak to peak, in kelvin of difference (the equator's cut by the cosine of the latitude); 0 on every ring of a world that turns once an orbit. |
| `band1_day_hours` to `band7_day_hours`, `band1_night_hours` to `band7_night_hours` | DOUBLE | nullable | Added in v196 (GEN.1108). The hours of day and of night in a solar day at the band's latitudes (the share of one turn the upper limb of the star is up, `daylight.day_share_at`, averaged over the band). NULL on a ring and where the star never sets or never rises. |
| `band1_longest_day_days` to `band7_longest_day_days`, `band1_longest_night_days` to `band7_longest_night_days` | DOUBLE | nullable | Added in v196 (GEN.1108). The longest stretch of the year, in days, in which the star stays up (or down) at the poleward edge of the band; 0 where it does not happen. NULL on a ring. |
| `band1_state` to `band7_state` | VARCHAR(24) | nullable | Added in v196 (GEN.1108). The state of the band's day and night: `day and night`, `polar day and night` or `non-rotating`; on a world that turns once an orbit `permanent day side`, `terminator band` or `permanent night side`. |
| `ice_co2_mass_fraction`, `ice_ch4_mass_fraction`, `ice_nh3_mass_fraction`, `ice_n2_mass_fraction`, `ice_so2_mass_fraction`, `ice_co_mass_fraction`, `ice_hcn_mass_fraction`, `ice_other_mass_fraction` | DOUBLE | nullable | Added in v183 (GEN.834). The share of the whole body's mass that each non-water ice holds (carbon dioxide, methane, ammonia, nitrogen, sulfur dioxide, carbon monoxide, hydrogen cyanide; `ice_other_mass_fraction` is ethane, methanol, formaldehyde and hydrogen sulfide), adding up to `ice_mass_fraction`. 0 for a body with no ice; NULL for a row made before v183. |
| `outside_hours_day`, `outside_hours_limit` | DOUBLE, VARCHAR(16) | nullable | Added in v187 (GEN.768). The hours of one local day a person can be outside in the equipment the survival rules list for the world, and the limit that sets them: `daylight`, `dose`, `heat`, `cold`, `air`, `uv` or `vehicle` (`physics/outside_hours.py`). A tidally locked world is judged in its terminator belt over a 24-hour shift. NULL for a gas giant, an unscored body, a world whose daylight is not stored and a row made before v187. |
| `day_night_state` | VARCHAR(24) | nullable | Added in v193 (GEN.1097). The state of the body's day and night (`physics/day_night.py`): `day and night`, `tidally locked`, `spin-orbit resonance`, `non-rotating`, `polar day and night` or `no solid surface` (a giant). NULL for a body whose spin is not stored and a row made before v188. |
| `night_hours`, `days_per_year` | DOUBLE | nullable | Added in v193 (GEN.1097). The hours of night in one solar day at the equator (the solar day less `daylight_hours`) and the number of solar days in one orbit round the star (365.25 for Earth, 0.5 for Mercury). NULL for a world that turns once an orbit and for a row made before v188. |
| `roundness_ratio`, `yield_strength_pa`, `binding_energy_j` | DOUBLE | nullable | Added in v197 (GEN.1159, GEN.1178). The mean radius over the round radius of the body (1 or more is round, less is not; `physics/roundness.py`), the whole-body strength of its material mix in Pa, drawn from the range of each material by a quantile from a stream of the body's own seed, and the energy that assembling it releases in J (`physics/binding_energy.py`: two uniform shells for a rocky body with a core, the n = 1 polytrope for a gas giant). `physics/body_strength.py` sets them. A gas giant keeps only `binding_energy_j`; the other two are NULL for it. NULL for a row made before v195; nothing is backfilled. |
| `temp_mean_swing_k` | DOUBLE | nullable | Added in v182 (GEN.848), kept in v196 (GEN.1108). The day and night swing of the whole surface, peak to peak, in kelvin of difference: the area-weighted mean of the swings of the seven bands (`band1_swing_k` to `band7_swing_k`). NULL for a gas giant, an unscored body and a row made before v195. |
| `l_solv`, `l_chem`, `l_ener`, `l_rad` | DOUBLE | nullable | Added in v76 (GEN.89). The four likelihoods behind microbial life: `l_solv` solvent, `l_chem` chemistry, `l_ener` energy and `l_rad` radiation, each 0 to 1. (v76 also stored `phi_bio`, `phi_cpx` and `phi_tech`; v136 dropped them for the four life scores below.) |
| `life_microbial`, `life_flora`, `life_fauna`, `life_sentient` | TINYINT UNSIGNED | nullable | Added in v150 (GEN.649). The four life scores, 0 (none) to 9 (abundant), from the evolutionary path the body reached and its surface conditions ([habitability-index.md](../design/habitability/habitability-index.md) section 11). NULL for a gas giant or a body with no score. |
| `energy_flux_w_m2` | DOUBLE | nullable | Added in v76 (GEN.89). The light or chemical power life can use: 6% of the starlight, or 1% of the body's own heat flow, whichever is more (the heat only, under an ice lid). The one input the score keeps, so the dose refresh can re-score from the row. |
| `flux_variation` | DOUBLE | nullable | Added in v87 (DB.26, GEN.148). The relative swing of the stellar flux over the orbit or pair cycle; 0 is steady. NULL until regenerated. |
| `induced_moment_a_m2` | DOUBLE | nullable | Added in v87 (DB.26, GEN.178), placed after `dipole_class`. The induced magnetic moment in A m2; NULL when not computed or none. Since GEN.981 a moon's holds its core's answer as well as its ocean's, and a rocky planet's holds the answer of its core and ocean to its star's field (a giant's stays NULL); a row made before that keeps the ocean-only value, and a planet's is NULL until its galaxy is regenerated. |
| `phi4_water`, `phi4_suitability` | DOUBLE | nullable | Added in v87 (DB.26, GEN.207, PHI-6), beside `phi4_air`. The water score and the suitability (the lowest of water, food and shelter). NULL until regenerated. |
| `tier_water`, `tier_suitability` | TINYINT | nullable | Added in v87 (DB.26, GEN.207). The colour of each, 0 Blue to 3 Red like the other tiers. |
| `metal_mass_fraction`, `rock_mass_fraction`, `ice_mass_fraction`, `envelope_mass_fraction` | DOUBLE | nullable | Added in v144 (DB.69, GEN.428). The body's bulk composition as shares of the whole mass: the iron-rich core (light elements included), the rock, every ice other than water, and the hydrogen and helium envelope. With `water_mass_fraction` they add up to 1. NULL on a row saved before v135. Read through `planetgen.physics.bulk_composition.bulk_fractions`. |
| `rock_feo_fraction`, `rock_al2o3_fraction`, `rock_mg_si_ratio` | DOUBLE | nullable | Added in v144. FeO and Al2O3 as mass fractions of the rock, and the rock's molar Mg/Si. NULL for a giant's rock-less composition and on a row saved before v135. |
| `envelope_helium_fraction`, `envelope_z_fraction` | DOUBLE | nullable | Added in v144. Y and the heavy elements as mass fractions of the envelope; NULL for a body without one. |
| `porosity` | DOUBLE | nullable | Added in v144. The void share of a small icy body's volume; NULL for a body that has none. |
| `carbon_rock`, `composition_origin` | BOOLEAN, VARCHAR(8) (CHECK drawn / inverse) | nullable | Added in v144. Whether the rock is silicon carbide and graphite (a star of C/O 0.8 or more, or class C), and whether the composition was drawn first (`drawn`) or worked back from a given mass and radius (`inverse`). |
| `interior_kind` | VARCHAR(16) (CHECK gas_giant / ice_giant / brown_dwarf) | nullable | Added in v154 (DB.89, GEN.642). The kind of interior a giant was drawn with. NULL for a body that has none (every rocky world, and a giant saved before v150). |
| `core_mass_earth`, `dilute_core_mass_earth` | DOUBLE | nullable | Added in v154. The compact core's (an ice giant's rock core's) mass and the dilute core's mass, in Earth masses; 0 with none. |
| `metallization_pressure_gpa`, `helium_depletion` | DOUBLE | nullable | Added in v154. The drawn pressure where hydrogen turns metallic, in GPa, and the molecular envelope's helium share against the natal share (1 with no helium rain). |
| `has_helium_shell` | BOOLEAN | nullable | Added in v154. Whether a helium-rich shell sits above the metallic envelope. Read with the layers through `planetgen.physics.giant_interior.interior_of`. |
| `ground_type`, `ground_grain` | VARCHAR(16), VARCHAR(8) | nullable | Added in v159 (GEN.441). The ground type with the largest areal share (`ocean`, `ice` or a rock type: `ultramafic`, `anorthositic`, `basaltic`, `felsic`, `intermediate`, `clay`, `oxidised`, `metal`, `carbon`, `evaporite`, `sulfur`) and the grain state (`slab`, `crushed`, `powder`). NULL for a giant and on a row saved before v152. |
| `surface_liquid`, `surface_liquid_fraction`, `surface_liquid_depth_km` | VARCHAR(32), DOUBLE, DOUBLE | nullable | Added in v188 (GEN.1035). The liquid the surface holds (an id of `physics/surface_liquids.py`: `water`, a hydrocarbon, a brine, a lava...), the share of the surface it covers and its mean depth over that share, km (NULL for a lava surface, whose depth is not known). NULL for a world with no liquid at its surface and for a row made before v188. |
| `ground_ocean_fraction`, `ground_ice_fraction`, `ground_ultramafic_fraction`, `ground_anorthositic_fraction`, `ground_basaltic_fraction`, `ground_felsic_fraction`, `ground_clay_fraction`, `ground_oxidised_fraction`, `ground_metal_fraction`, `ground_carbon_fraction`, `ground_evaporite_fraction`, `ground_sulfur_fraction`, `ground_intermediate_fraction` | DOUBLE | nullable | Added in v152 (`ground_intermediate_fraction` in v190, GEN.474). The areal share of the surface that is each ground type; the thirteen add up to 1. The minerals, ices and ocean they stand for are rows of the materials table ([surface-composition.md](../design/composition/surface-composition.md)). |
| `moon_origin` | VARCHAR(8) | nullable | Added in v180 (DB.70). How a moon came to orbit its parent: `disk`, `impact`, `captured`, `ring` or `binary` ([small-bodies-composition.md](../design/composition/small-bodies-composition.md)). NULL for a planet and for a moon saved before v179. |
| `moon_disk_lifetime_myr`, `moon_disk_alpha_over_f`, `moon_disk_f_inf`, `moon_disk_ice_line_km` | DOUBLE | nullable | Added in v180 (DB.70). A giant planet's circumplanetary disk, where its moons formed: the disk's lifetime in Myr, its alpha over f, the asymptotic ice fraction of its ramp and the distance of its ice line from the giant in km. NULL for a body with no disk. |

| `bond_albedo`, `greenhouse_factor` | DOUBLE | nullable | Added in v142 (GEN.412). The two draws the surface temperature is made of: the Bond albedo A and the greenhouse quantity G, with T_e = ((1 - A) S / (4 sigma))^(1/4) and T_s = T_e (1 + G)^(1/4) (S the starlight at the orbit; a moon's is its planet's; never below 2.725 K). G is 0 for a body with no air. |
| `albedo_basis` | TINYINT UNSIGNED | nullable | Added in v142 (GEN.412). Where the two came from: 0 recorded at generation, 1 the exact inverse of a row saved without them (from its stored scale height), 2 a seeded draw from the class range where the inverse had no solution. NULL for a row with no temperature or no starlight to work from. |
| `tidal_flux_w_m2` | DOUBLE | nullable | Added in v152 (GEN.700). The eccentricity-tide heat flow in W/m2 of an icy moon with a frozen surface (the same flow that thins its ice lid), which the chemical energy of its ocean's life counts: `energy_flux_w_m2` is 1% of the moon's own heat plus this. NULL for a planet, a moon with a liquid or dry surface, a rocky moon, and a body saved before v144. |

### `moon_reflection_spectrum`

Moon counterpart of `planet_reflection_spectrum` — identical shape,
`moon_id` FK -> `moons.id` (`ON DELETE CASCADE`) instead of `planet_id`.

### `planet_interior_layers`

Added in v154 (DB.89, GEN.642). One row per layer of a giant planet's interior, outside in. `planet_id` FK -> `planets.id` (`ON DELETE CASCADE`), `layer_index` (SMALLINT UNSIGNED, 0 at the top; unique with `planet_id`), `kind` (CHECK molecular_envelope / metallic_envelope / helium_gradient / helium_shell / dilute_core / compact_core / ice_mantle / rock_core), `state` (CHECK molecular_fluid / metallic_fluid / ionic_fluid / superionic / molten), `top_pressure_gpa` and `bottom_pressure_gpa` (GPa), `top_radius_fraction` and `bottom_radius_fraction` (a share of the body's radius), `mass_earth` (the layer's mass) and the mass of each material in it, in Earth masses (`hydrogen`, `helium`, `iron`, `rock`, `water_liquid`, `ice_x`, `nh3_ice`, `ch4_ice`; 0 where the layer has none). The layers add up to the body's mass, and the body's bulk composition columns are their totals ([giant-interiors.md](../design/composition/giant-interiors.md)). A body with no interior has no rows.

### `moon_interior_layers`

Moon counterpart of `planet_interior_layers` — identical shape, `moon_id` FK -> `moons.id` (`ON DELETE CASCADE`) instead of `planet_id`.

### `planet_gases` / `moon_gases`

Added in v155 (DB.71, GEN.424). The air of a planet or moon: one row per gas above zero, replacing the ten `p_<gas>_kpa` columns. A class's balance gas (helium, sodium, ammonia) is a row like the others, so the surface pressure is the sum of the rows. No row means no air (or a body saved before GEN.85).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `planet_id` / `moon_id` | BIGINT UNSIGNED | NOT NULL | FK -> `planets.id` / `moons.id` (`ON DELETE CASCADE`). |
| `gas` | VARCHAR(16) | NOT NULL | The formula in the gas table (`physics/data/gases.json`): `O2`, `CO2`, `N2`, `He`, `Na`, `NH3`. UNIQUE with the owner. |
| `partial_pa` | DOUBLE | NOT NULL, CHECK above 0 | The partial pressure, Pa. |

### `asteroid_belts`

One row per asteroid belt. Belts have no properties-dict data table
(their text was prose only), so their searchable
columns capture the facts the prose always states instead.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `star_system_id` | INTEGER | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | |
| `star_id` | INTEGER | FK -> `stars.id`, `ON DELETE SET NULL`, nullable | Added in v15. Same "which specific star this orbits" semantics as `planets.star_id` above — NULL for a single star's or a `'close'` binary's belt, set to the owning star for a `'wide'` binary's. |
| `orbital_index` | INTEGER | NOT NULL | |
| `distance_km` | DOUBLE | NOT NULL | Average distance from the star. |
| `lower_limit_km`, `upper_limit_km` | DOUBLE | NOT NULL | The belt's inner/outer boundary. |
| `density_g_cm3`, `ice_mass_fraction` | DOUBLE | nullable | Added in v197 (GEN.1159). The mean bulk density of the members in g/cm3 and the mean share of their mass that is ice of any kind, from the drawn makeup (`physics/small_body_composition.py`). NULL for a row made before v195. The roundness ratio, strength and binding energy belong to a member, not to the population. |
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
`docs/design/systems/comet-orbital-realism.md`. Contrast `interstellar_comets`
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
| `density_g_cm3`, `ice_mass_fraction`, `roundness_ratio`, `yield_strength_pa`, `binding_energy_j` | DOUBLE | nullable | Added in v197 (GEN.1159, GEN.1178). The bulk density of the nucleus in g/cm3, the share of its mass that is ice of any kind, and the same three values a planet keeps (`physics/body_strength.py`). NULL for a row made before v195. |
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
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the orbit normal (YORP: near 10 or 170 degrees); see `planetgen/physics/spin.py` and docs/design/orbits/orbital-updates.md section 6. NULL on a row generated before v68. |

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
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the galactic pole; see `planetgen/physics/spin.py` and docs/design/orbits/orbital-updates.md section 6. NULL on a row generated before v68. NULL for an anchored remnant: its spin is on its `stars` row. |

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
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v18 (`planetgen phenomenon --sector-id`): this nebula's own galaxy-frame center, in the same Cartesian space `sectors.center_x/y/z_pc` uses -- a sphere (this + `radius_ly`), not a sector-relative offset, since a nebula (up to 200 ly across) is frequently far larger than one sector (default edge 4 pc, ~13 ly) and may overlap several. NULL together: never placed in the galaxy (still the default -- `--sector-id` is optional). See `queryDb.phenomena_near_sector`/`galaxy_placed_phenomena` for how this is read back, and `docs/design/galaxy/galaxy-coordinate-system.md` for the coordinate system itself. |
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
| `phi4`, `phi4_pressure`, `phi4_temperature`, `phi4_air`, `phi4_water`, `phi4_radiation`, `phi4_suitability` | DOUBLE | nullable | Added in v153 (GEN.701). The PHI-6 score of the rogue planet's surface and its six domains, as on `planets`. NULL for a gas giant or brown dwarf. |
| `tier_pressure`, `tier_temperature`, `tier_air`, `tier_water`, `tier_radiation`, `tier_suitability` | TINYINT | nullable | Added in v153. Each domain's colour (0 Blue to 3 Red). |
| `life_microbial`, `life_flora`, `life_fauna`, `life_sentient` | TINYINT UNSIGNED | nullable | Added in v153. The four life scores. A rogue planet has no life draw, so a scored one holds 0 in each. |
| `l_solv`, `l_chem`, `l_ener`, `l_rad` | DOUBLE | nullable | Added in v153. The four likelihoods behind microbial life, as on `planets`. |
| `eq_body`, `eq_body_kind`, `eq_body_rec`, `eq_body_rec_kind`, `eq_rad_suit`, `eq_rad_suit_rec`, `eq_mobility`, `eq_dust`, `eq_gas_filter`, `eq_co2_scrub`, `eq_supply`, `eq_supply_rec`, `eq_shelter`, `eq_shelter_rec` | TINYINT | nullable | Added in v153. The survival equipment levels a visitor needs, as on `planets`. |
| `shelter_column_g_cm2` | DOUBLE | nullable | Added in v153. The shelter column, as on `planets` (NULL when no shelter applies). |
| `energy_flux_w_m2` | DOUBLE | nullable | Added in v153. The chemical power life could use: 1% of `internal_heat_flux_w_m2` (no starlight). |
| `galactic_orbital_period_gy`, `_phase_deg`, `_min_update_interval_years` | DOUBLE | NOT NULL | Added in v17. Always populated (a rogue planet is always standalone). Its speed column was dropped in v103 (GEN.336): the speed is the length of `velocity_x/y/z_kms` (a star and a close pair read their system's). |
| `center_x_pc`, `center_y_pc`, `center_z_pc`, `galactic_radius_pc` | DOUBLE | nullable, NULL together | Added in v28: this rogue planet's own galaxy-frame center, the same shape `nebulae` uses. NULL together: never placed in the galaxy. |
| `quadrant` | VARCHAR(4) | nullable, CHECK IN ('I'..'VIII') | Added in v41. The sector octant its center sits in, as `star_systems.quadrant`. NULL when unplaced. |
| `inside_nebula_id`, `inside_remnant_id` | BIGINT UNSIGNED | nullable, `ON DELETE SET NULL` | Added in v39. The innermost cloud holding it, as on `star_systems`. |
| `created_at`, `modified_at` | TIMESTAMP / TIMESTAMP(3) | NOT NULL | Added in v27. Row timestamps, `modified_at` indexed. |
| `rotation_period_hours` | DOUBLE | nullable | Added in v68 (GEN.104). A planet's day for its type. The spin vector is the axis times 2 pi over this period. |
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the galactic pole; see `planetgen/physics/spin.py` and docs/design/orbits/orbital-updates.md section 6. NULL on a row generated before v68. |
| `log_lx_lbol`, `l_xuv_w`, `xuv_saturated`, `flare_n33_per_yr`, `flare_alpha`, `xuv_fluence_j` | DOUBLE (`xuv_saturated` BOOLEAN) | nullable | Added in v72 (GEN.86). Stellar activity: log10 of the coronal X-ray share of the luminosity (saturated at 7.4e-4 for G, K and M stars, decaying after a mass-dependent saturation age; NULL for a white dwarf or compact remnant), the X-ray plus EUV luminosity in W (a hot photosphere's output above 13.6 eV included; a neutron star's thermal or spin-down X-rays; a black hole's disk), whether the corona is still saturated, flares above 1e33 erg a year and their energy slope, and the XUV given off over the star's age, J. See `planetgen/physics/activity.py` and docs/design/habitability/activity-magnetism-radiation-hydrosphere.md section 2. NULL on a row generated before v72. |

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
| `spin_axis_x`, `_y`, `_z`, `axial_tilt_deg` | DOUBLE | nullable | Added in v68 (GEN.104). The spin axis as a unit vector and its angle, degrees, from the galactic pole (YORP: near 10 or 170 degrees); see `planetgen/physics/spin.py` and docs/design/orbits/orbital-updates.md section 6. NULL on a row generated before v68. |

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
| `density_g_cm3`, `ice_mass_fraction` | DOUBLE | nullable | Added in v197 (GEN.1159). As on `asteroid_belts`: the mean bulk density of the members and the mean ice share of their mass. NULL for a row made before v195. |
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

Added in v89 (DB.28, POP.10); design in `docs/design/habitability/population-and-politics.md`, "Data model" of the facility section. A type refines a kind.

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
| `rotation_period_hours` | DOUBLE | nullable | Added in v185 (GEN.846). The star's spin period in hours, drawn at scatter time by the same cascade a built star uses (`physics/spin.py`: a white dwarf about a day, a giant slow, a cool dwarf by gyrochronology, a hot star by a drawn speed under breakup) on a stream of its own seed; the fill that builds its system keeps it. NULL on rows made before v184, nothing is backfilled. |
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
| `habitable_worlds` | INT | nullable | The planets and moons in its systems that are potentially habitable (overall PHI-6 Blue or Green); NULL before a fill (v135, MAP.143). The Galaxy Map's Habitable worlds color comes from it. |
| `total_mass_kg` | DOUBLE | nullable | The summed mass (kg) of the bodies in it (v162, DB.108); 0 for a generated sector with none, NULL before a fill or when a body has no stored mass. See "v168 (DB.108)". |
| `filled_at` | TIMESTAMP(3) | nullable | When the fill saved it. |
| `updated_at` | TIMESTAMP | NOT NULL | |

### `sector_paths`

Added in v62 (GEN.123). One row per star system, rogue planet or
interstellar comet: the path it takes from where it is now to where it
leaves its sector (`physics/sector_path.py`), recomputed by
`sector_paths.compute_sector_paths` for every sector holding one of them
at each orbit update, and for the sectors a `galaxy` run created and the
sectors around them as the last step of the run (`planetgen settle`).
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
| `generator_epoch` | SMALLINT UNSIGNED | nullable | Added in v131 (OPS.58). The generator epoch the run's code generated under, so a galaxy says which epoch each pass of it (a scatter, a fill, a redo) was made at. NULL on a run from before v131. |

The command line is in `generation_run_arguments` (v63, DB.13), which replaced a JSON text column of `generation_runs`: `run_id` (FK -> `generation_runs.id`, `ON DELETE CASCADE`), `position` (the argument's place, from 0) and `value` (VARCHAR(1024)), primary key `(run_id, position)`. It leaves out the `--mysql-*` and `--debug` options (`generate._run_argv`).

### `scatter_epochs`

Added in v161 (GEN.764). The generator epoch each kind of scatter was drawn under. A scatter is galaxy-wide and its rows carry no epoch of their own, so the draw that finishes a kind writes one row here and a clear of that kind deletes it. `kind` (VARCHAR(16), primary key; CHECK `mass`, `luminosity`, `band`, `phenomena` or `nebula-field`), `generator_epoch` (SMALLINT UNSIGNED, nullable: unknown), `run_id` (FK -> `generation_runs.id`, `ON DELETE SET NULL`, nullable) and `recorded_at` (TIMESTAMP(3)). The sectors (v121), the plan (v131) and the species (v156) carry their own epoch, so every stage of a galaxy says what made it; `store.scatter_epochs(conn)` reads this table. The migration gave a galaxy made earlier the plan's epoch for each kind that has rows. Not read by the fingerprint.

### `species` / `polities` / `system_owners` / `population_state`

Added in v44, filled only by the population pass
(`planetgen/population/model.py`); see "Population and politics (v44)"
above and `docs/design/habitability/population-and-politics.md`.

**`species`**: one dominant species per life world.
`generator_epoch` (v161, GEN.764) is the epoch the population pass made it under. Its indexes: the unique `name` and homeworld, `idx_species_system`, and since v177 `idx_species_spacefaring_name`, `idx_species_tech_level_name` and `idx_species_era_name`, each a sort column followed by the name, so the Species table pages by key in index order (PERF.205).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `name` | VARCHAR(255) | NOT NULL, UNIQUE | |
| `homeworld_planet_id` | BIGINT UNSIGNED | FK -> `planets.id`, `ON DELETE CASCADE`, NOT NULL, UNIQUE | |
| `star_system_id` | BIGINT UNSIGNED | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | |
| `life_class` | VARCHAR(2) | nullable | v175: the id of the life class of its homeworld (`1a`, `4c`, `6` ...); replaces `life_chemical`. |
| `life_stage` | VARCHAR(32) | NOT NULL, CHECK IN ('multicellularity','technological_civilization') | |
| `build`, `climate`, `size` | VARCHAR(16) | NOT NULL | |
| `civilization_age_years`, `era` | DOUBLE, VARCHAR(16) | nullable | NULL without a civilization. |
| `spacefaring` | TINYINT(1) | NOT NULL, default 0 | Era Interstellar or later. Indexed. |
| `energy_index`, `materials_index`, `information_index`, `medical_index`, `propulsion_index`, `defense_index` | TINYINT UNSIGNED | nullable, CHECK 0 to 7 | Added in v88 (POP.9). The six tech indexes; NULL until drawn. |
| `tech_level` | DECIMAL(3,2) | nullable | Added in v88 (POP.9). The weighted sum of the indexes. Indexed (`idx_species_tech_level_name`, with the name) for range search and the Species table's tech sort. |
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

Added in v82 (PERF.74); `system_density` and `quadrant` and a row for every sector from v181 (PERF.208). Each sector's system count, system density and Quadrant as of the last refresh (`query.refresh_sector_system_counts`), so the Sectors list reads these instead of counting every system per request and its sorts by systems, density and position page by key on the three indexes. Derived data: the refresh makes every row again, in the background and not more often than every 30 s; a sector made since has no row until it runs, and the list shows it after the others with a count of 0.

| Column | Type | Null | Notes |
|---|---|---|---|
| `sector_id` | BIGINT UNSIGNED | PK | FK-like to `sectors.id` (no constraint). |
| `system_count` | INT UNSIGNED | NOT NULL | The systems it holds. Indexed with `sector_id` (`idx_sector_system_counts_count`). |
| `system_density` | DOUBLE | NULL | v181. Systems per cubic milliparsec of the sector (`system_count` over `sectors.edge_mpc` cubed). NULL for a sector with no galaxy position. Indexed with `sector_id` (`idx_sector_system_counts_density`). |
| `quadrant` | VARCHAR(8) | NULL | v181. `'I'` to `'IV'`, the bands of `spherical.sector_quadrant`. NULL for a sector with no galaxy position. Indexed with `sector_id` (`idx_sector_system_counts_quadrant`). |

### `system_interests`

Added in v165 (MAP.305, MAP.303); replaces `habitable_systems` (v148). One row for each star system and reason it is worth a mark on the Galaxy Map, so the marker layer finds the systems with a reason inside a box without reading every planet (`query.galaxy_system_interests_in_box`). The reasons are the registry in `db/interests.py` (`inhabited` (v167): a planet or moon whose timeline reached a technological civilization or that has a terrestrial colony on it, as `query._with_life_fields` decides it, with the best score; `habitable`: a planet or moon whose overall PHI-6 score is Yellow or better, `habitability.YELLOW_FLOOR`; `uninhabitable` (v166): a scored terrestrial planet or moon under that floor; `gas-giant`: a planet with `body_type = 'g'`; `asteroid-field` (v166): a system with an asteroid belt; `civilization` (v173, MAP.311): a planet or moon whose timeline reached a technological civilization on a class that can carry life, with (v174, MAP.344) the technology level of its most advanced species as `best_score`, NULL until a species is drawn; `sentient` (v173): a sentient life score of 1 or more; `life` (v173): any of the four life scores of 1 or more, `sentient` and `life` keep no score; `class-A` to `class-V` (v172, MAP.220): a system with a planet or moon of that planet class, `n` how many, no score; `lifeclass-1a` to `lifeclass-9` (v176, MAP.351): a system with a planet or moon on which that life class lives, `n` how many, no score), in the order the map prefers them: a system's mark takes the symbol of its best world, an inhabited one before a habitable one before a barren one before a gas giant before a belt (MAP.321). Derived data: `store.record_system_interests` rebuilds a system's rows when its sector is filled, whenever its scores are worked out again and when a facility is added, changed or removed (a colony makes its world inhabited). A system with no reason, or with no galaxy position, has no row.

| Column | Type | Null | Notes |
|---|---|---|---|
| `star_system_id` | BIGINT UNSIGNED | PK | FK to `star_systems.id`, ON DELETE CASCADE. |
| `reason` | VARCHAR(24) | PK | The reason's id in `interests.REASONS`. |
| `n` | SMALLINT UNSIGNED | NOT NULL | The planets and moons of the system the reason counts. |
| `best_score` | DOUBLE | NULL | The largest habitability score among them; NULL for a reason without a score. |
| `galaxy_x_pc`, `galaxy_y_pc`, `galaxy_z_pc` | DOUBLE | NOT NULL | The system's `galaxy_*_pc`, copied (they never change). Indexed together (`idx_system_interests_position`). |

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

The discrete stellar components of a galaxy beyond the disk, one row each, of every kind in one table (v97, GEN.159): the globular clusters (GEN.159) and the satellites (GEN.237: the Magellanic Clouds, Sagittarius and eight dwarfs, with an `exponential` or `plummer` profile), and streams later. A component is a density, not a list of stars: a profile of its mass about a centre, from which a sector's stars are drawn when the sector is generated (`docs/design/galaxy/globular-clusters.md`, `docs/design/galaxy/intergalactic-and-halo-contents.md` section 5.5). The `king` profile is `planetgen/physics/king.py`: ten tables for central potentials 3 to 12, solved once and scaled by length; the `exponential` and `plummer` profiles are `planetgen/physics/component_profiles.py`, each cut at `tidal_radius_pc` (the column is set for every profile, `king_w0` for `king` alone). `planetgen/galaxy/components.py` has the `Component` record and `store.replace_galaxy_components` / `get_galaxy_components` write and read it. `planetgen plan` draws the clusters from the galaxy seed with `planetgen/generation/clusters.py` before it walks the outline (GEN.160), replacing the `cluster` rows, and places the Milky Way's eleven satellites with `planetgen/generation/satellites.py` (GEN.237), replacing the `satellite` rows; the read-back shape carries both as `GalaxyShape.components` (`planetgen/galaxy/component_field.py`), which adds their density to the sector gate.

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
                  │             ├──1───* planet_reflection_spectrum
                  │             ├──1───* planet_interior_layers
                  │             ├──1───* planet_gases
                  │             └──1───* moons
                  │                       │
                  │                       ├──1───* moon_reflection_spectrum
                  │                       ├──1───* moon_interior_layers
                  │                       └──1───* moon_gases
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

### v179 (GEN.908): the brine seep

Migration 0179 adds no column. `planets.hydrosphere` and `moons.hydrosphere` take a seventh value, `brine seep` (their CHECKs, `planets_chk_5` and `moons_chk_5`, list it), for a liquid surface whose ocean class is the chloride brine (surface water under 2e-5 of the mass): a puddle of concentrated salt water that melts below 273 K, often under thin or no air, not a sea. The migration relabels a stored `surface ocean` whose `ocean_class` is `chloride brine`; the generator gives the new label to bodies it draws from the new epoch on. The checks that ask for a `surface ocean` (a sunlit sea for life, plate tectonics, an induced field) do not take it; the ocean columns, the water solvent and the water source (`liquid`) do apply.

### v176 (MAP.351): the life class reasons

Migration 0176 adds no column. It fills `system_interests` with one reason for each life class, `lifeclass-1a` to `lifeclass-9` (`interests.LIFE_CLASS_REASONS`): a system has a row for each life class that lives on one of its planets or moons (`life_class` and `life_class_variant`), `n` how many, no score. It is filled from the planets and moons already stored; a body whose life class is NULL gives none. Nothing a generator draws changes, so there is no epoch. Since GEN.1001 a body counts only when it also has a `life_stage` (a world that suits a class but whose life has not started gives none); a database migrated before that keeps the rows of such bodies until its systems are saved again.

### v192 (GEN.1057): the age of a ring

Migration 0192 adds `ring_age_yr` (DOUBLE, NULL) to `planets`, `moons` and `rogue_planets`, after `ring_ice_fraction`: the years since the event that made the ring, a moon that fell inside the Roche limit or a comet that passed inside it (`physics/ring_origin.py`). NULL for a body with no ring, for the narrow ring of a small body's spin-orbit resonance, and for a row made before v187; nothing is backfilled. A column already there is skipped, so a stopped run can be started again.

### v191 (GEN.1034): the open-ocean tide of a world with a sea

Migration 0191 adds three nullable columns at the end of `planets` and `moons`, derived by `physics/ocean_tides.py` from the stored mass, radius, spin axis and rotation period, orbits, moons and hydrosphere: `tide_range_spring_m` and `tide_range_neap_m` (DOUBLE, metres: the equilibrium high-to-low range with every tide-raiser in line, and with the largest against the rest; `0` for an ice-covered ocean, NULL for a body without a free sea, that is any other hydrosphere, a gas giant, or a body with no spin on record) and `ocean_liquid` (VARCHAR(40): the liquid of that sea, `water` until the surface liquid inventory of GEN.1035 names others). They are the value a deep ocean would have if it followed the tidal potential freely; the range on a coast depends on the shape and depth of the basin and on its resonance and can exceed it many times. A row made before v187 keeps NULL; nothing is backfilled; a stopped run can be started again (a column already there is skipped).
### v197 (GEN.1159): the roundness ratio, the strength and the binding energy

Migration 0195 adds `roundness_ratio`, `yield_strength_pa` and `binding_energy_j` (DOUBLE, nullable) to `planets`, `moons` and `comets`, `density_g_cm3` and `ice_mass_fraction` to `comets`, `asteroid_belts` and `asteroid_fields`, skipping a column already there. All are NULL after the migration: they are worked out by `physics/body_strength.py` when a body is generated, so a row made before v195 gets them when the galaxy is regenerated. Nothing is backfilled.

### v189 (GEN.1009): the largest radius and mass of a member

Migration 0189 adds two nullable DOUBLE columns to `scatter_members`: `radius_max_au`, the largest radius the member has had (a white dwarf's red-giant or asymptotic-branch tip; the present radius for a member that has not left its main sequence), and `mass_max_solar`, the largest mass it has had (its birth mass). Both come from the member record (`generation/member_record.py`) that the tree and the build read. NULL on a row made before v187; nothing is backfilled, a row gets them when the galaxy is regenerated. A column already there is skipped, so a stopped run can be started again.

### v196 (GEN.1108): the climate of a body in seven bands of latitude

Migration 0196 adds 49 columns to `planets` and `moons` (`band1_mean_k`, `band1_swing_k`, `band1_day_hours`, `band1_night_hours`, `band1_longest_day_days`, `band1_longest_night_days`, `band1_state`, and the same for `band2` to `band7`; DOUBLE except the VARCHAR(24) states, all nullable) and drops the 14 columns they replace (`temp_polar_k`, `temp_north_k`, `temp_equator_k`, `temp_south_k`, `temp_swing_k`, `temp_day_side_k`, `temp_terminator_k`, `temp_night_side_k`, `temp_polar_swing_k`, `temp_north_swing_k`, `temp_south_swing_k`, `temp_day_side_swing_k`, `temp_terminator_swing_k`, `temp_night_side_swing_k`); `temp_mean_swing_k` stays. `planets` and `moons` carry a FULLTEXT key, so each table is rebuilt once with every change in one ALTER (a plain ALTER, not INSTANT). A column already added or already dropped is skipped, so a stopped run can be started again. All the new columns are NULL after the migration: they are worked out from the stored temperature, spin, orbit and pressure when a body is generated, so a regenerated galaxy fills them (`physics/regional_temperature.py`, `physics/latitude_bands.py`).

### v194 (GEN.1013): the members of a tree that are not stars

Migration 0194 creates `dark_members`, keyed by `star_system_id` and `member_no`: one row for each brown dwarf, neutron star, black hole or stripped helium star that hangs in a built system's orbit tree without being one of its stars. Columns: `parent_node_no` and `slot` (the numbering of `orbit_nodes`), `kind` (a CHECK on the four kinds), `mass_solar` (DOUBLE, not null) and `radius_au` (DOUBLE, nullable). A node's `mass_solar` in `orbit_nodes` is the sum of what hangs below it, these included. A row goes with its system (foreign key, ON DELETE CASCADE). A system with a dark member and no other star is a single star: its `stars` row has role `single`, `parent_node_id` the node it hangs in and its `slot`. No row is made for a system already stored; a system gets its rows when it is generated again. A table already there is kept, so a stopped run can be started again.

### v188 (GEN.1035): the liquid at the surface

Migration 0188 adds `surface_liquid` (VARCHAR(32)), `surface_liquid_fraction` and `surface_liquid_depth_km` (DOUBLE, all nullable) to `planets` and `moons`, skipping a column already there. All NULL after the migration: the liquid is found with the ground, so a regenerated galaxy fills them.
### v195 (GEN.983): the dose a moon takes from its parent's radiation belts

Migration 0195 adds `dose_belt_msv_yr` (DOUBLE, nullable) to `planets` and `moons`, skipping a column already there. It is the part of `surface_dose_msv_yr` that comes from the trapped particles of the parent planet's field at the moon's orbit (0 for a planet, NULL for a giant and on a row made before v187; `physics/radiation_belts.py`). NULL after the migration: a regenerated galaxy fills it.

### v186 (PERF.211): the Phenomena table's sorts by key

Migration 0186 adds indexes to the eight phenomenon tables (`nebulae`, `asteroid_fields`, `black_holes`, `neutron_stars`, `supernova_remnants`, `rogue_planets`, `interstellar_comets`, `quasars`), each ending in `name`: `idx_<table>_placed_name` on `(center_x_pc IS NOT NULL)`, `idx_<table>_descriptor_name` on the descriptor the list shows (`nebula_type`, `density`, `pulsar_type`, `morphology`, or the text a CASE works out for a black hole, rogue planet, interstellar comet and quasar, with the `utf8mb4_unicode_ci` collation), and, for the three tables with a radius (`nebulae`, `asteroid_fields`, `supernova_remnants`), `idx_<table>_radius_name` on `radius_ly`. The black hole and neutron star indexes lead with `star_id`, because the list leaves out the ones that belong to a star system. Some key parts are expressions (functional indexes); a query uses one only when it writes the expression the same way (`query._built_sort_key`), and `test_phenomenon_sort_keys` checks that each index serves the query's order without a sort. No column and no stored value changes; the indexes on columns build online (`LOCK=NONE`), the ones on an expression cannot (the server refuses), so while one builds reads carry on and writes to that table wait (`LOCK=SHARED`); a stopped run can be started again. With 1,000,000 rows in each of the eight tables a first page by type took 56 s and each sort now answers a page, the first or the 31st after a row, in 4 to 21 ms.

### v183 (GEN.834): the non-water ices by species

Migration 0183 adds eight DOUBLE columns (`ice_co2_mass_fraction`, `ice_ch4_mass_fraction`, `ice_nh3_mass_fraction`, `ice_n2_mass_fraction`, `ice_so2_mass_fraction`, `ice_co_mass_fraction`, `ice_hcn_mass_fraction`, `ice_other_mass_fraction`, all nullable) to `planets` and `moons`, skipping a column already there. All NULL after the migration: the species are drawn with the body's bulk composition, so a regenerated galaxy fills them.

### v187 (GEN.768): the hours outside in a day

Migration 0187 adds `outside_hours_day` (DOUBLE) and `outside_hours_limit` (VARCHAR(16)), both nullable, to `planets` and `moons`, skipping a column already there. Both are NULL after the migration: they are worked out when a body is generated or scored again, so a regenerated galaxy fills them.


### v193 (GEN.1097): the state of a day and night

Migration 0193 adds `day_night_state` (VARCHAR(24)), `night_hours` (DOUBLE) and `days_per_year` (DOUBLE), all nullable, to `planets` and `moons`, skipping a column already there. All three are NULL after the migration: they are worked out from the stored spin, orbit and polar circle when a body is generated, so a regenerated galaxy fills them.

### v182 (GEN.848): a day and night swing for every temperature

Migration 0182 adds seven DOUBLE columns (`temp_polar_swing_k`, `temp_north_swing_k`, `temp_south_swing_k`, `temp_mean_swing_k`, `temp_day_side_swing_k`, `temp_terminator_swing_k`, `temp_night_side_swing_k`, all nullable) to `planets` and `moons`, skipping a column already there. All NULL after the migration: they are worked out when a body is generated, so a regenerated galaxy fills them.

### v175 (GEN.889, GEN.890): the life class switch

Migration 0175 drops `life_chemical` (and its index) from `planets` and `moons` and adds `life_class` (TINYINT UNSIGNED, the class number), `life_class_variant` (CHAR(1)), `life_class_best` (VARCHAR(2)), `life_class_best_score` (TINYINT UNSIGNED), `life_animal_grade` (TINYINT UNSIGNED) and `life_complexity` (DOUBLE), all nullable, plus the child tables `planet_life_class_scores` and `moon_life_class_scores` (`planet_id` or `moon_id`, `life_class` VARCHAR(2), `score` 1 to 9; unique for a body and class; indexed on class and score). `species.life_chemical` becomes `species.life_class` (the class id). A row made before keeps the class its chemistry stood for (Chlorophyll a is 1a, Blue-Optimized Porphyrins 1b, Bacteriochlorophylls and Retinal 2, Melanin 3b, Zinc-Bacteriochlorophyll 5b) and nothing else: no scores, grade or complexity are worked out, so a galaxy shows the new life only after `planetgen bring-to-epoch` or a regeneration (epoch 120). The columns are added with one ALTER for each table (no `ALGORITHM=INSTANT`: both tables carry a FULLTEXT index), and a stopped run starts again where it stopped.

### v174 (MAP.344): the technology level of a civilization

Migration 0174 adds no column. It sets `best_score` of the `civilization` rows of `system_interests` to the technology level (`species.tech_level`, 0 to 7) of the most advanced species whose homeworld the reason counts, from the species already stored; NULL when none is drawn. The population pass sets it from now on (`interests.refresh_technology`), because a species is drawn after the planets are stored. Unlike the other reasons the score is not a habitability score: the node lists it as `tech` and `band` (`tuning.TECH_BANDS`) and no `tier`. Nothing a generator draws changes, so there is no epoch.

### v173 (MAP.311): the life reasons

Migration 0173 adds no column. It fills `system_interests` with the reasons `civilization`, `sentient` and `life` (`interests.CIVILIZATION`, `SENTIENT`, `LIFE`): the count of planets and moons that reached a technological civilization on a class that can carry life, that have a sentient life score of 1 or more, and that have any of the four life scores of 1 or more, from what is already stored; no score is kept. Nothing a generator draws changes, so there is no epoch. The map reads the stored scores as they are: how many worlds qualify moves with the generator.

### v172 (MAP.220): the planet class reasons

Migration 0172 adds no column. It fills `system_interests` with one reason for each planet class, `class-A` to `class-V` (`interests.CLASS_REASONS`): a system has a row for each class it has a planet or moon of, with how many, from what is already stored. Nothing a generator draws changes, so there is no epoch. The nodes list these reasons only for blocks of `nodes.CLASS_ABOVE` (3) shells or fewer, a limited number of systems for each class (`NODE_CLASS_SYSTEMS`, `NODE_SECTOR_CLASS_SYSTEMS`), so a node does not carry twenty classes for every system.

### v171 (GEN.836): the ring system of a rogue planet

Migration 0171 adds the seven ring columns of v160 (`ring_kind`, `ring_origin`, `ring_inner_km`, `ring_outer_km`, `ring_mass_kg`, `ring_optical_depth`, `ring_ice_fraction`) to `rogue_planets`, after `hp_ice_km`; NULL for a planet with no ring, a brown dwarf, and a row made before v171.

### v170 (GEN.767): the day and night shares and the hours of daylight

Migration 0170 adds six nullable DOUBLE columns at the end of `planets` and `moons`, derived by `physics/daylight.py` from the stored rotation, orbit, tilt and spin axis and the star's radius: `day_share_equator` (share of the time the star is up at the equator over a year, NULL for a world that turns once an orbit), `day_share_mean` (the same averaged over the surface), `parent_shadow_share` (a moon only: share of its time in its planet's shadow, which the two day shares already lose), `solar_day_hours` and `daylight_hours` (hours from noon to noon, and of daylight in it at the equator; NULL for a world that turns once an orbit) and `polar_circle_deg` (latitude beyond which the star stays up or down for at least a day). A row made before v170 keeps NULL; nothing is backfilled.

### v169 (GEN.721): `lid_regime`, the tectonic regime kept with the body

Migration 0169 adds `lid_regime` (VARCHAR(10), nullable) at the end of `planets` and `moons`: `mobile` (plate tectonics) or `stagnant` (one rigid plate). It was three separate coin flips on the same chance (`magnetism.mobile_lid_chance`: a liquid surface ocean and a mass of 0.5 to 5 Earth masses at 0.6, heavier at 0.3), one in the magnetic field and one in the ground dose; GEN.721 draws it once, after the hydrosphere, and the dynamo's lifetime, the ground dose's felsic crust and the interior's heat flux law all read it. A row made before v169 has NULL there and keeps the heat flux of a stagnant lid; nothing is backfilled.
### v168 (DB.108): the mass of each sector, block and region

Migration 0168 adds `sector_stats.total_mass_kg` (DOUBLE, nullable, kg: the summed mass of the sector's stars, planets, moons, asteroid belts and comets, reached through `star_systems.sector_id`, and of its black holes, neutron stars, quasars (their black hole), rogue planets, interstellar comets and asteroid fields; a black hole's, neutron star's and quasar's black hole's mass is stored in solar masses and multiplied by `constants.SOLAR_MASS_TO_KG`; nebulae and supernova remnants are diffuse gas, counted by volume under their own label slice, and are left out; a generated sector with no body holds 0; NULL before the fill or when a body of it has no stored mass), and `galaxy_block_stats.mass_sum_kg` and `region_stats.mass_sum_kg` (DOUBLE, NOT NULL, default 0, kg: the sum of the generated sectors of the block or region; a sector with a NULL mass counts as 0). The sector's sum is made exactly (the masses are added as decimals and rounded once to a double), the region's by `math.fsum`, so the same bodies give the same bits whatever order they were saved in. `db/sector_mass.py` makes the sector figures: `store.record_sector_stats` writes it with the rest of the sector's stats, `block_stats.refresh_chain` and `refresh_sectors` put it right before the blocks above add it up, and `block_stats.rebuild` (migration 0162, `planetgen fill-block-totals`) makes every sector's again and lists the sectors it could not sum (`on_skipped`). The three columns have unit entries in `db/column_units.py`. No generated value changes (the three tables are left out of the fingerprints), so no epoch.

### v167 (MAP.321): the inhabited reason

Migration 0167 adds no column. It fills `system_interests` with the reason `inhabited` (the count of planets and moons that reached a technological civilization on a class that can carry life, or have a terrestrial colony facility, and the best habitability score among them) from what is already stored. Nothing a generator draws changes, so there is no epoch.

### v166 (MAP.321): two more interest reasons

Migration 0166 adds no column. It fills `system_interests` with the reasons `uninhabitable` (a scored terrestrial planet or moon under `habitability.YELLOW_FLOOR`, with the count and the best score) and `asteroid-field` (the count of belts) from the planets, moons and belts already stored. Nothing a generator draws changes, so there is no epoch.

### v178 (PERF.205): name indexes on the last two phenomenon tables

Migration 0178 adds `idx_asteroid_fields_name` and `idx_interstellar_comets_name` (a `name` index, built online), so every phenomenon table has one and the Phenomena table pages by name in index order. No column and no stored value changes.

### v177 (PERF.205): the Species table's sort indexes

Migration 0177 adds `idx_species_spacefaring_name (spacefaring, name)`, `idx_species_tech_level_name (tech_level, name)` and `idx_species_era_name (era, name)` to `species`, online, and drops `idx_species_spacefaring` and `idx_species_tech_level`, which were their first columns. No column and no stored value changes.

### v165 (MAP.305, MAP.303): `system_interests` replaces `habitable_systems`

Migration 0165 makes `system_interests` (see its section) and fills it from the planets and moons already stored: a `habitable` row for each system with a planet or moon at or above `habitability.YELLOW_FLOOR`, a `gas-giant` row for each system with a planet of `body_type` `g`. It then drops `habitable_systems`. Nothing a generator draws changes, so there is no epoch. The Galaxy Map's nodes list its rows inside each block of 81 shells and smaller (`interests`, see the wire format).

### v164 (GEN.764): the epoch of each scatter and each species

Migration 0164 creates `scatter_epochs` (above) and adds `species.generator_epoch` (SMALLINT UNSIGNED, nullable): the epoch the population pass made the species under. A species made before v164 gets the epoch of the sector its homeworld's system is in.

### v162 (DB.96): the atmosphere colour of a planet or moon

Migration 0162 adds `atmosphere_colour` (INT UNSIGNED, nullable) to `planets` and `moons`: the colour of the body's air as one sRGB integer, `red * 65536 + green * 256 + blue` (unit key `rgb24` in `db/column_units.py`). It is worked out once with the body at the generator epoch (`planetgen.physics.body_colour`): a solid world gets the colour of the light its air scatters to the ground (GEN.418), a gas or ice giant the colour of its air at the reference level (GEN.419), from the air, the star's temperature, the ozone column and the stored haze, dust and cloud. A body with no air, no star temperature or no mixture, and a row made before v161, keeps NULL. The cloud cover is the existing `cloud_cover` (v138). The ground under the sky is the labelled stand-in reflectance until the ground's own is stored.

### v161 (DB.107): `region_stats`, the named regions of the Galaxy Map

Migration 0161 creates `region_stats` and fills it from `sectors`. A region is a 9-block of the block ladder (about 780 sectors, `galaxy/named_regions.py`) that holds at least one generated sector. A row is one such cell, `(shell_index, band_index, wedge_index)` (the level, 9, is implied; PK, all INT), with `sector_count` (its generated sectors) and `capacity` (the sectors its outline holds), so the share charted is `sector_count / capacity` (INT UNSIGNED; a standard fill of 109 sectors is about 14 percent of a cell); `clump_count` (INT UNSIGNED, the face-connected clumps those sectors make); and `extent_x_pc`, `extent_y_pc`, `extent_z_pc`, `extent_radius_pc` (DOUBLE, galaxy frame, parsecs: the bounding sphere of the largest clump, never under the region's framing sphere, which the map frames and draws the region with). A cell with no generated sector has no row, and the region does not exist. `db/region_stats.py` keeps it: `refresh_cells` recomputes a cell from its sectors (and drops the row of a cell with none), and `db/block_stats.py` calls `refresh_chain`, `refresh_sectors` and `rebuild` here whenever it puts right the blocks above a sector, so a save, a deletion and `planetgen fill-block-totals` keep the regions the same whatever order the sectors came or went in. `region_stats.region_of(conn, position_pc)`, `region_row(conn, cell)`, `charted_percent(conn, cell)` and `regions(conn, cells)` read it. No foreign key: a cell is not a row of its own. The version is detected from the table.

### v160 (GEN.732): the ring system of a planet

Migration 0160 adds `ring_kind` (`dusty`, `narrow` or `broad`), `ring_origin` (`moon_erosion`, `disrupted_moon` or `captured_body`), `ring_inner_km`, `ring_outer_km`, `ring_mass_kg`, `ring_optical_depth` (a mean over the area the ring covers) and `ring_ice_fraction` (DOUBLE, nullable) to `planets` and `moons`, after `has_helium_shell`. One summed ring system per body (`planetgen.physics.rings`, `docs/design/composition/planetary-rings.md`); NULL for a body with no ring, for every moon, and for a row made before v157, which gets its rings when the galaxy is regenerated.

### v158 (GEN.701): the habitability score of a rogue planet

Migration 0158 adds 37 score columns (`phi4` and its six domains, the six tiers, the four life scores, `l_solv`, `l_chem`, `l_ener`, `l_rad`, the fourteen `eq_*` levels with `shelter_column_g_cm2`, and `energy_flux_w_m2`) to `rogue_planets`, at the end of the table (no regional temperatures or hours outside: a rogue planet has no star), and fills them for every terrestrial row from the surface it already holds, by `physics/rogue_habitability.py`, in pages of ids (a row with `phi4` set is done). Nothing is drawn. A gas giant or brown dwarf keeps NULL.

### v157 (GEN.654, DB.92): the aquifer of a world

Migration 0157 adds `aquifer_depth_km` and `aquifer_thickness_km` (DOUBLE, nullable) to `planets` and `moons`: the depth to the top of the liquid-water layer in the pores of the crust of a world whose surface has no liquid (hydrosphere `dry`, `ice` or `vapour`), and its thickness, km ([activity-magnetism-radiation-hydrosphere.md](../design/habitability/activity-magnetism-radiation-hydrosphere.md) section 5.6). NULL for a world with none, a gas giant and a row saved earlier: the surface reservoir the layer is filled from was a draw that is not stored, so an older row is not worked out again and its stored scores stay as they were. One ALTER for each table with no `ALGORITHM=INSTANT` (both have a FULLTEXT index); a column already there is left alone. The version is detected from `planets.aquifer_depth_km`.

### v156 (DB.100): the milestone schedule of a world's pace

Migration 0156 adds `life_age_abiogenesis_gyr`, `life_age_photosynthesis_gyr`, `life_age_complex_cells_gyr`, `life_age_multicellularity_gyr` and `life_age_technological_civilization_gyr` (DOUBLE, Gyr, nullable) to `planets` and `moons`, after `life_system_age_gyr`. They are the age at which the world's `evolutionary_pace` reaches each of the five milestones (`tuning.LIFE_MILESTONES`), written with the pace by `life.apply_life_data` for every body that has a life chemistry (NULL without one). With `life_chemical`, `evolutionary_pace`, `life_stage`, `life_milestone_age_gyr` and `life_system_age_gyr` they are all the life and evolution data the life paragraph and the timeline paragraph need: `describe/life_facts.py` reads the row alone (the words of a pace are `tuning.EVOLUTIONARY_PACE_LABELS`, the text of a stage is a lookup by its key) and no draw or pace table is read when a page is written, so tuning the table later moves no saved world. A row saved before gets the ages of its stored pace (`evolution.schedule`); a row with no pace keeps NULL. The four life scores (v150) and the survival kit are stored as before. One ALTER for each table with no `ALGORITHM=INSTANT` (both have a FULLTEXT index), rows filled in pages of ids (a row with `life_age_abiogenesis_gyr` set is done). The version is detected from `planets.life_age_abiogenesis_gyr`.

### v155 (DB.71, GEN.424): the air moves to `planet_gases` and `moon_gases`

The ten `p_<gas>_kpa` columns of `planets` and `moons` are gone; each gas above zero is a row of `planet_gases` or `moon_gases` (formula, Pa), the balance gas included. The migration copies the ten partial pressures (kPa times 1000); a row made before this change has no balance-gas row until the galaxy is regenerated (generator epoch 89).

### v154 (DB.89, GEN.642): the interior of a giant

Migration 0154 adds to `planets` and `moons`: `interior_kind` (VARCHAR(16), nullable, CHECK gas_giant / ice_giant / brown_dwarf), `core_mass_earth`, `dilute_core_mass_earth`, `metallization_pressure_gpa`, `helium_depletion` (DOUBLE, nullable) and `has_helium_shell` (BOOLEAN, nullable), each at the end of the table, and creates `planet_interior_layers` and `moon_interior_layers` (the layers, outside in; see above). Nothing is backfilled: an interior needs the draws of the generator epoch, so a row made before the revision keeps NULL and no layer rows, and the regenerated galaxy replaces it. A stopped run can be started again; the version is detected from `planets.has_helium_shell`.
### v153 (GEN.699): the mantle's rheology and the heat of the interior

Migration 0153 adds `mantle_eta0_pa_s`, `mantle_activation_kj_mol` and `internal_heat_flux_w_m2` (DOUBLE, nullable) at the end of `planets` and `moons`. The first two are the draws of the mantle's viscosity law; the third is the heat that comes out of the interior over the body's surface: the radiogenic heat of the stored concentrations (v141) plus the primordial heat, the cooling of the mantle since the body was made, from a one-box convection model (`physics/interior_heat.py`, design `docs/design/composition/interior-heat.md`). `urey_ratio` keeps its column and is now computed from that history (the radiogenic share of the heat out) instead of drawn between 0.1 and 0.9. A row from before v153 has NULL in the three new columns and keeps the Urey ratio it was drawn with; nothing is backfilled. Generator epoch 72 marks it.

### v152 (GEN.700): the tide heat flow of an icy moon

Migration 0152 adds `tidal_flux_w_m2` (DOUBLE, nullable) to `planets` and `moons`, at the end of each table (a planet has no tide and keeps NULL). The generator drew the tide before and discarded it after thickening or thinning the ice lid; it is now kept and added to the heat flow the chemical energy of the ocean's life reads (`energy_flux_w_m2`, [habitability-index.md](../design/habitability/habitability-index.md) section 8). The forced eccentricity behind a tide is not stored, so a moon saved earlier keeps NULL and its stored energy flux; only a moon generated now has it. One ALTER for each table with no `ALGORITHM=INSTANT` (both have a FULLTEXT index); a column already there is left alone. The version is detected from `planets.tidal_flux_w_m2`.
### v181 (PERF.208): the Sectors list's sorts by key

Migration 0181 adds `system_density` and `quadrant` to `sector_system_counts`, each with an index ending in `sector_id`, and fills the table again so that every sector has a row (before, a sector with no systems had none). Nothing else moves, and the table is derived data. The Sectors list sorted by systems, density or position pages by key on these indexes: with 1,000,000 sectors a page took 2.5 to 8 s whatever its depth and now takes 2 to 40 ms.

### v180 (DB.70): a moon's origin and its giant's disk

Migration 0180 adds to `planets` and `moons`: `moon_origin` (VARCHAR(8)) and four DOUBLE columns, `moon_disk_lifetime_myr`, `moon_disk_alpha_over_f`, `moon_disk_f_inf` and `moon_disk_ice_line_km`, all nullable, each at the end of the table. Nothing is backfilled: the origin and the disk are draws of the generator epoch, so a row made before the revision keeps NULL and the regenerated galaxy replaces it. A stopped run can be started again (a column already there is skipped).

### v159 (GEN.441): the surface composition of every planet and moon

Migration 0160 adds to `planets` and `moons`: `ground_type` (VARCHAR(16)), `ground_grain` (VARCHAR(8)) and twelve DOUBLE columns, `ground_<type>_fraction` for `ocean`, `ice`, `ultramafic`, `anorthositic`, `basaltic`, `felsic`, `clay`, `oxidised`, `metal`, `carbon`, `evaporite` and `sulfur`, all nullable, each at the end of the table. Nothing is backfilled: a ground needs the draws of the generator epoch, so a row made before the revision keeps NULL and the regenerated galaxy replaces it. A stopped run can be started again (a column already there is skipped).

### v151 (GEN.461): the pressure under the water layer

Migration 0151 adds `water_base_pressure_pa` (DOUBLE, nullable) at the end of `planets` and `moons`: the pressure at the bottom of the water over the covered area, Pa. A row made before v151 has NULL there and keeps its water-equivalent `ocean_depth_km`, `ice_shell_km` and `hp_ice_km`; from v151 those are the real, compressed kilometres (746 km, not 1047 km, for a tenth of an Earth mass of water on an Earth mass). The generator epoch marks the old rows. Nothing is backfilled.

### v150 (GEN.649): the four life scores

Migration 0150 adds `life_microbial`, `life_flora`, `life_fauna` and `life_sentient` (TINYINT UNSIGNED, nullable) to `planets` and `moons`, fills them for every row that has a habitability score from the two likelihoods (`phi_bio`, `phi_cpx`) and the evolutionary path it stored (`life_chemical`, `life_stage`, `evolutionary_pace`, `life_milestone_age_gyr`, `life_system_age_gyr`) by `planetgen/physics/life_scores.py`, in pages of ids (a row with `life_microbial` set is done), and drops `phi_bio`, `phi_cpx` and `phi_tech`. The human operability score has no successor. A gas giant or an unscored body keeps NULL.

### v149 (OPS.69): the wiki links are gone

Migration 0149 drops `sectors.wiki_url`, `star_systems.wikijs_url` and `star_systems.mediawiki_url`. The wiki upload that wrote them was taken out of the site, and nothing read them any more.

### v148 (MAP.220): the `habitable_systems` table

Migration 0148 makes `habitable_systems` (`star_system_id` PK with a cascading FK, `habitable_worlds` SMALLINT UNSIGNED, the system's `galaxy_x_pc`, `galaxy_y_pc`, `galaxy_z_pc` and a key on the three) and fills it from the planets and moons already stored. It is a table of its own because `star_systems` has a FULLTEXT key, so adding a column there would copy the whole table. The Galaxy Map's nodes list its rows inside each block of 81 shells and smaller (`worlds` then, `interests` since v165, see the wire format).

### v147 (MAP.259): `galaxy_block_lights`, the light of the pre-placed objects in each block

Migration 0147 creates `galaxy_block_lights` and fills it from the stored scatter. The single-space map draws a block it has not opened as one node: its brightest objects and a glow made from the light of the rest, so a node needs the summed light of what is in a block without reading the objects. A row is one block `(m, shell_index, band_index, wedge_index)` of the eleven levels above a sector that holds at least one bright star no system was built around (`bright_stars.star_system_id IS NULL`) or scattered phenomenon not yet built (`phenomenon_scatter.built_at IS NULL`). `star_count`, `star_luminosity_sol`, `star_lpx`, `star_lpy`, `star_lpz` (luminosity times position, solar luminosities times parsecs, so a centre of light is each over `star_luminosity_sol`) and `star_lt` (luminosity times temperature, so the mean temperature of the light is `star_lt / star_luminosity_sol`) sum the stars; `point_count`, `point_luminosity_sol`, `point_lpx`, `point_lpy` and `point_lpz` sum the phenomena the same way (`map_luminosity_sol`, no temperature). Built objects are counted in `galaxy_block_stats` (the generated sector's total luminosity), so a block's two tables never count one light twice. `db/block_lights.py` keeps it, as `block_stats.py` keeps its table: a level-3 block is summed from the objects of its sectors, a bigger block from its children's rows; `refresh_chain` puts right the blocks above one sector (a deletion), `refresh_sectors` those above the sectors a `galaxy` run saved (once at its end), and `rebuild` makes every row, a unit of sector shells at a time (the shells of one level-729 block), committing each unit, then the levels above the unit from the stored rows. A scatter, a clear and the map-luminosity fill rebuild the table at their end, and `planetgen fill-block-totals` makes it (and `galaxy_block_stats`) again. No foreign keys: a block is not a row of its own. The version is detected from the table.

### v146 (MAP.258): the block ladder's root and `galaxy_block_stats` for every level

Migration 0146 widens `galaxy_block_stats.m` from SMALLINT UNSIGNED to INT UNSIGNED, fills the table again and makes the nebula block tables again. The block ladder (`galaxy/drill.py`) was 243, 81, 27, 9 and 3 shells deep and now runs up from 243 through 729, 2,187, 6,561, 19,683 and 59,049 to a root of 177,147 shells (708 kpc, past the grid's 262 kpc reach): the root's one block shell is one row of three wedges and holds the whole galaxy, so a walk from the root starts with three blocks. The single-space map's nodes are these blocks. A level's rows and wedges are its parent's counts times a rounded ratio, so putting levels above 243 moved the rows and wedges of every level below it: a block `(243, shell, band, wedge)` is another block than it was, and every stored per-block aggregate names the wrong block. The table gets a row for each of the eleven levels above a sector (a bigger block adds up the stored rows of its children, as before), and the 177,147 level does not fit SMALLINT UNSIGNED's 65,535. The migration empties the table and fills it from `sectors` and `sector_stats` a shell at a time, and, when `nebula_block_state` holds its row, makes `nebula_block_stats` and `nebula_block_classes` again from the nebulae (`db.nebula_blocks.rebuild`, the five stage levels they keep), so a stopped run can be started again. The stages and map windows still stop at the level-243 blocks (`drill.STAGE_LEVELS`). The version is detected from `galaxy_block_stats.m` being INT UNSIGNED.

### v145 (GEN.443): `stars.uv_c_share`

Migration 0145 adds `uv_c_share` (DOUBLE, nullable) after `snow_line_km` on `stars`: the share of the star's blackbody UV-C that it gives off, drawn once per ordinary star on its own labelled stream (`uv_c_share:` plus its mass, age and temperature) so no draw of the run moves. A row saved before the revision keeps NULL, which the radiation model (`physics/body_radiation.py`) reads as the blackbody; the regenerated galaxy (the generator epoch of this change) draws it. A column already there is skipped, so a stopped run can be started again.

### v144 (DB.69, GEN.428): the bulk composition of every planet and moon

Migration 0144 adds to `planets` and `moons`: `metal_mass_fraction`, `rock_mass_fraction`, `ice_mass_fraction`, `envelope_mass_fraction` (DOUBLE, nullable; with `water_mass_fraction` they add up to 1), `rock_feo_fraction`, `rock_al2o3_fraction`, `rock_mg_si_ratio`, `envelope_helium_fraction`, `envelope_z_fraction`, `porosity` (DOUBLE, nullable), `carbon_rock` (BOOLEAN, nullable) and `composition_origin` (VARCHAR(8), nullable, CHECK drawn / inverse), each at the end of the table. Nothing is backfilled: a composition needs the star's abundances and the draws, so a row made before the revision keeps NULL and the regenerated galaxy of the generator epoch replaces it. A stopped run can be started again; the version is detected from `planets.composition_origin`. The mass, radius and density of a body now follow from the composition ([bulk-composition.md](../design/composition/bulk-composition.md)); `water_mass_fraction` is drawn with the rest and `hydrosphere` reads it.

### v143 (GEN.666): `evolutionary_speed` is gone, the pace is stored once

Migration 0143 drops `planets.evolutionary_speed` and `moons.evolutionary_speed`. A world's pace (`fast`, `normal` or `slow`, a `tuning.EVOLUTIONARY_TIMELINES` key) was stored twice, once as the speed its life sentence named and once as the pace its timeline was built on, and the two were separate draws, so a page could name two paces. The generator now draws it once and stores it in `evolutionary_pace`, for every body that has a life chemistry (a moon or a world outside the habitable zone has a pace and no timeline: `life_system_age_gyr` is NULL). A row that had a speed and no pace takes the speed as its pace; a row with both keeps its pace, which its ages and milestone were built on. The version is detected from `planets.evolutionary_speed` being gone.

### v141 (GEN.459): the radiogenic content of a body

Migration 0141 adds `c_k0_ppm`, `c_th0_ppm`, `c_u0_ppm` (the element's mass in the bulk silicate at formation, ppm by mass), `al26_ratio` (26Al/27Al at formation) and `urey_ratio` (the radiogenic share of the heat flowing out, today) to `planets` and `moons`, DOUBLE, nullable, at the end of each table. A body's interior heat flux (`physics/radiogenic.heat_flow_w_m2`) and its ground dose (gamma and radon, `physics/radiation.ground_dose_msv_yr`) read these five and nothing else (the Turcotte and Schubert mantle values, the 0.67 Urey ratio and the `a_rad` draw are gone). A row from before v141 has NULL there and keeps the heat and the dose it was generated with; nothing is backfilled. The version is detected from `planets.c_k0_ppm`.

### v140 (DB.73): one parent link for every kind

Migration 0140 gives every orbited-body table the same two columns the planets and moons have carried since v102: `parent_type` (VARCHAR(16), not null; `star`, `planet`, `system`, `node`, `rogue_planet` or `free`, a CHECK) and `parent_id` (BIGINT UNSIGNED, nullable: the row of the table `parent_type` names, or for `node` the `orbit_nodes.node_no` of the row's own system; NULL for `free`), with an index `idx_<table>_parent` on both. They are added to `stars`, `asteroid_belts`, `comets`, `black_holes`, `neutron_stars`, `rogue_planets` and `interstellar_comets`, and `planets` and `moons` take the wider `parent_type`. A row saved before is backfilled: a star from `parent_node_id` (a star with none is `system`), a belt or comet from `star_id` (`system` when it has none), a black hole or neutron star from `star_id` (`free` when it has none), a rogue planet or interstellar comet as `free`.

`moons.planet_id` and `moons.star_system_id` are nullable now. A moon of a free planet has `parent_type` `rogue_planet`, the rogue planet's row in `parent_id` and in the new `moons.rogue_planet_id` (a foreign key, ON DELETE CASCADE, indexed), and no planet or system. Its ID is a body number under the rogue planet's own ID (`insert_moon_of_rogue_planet`, `load_moons_of_rogue_planet`), and the fingerprint files it with that planet. A planet that leaves its star keeps its row and its ID: `store.set_parent` rewrites only the link (`free` clears `parent_id`). There is no foreign key for a polymorphic `parent_id`; `parent_type` is the discriminator. The `sector_objects` view lists the bodies of systems only, so a moon of a free planet is not in it; the object index, `resolve_object` and the search panels do not read these moons yet (follow-up items). Every step is skipped when its column is already there, so a stopped run can be started again. The version is detected from `interstellar_comets.parent_id`, the last column the migration adds.

### v138 (GEN.414, DB.68): the haze, dust and cloud depths of a planet or moon

Migration 0138 adds `tau_haze`, `tau_dust`, `tau_cloud` (the optical depths at 550 nm), `cloud_cover` (the share of the sky under cloud, 0 to 1) and `aerosol_draw` (the number in [0, 1) that picks the optics inside their brackets) to `planets` and `moons`, all DOUBLE and nullable. They are drawn from the body's own labelled stream (`aerosols:` plus its class, mass, radius, distance and temperature), so no draw of the run moves; SQL cannot make them from a stored row, so a row saved before the revision keeps NULL, which the radiation model reads as no aerosol, and the regenerated galaxy (the generator epoch of this change) draws them. A column already there is skipped, so a stopped run can be started again; the version is detected from `planets.tau_cloud`.

### v137 (GEN.656): regional temperatures and the day and night swing

Migration 0137 adds eight DOUBLE columns (`temp_polar_k`, `temp_north_k`, `temp_equator_k`, `temp_south_k`, `temp_swing_k`, `temp_day_side_k`, `temp_terminator_k`, `temp_night_side_k`, all nullable) to `planets` and `moons`, one `ALGORITHM=COPY` statement per table. They are a pure function of stored columns, so there is no backfill and no generator epoch: `planetgen rescore-habitability` fills them for a galaxy made earlier. See docs/design/habitability/survival-equipment-fields.md section 5.5.

### v136 (GEN.655, DB.93): the survival equipment levels replace `equipment_tier`

Migration 0136 adds the `eq_*` TINYINT columns, `shelter_column_g_cm2` and `outside_hours` to `planets` and `moons` (all nullable) and drops `equipment_tier`, in one `ALGORITHM=COPY` statement per table (a FULLTEXT index on those tables refuses an instant drop). The levels are a pure function of the stored conditions, so there is no backfill and no generator epoch: `planetgen rescore-habitability` fills them for a galaxy made earlier. See docs/design/survival-equipment-fields.md.

### v135 (GEN.413): `planets.ozone_column_du` and `moons.ozone_column_du`

Migration 0135 adds `ozone_column_du` (DOUBLE, nullable) at the end of `planets` and `moons`: the ozone column in Dobson units (1 DU = 2.687e16 molecules/cm2), drawn with the body's surface dose from its oxygen partial pressure, gravity and its star's temperature (`planetgen/physics/ozone.py`: the Kozakis et al. 2022 fit by oxygen level and host, a structure factor of 0.85 to 1.71, the abiotic branch below 0.01 percent of the present oxygen; docs/design/atmosphere-colour.md section 9). 0 for a world with no oxygen, NULL for a gas giant and for a row made before the revision (nothing was drawn for it; no backfill, the regenerated galaxy has them). `uv_surface_index` and `ozone_loss_flag` read it. The version is detected from `planets.ozone_column_du`.

### v139 (MAP.143): `habitable_worlds` on `sector_stats` and `galaxy_block_stats`

Migration 0139 adds `sector_stats.habitable_worlds` (INT, nullable: the planets and moons in the sector's systems that are potentially habitable: an overall PHI-6 score that is Blue or Green, `habitability.POTENTIALLY_HABITABLE_MIN_SCORE`, read from the stored `phi4`) and `galaxy_block_stats.habitable_worlds` (BIGINT UNSIGNED, not null, default 0: their sum over the block's generated sectors). `store.record_sector_stats` writes the sector's count when it is filled (`store.sector_habitable_worlds`), `forget_sector_fill` clears it, and `db/block_stats.py` adds it up the ladder like the other sums; the stage API (`query.galaxy_stage`, `galaxy_block_stats`) reports it as `stats.habitable_worlds`, which the Galaxy Map's Color by Habitable worlds reads. The migration counts every filled sector from the scores already stored, a shell at a time, and rebuilds the block totals. The version is detected from `sector_stats.habitable_worlds`.

### v142 (GEN.412): the albedo and the greenhouse quantity on a planet or moon

Migration 0142 adds `bond_albedo`, `greenhouse_factor` (both DOUBLE, nullable) and `albedo_basis` (TINYINT UNSIGNED, nullable) to `planets` and `moons`, at the end of each table. The generator keeps the two numbers it draws for the temperature (`planetgen/physics/climate_draws.py`; [atmosphere-colour.md](../design/composition/atmosphere-colour.md) sections 5 and 8). A row saved before the revision gets them back with no draw: its stored scale height holds the temperature it would have with no air, so `1 - A = 4 sigma T_e^4 / S` and `G = (T_s / T_e)^4 - 1` follow (`albedo_basis` 1). Where there is no solution (no starlight, a temperature on the 2.725 K floor, an albedo outside 0 to 1, a pressure above the barometric law) the albedo is a seeded draw from the class range and G is worked out from the stored temperature, which never moves (`albedo_basis` 2). A stopped run can be started again; the version is detected from `planets.bond_albedo`.


### v134 (GEN.267): the star count, its cap and the multi-star prevalence on a stored config

Migration 0134 adds `system_configs.stars` (TINYINT UNSIGNED, nullable: the star count the system was forced to, `--stars N`), `max_stars` (TINYINT UNSIGNED, not null, default 12: the cap on the star count, `--max-stars`) and `prevalence_multi_star` (DOUBLE, not null, default 0: the percent from the usual chance of five or more stars, as the other prevalence columns), each at the end of the table. The defaults are what every config before the revision meant, so no row is touched or backfilled. The version is detected from `system_configs.max_stars`.

### v133 (GEN.630): `stars.role` 'member' and `stars.slot`

Migration 0133 lets a system keep a `stars` row for every one of its stars ([multistar-and-compact-systems.md](../design/systems/multistar-and-compact-systems.md), section 6.1). `stars.role` takes a fourth value, `member`, for each star beyond the pair (`primary`, `secondary`); a single star is still `single`. `stars.slot` (TINYINT UNSIGNED, nullable) is which child of its parent node the star is, 0 or 1, the numbering of `orbit_nodes.slot`: `parent_node_id` alone cannot say it for a node whose two children are both stars. A pair already stored gets `slot` 0 for its primary and 1 for its secondary where it hangs on a node. The fill still builds single stars and pairs, so every row it writes today is one of those; a system with a third star is written and read back (`StarSystem.member_stars`) but nothing generates one until GEN.631.

### v132 (GEN.458): each star's own snow line, and whether a body formed beyond it

Migration 0132 adds `stars.snow_line_km` (DOUBLE, not null): the star's own snow line, drawn once around the 170 K blackbody position for its luminosity with 0.2 dex of scatter (cut at 3 sigma) and stored with the star, so the water share, the small-body zones and the moons read one value (the fixed 2.7 AU constant is gone). A close pair's merged proxy has no column of its own: its line is the square root of the sum of the squares of its two stars' lines, computed from the stored pair. `planets.born_beyond_snow_line` and `moons.born_beyond_snow_line` (BOOLEAN, nullable) record whether the body formed at or beyond its star's line (a circumbinary planet is compared with the pair's line, a moon takes its planet's flag); belts and comets read the stored star line when GEN.442 lands. A row made before the revision gets the blackbody position and the flag computed from its distance; the regenerated galaxy replaces them with drawn values. A stopped run can be started again; the version is detected from `stars.snow_line_km`.

### v131 (OPS.58): the generator epoch, the fingerprint format and the battery on the galaxy record

Migration 0131 adds `generator_epoch`, `fp_spec`, `battery_digest` and `battery_cases` to `galaxy_shape` and `generator_epoch` to `generation_runs`, each at the end of its table and nullable. A plan that draws a new galaxy seed writes the first four with the code's version key; every run writes the epoch it generates under. `planetgen.galaxy.version_check` reads them: a different epoch is an error, a different battery digest at the same epoch a warning that names the first case that differs ([generation-determinism.md](../design/generation/generation-determinism.md) section 5.2). Rows from before the migration keep NULL, since the epoch they were made under was never recorded, and nothing is backfilled. Both tables are small (one row, and one row a run), so each change is one statement; a column already there is skipped, so a stopped run can be started again. The version is detected from `generation_runs.generator_epoch`.

### v130 (DB.85): globular clusters and superbubbles stored as containers

Migration 0130 adds the table `galaxy_containers` ([globular-clusters.md](../design/galaxy/globular-clusters.md), [superbubble-density.md](../design/galaxy/superbubble-density.md), UX.144). One row per globular cluster or superbubble: `uid` (its object ID, unique), `kind` (`globular_cluster` or `superbubble`), `name`, `component_id` (the `galaxy_components` row of a cluster; removing the component removes the container) or `bubble_col`, `bubble_row` and `bubble_n` (the column and draw of the bubble in the superbubble field; unique together), `center_x_pc`, `center_y_pc`, `center_z_pc` (the galaxy frame), `extent_pc` (a cluster's tidal radius or a bubble's radius, at least 1) and `shell_index`, `band_index`, `wedge_index` (the sector its centre is born in, the region it belongs to). A plan writes them from the galaxy seed, so a replan replaces them; `db/containers.py` finds the containers that hold a point and what a container holds. No other table changes.

### v129 (DB.63): ships and stations

Migration 0129 adds the table `vessels` for the ships and stations of 9.3 (USR.10). Nothing writes or reads it yet and no other table changes. One row per vessel: `uid` (its object ID, unique; the object-kind code is added to `names/object_id.py` with the code that writes the table), `kind` (`ship` or `station`), `name`, `origin` (`generated` or `user`), `owner_id` and `owner_name` (the account's id in the control database and its name as it was, NULL for a generated vessel), the position (`center_x_pc`, `center_y_pc`, `center_z_pc`), `mass_kg` and the movement vector (`velocity_x_kms`, `velocity_y_kms`, `velocity_z_kms`) it had at `epoch_unix`, which are the fields of the SpatialPosition3D of GEN.294, `course` (JSON, the route plan, NULL while the vessel holds still), `seed` and `version_key` (what a generated vessel was made under) and `created_at`; `idx_vessels_owner` lists one owner's vessels. A table already there is left alone, so a stopped run can be started again. The version is detected from the table.

### v128 (DB.64): stored constellations

Migration 0128 adds the table `constellations` for the view from a planet ([sky-view.md](../design/orbits/sky-view.md), VIEW.3): the constellations of a system's sky, stored so a planet keeps the same ones each time. Nothing writes or reads it yet and no other table changes. One row per constellation: `system_uid` (the system's object ID) and `constellation_index` (its number in that system, unique together), `name_id` (the object ID of its name, the `constellation` kind of VIEW.4; NULL until it is named), `members` (JSON, the object IDs of the member stars), `edges` (JSON, the pairs of members the lines join), `seed` and `version_key` (what it was built under, so it can be rebuilt and told stale) and `created_at`. A table already there is left alone, so a stopped run can be started again. The version is detected from the table.

### v127 (DB.61): unnamed clusters and the star-to-cluster link

Migration 0127 makes `galaxy_components.name` nullable, so a cluster nobody has named can be a row; the unique key over `(kind, name)` stays, and since MySQL lets any number of rows hold a NULL in a unique key, it still holds over the named rows. It also adds `cluster_id` (BIGINT UNSIGNED, NULL) to `stars` and `bright_stars`, at the end of each table: a reference to `galaxy_components.id` (`fk_stars_cluster`, `fk_bright_stars_cluster`, cleared if the cluster row is removed) with its own index (`idx_stars_cluster`, `idx_bright_stars_cluster`), so a star can point to its cluster as its parent ([globular-clusters.md](../design/galaxy/globular-clusters.md) has the cluster components). Every existing cluster keeps its name, every star's link is NULL, and nothing writes or reads the link yet. Each star table is rebuilt once in one `ALGORITHM=COPY` statement (a foreign key cannot be added in place with checks on, and `stars` has a FULLTEXT index), which takes time in proportion to its rows; a step already done is skipped, so a stopped run can be started again. The version is detected from `stars.cluster_id`.

### v126 (DB.39): edit request comments and notices

Migration 0126 adds two tables for the 10.0 edit requests ([roles-and-requests.md](../design/interface/roles-and-requests.md)); a request is a `change_plans` row (DB.36), so the requester, the approver role and the decision are not repeated here. Nothing reads them yet and no other table changes. `request_comments` is the conversation on a request: `plan_id`, `author_id` (the account's id in the control database, NULL when none), `author_name` and `author_role` (`user`, `moderator` or `admin`) as they were, `body` and `created_at`, listed in order by `idx_request_comments_plan`. `request_notifications` is one row per notice owed to an account: `plan_id`, `recipient_id`, `kind` (what the notice is about, such as a comment or a decision), `created_at` and `sent_at` (NULL until it is sent); `idx_request_notifications_recipient` finds a recipient's unsent notices. Both tables go with their plan. A table already there is left alone, so a stopped run can be started again. The version is detected from `request_comments`.

### v125 (DB.60): the density a sector was filled under

Migration 0125 adds two columns to `sectors`, at the end of the table: `density_multiplier_used` (DECIMAL(10,4), NULL) and `override_id` (BIGINT UNSIGNED, NULL). They record the multiplier of the region density override (`region_density_overrides`, v116) that applied when the sector was filled, and that override's `id`, so a rebuilt galaxy can tell which setting a sector was filled under ([map-ui-and-frontend-libraries.md](../design/interface/map-ui-and-frontend-libraries.md) section 8). Both are NULL for a sector filled with no override and for every sector saved before the revision. `override_id` has no foreign key: an override row is never deleted, and `sectors` is made before that table in a new galaxy. Nothing writes or reads the columns yet. The columns are added in one `ALGORITHM=COPY` rebuild (`sectors` has a FULLTEXT index, so an instant change is refused), and a column already there is skipped, so a stopped run can be started again. The version is detected from `override_id`.

### v124 (DB.38): region density overrides

Migration 0124 adds the table `region_density_overrides`, the append-only list of an admin's density multipliers for a region ([map-ui-and-frontend-libraries.md](../design/interface/map-ui-and-frontend-libraries.md) section 8, [edit-everything.md](../design/interface/edit-everything.md)); no other table changes and nothing reads it yet. One row per setting: `shape` (`sphere`, `box`, `shell`, `band`, `wedge` or `sectors`) and `params` (JSON, the parameters of that shape, the region object of GEN.302), `multiplier` (DECIMAL(10,4), at least 0.0001, so above zero: the 0.1 to 10 range belongs to the interface, so a typed value beyond it can carry a warning), `set_by_id` (the account's id in the control database, NULL when none), `set_by_name`, `set_at`, an optional `reason` (up to 500 characters) and `undone_at`, set once when the setting is undone and the only change a row gets. `idx_region_density_overrides_active` lists the settings still in force (`undone_at IS NULL`) in order. A table already there is left alone, so a stopped run can be started again. The version is detected from the table.

### v123 (DB.37): lock marks

Migration 0123 adds the table `object_locks`, a standing keep mark on an object or on one aspect of it ([edit-everything.md](../design/interface/edit-everything.md)); no other table changes and nothing reads it yet. The primary key is `(object_uid, aspect)`: `object_uid` is the 96-bit object ID (BINARY(12)) and `aspect` a name from the aspect catalogue (ADM.62) or `*` for the whole object (the default), so an object has at most one lock per aspect. `locked_by_id` is the locking account's id in the control database (NULL when none), `locked_by_name` its name when it locked, `locked_at` when, and `reason` an optional note of up to 500 characters; `idx_object_locks_locked_at` lists the newest locks. A table already there is left alone, so a stopped run can be started again. The version is detected from the table.

### v122 (DB.36): change plans, their operations and the edit journal

Migration 0122 adds three tables for the 8.6 edit-everything plans ([edit-everything.md](../design/interface/edit-everything.md)); no row of any other table changes and nothing reads them yet.

`change_plans`, one row per plan: `requester_id` (the account's id in the control database, NULL when none), `requester_name` and `requester_role` (`user`, `moderator` or `admin`) as they were when the plan was asked for, `request_text` (the request in the requester's words), `state` (`draft`, `requested`, `proposed`, `approved`, `applied`, `undone`, `stale`, `rejected` or `withdrawn`), `approver_role` (`moderator` or `admin`, the role that must approve it; NULL when the admin confirms their own plan), `decided_by_id`, `decided_by_name`, `decided_at` and `decision_note` (who approved or rejected it, when, and why), `apply_job_id` (the control database job that applied it) and `created_at` and `updated_at`. The requester, approver and state fields are here from the start, so the 10.0 roles (a request is a plan, DB.39) need no rebuild. `change_plan_ops`, one row per operation: `plan_id` (removed with its plan), `op_order` (unique within the plan), `kind`, `target_uid` (the object's 96-bit ID, BINARY(12)), `before_json` and `after_json` (the values as JSON) and `status` (`pending`, `applied`, `undone`, `skipped` or `failed`). `edit_journal`, the append-only list of every applied operation: `plan_id` (a plan with journal rows cannot be deleted), `op_id`, `object_uid`, `sector_id`, `aspect`, `before_json`, `after_json`, `generator_epoch` (the epoch the sector was under, so a row that no longer applies after an epoch change is listed and not dropped) and `applied_at`, with keys by object (`idx_edit_journal_object`), by sector and by plan. A table already there is left alone, so a stopped run can be started again. The version is detected from `change_plans`.

| `generator_epoch`, `run_id` | SMALLINT UNSIGNED, BIGINT UNSIGNED | nullable, `run_id` indexed, foreign key to `generation_runs.id` | Added in v117 (DB.16), replacing the four text columns of v67 (DB.7: `version_key`, `planetgen_version`, `python_version`, `platform`, dropped). `generator_epoch` is the epoch of `util/generator_epoch.py` the code ran under (NULL on a sector saved before v117, when it was not recorded); `run_id` is the `generation_runs` row of the run that saved the sector, which holds the packed version key (`planetgen/galaxy/version_key.py`) and the release, Python and platform strings. A sector saved outside a run points at one shared run row (command `save`) per version of the code. NULL `run_id` on a sector generated before v67. The Generate page and `planetgen galaxy` warn when the running version differs from a sector's (`galaxy/version_check.py`, `store.sector_versions`). |
### v121 (DB.16): the epoch and run on each sector

Migration 0121 replaces the four version text columns of `sectors` (v67, DB.7: `version_key`, `planetgen_version`, `python_version`, `platform`) with `generator_epoch` (SMALLINT UNSIGNED, NULL) and `run_id` (BIGINT UNSIGNED, NULL, foreign key `fk_sectors_run` to `generation_runs.id`, indexed by `idx_sectors_run`), added at the end of the table. `generation_runs` already holds the key, release, Python and platform of each run, so a sector no longer repeats them: `store.sector_versions` and the version-key check (`db/check.py`) read them through the run. `generator_epoch` is the epoch of [generation-determinism.md](../design/generation/generation-determinism.md) section 5 (`util/generator_epoch.py`) the code ran under; a sector saved before v121 has NULL, since the epoch was not recorded. A sector saved by a run carries that run's id; one saved outside a run (a test, the admin pages) points at one shared `generation_runs` row per version of the code (command `save`, `store.sector_run_id`). The migration gives each distinct set of the four old values one run row (command `migrated`, outcome `ok`, started and finished at the first such sector's `created_at`), points those sectors at it in batches, then drops the four columns and `idx_sectors_version_key`; a sector with no recorded version (made before v67) keeps a NULL `run_id`. `sectors` has a FULLTEXT index, so each of its two changes (adding the columns, dropping the old ones) is one `ALGORITHM=COPY` rebuild, and a step already done is skipped, so a stopped run can be started again. The version is detected from `sectors.run_id`.

### v120 (DB.49): the leaf digest and the shell-band tree

Migration 0120 adds `sectors.leaf_digest` (BINARY(32), NULL) and makes the table `sector_tree_nodes`. The leaf digest is the sector's fingerprint (the canonical encoder of GEN.136, SHA-256) stored at save time; it is NULL for a sector saved before it. `sector_tree_nodes` holds the digests above the leaves ([generation-determinism.md](../design/generation/generation-determinism.md) section 4.3): `level` 1 is a band node (one per `shell_index` and `band_index`, over the leaves of its wedges) and `level` 2 a shell node (one per `shell_index`, `band_index` 0, over its band nodes); the galaxy root is worked out on demand and not stored. The primary key is `(level, shell_index, band_index)`. `digest` is the node's SHA-256 (NULL until it is first made), `child_count` the number of leaves (level 1) or band nodes (level 2) under it, so a region's digest can include the count, `dirty` is 1 while a save or an edit has changed something under it since `digest` was made, and `updated_at` says when the row last changed; `idx_sector_tree_nodes_dirty` finds the nodes to rehash. No row is written and no code reads either one until GEN.305. The column is added in one `ALGORITHM=COPY` rebuild (`sectors` has a FULLTEXT index, so an instant change is refused), and a table already there is left alone, so a stopped run can be started again. The version is detected from the table.

### v119 (DB.30): the per-species environment limits

Migration 0119 adds the table `species_limit_values`, the environment limits a species does not share with the human visitor ([habitability-index.md](../design/habitability/habitability-index.md) section 10). The primary key is `(species_id, limit_name)`: one row for each limit the species replaces, `limit_name` the name of a row of `planetgen.physics.protection_limits.LIMIT_ROWS` (checked in code, by `species_limits.validate`, not by a constraint, so a new limit needs no migration) and `value` the species' value in that limit's unit (kPa, deg C or Sv/yr). A limit with no row is the human value, so a species with the human limits has no rows. `fk_species_limit_values_species` removes a species' rows with the species. The table is empty until the species editor (POP.13) writes it; `species_fit.stored_overrides` reads it and a stored value wins over the one generated from the species' homeworld. No other table changes. A table already there is left alone, so a stopped run can be started again. The version is detected from the table.

### v118 (DB.56): the orbital event record

Migration 0118 adds two tables; `planetgen/db/orbital_events.py` (GEN.197: `record_event`) writes them from a `collision_outcome.Outcome`, the participants first and then the products, and nothing reads them yet (GEN.143's report does); the orbital update that finds the pairs (`find_contact`) is not built, so no row is made in a run yet. `orbital_events` has one row per event: `kind` (`collision`, `merger`, `split`, `ejection`, `capture`, `unbind` or `warning`, checked by `chk_orbital_events_kind`), `outcome_code` (the Record column of the outcome table in docs/design/orbits/collisions-and-mergers.md section 5.1, such as `merger`, `hit_and_run`, `giant_merger` or `star_created`) and `regime`, `event_time_unix` (when it happened) and `detected_time_unix` (when the update found it), `seed_string` (what the event's draws were seeded with, so the outcome can be replayed) and `run_id` (the `generation_runs` row, with no foreign key so a purged run keeps its events). The numbers the outcome was decided from are nullable DOUBLE columns in km, km/s, J/kg and seconds: `separation_start_km`, `separation_impact_km`, `impact_parameter_km`, `critical_impact_parameter_km`, `v_inf_kms`, `v_esc_kms`, `v_imp_kms`, `q_r_j_per_kg`, `q_star_rd_j_per_kg`, `erosion_x`, `interacting_mass_fraction`, `catastrophic_velocity_kms`, `window_s`, with `window_capped` and `focusing_radius_used` as 0/1 flags. The conservation check is `mass_error` and `momentum_error` (relative), `conservation_tolerance`, `impact_energy_j` and `heat_share`. The products are `lost_mass_kg`, `merged_mass_kg`, `new_radius_km`, `spin_breakup_ratio`, the fusion lines used (`deuterium_limit_mj`, `hydrogen_limit_sol`) and `borderline` (0/1, within 3 percent of a line). The context is `containing_uid` (BINARY(12), the nebula or remnant the event happened in), `influence_body_count`, `sectors_dirtied`, `suppressed_duplicates`, `already_overlapping`, `numeric_warning` (0/1, the position resolution exceeds 1 percent of the summed radii) and `flags` (a comma-separated list), then `recorded_at`. Indexes on the event time, the kind with the event time, and the run.

`orbital_event_bodies` has one row per body of an event, keyed `(event_id, role, ordinal)` with a cascading foreign key to `orbital_events.id`: `role` is `participant` or `product`, `uid` the 96-bit object ID (`BINARY(12)`, no foreign key to an object table because a destroyed body's ID is retired and never reused), `object_kind`, `fate` (nullable; `destroyed`, `survived`, `absorbed` or `created`), and the body's state at the event: `name`, `object_class`, `mass_kg`, `radius_km`, `composition`, `sector_id`, `x_pc`, `y_pc`, `z_pc`, `vx_kms`, `vy_kms`, `vz_kms`, `admin_edited` (0/1) and a `note`. Both are left out of a fingerprint (a record of what the update did, not content). The version is detected from the table `orbital_events`.

### v117 (DB.66): the galaxy family and the numbers of its components

Migration 0117 adds fifteen columns to `galaxy_shape` in one statement, after `pinned_stellar_mass_solar` (docs/design/galaxy/galaxy-types.md section 10). `family` (VARCHAR(16) NOT NULL, default `spiral`; one of `spiral`, `lenticular`, `elliptical`, `irregular`, `dwarf_spheroid`, checked by `chk_galaxy_shape_family`). The bar: `bar_present` (TINYINT), `bar_axis_ratio_b`, `bar_axis_ratio_c` (DOUBLE, in (0, 1]) and `bar_angle_rad`, and the thick disk's `thick_disk_amplitude` (DOUBLE, at least 0), all NULL. The Sersic spheroid, all NOT NULL: `spheroid_amplitude` (default 0, and 0 removes the component exactly), `spheroid_sersic_n` and `spheroid_effective_radius_pc` (default 0, unused while the amplitude is 0), `spheroid_axis_ratio_b` and `spheroid_axis_ratio_c` (default 1, in (0, 1]) and `spheroid_orientation_1_rad`, `spheroid_orientation_2_rad`, `spheroid_orientation_3_rad` (default 0). The irregular family's `body_clump_share` (NOT NULL, default 0, between 0 and 1). One CHECK per number holds the ranges (`chk_galaxy_shape_bar_axis_ratio_b`, `_c`, `chk_galaxy_shape_thick_disk_amplitude`, `chk_galaxy_shape_spheroid_amplitude`, `_sersic_n`, `_effective_radius`, `_axis_ratio_b`, `_axis_ratio_c`, `chk_galaxy_shape_body_clump_share`); a NULL passes, which is how the bar and thick disk stay unset. A galaxy planned before v117 becomes a `spiral` with no spheroid and no clump share, so no stored density changes. The bar and the thick disk were never stored: NULL means "as the structure draw gives them from the seed", which is what every existing galaxy has; a family that has none stores amplitude 0. `get_galaxy_shape` reads them into the `GalaxyShape` fields of the same names (GEN.379; a NULL stays `None` and the structure draw's value stands, a number pins the shape's `structure` too) and `save_galaxy_shape` writes them, every plan replacing what an earlier one stored; a stored spheroid amplitude above 0 is refused until the spheroid component exists (GEN.381). The migration run again leaves them as they are.

### v116 (MAP.172): the nebula totals of the block ladder

Migration 0116 adds three tables. `nebula_block_stats` has one row per block `(m, shell_index, band_index, wedge_index)` of the Galaxy Map's block ladder (`galaxy/drill.py`, `m` 243, 81, 27, 9 and 3) that a nebula's centre is in or its sphere reaches: `nebula_count` (the nebulae whose centre lies in the block, 0 for a block a cloud only reaches), `volume_pc3` (their summed sphere volume), `centroid_x/y/z_pc` (volume weighted) and `bound_radius_pc` (every member's sphere lies within it of the centroid; a catch-up keeps it as an upper bound, a rebuild makes it exact), the volume of each family in `volume_diffuse_pc3`, `volume_emission_pc3`, `volume_reflection_pc3`, `volume_planetary_pc3`, `volume_dark_pc3` and `volume_remnant_pc3`, `dominant_class` (the class with the most volume) and the cover: `cover_volume_pc3`, the volume of the part of each sphere that lies inside the block, summed, and `cover` = `1 - exp(-cover_volume_pc3 / the block's volume)`, the share of the block inside a nebula. A sphere's volume is split among the level-3 blocks by 64, 256 or 1,024 evenly spread points inside it (by its size), and a bigger block's cover is the sum of the level-3 blocks inside it, so the figures add up through the ladder exactly. `nebula_block_classes` holds the volume of each class in a block (by centre), so a nebula can be added to a block without reading the others. `nebula_block_state` has one row, `nebulae_through_id`, the highest `nebulae.id` the two tables have folded in. The nebulae are the rows of `nebula_field` and the stored placed `nebulae` that no `nebula_field` row holds (matched by `uid`). The migration fills the tables from them; `planetgen plan` makes them again with the nebula field, `planetgen plan --fill-nebula-blocks` on its own, and a `galaxy` run adds the nebulae it stored at its end (`db/nebula_blocks.py`). No foreign keys, and none of the three is part of a fingerprint. The version is detected from the table `nebula_block_stats`.

### v115 (DB.44): the sector storage checksum and the parity dirty flag

Migration 0115 adds two columns to `sectors`, at the end of the table: `storage_checksum` (BINARY(32), NULL) and `parity_dirty` (TINYINT NOT NULL, default 1). They are the two columns the rest of DB.9 reads: the checksum is the digest of the sector's stored rows once DB.9's code writes one, and the flag says the sector's parity record is not current, so every sector that existed before the revision counts as not yet protected (default 1). Both are added in one statement with `ALGORITHM=COPY`: `sectors` has a FULLTEXT index (`ft_sectors_name`), and InnoDB refuses an instant or in-place change of such a table (error 1846), so the table is rebuilt once (it holds one row per sector). No row is touched and no code reads the columns yet. A column that is already there is skipped, so a stopped run can be started again. The version is detected from `sectors.parity_dirty`.

The columns are read and written by `planetgen/db/storage_checksum.py` (DB.45; no schema change). The checksum is the SHA-256 of a lossless dump of everything a sector holds as the database stores it: its `sectors` row, the star systems, stars, planets, moons, belts, comets, phenomena and facilities filed in it, and their child rows, every column of every row with `id`, foreign keys and `uid` included, floats by their 64 bits, rows in primary key order. Left out are what the daily orbit update rewrites (`VOLATILE_COLUMNS`: positions, velocities, phases, reflex offsets, binary positions, comet anomalies, `location`, the `galaxy_*_pc` of a system, and each row's clocks), the generated `cell_key` columns, `leaf_digest`, `modified_at` and the two columns themselves; the address-keyed rows (`bright_stars`, `phenomenon_scatter`) and the rebuilt tables are not part of any sector. Every writer of a covered column sets `parity_dirty` to 1 in its own transaction without moving `modified_at`: an INSERT of a sector- or system-filed row through `Connection._flag_insert`, the other writers by name (`tests/test_storage_checksum.py` lists them and fails for one it does not know). `seal` stores the checksum worked out from the rows as they read back and clears the flag; the parity writer (DB.46/48) calls it. `planetgen check-db storage` hashes every sector whose flag is clear against its stored checksum: a mismatch is damage, and the sectors with the flag set, or never sealed, are one warning, not damage.

### v114 (DB.81): the plan stores which shape levers it was told

Migration 0114 adds `galaxy_shape.pinned_stellar_mass_solar` (DOUBLE, nullable, after `halo_min_luminosity_sol`) and the table `galaxy_shape_pins`. A plan writes both with the shape: the column holds the stellar mass the Scale dial fixed (a pinned mass, or the one a pinned size and density hold; NULL when the seed drew it), which `store.get_galaxy_shape` draws the mass chain and structure again at, and the table holds one row per lever the plan was told (`galaxy.dials` lever name, value in the table's units, lengths in parsecs), so the galaxy page can say "mass pinned at X, the rest drawn". A galaxy planned before this has no pins and no pinned mass: everything reads as drawn, and nothing is touched. A step already done is skipped, so a stopped run can be started again. The version is detected from the table `galaxy_shape_pins`.

### v113 (DB.79): the object index

Migration 0113 adds the table `object_index`, one row per object of every kind the site lists. It is empty until the index builder (PERF.140) fills and refreshes it; no other table changes. The primary key is `(kind, uid)`: `kind` is `system`, `star`, `planet`, `moon`, `belt`, `comet`, `facility` or a phenomenon kind (`planetgen.galaxy.objectref.PHENOMENON_KINDS`), `uid` the 96-bit object ID (BINARY(12)). `source_id` is the row id in the object's own table, `name` its name, `sector_id` and `system_id` where it is (NULL for what has none), `parent_uid` the object ID of what it orbits or belongs to, `object_class` its class or type (a star type, a planet class, a phenomenon descriptor), `mass_kg` and `radius_km` (NULL where the kind has none), `placed` (1 when it has a place on the galaxy) and `indexed_at` (when the builder wrote the row). The keyset keys put the kind first, then the sort column, then the object ID: `idx_object_index_name`, `_class`, `_mass` and `_radius`, so a page after a row is one range read; `idx_object_index_source` finds a row from its own table's id, `idx_object_index_sector` and `_system` list what is in a sector or a system, and `ft_object_index_name` is the FULLTEXT name index. There is no foreign key (the builder keeps the table in step, and a key on a table this size costs every write). A table already there is left alone, so a stopped run can be started again. The version is detected from the table.

### GEN.311: a built system stores and loads its orbit tree

No revision: `insert_star_system` writes the tree into the `orbit_nodes` of v104 and `load_star_system` reads it back (`planetgen/db/orbit_trees.py`). A generated binary is one node, number 0, made from the pair it was generated with: `mass_solar` the two stars' mass, `separation_km` the semi-major axis, `eccentricity` (0 for a close pair, the pair's own for a wide one), the plane, `mean_anomaly_deg` the stored phase, `epoch_unix` the database's orbit epoch (the time of the insert before any update has run), `kind` `close` or `wide`, and `interaction`, the Roche classification of the pair at its closest approach (`multistar.classify_pair`: `detached`, `semi-detached`, `contact`, `accreting` or `post-common-envelope`; the classifier the member-tree generator uses). The old `binary_*` columns stay filled and equal the node. Each star row now says where it hangs: `parent_node_id` (0, or NULL for a single star) and its offset from the system's centre of mass (`offset_x/y/z_au`, the stored binary position, zero for a single star). A system that carries its own `orbit_nodes` (a tree drawn by `generation.multistar`) has them stored as they are. A star takes the free slot of its node, in the order of its id (the primary is slot 0). The orbital update still advances the old binary columns and leaves a node's phase alone: a node's place at any time is its mean anomaly advanced from its own epoch (`physics/system_tree.py`).

### v112 (GEN.244): a generated galaxy's young open clusters

Migration 0112 lets `bright_stars.population` take `open_cluster_00` to `open_cluster_11`, the twelve age bins (a quarter dex wide, 1 Myr to 1 Gyr) of the stars of a generated galaxy's young open clusters (`tuning.OPEN_CLUSTER_POPULATIONS`): every star of a cluster shares its age, so the cluster's stars take the bin its age falls in. The clusters themselves are `galaxy_components` rows of kind `cluster` with the `plummer` profile. No row is touched. The version is detected from the check constraint `chk_bright_stars_population` listing `open_cluster_00`.

### v111 (DB.52): the prose is dropped

Migration 0111 drops what v110 and DB.51 made unused: the tables `planet_evolutionary_paragraphs` and `moon_evolutionary_paragraphs`, the columns `flavor_text` and `hab_note` of `planets` and `moons`, and `system_configs.markdown`. The life timeline is the `life_stage`, `evolutionary_pace` and age columns, the flavour line is `flavor_key`, and the habitability note is worked out from the `tier_<domain>` columns, so no information is lost. `flavor_text_count` stays, because the generator still counts the line it draws. A table or column already gone is skipped, so a stopped run can be started again. The version is detected from `planets.flavor_text` being gone.

### v110 (DB.50): structure columns for life and flavour

Migration 0110 adds five nullable columns to `planets` and `moons`, after `flavor_text_count`: `life_stage` (VARCHAR(32), a `generation.evolution.MILESTONE_KEYS` key, NULL before any milestone), `evolutionary_pace` (VARCHAR(16), the timeline scale `fast`, `normal` or `slow`), `life_milestone_age_gyr` and `life_system_age_gyr` (DOUBLE, the age of the milestone and the age of the system the timeline was written at, Gy) and `flavor_key` (VARCHAR(32), which flavour line, replacing `flavor_text`). They replace what is only recoverable now by reading the evolutionary paragraphs and `flavor_text`. Nothing is dropped and no row is touched. Generation writes the columns from DB.51 on (`life.apply_life_data` and `decide_flavor_text`, which draw as before), the loaders and `Planet.evolutionary_data`/`flavor_text` read them back, and the prose columns and tables are no longer written or read; v111 (DB.52) drops them. `flavor_key` is `category:index` into the `tuning` flavour list the line was drawn from (`habitable`, `planet` or `orbital`), so reordering a list changes the lines of stored bodies. `system_configs.markdown` is no longer read or written either. A column that is already there is skipped, so a stopped run can be started again. The version is detected from `planets.life_stage`.

### v109 (GEN.289): a planet or moon is its state vector

Migration 0109 drops `orbital_phase_deg` from `planets` and `moons`. The position, velocity, mass, epoch, parent and dynamical state of v102 hold the body whole: `advance_orbital_phases` moves a due planet or moon with `kepler.universal_step` (a compiled batch step, `physics/state_kernel.py`, written back in one `UPDATE ... JOIN VALUES` per batch), taking `mu` from the stored `distance_km` and `period_years`, and the next due time from the speed the new state has. The phase a page or the generator needs is read from the position (`body_state.argument_of_latitude_deg`). `insert_planet` and `insert_moon` now fill `parent_type` and `parent_id` (a circumbinary planet's parent is its system, and the migration sets that for rows saved before). The system scene carries each body's `position_km`, `velocity_kms` and `epoch_unix` in place of `phase_deg`, and `orbitpositions.js` steps that state with the same algorithm. A column that is already gone is skipped, so a stopped run can be started again. The version is detected from `planets.orbital_phase_deg` being absent.

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

Migration 0102 adds, nullable and backfilled, the columns the layout in `planetgen/db/point_columns.py` says an object lacks, so the code items that fill them (GEN.336, GEN.337, GEN.289) need no revision of their own. A movement vector (`velocity_x/y/z_kms`, galactic axes) goes on black holes, neutron stars, rogue planets, interstellar comets, asteroid fields, nebulae, supernova remnants, asteroid belts and the quasar: a free object takes the tangent of the galaxy's rotation at its place at its old `galactic_orbital_speed_kms`, and the quasar and a belt take a stated zero. A mass (`mass_kg`) goes on interstellar comets, asteroid fields, nebulae, remnants, belts, comets and facilities: a comet nucleus takes its sphere's mass at the density of 67P (533 kg/m3), a nebula or remnant the mass of its bounding sphere at its particle density (an upper bound for a cloud that fills part of its sphere, `physics/bulk_mass.py`), a field or belt a stated zero that the gravity sum skips, and a facility the mass of what it sits on. Belts also get `position_x/y/z_km` (the system frame, at the centre), belts and the quasar get `epoch_unix` and `next_update_due`, and scatter members get `epoch_unix`. Asteroid fields get `formed_at` and `origin_event_id` (GEN.110; no foreign key until the event table, DB.56). Planets and moons get `parent_type` (`star`, `planet` or `system`, with checks), `parent_id` and `dynamical_state` (`bound` or `unbound`), backfilled from `star_id` and `planet_id`. `sectors` gets `max_mass_kg` and `max_mass_id` (the heaviest object and its uid, for GEN.109), NULL until worked out. A star has no columns of its own: it moves with its system, whose point, vector and epoch it reads. The scalar `galactic_orbital_speed_kms` stays until GEN.336's drop revision. Every step is skipped when done, and the version is detected from `sectors.max_mass_id`.

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
- **Companions of a scattered system** (GEN.260): `member_count` and `system_mass_solar` on both scatters (the row's point is the barycentre), and `scatter_members` / `scatter_nodes` (docs/design/systems/multi-star-systems.md sections 5, 13, 14), written by the scatter passes (GEN.262, `planetgen/generation/scatter_companions.py`: a star row always draws its tree, a black hole, neutron star or planetary nebula row at the table's odds, and deleting a scatter row or clearing a scatter deletes its members and nodes): a member's `kind`, `state`, `flags`, masses, seed, and for a young group its offset, velocity and `decay_age_yr`; a node's orbit, `kind` (`close`, `wide`, `compact`, `group`), `interaction` and `flags`; `slot` is SMALLINT on both. No foreign keys (the parent is one of two tables). The nightly update (GEN.265) moves only the row's point, the barycentre; it never writes a member or node row, and a member's place is computed on read from the node phase and epoch (`planetgen/physics/system_tree.py`, `planetgen/db/scatter_systems.py`). GEN.952 (v184): every member also keeps `radius_au` (DOUBLE, nullable, the extent of the member in AU), so the X-ray luminosity of an accreting pair is worked out on read (`planetgen/db/compact_pairs.py`) and not stored; a row made before v184 has NULL and no luminosity, nothing is backfilled.
- **Cell key** (DB.34, DB.35; docs/design/galaxy/spherical-region-queries.md): `cell_key` BIGINT UNSIGNED, a stored generated column holding the Morton key of the object's 0.25 pc cube over plus or minus 262,144 pc (NULL outside), and the covering index `idx_<table>_cell (cell_key, x, y, z)`, on `star_systems`, the eight phenomenon tables, `facilities`, `nebula_field`, `bright_stars` and `phenomenon_scatter`. MySQL computes it from the position columns (`galaxy/cellkey.py` has the same function, tested against it), so it cannot disagree with them. `star_systems` keeps its position relative to its sector, so it also has `galaxy_x/y/z_pc`, its point in the galaxy frame, written by every save and by the orbit update.

### v92 (MAP.205): `galaxy_block_stats` on shells

Migration 0092 makes `galaxy_block_stats` again, for the Galaxy Map's block ladder on shells (`galaxy/drill.py`, `docs/design/galaxy/spherical-block-ladder.md`). It holds one row per block of 243, 81, 27, 9 or 3 shells that holds a generated sector: `m` (SMALLINT UNSIGNED), `shell_index`, `band_index` and `wedge_index` (INT, the block's own address; together the primary key), `sector_count`, `systems` (BIGINT UNSIGNED), `expected_systems` (DOUBLE), `stars` (BIGINT UNSIGNED), `age_sum` (DOUBLE, each sector's mean age times its stars, so a mean is `age_sum / stars`) and `luminosity_sol` (DOUBLE). There is no `ladder_step` column any more: the ladder is 3-ary and a block is named by its shell, band and wedge. The figures come from `sectors` and `sector_stats`; a block with no generated sector has no row. `query.galaxy_block_stats(conn, m, cells)` reads a window's cells in one statement. `db/block_stats.py` keeps the table: a sector saved on its own (a visit) or deleted updates the blocks above it (`drill_chain_of`), and a `galaxy` run updates the blocks above all the sectors it created once at its end. The migration fills the table from the sectors already there, a shell at a time, and skips itself when the table exists, so a stopped run can be started again.

### v91 (DB.32): sectors on spherical shells

The sector grid moved from the cylinder's ring, layer and slot to the sphere's shell, band and wedge (`docs/design/galaxy/spherical-sector-coordinates.md`). Migration 0091 renames `ring_index`, `layer_index` and `ring_slot_index` to `shell_index`, `band_index` and `wedge_index` on `sectors`, `bright_stars`, `sector_stats` and `phenomenon_scatter`, renames `galaxy_shape.outer_ring_index` to `outer_shell_index`, replaces `galaxy_layer` and `galaxy_column` with `galaxy_shell` and drops `galaxy_block_stats` (the block ladder is rebuilt on shells with MAP.205 and MAP.206).

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

`phenomenon_scatter.map_luminosity_sol` (FLOAT, nullable) is the luminosity, in solar luminosities, of a main-sequence star of the row's `mass_solar` by the Galaxy Map's relation (`stellar_evolution.map_luminosity_sol`, GEN.200), so a faint but massive neutron star or black hole is shown where a star of its mass would be. It is NULL for a row with no mass (planetary nebulae, supernova remnants, hypervelocity stars, the quasar). A scatter writes it with each row; `planetgen fill-map-luminosity` fills it, in id-range batches of 200,000 that readers never wait for, in a galaxy scattered before v84 (resumable: it only fills rows still empty). The migration adds the column empty. No index. FLOAT keeps about seven digits and the server prints six, so a value read back agrees with the function to 1e-5 relative.

### v83 (DB.24): `nebula_field`

`nebula_field` holds the galaxy's field clouds of classes M and N (giant molecular clouds and dark clouds; P and Q stay implicit in the seed), one row each, written by `planetgen scatter --phenomena` (or `--nebula-field` alone) from the galaxy seed over the 50 pc field cells (`galaxy/nebula_field.py`). Primary key `(cell_i, cell_j, cell_k, n)`: the field cell and the cloud's place in that cell's draw, so a row is found again from the seed alone. `nebula_class` (CHAR(1)), `radius_ly` and `center_x_pc`, `center_y_pc`, `center_z_pc` (DOUBLE, the exact numbers the seed gives) describe the cloud. `uid` (BINARY(12), nullable, indexed) and `is_stored` (TINYINT, default 0) stay empty until a sector's save stores the cloud as a `nebulae` row, which gives it its object ID: working an ID out for every cloud would cost a nebula shape each (about 0.1 s). The table is seed-fixed, emptied by every plan and rebuilt by every phenomena scatter and by `planetgen scatter --nebula-field`; `planetgen clean-nebula-field --apply` drops the clouds the keep rates now reject. The migration only creates the empty table.

### v82 (PERF.68, PERF.70, PERF.74): indexes and `sector_system_counts` for the pages that timed out

`phenomenon_scatter` gains `idx_phenomenon_scatter_class (kind, subtype, mass_solar)`, which the Galaxy Map's coarse tiles use for the few biggest classes instead of reading an address range (3 million scattered rows: 12.7 s to 0.25 s). `star_systems.idx_star_systems_quadrant` becomes `(quadrant, name, id)` and `idx_star_systems_binary_name (is_binary, name, id)` is new: the Systems list sorted by Octant, Binary or Sector counts group after group and reads the page from the group it lands in (`query._systems_by_group`) instead of sorting every system. `sector_system_counts (sector_id PK, system_count)` stores each sector's system count, with an index on the count; the migration fills it, and `query.refresh_sector_system_counts` (one `GROUP BY` on the sector index, run in the background through `countcache`, not more often than every 30 s) keeps it. A sector with no systems has no row. The Sectors list and its sorts by systems and density read it.

### v81 (ADM.49): `galaxy_shape.arm_level`, `galaxy_shape.core_amplitude`

`galaxy_shape.arm_level` (DOUBLE, NOT NULL, default 1) is the thin disk's level between its arm crest and its inter-arm trough: the arm factor is `arm_level * (1 + arm_amplitude * cos)`, so the crest density is `arm_level * (1 + arm_amplitude)` and the trough `arm_level * (1 - arm_amplitude)`. `galaxy_shape.core_amplitude` (DOUBLE, NOT NULL, default 0) was the galactic core's centre density, a Gaussian added to the bulge; v99 (GEN.280) replaced it with `nuclear_amplitude`. The defaults are the model before v81, so an upgraded galaxy keeps its density. The command line and the Generate page take one arm contrast (GEN.611; `planetgen plan --arm-contrast 0.4`), stored as `arm_amplitude`, and the nuclear amplitude (since v99); a new plan always stores `arm_level` 1.

### v80 (MAP.165): `phenomenon_scatter.mass_solar`

`phenomenon_scatter.mass_solar` (DOUBLE, nullable): the solar mass a scattered black hole (the central one too) or neutron star was drawn with, inside its class range and above the scatter's cut. The sector build gives the object exactly this mass, and the Galaxy Map sizes the point from it (a straight line in log mass from 0.2 at one solar mass to 1.0 at 1e8), so a heavier object is never drawn smaller. NULL for quasars, hypervelocity stars, nebulae and remnants, and on rows scattered before v80 (the map then sizes those by mass class).

### v71 (GEN.167): `galaxy_shape.phenomenon_min_mass_solar`

`galaxy_shape.phenomenon_min_mass_solar` (DOUBLE, nullable): the phenomenon scatter's mass cut. The scatter places only the neutron stars and black holes at or above it (with 20, only the intermediate-mass black holes: about 2.7e5 rows at a 0.1% share, about 800 at the 3e-6 of GEN.259, instead of 1.17e9); a sector fill draws the share below it on the stream `{seed}:phenomena-fill:{ring}:{layer}:{slot}`, so the two add up to the whole population. `planetgen clear-scatter --phenomena` and `planetgen scatter --phenomena` re-scatter at a new cut (`--mass --phenomenon-min-mass M` stores the stellar limit); filled sectors already hold every mass and are left out. See docs/design/systems/phenomenon-scatter-mass-cut.md.

### v72 (GEN.86): stellar activity and magnetic fields

`stars` gains `log_lx_lbol`, `l_xuv_w`, `xuv_saturated`, `flare_n33_per_yr`, `flare_alpha` and `xuv_fluence_j`; `planets` and `moons` gain `magnetic_moment_a_m2`, `dipole_class`, `magnetopause_rp`, `xuv_flux_earth`, `xuv_exposure_index` and `flare_irradiation_index`. They are the inputs of the surface radiation dose (GEN.87) and the habitability index (GEN.83, GEN.89).

### v74 (GEN.88): hydrosphere and ocean chemistry

`planets` and `moons` gain `water_mass_fraction`, `hydrosphere`, `ocean_fraction`, `land_fraction`, `ocean_depth_km`, `ice_shell_km`, `hp_ice_km`, `ocean_class`, `ocean_ph`, `water_activity` and `phosphorus`; `rogue_planets` gains `hp_ice_km` and the `hycean` regime. Rogue oceans are now capped by the liquidus, and their ice lid counts as 1.09 km of ice per km of water with a pressure-lowered base. See `planetgen/physics/hydrosphere.py`.

### v75 (GEN.87): surface radiation dose, UV and the galactic hazard

`planets` and `moons` gain `surface_dose_msv_yr`, `dose_gcr_msv_yr`, `dose_sep_msv_yr`, `dose_ground_msv_yr`, `dose_helio_mult`, `uv_surface_index` and `ozone_loss_flag`; `stars` gains `lethal_event_rate_per_gyr`. (v135 adds `ozone_column_du`, below.) Generated with the planet's spin (it needs the air, water and magnetic field) and updated by `refresh_containment` when a system enters or leaves a nebula or remnant. See `planetgen/physics/radiation.py` and docs/design/habitability/activity-magnetism-radiation-hydrosphere.md section 4.

### v76 (GEN.89): the habitability score

`planets` and `moons` gain `phi4`, `phi4_pressure`, `phi4_temperature`, `phi4_chemistry`, `phi4_radiation`, `tier_pressure`, `tier_temperature`, `tier_chemistry`, `tier_radiation`, `phi_bio`, `phi_cpx`, `phi_tech`, `l_solv`, `l_chem`, `l_ener`, `l_rad`, `equipment_tier`, `hab_note` and `energy_flux_w_m2`, all NULL for a gas giant and on a row generated before v76. Computed last in a body's generation (after its air, water, dose and light) and again by `refresh_containment` when the dose changes. See `planetgen/physics/habitability_world.py` and docs/design/habitability/habitability-index.md section 8.

### v77: `phenomenon_scatter_classes`

`phenomenon_scatter_classes` (`kind`, `subtype`, `placed`, `built`; `subtype` is `''` where a class has none) holds how many scattered phenomena of each class the plan placed and how many a sector's save has built since. The Phenomena table counts and pages the unbuilt ones from it instead of scanning `phenomenon_scatter`, which holds over a hundred million rows in a real galaxy. Written when a phenomena scatter ends, kept by `mark_phenomena_built`, emptied with the scatter.

### v78 (DB.20): object IDs

`uid` is `BINARY(12)` (`BINARY(10)` before v93) and UNIQUE on its own on `star_systems`, `stars`, `planets`, `moons`, `asteroid_belts`, `comets`, the eight phenomenon tables (`black_holes`, `neutron_stars`, `nebulae`, `supernova_remnants`, `quasars`, `rogue_planets`, `interstellar_comets`, `asteroid_fields`) and `facilities`; the old `UNIQUE (star_system_id, uid)` and the 96- and 64-bit hashes are gone. An ID is 96 bits since v93 (80 before): the birth sector's address (shell 19, biased band 18, wedge 19 bits; ring 12, layer 12, slot 16 in the 80-bit form), a serial in that sector (28 bits; the top two say generated, run-time or field-drawn) and a body number (12 bits; 0 for the top-level object, 1 and up for the stars, planets, moons, belts and comets born in a system, one counter across the five kinds). Printed `FE81000A2B-0000005-000`; the layout is `planetgen/galaxy/object_uid.py` (GEN.170). The sector fill gives a sector's systems, then its phenomena, generated serials 0, 1, 2 ... in the order they are inserted (`store._UidIssuer`); a row saved later (an admin-added body, a system saved on its own) takes a run-time serial from the sector's `id_counters` row, or for a body the system's next body number (`store.assign_uids`). `sectors.uid` is unchanged (the designation).

`id_counters (kind, scope, next_value)`, primary key `(kind, scope)`: `kind = 'sector'` with the 5-byte birth sector as `scope` holds the next run-time serial there; `kind = 'system'` with the system's 10-byte ID holds the next body number in it (a system with no row starts past the bodies it holds). A row only moves up, on its own autocommitted connection (`store._reserve_counter`, like `id_blocks`), and `planetgen reset` keeps it, so a number is never given twice. The migration numbers existing rows by row order (a sector's systems, then its phenomena; a system's stars, planets, moons, belts, comets) and gives rows with no sector address run-time serials at the no-sector address (ring 0, layer 0, slot 0). `facilities.uid` is added NULL; GEN.172 fills it.

### v190 (GEN.474, GEN.475): intermediate ground and the chondrite group

Migration 0190 adds `ground_intermediate_fraction` (DOUBLE, nullable, at the end of the table) to `planets` and `moons`: the areal share of the surface that is intermediate (andesitic) ground, the volcanic crust of a body whose mantle melt is 52 to 63 percent silica. Nothing is backfilled: a row made before the revision keeps NULL and the regenerated galaxy replaces it. Migration 0187 also adds `chondrite_group` (VARCHAR(4), nullable, at the end of the table) to `asteroid_belts` and `asteroid_fields`: the meteorite group the body's mean composition is nearest (`CI`, `CM`, `H`, `L` and so on, or `none` far from every group; `physics/meteorite_groups.py`), NULL on a row made before and left out of the sheet. A stopped run can be started again (a column already there is skipped).
