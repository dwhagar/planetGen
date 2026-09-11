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
tables (`sectorGen.py` calls it automatically on every run), and
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
  between (today, always a no-op past the first connection, since every
  MySQL database this project creates already starts at the current
  schema — see "Versioning" below). An existing pre-MySQL-port SQLite
  database is brought in with the separate, one-time
  `src/migrateSqliteToMysql.py` instead (a straight column-preserving copy,
  documented in its own module docstring), not this function.
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
- Every `table_*`/`composition_summary` column is read from a
  `get_table_properties()` method on `Star`/`BinaryStarProxy`/`Planet` (and
  `get_composition_summary()` on `AsteroidBelt`), extracted from each
  class's own `to_paragraph_list()` so the database and the rendered wiki
  page can never drift apart — one formatting implementation, two
  consumers.

## How to read this document

- Every table is listed with its columns, each column's type, nullability,
  and unit (where relevant), and a short note on where the value comes from
  in the generator.
- `FK -> table.column` marks a foreign key.
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

`table_*` columns (see below) are the one exception to both conventions —
they're copies of already-formatted display text (e.g. `"1.2 R☉"`, still
whatever unit the wiki page itself shows: AU, ly, or km), independent of
the raw column's storage unit.

### The searchable-field principle

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
  `9`). Replaces SQLite's `PRAGMA user_version`, which has no MySQL
  equivalent — see `schema.sql`'s "MySQL port" header note.
- `star_systems.schema_version` (per row) — the version of the serialized
  object-graph shape (Phase 1's `to_dict()`) that produced that row.
  Independent of the DDL version because a JSON export/import could bring
  an older object-graph shape into a newer database.

The schema evolved through several versions while still SQLite-backed;
each version's structural change is recorded in `schema.sql`'s own header
comment ("v2" through "v8" notes) rather than duplicated here, since that
file is the one place both the current column list and the historical
rationale for it live together. In brief: v1→v2 split moons out of the
shared `planets` table into their own `moons` table; v2→v3 added
`star_systems.location`; v3→v4 added `sectors`' galaxy-frame placement
columns; v4→v5 dropped every pre-rendered `table_*`/`binary_table_*`
display-string column (superseded by computing display formatting on
demand from the underlying data columns, e.g. `html/lib/tabledisplay.py`);
v6/v7 gave every galaxy-placed sector exact vertices (`sector_vertices`,
built from an exact local spherical Voronoi tessellation among its
same-shell neighbors — see `stellarObjects/sectorGeometry.py`); v8 added
the galaxy-wide density "skeleton" (`galaxy_shape`/`galaxy_shell_band`,
built by `galaxyPlan.py` — see `stellarObjects/galaxyDensity.py`/
`galaxySkeleton.py`) plus a `UNIQUE (shell_index, shell_slot_index)`
constraint on `sectors`, turning a concurrent lazy-generation race
(`galaxyGen.ensure_sector_generated`) into a recoverable `IntegrityError`
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
for the v11 planet/moon position columns, and `_migrate_v11_to_v12` for
the v12 noticeable-motion-interval columns. `migrate_database` applies
whatever steps are needed to reach `SCHEMA_VERSION`, one call `migrateDb.py`
wraps as a CLI (also run automatically by `install.sh`/`update.sh` on
every deploy). A pre-existing SQLite database from before the MySQL port
itself is brought in with the separate, one-time
`src/migrateSqliteToMysql.py` script instead (see its module docstring)
— it only accepts a source already at the database's current
`SCHEMA_VERSION` (today, v12), so a database still on an older SQLite
schema needs a pre-MySQL-port release of this project first.

**This versioning is independent of the control schema's own.** Admin
logins/sessions/API keys/the write-action audit log live in a separate
MySQL schema entirely (`stellarObjects/control_schema.sql`,
`control_schema_migrations`, currently version 1) — see "The control
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

Four tables, versioned independently via `control_schema_migrations`
(currently version 1, mirroring `schema_migrations`'s own shape):

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
  routes (`html/api/authz.audit`) after each one actually succeeds.

`stellarObjects/adminAuth.py` is the only code that reads/writes these
tables directly — `bootstrap_control_schema` creates the schema and seeds
the default `admin`/`password` row (`migrateDb.py` calls this
automatically, alongside its usual content-schema migration), and every
other function there implements one piece of the login/session/API-key/
audit lifecycle `html/api/auth.py`'s routes expose.

### Booleans and tri-state flags

MySQL has no dedicated boolean type either. Plain booleans are `TINYINT(1)`
`0`/`1` with a `CHECK` constraint. `SystemConfig`'s tri-state flags (`True`/
`False`/`None` in Python — force-on / force-off / random) are nullable
`INTEGER` `0`/`1`/`NULL`.

### Rendered wiki text and URLs

`star_systems` holds the complete rendered page for a system twice —
`wikitext_content` (MediaWiki markup) and `markdown_content` (Wiki.js-style
Markdown) — both produced from the *same* generated `StarSystem` object,
rendered back-to-back at save time. They can't be independently
regenerated later and still match, because generation mixes the unseedable
`secrets` module with the seedable `random` module. Alongside them,
`mediawiki_url`/`wikijs_url` record where that page is expected to live (or
does live, once uploaded) on each wiki — one system is one wiki page;
individual stars/planets/moons are sections within that one page, not
separate pages.

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
| `center_x_pc`, `center_y_pc`, `center_z_pc` | DOUBLE | nullable | The sector's center, in a galaxy-frame Cartesian coordinate system whose origin is the galactic center (parsecs — see `docs/design/galaxy-coordinate-system.md`). NULL together iff this sector has never been placed in a galaxy (`sectorGen.py`'s own standalone CLI, or a sector migrated from a pre-v4 database). |
| `galactic_radius_pc` | DOUBLE | nullable | `sqrt(x^2+y^2+z^2)`, persisted (not just derivable) so "sectors within radius R of the core" is a plain indexed range scan — same treatment `star_systems.quadrant` gets. NULL iff the center columns are NULL. |
| `shell_index`, `shell_slot_index` | INTEGER | nullable | This sector's stable address within the shell/Fibonacci-sphere radial tiling scheme (`galaxyGen.py`) — `shell_index` is the radial shell, `shell_slot_index` its placement index within that shell's deterministic ordering. Independently nullable from the center/radius columns above (not part of the same CHECK) — a sector could in principle have a hand-authored galaxy position without this particular placement algorithm's own addressing. |

A `CHECK` constraint enforces `center_x_pc`/`center_y_pc`/`center_z_pc`/
`galactic_radius_pc` being NULL together (see "v3 → v4" in "Schema
history" above for why that addition needed a wholesale table rewrite
rather than an incremental `ALTER TABLE`). This sector's vertices live in
the separate `sector_vertices` table below, present iff this sector has
been placed in a galaxy — neither SQLite nor MySQL can express "rows
exist in another table" as a `CHECK` constraint, so that condition is
enforced at the application level (`stellarObjects._db.insert_sector`)
rather than by the schema itself. A `UNIQUE (shell_index, shell_slot_index)` constraint (v8)
guarantees at most one sector per galaxy address — NULL-together rows
(never placed in a galaxy) don't collide with each other or with a placed
sector, ordinary SQL `NULL` semantics for `UNIQUE`. This is what turns a
lazy-generation race (`galaxyGen.ensure_sector_generated`, see below) into
a clear `IntegrityError` its caller recovers from, instead of a silent
duplicate row at the same address.

### `sector_vertices`

This sector's own exact vertices — one row per vertex, a normalized child
table rather than a JSON column (this schema has no JSON-blob columns
anywhere; see `planet_reflection_spectrum` below for the same treatment
applied to another variable-length list). Built from an exact local
spherical Voronoi tessellation among a sector's same-shell neighbors,
extruded radially between its shell's inner and outer bounding spheres
(`stellarObjects/sectorGeometry.prism_vertices`) — genuinely gap-free
laterally (against same-shell neighbors, not merely reduced — see that
module's own docstring for why exact circumcenters make this possible) and
area-matched (not vertex-matched) radially against the shells in front of
and behind it. No rows exist for a sector never placed in a galaxy.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `sector_id` | BIGINT UNSIGNED | FK -> `sectors.id`, `ON DELETE CASCADE`, NOT NULL | |
| `ring` | VARCHAR(8) | NOT NULL | `'inner'` (on the shell's inner bounding sphere) or `'outer'` (its outer one). |
| `vertex_index` | INT | NOT NULL | This vertex's cyclic position (0-based) within its own ring — pairing `(sector_id, vertex_index)` across the two rings gives the lateral edge each inner/outer vertex pair spans. Not a global ordering across rings. |
| `x_pc`, `y_pc`, `z_pc` | DOUBLE | NOT NULL | Galaxy-frame Cartesian position, parsecs (same frame as `sectors.center_x/y/z_pc`). |

`UNIQUE (sector_id, ring, vertex_index)`; indexed on `sector_id` for the
"every vertex of this sector" query pattern.

### `galaxy_shape`

The galaxy-wide density "skeleton" (v8) — a singleton row (`id` pinned to
`1`) holding everything needed to recompute any sector's exact position
and density on demand, built by `galaxyPlan.py`. Deliberately **not** one
row per sector: a sector's position (`stellarObjects.galaxyGeometry.
sector_position_pc`), density (`stellarObjects.galaxyDensity.
relative_density`), and vertices (`sectorGeometry.prism_vertices`) are all
pure deterministic functions of its `(shell_index, shell_slot_index)`
address plus this handful of galaxy-wide numbers — cheaper to recompute
on demand (sub-millisecond per sector) than to look up, so none of it is
stored per address. What genuinely needs precomputing — where the galaxy
has any content at all — lives in `galaxy_shell_band` below instead.
Building or rebuilding the skeleton replaces this row (and every
`galaxy_shell_band` row) wholesale; there is no partial update, since a
full build is well under a minute even at real Milky-Way scale (the
per-shell work is a closed-form calculation over ~4,100 shells, not a
per-sector scan over billions of candidates).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK, `CHECK (id = 1)` | Pinned to `1` — there is exactly one galaxy. |
| `disk_scale_length_pc`, `disk_scale_height_pc`, `bulge_scale_radius_pc`, `bulge_amplitude`, `pitch_angle_rad`, `arm_amplitude`, `spiral_reference_radius_pc`, `spiral_reference_angle_rad`, `k_norm` | DOUBLE | NOT NULL | `stellarObjects.galaxyDensity.GalaxyShape`'s own fields, verbatim — see that module for what each means and how `k_norm` is calibrated. |
| `arm_count` | INT | NOT NULL | Same source. |
| `edge_pc` | DOUBLE | NOT NULL | The sector edge length this skeleton was built at, parsecs. |
| `expected_system_count_at_density_1` | DOUBLE | NOT NULL | `SpaceSector(edge_ly=...).expected_system_count()` at `relative_density = 1` — cached since every qualification check needs it. |
| `outer_shell_index` | INT | NOT NULL | The last shell index with any qualifying content — this galaxy's real edge, discovered by `galaxyPlan.py` (a run of consecutive empty shells beyond it), not an arbitrary radius. |

### `galaxy_shell_band`

One row per contiguous *candidate* slot-index band per shell (v8) — a
safe, cheap-to-compute superset of where a shell's qualifying sectors
could be (`stellarObjects.galaxySkeleton.find_shell_bands`), not an exact
per-sector list: found via an exact upper bound over every possible spiral-
arm azimuth at a given polar angle, so a slot inside the band isn't
guaranteed to individually qualify, but a slot outside every band is
guaranteed *not* to (the bound can only overstate density, never
understate it). The exact per-slot answer — a single `relative_density`
evaluation — is deferred to the moment that slot is actually visited
(`galaxyGen.ensure_sector_generated`), not computed or stored here.
Almost every shell has exactly one band (this density model is symmetric
about and peaks at the galactic plane for any realistic parameter choice —
verified directly across a full real-Milky-Way-scale build: ~4,100 rows
total, all singletons), but nothing stops a shell from having more than
one (`band_index` orders them), and a shell with no qualifying content at
all simply has no rows here.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | BIGINT UNSIGNED | PK | |
| `shell_index` | INT | NOT NULL | The shell this band belongs to. |
| `band_index` | INT | NOT NULL | 0-based position of this band within its own shell (almost always just `0`). |
| `slot_index_min`, `slot_index_max` | INT | NOT NULL | The candidate slot-index range, inclusive. |

`UNIQUE (shell_index, band_index)`; indexed on `shell_index` for the
"every band of this shell" query `ensure_sector_generated` makes on every
visit to a not-yet-generated address.

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
sudo ../examples/maintenance/install-maintenance-timer.sh [database ...]
```

`updateOrbits.py` mutates rows, so it needs the same read-write database
account `sectorGen.py`/`systemGen.py` use, not `queryDb.py`'s read-only one.

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
| `large_star` | INTEGER (0/1) | tri-state | Force/forbid/random a larger star. |
| `moons` | INTEGER (0/1) | tri-state | Force/forbid/random moon generation. |
| `max_planets` | INTEGER (0/1) | tri-state | Force max vs. min planet count. |
| `planets` | INTEGER (0/1) | tri-state | Force/forbid at least one planet or belt. |
| `star_type` | TEXT | nullable | Explicit spectral type, e.g. `"G2V"`. |
| `name` | TEXT | nullable | Forced system name, if any. |
| `age` | TEXT | nullable, CHECK IN ('young','old') | |
| `intelligent_life` | INTEGER (0/1) | tri-state | |
| `binary_system` | INTEGER (0/1) | tri-state | |
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
| `system_config_id` | INTEGER | FK -> `system_configs.id`, NOT NULL | The recipe this system was generated from. For binaries, this is the *shared* config the primary/proxy/planets all use — the secondary star's transient `LARGE_STAR=False` deep copy (`systemData.py:97-100`) has no field this schema captures, so there's no second config row. |
| `name` | TEXT | NOT NULL | The system's display name — the primary star's name, or `"X Binary System"` for a binary. |
| `position_x_mpc`, `position_y_mpc`, `position_z_mpc` | DOUBLE | nullable | Position relative to the sector's cubic center. NULL iff not placed in a sector. |
| `quadrant` | TEXT | nullable, CHECK IN ('I'..'VIII') | The sector octant label derived from the position above (see "Quadrant labeling" below). NULL iff position is NULL. |
| `location` | TEXT | nullable | Human-readable "sector name + nearest neighbors" summary, e.g. `"Voranthis Kelmoor — nearest: Alpha Prime (4.2 ly), Beta Cerise (7.8 ly), Gamma Ost (9.1 ly)"` — up to 3 neighbors, nearest first, computed once at write time from `SpaceSector.nearest_neighbors` (see "v3" in "Schema history" above). NULL iff position is NULL. |
| `is_binary` | INTEGER (0/1) | NOT NULL, default 0 | |
| `binary_separation_km` | DOUBLE | nullable | Orbital separation between the two stars. NULL for single-star systems. |
| `binary_type` | TEXT | nullable | e.g. `"Binary (G/K)"`. |
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
| `binary_mutual_position_x_km`, `_y_km`, `_z_km` | DOUBLE | nullable | Added in v14. The secondary's Cartesian position relative to the primary — same "position relative to whatever this orbit is around" convention as `planets.position_x/y/z_km` (`utils.orbital_position_au`, applied to `binary_separation_km` and the mutual-orbit orientation columns above). Recomputed by `_db.advance_orbital_phases` in lockstep every time `binary_mutual_orbital_phase_deg` advances. |
| `binary_table_type`, `_mass`, `_lum`, `_hab`, `_separation`, `_loc` | TEXT | nullable | The "Binary System Data" table (`doubleStar.py:158-170`), one column per key. This is the *only* properties table with no owning row elsewhere — `BinaryStarProxy` is never itself stored as a `stars` row (see below). All NULL unless `is_binary`. |
| `system_flavor_text` | TEXT | nullable | Decided once at generation time (Phase 0 fix). |
| `schema_version` | INTEGER | NOT NULL, default 1 | See "Versioning" above. |
| `wikitext_content` | TEXT | nullable | Full rendered page, MediaWiki markup. |
| `markdown_content` | TEXT | nullable | Full rendered page, Markdown. |
| `mediawiki_url` | TEXT | nullable | Where this system's page lives (or should live) on MediaWiki. |
| `wikijs_url` | TEXT | nullable | Where this system's page lives (or should live) on Wiki.js. |
| `created_at` | TEXT | NOT NULL, default `CURRENT_TIMESTAMP` | |

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

### `stars`

One row per individual star: one row for a single-star system, two rows
(`primary`/`secondary`) for a binary. There is **never** a row for the
`BinaryStarProxy` itself — its combined-pair values live on `star_systems`
above (the `binary_*` columns), not here.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `star_system_id` | INTEGER | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | |
| `role` | TEXT | NOT NULL, CHECK IN ('primary','secondary','single') | |
| `name` | TEXT | NOT NULL | |
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

### `planets`

One row per **top-level** planet (as of schema v2 — moons live in their
own `moons` table below, not here; see "Schema history" above). Covers
both terrestrial and gas-giant bodies (`body_type`).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `star_system_id` | INTEGER | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | Always the reliable owning link, regardless of `star_id`. |
| `star_id` | INTEGER | FK -> `stars.id`, `ON DELETE SET NULL`, nullable | The specific star this planet orbits, when that's a real stored `stars` row — true for every single-star system. **NULL for a binary's planets**: the generator always builds planets against `self.star`, which for a binary is the `BinaryStarProxy` (never one individual constituent star — there's no S-type/circumbinary distinction in the current generator), and the proxy has no `stars` row to point at. |
| `orbital_index` | INTEGER | NOT NULL | Position in the star's ordered `planets` list. |
| `body_type` | TEXT | NOT NULL, CHECK IN ('t','g') | Terrestrial or gas giant. Unrelated to `stars.star_type`. |
| `name` | TEXT | NOT NULL | |
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
(`Planet.__init__` only calls `generate_moons` `if not self.is_moon`).

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

### `sector_objects` (view, not a table)

`UNION ALL` across `stars`, `planets`, `moons`, and `asteroid_belts` (each
joined to `star_systems` for `sector_id`), for "every stellar object in
this sector" queries without hand-writing the union each time:

```sql
SELECT * FROM sector_objects WHERE sector_id = ?;
```

Columns: `object_type` (`'star'`/`'planet'`/`'moon'`/`'asteroid_belt'`),
`object_id` (the row's real id in its own table), `star_system_id`,
`sector_id`, `name`, `summary` (a short type-appropriate label — a star's
`table_type`, a planet's/moon's class or body type, a belt's density),
`orbital_index` (NULL for stars).

This is a view rather than a fifth physical table specifically to avoid
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

planets.star_id ──────────> stars.id   (nullable; NULL for binary systems)
moons.star_id   ──────────> stars.id   (nullable; NULL for binary systems)
moons.star_system_id ─────> star_systems.id   (redundant with moons.planet_id's own owner)
```

`galaxy_shape`/`galaxy_shell_band` (v8) stand apart from the tree above —
neither has a foreign key to `sectors` or anything else. They describe the
galaxy as a whole (its shape parameters, and which addresses could hold
content), not any individual sector; `sectors.shell_index`/
`shell_slot_index` is the address both an actual generated sector and a
`galaxy_shell_band` row are independently expressed in, not an FK
relationship.
