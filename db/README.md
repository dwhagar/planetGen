# planetGen Database Format

This document describes the SQLite database schema defined in
[`stellarObjects/schema.sql`](../stellarObjects/schema.sql). It's the
reference for anyone reading, querying, or extending the database — table by
table, every column's meaning and unit, and the conventions that hold the
schema together.

**Status**: the schema is implemented and being written to.
[`stellarObjects/_db.py`](../stellarObjects/_db.py) (private — leading
underscore, not part of the package's public generation API) writes
already-generated `StarSystem`/`SpaceSector` objects straight into these
tables; `sectorGen.py` calls it automatically on every run. There is no
read path yet — nothing reconstructs a live object from database rows, so
a "list what's stored" tool is still future work (`TODO.md` Phase 3). The
actual `.db` file lives in this `db/` directory (gitignored — see the
repo's `.gitignore`, `*.db`/`*.db-journal`) and is created automatically
the first time something writes to it (default path
`db/planetgen.db`, overridable via `sectorGen.py --db-path`). See
`TODO.md` for the full roadmap.

## Persistence layer

`stellarObjects/_db.py` owns every unit conversion at the point of writing
a value into the database; nothing else in the package imports it, and it
never mutates the generation/physics code's own native units. Its shape:

- `get_connection(db_path=None)` / `save_sector(sector, db_path=None)` —
  the two entry points most callers need. `save_sector` opens (or creates)
  the database, applies the schema if needed (idempotent —
  `CREATE ... IF NOT EXISTS` throughout `schema.sql`), and writes the whole
  sector in one transaction.
- `insert_sector` / `insert_star_system` / `insert_star` / `insert_planet`
  (recurses over moons) / `insert_asteroid_belt` / `insert_system_config` —
  the per-table building blocks, callable individually for a standalone
  system with no sector (`insert_star_system(conn, star_system,
  system_config, sector_id=None, position=None)`).
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
`stellarObjects/starData.py`, `doubleStar.py`, `planetData.py`,
`asteroidData.py`, and `spaceSector.py` keeps its own native unit (km, AU,
or ly) exactly as today. Conversion only happens at the persistence
boundary, once it's built: `stellarObjects/utils.py` provides
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

- `PRAGMA user_version` on the database file — the DDL structure version
  (this schema is version `1`).
- `star_systems.schema_version` (per row) — the version of the serialized
  object-graph shape (from the future Phase 1 `to_dict()`) that produced
  that row. Independent of the DDL version because a JSON export/import
  could bring an older object-graph shape into a newer database file.

### Booleans and tri-state flags

SQLite has no native boolean type. Plain booleans are `INTEGER` `0`/`1`
with a `CHECK` constraint. `SystemConfig`'s tri-state flags (`True`/
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

### `sectors`

One row per generated sector.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `name` | TEXT | NOT NULL | e.g. `"Voranthis Sector"` |
| `edge_mpc` | REAL | NOT NULL | Cube edge length, milliparsecs. Native generator value is `SpaceSector.edge_ly` (light-years). |

### `system_configs`

One row per `SystemConfig` "recipe" — the generation parameters a
`StarSystem` was built from (`stellarObjects/config.py`,
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
| `position_x_mpc`, `position_y_mpc`, `position_z_mpc` | REAL | nullable | Position relative to the sector's cubic center. NULL iff not placed in a sector. |
| `quadrant` | TEXT | nullable, CHECK IN ('I'..'VIII') | The sector octant label derived from the position above (see "Quadrant labeling" below). NULL iff position is NULL. |
| `is_binary` | INTEGER (0/1) | NOT NULL, default 0 | |
| `binary_separation_km` | REAL | nullable | Orbital separation between the two stars. NULL for single-star systems. |
| `binary_type` | TEXT | nullable | e.g. `"Binary (G/K)"`. |
| `binary_temperature_k` | REAL | nullable | Average of the two stars' temperatures. |
| `binary_radius_km` | REAL | nullable | The larger constituent star's radius (used as an approximation). |
| `binary_effective_mass_kg` | REAL | nullable | Sum of both stars' masses. |
| `binary_effective_luminosity_w` | REAL | nullable | Sum of both stars' luminosities. |
| `binary_age_gy`, `binary_lifespan_gy` | REAL | nullable | Max of the two stars' age/lifespan. `binary_lifespan_gy` NULL = infinite (a white-dwarf constituent). |
| `binary_habitable_zone_inner_km`, `_outer_km` | REAL | nullable | Computed from the pair's combined luminosity. |
| `binary_system_perimeter_km` | REAL | nullable | Hill sphere, combined mass. |
| `binary_heliosphere_radius_km` | REAL | nullable | |
| `binary_table_type`, `_mass`, `_lum`, `_hab`, `_separation`, `_loc` | TEXT | nullable | The "Binary System Data" table (`doubleStar.py:158-170`), one column per key. This is the *only* properties table with no owning row elsewhere — `BinaryStarProxy` is never itself stored as a `stars` row (see below). All NULL unless `is_binary`. |
| `system_flavor_text` | TEXT | nullable | Decided once at generation time (Phase 0 fix). |
| `schema_version` | INTEGER | NOT NULL, default 1 | See "Versioning" above. |
| `wikitext_content` | TEXT | nullable | Full rendered page, MediaWiki markup. |
| `markdown_content` | TEXT | nullable | Full rendered page, Markdown. |
| `mediawiki_url` | TEXT | nullable | Where this system's page lives (or should live) on MediaWiki. |
| `wikijs_url` | TEXT | nullable | Where this system's page lives (or should live) on Wiki.js. |
| `created_at` | TEXT | NOT NULL, default `CURRENT_TIMESTAMP` | |

**Quadrant labeling.** `quadrant` reuses the generator's own octant scheme
(`stellarObjects/spaceSector.py`'s `classify_octant`, backed by
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
| `mass_kg`, `radius_km`, `luminosity_w` | REAL | NOT NULL | |
| `temperature_k` | REAL | NOT NULL | |
| `age_gy` | REAL | NOT NULL | |
| `lifespan_gy` | REAL | nullable | **NULL = `float('inf')`** — white dwarfs (`yerkes_class` `'VII'`/`'D'`). Never the JSON `Infinity` token. |
| `habitable_zone_inner_km`, `_outer_km` | REAL | NOT NULL | |
| `system_perimeter_km` | REAL | NOT NULL | Hill sphere relative to the galaxy. |
| `heliosphere_radius_km` | REAL | NOT NULL | |
| `table_type`, `table_radius`, `table_mass`, `table_temp`, `table_lum`, `table_hab`, `table_loc` | TEXT | NOT NULL | The "Star Data" table (`starData.py:488-509`), one column per key. Always present — every constituent star renders its own individual table even inside a binary (in addition to the combined `star_systems.binary_table_*` block). |

### `planets`

One row per planet *or* moon (moons self-reference via `parent_planet_id`
— there is no separate moons table). Covers both terrestrial and gas-giant
bodies (`body_type`).

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `star_system_id` | INTEGER | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | Always the reliable owning link, regardless of `star_id`. |
| `star_id` | INTEGER | FK -> `stars.id`, `ON DELETE SET NULL`, nullable | The specific star this planet orbits, when that's a real stored `stars` row — true for every single-star system. **NULL for a binary's planets**: the generator always builds planets against `self.star`, which for a binary is the `BinaryStarProxy` (never one individual constituent star — there's no S-type/circumbinary distinction in the current generator), and the proxy has no `stars` row to point at. |
| `parent_planet_id` | INTEGER | FK -> `planets.id` (self), `ON DELETE CASCADE`, nullable | Set for moons; NULL for top-level planets/belts. |
| `orbital_index` | INTEGER | NOT NULL | Position in the parent's ordered list (star's planets, or a planet's moons). |
| `is_moon` | INTEGER (0/1) | NOT NULL, default 0 | |
| `body_type` | TEXT | NOT NULL, CHECK IN ('t','g') | Terrestrial or gas giant. Unrelated to `stars.star_type`. |
| `name` | TEXT | NOT NULL | |
| `planet_class` | TEXT | nullable | e.g. `"M"`. |
| `distance_km` | REAL | NOT NULL | From the star (top-level) or the parent planet (moon). |
| `radius_km`, `mass_kg` | REAL | NOT NULL | |
| `volume_km3` | REAL | NOT NULL | Derived, stored not recomputed. |
| `period_years` | REAL | NOT NULL | Derived, stored not recomputed. |
| `zone` | TEXT | nullable, CHECK IN ('h','e','c') | Hot / ecosphere / cold. |
| `description` | TEXT | nullable | |
| `gravity_g` | REAL | nullable | Surface gravity, g's. |
| `surface_temperature_k` | REAL | nullable | |
| `density_g_cm3` | REAL | nullable | |
| `atmosphere` | TEXT | nullable | Composition description. |
| `atm_density`, `atm_molar_density`, `atmospheric_pressure_pa` | REAL | nullable | |
| `composition` | TEXT | nullable | Descriptive string — contrast `asteroid_belt_composition`'s structured (component, concentration) rows. |
| `scale_height_km`, `hill_radius_km`, `min_orbit_distance_km` | REAL | nullable | |
| `habitable_zone_inner_km`, `_outer_km` | REAL | NOT NULL | Copied from the host star at generation time, never recomputed. |
| `life_chemical`, `evolutionary_speed` | TEXT | nullable | Set by life-data generation, if any. |
| `flavor_text` | TEXT | nullable | |
| `flavor_text_count` | INTEGER | NOT NULL, default 0 | |
| `table_class`, `table_distance`, `table_period`, `table_radius`, `table_gravity` | TEXT | `table_distance`/`table_period`/`table_radius` NOT NULL, `table_class`/`table_gravity` nullable | The "Planet Data" (or "Class Data" for moons) table (`planetData.py:291-302`), one column per key — same dict shape for both. |

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

### `asteroid_belts`

One row per asteroid belt. Belts have no properties-dict data table
(`AsteroidBelt.to_paragraph_list()` is prose only), so their searchable
columns capture the facts the prose always states instead.

| Column | Type | Null | Notes |
|---|---|---|---|
| `id` | INTEGER | PK | |
| `star_system_id` | INTEGER | FK -> `star_systems.id`, `ON DELETE CASCADE`, NOT NULL | |
| `orbital_index` | INTEGER | NOT NULL | |
| `distance_km` | REAL | NOT NULL | Average distance from the star. |
| `lower_limit_km`, `upper_limit_km` | REAL | NOT NULL | The belt's inner/outer boundary. |
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

`UNION ALL` across `stars`, `planets`, and `asteroid_belts` (each joined to
`star_systems` for `sector_id`), for "every stellar object in this sector"
queries without hand-writing the union each time:

```sql
SELECT * FROM sector_objects WHERE sector_id = ?;
```

Columns: `object_type` (`'star'`/`'planet'`/`'asteroid_belt'`), `object_id`
(the row's real id in its own table), `star_system_id`, `sector_id`,
`name`, `summary` (a short type-appropriate label — a star's `table_type`,
a planet's class or body type, a belt's density), `orbital_index` (NULL for
stars).

This is a view rather than a fourth physical table specifically to avoid
write-side upkeep: a real table would need to be kept in sync on every
insert/update/delete to the tables it mirrors, or drift out of sync. A view
has no storage and resolves against current data on every query.

## Relationships at a glance

```
sectors 1───* star_systems *───1 system_configs 1───* system_config_slots
                  │
                  ├──1───* stars
                  │
                  ├──1───* planets (self-referential: parent_planet_id for moons)
                  │             │
                  │             ├──1───* planet_evolutionary_paragraphs
                  │             └──1───* planet_reflection_spectrum
                  │
                  └──1───* asteroid_belts 1───* asteroid_belt_composition

planets.star_id ──────────> stars.id   (nullable; NULL for binary systems)
```
