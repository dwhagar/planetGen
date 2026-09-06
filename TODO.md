# planetGen Roadmap: Full Database Storage + Web Interface

## Vision

Long-term goal: a fully populated galaxy (every sector, every star system,
every planet/moon/asteroid belt, with full generated detail preserved)
stored in a real relational database, eventually served through a web
interface for browsing and search. This is a multi-phase, multi-session
effort — phases below are ordered by dependency, not necessarily by when
they'll be tackled.

## How to use this document

Work top to bottom — each phase depends on the ones above it (Phase 2's
database needs Phase 1's serialization; Phase 1's serialization needs
Phase 0's bug fixed first, or the frozen values it captures would be wrong).
Before starting **any** implementation, resolve the "Open questions still
to resolve" section at the bottom — several of those are decisions
(dependencies, storage shape for a couple of fields) that change what gets
built, not just how. All file:line references throughout have been checked
against the current source as of the commit that added this document; if
the code has moved since, treat the referenced method/attribute names as
authoritative over the exact line numbers.

## Phase 0 — Prerequisite correctness fix (blocks faithful persistence)

Root cause: `Planet._generate_life_and_flavor_paragraphs`
(`stellarObjects/planetData.py:181-267`) rolls flavor text and mutates
`system_config.system_flavor_count`/`recent_flavor_texts` **every time it's
called**, and it's only ever called from `to_paragraph_list()` — i.e. render
time, not generation time (confirmed: it's called from nowhere else). Likewise
`StarSystem.__str__` (`systemData.py:592-595` and `:626-629`, duplicated in
both the binary and single-star branches) rolls
`random.choice(program_constants.SYSTEM_FLAVOR)`
into a local variable never assigned to `self`. Both must move to
generation time so `to_paragraph_list()`/`__str__()` become pure reads —
required for faithful persistence (there's no attribute to persist
otherwise) and a real correctness bug independent of persistence
(double-rendering currently double-rolls and double-mutates shared
counters).

- [x] **Planet-level fix**: add `planetLife.decide_flavor_text(planet)` to
  `stellarObjects/planetLife.py`, containing exactly the roll/selection
  logic currently at `planetData.py:222-265` (chance gate, recent-list
  filtering/reset, habitability/multicellular check reading
  `planet.evolutionary_data`, `secrets.choice`, setting
  `planet.flavor_text`/`flavor_text_count`, updating
  `system_config.system_flavor_count`/`recent_flavor_texts`) — minus the
  "Sensors show ..." string construction, which stays at render time
  reading the now-already-decided `self.flavor_text`. Must run **after**
  `planetLife.apply_life_data(planet)` for the same planet (the
  habitability check reads `evolutionary_data`, set by `apply_life_data`).
  Call it from `StarSystem.__init__`'s existing life-data loop
  (`systemData.py:243-248`), right after `apply_life_data`, for each planet
  and each moon. Trim `_generate_life_and_flavor_paragraphs` down to just
  appending the sentence from the already-set `self.flavor_text` — no more
  `system_config` mutation at render time.
- [x] **System-level fix**: add `self.system_flavor_text = None` in
  `StarSystem.__init__`, decided once, right after
  `self.star.adjust_age_for_planets(self.planets)` and before the life-data
  loop (preserves current ordering). Move the identical roll block here
  once. In `__str__`, replace both duplicated roll-and-append blocks with a
  read of `self.system_flavor_text`.
- [x] Confirmed OK, no action needed: `evolution.py`'s
  `get_evolutionary_timeline`/`planetLife.apply_life_data` are already
  correctly generation-time-only (run exactly once per planet/moon, output
  never touched again) — only the two flavor-text sites have the bug.
- [x] Add regression tests proving `to_paragraph_list()`/`__str__()` are now
  idempotent — calling twice on the same object gives identical text and
  doesn't further mutate `system_config.system_flavor_count`.

## Naming system — ASCII-7-bit-printable audit (independent of the DB phases)

Every generated name ultimately gets pasted into a wiki page and stored as
`TEXT` in the database — it should never depend on non-ASCII characters
rendering correctly somewhere downstream. `stellarObjects/names.py`'s new
`UNIVERSAL_PHONEMES` pool (added for cross-cultural name-salad diversity,
used by `generate_phoneme_salad_name` for stars/planets/moons/sectors
alike) was built and verified 100% ASCII-7-bit-printable, but the rest of
the file predates that constraint and hasn't been audited against it.

- [ ] **Known violation, not yet fixed**: `MOON_PREFIXES` contains `"Rêv"`
  (French, `ê` = U+00EA) and `"Sueñ"` (Spanish, `ñ` = U+00F1) — both
  outside 7-bit ASCII. Replace with plain-ASCII transliterations (e.g.
  `"Rev"`, `"Suen"`) that preserve the intended flavor.
- [ ] Audit every other list in `stellarObjects/names.py`
  (`STAR_NAMES`/`STAR_PREFIXES`/`STAR_SUFFIXES`,
  `PLANET_NAMES`/`PLANET_PREFIXES`/`PLANET_SUFFIXES`,
  `MOON_NAMES`/`MOON_SUFFIXES`, `SECTOR_NAMES`/`SECTOR_PREFIXES`/
  `SECTOR_SUFFIXES`) for the same problem — none currently enforce it.
- [ ] Add a regression test (e.g. `tests/test_names.py`) that walks every
  public list-of-strings constant in `stellarObjects/names.py` and asserts
  each string is ASCII and printable (`all(32 <= ord(c) < 127 for c in s)`)
  — catches any future addition automatically instead of relying on manual
  review each time a list is extended.

## Phase 1 — Full object-graph serialization

Every generated object needs a canonical, storage-agnostic representation
(feeds the database in Phase 2). Nothing except `SystemConfig` and
`SpaceSector`/`SectorSystemEntry` has any serialization today. Scale:
roughly 300-350 distinct fields for one modest single-star system (1 star +
5 planets + 3 moons + 1 asteroid belt + 1 habitable world), growing
roughly linearly with body count.

- [x] **Serialization mechanism — decided**: a small shared helper module,
  `stellarObjects/serialization.py`, with `fields_to_dict(obj, fields)`
  (`{f: getattr(obj, f) for f in fields}`) and `fields_from_dict(obj, data,
  fields)` (`for f in fields: setattr(obj, f, data.get(f))`), used by each
  class's own `to_dict`/`from_dict`, each owning an explicit
  `SERIALIZABLE_FIELDS` allowlist constant — exactly the pattern
  `SystemConfig` already uses (`config.py`'s `SERIALIZABLE_FIELDS`),
  applied to the other four classes. Rejected: plain `pickle` (not
  readable/diffable/queryable, no versioning story, code-exec risk if
  shared); a blanket `vars(obj)`-minus-exclude-list (fails *unsafe* — a
  future new attribute gets serialized by default until someone remembers
  to blacklist it; an allowlist fails *safe* instead). Reconstruction uses
  `object.__new__(cls)` + attribute assignment, bypassing `__init__`
  (`__init__` re-runs full random generation, not just assignment).
  Caveat: `BinaryStarProxy.mass`/`.luminosity` are read-only properties
  backed by `_effective_mass`/`_effective_luminosity` — its allowlist must
  name the real backing attributes; `from_dict` can't `setattr` the
  property names.
- [x] **Shared back-references — decided**: `star`/`system_config` are
  serialized exactly ONCE, at the `StarSystem` level, and the same live
  instance is threaded down as a parameter into every child's `from_dict`
  call (never stored redundantly in each child's own dict).
- [ ] `Star` (`starData.py`) — `Star.SERIALIZABLE_FIELDS = ["name", "type",
  "yerkes_class", "mass", "radius", "temperature", "luminosity", "age",
  "lifespan", "habitable_zone", "system_perimeter", "heliosphere_radius"]`.
  `to_dict`/`from_dict(cls, data, system_config)` per the pattern above.
  Special-case `lifespan == float('inf')` (white dwarfs) — see risk below.
- [ ] `BinaryStarProxy` (`doubleStar.py`) —
  `SERIALIZABLE_FIELDS = ["name", "type", "temperature", "radius", "age",
  "lifespan", "habitable_zone", "system_perimeter", "heliosphere_radius",
  "_binary_separation_au"]` (mirrors `Star`'s list minus the two
  properties). `to_dict` nests `primary`/`secondary` inline as full
  `Star.to_dict()` calls (never shared outside this one proxy, so no
  dedup needed) plus `effective_mass`/`effective_luminosity`. Store the
  derived values too rather than re-deriving on load — cheap here, and
  avoids having the derivation formulas live in two places.
- [ ] `Planet` (`planetData.py`; covers moons too — nesting confirmed
  exactly 2 levels, a moon's own `.moons` is always empty) —
  `SERIALIZABLE_FIELDS` covering every attribute set in `__init__`
  (`planetData.py:124-179`, ~34 names — re-derive the exact list from the
  current `__init__` body at implementation time rather than trusting a
  recount here) plus `volume`/`period`, excluding `star`/`system_config`
  (back-refs) and `moons` (nested list, serialized recursively as
  `[m.to_dict() for m in self.moons]`). `from_dict(cls, data, star,
  system_config)` threads the same `star`/`system_config` into every moon.
- [ ] `AsteroidBelt` (`asteroidData.py`) —
  `SERIALIZABLE_FIELDS = ["distance", "lower_limit", "upper_limit", "type",
  "density"]`; `composition` handled separately as
  `[list(pair) for pair in self.composition]` on save (JSON has no tuple
  type), `[tuple(pair) for pair in data["composition"]]` on load.
- [ ] `StarSystem` (`systemData.py`) — `to_dict`/`from_dict` live directly
  on `StarSystem` (matching how `SystemConfig`/`SpaceSector` already own
  theirs). `to_dict` emits `schema_version`, `system_config.to_dict()`,
  `star.to_dict()` (polymorphic — includes nested primary/secondary if
  binary), `planets` (list, order = orbital sequence), `system_flavor_text`
  (Phase 0). Deliberately omits `planet_count`/`belt_count`/`moon_count`/
  `hab_count`/`m_count` (recomputed on load via existing `count_objects`/
  `count_habitable`, never trusted from disk) and `stars`/`primary_star`/
  `secondary_star` (resolvable from `star`). `from_dict` is the single
  place that resolves both shared back-references once and re-attaches the
  same instances everywhere (builds `system_config` once; builds `star`
  once via a `class` discriminator key choosing `Star.from_dict` vs.
  `BinaryStarProxy.from_dict`; passes that same `star` into every top-level
  `Planet.from_dict`/`AsteroidBelt.from_dict` call, which passes it into
  every moon too).
- [ ] New tests (`tests/test_serialization.py`, new file, mirroring
  `test_space_sector.py`'s style): round-trip each class and assert every
  field matches; `test_star_from_dict_does_not_rerun_generation`; binary
  round-trip preserves primary/secondary/separation; planet round-trip
  with `HABITABLE_WORLD=True` preserves frozen life data exactly;
  moon-nesting invariant (`is_moon is True`, `moons == []`) verified, not
  assumed; asteroid composition round-trips as tuples
  (`isinstance(..., tuple)`); identity checks (`is`, not `==`) that a
  reloaded system's star is the SAME object across every planet and moon,
  and the SAME `system_config` everywhere; bookkeeping is recomputed, not
  trusted; orbital order preserved; flavor text is already decided before
  first render; rendering is idempotent (call `str(system)` twice, assert
  identical output and unchanged `system_config.system_flavor_count`).
- [ ] Not yet needed: confirmed `evolution.py`'s `get_evolutionary_timeline`
  holds no additional state beyond what's already in
  `Planet.evolutionary_data`.

## Phase 1.5 — Interim JSON extension (can land before the database exists)

- [ ] `SectorSystemEntry.to_dict()` keeps its existing `config` key
  unchanged and adds a new `"generated": self.star_system.to_dict()` key
  (Phase 1's output). `SpaceSector.from_dict()` prefers `generated` when
  present (`StarSystem.from_dict`, faithful) and falls back to today's
  regenerate-from-`config` path when absent — fully backward compatible;
  every existing test in `test_space_sector.py` keeps exercising the
  unchanged fallback branch since old files only ever wrote `config`.
- [ ] Update `spaceSector.py`'s module docstring "Persistence" section
  (currently states reload is *never* byte-identical as a blanket fact) to
  describe the new two-tier behavior.
- [ ] New tests: round-trip with the `generated` key reproduces the exact
  system; `from_dict` falls back to recipe-regeneration when `generated` is
  absent (explicitly delete the key and confirm the old path still works);
  file save/load round-trip preserves the generated system.
- [ ] This is what Phase 2 eventually reads from/writes to instead of flat
  files — a natural stepping stone, not wasted work.

## Phase 2 — Relational database schema + persistence layer

Per explicit decision: a full, real, normalized relational database — not a
JSON-blob-in-a-column hybrid. Scope has grown beyond re-hydrating the raw
generative object graph: the schema also has to support the "publish"
workflow — searching stored systems via DB queries, and uploading a
system's rendered page to a wiki (MediaWiki and/or Wiki.js) by copy/paste or
API. That means storing, per system: its position within its sector, a
queryable inventory of its contents, the exact "finished" property values
that appear in each body's data table (as actually rendered, not
re-derived), the fully rendered page text ready to paste/upload, and the two
URLs the page is expected to live at.

- [x] **Tooling — decided**: raw `sqlite3` (stdlib), plain SQL DDL. No
  SQLAlchemy/ORM, no Alembic — the project has exactly one dependency
  (`nltk`) today and stays that way. `PRAGMA user_version` tracks the DDL
  structure version (see schema-versioning bullet below); any future
  structural change gets a small sequential `_migrate_vN_to_vN+1()` function
  in the persistence module rather than a migration framework.
- [x] `.gitignore` now has a `*.db`/`*.db-journal` glob covering whatever
  the actual database filename ends up being — the database itself is
  generated data and should never be committed, same as any other build
  artifact.
- [x] `db/` scaffolded at the project root (empty for now, and effectively
  untracked by git since it holds nothing but the eventual `*.db` file,
  already ignored) — where the persistence layer will create/open the
  actual SQLite database file once it exists.
- [x] **Distance-unit convention — resolved**: TWO standard units, by
  scale. Sector-scale placement (`star_systems.position_x/y/z`,
  `sectors.edge`) is stored in **milliparsecs** (`_mpc` suffix) — the one
  place kilometers' numbers get unwieldy (an 11.5-ly sector edge is
  ~1.09×10^14 km vs. ~3526 mpc), and sector geometry has no other unit
  competing for consistency the way orbital distances do. Every other
  distance/length-shaped column — orbital distances, habitable zones,
  system perimeter, heliosphere radius, hill radius, scale height, binary
  separation, AND radius/volume — is stored in **kilometers** (`_km`/`_km3`
  suffix), one standard unit instead of the generator's native per-attribute
  mix of km/AU, so search/comparison across system-scale quantities never
  needs unit-aware query logic (this also means no radius/volume
  "exception" is needed at that scale, since those were already km).
  `table_*` columns (the rendered wiki-table snapshot strings) are
  unaffected by either convention — copies of already-formatted display
  text (still AU/ly/km as the wiki shows today), independent of the raw
  column's storage unit. Conversion helpers for the future persistence
  layer live in `stellarObjects/utils.py`: `ly_to_milliparsecs`/
  `milliparsecs_to_ly` for the sector-scale columns; AU-to-km needs no
  helper for the km columns, it's a single multiply by the existing
  `physical_constants.AU_TO_KM`. Generation/physics code is untouched
  either way, confirmed unnecessary by an audit of every distance attribute
  in `starData.py`/`doubleStar.py`/`planetData.py`/`asteroidData.py`/
  `spaceSector.py`/`planetPhysics.py`: each already has exactly one fixed
  native unit, and every consumption boundary already does an explicit
  `physical_constants`-based conversion, so nothing implicit could break.
  Two unrelated pre-existing bugs surfaced during that audit (`Planet.volume`
  silently ending up in m³ despite being documented/computed as km³,
  `planetData.py:192` vs. the overwrite at `:313`; `Planet.scale_height`'s
  formula at `planetPhysics.py:380-382` being dimensionally meters but
  treated as km downstream) are tracked separately (background task
  `task_8913d9bf`), out of scope here.
- [x] **`star_systems.quadrant` — added**: the sector octant label (Roman
  numeral I-VIII, one of the 8 sign(x)/sign(y)/sign(z) combinations — see
  `spaceSector.py`'s `classify_octant`/`program_constants.SECTOR_OCTANT_LABELS`;
  the code's own naming calls these "quadrants" despite being a 3D octant
  scheme), stored alongside the raw `position_x/y/z_mpc` rather than only
  derived on read — it's purely a function of position, but worth
  persisting so it's directly queryable/indexable ("every system in
  Quadrant III") without a UDF or generated-column expression. NULL iff
  position is NULL (not placed in a sector).
- [x] **Schema — written**: `stellarObjects/schema.sql`, verified to load
  cleanly via `sqlite3.executescript()` and column lists cross-checked
  directly against the real `__init__` bodies (`config.py`, `starData.py`,
  `doubleStar.py`, `planetData.py`, `asteroidData.py`) rather than the
  estimates below (e.g. `Planet` turned out to have 27 scalar columns, not
  ~34 — the estimate had folded in the `moons`/`evolutionary_data` lists,
  which are child tables, not columns):
  - `sectors` (id, name, edge_mpc)
  - `system_configs` (id, the ~13 `SystemConfig.SERIALIZABLE_FIELDS` as real
    columns — tri-state fields as nullable INTEGER 0/1/NULL, `markdown` as
    non-nullable INTEGER 0/1)
  - `system_config_slots` (id, config_id FK, orbit_index, type, planet_class,
    moons) — child table for the variable-length `slots` recipe list
  - `star_systems` (id, sector_id FK nullable [standalone systems allowed],
    system_config_id FK, name, position_x_mpc/position_y_mpc/position_z_mpc
    nullable [NULL iff not placed in a sector], quadrant nullable [Roman
    numeral I-VIII octant label derived from position, see the resolved
    bullet above], is_binary, binary_separation_km nullable, plus one column
    per `BinaryStarProxy`
    derived field — binary_type, binary_temperature_k, binary_radius_km,
    binary_effective_mass_kg, binary_effective_luminosity_w, binary_age_gy,
    binary_lifespan_gy [nullable = infinite], binary_habitable_zone_inner_km/
    _outer_km, binary_system_perimeter_km, binary_heliosphere_radius_km —
    all nullable, populated only for binaries, stored rather than re-derived
    since `_effective_mass`/`_effective_luminosity` are computed once at
    generation time; plus one column per key of the "Binary System Data"
    table (`doubleStar.py:158-170`) — binary_table_type, binary_table_mass,
    binary_table_lum, binary_table_hab, binary_table_separation,
    binary_table_loc, all TEXT nullable, populated only for binaries [the
    one properties table with no owning row elsewhere, since
    `BinaryStarProxy` is never itself stored as a `stars` row];
    system_flavor_text nullable [Phase 0]; schema_version; wikitext_content
    nullable; markdown_content nullable [both renderings of the *same*
    generated system — see note below]; mediawiki_url nullable; wikijs_url
    nullable [one MediaWiki + one Wiki.js URL per system page — a system is
    one wiki page, stars/planets/moons are sections within it, per
    `StarSystem.__str__`]; created_at)
  - `stars` (id, star_system_id FK, role: primary/secondary/single, + the
    ~12 `Star` fields as columns [lifespan_gy nullable = infinite, white
    dwarfs], plus one column per key of the "Star Data" table
    (`starData.py:488-509`) — table_type, table_radius, table_mass,
    table_temp, table_lum, table_hab, table_loc, all TEXT, always present
    [every constituent `Star` still renders its own individual table even
    inside a binary])
  - `planets` (id, star_system_id FK, star_id FK→stars nullable
    [the specific star this planet orbits, when that's a real stored
    `stars` row — always true for single-star systems; NULL for a
    binary's planets, since `systemData.py` always generates planets
    against `self.star`, which for binaries is the `BinaryStarProxy`
    (never one individual constituent star — there's no S-type/
    circumbinary choice in the current generator), and the proxy is
    deliberately not stored as its own `stars` row. `star_system_id`
    remains the reliable owning link regardless], parent_planet_id FK
    nullable **self-referential** for moons, orbital_index, + the 27
    scalar `Planet` fields as columns [excluding reflection_spectrum_visible/
    _non_visible — see child table below], plus one column per key of the
    "Planet Data"/"Class Data" table (`planetData.py:291-302`) —
    table_class, table_distance, table_period, table_radius, table_gravity,
    all TEXT [same dict shape for planets and moons])
  - `planet_evolutionary_paragraphs` (planet_id FK, position, paragraph) —
    child table for the variable-length narrative list
  - `planet_reflection_spectrum` (id, planet_id FK, spectrum_type TEXT CHECK
    IN ('visible','non_visible'), position INTEGER, value TEXT) — child
    table replacing the `reflection_spectrum_visible`/`non_visible` scalar
    fields, same position+value shape as `planet_evolutionary_paragraphs`/
    `asteroid_belt_composition` below (**resolved**, see note below: no JSON
    column, this is a real normalized table)
  - `asteroid_belts` (id, star_system_id FK, orbital_index, distance,
    lower_limit, upper_limit, density, composition_summary) — belts have no
    properties-dict data table, prose only (`asteroidData.py:93-141`), so
    their searchable columns are instead the facts that prose always
    states: density, the distance range (lower/upper limit), and a
    composition summary string built the same way as the prose sentence
    (`asteroidData.py:119-137`) — kept alongside the structured
    per-component `asteroid_belt_composition` child table below, not
    instead of it, for queries that need one specific component
  - `sector_objects` — a VIEW, not a table, `UNION ALL`-ing `stars`,
    `planets`, and `asteroid_belts` (each joined up to `star_systems` for
    `sector_id`) into one `(object_type, object_id, star_system_id,
    sector_id, name, summary, orbital_index)` shape, for "every stellar
    object in this sector" queries. Deliberately a view, not a fourth
    physical table: a real table would need to be kept in sync on every
    write to the tables it mirrors (or drift), while a view has no storage
    and resolves against current data on every query — the standard
    relational answer to "I need to search across several typed tables
    uniformly" is a view/union query over normalized per-type tables, not
    an EAV table or a duplicated denormalized one.
  - `asteroid_belt_composition` (belt_id FK, position, component,
    concentration) — child table
  - "Contents of the system" (a queryable inventory of what's in it) is
    deliberately *not* a separate table — the `stars`/`planets`/
    `asteroid_belts` rows above, each carrying `star_system_id` +
    name/type/class, already serve that (e.g. "find all systems with a
    habitable planet" is a query across `planets`, no extra table needed).
  - Design principle behind every `table_*` column: `properties_to_string`
    (`stellarObjects/utils.py:480-524`) builds each data table from a plain
    dict of already-formatted value strings *before* branching on
    `system_config.MARKDOWN` — the dict itself is identical between the
    wikitext and Markdown renderings, only the header labels differ. Since
    each dict's key set is small and fixed per class (Star: 7 keys, Binary:
    6, Planet/Class: 5), it maps directly onto named typed columns rather
    than a JSON blob or a generic key-value table — real columns, fully
    queryable, no `json_extract()` needed, and matching this document's own
    "not a JSON-blob-in-a-column hybrid" decision at the top of this phase.
    They hold the "finished" property values exactly as published,
    independent of rendering format, and don't need to be re-derived from
    the raw scalar columns (whose formatting logic could change across
    versions). Variable-*length* data (the reflection spectrum lists, the
    evolutionary paragraphs, the asteroid composition) is the one case that
    still needs a child table rather than a fixed column set — same
    reasoning, different shape.
- [ ] **Two renderings of one system**: because generation mixes unseedable
  `secrets` with seedable `random` (`spaceSector.py`'s own module
  docstring), the *other* format can never be regenerated later from a
  `StarSystem` and match the first (flavor text, etc. would differ). Both
  `wikitext_content` and `markdown_content` must be produced from the
  **same already-generated** `StarSystem` object, rendered twice back-to-back
  at save time (toggle `system_config.MARKDOWN`, call `str(system)`, toggle
  back, call again) — never regenerated independently after the fact.
- [x] **Schema versioning — resolved**: two independent version numbers,
  neither tied to `_version.py`'s software semver (which tracks CLI/package
  releases, not data shape). `star_systems.schema_version` is row-level —
  tracks the shape of the serialized object graph (Phase 1's `to_dict()`
  output) that produced that row, so a JSON export/import can carry an
  older shape into a newer database file; `from_dict`/DB-read raises a clear
  error if a row's version is newer than the code understands. SQLite's own
  `PRAGMA user_version` is DDL-level — tracks the table/column structure
  itself, checked before applying any future `_migrate_vN_to_vN+1()` step.
  Both start at `1`.
- [x] **Write path — built**: `stellarObjects/_db.py` (leading underscore —
  private, not part of the package's public generation API). Deviates from
  this phase's original plan in one deliberate way: rather than routing
  through a separate Phase 1 `to_dict()`/`from_dict()` object-graph
  serialization module (still not built — see Phase 1's checklist above,
  still all unchecked), `_db.py` reads live `StarSystem`/`Star`/
  `BinaryStarProxy`/`Planet`/`AsteroidBelt` attributes directly and writes
  them straight into the tables above in one transaction
  (`insert_sector`/`insert_star_system`/`insert_star`/`insert_planet`/
  `insert_asteroid_belt`, plus a `save_sector(sector, db_path=None)`
  convenience wrapper). Every `table_*`/`composition_summary` column reads
  from a new `get_table_properties()` method added to `Star`,
  `BinaryStarProxy`, and `Planet` (and `get_composition_summary()` on
  `AsteroidBelt`) — each extracted from that class's existing
  `to_paragraph_list()` so the exact same formatting logic backs both the
  rendered wiki text and the database, with no duplication. Both
  `wikitext_content`/`markdown_content` are rendered from the same
  already-generated object, toggling `system_config.MARKDOWN` and
  restoring it afterward, exactly as planned above. `schema.sql` itself
  was made idempotent (every `CREATE TABLE`/`INDEX`/`VIEW` now says
  `IF NOT EXISTS`) so `_db.get_connection` can safely apply it on every
  connection rather than needing a first-run/migration check.
- [ ] **Read path — not built**: nothing yet reconstructs a live
  `StarSystem`/`SpaceSector` from database rows. This is where Phase 1's
  `to_dict()`/`from_dict()` allowlists (or an equivalent direct
  row-to-object mapping in `_db.py`) still matter — needed for Phase 3's
  "list/query what's already stored" tooling to return anything richer
  than raw rows, and for any future re-upload/re-render workflow.
- [ ] Tests: round-trip through the actual database (not just in-memory
  dicts) — single star, binary, multi-moon, asteroid-belt, and
  sector-with-position cases — verifying both text columns are populated
  and mutually consistent, `lifespan_gy IS NULL` round-trips to
  `float('inf')`, and `PRAGMA foreign_key_check` is clean. Partially done:
  `tests/test_db_persistence.py` now covers the multi-moon/asteroid-belt
  save path end-to-end (added alongside the v2 schema change below, since
  hand-testing that change by running `sectorGen.py` and eyeballing the
  result is exactly what missed a genuine bug — an `INSERT`'s column list
  and value tuple silently drifting one apart in `insert_moon` — that a
  test would have caught immediately). Binary systems, `lifespan_gy`
  round-tripping, and `PRAGMA foreign_key_check` are still untested.
- [x] **Schema v2 — moons split into their own `moons` table**: a `Class D`
  search tag was silently mixing planets and moons together, because both
  shared one `planets` table discriminated only by `is_moon` — a search
  facet querying `planets` had no way to mean "planet" without also
  matching every `is_moon=1` row of the same class. Fixed the same way
  `star_systems`/`planets` already get their own tables: moons now live
  in `moons` (`planet_id` FK to the planet they orbit, plus their own
  `moon_evolutionary_paragraphs`/`moon_reflection_spectrum` child tables),
  not self-referencing via `planets.parent_planet_id`/`is_moon`. No
  self-reference on `moons` itself — confirmed via
  `tests/test_moons.py::test_moons_cannot_themselves_have_moons` that a
  moon never generates its own moons, since `Planet.__init__` only calls
  `generate_moons` `if not self.is_moon`. `PRAGMA user_version` bumped to
  `2`; `stellarObjects/_db.py`'s `migrate_database` converts an existing
  v1 database in place (backup first), and `migrateDb.py` runs it over
  every `*.db` in a directory — wired into `install.sh` (and so
  `update.sh`, which calls it) so a deployed database keeps working after
  a pull brings in the new schema. `html/search.py`'s tag facets and name
  search follow the split: separate "Planet Class"/"Moon Class" (etc.)
  tag groups and a separate Moon name field, instead of one set covering
  both.

## Phase 3 — Migrate the CLI tools to the database

- [x] `sectorGen.py`: every run now builds a `SpaceSector` (via
  `SpaceSector.add_system`, so each system gets a real Hill-sphere-placed
  position) alongside its existing text generation, and saves it to the
  database via `_db.save_sector` — unconditionally, not behind an opt-in
  flag, since populating the database *is* the point now. `--db-path`
  overrides the default `db/planetgen.db` location. Existing stdout/
  `--output` text behavior is unchanged.
- [ ] `systemGen.py`: still text-only, no database option. A standalone
  (non-sector) system has no natural `sector_id`/position, but
  `_db.insert_star_system` already accepts `sector_id=None, position=None`
  for exactly this case — wiring it in is mechanical once wanted.
- [ ] A way to list/query what's already stored (e.g. a new small CLI or
  script: "every G-type system", "everything within 50 ly of Sol", listing
  sectors) — blocked on the read path above for anything beyond raw SQL.

## Phase 4 — Galaxy-scale generation

- [ ] Tooling to generate and store MANY sectors at once, eventually
  covering "every sector of every space in the galaxy" per the long-term
  vision — needs its own design pass once Phases 1-3 exist (a galaxy-scale
  coordinate system beyond a single sector's local cube, how sectors
  tile/connect to each other, batch-generation performance at that scale).
- [ ] `galaxyGen.py` — a new top-level CLI script (sibling to
  `sectorGen.py`/`systemGen.py`) that generates and persists many sectors
  as one galaxy, reusing `sectorGen.py`'s own per-sector generation/save
  logic for each one (the same way `sectorGen.py` itself reuses
  `systemGen.py`'s per-system logic). Explicitly gated on sector
  generation "working properly" first — at minimum, the Phase 2 read path
  and Phase 3's list/query tooling (both still unbuilt above) should land
  before this, since a galaxy-scale tool multiplies whatever gaps exist in
  the single-sector path rather than surfacing them any more gently.
  Needs the galaxy-scale coordinate system this bullet already calls out
  (how sectors tile/connect) decided first.

## Phase 5 — Web interface (long-term; needs its own dedicated planning pass)

**Superseded in part**: an interim, dependency-free read-only browser
already exists at [`html/`](html/README.md) (plain Python CGI scripts, no
framework) plus [`apache/`](apache/README.md) (example vhost config +
`set-permissions.sh`), meant for a single-user/small-scale Apache2
deployment today rather than the full multi-user vision below. The bullets
below split between near-term enhancements to that interim tool and the
still-undecided long-term Flask/FastAPI+Postgres rebuild.

- [ ] Backend API (framework TBD — Flask/FastAPI are natural fits given the
  Python-native stack) serving the database from Phase 2.
- [ ] Frontend for browsing/searching the galaxy (sector maps, system detail
  pages, search/filter UI).
- [x] Deployment target (where does this actually run/get hosted?) —
  undecided, needs its own decision. Web server will be the same one the
  mediawiki site is hosted on.
- [ ] Almost certainly means moving off SQLite to PostgreSQL (or similar)
  for real concurrent multi-user access — Phase 2 deliberately stuck to raw
  `sqlite3`/plain SQL for now, so this move will mean hand-porting the DDL
  and persistence layer rather than a drop-in config change; revisit
  tooling (e.g. an ORM) at that point if the port proves painful.

### Near-term: interim `html/` browser enhancements

- [ ] **Wiki-URL reachability check + clipboard fallback**, on
  `html/system.py`: for each of `star_systems.mediawiki_url`/`wikijs_url`,
  check (server-side, short timeout, some form of caching so a page load
  doesn't stall on two live HTTP round-trips every time) whether the URL
  is set *and* actually resolves; if so, render it as a normal clickable
  link. If it's unset or unreachable, show a "Copy to clipboard" button
  next to the matching content instead (`wikitext_content` for
  `mediawiki_url`, `markdown_content` for `wikijs_url`) — needs a small
  inline JS snippet (`navigator.clipboard.writeText`) since the current
  plain `<textarea>` select-all-and-copy affordance is manual. Requires
  deciding how "does it exist" is actually checked (HTTP HEAD vs. GET,
  what counts as a hit for a wiki that 200s its own "page doesn't exist
  yet" placeholder) and whether the check result is cached anywhere or
  re-checked on every view.
- [ ] **Sprite-based graphical system view**: render a system's star,
  planets, and moons as small icon sprites sized relative to each other
  (from `radius_km`) for an at-a-glance size comparison, the way real
  solar-system scale charts do. Open design questions: where the sprite
  art comes from (hand-authored set keyed by `body_type`/`planet_class`/
  star spectral type vs. some procedural generation), linear vs.
  logarithmic size scaling (a gas giant vs. a moon differ by 2-3 orders of
  magnitude in radius — linear scaling would render most bodies as
  invisible dots), and whether this is a simple side-by-side comparison
  row (answers "size comparison" directly, much simpler) or a full
  scaled-orbit diagram (answers more, substantially harder — distances and
  radii can't share one linear scale without one becoming unreadable).
- [ ] **Distance calculator**: given two systems (same sector, since
  `position_x/y/z_mpc` is stored relative to that sector's own center —
  see `db/README.md`), compute the straight-line distance between their
  stored positions and display it in both light-years and milliparsecs
  (`stellarObjects.utils.milliparsecs_to_ly` already exists for the
  conversion).
- [ ] **"Nearest N systems within radius R"**, surfaced from
  `html/sector.py`: given a chosen system (or the sector's own center) plus
  a radius (accepted in either ly or mpc, converted via the same helper as
  above) and a count N, return the N closest other systems in that sector
  within the radius — a straightforward Euclidean-distance query over the
  existing `position_x/y/z_mpc` columns, no schema change needed. Doesn't
  extend across sectors — there's no cross-sector coordinate system yet
  (see Phase 4's galaxy-scale coordinate bullet above).

### Deployment bugs found in production (interim `html/` browser)

Found live on the deployed Ubuntu/Apache2 VPS (`starmap.moltenaether.com`)
shortly after first rollout — root-caused against the actual server
logs/config, not just reproduced locally. **Resolved** by adding
[`install.sh`](../install.sh) at the repo root as the answer to the open
design question this section originally ended with: `setup.py` stays
scoped to the Python package only (`stellarObjects` + the
`sectorgen`/`systemgen` entry points, installable on any OS); a separate,
explicit `install.sh` (Linux-only, run once per deploy as root) owns
everything Apache/CGI-specific, calling `apache/set-permissions.sh` as
one of its steps.

- [x] **CGI scripts deployed non-executable — root cause found and
  fixed**: Apache logged `AH01215: (13)Permission denied: exec of
  '/var/lib/planetGen/html/index.py' failed` for every request. Root
  cause: this repo's local (Windows) git config had `core.fileMode=false`,
  which silently drops any `chmod +x` before it reaches a commit —
  confirmed via `git ls-files -s html/*.py apache/*.sh`, which showed all
  of them committed as mode `100644` despite being `chmod +x`-ed in the
  working tree when authored. Fixed by no longer trusting git to carry
  the executable bit at all: `install.sh` unconditionally
  `chmod +x`-es `html/*.py` and `apache/set-permissions.sh` on every run,
  independent of whatever mode git stored.
  - **Follow-up bug, also fixed**: that first fix was itself incomplete
    for anything under a subdirectory — `install.sh`'s own
    `chmod +x` line used `find "$HTML_DIR" -maxdepth 1 -name '*.py'`,
    which only ever reached the top-level scripts, silently skipping
    `html/lib/*.py`. `apache/set-permissions.sh`'s own pass had no such
    depth limit and should have caught it regardless, but it only did
    `chgrp` (group ownership), not `chown` (user *and* group) — reported
    as "still doesn't work" for ensuring `html/lib/*.py` was executable
    by `www-data:www-data` specifically. Fixed both: `install.sh`'s
    `chmod +x` line dropped `-maxdepth 1` entirely, and
    `set-permissions.sh` now `chown -R`s (not just `chgrp -R`s) to the
    detected `user:group`, and prints a count of how many `.py` files it
    touched so a wrong path is obvious immediately instead of silently
    matching nothing.
- [ ] **CRLF line endings on the deployed scripts** — not yet fixed. The
  same production session also had to run `sed -i 's/\r$//' index.py` to
  get a clean run. The blobs currently committed at `HEAD` are confirmed
  LF-only, so this was very likely a one-time artifact of exactly when/how
  that particular commit was made from the Windows authoring machine
  (`core.autocrlf=true` locally) — but nothing today actually *prevents* a
  future Windows-side commit from reintroducing CRLF. Add a
  `.gitattributes` entry pinning `html/**/*.py` and `apache/*.sh` to
  `text eol=lf` so this is enforced by the repo itself rather than by
  convention. (Left open — not part of `install.sh`, since it's a repo
  hygiene fix, not a deploy-time one.)
- [x] **`nltk.download('words', quiet=True)` fails under the Apache
  worker user — root cause found and fixed**: `stellarObjects/names.py`
  called this unconditionally at import time. `nltk`'s `download()`
  always targets its *own* default download directory (resolved against
  the current user's home dir) and always attempts `os.makedirs()` there
  first, regardless of whether the corpus already exists elsewhere on
  `nltk.data.path` — confirmed by reading `Downloader._download_package`
  directly. Run as `www-data` under `mod_cgi` (Debian/Ubuntu default home
  `/var/www`), that resolved to `/var/www/nltk_data`, which doesn't exist
  and `www-data` can't create (`PermissionError: [Errno 13] Permission
  denied: '/var/www/nltk_data'`, confirmed from `planetgen.error.log`) —
  breaking every single page of the web interface, since every CGI
  request re-imports the package fresh and re-triggers the same check. It
  had "worked" for the CLI tools only because those happened to run as
  `root` (writable home directory), masking the bug until the web
  interface hit it as a different user. Fixed two ways together (verified
  end-to-end: pre-fetch to a directory on `nltk.data.path`, then confirm
  `nltk.data.find` locates it without touching `download()` at all):
  - `install.sh` runs `python3 -m nltk.downloader -d
    /usr/local/share/nltk_data words` (one of nltk's own default
    system-wide search paths, so no code needs to know the location
    explicitly) and `chmod -R a+rX` on it — every consumer (CLI as any
    user, or the web interface as `www-data`) now reads the same
    pre-populated, read-only copy.
  - `stellarObjects/names.py` now calls `nltk.data.find('corpora/words')`
    first, only falling back to `download()` on `LookupError` — removing
    the unconditional download attempt (and its per-CGI-request
    filesystem overhead) entirely once the corpus is pre-fetched.
- [x] **`apache/set-permissions.sh`: no safe fallback when apache2 isn't
  detectable — fixed**: previously, if `/etc/apache2/envvars` wasn't
  readable *and* no running `apache2` process was found,
  `detect_apache_group` failed the whole script. That fired during
  exactly the install-time window `install.sh` now runs it in — apache2
  may be installed but not yet started/enabled at that point. Changed to
  warn (to stderr) and default to `www-data:www-data` (the standard
  Debian/Ubuntu identity) instead of hard-failing.
- [x] **`update.sh` added** (repo root, alongside `install.sh`): a plain
  `git pull` on a deployment isn't enough on its own — pulling a changed
  file rewrites it with whatever mode is tracked in the repo, silently
  re-breaking the executable-bit fix above for any `.py` file that
  happened to change upstream. `update.sh` refuses to run over
  uncommitted local changes (checks `git status --porcelain` first, so a
  hand-edited deployment isn't silently clobbered), pulls with
  `--ff-only` (fails loudly rather than creating a surprise merge commit
  if history has diverged), and then re-runs `install.sh` unconditionally
  so permissions are guaranteed correct again regardless of what changed.
  Verified all four paths (clean fast-forward, already-up-to-date, dirty
  working tree, diverged history) against a local throwaway git
  remote/clone pair.
- [x] **`install.sh` step 1 (`pip install --upgrade setuptools`) broke on
  the deployed Ubuntu 20.04/Python 3.8 VPS — root cause found and fixed**:
  `sudo ./update.sh` pulled the latest commit fine, but its re-run of
  `install.sh` crashed on `python3 setup.py install` with `AttributeError:
  module 'importlib_metadata' has no attribute 'EntryPoints'`, right after
  pip reported `Successfully installed setuptools-75.3.4`. Root cause:
  setuptools' own `extern` shim (used to reach its bundled, newer copy of
  `importlib_metadata`) prefers a *real*, already-installed
  `importlib_metadata` package over that bundled copy when one is present
  on the import path — and Ubuntu 20.04's apt-provided
  `python3-importlib-metadata` (~1.5.0) is old enough to be missing the
  `EntryPoints` class current setuptools calls, so it silently shadows the
  working vendored copy instead of falling through to it. This wasn't a
  bug in this repo's own code — `pip install --upgrade setuptools` with no
  version pin just started resolving to a newer PyPI release than it used
  to, on a system whose apt-managed `importlib_metadata` never gets
  upgraded, so the exact same command could re-break again in the future
  as PyPI's "latest setuptools" keeps advancing. Fixed by having
  `install.sh` upgrade `importlib_metadata` in the same `pip install`
  call as `setuptools` — installing a current version ahead of the stale
  apt one on the import path resolves the shadowing without needing to
  pin `setuptools` to some specific past version.
- [x] **`update.sh`/`install.sh`: executable-bit fix only ever covered
  `html/*.py`, not `.sh` scripts themselves — fixed**: hit on the same
  session as the `install.sh` permission-denied incident above (`sudo
  ./update.sh` failing on its own `install.sh` invocation with
  "Permission denied" until a manual `chmod +x install.sh`). Root cause:
  the original executable-bit fix (see the `install.sh` bullet earlier in
  this section) only ever re-`chmod +x`ed `html/*.py` plus one
  hardcoded `apache/set-permissions.sh` path — it never covered `.sh`
  scripts in general, including `install.sh` and `update.sh` themselves,
  which are just as exposed to git dropping the executable bit on a pull
  (the same `core.fileMode=false` root cause, not specific to `.py`
  files). Worse, `update.sh` invokes `install.sh` directly
  (`"$SCRIPT_DIR/install.sh"`, not `bash install.sh`), so if a pull
  dropped *its* executable bit, install.sh's own fix never got a chance
  to run at all — the shell refuses to exec a non-executable file before
  any of its own logic runs. Fixed two ways together: `install.sh`'s step
  3 now also does `find "$SCRIPT_DIR" -name '*.sh' -exec chmod +x {} +`
  (every `.sh` in the repo, not a hardcoded list, so a future script
  added anywhere doesn't need this updated to be covered); `update.sh`
  now runs that same `find`/`chmod` itself, immediately before invoking
  `install.sh`, so a dropped executable bit on `install.sh` (or on
  `update.sh` for its *next* run) is fixed before the exec attempt that
  would otherwise fail. One remaining wrinkle, inherent to a
  self-updating script and not fixable further: the *specific* `update.sh`
  process that pulls this fix for the first time is still running the old
  (pre-fix) script text it already had open when `git pull` replaced the
  file on disk out from under it, so that one migration run can still
  need a manual `chmod +x install.sh`. Every run after that first one
  starts by reading the fixed `update.sh` from disk, so the self-heal
  applies from the very start of its own execution and no further manual
  intervention should be needed.
- [x] **`install.sh`: `setup.py install` broke a *second* time, same root
  cause as the `importlib_metadata` incident above, this time for
  `packaging` — fixed by dropping the legacy `setup.py install` path
  entirely**: right after the `importlib_metadata` fix landed, the next
  deploy hit `TypeError: canonicalize_version() got an unexpected keyword
  argument 'strip_trailing_zero'`, from the same `setup.py install`
  invocation. Same shape of bug, different dependency: setuptools'
  `extern` vendoring shim prefers a real, already-installed copy of a
  dependency it vendors over its own newer bundled copy whenever a real
  one is importable, and this time it was Ubuntu 20.04's old apt-provided
  `packaging` doing the shadowing instead of `importlib_metadata`.
  Patching this dependency-by-dependency (add `packaging` to the `pip
  install --upgrade` line, wait for the next one to break, repeat) was
  recognized as a whack-a-mole fix, not a real one — there was no reason
  to believe `packaging` and `importlib_metadata` were the only two
  vendored dependencies old enough system packages could shadow. Root
  fix: stop running the deprecated `setup.py install` at all (setuptools
  itself prints "Please avoid running setup.py directly" for this exact
  invocation) and switch to `pip install --upgrade --force-reinstall
  "$SCRIPT_DIR"` (a proper PEP 517 build-isolated install) instead. Build
  isolation builds the package in a throwaway environment that can't see
  this system's site-packages at all — only the stdlib and pip's own
  freshly fetched build dependencies — so the shadowing can't happen
  there for *any* vendored dependency, present or future, closing the
  whole bug class rather than the one instance of it that happened to
  surface. `--force-reinstall` was added deliberately (a plain `pip
  install .` skips reinstalling when it thinks the same version is
  already installed, which would silently stop `update.sh` from
  redeploying freshly pulled source between version bumps in
  `stellarObjects/_version.py` — the same "stale copy shadows the fresh
  one" failure mode this whole fix is otherwise about avoiding). The
  explicit `pip install --upgrade setuptools importlib_metadata` line
  from the previous fix was removed as no longer needed: this system's
  global `setuptools` no longer needs upgrading at all for this step,
  since the isolated build environment resolves its own build
  dependencies independent of whatever is already installed globally.

## Open questions still to resolve

Resolved since the last pass (kept here only as a pointer, full reasoning
lives where the decision is used): schema versioning, `float('inf')`
white-dwarf lifespan storage, and `reflection_spectrum_visible`/
`non_visible` storage shape (a normalized `planet_reflection_spectrum` child
table, not a TEXT/JSON column — no JSON-blob columns anywhere in this
schema, per the pure-relational-DB decision) — see the Phase 2 schema
section above. `SQLAlchemy`/Alembic vs. raw `sqlite3` — resolved in favor of
raw `sqlite3` (Phase 2 tooling bullet above).

- [ ] **Secondary-star `system_config` asymmetry**: `StarSystem.__init__`
  (`systemData.py:97-100`) currently `copy.deepcopy`s the config for a
  binary's secondary star (forcing `LARGE_STAR=False`), so today the
  primary/proxy/planets share one `SystemConfig` while the secondary
  privately holds its own deep copy. The `from_dict` design (one shared
  `system_config` passed to both `Star.from_dict` calls), and the Phase 2
  schema's single `star_systems.system_config_id`, would both silently
  collapse that asymmetry on reload/persist. `LARGE_STAR` is only ever
  consulted *during* `generate_star()` (pre-persistence-relevant), so this
  is very likely inconsequential post-generation — but confirm with the
  user rather than assume before Phase 1 locks in the single-shared-config
  design.
- [ ] Whether to recompute `BinaryStarProxy`'s derived fields on load vs.
  snapshot them redundantly — leaning snapshot (see Phase 1 and the Phase 2
  `star_systems.binary_*` columns), confirm.
- [ ] By design, there is no seed-based replay anywhere in this plan —
  generation mixes the unseedable `secrets` module with the seedable
  `random` module (already documented in `spaceSector.py`'s own
  docstring), so results are stored directly rather than via a seed. State
  this explicitly in the new `to_dict`/`from_dict` docstrings so a future
  contributor doesn't try to "optimize" this into seed replay.
- [ ] `examples/*.json` are unaffected by any of this — those are
  `systemGen.py --system-file` recipe *inputs*, a wholly separate format
  from the generated-*result* persistence this plan adds.
