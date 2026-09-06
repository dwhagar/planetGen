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

- [ ] **Planet-level fix**: add `planetLife.decide_flavor_text(planet)` to
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
- [ ] **System-level fix**: add `self.system_flavor_text = None` in
  `StarSystem.__init__`, decided once, right after
  `self.star.adjust_age_for_planets(self.planets)` and before the life-data
  loop (preserves current ordering). Move the identical roll block here
  once. In `__str__`, replace both duplicated roll-and-append blocks with a
  read of `self.system_flavor_text`.
- [ ] Confirmed OK, no action needed: `evolution.py`'s
  `get_evolutionary_timeline`/`planetLife.apply_life_data` are already
  correctly generation-time-only (run exactly once per planet/moon, output
  never touched again) — only the two flavor-text sites have the bug.
- [ ] Add regression tests proving `to_paragraph_list()`/`__str__()` are now
  idempotent — calling twice on the same object gives identical text and
  doesn't further mutate `system_config.system_flavor_count`.

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
JSON-blob-in-a-column hybrid.

- [ ] **Pick tooling** (open decision, needs confirmation before
  committing): recommend SQLAlchemy (ORM) + Alembic (migrations), targeting
  SQLite for now with a documented, low-friction path to PostgreSQL later
  for a real web deployment (SQLAlchemy abstracts the SQL dialect; Alembic
  tracks schema changes). This adds two new dependencies beyond the
  project's current sole dependency (`nltk`) — flag and confirm explicitly.
- [ ] `.gitignore` already has stray `db.sqlite3`/`db.sqlite3-journal`
  entries (leftover boilerplate, not added for this project). Whatever the
  actual database filename ends up being needs its own `.gitignore` entry
  (or a glob covering it, e.g. `*.db`) — the database itself is generated
  data and should never be committed, same as any other build artifact.
- [ ] **Draft schema** (first pass, expect revision once implementation
  starts):
  - `sectors` (id, name, edge_ly)
  - `system_configs` (id, the ~13 `SystemConfig.SERIALIZABLE_FIELDS` as real
    columns)
  - `system_config_slots` (config_id FK, orbit_index, type, planet_class,
    moons) — child table for the variable-length `slots` recipe list
  - `star_systems` (id, sector_id FK nullable, system_config_id FK, name,
    position_x/y/z nullable, binary_separation_au nullable,
    system_flavor_text nullable — from Phase 0, schema_version)
  - `stars` (id, star_system_id FK, role: primary/secondary/single, + the
    ~12 `Star` fields as columns)
  - `planets` (id, star_system_id FK, parent_planet_id FK nullable
    **self-referential** for moons, orbital_index, + the ~34 scalar
    `Planet` fields as columns)
  - `planet_evolutionary_paragraphs` (planet_id FK, position, paragraph) —
    child table for the variable-length narrative list
  - `asteroid_belts` (id, star_system_id FK, orbital_index, distance,
    lower_limit, upper_limit, density)
  - `asteroid_belt_composition` (belt_id FK, position, component,
    concentration) — child table
  - Open question to settle during implementation:
    `reflection_spectrum_visible`/`non_visible` (descriptive string lists
    with no real per-item query value) — own child table for full
    normalization, or a single delimited/TEXT column? Leaning TEXT column
    since nobody will query "find planets with wavelength X" specifically —
    flag for an explicit decision, don't silently pick.
- [ ] Build the persistence layer: functions/classes to write a generated
  `StarSystem`/`SpaceSector` (using Phase 1's serialization) into these
  tables, and to read them back out into live `Star`/`Planet`/etc. objects.
- [ ] Tests: round-trip through the actual database (not just in-memory
  dicts), covering the same cases as Phase 1's tests plus sector-level
  cases (multiple systems, positions, a home system).

## Phase 3 — Migrate the CLI tools to the database

- [ ] `systemGen.py`/`sectorGen.py`: add options to persist generated output
  into the database instead of (or alongside) printing wikitext/markdown.
- [ ] A way to list/query what's already stored (e.g. a new small CLI or
  script: "every G-type system", "everything within 50 ly of Sol", listing
  sectors).

## Phase 4 — Galaxy-scale generation

- [ ] Tooling to generate and store MANY sectors at once, eventually
  covering "every sector of every space in the galaxy" per the long-term
  vision — needs its own design pass once Phases 1-3 exist (a galaxy-scale
  coordinate system beyond a single sector's local cube, how sectors
  tile/connect to each other, batch-generation performance at that scale).

## Phase 5 — Web interface (long-term; needs its own dedicated planning pass)

- [ ] Backend API (framework TBD — Flask/FastAPI are natural fits given the
  Python-native stack) serving the database from Phase 2.
- [ ] Frontend for browsing/searching the galaxy (sector maps, system detail
  pages, search/filter UI).
- [ ] Deployment target (where does this actually run/get hosted?) —
  undecided, needs its own decision.
- [ ] Almost certainly means moving off SQLite to PostgreSQL (or similar)
  for real concurrent multi-user access — exactly why Phase 2 recommends
  SQLAlchemy now, to make that switch mostly a configuration change later.

## Open questions still to resolve

- [ ] **Schema versioning**: add `schema_version`/`SCHEMA_VERSION` to
  `StarSystem.to_dict()` from day one (retrofitting onto files that never
  had it is the harder direction). `from_dict` should raise a clear error
  if a file's version is newer than the code understands; older versions
  are the hook point for future migrations, none needed yet.
- [ ] **`float('inf')` for white-dwarf `Star.lifespan`**: Python's `json`
  round-trips this via the non-standard `Infinity` token, which works
  today but isn't strictly valid JSON and may break future SQLite
  `json_extract()`/`json_valid()` querying. Recommend `Star.to_dict`/
  `from_dict` special-case it (`null`/omitted when infinite, reconstruct
  `float('inf')` on load) before any real data exists in this format.
- [ ] **Secondary-star `system_config` asymmetry**: `StarSystem.__init__`
  (`systemData.py:97-100`) currently `copy.deepcopy`s the config for a
  binary's secondary star (forcing `LARGE_STAR=False`), so today the
  primary/proxy/planets share one `SystemConfig` while the secondary
  privately holds its own deep copy. The `from_dict` design (one shared
  `system_config` passed to both `Star.from_dict` calls) would silently
  collapse that asymmetry on reload. `LARGE_STAR` is only ever consulted
  *during* `generate_star()` (pre-persistence-relevant), so this is very
  likely inconsequential post-generation — but confirm with the user
  rather than assume.
- [ ] `reflection_spectrum_visible`/`non_visible` storage shape in the DB
  (child table vs. TEXT column) — see Phase 2; leaning TEXT column since
  nobody will query "find planets with wavelength X" specifically.
- [ ] Whether to recompute `BinaryStarProxy`'s derived fields on load vs.
  snapshot them redundantly — leaning snapshot (see Phase 1), confirm.
- [ ] SQLAlchemy + Alembic as new dependencies for Phase 2 — confirm before
  adopting (project currently has exactly one dependency, `nltk`).
- [ ] By design, there is no seed-based replay anywhere in this plan —
  generation mixes the unseedable `secrets` module with the seedable
  `random` module (already documented in `spaceSector.py`'s own
  docstring), so results are stored directly rather than via a seed. State
  this explicitly in the new `to_dict`/`from_dict` docstrings so a future
  contributor doesn't try to "optimize" this into seed replay.
- [ ] `examples/*.json` are unaffected by any of this — those are
  `systemGen.py --system-file` recipe *inputs*, a wholly separate format
  from the generated-*result* persistence this plan adds.
