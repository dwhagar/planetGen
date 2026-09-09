# planetGen Roadmap: Full Database Storage + Web Interface

## Vision

Long-term goal: a fully populated galaxy (every sector, every star system,
every planet/moon/asteroid belt, with full generated detail preserved)
stored in a real relational database, eventually served through a web
interface for browsing and search. This is a multi-phase, multi-session
effort — phases below are ordered by dependency, not necessarily by when
they'll be tackled.

## How to use this document

Work top to bottom. Completed work is condensed into short summary bullets
(full historical rationale lives in git history and `CHANGELOG.md`, not
here) — read `src/stellarObjects/schema.sql` and `docs/database-schema.md`
for the actual current schema rather than this document. Trim finished
items down to a pointer whenever a section is revisited; only genuinely
open items need working detail. Phase 4 and parts of Phase 5 are still
open.

## Phases 0-3 — Complete: serialization, database, and CLI persistence

- **Phase 0 (flavor text bug)**: flavor text is decided once at generation
  time instead of being re-rolled at render time.
- **Naming ASCII audit**: every string constant in `stellarObjects/names.py`
  is verified 7-bit-ASCII-printable, with a regression test guarding future
  additions.
- **Phase 1 (object-graph serialization)**: every generated class has a
  `to_dict`/`from_dict` pair (`stellarObjects/serialization.py`);
  reconstruction bypasses `__init__` so loading never re-runs random
  generation.
- **Phase 2 (relational database)**: a normalized SQLite schema
  (`stellarObjects/schema.sql`) and persistence layer (`stellarObjects/_db.py`)
  with versioned migrations (`PRAGMA user_version`) — see
  `docs/database-schema.md` for schema history.
- **Phase 3 (CLI tools use the database)**: `sectorGen.py`/`systemGen.py`
  save every run to the database unconditionally; `src/queryDb.py` is a
  read-only search CLI.

## Completed work log (through 2026-09-07)

Pointer index only — full rationale/detail for each is in `CHANGELOG.md`
and git history.

- Backend API (Flask), mounted at `/api/` alongside the interim `html/`
  browser (commit `ee8daab`); see `docs/api.md`.
- Evolved-star pre-Big-Bang mass fix via reject-and-resample
  (`starData.py`) — CHANGELOG [5.3.1].
- Galaxy coordinate system: design doc
  (`docs/design/galaxy-coordinate-system.md`) with all 8 open questions
  decided, then implemented and merged as Phase 4 Track C (see Phase 4
  below) — CHANGELOG [5.3.5]+.
- Habitable/atmosphere sanity review
  (`docs/analysis/habitability-atmosphere-sanity-review.md`) found Class
  M/P statistically indistinguishable, an inverted greenhouse-factor
  formula, pressure mathematically decoupled from gravity, and
  implausibly-low gas-giant densities — all fixed via the
  physical-plausibility anomaly finder (`stellarObjects/plausibility.py`,
  `src/tests/physical_plausibility_cli.py`) and Track A physics fixes
  (mass-weighted harmonic-mean gas-giant density blend, gravity-coupled
  atmospheric retention, Class P's own `albedo_range`) — CHANGELOG [5.3.5].
- Greenhouse-formula fix + per-class climate tuning for M/O/H/K/L/N/E/F/G/V
  (`greenhouse_multiplier_range`/`atm_molar_density_range`/
  `atm_density_range` overrides, `src/tests/climate_tuning_cli.py`) —
  CHANGELOG [5.3.7].
- Life-bearing classes (Q, W) restricted to the ecosphere zone only, plus a
  description/atmosphere-text cleanup pass across `PLANET_CLASSES` —
  CHANGELOG [5.3.8].
- New-class investigation: gas giants excluded from zones `h`/`e`, two
  physically-impossible small hot-zone rocky classes merged into A/B, two
  redundant brown-dwarf-scale classes cut entirely, per-class `size_mode`
  bell-curve radius distributions added — CHANGELOG [5.3.9]/[5.4.0].
- `sectorGen.py --density`: controllable sector density, shared with
  `galaxyGen.py` — CHANGELOG [5.4.1].
- Class R (unreachable — `h`/`e`/`c` all `False`) removed entirely rather
  than left as dead weight; naming's triple-consonant validation bug
  fixed; single-database navbar link fixed; system-page TOC made
  collapsible (`<details>`/`<summary>`, collapsed by default) — CHANGELOG
  [5.4.2].
- System-page TOC moved out of the inline description flow to a fixed
  right-margin rail (mirroring the left `.sidenav`), shown only above
  `min-width: 90rem` — CHANGELOG [5.4.3].

## Investigate Further

- [ ] Class W's "tidally locked world with extreme temperature variations"
  identity is a day/night split that no per-class range (albedo, molar
  density, greenhouse multiplier, or atmosphere density) can produce from a
  single global `surface_temperature` scalar -- would need an actual
  dayside/nightside model, flagged during the greenhouse-formula/per-class
  climate tuning pass (CHANGELOG.md [5.3.7]) but out of scope for it.
- [ ] Class K (Mars analog) and, by construction, every other ecosphere-zone
  class are generated at the same zone-midpoint orbital distance as Class M
  -- this generator doesn't place different terrestrial classes at
  different distances within (or beyond) the habitable zone the way real
  Mars sits much farther from the Sun than Earth. K's tuned values get as
  close to real Mars' absolute temperature/pressure as achievable under
  that constraint (~231K/~0.57kPa vs real ~210K/~610Pa) but can't fully
  close the gap without a zone/distance-placement change, which is a larger
  design question than per-class range tuning.
- [ ] Classes P and W could still receive the same
  `atm_molar_density_range`/`atm_density_range`/`greenhouse_multiplier_range`
  treatment the M/O/H/K/L/N/E/F/G/V pass (CHANGELOG.md [5.3.7]) gave the
  other habitable classes -- P already has a working `albedo_range` from an
  earlier pass and wasn't part of this round's ordered list; W has the
  day/night structural gap noted above.
- [ ] Render an image or web interface to visualize the location of 2 points in galactic space -- see TODO in src/api/routes.py near `systems_near`.
- [ ] Introducing realistic orbital paths and speeds to all bodies in space, would need a dedicated update script to update like once a month or something to adjust all of the coordinates.
- [ ] Search parameter for searching by not only planet class but planet size, or sort by size in the tagged search field -- see TODO in src/queryDb.py near `process_args` and src/html/search.py near `_planets_panel`.
- [ ] Still open from the file-system cleanup (5.3.2/5.3.3): consider
  moving `src/api/` into `../src/html/` to expose the API endpoint from the same
  served tree. Deliberately not done -- it was posed as an open question,
  not a decision, and moving a Flask package into `../src/html/`'s Apache
  `DocumentRoot` needs its own look at exposure/routing implications
  first. The rest of the file-system cleanup this pointed at is done as
  of 5.3.3 (see CHANGELOG.md [5.3.3]); short pointer note in
  `src/api/__init__.py`.

## File Management

Both done: schema-migration backups are gzip-compressed and excluded from
the web database picker and from a subsequent migration run
(`stellarObjects._db.migrate_database`'s `BACKUP_MARKER`, honored by
`../src/html/lib/dbutil.py` and `../src/migrateDb.py`). See CHANGELOG.md.

## Population and Politics

- [ ] Need to start thinking about assigning government ownership to a particular star system in the database so that together the systems make territories that are mapped out in 3D space by the star systems.
- [ ] Worlds with life on them can be flagged for generated names of dominant races.
- [ ] Probably going to need a space fairing species database.
- [ ] Need to think about under-developed / older civilizations and the differences and how to store and present that data based on society age.

## Phase 4 — Galaxy-scale generation

- [x] **Galaxy coordinate system (Track C) — implemented and merged.**
  Single galaxy per database, `GALACTIC_CENTER_DISTANCE_LY` computed from a
  sector's real position, fixed orientation convention, deterministic
  Fibonacci-sphere shell placement, computed (not stored) adjacency, and a
  radius-based generation-unit primitive covering both whole-shell batch
  generation and a localized neighborhood mode — see
  `docs/design/galaxy-coordinate-system.md` section 8. Implemented across
  the v3->v4 `sectors` schema migration, `src/stellarObjects/galaxyGeometry.py`,
  and `galaxyGen.py` (`--shell K` / `--center-sector ID --radius-pc R`),
  with end-to-end CLI coverage in `src/tests/test_galaxy_gen.py`.
- [x] **Sector cube vertices with neighbor relaxation — implemented and
  merged.** Every galaxy-placed sector's cube is now built as 8 explicit
  vertices (`sectors.vertices_pc`, v6) from Track C's fixed orientation
  convention, then nudged toward nearby sectors' matching corners
  (`stellarObjects/sectorGeometry.relax_vertices`) to shrink -- not
  eliminate, which isn't geometrically possible for a cube tiling of a
  sphere -- the seams between them. See
  `docs/design/galaxy-coordinate-system.md` section 9.
- [ ] **Galaxy thickness / shape (disk-density envelope) — design drafted
  (revision 2), not implemented.** Milky-Way-scale exponential-disk-plus-
  bulge-plus-spiral-arm density model, evaluated and persisted for every
  possible sector up front in a new `galaxy_sector_plan` table (built by a
  new `galaxyPlan.py` batch tool) rather than decided at generation time —
  gated by a deterministic "predicted less than 1 star per sector" cutoff
  that stops the radial scan outward, with each planned sector getting a
  stable sequential index. `galaxyGen.py` consumes plan rows instead of
  computing occupancy itself. A feasibility investigation during this
  design pass measured ~10.5 billion qualifying sectors at real Milky-Way
  scale (~1 TB, several hours to build with pruning + NumPy vectorization,
  down from a naive ~12+ days) — see
  `docs/design/galaxy-disk-density.md` for the full model, schema, and
  benchmarked numbers (the design pass this item and
  `docs/design/galaxy-coordinate-system.md` section 7 question 2 called
  for).

## Phase 5 — Web interface (long-term; needs its own dedicated planning pass)

An interim, dependency-free read-only browser already exists at
[`../src/html/`](html-interface.md) (plain Python CGI scripts, no framework) plus
[`examples/apache/`](apache-deployment.md) (example vhost config +
`set-permissions.sh`), meant for a single-user/small-scale Apache2
deployment today rather than the full multi-user vision below.

- [x] Backend API framework: **Flask**, chosen over FastAPI/Django REST
  Framework — no ORM opinion, deploys via `mod_wsgi` in the same Apache
  process model the interim browser already uses. See
  [`docs/api.md`](api.md).
- [ ] Frontend for browsing/searching the galaxy (sector maps, system detail
  pages, search/filter UI).
- [ ] Almost certainly means moving off SQLite to **MySQL** (the database
  server already present on the deployment host, per user decision — not
  PostgreSQL as this section previously speculated) for real concurrent
  multi-user access. Phase 2 deliberately stuck to raw `sqlite3`/plain SQL
  for now, so this is hand-porting the DDL and persistence layer rather
  than a drop-in config change. Major changes this will require, not yet
  started:
    - **`stellarObjects/schema.sql` porting**: `PRAGMA user_version`
      (the current schema-version mechanism, see `_db.py`'s
      `SCHEMA_VERSION`/`migrate_database`) has no MySQL equivalent — needs
      a real `schema_migrations` tracking table instead. `INTEGER PRIMARY
      KEY` autoincrement semantics, `CHECK` constraint support, and
      `TEXT`/`REAL` column types all differ between SQLite and MySQL and
      need explicit type mapping (`REAL` → `DOUBLE`, `TEXT` → `VARCHAR`/
      `TEXT` per-column, etc.); pick an explicit storage engine (InnoDB,
      for real foreign-key enforcement matching `PRAGMA foreign_keys = ON`
      today).
    - **`stellarObjects/_db.py` porting**: swap the `sqlite3` stdlib
      module for a MySQL driver (`PyMySQL` or `mysqlclient` — pick one and
      justify it, similar to how this section already justified Flask
      over FastAPI); every `?` positional placeholder becomes `%s`; add
      real connection pooling (SQLite's single-file-lock model doesn't
      carry over — MySQL wants a proper pool for concurrent access, which
      is the whole point of this migration).
    - **`src/api/config.py` / `src/queryDb.py` porting**: both currently assume a
      SQLite file path (`DB_PATH`/`--db-path`) rather than a connection
      string/host+credentials pair — see TODO in src/api/config.py near
      `Config.DB_PATH` and src/queryDb.py near `open_readonly`/`--db-path`.
    - **Data migration**: existing `.db` files need a one-time export/
      import into MySQL; no tooling for this exists yet.
    - Revisit tooling (e.g. an ORM/migration framework) at that point if
      the hand-ported approach proves painful, per the original note here.

### Near-term: interim `../src/html/` browser enhancements

Done: distance calculator superseded by stored `star_systems.location`
(schema v3, CHANGELOG [5.3.0]); Search reachable from every page's shared
header; database picker skipped when only one database exists (and, as of
[5.4.2], reachable again via the navbar's `index.py?all=1` link even in
that case); site configuration via `webconfig.json`; wiki-URL reachability
check with a clipboard-copy fallback (commit `ee8daab`).

- [ ] **Sprite-based graphical system view**: render a system's star,
  planets, and moons as small icon sprites sized relative to each other --
  see TODO in src/html/system.py near `_bodies_html`.
- [x] **Sector Map: interactive 3D (drag-to-rotate, scroll-to-zoom)** --
  done. `src/html/lib/starmap.py` now emits a real CSS 3D scene instead of
  a fixed-projection SVG: 6 bordered `<div>` cube faces (the standard
  "CSS 3D cube" recipe) plus one billboarded `<div>` per star, positioned
  via plain layout (`left`/`top`) for x/y and `transform: translateZ()`
  for z. `src/html/static/sectormap.js` tracks two rotation angles and a
  zoom factor from pointer drag / wheel / +/- buttons and feeds them to
  `#starmap-scene`'s CSS transform every frame; the browser's own
  `preserve-3d` compositor handles rotation and occlusion, no manual
  painter's-algorithm re-sort needed. Two things the original plan didn't
  anticipate, both resolved:
    - **No `perspective`.** Billboarding each dot (counter-rotating it so
      it keeps facing the camera instead of going edge-on as the scene
      turns) needs the plain algebraic inverse of the scene's rotation --
      correct only for a pure-rotation (orthographic) composition. With
      `perspective` on `.starmap-stage`, the composition becomes genuinely
      projective and that inverse stops cancelling correctly (confirmed
      directly: it only worked at the one angle it was tested at). Fixed
      by dropping `perspective` entirely -- `preserve-3d` occlusion is
      unaffected either way, and a schematic sector map has no real need
      for vanishing-point foreshortening.
    - **Dot-click detection is resolved by geometry, not native
      hit-testing.** `elementFromPoint`/real click dispatch turned out
      unreliable for an element nested this deep in a rotated
      `preserve-3d` hierarchy (confirmed directly: `elementFromPoint` and
      `elementsFromPoint()[0]` disagreed for the identical coordinate).
      `sectormap.js` instead resolves a click against every dot's own
      `getBoundingClientRect()` (which stays reliable regardless), ties
      going to whichever center is nearest the click.
  See CHANGELOG.md [5.4.4].

### Deployment history (interim `../src/html/` browser)

Several deployment bugs surfaced on first production rollout to the
Ubuntu/Apache2 VPS (`starmap.moltenaether.com`) and are now resolved by
[`install.sh`](../install.sh)/[`update.sh`](../update.sh) (repo root): CGI
scripts deployed non-executable, CRLF line endings, `nltk`'s corpus
download failing under `www-data`'s unwritable home directory, and
`setup.py install` breaking against old apt-provided packaging shadowing
setuptools' own vendored copies. Full root-cause detail is in git history
around `install.sh`/`update.sh`/`examples/apache/set-permissions.sh`.

## Future ideas — not scheduled, just parking so it isn't lost

- [ ] **Subdwarf (Yerkes VI) progenitor mass/age modeling is a known gap**:
  fixing the pre-Big-Bang evolved-star mass issue required excluding
  Yerkes class VI from the new mass-sampling check, because its *entire*
  allowed mass range (0.1-0.8 Msun,
  `physical_constants.YERKES_MASS_CONSTRAINTS["VI"]`) sits below the
  ~0.88 Msun cutoff where a progenitor's own main-sequence lifespan would
  already exceed `UNIVERSE_AGE_GY` — every possible mass would be
  rejected. `Star._calculate_initial_star_age_and_lifespan` still runs
  subdwarfs through the same "derive lifespan from progenitor mass"
  evolved-star logic as giants/supergiants, so a generated subdwarf's age
  can still come out older than the universe. Real subdwarfs (sdB/sdO) are
  thought to form via binary mass-stripping rather than single-star
  post-main-sequence evolution, so the single-star progenitor-lifespan
  model may just be the wrong model for this class entirely — needs its
  own design pass (a different age-generation path for VI, rather than
  another mass-sampling tweak).
