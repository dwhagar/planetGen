# planetGen Roadmap: Full Database Storage + Web Interface

## Vision

Long-term goal: a fully populated galaxy (every sector, every star system,
every planet/moon/asteroid belt, with full generated detail preserved)
stored in a real relational database, eventually served through a web
interface for browsing and search. This is a multi-phase, multi-session
effort — phases below are ordered by dependency, not necessarily by when
they'll be tackled.

## How to use this document

Work top to bottom. Phases 0-3 are complete and condensed below into short
summaries (full historical rationale lives in git history / the commits
that did the work, not here) — read `src/stellarObjects/schema.sql` and
`docs/database-schema.md` for the actual current schema rather than this
document.
Phase 4 and parts of Phase 5 are still open.

## Phases 0-3 — Complete: serialization, database, and CLI persistence

- **Phase 0 (flavor text bug)**: fixed a bug where flavor text was rolled
  (and shared counters mutated) at *render* time instead of generation
  time, so rendering an object twice would double-roll/double-mutate.
  Flavor text is now decided once during generation; rendering is a pure
  read.
- **Naming ASCII audit**: every string constant in `stellarObjects/names.py`
  is now verified 7-bit-ASCII-printable (two non-ASCII `MOON_PREFIXES`
  entries were found and transliterated), with a regression test that
  checks any future addition automatically.
- **Phase 1 (object-graph serialization)**: every generated class
  (`Star`/`BinaryStarProxy`/`Planet`/`AsteroidBelt`/`StarSystem`, plus the
  pre-existing `SystemConfig`/`SpaceSector`) has a `to_dict`/`from_dict`
  pair, built on a shared `fields_to_dict`/`fields_from_dict` helper
  (`stellarObjects/serialization.py`) with each class owning an explicit
  `SERIALIZABLE_FIELDS` allowlist. Reconstruction bypasses `__init__`
  entirely (`object.__new__` + attribute assignment) so loading never
  re-runs random generation. Round-trip tests assert byte-for-byte
  identical rendering and correct shared-object identity (same `star`/
  `system_config` instance across every planet/moon).
- **Phase 2 (relational database)**: a real, normalized SQLite schema
  (`stellarObjects/schema.sql`, plain SQL DDL, no ORM) and persistence
  layer (`stellarObjects/_db.py`) — write path (`insert_sector`/
  `insert_star_system`/etc.), read path (`load_star_system`/`load_sector`),
  and versioned migrations (`PRAGMA user_version`, currently at 3 — see
  `docs/database-schema.md`'s schema history for what changed at each
  version).
  Distances are stored in milliparsecs for sector-scale placement and
  kilometers for everything else (see `schema.sql`'s header comment for
  the full convention). Round-trip tests verify byte-for-byte render
  fidelity through the actual database, not just in-memory dicts.
- **Phase 3 (CLI tools use the database)**: `sectorGen.py`/`systemGen.py`
  save every run to the database unconditionally (`--db-path` to
  override the default `db/planetgen.db`); `src/queryDb.py` is a read-only
  CLI for listing/searching what's stored (`sectors`, `systems`, `near`
  subcommands).

## Session status (2026-09-06/07)

- [x] Backend API (Flask) and the wiki-URL reachability check — merged in
  commit `ee8daab`.
- [x] Evolved-star mass fix, galaxy coordinate system design, and the
  habitable/atmosphere sanity review — all three merged into `main`
  (8592 tests passing). Detail on each:
    - **Evolved-star mass fix** (`stellarObjects/starData.py`): fixes the
      pre-Big-Bang evolved-star mass issue via reject-and-resample. Bumped
      the version to 5.3.1 with a CHANGELOG entry. See "Future ideas"
      below for the one narrow follow-up it surfaced (Yerkes class VI).
    - **Galaxy coordinate system design**
      (`docs/design/galaxy-coordinate-system.md`, proposal only, no code):
      spherical-to-XYZ sector placement, **parsecs** recommended for
      galaxy-scale distance, a shell-based radial tiling scheme
      (Fibonacci-sphere sequence per shell), and a v3->v4 schema migration
      sketch for the `sectors` table. Surfaces one finding needing a
      decision: `physical_constants.GALACTIC_CENTER_DISTANCE_LY` is
      currently a fixed constant assumed by every star's Hill-sphere math,
      which becomes wrong once sectors have real galactic positions. Has 8
      open questions flagged for review before `galaxyGen.py`
      implementation starts (this doc is the prerequisite for Phase 4
      below — still needs a read-through and a decision on those 8
      questions, not yet acted on).
    - **Habitable/atmosphere sanity review**
      (`docs/analysis/habitability-atmosphere-sanity-review.md`, analysis,
      not a fix): generated 300+ samples per class and found real problems
      (not just "the clamps are fine to leave removed" as hoped):
        - Class P is no longer meaningfully colder than Class M or any
          other ecosphere terrestrial class — the disabled clamp was the
          *only* thing giving P its cold identity; it has zero built-in
          bias now.
        - Class M's atmospheric pressure never approaches 1 atm across 600
          samples (mean ~1/3 atm) — contradicts the assumption that the
          fixed pressure formula "lands close to realistic ranges on its
          own."
        - Likely root cause for both: the `greenhouse_factor` formula in
          `planetPhysics.py` (~L402) looks physically inverted — rewards
          an atmosphere for being *far* from CO2's molar density rather
          than for having more CO2 in it.
        - Atmospheric pressure is mathematically independent of gravity
          for every class (proven algebraically and numerically) — the
          `scale_height`/gravity terms cancel exactly, which is why gas
          giants spanning a >1000x gravity range all produce nearly
          identical pressure.
        - Gas-giant density blending (`planetPhysics.py` L308-311) can
          produce implausibly "fluffy" planets (densities as low as 0.026
          g/cm^3) — likely a real bug, fix depends on clarifying what the
          blended ratio should represent.
        - Most ecosphere terrestrial classes are statistically
          indistinguishable in pressure/temperature since they share the
          same `"t"`-keyed parameter ranges regardless of class-specific
          flavor text.
      One trivial doc-comment fix already applied and tested (mislabeled
      `ATMOSPHERE_DENSITY["t"]` range comment). Everything else above is
      **still a recommendation, not yet implemented or decided on** —
      this should probably be resolved before further physics tuning,
      since it affects M/P/gas-giant generation broadly, not just the
      originally-scoped clamps question.
- [x] **Physical-plausibility anomaly finder** — merged (`src/tests/physical_plausibility_cli.py`
  CLI, `stellarObjects/plausibility.py` engine, `src/tests/test_physical_plausibility.py`,
  8669 tests passing). Batch-generates across every (class, zone) pair and
  a broad host-star spectral grid; hard-invariant checks (gravity bounds,
  finiteness/sign) run as always-on tests, statistical outlier detection
  (Tukey's fences, including density now) stays a human-report-only CLI.
  A real run independently **corroborated** the habitability sanity review
  above: Class M/P statistically indistinguishable (277.5K vs 278.4K mean),
  Class M pressure maxes out at 77.7 kPa (never near 101,325 Pa/1 atm),
  gas-giant density as low as 0.027 g/cm^3, and pressure/gravity
  correlation ~= -0.11 (confirms the algebraic cancellation). 0 hard
  invariant violations across 17,294 generated bodies — nothing here rises
  to "generator is crashing/producing nonsense," it's specifically the
  greenhouse/pressure-formula and gas-giant-density concerns already
  flagged above that need a decision.

- [x] **Gas-giant density blend and inverted greenhouse factor** — fixed
  (Track A, partial). `planetPhysics.py`'s gas-giant density
  blend now uses a mass-weighted harmonic mean instead of an arithmetic
  mean over a mass fraction (was producing densities as low as
  0.026 g/cm^3); `plausibility.py`'s mirrored `theoretical_gravity_bounds_g`
  updated in lockstep, docstring corrected from "multilinear" to
  "monotonic in each argument." The greenhouse-factor formula no longer
  rewards atmospheres for being *far* from CO2's molar density — it now
  scales with `atm_molar_density` directly, the only atmosphere-composition
  signal that exists in the data model today.
- [x] **Atmospheric pressure/gravity decoupling and Class M/P
  indistinguishability — Track A completed.** `_atmosphere_retention_factor`
  (linear, normalized at Earth gravity) now scales an effective atmospheric
  density used only in the pressure calculation, reintroducing a real
  gravity/pressure relationship (Spearman correlation > 0.5 on a mixed
  terrestrial/gas-giant sample, vs. ~-0.11 before). Class P now has its own
  `albedo_range` (0.5-0.7, real ice/snow Bond albedo) instead of sharing the
  default (0.12-0.35) with every other terrestrial class, so it's
  meaningfully colder than Class M again from the physics itself rather
  than a post-hoc clamp. 8682 tests passing (10 new in
  `src/tests/test_planet_physics_fixes.py`). See CHANGELOG.md [5.3.5].
  Track A is now fully done.
- [x] **Galaxy coordinate system: the 8 open questions were decided this
  session** (single galaxy per database, fix `GALACTIC_CENTER_DISTANCE_LY`
  now rather than defer, no stored per-sector roll angle, and a
  radius-based "sectors within R of a point" generation primitive that
  covers both whole-shell batch generation and a localized
  "observable-region"/cluster mode around an existing sector), **and
  implementation (Track C) is now complete and merged** — see the Phase 4
  entry below for what was written and what closed out the two
  previously-open gaps.
- [x] **Greenhouse-formula fix and per-class climate tuning (M/O/H/K/L/N/E/F/G/V).**
  The `greenhouse_factor` formula's only lever (`atm_molar_density`, scaled
  by a single shared `CO2_MAX_GREENHOUSE_FACTOR` cap) put every terrestrial
  class's warming in the same narrow band regardless of composition — real
  Earth air's own molar mass already produced `greenhouse_factor ≈ 3.33`
  under the old formula (cap=5), driving Class M to a mean 362K instead of
  ~288K; separately, Mars and Venus's near-identical real molar mass
  (~43.3 vs 43.45 g/mol) despite ~100x different real greenhouse forcing
  meant molar mass alone could never tell a "thin, weak" class from a
  "dense, powerful" one. Fixed by raising `CO2_MAX_GREENHOUSE_FACTOR` to a
  generous safety ceiling (500) and adding a per-class
  `greenhouse_multiplier_range` as the real calibration knob, alongside new
  per-class `atm_molar_density_range` and `atm_density_range` overrides
  (extending the override pattern Class P's `albedo_range` already
  established) — see CHANGELOG.md [5.3.7] for the full list of what was
  tuned and to what real-world/relative targets, and the new
  `src/tests/climate_tuning_cli.py` (human-driven iteration tool) /
  `src/tests/test_climate_tuning.py` (regression suite) this introduced.
  Classes P and W intentionally untouched this pass (P already tuned; W's
  day/night "extreme temperature variations" identity needs a model this
  generator doesn't have, not just range tuning — see "Investigate
  Further" below).
- [x] **Habitable/life-bearing classes restricted to the ecosphere zone,
  and a description/atmosphere-text cleanup pass across `PLANET_CLASSES`.**
  Classes Q and W both carried a `life_chemical` while still being valid
  outside zone `e` (Q in all three zones, W in `h` too) — fixed, and locked
  in going forward by a new regression test
  (`test_life_bearing_classes_are_ecosphere_only`,
  `src/tests/test_planets.py`) covering every `life_chemical`-bearing class,
  not just the ones on `HABITABLE_PLANET_CLASSES`. Separately, several
  classes' `description` text repeated a word the render template already
  supplies on its own, producing broken rendered sentences (e.g. "...with a
  thin atmosphere with an atmosphere of...") — fixed for B, E, K, N, Q, W,
  X, and Y. See CHANGELOG.md [5.3.8].

## Investigate Further

- [x] **New-class investigation, acted on.** The prior session's findings
  (gas giants never valid in zone `h`/`e`; Classes S/T/U's radius ranges
  physically impossible for real sub-stellar objects; Classes X/Y redundant
  with A/B/C) were implemented — see CHANGELOG.md [5.3.9] for the full
  real-science-grounded rework (gas-giant zone flags, S/T/U radius +
  density_range correction, X/Y removal/merge, and the gas-giant
  density-blend bug this surfaced and fixed).
- [ ] Class R ("an ejected, geologically active world") still has `h`/`e`/`c`
  all `False` -- zero probability weight, unreachable outside a manual
  `zone_override` (caught by
  `test_known_issue_class_with_no_valid_zone_is_unreachable`,
  `src/tests/test_planets.py`). A genuinely free-floating/rogue planet (no
  host star at all) is a real, increasingly-studied exoplanet category, but
  doesn't fit this generator's star-centric `h`/`e`/`c` zone model at all --
  "which zone" is the wrong question for an object with no star to be zoned
  relative to. Fixing this is a structural question (a rogue-planet
  generation path independent of `StarSystem`), not a one-line zone-flag
  change -- not attempted in the [5.3.9] rework, which stayed within the
  existing star-centric zone model.
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
  of 5.3.3: markdown consolidated into `docs/` (only `README.md`/
  `LICENSE.md`/`CHANGELOG.md` remain at the repo root), `wsgi.py`/
  `queryDb.py`/`migrateDb.py` moved into `src/`,
  `physicalPlausibility.py` renamed into `src/tests/physical_plausibility_cli.py`,
  `apache/` moved into `examples/apache/`, and `examples/*.json` moved
  into `examples/systems/`. `setup.py` and `pytest.ini` were evaluated
  and can't move into `src/` without breaking (see CHANGELOG.md [5.3.3]
  for why, verified empirically not just assumed). See the short pointer
  note in `src/api/__init__.py`.

## File Management

- [x] **Database backups during schema's should be gzipped and not available on the database picker.**
  `stellarObjects._db.migrate_database` now writes a gzip-compressed backup
  named `<db>.v<sourceVersion>-backup-<timestamp>.db.gz` (the `-backup-`
  marker is `_db.BACKUP_MARKER`) instead of a plain `.db` copy.
  `../src/html/lib/dbutil.py`'s `list_databases`/`resolve_db_path` and
  `../src/migrateDb.py`'s directory scan all explicitly exclude anything
  matching `BACKUP_MARKER`, so a backup is never offered by the web
  database picker.
- [x] **Don't upgrade backups, if a backup database is present from an upgrade don't upgrade the backup.**
  Same `BACKUP_MARKER` exclusion in `../src/migrateDb.py` also keeps a
  backup from ever being handed back into `migrate_database` on a
  subsequent run, so re-running the migration CLI can't re-migrate (and
  further nest-backup) its own prior backup.

## Population and Politics

- [ ] Need to start thinking about assigning government ownership to a particular star system in the database so that together the systems make territories that are mapped out in 3D space by the star systems.
- [ ] Worlds with life on them can be flagged for generated names of dominant races.
- [ ] Probably going to need a space fairing species database.
- [ ] Need to think about under-developed / older civilizations and the differences and how to store and present that data based on society age.

## Phase 4 — Galaxy-scale generation

- [x] **Galaxy coordinate system (Track C) — implemented and merged.** All
  8 open questions in `docs/design/galaxy-coordinate-system.md` are
  decided (single galaxy per database, fix `GALACTIC_CENTER_DISTANCE_LY`
  to use a sector's real position rather than defer it, fixed orientation
  convention with no stored roll angle, deterministic Fibonacci
  placement, uniform sector size per shell, computed rather than stored
  adjacency, and a radius-based generation-unit primitive — see that
  doc's section 8). Rescued from the paused worktree
  (`agent-a36f801e275fb2b71`, cut before Track A's completion) and
  finished, then merged into `main` after Track A's completion (5.3.5)
  and the TODO/FIXME-comment migration — the only real conflict was
  `_db.py` composing this schema's v3->v4 migration with Track B's
  already-merged gzip-compressed backups, which merged cleanly with both
  behaviors intact:
    - The v3->v4 `sectors` schema migration
      (`src/stellarObjects/schema.sql`, `_db.py`'s `_migrate_v3_to_v4`).
    - `src/stellarObjects/galaxyGeometry.py` (shell/Fibonacci-sphere
      tiling + the neighborhood-enumeration primitive).
    - `GALACTIC_CENTER_DISTANCE_LY` threaded per-sector through
      `starData.py`/`doubleStar.py` (falling back to the old constant for
      unplaced/standalone sectors).
    - `galaxyGen.py` (`--shell K` batch mode, `--center-sector ID
      --radius-pc R` local-neighborhood mode).
    - **The two previously-open gaps are now closed**:
      `src/tests/test_galaxy_gen.py` adds real end-to-end coverage for
      both `galaxyGen.py` modes (running the actual CLI entry point
      against a temporary database — correct shell/slot addresses,
      correct stored positions, no duplicate slots, already-occupied
      slots skipped on re-run, the large-shell guard, and the
      no-galaxy-position rejection for `--center-sector`), and a full
      `python -m pytest -q` run is green (8699 passed). One regression
      the rebase itself introduced was found and fixed:
      `test_migrate_v3_to_v4_adds_null_galaxy_columns_to_existing_sectors`
      still tried to open the migration backup as a plain SQLite file,
      but Track B's gzip-compressed-backup change (already on `main`)
      made that backup a `.db.gz` — fixed to decompress first, the same
      way the existing v1->v2 backup test already did.
- [ ] Galaxy thickness / shape (disk-density envelope, and support for
  different overall galaxy shapes — spiral, elliptical, irregular, etc.)
  — still needs its own dedicated design pass, deliberately not tackled
  as part of Track C above. Consolidates the design doc's own
  "disk-density envelope" open question (section 7, question 2) with an
  explicit user request this session to flag it for later.

## Phase 5 — Web interface (long-term; needs its own dedicated planning pass)

An interim, dependency-free read-only browser already exists at
[`../src/html/`](html-interface.md) (plain Python CGI scripts, no framework) plus
[`examples/apache/`](apache-deployment.md) (example vhost config +
`set-permissions.sh`), meant for a single-user/small-scale Apache2
deployment today rather than the full multi-user vision below.

- [x] Backend API framework: **Flask**, chosen over FastAPI/Django REST
  Framework — no ORM opinion (fits the existing raw-`sqlite3` persistence
  layer with zero glue), deploys via `mod_wsgi` in the same Apache process
  model the interim `../src/html/` CGI browser already uses, and mounts at `/api/`
  alongside `../src/html/` for an incremental rollout. FastAPI's async/auto-docs
  advantages don't pay for themselves yet (no separate frontend consuming
  the API, and `sqlite3`'s driver is synchronous regardless of framework);
  revisit if a dedicated frontend makes API-contract docs valuable. A
  read-only scaffold now exists — see [`docs/api.md`](api.md) for
  the endpoints, how to run it, and the `mod_wsgi` deployment story.
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

- [x] **Distance calculator** / **"Nearest N systems within radius R"**:
  superseded by `star_systems.location` (schema v3) — every system placed
  in a sector now stores its sector name plus distance (in ly) to its 3
  nearest neighbors, computed at write time from the existing
  `position_x/y/z_mpc` columns. See `CHANGELOG.md` [5.3.0].
- [x] **Search reachable from every page**: the shared page header
  (`../src/html/lib/page.py`) now includes a Search link whenever a database is
  selected, instead of requiring a detour back through `browse.py`.
- [x] **Skip the database picker when only one database exists**:
  `../src/html/index.py` now redirects straight to `browse.py` if the configured
  database directory contains exactly one `.db` file.
- [x] **Site configuration file**: `webconfig.json` (repo root, gitignored;
  `webconfig.json.example` committed) holds `site_name`/`base_url` today,
  plus unused placeholder fields (`db_username`/`db_password`/`db_name`)
  for a possible future non-SQLite backend. See [`webconfig.md`](webconfig.md).
- [x] **Wiki-URL reachability check + clipboard fallback**: implemented in
  `../src/html/system.py` (merged in commit `ee8daab`), deliberately the bare
  minimum rather than the fuller version this item originally speculated
  about below — a plain HEAD-then-GET check treating any non-200 status,
  timeout, or connection error as unreachable, and no caching layer (a CGI
  script re-executed fresh per request has nothing worth memoizing across
  one page load). Falls back to a "Copy link" button with a clipboard-API
  script when unreachable; renders nothing when the URL field is NULL.
- [ ] **Sprite-based graphical system view**: render a system's star,
  planets, and moons as small icon sprites sized relative to each other --
  see TODO in src/html/system.py near `_bodies_html`.

### Deployment history (interim `../src/html/` browser)

Several deployment bugs surfaced on first production rollout to the
Ubuntu/Apache2 VPS (`starmap.moltenaether.com`) and are now resolved by
[`install.sh`](../install.sh)/[`update.sh`](../update.sh) (repo root):
CGI scripts deployed non-executable (git's `core.fileMode=false` silently
drops the executable bit — both scripts now `chmod +x` unconditionally on
every run, regardless of what mode git stored), CRLF line endings
(`.gitattributes` now pins `../src/html/**/*.py`/`examples/apache/*.sh` to `text eol=lf`),
`nltk`'s corpus download failing under the `www-data` user's unwritable
home directory (`install.sh` now pre-fetches the corpus system-wide;
`names.py` checks `nltk.data.find` before ever attempting a download), and
two rounds of `setup.py install` breaking against old apt-provided
`importlib_metadata`/`packaging` shadowing setuptools' own vendored
copies (resolved by dropping the deprecated `setup.py install` entirely in
favor of a build-isolated `pip install --upgrade --force-reinstall`).
Full root-cause detail for each is in the git history around
`install.sh`/`update.sh`/`examples/apache/set-permissions.sh`.

## Future ideas — not scheduled, just parking so it isn't lost

The physical-plausibility anomaly finder and the original evolved-star
pre-Big-Bang mass-sampling issue moved to "Session status" above once work
on them actually started this session — see that section for current
status instead of here.

- [ ] **Subdwarf (Yerkes VI) progenitor mass/age modeling is a known gap**:
  fixing the pre-Big-Bang evolved-star mass issue (above) required
  excluding Yerkes class VI from the new mass-sampling check, because its
  *entire* allowed mass range (0.1-0.8 Msun,
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