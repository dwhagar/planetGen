# Changelog

## [5.3.5] - 2026-09-07

### Added
- **Galaxy-scale coordinate system (Track C), merged.** Sectors can now
  be placed on a galaxy-wide, radial shell/Fibonacci-sphere tiling
  (`docs/design/galaxy-coordinate-system.md` sections 0-8) instead of
  existing only in isolation:
  - **Schema v3 -> v4** (`stellarObjects/schema.sql`): `sectors` gains six
    nullable galaxy-frame columns (`center_x/y/z_pc`, `galactic_radius_pc`,
    `shell_index`, `shell_slot_index`), NULL together for a sector never
    placed in a galaxy. `_db.migrate_database`'s new `_migrate_v3_to_v4`
    handles the upgrade (existing `_migrate_v1_to_v2`/`_migrate_v2_to_v3`
    updated in lockstep, per that function's "each hop maps straight to
    the current schema" design); `docs/database-schema.md`'s schema
    history documents the change.
  - **`stellarObjects/galaxyGeometry.py`** (new): the shell/Fibonacci-sphere
    tiling primitives (`shell_sector_count`, `shell_radius_pc`,
    `sector_position_pc`) and `enumerate_sectors_within_radius` — an
    exact, two-prune generation-unit primitive that finds every sector
    address within a radius of an arbitrary galaxy-space point without
    ever scanning a whole shell's slot count (design doc section 8).
  - **`galaxyGen.py`** (new, repo root): a CLI generating many sectors as
    one galaxy, reusing `sectorGen.py`'s own per-sector generation/save
    path. `--shell K` batch-generates a whole radial shell (guarded by
    `LARGE_SHELL_WARNING_THRESHOLD`, needing `--limit`/`--yes` above it);
    `--center-sector ID --radius-pc R` generates a local neighborhood
    around an already galaxy-placed sector. Either mode skips slots a
    sector already occupies.
  - **`GALACTIC_CENTER_DISTANCE_LY` is now per-sector**, not a single
    fixed constant: `Star.calculate_system_perimeter`/
    `BinaryStarProxy._calculate_system_perimeter_static` accept a
    `galactic_center_dist_ly` override, threaded from `galaxyGen.py`
    through `StarSystem`/`Star`/`BinaryStarProxy`, falling back to the
    old fixed constant for unplaced/standalone sectors
    (`sectorGen.py`'s own CLI unaffected).
  - `stellarObjects/utils.py` gains `mpc_to_pc`/`pc_to_mpc` (exact) and
    `pc_to_ly`/`ly_to_pc` (display-string conversions), per the design
    doc's unit-choice section.
  - New end-to-end coverage: `src/tests/test_galaxy_gen.py` runs
    `galaxyGen.py`'s actual CLI entry point (both `--shell` batch mode
    and `--center-sector` local-neighborhood mode) against a real
    temporary database, asserting correct shell/slot addresses, correct
    stored positions, no duplicate slots, and that already-occupied
    slots are skipped on re-run — the gap `docs/TODO.md`'s Phase 4 entry
    flagged as not yet done when this was paused mid-session. Rebased
    onto current `main` (past Track A's physics fixes and Track B's
    gzip-compressed migration backups) with no functional changes beyond
    a `test_db_migration.py` assertion that needed to decompress the
    (now gzip-compressed) v3->v4 migration backup the same way the
    existing v1->v2 backup test already did.

## [5.3.4] - 2026-09-07

### Fixed
- **Gas-giant density blend** (`stellarObjects/planetPhysics.py`): the
  core/atmosphere blend used a mass fraction as an arithmetic-mean weight
  between two densities, which is dimensionally wrong and could produce
  gas giants as low as 0.026 g/cm^3. Replaced with the mass-weighted
  harmonic mean, the physically correct way to combine two component
  densities via a mass fraction. `plausibility.py`'s independently
  reimplemented copy of this formula (used to derive analytical
  hard-invariant gravity bounds) was updated in lockstep, including its
  docstring's justification for corner-evaluation (still valid: the new
  formula is monotonic in each argument, just no longer multilinear).
- **Inverted greenhouse factor** (`stellarObjects/planetPhysics.py`): the
  formula rewarded an atmosphere for being *far* from CO2's molar density
  rather than for actually containing more CO2 — backwards from physical
  reality. Now scales directly with `atm_molar_density`, the only
  atmosphere-composition signal the data model has today.

### Changed
- **Database schema-migration backups are now gzip-compressed** and
  excluded from the web database picker and from a subsequent migration
  run (previously a plain `.db` copy that a naive `*.db` glob would both
  surface in the picker and silently re-migrate on the next run). See
  `docs/TODO.md`'s "File Management" section for detail.

### In progress, not yet merged (see `docs/TODO.md` for exact state)
- Atmospheric pressure is still algebraically independent of gravity, and
  Class M/Class P remain statistically indistinguishable — both scoped
  and partially started, paused mid-session in worktree
  `agent-a8acb02b98bed5b8d`.
- The galaxy-scale coordinate system's 8 open design questions were
  decided this session, and a schema v3->v4 migration,
  `GALACTIC_CENTER_DISTANCE_LY` fix, and a first `galaxyGen.py` were
  written but paused uncommitted in worktree `agent-a36f801e275fb2b71`
  before a full test pass — not part of this release.

## [5.3.3] - 2026-09-07

### Changed
- **Markdown consolidated into `docs/`, renamed for clarity.** Only
  `README.md`/`LICENSE.md`/`CHANGELOG.md` remain at the repo root; every
  other README and loose doc moved into `docs/` with a descriptive name:
  `TODO.md` -> `docs/TODO.md`, `src/html/README.md` -> `docs/html-interface.md`,
  `db/README.md` -> `docs/database-schema.md`,
  `src/api/README.md` -> `docs/api.md`, `apache/README.md` ->
  `docs/apache-deployment.md`, `examples/EXAMPLES.md` ->
  `docs/example-systems.md`, `examples/JSON.md` ->
  `docs/system-file-format.md`. Every cross-reference between them (and
  from code/scripts) was updated to match; a couple of pre-existing
  broken/mismatched links were caught and fixed along the way
  (`docs/system-file-format.md`'s "full command-line reference" link text
  didn't match its own target; `src/html/lib/dbutil.py`'s docstring still
  pointed at the pre-5.3.2 `stellarObjects/` path instead of
  `src/stellarObjects/`).
- **`wsgi.py`, `queryDb.py`, `migrateDb.py` moved into `src/`**, alongside
  `stellarObjects`/`api`/`tests`, so only the `*Gen.py` scripts
  (`sectorGen.py`/`systemGen.py`) are visible at the repo root as CLI
  entry points. Since these three now sit as direct siblings of
  `stellarObjects`/`api` under `src/`, Python's own sys.path[0] (the
  running script's own directory) already makes those packages
  importable -- the sys.path shims 5.3.2 added to them are gone, no
  longer needed (unlike `sectorGen.py`/`systemGen.py`, which stay one
  directory further away at the repo root and keep theirs).
  `migrateDb.py`'s `DEFAULT_DB_DIR` needed an extra `os.path.dirname()`
  level to still resolve to the repo-root `db/`, one directory deeper
  than before.
- **`physicalPlausibility.py` moved to `src/tests/physical_plausibility_cli.py`**,
  matching that directory's naming scheme, without a `test_` prefix so
  pytest doesn't try to collect it as a test module (it's a human-facing
  report generator, not a pass/fail check -- `src/tests/test_physical_plausibility.py`
  remains the actual automated test).
- **`apache/` moved into `examples/apache/`** (its `README.md` moved to
  `docs/apache-deployment.md` per the markdown-consolidation rule above);
  **`examples/*.json` moved into `examples/systems/`**, so `examples/`
  now holds two clearly-separated subfolders (`apache/` deployment
  config, `systems/` recipe files) instead of a flat mix of both kinds of
  example content. `.gitattributes`' LF-pinning rule for
  `apache/*.sh` was updated to `examples/apache/*.sh` to keep matching
  the actual file.
- `setup.py` and `pytest.ini` were **evaluated for a move into `src/` and
  kept at the repo root** -- both are genuinely not possible without
  breaking things, verified empirically rather than assumed:
  `pip install -e .` with `setup.py` moved silently built a bogus,
  empty `UNKNOWN-0.0.0` package instead of erroring (pip's PEP 517
  build only looks for `setup.py`/a full `pyproject.toml` project table
  at the invocation root, and this repo's `pyproject.toml` only declares
  a build backend, no project metadata of its own to fall back on); with
  `pytest.ini` moved, `pytest`'s own config-file discovery (which only
  searches the invocation directory and its parents, never a
  subdirectory) silently fell back to `pyproject.toml` and ignored
  `testpaths`/`pythonpath` entirely -- tests still happened to pass
  either way (coincidentally, via `pytest`'s own unrelated `__init__.py`-walkup
  sys.path behavior and an editable install's global registration of the
  root-level `py_modules`), which is exactly the kind of silent,
  environment-dependent fragility not worth introducing on purpose.

## [5.3.2] - 2026-09-07

### Changed
- **Repo layout: `stellarObjects`/`api`/`tests` moved under a new `src/`
  directory** (`src/stellarObjects/`, `src/api/`, `src/tests/`), so only
  the top-level CLI entry points (`sectorGen.py`, `systemGen.py`,
  `queryDb.py`, `migrateDb.py`, `physicalPlausibility.py`, `wsgi.py`) are
  visible at the repo root. `setup.py` now declares an explicit
  `package_dir` per discovered package rather than a blanket
  `package_dir={'': 'src'}`, since that would have also redirected the
  root-level `py_modules` lookups (`systemGen`/`sectorGen`) into `src/`,
  where they don't live. Every root entry script gained a small
  `sys.path` shim (inserting `src/` before its `stellarObjects`/`api`
  imports) so they keep working without requiring `pip install .` first,
  matching the no-install fallback `src/html/`'s CGI scripts already relied
  on -- those fallbacks (`src/html/lib/dbutil.py`, `src/html/sector.py`,
  `src/html/search.py`) were updated the same way. Two internal
  repo-root-relative path computations (`stellarObjects/webconfig.py`'s
  `_PROJECT_ROOT`, `stellarObjects/_db.py`'s `DEFAULT_DB_PATH`) needed an
  extra `os.path.dirname()` level to still resolve correctly one
  directory deeper; `pytest.ini` gained an explicit `pythonpath = . src`
  so both the entry scripts and the moved packages resolve during tests
  regardless of pytest's own import-mode heuristics.
- **`webconfig.json.example` moved into `src/html/`**; the real, gitignored
  `webconfig.json` stays at the repo root, outside Apache's `src/html/`
  `DocumentRoot`, for the same security reason `db/` already lives there
  (see `docs/WEBCONFIG.md`).
- **`WEBCONFIG.md` moved into `docs/`**, alongside this session's
  `docs/design/`/`docs/analysis/` additions, consolidating loose
  documentation in one place (`README.md`/`LICENSE.md`/`TODO.md`/
  `CHANGELOG.md` stay at the repo root, and per-directory READMEs
  `src/html/README.md`/`apache/README.md`/`db/README.md`/`src/api/README.md`
  stay next to the code they document).

## [5.3.1] - 2026-09-06

### Fixed
- **Evolved-star mass sampling could imply a pre-Big-Bang star.** An
  evolved-class star (`Yerkes != V`, e.g. a giant or supergiant) derives
  its required main-sequence lifespan from its own already-generated mass
  (`Star._calculate_initial_star_age_and_lifespan`'s evolved-star branch);
  for a sub-solar-mass progenitor (roughly under ~0.88 Msun), that implied
  lifespan alone already exceeded `UNIVERSE_AGE_GY` (13.8 Gy) -- meaning
  such a star couldn't actually have finished its main-sequence phase yet
  in the real universe. The [5.3.0] universe-age fix deliberately didn't
  paper over this by capping age below its own required floor, since that
  would produce a self-contradictory star (e.g. a red giant younger than
  its own progenitor's main-sequence lifespan); this is the deeper fix it
  called for. `Star.generate_star`'s evolved-star mass sampling now uses a
  new `_sample_evolved_star_mass_sol` helper (`stellarObjects/starData.py`)
  that rejects and resamples (not clamps, which would just pile an
  artificial spike at the cutoff) any candidate mass whose implied
  main-sequence lifespan would exceed `UNIVERSE_AGE_GY`, capped at
  `program_constants.EVOLVED_STAR_MASS_MAX_RESAMPLE_ATTEMPTS` (100)
  attempts before raising `ValueError` -- in practice a no-op resample for
  every class but III (Giant), whose 0.8-8 Msun range straddles the
  cutoff. Yerkes class VI (subdwarf) is deliberately excluded, since its
  entire allowed mass range (0.1-0.8 Msun) sits below the cutoff and would
  reject every draw; that's tracked as a separate, still-open modeling
  question in `TODO.md`.

## [5.3.0] - 2026-09-06

### Fixed
- **Atmospheric pressure was ~5 orders of magnitude too low for every
  planet class except M.** `planetPhysics.calculate_atmospheric_conditions`
  summed "shell" volumes derived from `planet.radius`/`scale_height`
  (stored in km) directly against `planet.atm_density` (kg/m³) without
  converting km³→m³, undercounting `atmospheric_mass` by ~10⁹×; an
  unexplained `* 7500` fudge factor only clawed back about 4 of those ~9
  orders of magnitude. This rendered as "0.0 kPa" for every atmosphere-
  bearing class but M (Class M never showed it, because a hardcoded clamp
  force-overrode its pressure into an Earth-like range regardless of what
  was computed). Replaced the whole shell-integration loop with the
  closed-form barometric formula for an isothermal, hydrostatic atmosphere
  (`P_surface = ρ_surface · g · H`), which needs no arbitrary zone count or
  fudge factor and lands within the right order of magnitude for both
  Earth-like and Mars-like test cases. `tests/test_planets.py`/
  `tests/test_full_matrix.py`'s pressure assertions were tightened from a
  no-op `>= 0` check into real physical sanity bounds (1 Pa – 10 MPa) so a
  regression like this can't silently pass again.
- **Stars could be reported as hundreds of billions of years old** (e.g.
  "918.77 Billion Years old"), which is impossible since the universe
  itself is only ~13.8 billion years old — even though the star's age never
  actually exceeded its own (very long, and realistically so: real M dwarfs
  are predicted to live trillions of years) lifespan. Added a
  `UNIVERSE_AGE_GY = 13.8` ceiling, applied in both the main-sequence and
  evolved-star branches of `Star._calculate_initial_star_age_and_lifespan`
  and in `adjust_age_for_planets`'s age-adjustment logic. (White dwarfs
  already had their own bounded 0.1–12.0 Gy cooling-age range and needed no
  change.) A sub-solar-mass evolved-star progenitor whose own
  main-sequence lifespan already exceeds the universe's age is a separate,
  deeper mass-sampling issue, noted in `TODO.md` rather than papered over
  here.
- **Rare `IndexError: string index out of range` crash in name
  generation.** `generate_phoneme_salad_name` (used for star, planet, and
  moon names) could crash when a base name containing a literal apostrophe
  (e.g. `PLANET_NAMES`'s "Hi'iaka") landed adjacent to a `UNIVERSAL_PHONEMES`
  chunk also ending in an apostrophe (e.g. "ch'"), producing a literal `''`
  in the assembled name; the final capitalization step split on `'` and
  indexed `part[0]` on every piece, crashing on the empty piece between the
  two apostrophes. Not specific to binary systems — any name generation
  call could hit it given enough attempts, which is why it only surfaced
  rarely across a large batch of generated systems.

### Changed
- Class M's (and Class P's) hardcoded gravity/pressure/temperature
  "forcing" clamps in `planetPhysics.py` are commented out, not deleted —
  the fixed atmospheric-pressure formula already lands close to realistic
  ranges on its own, so the band-aid that had been masking the pressure bug
  for Class M is no longer needed. May be restored later if Class M/P need
  tighter narrative guarantees again.
- `star_systems` gained a `location` column (schema bumped to v3): for any
  system placed in a sector, its sector's name plus distance (in
  light-years) to its 3 nearest neighboring systems, e.g. `"Voranthis
  Kelmoor -- nearest: Aldenar (4.2 ly), Brekthos (7.8 ly), Corvane (9.1
  ly)"`. Computed once at write time (mirroring how `quadrant` is already
  derived and persisted) via the existing `SpaceSector.nearest_neighbors`/
  `distance_between` helpers; `NULL` for systems never placed in a sector.
  Existing v1/v2 databases migrate automatically (backed up first, as
  usual) and have their `location` backfilled from already-stored position
  data.
- The web interface's Search link is now reachable from every page (it
  previously required detouring back through the database-browse page) —
  `src/html/lib/page.py`'s shared page header now includes it whenever a
  database is selected.
- The database-picker landing page (`src/html/index.py`) now redirects straight
  to browsing the one database present, if the configured database
  directory contains exactly one `.db` file, instead of always showing the
  picker.

### Added
- `webconfig.json` (repo root, gitignored — `webconfig.json.example` is the
  committed template): site-level configuration, currently `site_name` and
  `base_url`, plus unused placeholder fields (`db_username`, `db_password`,
  `db_name`) reserved for a possible future non-SQLite backend. Kept
  outside `src/html/`'s served document root, the same way `db/` already is.
  See [`WEBCONFIG.md`](WEBCONFIG.md) for full documentation.

## [5.2.5] - 2026-09-06

### Changed
- **Database schema bumped to v2**: moons split out of the shared
  `planets` table into their own `moons` table (`planet_id` FK to the
  planet they orbit, plus their own `moon_evolutionary_paragraphs`/
  `moon_reflection_spectrum` child tables), instead of self-referencing
  via `planets.parent_planet_id`/`is_moon`. This is what `src/html/search.py`'s
  attribute tags needed to actually tell a planet from a moon: previously
  a "Class D Planet" tag queried `planets` with no way to exclude
  `is_moon=1` rows of the same class, so it silently listed moons too.
  `stellarObjects/_db.py` gained a dedicated `insert_moon` (mirroring
  `insert_planet`, but writing to the new table); `insert_planet` no
  longer recurses into itself for moons.
- `src/html/search.py`'s tag facets and name search now follow the same
  split: "Planet Class"/"Planet Body Type"/"Planet Supported Life
  Chemistry" only ever match top-level planets, with an identically-
  shaped "Moon Class"/"Moon Body Type"/"Moon Supported Life Chemistry"
  set of tags (and a separate autocompleting "Moon name" field) for
  moons -- each still only rendering the tag buttons for values actually
  present in the chosen database. `src/html/system.py`'s planet/moon table
  now reads moons from the new table instead of a recursive
  self-join, and no longer needs to recurse (moons never generate their
  own moons -- confirmed by the existing
  `tests/test_moons.py::test_moons_cannot_themselves_have_moons`).

### Added
- `migrateDb.py` (repo root): converts every `*.db` file in a directory
  from schema v1 to v2, backing up each original first
  (`<name>.db.v1-backup-<timestamp>.db`) -- a no-op for a database
  that's already current. Backed by a new `stellarObjects._db.
  migrate_database`/`_migrate_v1_to_v2`, which builds the converted
  database at a temporary path and only atomically swaps it into place
  at the very end, so a crash or error partway through leaves the
  original file (and the backup already made) untouched. `install.sh`
  now runs this automatically (a new step, right after installing the
  Python package) over every database in `db/`, so `sudo ./update.sh`
  keeps an existing deployment's database working across the schema
  change with no manual step.
- `tests/test_db_migration.py` and `tests/test_db_persistence.py`: the
  first automated tests of `stellarObjects/_db.py`'s save/migrate paths
  (previously "manually smoke-tested via `sectorGen.py`" per `TODO.md`).
  Between the two, they cover the v1->v2 migration (including the
  moon-owned child-row split, which the real generated data used for
  manual smoke-testing happened not to exercise) and a genuine
  generate-then-save-then-query round trip through `insert_star_system` --
  the latter caught a real bug during development (`insert_moon`'s
  `INSERT` had one more column than value placeholder).

## [5.2.4] - 2026-09-06

### Added
- `src/html/search.py`: a faceted/name search page for the web interface.
  Two complementary ways to find an object in the chosen database:
  - Click-to-filter tag buttons for object type (star/planet/moon/
    asteroid belt), star spectral class and Yerkes luminosity class,
    planet class and body type, and supported life chemistry -- each
    group renders a button only for values actually present in that
    specific database (e.g. no Yerkes-luminosity buttons beyond "Main
    Sequence" if nothing else was generated), per `SELECT DISTINCT ...
    GROUP BY` queries against `stars`/`planets`/`asteroid_belts`, not a
    fixed enumeration. Multiple tags within one group OR together (any
    matching value); tags across groups AND together, except that
    selecting an explicit Object Type tag acts as a master filter (e.g.
    selecting only "Stars" hides the Planets panel even if a
    planet-class tag also happens to be selected). Clicking a tag
    toggles it via a plain link that rewrites the query string, so this
    works with JavaScript entirely disabled, same as every other page in
    `src/html/`.
  - A name search with one field each for sector, star system, star, and
    planet/moon names, each with its own HTML5 `<datalist>` for
    autocomplete (native browser suggestions, no JavaScript, populated
    from that entity's own distinct names). Asteroid belts have no name
    of their own (per `schema.sql`'s own note on this), so they're only
    reachable via the "Asteroid Belt" object-type tag.

  Every filter is a parameterized SQL query (`IN (...)`/`LIKE ... ESCAPE
  '\'`), and each result panel is capped (300 rows, with a "showing the
  first N" note) to keep a broad tag click from dumping an entire large
  database into one page. Linked from `browse.py`'s breadcrumb.

## [5.2.3] - 2026-09-06

### Changed
- `install.sh` no longer installs the Python package via the deprecated
  `python3 setup.py install`; it now runs a proper, build-isolated `pip
  install --upgrade --force-reinstall "$SCRIPT_DIR"` instead.
  `setup.py install` broke on a deployed Ubuntu 20.04/Python 3.8 host
  *twice* in a row, both times because setuptools' own vendoring shim
  (`extern`) prefers a real, already-installed copy of a dependency it
  vendors over its own newer bundled copy whenever one is importable —
  first `importlib_metadata` (`AttributeError: ... no attribute
  'EntryPoints'`), then, after that was patched around, `packaging`
  (`TypeError: canonicalize_version() got an unexpected keyword argument
  'strip_trailing_zero'`) — because that distribution's old apt-provided
  copies of each in turn shadowed the working vendored one. Chasing each
  shadowed dependency individually only fixes the one that broke that
  day, not the next one; `pip install`'s build isolation builds the
  package in a throwaway environment that can't see the system's
  site-packages at all, closing the whole bug class instead of patching
  it dependency-by-dependency, and means `install.sh` no longer needs to
  upgrade the system's global `setuptools`/`importlib_metadata` at all
  for this step. `--force-reinstall` is deliberate: a plain `pip install
  .` skips reinstalling when pip thinks the same version is already
  installed — true on every `update.sh` run between version bumps in
  `stellarObjects/_version.py` — which would otherwise silently leave the
  previous run's install in place instead of the source `update.sh` just
  pulled.
- `install.sh` and `update.sh` now re-`chmod +x` every `*.sh` file in the
  repo (not just `src/html/*.py` and one hardcoded `apache/set-permissions.sh`
  path), including themselves. A `core.fileMode=false` git config on the
  authoring machine drops the executable bit on any file type on
  checkout, not just `src/html/`'s `.py` scripts, and `update.sh` invokes
  `install.sh` directly (`"$SCRIPT_DIR/install.sh"`, not `bash
  install.sh`) — if a pull had dropped *its* executable bit, the shell
  would refuse to exec it before `install.sh`'s own permission fix ever
  got a chance to run. `update.sh` now fixes this itself, immediately
  before invoking `install.sh`, so a dropped bit on `install.sh` (or on
  `update.sh` for its own next run) no longer requires a manual `chmod
  +x` to recover from. (One unavoidable exception: the very first
  `update.sh` run to pull this fix is still executing the old script
  text it already had open when `git pull` replaced the file on disk out
  from under it, so that one migration may still need a manual `chmod +x
  install.sh` — every run after that starts from the fixed script.)

## [5.2.2] - 2026-09-06

### Added
- `src/html/system.py` now renders a system's description as actual HTML by
  default (`?view=rendered`), via a new small, purpose-built Markdown-to-
  HTML converter (`src/html/lib/mdconvert.py`) targeting exactly the narrow
  Markdown subset `StarSystem.__str__` generates (headers, pipe tables,
  paragraphs, `<sup>` exponents) — not a general-purpose parser, and no
  new dependency. The original raw wikitext/Markdown source is still one
  click away (`?view=source&format=...`) for copy-pasting into a wiki.
  Escapes every block in full before emitting markup, then narrowly
  re-enables only the one legitimate raw-HTML pattern generated content
  contains, so a mischievous `--name`/`--star-type` value can't inject
  live HTML into a rendered page (covered by new tests in
  `tests/test_mdconvert.py`).
- `src/html/static/style.css` rewritten as a small design system: CSS custom
  properties, automatic light/dark via `prefers-color-scheme`, card-style
  panels, badges, breadcrumbs, a proper type scale, hover states,
  focus-visible outlines, and a responsive breakpoint. Applied
  consistently across every page (`index.py`, `browse.py`, `sector.py`,
  `system.py`, and the shared error page in `lib/page.py`), not just the
  system page.
- `update.sh` (repo root, alongside `install.sh`): pulls the latest
  changes from git and re-runs `install.sh` so permissions stay correct
  afterward. Refuses to run over uncommitted local changes, and pulls
  with `--ff-only` (fails loudly rather than creating a surprise merge
  commit if history has diverged) instead of a plain `git pull`.

### Fixed
- `apache/set-permissions.sh` and `install.sh` both failed to make every
  `*.py` file under `src/html/` executable by `www-data` — `install.sh`'s own
  `chmod +x` step used `find -maxdepth 1`, silently skipping
  `src/html/lib/*.py`, and `set-permissions.sh` only `chgrp`'d (group
  ownership) rather than `chown`'d (user *and* group) the deployed
  directories. Both fixed: the `-maxdepth 1` restriction is gone, and
  `set-permissions.sh` now `chown -R`s to the detected Apache user:group
  and reports how many `.py` files it made executable, so a wrong path is
  obvious rather than silently matching nothing.

## [5.2.1] - 2026-09-06

### Added
- `src/html/`: a small, dependency-free web interface (plain Python CGI
  scripts, standard library only) for browsing the SQLite databases from
  `db/README.md` — pick a database, drill into its sectors and star
  systems, and view (or copy, via a `<textarea>`) the rendered
  wikitext/Markdown page saved for each one. See `src/html/README.md`.
- `apache/`: deployment tooling for the web interface —
  `planetgen.conf.example` (an example Apache2 virtual host) and
  `set-permissions.sh` (detects the user/group Apache2 actually runs as
  and sets ownership/permissions on the deployed `src/html/`/`db/`
  directories accordingly).
- `install.sh`: a one-shot Linux installer tying the above together —
  runs `setup.py install`, pre-fetches the NLTK `words` corpus into a
  shared, world-readable location, makes the CGI scripts executable,
  enables Apache's `cgid` module, and runs `apache/set-permissions.sh`;
  prints the one remaining manual step (copying/enabling the example
  vhost) rather than touching Apache's site configuration itself.
- `stellarObjects/names.py` gained `UNIVERSAL_PHONEMES`: a pool of ~100
  short, ASCII-7-bit-printable phoneme chunks romanized from roughly 18
  language families (Romance, Germanic, Slavic, Arabic/Hebrew/Persian,
  Turkish, South Asian, Mandarin, Japanese, Korean, Vietnamese,
  Austronesian, Polynesian, Bantu, Mesoamerican, Andean, Celtic,
  Finno-Ugric, Caucasian). `generate_phoneme_salad_name` now splices one
  of these into every generated name — star, planet, moon, and sector
  alike — with `UNIVERSAL_PHONEME_CHANCE` (40%) odds, widening the
  cultural range generated names are drawn from beyond each type's own
  base name list.

### Changed
- Removed `SECTOR_DESIGNATORS` from `stellarObjects/names.py` (dead code
  — defined but never referenced anywhere).
- `generate_phoneme_salad_name` gained an `allow_split` parameter
  (default `True`, unchanged for stars/planets/moons); see the sector
  name bug fix below.

### Fixed
- Sector names could come out as 3-4 words instead of the intended 2 (in
  a 500-name stress test, this happened 95% of the time).
  `generate_phoneme_salad_name` can split a long result into two
  space-separated words on its own (`split_long_word`), and
  `sectorGen.generate_sector_name` already joins two independent calls
  into one name — an internal split on either half silently produced
  3-4 words in the final result. `generate_sector_name` now passes the
  new `allow_split=False` for both halves.
- `stellarObjects/names.py` unconditionally called
  `nltk.download('words', quiet=True)` at import time; `nltk`'s
  `download()` always targets the *current user's* default download
  directory and attempts to create it, regardless of whether the corpus
  already exists elsewhere on `nltk.data.path` — this broke the web
  interface entirely under Apache's locked-down `www-data` account
  (`PermissionError: ... '/var/www/nltk_data'`), since that account has
  no writable home directory. It had only ever "worked" for the CLI
  tools because those ran as `root`. Now checks
  `nltk.data.find('corpora/words')` first and only downloads on
  `LookupError`; `install.sh` pre-fetches the corpus into a shared,
  world-readable path so that check always succeeds once installed.
- `apache/set-permissions.sh` hard-failed when it couldn't detect
  Apache's user/group (e.g. because apache2 isn't started/enabled yet,
  which is exactly the case during a fresh `install.sh` run) — now warns
  and defaults to the standard Debian/Ubuntu `www-data:www-data` instead.
- CGI scripts could be deployed non-executable regardless of a local
  `chmod +x`, because a `core.fileMode=false` git config on the
  authoring machine silently drops the executable bit before it reaches
  a commit. `install.sh` now `chmod +x`s the CGI scripts and
  `apache/set-permissions.sh` directly on every install, independent of
  whatever mode git happened to store.

## [5.2.0] - 2026-09-06

### Added
- Full SQLite database persistence: `stellarObjects/schema.sql` defines the
  schema, and the new `stellarObjects/_db.py` writes an already-generated
  `SpaceSector` (every system, star, planet, moon, and asteroid belt,
  plus a rendered copy of the wiki page in both wikitext and Markdown)
  into it in a single transaction. Documented column-by-column in the new
  `db/README.md`. `sectorGen.py` now calls this automatically on every
  run, saving to `db/planetgen.db` by default (overridable via the new
  `--db-path` option); the database file itself is gitignored.
- Sector name generation: sectors now get a random two-word name (e.g.
  "Voranthis Kelmoor") drawn from a new sector-flavored name list
  (`SECTOR_NAMES`/`SECTOR_PREFIXES`/`SECTOR_SUFFIXES` in
  `stellarObjects/names.py`, pulling from real galactic structures and
  well-known science-fiction sector names) instead of reusing star names
  with "Sector" appended.
- `Star`, `Planet`, `BinaryStarProxy`, and `AsteroidBelt` each gained a
  `get_table_properties()`/`get_composition_summary()` method that returns
  the same already-formatted values their `to_paragraph_list()` renders
  into text, so the new database layer can store exactly what was
  published without duplicating any formatting logic.

### Changed
- **Breaking:** `sectorGen.py`'s `-n` short flag now means `--name`
  (hard-sets the sector's own name) instead of `--num-systems`, which no
  longer has a short flag; the old `--sector-name` option was renamed to
  `--name`.
- Distances stored in the database use a two-tier unit convention
  (milliparsecs for sector-scale position/geometry, kilometers everywhere
  else); `stellarObjects/utils.py` gained `ly_to_milliparsecs`/
  `milliparsecs_to_ly` and `physical_constants.py` gained
  `AU_PER_PARSEC`/`AU_PER_MILLIPARSEC` to support the conversion, used only
  by the persistence layer.

### Fixed
- Atmospheric scale height (`planetPhysics.calculate_atmospheric_conditions`)
  was computed in meters but combined directly with `planet.radius` (in
  km) without converting units first, throwing off atmosphere thickness
  and volume for every planet with an atmosphere.
- `utils.calculate_object_mass` and `Planet.__init__` each independently
  (and inconsistently) converted radius-in-km to a volume, one of them
  mixing up its km/m factor; `Planet` now takes its volume solely from
  `calculate_object_mass`'s corrected calculation instead of recomputing
  it a second time.

## [5.1.0] - 2026-09-06

### Added
- `--version` option on both `systemGen.py` and `sectorGen.py`, printing the
  program's version, this repository's URL
  (https://github.com/dwhagar/planetGen), and a license summary, then
  exiting immediately.
- `stellarObjects/_version.py`: a single, dependency-free source of truth
  for the project's version number, shared by both CLI scripts and
  `setup.py`.

### Changed
- `setup.py` now reads its `version` from `stellarObjects/_version.py`
  instead of a hardcoded, never-updated placeholder.
- README updated with current version information and a link to this
  changelog.

## [5.0.1] - 2026-09-06

### Fixed
- Flavor text (both system-level and planet/moon-level) was being rolled —
  and shared `system_config` counters mutated — every time a system or planet
  was *rendered* (`__str__`/`to_paragraph_list()`) instead of once at
  generation time, so rendering the same object twice would silently re-roll
  and double-count it. It is now decided exactly once per object during
  generation (`StarSystem.__init__`), and rendering is a pure, idempotent
  read of the already-decided text.

### Changed
- Internal property cleanup across `Planet`, `Star`, and `AsteroidBelt`:
  removed five redundant star-property snapshots from `Planet` that
  duplicated data already reachable via its `star` reference (one of
  which — `star_radius` — had a latent double-unit-conversion bug whenever a
  moon was generated); renamed ambiguous or colliding attributes
  (`evolution` → `evolutionary_speed`, `hab` → `habitable_zone`, and
  `type` → `body_type` on `Planet`/`AsteroidBelt` so it no longer shares a
  name with the unrelated `Star.type`); filled in missing/stale docstring
  documentation.

### Added
- Regression tests proving rendering is idempotent (calling
  `to_paragraph_list()`/`__str__()` twice produces identical output and
  performs no further mutation).

## [5.0.0] - 2026-09-05 (evening)

### Added
- `spaceSector.py`: a new `SpaceSector`/`SectorSystemEntry` layer that places
  generated `StarSystem`s at `(x, y, z)` positions within a cubic sector, with
  save/load to JSON (`SystemConfig` gains `to_dict`/`from_dict` to support
  this).
- Poisson-disk-based sector growth for realistic inter-system spacing, plus
  named-location ("quadrant") formatting for a system's position within a
  sector.
- `TODO.md`: a long-term roadmap toward full database storage and a web
  interface.

### Changed
- Reworked sector growth to sample candidate positions and then fine-tune
  each one against every existing neighbor, rather than a single-pass
  placement.
- Sector density and the minimum allowed separation between two systems are
  now grounded in real astronomical data — density from real local stellar
  surveys, and minimum separation from each system's own gravitationally
  derived Hill-sphere radius — rather than arbitrary constants.

## [4.1.0] - 2026-09-05 (afternoon)

### Added
- JSON system-file support: a system's exact contents (star type, forced
  features, per-orbit slots) can now be specified via a JSON recipe file,
  with several example systems included under `examples/`.
- The project's first automated test suite (pytest), covering example
  systems, moon generation, and planet generation.
- `sectorGen.py` CLI entry point, ahead of the sector-generation logic it
  would call the next day.

### Changed
- Major standardization pass across CLI options and internal call
  signatures.
- Renamed the main script from `planetGen.py` to `systemGen.py`, reflecting
  that it now generates full star systems rather than a single planet.

## [4.0.0] - 2026-09-03

### Changed
- Major internal refactor: generation logic was split out of the `Planet`
  class into dedicated `planetPhysics.py` (physical/orbital generation) and
  `planetLife.py` (life chemistry and evolution) modules, so `Planet` itself
  holds state and presentation while generation logic lives alongside it as
  free functions.
- The monolithic `constants.py` was split into `physical_constants.py`
  (real-world physical constants) and `program_constants.py`
  (generation/tuning constants), for readability and maintainability.
- Import structure cleaned up throughout the package; install script
  (`setup.py`) updated.

## [3.0.0] - 2026-07-07

### Added
- Binary star systems: a new `doubleStar.py` module (`BinaryStarProxy`)
  represents a double star system as a single effective star (combined
  mass/luminosity/habitable zone) for the purposes of planet placement.

### Fixed
- Primary/secondary star identification and overall system age calculation
  for binary systems.
- An age/lifespan bug affecting planets.

### Changed
- Binary star output formatting refined; the large-star-forcing option
  updated to work correctly alongside binary generation; README updated.

## [2.2.3] - 2026-06-27

### Added
- More flavor text variants across different planet classes.

### Changed
- Refined repeat-prevention and selection logic for flavor text so the same
  text is less likely to recur in quick succession.

## [2.2.2] - 2026-06-26

### Fixed
- Flavor text selection options.

### Changed
- Switched random-number generation to Python's `secrets` module for
  higher-quality entropy, avoiding similar successive sequences.

## [2.2.1] - 2026-06-25

### Changed
- Wording and grammar adjustments throughout the generated text output.

## [2.2.0] - 2026-06-24

### Added
- Options to specify a system's age and to force (or forbid) intelligent
  life.
- Flavor text: a random chance of extra descriptive "sensor" text being
  appended to a system or a planet's description.

### Changed
- Centralized magic numbers into `constants.py` as named constants, making
  the physics/generation formulas easier to read.
- Unified the set of properties referenced across all object types onto a
  shared convention.
- README updated for the new CLI options.

## [2.1.1] - 2026-06-23

### Changed
- Asteroid belt data and logic split out into its own module
  (`asteroidData.py`).
- Expanded asteroid belt composition descriptions and grammar.

### Added
- Option to specify a custom name for the star system.

## [2.1.0] - 2026-06-21

### Added
- Stellar age is now generated and reported for every system.
- Evolutionary timeline narratives: systems now estimate the plausibility of
  life at every stage, from simple single cells through technological
  civilizations.

### Fixed
- Assorted output-formatting issues.

## [2.0.0] - 2026-06-20

### Added
- Life chemistry system: planets and star types now carry information about
  the chemical processes that could plausibly give rise to life, as a
  foundation for the evolutionary modeling that follows.

### Changed
- Command-line options reworked so multiple flags combine correctly together.

## [1.5.0] - 2026-06-12

### Added
- Option to generate a star system with no planets.
- Option to specify a star's exact spectral type from the command line.

## [1.4.0] - 2026-06-06

### Added
- Markdown export support, and the ability to write generated output to a
  file (in addition to the console).
- Heliosphere radius and stellar gravitational-influence ("system perimeter")
  calculations, describing the outer boundaries of a system.

### Fixed
- Name-generation issues.
- Scientific-notation formatting issues.

## [1.3.1] - 2026-06-05

### Changed
- Further cleanup of generated text output for clarity.
- Refinements to procedural naming.

## [1.3.0] - 2026-06-04

### Added
- Procedural name generation for stars, planets, and moons.
- Project README.

### Changed
- General code cleanup pass.

*(~23-month gap in development between 2024-07-15 and 2026-06-04.)*

## [1.2.2] - 2024-07-15

### Fixed
- Star age/temperature calculations corrected to consistently use Kelvin
  throughout.

## [1.2.1] - 2024-07-14

### Fixed
- Class M worlds are now reliably clamped to Earth-like gravity and
  atmospheric pressure.
- Class P (and other) planet classes now get appropriate temperature
  treatment for their class.

### Changed
- Asteroid belt sizing now makes better use of the habitable zone.

## [1.2.0] - 2024-07-13

### Added
- Options to force a habitable world, a large star, and/or an asteroid belt
  (forcing both a habitable world and an asteroid belt together
  automatically forces a large star).
- Options to control the overall size (object count) of a generated system.
- Expanded command-line help text.

### Fixed
- A long-standing infinite-loop bug in star generation.
- Several edge cases in generation logic.

### Changed
- Capped the maximum number of system objects at 500 to prevent runaway
  generations.
- Refined the system description text and radius output formatting.

## [1.1.1] - 2024-07-11

### Fixed
- Completed and debugged the moon system: moon orbital distances are now
  tracked and calculated correctly, and remaining orbital overlap between
  planets and asteroid belts was removed.

### Added
- Weighted probabilities for moon class selection.

### Changed
- Further wording, grammar, and formatting refinements.

## [1.1.0] - 2024-06-29

### Added
- First implementation of the moon-generation function (not yet wired into
  planet creation or tested at commit time).

## [1.0.2] - 2024-06-28

### Added
- Mass validation for a specified planet, laying the groundwork for moon
  generation.

## [1.0.1] - 2024-06-26

### Fixed
- Orbital placement now accounts for planetary position and Hill radius, so
  systems no longer generate overlapping or too-closely-spaced orbits.
- Asteroid belt inner/outer bounds no longer overlap neighboring planets.

### Changed
- Further readability tweaks to the generated text output.

## [1.0.0] - 2024-06-25

### Added
- First fully working end-to-end planet generator, validated against real
  Earth reference values ("I think it finally works!").
- Special-cased handling for Class N worlds (boosted atmospheric density and
  molar mass).

### Changed
- Multiple formatting and descriptive-text passes on the generated output.

## [0.2.0] - 2024-06-23

### Added
- Core stellar physics: working star generation and planet-count estimation.
- Initial atmospheric modeling: mass, gravity, and atmosphere calculations.

### Fixed
- Numerous early bugs in surface pressure and temperature calculations.

## [0.1.1] - 2024-06-22

### Fixed
- Continued debugging to get the initial prototype running end-to-end (moved
  development from VS Code to PyCharm along the way).

## [0.1.0] - 2024-06-21

### Added
- Initial project scaffolding: the `stellarObjects` package with the first
  `Planet`, `Star`, and `StarSystem` classes, and a `main.py` entry point.
- First rough pass at generating a single planet's basic properties.