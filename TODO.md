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
that did the work, not here) — read `stellarObjects/schema.sql` and
`db/README.md` for the actual current schema rather than this document.
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
  `db/README.md`'s schema history for what changed at each version).
  Distances are stored in milliparsecs for sector-scale placement and
  kilometers for everything else (see `schema.sql`'s header comment for
  the full convention). Round-trip tests verify byte-for-byte render
  fidelity through the actual database, not just in-memory dicts.
- **Phase 3 (CLI tools use the database)**: `sectorGen.py`/`systemGen.py`
  save every run to the database unconditionally (`--db-path` to
  override the default `db/planetgen.db`); `queryDb.py` is a read-only
  CLI for listing/searching what's stored (`sectors`, `systems`, `near`
  subcommands).

## Unmerged work from the 2026-09-06 session (review before continuing)

Five background agents ran concurrently that night alongside the API/wiki
work already merged into `main` (commit `ee8daab`). None of the below are
merged — each is a separate git worktree/branch under `.claude/worktrees/`
still sitting on top of the pre-session commit. Review, cherry-pick or
merge, then delete the worktree (`git worktree remove`) once done.

- [ ] **Evolved-star mass fix** (branch `worktree-agent-a4307a26b82478c8d`)
  — complete and tested (8590 tests passing, 3 new regression tests). Fixes
  the pre-Big-Bang evolved-star mass issue via reject-and-resample in
  `stellarObjects/starData.py`. Also bumps version to 5.3.1 and adds a
  CHANGELOG entry — decide whether to keep that version bump given other
  unreleased changes already on `main`. Flags one new narrow follow-up:
  Yerkes class VI (subdwarf) is exempted from the check because its entire
  mass range (0.1-0.8 Msun) is below the universe-age cutoff, which needs
  its own look eventually.
- [ ] **Galaxy coordinate system design** (branch
  `worktree-agent-afe9fc7f7137bafe8`) — a full written proposal at
  `docs/design/galaxy-coordinate-system.md` (not code): spherical-to-XYZ
  sector placement, **parsecs** recommended for galaxy-scale distance, a
  shell-based radial tiling scheme (Fibonacci-sphere sequence per shell),
  and a v3->v4 schema migration sketch for the `sectors` table. Surfaces one
  finding needing a decision: `physical_constants.GALACTIC_CENTER_DISTANCE_LY`
  is currently a fixed constant assumed by every star's Hill-sphere math,
  which becomes wrong once sectors have real galactic positions. Has 8 open
  questions flagged for review before `galaxyGen.py` implementation starts
  (this doc is the prerequisite for Phase 4 below).
- [ ] **Habitable/atmosphere surface-condition sanity review** (branch
  `worktree-agent-a75772018725a8694`) — analysis at
  `docs/analysis/habitability-atmosphere-sanity-review.md`, not a fix.
  Generated 300+ samples per class and found real problems (not just "the
  clamps are fine to leave removed" as hoped):
    - Class P is no longer meaningfully colder than Class M or any other
      ecosphere terrestrial class — the disabled clamp was the *only* thing
      giving P its cold identity; it has zero built-in bias now.
    - Class M's atmospheric pressure never approaches 1 atm across 600
      samples (mean ~1/3 atm) — contradicts the assumption that the fixed
      pressure formula "lands close to realistic ranges on its own."
    - Likely root cause for both: the `greenhouse_factor` formula in
      `planetPhysics.py` (~L402) looks physically inverted — rewards an
      atmosphere for being *far* from CO2's molar density rather than for
      having more CO2 in it.
    - Atmospheric pressure is mathematically independent of gravity for
      every class (proven algebraically and numerically) — the
      `scale_height`/gravity terms cancel exactly, which is why gas giants
      spanning a >1000x gravity range all produce nearly identical pressure.
    - Gas-giant density blending (`planetPhysics.py` L308-311) can produce
      implausibly "fluffy" planets (densities as low as 0.026 g/cm^3) —
      likely a real bug, fix depends on clarifying what the blended ratio
      should represent.
    - Most ecosphere terrestrial classes are statistically indistinguishable
      in pressure/temperature since they share the same `"t"`-keyed
      parameter ranges regardless of class-specific flavor text.
  One trivial doc-comment fix already applied and tested (mislabeled
  `ATMOSPHERE_DENSITY["t"]` range comment). Everything else above is a
  recommendation, not yet implemented — **this should probably be resolved
  before further physics tuning**, since it affects M/P/gas-giant
  generation broadly, not just the originally-scoped clamps question.
- [ ] **Physical-plausibility anomaly finder** (branch
  `worktree-agent-afdb05a2b331ff123`) — was cut off mid-run by a
  usage-limit error before it could report results; **resumed 2026-09-07**
  in a fresh worktree that reads the original's `physicalPlausibility.py`
  (CLI, ~100 lines) and `stellarObjects/plausibility.py` (design/logic,
  ~407 lines — hard invariant checks plus Tukey's-fences statistical
  outlier detection per class) as a starting point, finishes verifying it
  actually runs correctly, adds test coverage, and reports real findings.
  Given the sanity-review findings above, this tool should be immediately
  useful for confirming how far M/gas-giant pressure and Class P
  temperature actually drift.
- [x] ~~Wiki-URL reachability check~~ and ~~Backend API~~ (Flask) — both
  already merged into `main` in commit `ee8daab`.

## Investigate Further

- [ ] Ways to render an image or maybe provide a web interface to visualize the location of 2 points within the galactic space.
- [ ] Introducing realistic orbital paths and speeds to all bodies in space, would need a dedicated update script to update like once a month or something to adjust all of the coordinates.
- [ ] Search parameter for searching by not only planet class but planet size or in the tagged search field sort by planet size.

## Organize File System

- [ ] Organize all documentation into a dedicated docs/ directory under the
  repo root.
- [ ] Do the same for most of the python code that is not part of the HTML interface
  into an src directory so that only the Gen python entry points are visible from
  repo root.
- [ ] All test code should go in src as well.
- [ ] Consider moving api directory into html to expose API end point?
- [ ] Move webconfig file into /html and change the name of examples.

## File Management

- [ ] Database backups during schema's should be gzipped and not available on the database picker.
- [ ] Don't upgrade backups, if a backup database is present from an upgrade don't ugprade the backup.

## Population and Politics

- [ ] Need to start thinking about assigning government ownership to a particular star system in the database so that together the systems make territories that are mapped out in 3D space by the star systems.
- [ ] Worlds with life on them can be flagged for generated names of dominant races.
- [ ] Probably going to need a space fairing species database.
- [ ] Need to think about under-developed / older civilizations and the differences and how to store and present that data based on society age.

## Phase 4 — Galaxy-scale generation

- [ ] Galaxy-scale coordinate system: no longer a blank design slate — a
  full proposal exists (see "Unmerged work" above) at
  `docs/design/galaxy-coordinate-system.md`: spherical-to-XYZ sector
  placement, **parsecs** for galaxy-scale distance, a shell-based radial
  tiling scheme (Fibonacci-sphere sequence per shell), and a v3->v4
  `sectors` schema migration sketch. Not yet reviewed/approved — the doc
  itself flags 8 open questions (eager vs. lazy shell generation, the
  disk-density envelope, the `GALACTIC_CENTER_DISTANCE_LY` fixed-constant
  problem it surfaced, and others) that need a decision before this is
  considered settled.
- [ ] `galaxyGen.py` — a new top-level CLI script (sibling to
  `sectorGen.py`/`systemGen.py`) that generates and persists many sectors
  as one galaxy, reusing `sectorGen.py`'s own per-sector generation/save
  logic for each one. Gated on the coordinate-system design above being
  reviewed and merged first.

## Phase 5 — Web interface (long-term; needs its own dedicated planning pass)

An interim, dependency-free read-only browser already exists at
[`html/`](html/README.md) (plain Python CGI scripts, no framework) plus
[`apache/`](apache/README.md) (example vhost config + `set-permissions.sh`),
meant for a single-user/small-scale Apache2 deployment today rather than
the full multi-user vision below.

- [x] Backend API framework: **Flask**, chosen over FastAPI/Django REST
  Framework — no ORM opinion (fits the existing raw-`sqlite3` persistence
  layer with zero glue), deploys via `mod_wsgi` in the same Apache process
  model the interim `html/` CGI browser already uses, and mounts at `/api/`
  alongside `html/` for an incremental rollout. FastAPI's async/auto-docs
  advantages don't pay for themselves yet (no separate frontend consuming
  the API, and `sqlite3`'s driver is synchronous regardless of framework);
  revisit if a dedicated frontend makes API-contract docs valuable. A
  read-only scaffold now exists — see [`api/README.md`](api/README.md) for
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
    - **`api/config.py` / `queryDb.py` porting**: both currently assume a
      SQLite file path (`DB_PATH`/`--db-path`) — becomes a connection
      string/host+credentials, need a secrets-handling story (env vars at
      minimum) rather than a bare file path.
    - **Data migration**: existing `.db` files need a one-time export/
      import into MySQL; no tooling for this exists yet.
    - Revisit tooling (e.g. an ORM/migration framework) at that point if
      the hand-ported approach proves painful, per the original note here.

### Near-term: interim `html/` browser enhancements

- [x] **Distance calculator** / **"Nearest N systems within radius R"**:
  superseded by `star_systems.location` (schema v3) — every system placed
  in a sector now stores its sector name plus distance (in ly) to its 3
  nearest neighbors, computed at write time from the existing
  `position_x/y/z_mpc` columns. See `CHANGELOG.md` [5.3.0].
- [x] **Search reachable from every page**: the shared page header
  (`html/lib/page.py`) now includes a Search link whenever a database is
  selected, instead of requiring a detour back through `browse.py`.
- [x] **Skip the database picker when only one database exists**:
  `html/index.py` now redirects straight to `browse.py` if the configured
  database directory contains exactly one `.db` file.
- [x] **Site configuration file**: `webconfig.json` (repo root, gitignored;
  `webconfig.json.example` committed) holds `site_name`/`base_url` today,
  plus unused placeholder fields (`db_username`/`db_password`/`db_name`)
  for a possible future non-SQLite backend. See [`WEBCONFIG.md`](WEBCONFIG.md).
- [x] **Wiki-URL reachability check + clipboard fallback**: implemented in
  `html/system.py` (merged in commit `ee8daab`), deliberately the bare
  minimum rather than the fuller version this item originally speculated
  about below — a plain HEAD-then-GET check treating any non-200 status,
  timeout, or connection error as unreachable, and no caching layer (a CGI
  script re-executed fresh per request has nothing worth memoizing across
  one page load). Falls back to a "Copy link" button with a clipboard-API
  script when unreachable; renders nothing when the URL field is NULL.
- [ ] **Sprite-based graphical system view**: render a system's star,
  planets, and moons as small icon sprites sized relative to each other
  (from `radius_km`) for an at-a-glance size comparison. Open questions:
  where the sprite art comes from, linear vs. logarithmic size scaling (a
  gas giant vs. a moon differ by 2-3 orders of magnitude in radius), and
  whether this is a simple size-comparison row or a full scaled-orbit
  diagram.

### Deployment history (interim `html/` browser)

Several deployment bugs surfaced on first production rollout to the
Ubuntu/Apache2 VPS (`starmap.moltenaether.com`) and are now resolved by
[`install.sh`](../install.sh)/[`update.sh`](../update.sh) (repo root):
CGI scripts deployed non-executable (git's `core.fileMode=false` silently
drops the executable bit — both scripts now `chmod +x` unconditionally on
every run, regardless of what mode git stored), CRLF line endings
(`.gitattributes` now pins `html/**/*.py`/`apache/*.sh` to `text eol=lf`),
`nltk`'s corpus download failing under the `www-data` user's unwritable
home directory (`install.sh` now pre-fetches the corpus system-wide;
`names.py` checks `nltk.data.find` before ever attempting a download), and
two rounds of `setup.py install` breaking against old apt-provided
`importlib_metadata`/`packaging` shadowing setuptools' own vendored
copies (resolved by dropping the deprecated `setup.py install` entirely in
favor of a build-isolated `pip install --upgrade --force-reinstall`).
Full root-cause detail for each is in the git history around
`install.sh`/`update.sh`/`apache/set-permissions.sh`.

## Future ideas — not scheduled, just parking so it isn't lost

Both items formerly here (the physical-plausibility anomaly finder, and the
evolved-star pre-Big-Bang mass-sampling issue) moved to "Unmerged work"
above once work on them actually started this session — see that section
for current status instead of here.