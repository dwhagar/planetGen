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

## Phase 4 — Galaxy-scale generation

- [ ] Tooling to generate and store MANY sectors at once, eventually
  covering "every sector of every space in the galaxy" per the long-term
  vision — needs its own design pass (a galaxy-scale coordinate system
  beyond a single sector's local cube, how sectors tile/connect to each
  other, batch-generation performance at that scale).
- [ ] `galaxyGen.py` — a new top-level CLI script (sibling to
  `sectorGen.py`/`systemGen.py`) that generates and persists many sectors
  as one galaxy, reusing `sectorGen.py`'s own per-sector generation/save
  logic for each one. Gated on the galaxy-scale coordinate system above
  being decided first.

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
- [ ] **Wiki-URL reachability check + clipboard fallback**, on
  `html/system.py`: for each of `star_systems.mediawiki_url`/`wikijs_url`,
  check (server-side, short timeout, cached so a page load doesn't stall on
  live HTTP round-trips) whether the URL is set *and* resolves; if so,
  render it as a normal link, otherwise show a "Copy to clipboard" button
  next to the matching content instead. Open questions: how "does it
  exist" is checked (HEAD vs. GET, what counts as a hit for a wiki that
  200s its own placeholder page), and whether/how the check result is
  cached.
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

- [ ] **Physical-plausibility test suite (anomaly finder)**: a testing tool
  that builds "sane" reference ranges for pressure, temperature, gravity, and
  the other generated physical factors of a planet/star (by class/type), then
  generates a large batch of bodies and flags any whose values fall outside
  those ranges — a way to catch generator bugs like the atmospheric-pressure
  unit bug (fixed in [5.3.0] — see CHANGELOG.md) by statistical outlier
  detection rather than hand-noticing a bad value in one rendered page.
  Needs its own design pass later (where the "sane" ranges come from —
  hand-authored per class vs. derived from a large sample; how many bodies
  to batch-generate per run; whether this lives in `tests/` as a slow/opt-in
  suite or as a separate standalone script). Deliberately low priority.
- [ ] **Evolved-star mass sampling can imply a pre-Big-Bang star**: an
  evolved-class star (`Yerkes != V`) derives its required main-sequence
  lifespan from its own already-generated mass
  (`Star._calculate_initial_star_age_and_lifespan`'s evolved-star branch,
  `stellarObjects/starData.py`); for a sub-solar-mass progenitor (roughly
  under ~0.88 Msun), that implied main-sequence lifespan alone already
  exceeds `UNIVERSE_AGE_GY` (13.8 Gy) — meaning such a star couldn't
  actually have evolved off the main sequence yet in the real universe.
  The [5.3.0] universe-age fix deliberately does *not* paper over this by
  capping age below its own required floor (that would produce a
  self-contradictory star, e.g. a red giant younger than its own
  progenitor's main-sequence lifespan) — it's a deeper issue with how
  evolved-star masses are sampled in the first place, and needs its own
  design pass (e.g. resampling/rejecting masses whose implied
  main-sequence lifespan exceeds `UNIVERSE_AGE_GY` before generating an
  evolved-class star at all).
- [ ] **Class M/P forcing clamps**: commented out (not deleted) in
  [5.3.0]'s `stellarObjects/planetPhysics.py` now that the corrected
  atmospheric-pressure formula lands close to realistic ranges on its own.
  Revisit whether they're still needed at all, or should be restored in a
  narrower form, once more real-world output has been reviewed.

## Open questions still to resolve

Resolved (kept here only as a pointer, full reasoning lives where the
decision is used): schema versioning, `float('inf')` white-dwarf lifespan
storage, `reflection_spectrum_visible`/`non_visible` storage shape (a
normalized child table, no JSON columns anywhere in this schema),
`SQLAlchemy`/Alembic vs. raw `sqlite3` (raw `sqlite3`), secondary-star
`system_config` asymmetry (collapsed), and `BinaryStarProxy` derived-field
recompute-vs-snapshot (snapshot) — see `schema.sql` and git history for
where each is used.

- [x] By design, there is no seed-based replay anywhere in this plan —
  generation mixes the unseedable `secrets` module with the seedable
  `random` module (documented in `spaceSector.py`'s own docstring), so
  results are stored directly rather than via a seed.
- [ ] `examples/*.json` are unaffected by any of this — those are
  `systemGen.py --system-file` recipe *inputs*, a wholly separate format
  from the generated-*result* persistence this plan adds.
