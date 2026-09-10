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

## Completed work log (through 2026-09-09)

Pointer index only — full rationale/detail for each is in `CHANGELOG.md`
and git history.

- **MySQL migration (Phase 5) and Flask API pagination/validation/health
  check/JSON error handling** — CHANGELOG [5.5.0]; see "Phase 5" below
  and `docs/api.md`/`docs/database-schema.md`.
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
- [x] Render an image or web interface to visualize the location of 2 points in galactic space -- see TODO in src/api/routes.py near `systems_near`. Done via the NAV feature (CHANGELOG.md [5.8.0]): `GET /api/nav`/`src/html/nav.py` give course/distance/route between two systems. The rendered-image gap that first pass left open is closed too (CHANGELOG.md [5.8.1]): `src/html/lib/navmap.py`'s "NAV Map" panel plots the origin, destination, and route hops as a flat, top-down SVG in the galactic X-Y plane (deliberately blind to altitude, same as the Galaxy Map's Quadrant view -- the course panel's own Altitude figure already covers that axis).
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
- [x] **Sector prism vertices with exact per-shell Voronoi tessellation —
  implemented and merged.** Every galaxy-placed sector's vertices
  (`sector_vertices` table, v7 -- one row per vertex, no JSON blobs) are
  now built from an exact local spherical Voronoi cell among its
  same-shell neighbors (genuinely gap-free
  laterally, not just reduced -- vertex count varies per sector, typically
  5-7, since a fixed-shape cube can't reconcile more neighbors than it has
  faces), extruded radially between the shell's inner/outer bounding
  spheres (area-matched, not vertex-matched, against adjacent shells).
  Same-shell neighbor search exploits this placement's Fibonacci-lattice
  structure for ~1000x speedup on large outer shells versus a naive radius
  search. Supersedes an earlier, never-released fixed-8-vertex corner-
  relaxation approach. See `docs/design/galaxy-coordinate-system.md`
  section 9.
- [x] **Galaxy thickness / shape (disk-density envelope) — implemented,
  as a compact "skeleton" rather than the originally-designed
  `galaxy_sector_plan` table.** The Milky-Way-scale exponential-disk-plus-
  bulge-plus-spiral-arm density model (`stellarObjects/galaxyDensity.py`,
  `relative_density`/`predicted_star_count`) is implemented and tested.
  `docs/design/galaxy-disk-density.md`'s revision-2 plan -- persisting
  every one of ~10.5 billion qualifying sectors' position/density up front
  in a `galaxy_sector_plan` table (~1 TB) -- was superseded before being
  built: a sector's position, density, and vertices are all pure
  deterministic functions of its `(shell_index, shell_slot_index)` address
  and a handful of galaxy-wide shape parameters, so none of that needs
  storing per sector at all -- it's cheaper to recompute on demand
  (sub-millisecond) than to look up. What's actually implemented instead
  (`stellarObjects/galaxySkeleton.py`, `galaxyPlan.py`, schema v8's
  `galaxy_shape`/`galaxy_shell_band`): one singleton row for the galaxy's
  shape parameters, plus one row per shell (almost always exactly one)
  recording the *candidate* slot-index band that shell's qualifying
  sectors could fall in -- a safe superset (an exact upper bound over
  every possible spiral-arm azimuth), not a per-sector list. Built in
  parallel across shells (`multiprocessing.Pool`, since each shell's
  band-finding is independent) by `galaxyPlan.py`; a full real-Milky-Way-
  scale build takes well under a second of actual compute and produces
  roughly 4,000 rows (~350 KB total), a many-orders-of-magnitude reduction
  from the ~1 TB `galaxy_sector_plan` design would have needed for the
  same information. Individual sectors are never batch-generated from the
  skeleton -- `galaxyGen.ensure_sector_generated(shell_index,
  shell_slot_index)` is the actual per-address entry point, called lazily
  the moment an address is visited: it consults the stored skeleton to
  decide, cheaply and exactly, whether that address holds anything at all,
  and if so generates and persists it on the spot (with that position's
  own `relative_density` driving the actual system count, not a uniform
  default) -- most of the galaxy is never visited, so it's never
  generated. A `sectors.UNIQUE (shell_index, shell_slot_index)` constraint
  (schema v8) turns a concurrent visit to the same never-before-generated
  address into a recoverable `IntegrityError` rather than a duplicate row.
  See `docs/design/galaxy-coordinate-system.md` section 9's storage-
  analysis addendum for the full reasoning and the real numbers measured
  against a real Milky-Way-scale build.

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
- [x] **Moved off SQLite to MySQL** for real concurrent multi-user access
  — implemented in full: `schema.sql` ported to MySQL/InnoDB DDL (a real
  `schema_migrations` table replacing `PRAGMA user_version`, explicit
  type mapping, every index/FK declared inline per table); `_db.py` now
  uses `pymysql` (pure Python, no build-time system libraries needed on
  the deployment host — chosen over `mysqlclient` for that reason) behind
  a small `Connection` compatibility wrapper, with real pooling via
  `DBUtils.PooledDB`; every entry point (`sectorGen.py`/`systemGen.py`/
  `galaxyGen.py`/`queryDb.py`/`migrateDb.py`/`src/api/`/`src/html/`) takes
  `--mysql-*`/`PLANETGEN_MYSQL_*` connection config instead of a SQLite
  path; a one-time `src/migrateSqliteToMysql.py` imports an existing
  pre-port database. See `docs/database-schema.md` and
  `schema.sql`'s "MySQL port" header note for the full rationale;
  CHANGELOG [5.5.0].

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
- [x] **Sector Map: real on-shell wedge shape, compass arrow, scale
  bar** -- done. A galaxy-placed sector (`shell_index`/`shell_slot_index`
  set) now draws its actual approximate on-shell wedge cell (new
  `sector_wedge_vertices_pc` in `src/stellarObjects/galaxyGeometry.py`)
  as a 12-edge wireframe instead of a generic cube; a sector with no
  placement still gets the plain cube. Also added a "Galactic Center"
  compass arrow (computed from the sector's own `center_x/y/z_pc`) and a
  live distance scale bar. Both the wedge and the compass arrow assume
  the sector's own local (x, y, z) axes run parallel to the galaxy
  frame's axes -- `galaxyGen.py` never actually rotates a sector's local
  star positions to the "Cube orientation" convention
  docs/design/galaxy-coordinate-system.md describes (that section only
  ever proposes it as a default, never wires it up), so this is the only
  assumption consistent with how every other galaxy-frame quantity on
  this map already renders. Implementing that orientation convention for
  real (rotating stored positions, or rotating only at render time) would
  let both drop this assumption -- open follow-up, not done here.
  See CHANGELOG.md [5.4.7].

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
