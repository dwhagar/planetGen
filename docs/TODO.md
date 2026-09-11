# planetGen Roadmap: Full Database Storage + Web Interface

## Vision

Long-term goal: a fully populated galaxy (every sector, every star system,
every planet/moon/asteroid belt, with full generated detail preserved)
stored in a real relational database, served through a web interface for
browsing and search. Phases 0-5 below are all substantially complete —
what's left is a short list of open design questions and follow-ups (see
"Open items"), plus the open-ended "Population and Politics" ideas that
haven't had a design pass yet.

## How to use this document

Work top to bottom in "Open items" — that's the only section with anything
actionable. Completed work is condensed into short pointer bullets (full
historical rationale lives in git history and `CHANGELOG.md`, not here) —
read `src/stellarObjects/schema.sql` and `docs/database-schema.md` for the
actual current schema, and `docs/html-interface.md`/`docs/api.md` for the
actual current web interface, rather than this document. Trim finished
items down to a pointer whenever a section is revisited; only genuinely
open items need working detail.

## Phases 0-5 — Complete

- **Phase 0 (flavor text bug)**: flavor text is decided once at generation
  time instead of being re-rolled at render time.
- **Naming ASCII audit**: every string constant in `stellarObjects/names.py`
  is verified 7-bit-ASCII-printable, with a regression test guarding future
  additions.
- **Phase 1 (object-graph serialization)**: every generated class has a
  `to_dict`/`from_dict` pair (`stellarObjects/serialization.py`);
  reconstruction bypasses `__init__` so loading never re-runs random
  generation.
- **Phase 2 (relational database)**: a normalized schema
  (`stellarObjects/schema.sql`) and persistence layer (`stellarObjects/_db.py`)
  with versioned migrations — see `docs/database-schema.md` for schema
  history.
- **Phase 3 (CLI tools use the database)**: `sectorGen.py`/`systemGen.py`
  save every run to the database unconditionally; `src/queryDb.py` is a
  read-only search CLI.
- **Phase 4 (galaxy-scale generation)**: sectors placed in real 3D
  galaxy-frame space — deterministic Fibonacci-sphere shell placement,
  exact per-shell Voronoi sector vertices, and a compact density
  "skeleton" (`galaxyPlan.py`, `galaxy_shape`/`galaxy_shell_band`) that
  drives lazy per-address generation instead of pre-planning ~10.5 billion
  sectors. See `docs/design/galaxy-coordinate-system.md` (sections 8-9),
  `stellarObjects/galaxyGeometry.py`/`galaxySkeleton.py`/`galaxyDensity.py`,
  and `galaxyGen.py`/`galaxyPlan.py`.
- **Phase 5 (web interface + MySQL backend)**: moved off SQLite to MySQL
  for real concurrent multi-user access (`schema.sql` ported to
  MySQL/InnoDB, `_db.py` on `pymysql`/`DBUtils.PooledDB`) — CHANGELOG
  [5.5.0]. A Flask JSON API (`src/html/api/`, see `docs/api.md`) backs a
  thin CGI browser (`src/html/`, see `docs/html-interface.md`) covering
  database/sector/system browsing, an interactive 3D Sector Map, a scaled
  System Map orbit diagram, a galaxy-scale Galaxy Map, faceted search, and
  NAV (course/distance/route between two systems) — CHANGELOG [5.8.0],
  [5.8.1], [5.10.0], [5.12.0]. Deployed via `install.sh`/`update.sh` to
  Apache2 (`examples/apache/`, see `docs/apache-deployment.md`). Both the
  read side and the write endpoints (admin-authenticated `POST`/`PATCH`/
  `DELETE` on sectors/systems) are real now — see "Open items" below for
  what those write endpoints still don't cover.

## Open items

### Search

- [ ] `queryDb.py`'s CLI (`sectors`/`systems`/`near` subcommands) still has
  no `planets`/`moons` equivalent at all -- `systems` only filters by
  `star_type_prefix`/`sector_id`, so "every Class D planet smaller than
  Earth" can't be asked of this particular tool (the web/API faceted
  search now supports a planet/moon/star size range plus class/body/life
  tags -- CHANGELOG [5.23.0] -- this item is specifically about the
  separate `queryDb.py` CLI never getting its own `planets` subcommand).
  See the TODO comment in `src/queryDb.py` near `process_args` for the
  concrete shape (`--class`, `--min-radius-km`/`--max-radius-km`) this
  would need.

### Web API/frontend

- [ ] **Sector-attached system creation isn't supported via the API.**
  `POST /api/systems` (real now — see `docs/api.md`'s "Systems — request
  body") only ever creates a **standalone** system (`sector_id = NULL`).
  Attaching a newly generated system to an existing sector needs that
  sector's own placement/Hill-sphere separation logic
  (`SpaceSector.add_system`), deliberately not wired into the write API
  in this pass (kept the validation/generation surface smaller for the
  admin-auth work that landed alongside it — see `docs/api.md`).
- [ ] **Editing a system's generated content isn't supported via the
  API.** `PATCH /api/systems/<id>` only renames a system — no way to
  modify its stars/planets/moons/belts short of `DELETE` + `POST`
  (regenerate). Not clear this needs solving at all (vs. "just
  regenerate"), but flagged in case it does.

## Population and Politics

Exploratory ideas, not yet scoped or designed:

- [ ] Need to start thinking about assigning government ownership to a particular star system in the database so that together the systems make territories that are mapped out in 3D space by the star systems.
- [ ] Worlds with life on them can be flagged for generated names of dominant races.
- [ ] Probably going to need a space fairing species database.
- [ ] Need to think about under-developed / older civilizations and the differences and how to store and present that data based on society age.

## Completed work log (through 2026-09-11)

Pointer index only — full rationale/detail for each is in `CHANGELOG.md`
and git history.

- **`BINARY_SYSTEM=None` now rolls real stellar-multiplicity chance
  instead of always producing a single star**, keyed by the primary's
  own spectral letter (Duchene & Kraus 2013; Raghavan et al. 2010; Moe &
  Di Stefano 2017 — `program_constants.
  BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS`,
  `StarSystem._should_generate_binary`) — the same tri-state contract
  (`True`/`False` force, `None` rolls chance) every other `SystemConfig`
  flag already follows — CHANGELOG [5.25.0].
- **`estimate_num_objects` now derives its planet/belt ceiling from real
  protoplanetary-disk physics** (Minimum Mass Solar Nebula surface
  density, Hayashi 1981; oligarchic-growth isolation mass, Lissauer 1993/
  Kokubo & Ida 2000-2002; real disk-mass-vs-stellar-mass scaling,
  Andrews 2013/Pascucci 2016) instead of an arbitrary curve fit to
  stellar mass, and reuses the exact mutual-Hill-radius spacing rule
  (`MUTUAL_HILL_RADII_SEPARATION`, [5.23.0]) that will later constrain
  actual placement, so the two are provably consistent
  (`StarSystem._estimate_max_objects_from_disk_physics`,
  `utils.snow_line_au`/`disk_surface_density_scale`/
  `mmsn_surface_density_gcm2`/`isolation_mass_kg`) — CHANGELOG [5.24.0].
- **Orbital spacing now uses the mutual Hill radius, not either planet's
  own individual one.** Real stability criteria (Gladman 1993; Chambers,
  Wetherill & Boslough 1996; Smith & Lissauer 2009) express minimum
  planet-planet separation in units of the *pair's* combined mass and
  average distance, not one body's own Hill sphere -- `StarSystem.
  validate_system`'s planet-planet spacing check now does the same
  (`_mutual_min_distance_au`, `utils.mutual_hill_radius_m`,
  `program_constants.MUTUAL_HILL_RADII_SEPARATION = 10`), solved in
  closed form for the exact minimum distance rather than approximated at
  the pre-correction position (the naive version measurably undershot,
  since the mutual radius depends on the average of both distances and
  grows once the correction moves one of them) — CHANGELOG [5.23.0].
- **`GET /api/health` no longer crashes on an unreachable database.**
  `routes.get_db()` now catches `open_readonly`'s `SystemExit` and
  re-raises it as an `ApiError` (503) -- fixes `/health` itself and,
  since every other read route shares `get_db()`, every one of them too
  (previously an uncaught `SystemExit`, a `BaseException`, would
  propagate straight through Flask's dispatch instead of becoming any
  HTTP response) — CHANGELOG [5.23.0].
- **Faceted search gained a min/max size (`radius_km`) range for
  stars/planets/moons.** `queryDb.search`'s new `sizes` argument (`GET
  /api/search`'s `star_min_radius_km`/`_max_radius_km`, likewise
  `planet_`/`moon_`) filters alongside the existing class/body/life tags;
  `html/search.py` gained matching form fields, active-filter chips, and
  a Radius column on the Stars/Planets/Moons result panels — CHANGELOG
  [5.23.0]. `queryDb.py`'s own CLI still has no `planets`/`moons`
  subcommand at all -- see "Open items" > "Search" above, a narrower,
  still-open item.
- **Sector Map star dots now rotate into the galaxy frame at render
  time.** The wedge outline and "Galactic Center" compass arrow were
  always computed directly from galaxy-frame quantities and so were
  always correct; star dots (`star_systems.position_x/y/z_mpc`) were
  plotted as if their own sector-local axes already ran parallel to the
  galaxy frame's, which generation never actually guarantees.
  `starmap.py`'s new `_rotate_to_galaxy_frame` applies the "Cube
  orientation" convention (`docs/design/galaxy-coordinate-system.md`
  section 3: radial-outward local `+Z`, projected-galactic-north local
  `+X`) from the sector's own stored `center_x/y/z_pc` alone, reusing
  `stellarObjects.sectorGeometry.cube_orientation` -- the same basis that
  module already computes for this sector's own wedge vertices -- so no
  new stored orientation column was needed — CHANGELOG [5.23.0].
- `validate_system`'s orbital-overlap correction could strand a planet
  outside the zone its class needs (e.g. an "Earth-like" Class M at
  ~30K, hundreds of thousands of AU out — 10% of a 400-system sample
  affected, almost entirely O/B-type stars), and the sequential
  placement loop had no outer bound to begin with. Fixed via
  `planetPhysics.reconcile_zone_and_class` (re-derives zone/class after
  any post-hoc distance correction, for a planet and its moons), a
  `star.system_perimeter` cap on the placement loop
  (`StarSystem._generate_planets`), and a generation-retry loop
  (`MAX_SYSTEM_GENERATION_ATTEMPTS`) so reconciliation can't silently
  break a requested `HABITABLE_WORLD`/`ASTEROID_BELT` guarantee —
  CHANGELOG [5.22.0].
- Ecosphere-zone classes (E/F/G/H/K/L/M/N/O/P/V) now generate at a
  class-appropriate distance within the zone via per-class
  `zone_position_mode`, instead of every class sharing the same
  distance-blind draw; K (Mars) and N (Venus) retuned to real-world
  targets now that real insolation does the work instead of a
  compensating albedo/greenhouse hack — CHANGELOG [5.21.0].
- **Write-capable API + admin auth.** `POST`/`PATCH`/`DELETE` on
  `/api/sectors`/`/api/systems` do real inserts/updates/deletes now,
  gated behind admin login (session cookie, `HttpOnly`/`Secure`/
  `SameSite=Strict`) or an API key (`Authorization: Bearer`) —
  `stellarObjects/adminAuth.py` + a new, deployment-global control schema
  (`stellarObjects/control_schema.sql`: `admin_users`/`admin_sessions`/
  `admin_api_keys`/`admin_audit_log`), seeded with a default `admin`/
  `password` login that's blocked from doing anything else until its
  credentials are changed. Writes run against a separate,
  less-privileged `PLANETGEN_MYSQL_WRITE_*` account, never the
  `SELECT`-only one every read endpoint uses. New admin web pages
  (`html/login.py`/`changecreds.py`/`admin.py`) for logging in, the
  forced credential change, and API key management. See `docs/api.md`'s
  "Authentication"/"Write endpoints" sections,
  `docs/database-schema.md`'s "The control schema", and
  `docs/apache-deployment.md`'s "MySQL accounts". Deliberately **not**
  included (per explicit direction): any database-file-management
  surface (create/rename/duplicate/delete a whole MySQL schema) — this
  deployment doesn't expose that over the Internet.
- Backend API (Flask) mounted at `/api/` alongside the interim `html/`
  browser, with pagination/validation/health-check/JSON-error-handling —
  CHANGELOG [5.5.0]; see `docs/api.md`.
- Evolved-star pre-Big-Bang mass fix via reject-and-resample
  (`starData.py`) — CHANGELOG [5.3.1].
- Habitable/atmosphere sanity review
  (`docs/analysis/habitability-atmosphere-sanity-review.md`) found Class
  M/P statistically indistinguishable, an inverted greenhouse-factor
  formula, pressure mathematically decoupled from gravity, and
  implausibly-low gas-giant densities — all fixed via the
  physical-plausibility anomaly finder (`stellarObjects/plausibility.py`,
  `src/tests/physical_plausibility_cli.py`) and Track A physics fixes —
  CHANGELOG [5.3.5].
- Greenhouse-formula fix + per-class climate tuning for
  M/O/H/K/L/N/E/F/G/V/P (`src/tests/climate_tuning_cli.py`) — CHANGELOG
  [5.3.7]/[5.9.1].
- Life-bearing classes (Q, formerly also W) restricted to the ecosphere
  zone only — CHANGELOG [5.3.8].
- Class cleanup: gas giants excluded from zones `h`/`e`; classes R
  (unreachable), S/U (folded into gas-giant classes), W (its day/night
  identity needed a dayside/nightside model out of scope for this
  generator), X/Y (folded into B/A) all removed or merged, leaving 19
  classes (A-Q, T, V); per-class `size_mode` bell-curve radius
  distributions added — CHANGELOG [5.3.9], [5.4.0], [5.4.2], [5.9.0].
- `sectorGen.py --density`: controllable sector density, shared with
  `galaxyGen.py` — CHANGELOG [5.4.1].
- System-page TOC made collapsible, then moved to a fixed right-margin
  rail above `min-width: 90rem` — CHANGELOG [5.4.2]/[5.4.3].
- NAV feature: `GET /api/nav`/`src/html/nav.py` give course/distance/route
  between two systems, with `src/html/lib/navmap.py`'s flat top-down SVG
  "NAV Map" plotting origin/destination/route hops — CHANGELOG [5.8.0]/[5.8.1].
- Sector Map: interactive 3D (drag-to-rotate, scroll-to-zoom, CSS 3D cube
  + billboarded star dots), then real on-shell wedge shape/compass
  arrow/scale bar for galaxy-placed sectors — CHANGELOG [5.4.4]/[5.4.7]
  (wedge/compass orientation-convention assumption tracked as an open item
  above).
- `src/api/` moved to `src/html/api/` (with `src/wsgi.py` ->
  `src/html/wsgi.py`), served from the same tree as the CGI browser via
  its own `WSGIDaemonProcess` — CHANGELOG [5.10.0].
- Subdwarf (Yerkes VI) age modeling given its own age-generation branch
  (binary mass-stripping origin, not single-star progenitor lifespan) —
  CHANGELOG [5.10.1].
- Orbital and rotational motion: every planet/moon has a fixed 3D orbital
  orientation and rotation period, plus a live `orbital_phase_deg`
  advanced by the new `src/updateOrbits.py` script (meant to run
  periodically, e.g. via cron) — schema v8 -> v9, CHANGELOG [5.11.0].
  Moon tidal locking replaced with a real tidal-despinning timescale
  estimate, plus a log-uniform (not linear-uniform) moon distance draw —
  CHANGELOG [5.11.1]. Extended in the same vein since: galactic orbital
  speed/period for every star system (schema v10, CHANGELOG [5.13.0]);
  real 3D Cartesian position + orbital speed for every planet/moon
  (schema v11, [5.14.0]); a floating-point update guard so
  `advance_orbital_phases` skips a body once `elapsed_years` can no
  longer move its phase (schema v12, [5.15.0]); stars themselves now
  advance a galactic orbital phase, and binary pairs get a real mutual
  orbit around their barycenter (schema v13, [5.16.0]), including that
  mutual orbit's own Cartesian position (schema v14, [5.17.0]) — all
  advanced by the same `updateOrbits.py` run. `examples/maintenance/`
  gained systemd timer units (`planetgen-orbits@.timer`, per-database
  template instance) as an Ubuntu/Debian-native alternative to a raw
  crontab line for running `updateOrbits.py` — [5.18.0] — plus a second
  timer (`planetgen-update.timer`) that runs `update.sh` on a schedule
  30 minutes ahead of it so code updates land before the orbit run picks
  them up, with `update.sh` itself now skipping its full reinstall step
  when `git pull` found nothing new — [5.20.0].
- `../src/html/` rearchitected as a thin frontend for the Flask API
  instead of a direct MySQL client (`html/lib/apiclient.py`,
  `PLANETGEN_API_BASE_URL`); API gained `?db=` on every route,
  `GET /api/databases`, `GET /api/galaxy/sectors`, `GET /api/search`, and
  richer sector/system detail responses. System Map gained a green
  "supports life" badge — CHANGELOG [5.12.0].
- Schema-migration backups are gzip-compressed and excluded from the web
  database picker and from a subsequent migration run
  (`stellarObjects._db.migrate_database`'s `BACKUP_MARKER`).
- Deployment bugs from first production rollout (CGI scripts deployed
  non-executable, CRLF line endings, `nltk` corpus download failing under
  `www-data`'s unwritable home directory, `setup.py install` breaking
  against old apt-provided packaging) resolved by
  `install.sh`/`update.sh`/`examples/apache/set-permissions.sh`.
