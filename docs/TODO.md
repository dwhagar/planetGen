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
  Apache2 (`examples/apache/`, see `docs/apache-deployment.md`). The read
  side is done; the API's write endpoints are still stubs — see "Open
  items" below.

## Open items

### Simulation / world-generation design questions

- [ ] Class K (Mars analog) and, by construction, every other ecosphere-zone
  class are generated at the same zone-midpoint orbital distance as Class M
  -- this generator doesn't place different terrestrial classes at
  different distances within (or beyond) the habitable zone the way real
  Mars sits much farther from the Sun than Earth. K's tuned values get as
  close to real Mars' absolute temperature/pressure as achievable under
  that constraint (~231K/~0.57kPa vs real ~210K/~610Pa) but can't fully
  close the gap without a zone/distance-placement change, which is a larger
  design question than per-class range tuning.

### Search

- [ ] Search parameter for searching by not only planet class but planet size, or sort by size in the tagged search field -- see TODO in src/queryDb.py near `process_args`/`_search_result_planets` and src/html/search.py near `_planets_panel`.

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
- [ ] **`GET /api/health` can crash instead of returning `503`.**
  `routes.get_db()` -> `queryDb.open_readonly` -> `_db.get_connection`
  raises a bare `SystemExit` (not `Exception`) when the configured MySQL
  server is unreachable -- `health()`'s own `except Exception` doesn't
  catch it, so the exact "database unreachable" case this endpoint exists
  to report as a clean `503` instead propagates out of the request
  entirely (confirmed via `python -c` against an unreachable MySQL host
  while building the admin-auth/write-API work below; pre-existing, not
  introduced by that work, and not fixed here since it touches the
  read-only path this pass otherwise left alone). Likely fix:
  `open_readonly`'s `SystemExit` was designed for `queryDb.py`'s CLI
  use, not for reuse inside a long-running Flask process -- `get_db()`
  probably wants its own `except (Exception, SystemExit)` (or a version
  of `open_readonly` that raises an ordinary exception instead).
- [ ] Sector Map's on-shell wedge shape and "Galactic Center" compass arrow
  (`docs/html-interface.md`'s `starmap.py` entry, CHANGELOG [5.4.7]) both
  assume a sector's own local (x, y, z) axes run parallel to the galaxy
  frame's axes, since `galaxyGen.py` never actually rotates a sector's
  local star positions to the "Cube orientation" convention
  `docs/design/galaxy-coordinate-system.md` describes (that section only
  ever proposes it as a default, never wires it up). Implementing that
  orientation convention for real (rotating stored positions, or rotating
  only at render time) would let both drop this assumption.

## Population and Politics

Exploratory ideas, not yet scoped or designed:

- [ ] Need to start thinking about assigning government ownership to a particular star system in the database so that together the systems make territories that are mapped out in 3D space by the star systems.
- [ ] Worlds with life on them can be flagged for generated names of dominant races.
- [ ] Probably going to need a space fairing species database.
- [ ] Need to think about under-developed / older civilizations and the differences and how to store and present that data based on society age.

## Completed work log (through 2026-09-10)

Pointer index only — full rationale/detail for each is in `CHANGELOG.md`
and git history.

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
  CHANGELOG [5.11.1].
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
