# Library migration and code layout

The plan of record for moving planetGen off its zero-dependency model and
onto third-party open-source libraries, and for reorganizing the code into
importable packages. It condenses Boss's research documents in this folder
("Library Migration Workflow.md", "Web UX Development Notes.md" and "Web UX
and Job Management Guide.md") and records where this plan departs from them.

## 1. What Boss asked for

Boss (2026-10-03 05:38Z): "Transition the code to open source 3rd party
libraries and simplify the code base, also in this step we'll use a redis
server for use in managing and executing the work queue. ... To be clear,
we have an explicit directive to move the code base from a 0-dependency
model into using 3rd party open source libraries to simplify our own code
deployment."

And: "Code reorganize to put things into distinct python modules as
different importable objects for modularity of code and easier
maintenance. Maximize the use of utility libraries for shared functions
across all modules."

Both are phase 0 groundwork: several bugs (the job log bugs, the 2FA and
rate-limit test flakes, the menu and button bugs, slow map zoom) are fixed
by the new libraries rather than patched in the old code.

## 2. Order

1. **Pins and Redis first.** The new packages go into the hash-pinned lock
   files, and install, update and CI provide a Redis server.
2. **Layout plan, then the move.** A package layout is written down and
   approved, then the code moves one package per PR with no behavior
   change. Doing this before the swaps means each file moves once.
3. **Library swaps**, each its own PR, deleting the code it replaces:

| Area | Today | Library | Notes |
|---|---|---|---|
| Two-step sign-in | `totp.py` (100 lines), `qrcodegen.py` (900 lines) | pyotp, segno | Same secrets, step and window, so enrolled users keep working. segno renders SVG with no further dependencies. |
| Markdown | `mdconvert.py` | markdown | Same output for pages and the wiki export. |
| Rate limits | `loginThrottle.py`, `api/limiter.py`, `api/loginguard.py` | Flask-Limiter (Redis storage) | The lockout rules in login-brute-force-protection.md are kept. |
| Work queue and web jobs | `workQueue.py` (1,600 lines), `planetgen.cli.job`, the `work_lease` table | RQ on Redis | See section 3. |
| Caches | `pagecache.py`, `tilecache.py` | cachetools for `pagecache.py` only (done: a `TTLCache` sized by body length, `max_entries` held on top; PERF.25) | Boss agreed (2026-10-07 13:27Z): `tilecache.py` stays (JSON only, never unpickles, prunes by size, checks its folder is private). Out: diskcache and sqlitedict (unfixed advisories PYSEC-2026-2447 and PYSEC-2026-1939, rejected by pip-audit), cachelib and Flask-Caching (pickle by default, count items not bytes). Fallback if disk caching proves slow: an optional Redis backend through redis-py storing JSON, on its own instance or with TTL'd keys so tile eviction can't evict rate-limit counters. |
| Database | `store.py` (9,800 lines), `planetgen.cli.migrate` with `schema_vNN.sql.gz` fixtures | SQLAlchemy, Alembic | Alembic starts from a baseline that recognises existing databases at the current schema. |
| Validation | `validation.py` | Pydantic | The same limits; errors list every field. |
| Physics | `keplerMotion.py`, `physical_constants.py` | scipy, astropy | Results checked against today's within stated tolerances. |
| Job logs and progress | `generatejobs.js`, `generatefolds.js`, `jobs.py`, `progressRate.py` | Xterm.js over Server-Sent Events, native `<progress>` | See section 4. |
| Buttons, menus, dialogs | hand-written HTML and CSS | Shoelace web components | Uses the approved icon set. |
| Tables | `tabledisplay.py` | TanStack Table and TanStack Virtual | The server still pages 50 rows. |
| Galaxy Map picking and streaming | `galaxymap3d.js`, `galaxyprisms.js` | three-mesh-bvh, 3D tiles, camera-relative rendering | Level of detail per tile. |
| Nebula meshes | none | scikit-image (marching cubes) | For the nebula shape (nebula-and-asteroid-field-classes.md). |

## 3. The work queue: Redis and RQ

The job guide proposes a broker-less `ProcessPoolExecutor` and mentions
SQLite. Boss chose Redis for the queue, and the store is MySQL or MariaDB,
so this plan uses **RQ on Redis**:

- Generation units (a sector, a layer scatter, a backfill block) and web
  jobs (Generate page runs, map-menu fills, API jobs) are RQ jobs.
- Worker count, retries, timeouts and status come from RQ. A worker that
  dies has its job retried instead of hanging the run.
- Each job publishes its progress and log lines to Redis, which the web
  log window and progress bars read.
- Results must match the old queue for the same seed at any worker count
  (reproducible-galaxies.md, section 5).
- The jobs folder keeps only downloadable logs.

Huey and Celery were looked at: Huey would avoid Redis but Boss chose
Redis; Celery is heavier than a single-host deployment needs.

**Windows**: Redis has no supported native Windows build. The installer
points at Redis in WSL2 (Boss, 2026-10-07 17:11Z: "Let's say Redis in WSL"; OPS.27).

## 4. Logs and progress in the browser

A running job's log streams over Server-Sent Events into an Xterm.js
terminal, and its progress into native `<progress>` bars. Under Apache,
SSE needs threaded mod_wsgi daemon processes (or a gevent worker behind a
proxy) so one long-lived stream doesn't hold a whole process; the
deployment docs give the settings. Log messages are single lines (the
terminal wraps them), a failed job's log stays open until the user
continues, and tracebacks appear in full with a Copy button.

## 5. Front-end libraries

Shoelace, TanStack and three-mesh-bvh are vendored as ES module builds and
served by Flask, with no bundler and no CDN at runtime (default until Boss
decides otherwise). Both themes, keyboard use and reduced motion are kept.

## 6. Package layout (OPS.23)

The plan for OPS.24's move. Boss approves it before the move starts.

### 6.1 Where things live today

- `src/stellarObjects/`: 75 flat modules, 46,000 lines, among them
  `store.py` (9,800), `program_constants.py` (2,900), `systemData.py`
  (2,100) and `utils.py` (1,900, a mix of formatting, unit conversion,
  random sampling, orbital formulas and word-salad names).
- `generate.py` (4,600 lines) at the repo root, and nine scripts in
  `src/` (`planetgen.db.query` alone is 4,300).
- `src/html/`: `api/` and `web/` import as top-level packages and the 17
  modules in `lib/` as top-level modules, because `wsgi.py`,
  `web/__init__.py`, `routes.py` and `apiclient.py` each push a folder
  onto `sys.path`. `fmt.py` repeats `utils.py`'s distance formatting and
  `galaxymap.py` repeats `galaxyGeometry.py`'s quadrant and zone helpers.
- `src/wikiClient/`: the wiki export client.

### 6.2 No compatibility layer

Boss (2026-10-07 13:04Z): "In the end I don't want shims or wrappers or
such, so far no one uses this but me, so I don't want to worry yet about
backward compatibility." So:

- A moved module leaves nothing at its old path. Each move PR changes
  every caller in the same PR: code, tests, scripts, install and update,
  CI, the Generate page's job commands and the docs.
- The command-line scripts don't stay at their old paths as wrappers.
  Each becomes a module in `planetgen.cli` run as
  `python3 -m planetgen.cli.<name>` (and `planetgen` stays the console
  command for generation).
- install and update install the checkout as an editable package
  (`pip install -e .`), so `planetgen` imports from anywhere with no
  `sys.path` pushes, and a `git pull` takes effect without a reinstall.
- `stellarObjects`, `src/html/api`, `src/html/web`, `src/html/lib`,
  `src/wikiClient`, `generate.py` and the `src/*.py` scripts are gone
  when the move is done.
- `src/html/` keeps only `wsgi.py` (the WSGI entry point Apache, gunicorn
  and waitress load) and `static/`, so the server's Apache config needs
  no change. wsgi.py imports `planetgen.web.app` and nothing else.

A branch open during the move (the Bugfixes lane's) merges main after
each move PR. Git follows a renamed file, so an edit to a moved module
usually merges cleanly; an import of an old path fails loudly and is
fixed to the new one.

### 6.3 The packages

One top-level package, `planetgen`, in `src/planetgen/`. Module names
become `snake_case`. Every package's `__init__.py` stays empty (no
imports), so no package pulls in another just by being imported and no
import cycle can form between them.

| Package | Modules (today's name in brackets) |
|---|---|
| `planetgen.util` | `log` (log), `appconfig` (appconfig), `serialization` (serialization), `checks` (`finite_domain` from utils), `format` (number, distance, speed, duration, temperature, pressure and age text from utils, merged with the same helpers in html/lib/fmt), `random` (seeded sampling: power law, bounded bell, the five copies of `_log_uniform`) |
| `planetgen.tuning` | one module (program_constants): every package reads it, so it sits outside them all |
| `planetgen.physics` | `constants` (physical_constants), `units` (ly, pc, mpc and AU conversions from utils), `orbits` (habitable zone, Hill sphere, Holman-Wiegert, reflex offset, orbital position and update interval from utils), `formation` (snow line, MMSN, isolation mass from utils), `kepler` (keplerMotion), `planets` (planetPhysics), `rogue_surface` (rogueSurface), `stellar_evolution` (stellarEvolution), `mathcheck` (mathCheck) |
| `planetgen.galaxy` | `geometry` (galaxyGeometry), `density` (galaxyDensity), `skeleton` (galaxySkeleton), `drill` (galaxyDrill), `viewport` (galaxyViewport), `seed` (galaxySeed), `version_key` (versionKey), `sector` (spaceSector), `nebula_field` (nebulaField), `navigation` (navigation), `nav_graph` (navGraph), `galactic_orbit` (the galactic orbit helpers from utils) |
| `planetgen.names` | `wordlists` (names, with offensive_words.txt), `wordsalad` (the syllable and phoneme-salad helpers from utils), `uniqueness` (nameUniqueness), `bodies` (bodyNames), `object_id` (objectId). GEN.67 adds the codec here (today's root `gatedPhonemeCodec.py`); GEN.71 deletes `wordlists` and `wordsalad`. |
| `planetgen.generation` | `config` (config: SystemConfig), `star` (starData), `planet` (planetData), `system` (systemData), `binary` (doubleStar), `wide_binary` (wideBinary), `belt` (asteroidData), `comet` (cometData), `life` (planetLife), `evolution` (evolution), `star_population` (stellarPopulation), `bright_stars` (brightStars), `limits` (generationLimits), `stats` (generationStats), `validation` (validation), `plausibility` (plausibility), `phenomena_plausibility` (phenomenaPlausibility); and `run_system`, `run_sector`, `run_galaxy`, `run_plan`, `run_phenomenon`, `run_population` (generate.py's sections) |
| `planetgen.generation.phenomena` | `asteroid_field`, `compact_remnant`, `nebula`, `quasar`, `rogue` (rogue planets and interstellar comets), `supernova_remnant` (the six phenomenon `*Data.py` modules) |
| `planetgen.population` | `model` (population), `facilities` (facilities) |
| `planetgen.db` | `store` (_db, with schema.sql and control_schema.sql beside it), `edits` (editStore), `render` (systemRender), `query` (queryDb's queries), `stats` (adminStats) |
| `planetgen.admin` | `auth` (adminAuth, with the common-password list), `throttle` (loginThrottle), `totp` (totp), `qrcode` (qrcodegen), `activity_log` (activitylog), `edits` (adminEdits) |
| `planetgen.queue` | `work` (workQueue), `progress_file` (progressFile), `progress_rate` (progressRate), `load` (systemLoad) |
| `planetgen.api` | everything in `src/html/api/` under the same module names |
| `planetgen.web` | everything in `src/html/web/` under the same names, with `templates/`; `app` (the Flask app factory, from api/app); `planetgen.web.lib` for html/lib's shared modules (apiclient, fmt's HTML helpers, pagination, pagecache, tilecache, classref, tabledisplay, mdconvert, privatedir, systempage); `planetgen.web.maps` for the map renderers (starmap, systemmap, navmap, galaxymap, galaxymap3d, phenomenonmap, phenomenonrender) |
| `planetgen.wiki` | the wiki client (`src/wikiClient/`) |
| `planetgen.cli` | `generate` (generate.py's argument parsing and dispatch), `query` (queryDb's command line), `migrate` (migrateDb), `reset` (resetDb), `orbits` (updateOrbits), `dedupe` (dedupeNames), `lockouts` (loginLockouts), `render_parity` (checkRenderParity), `job` (jobRunner, moved whole; RQ replaces it, section 3). Each of the six small scripts moved whole, logic and argument parsing together: migrate's logic is replaced by Alembic (DB.11) and the rest are a page or two each, so splitting them out into `planetgen.db` buys nothing. |

Large modules move whole. `store.py` is split by DB.11 (SQLAlchemy), and
`planetgen.db.query` with it; splitting them during the move would make every
open branch conflict twice. `utils.py`, `fmt.py` and `generate.py` are
the exceptions: they are split as listed, because nothing replaces them
later.

`src/tests/` stays where it is; each move PR changes the tests' imports
and string patch targets (`monkeypatch.setattr("stellarObjects.utils.x",
...)`) along with everything else. `_version.py` moves to
`planetgen/_version.py` in the first PR, with `setup.py`,
`scripts/bump_version.py` and the release workflow.

### 6.4 Shared helpers merged into utilities

Only where the merged helper gives identical results, so a seeded galaxy
is unchanged (reproducible-galaxies.md):

- Number and distance formatting: `utils.format_distance_*` and
  `fmt.format_distance_*` into `planetgen.util.format`; `planetgen.web.lib.fmt`
  keeps only the HTML helpers.
- `sector_quadrant` and `sector_zone`: galaxymap calls
  `planetgen.galaxy.geometry`'s.
- `_log_uniform` (generate.py, compactRemnant, quasarData,
  stellarEvolution, rogueSurface): one `planetgen.util.random.log_uniform`
  where the draws match; a copy whose draw order differs stays local,
  with a note, until a seeded test proves the swap is safe. Done in step
  13: every copy made the same single `uniform` draw, so all of them
  (and nebula's `_draw_in_range` and the wide-binary separation) now call
  `log_uniform`; stellarEvolution's copy was unused and went.
- utils.py's helpers the table above doesn't name went to their one
  caller's module: the star-profile helpers (`get_star_spectral_class`,
  `get_star_evolutionary_profile`) to `planetgen.generation.star`,
  `calculate_object_mass` to `planetgen.physics.planets`, the wide-binary
  samplers to `planetgen.generation.wide_binary`, and `to_paragraph` and
  `properties_to_string` to `planetgen.util.format`.
- `sys.path` pushes: all four go (the editable install makes them
  unnecessary).

### 6.5 The order of the move (OPS.24)

One package per PR, each mechanical (`git mv`, import lines, callers),
each with the full suite green and no change to behavior. Leaves first,
so each PR's modules import ones that already moved:

1. Scaffold `src/planetgen/` with `_version.py` and `planetgen.util`
   (log, appconfig, serialization). Every entry point already puts `src/`
   on `sys.path`, so `planetgen` imports the way `stellarObjects` does;
   the editable install comes with step 14, when those `sys.path` pushes
   go.
2. `planetgen.tuning` and `planetgen.physics`.
3. `planetgen.galaxy`.
4. `planetgen.names`.
5. `planetgen.generation` and its phenomena.
6. `planetgen.population`.
7. `planetgen.admin`.
8. `planetgen.db`, with the `src/` scripts into `planetgen.cli`. Until
   step 14's editable install, they run as `python3 -m planetgen.cli.<name>`
   from the checkout's `src/`; install and update do that, and the job
   runner puts `src/` on its steps' `PYTHONPATH`.
9. `planetgen.queue`, with jobRunner into `planetgen.cli.job` (the web
   app starts it with `python -m` from `src/`).
10. `planetgen.web.lib` and `planetgen.web.maps`.
11. `planetgen.api`.
12. `planetgen.web` with its templates, and wsgi.py; `src/html/` is down
    to `wsgi.py` and `static/`.
13. The split of `utils.py` and `fmt.py` into the utility modules
    (section 6.4); `stellarObjects` is gone.
14. `generate.py` into `planetgen.generation.run_*` and
    `planetgen.cli.generate`; the root file is gone, and `pytest.ini`'s
    `pythonpath` is gone. Done in two PRs: the split first, then the
    editable install with the `sys.path` pushes and `pythonpath`. What
    every command shares (the progress bar, the `--strict`
    refuse-or-warn rule, the run's counts, the work queue, the
    statistics and estimate) went to `planetgen.generation.run_common`,
    so the run modules import downward only (system, sector, galaxy, then
    plan); the bright-star band drawing moved from the galaxy section to
    `run_plan`, which owns the scatter. Command options stay in
    `planetgen.cli.generate`; `run_galaxy` imports the shared generation
    options from it inside `_default_generation_args`, the one place a
    run module needs them.
    The second PR moved the last stray package, `src/wikiClient`, to
    `planetgen.wiki` (no step had named it), limited `setup.py` to the
    `planetgen` packages, and replaced the `/usr/local/bin/planetgen`
    wrapper with pip's console script from the editable install (linked
    there from the venv on macOS). install and update make the editable
    install on every path (`--break-system-packages` on a PEP 668
    Python, as the libraries already were) and skip it when `planetgen`
    already imports from the checkout. The job runner no longer sets
    `PYTHONPATH`, and every command runs from the checkout's root.

### 6.6 Deployment

The server needs `update.sh` after each step that renames a script
install or update runs (steps 8, 9 and 14; step 14 also brings the
editable install). The Apache config doesn't change. The
web server's user needs read access to `src/planetgen/`;
`set-permissions.sh` covers it from step 1.

## 7. Why it works this way

- **Libraries over hand-rolled code**: Boss's directive. The research
  estimates over 10,000 lines removed, and the libraries are audited and
  tested where the hand-rolled code is not (TOTP, QR, SQL building).
- **Move first, swap second**: a swap inside a file that is about to move
  makes every open branch conflict twice.
- **One swap per PR**: each swap is checked against the old behavior on
  its own, so a regression points at one library.
