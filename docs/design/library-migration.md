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
| Work queue and web jobs | `workQueue.py` (1,600 lines), `jobRunner.py`, the `work_lease` table | RQ on Redis | See section 3. |
| Caches | `pagecache.py`, `tilecache.py` | cachetools, diskcache | Same keys and invalidation. |
| Database | `_db.py` (9,800 lines), `migrateDb.py` with `schema_vNN.sql.gz` fixtures | SQLAlchemy, Alembic | Alembic starts from a baseline that recognises existing databases at the current schema. |
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
points at Memurai or a Redis in WSL (default until Boss decides otherwise).

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

## 6. Package layout

The layout plan (to be approved before the move) groups the code roughly
as follows; names may change in the plan:

| Package | Holds |
|---|---|
| `planetgen.physics` | constants, units, orbits, the point-in-space object, planet physics |
| `planetgen.generation` | systems, stars, planets, phenomena, the galaxy model, backfill |
| `planetgen.galaxy` | geometry, skeleton, sectors, navigation |
| `planetgen.names` | the codec-based names from IDs (object-ids.md) |
| `planetgen.db` | SQLAlchemy models, Alembic migrations, queries |
| `planetgen.queue` | RQ jobs and progress publishing |
| `planetgen.population` | species, civilizations, polities, tech levels |
| `planetgen.web` | the Flask app, API and pages |
| `planetgen.util` | shared helpers: number formatting, seeded random helpers, logging setup |

Old import paths keep working through thin re-export modules until every
caller has moved.

## 7. Why it works this way

- **Libraries over hand-rolled code**: Boss's directive. The research
  estimates over 10,000 lines removed, and the libraries are audited and
  tested where the hand-rolled code is not (TOTP, QR, SQL building).
- **Move first, swap second**: a swap inside a file that is about to move
  makes every open branch conflict twice.
- **One swap per PR**: each swap is checked against the old behavior on
  its own, so a regression points at one library.
