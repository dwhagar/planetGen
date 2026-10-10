# Major design decisions

Short records of the choices that shape planetGen as a whole. Each says what
was chosen, when, why, and what was turned down. Versions are from
`CHANGELOG.md`; PR numbers are from the merge commits. Where a reason is not
written down anywhere in the repository, this file says so instead of
guessing.

The git history in this repository starts on 2026-09-23 (PR #66). Earlier
choices are known only from `CHANGELOG.md` and old copies of `docs/TODO.md`.

---

## 1. The cylindrical sector grid (rings, layers, slots; 4 pc sectors)

**Chosen:** every galaxy-placed sector is one cell of a cylindrical grid
around the galactic axis: ring `i` (radius), layer `j` (height, layer 0 on
the plane) and slot `k` (angle), each about one edge long. The edge is 4 pc
(13.05 ly). Rings hold a multiple of a master wedge count (3 at the center,
doubling outward), so slot boundaries line up from the center out. Details:
`galaxy-coordinate-system.md`.

**When:**
- 6.0.0 (2026-09-30, schema v32, PR #98): cylinders replace spherical
  shells; 11.5 ly cells.
- 7.0.0 (2026-09-30, schema v33): 4 pc edge, the same slot count on every
  layer (aligned columns), per-layer skeleton, bounds checked before
  generating.
- 7.13.0 (2026-09-30, schema v35, MAP.33): hybrid master-wedge slot counts.

**Why:**
- Sectors should "follow the flat disk instead of a ball" (CHANGELOG 6.0.0).
  The shell design spent most of its addresses far off the disk.
- The cells tile space exactly and every lookup (point to cell, neighbors,
  cells near a point) is closed form, so the Voronoi prism table
  (`sector_vertices`) and the per-shell band table (`galaxy_shell_band`)
  could be dropped (CHANGELOG 6.0.0).
- 4 pc is a whole number of parsecs and puts the galaxy's edge at about
  50,000 ly at the default shape, the real disk's radius (CHANGELOG 7.0.0).
- Master wedges make every master line a slot boundary in every ring
  outward, so the Galaxy Map's blocks and wedge lines cut on shared lines.
  Boss chose this hybrid on 2026-09-30.

**Rejected:**
- Spherical shells with Fibonacci-sphere slots (the 5.3.6 "Track C"
  design): cubes cannot tile a sphere, so it needed per-sector Voronoi
  prisms and a Fibonacci-number neighbor search.
- A stored plan of every qualifying sector (`galaxy_sector_plan`, about
  10.5 billion rows and 1 TB): never built, because position and density
  are pure functions of the address.
- 2 pc (a local sector would average under one star), 3 pc and 5 pc edges.
- A staggered brick pattern of layers: rejected for aligned columns; no
  reason beyond alignment is recorded.

**Cost:** each of v32, v33 and v35 deleted the affected galaxy-placed
sectors with their systems and phenomena.

---

## 2. Rendering system pages from the database instead of storing text

**Chosen:** a system's wikitext and Markdown pages are rendered on demand
from its rows (`planetgen/db/render.py`). The stored copies
(`star_systems.wikitext_content`, `markdown_content`) are gone.

**When:** 5.52.0 (2026-09-24, PR #85, schema v29; the commit message says
v28, the number before a merge renumbered it).

**Why:** pages now "follow renames, names made unique after generation, and
orbit ticks" (CHANGELOG 5.52.0). Stored text went stale whenever a row
changed after generation. Wiki upload uses the fresh render.

**Rejected:** keeping the stored text and regenerating it on every change.
Not discussed in the record beyond the reasons above. `planetgen.cli.render_parity`
was written to compare stored and rendered text on an old database before
upgrading, and to export the stored copies.

**Cost:** the migration deletes the stored copies; the release note says to
back up first. The rendered page must match generation exactly, which is
why loading a binary now re-derives which star is heavier (same release).

---

## 3. three.js for the Galaxy Map (and the other maps)

**Chosen:** the Sector, System and Galaxy Maps draw with three.js (r186),
vendored as one minified file in `src/html/static/vendor/`, not loaded from
a CDN.

**When:**
- 5.41.0 (2026-09-19): the Sector Map moves from a CSS 3D illusion to a
  three.js WebGL scene.
- 5.46.13 (2026-09-23): the 3D Galaxy Map, with live viewport queries.
- 2026-10-01 (PR #153, MAP.40): the choice re-checked for the Galaxy Map
  and recorded in `docs/html-interface.md`, "Why the maps use three.js".

**Why:**
- A real perspective camera gives correct occlusion, and sprites shrink with
  distance for free, which fixed the flat map's "star icon doesn't shrink as
  you zoom in" problem (CHANGELOG 5.46.13).
- Vendored, so the pages' `Content-Security-Policy: default-src 'self'`
  needs no exception (CHANGELOG 5.41.0).
- The project has no build step, and one vendored file suits that.
- The Galaxy Map's slow part was the JavaScript that lists blocks, not the
  drawing; that work moved to a Web Worker (`galaxyblocks.js`, 7.32.0).

**Rejected (2026-09-30 review):** Babylon.js (several megabytes); deck.gl
(needs a bundler); regl or raw WebGPU (picking, sprites and lighting by
hand). If WebGL runs out, three's own `InstancedMesh`, `BatchedMesh` and
`WebGPURenderer` come first.

---

## 4. A Flask-only website (CGI pages removed)

**Chosen:** every page is served by the same Flask app as the API
(`src/planetgen/web/`, Jinja2 templates, mounted at `/` through mod_wsgi or
gunicorn). Old `/<name>.py` URLs answer with a 301 to their
replacement (`web/old_urls.py`).

**When:**
- 5.2.1 (2026-09-06): the first web interface, plain Python CGI scripts
  with the standard library only.
- 5.54.0 (2026-09-24, PR #92): the Flask page foundation; the home page
  moves.
- PRs #93, #96, #97 and #100 (2026-09-24 to 2026-09-30): the sector, NAV,
  Galaxy Map, admin, system and phenomenon pages move.
- 7.1.0 (2026-09-30, PR #104): the CGI pages, shims and Apache CGI rules are
  removed.

**Why (from the 5.54.0 and 7.1.0 notes and commit d66fbed):**
- Pages are served "without a process start or an HTTP call back to the
  API": inside a request, page code calls the API routes in-process.
- Plain, bookmarkable GET URLs with no database name in them; the database
  comes from `config.json`.
- Autoescaped Jinja2 templates, an app-wide CSRF check, and HTML error pages
  that never show a traceback.
- One shared page shell (header, skip link, breadcrumbs, theme), checked in
  CI by the browser accessibility test (5.59.0).

**Rejected:** keeping the CGI scripts. They were the "interim CGI browser"
(5.10.0) from the days before the API existed. Why Flask was picked over
another framework is not recorded beyond the JSON API already being
Flask (it existed by 5.3.2).

---

## 5. MySQL as the store

**Chosen:** MySQL 8.0.16+ or a compatible MariaDB, InnoDB tables, through
`pymysql` with a `DBUtils.PooledDB` connection pool. Schema changes are
numbered migrations recorded in `schema_migrations` (now at v43).

**When:** 5.5.0 (2026-09-09), "TODO.md Phase 5", replacing SQLite entirely.
This predates the repository's git history, so there is no PR.

**Why:** the old `docs/TODO.md` Phase 5 entry says the move was "for real
concurrent multi-user access", together with "real connection pooling" for
the web interface. The `schema.sql` port notes add that `BIGINT` keys give
headroom at galaxy scale (tens of millions of bodies).

**Rejected:** SQLite (one writer at a time, which the shell-era plan table
also hit; `docs/design/archive/galaxy-disk-density-rev2.md`, section 5).
Other servers such as PostgreSQL are not discussed anywhere in the
repository.

**Cost:** SQLite's in-place v1 to v5 migrations were removed. The one-time
SQLite import script (`src/migrateSqliteToMysql.py`) was retired too
(TEST.61): it only took a file at the current schema version, which no
SQLite database ever reached (SQLite stopped at v5, and MySQL migrations
start at v8), so it could never run.

---

## 6. Release notes in `changes/`, versions stamped after merge

**Chosen:** a PR does not edit the version. It adds one note,
`changes/<name>.<patch|minor|major>.md`. After the merge,
`.github/workflows/stamp-version.yml` runs `scripts/bump_version.py`, which
gives each note its own version, writes `_version.py`, the README badge and
`CHANGELOG.md` together, and commits `Release x.y.z`. A PR check rejects a
hand-edited version or a missing note (unless labelled `no-release`).

**When:** 5.46.33 (2026-09-24, PR #68, commit 7ed99d6).

**Why:** parallel PRs kept claiming the same next version number and needed
hand-renumbering merges (5.41.0/5.41.1 and 5.46.15 are named in the release
note). A uniquely named note per PR cannot collide.

**Rejected:** bumping the version inside each PR (the old way).

**Later (2026-10-01, OPS.1):** the version became MAJOR.REVISION.BUILD.
A `major` note bumps MAJOR, any other note bumps REVISION, and BUILD was the
sum of the TODO category counters in `docs/design/todo-number-map.md`
(Boss: "major feature set.revision.build"). **Changed 2026-10-10:** Boss
wanted BUILD to change on every stamp, so BUILD is now the previous BUILD
plus one and never resets; the counters only allocate TODO IDs. See
`changes/README.md`.

---

## 7. System Python with apt on Linux, venvs elsewhere

**Chosen:** on Linux, `install.sh` and `update.sh` put the libraries into the
system Python. On an externally managed Python (PEP 668: Debian 12+,
Ubuntu 23.04+), each library comes from apt when apt has it at or above
`setup.py`'s floor; only what apt lacks or ships too old is pip-installed
system-wide into `/usr/local`, alongside apt's copy, from the hash-checked
`requirements.lock`. On macOS the installer builds a venv
(with gunicorn).

**When:**
- 5.58.0 (2026-09-30, PR #99): first support for externally managed Python,
  with apt packages plus a fallback venv joined by a `.pth` file.
- 7.2.0 (2026-09-30, PR #105, commit b3d533e): "No more venv: libraries go
  into the system Python."
- Commit 38a67c5 then moved everything into a dedicated venv at
  `/opt/planetgen/venv`, and commit 498df96, in the same PR, reverted it.
- 7.7.0: pip installs only locked, hash-checked files.
- 7.16.0: macOS installer, with a venv.

**Why:**
- Boss "asked in words for the system Python with apt first" (commit
  498df96).
- pip must never overwrite apt's files: a plain
  `pip install --upgrade --break-system-packages` deletes apt's copy of
  some packages (python3-pymysql) and breaks dpkg, so pip resolves first and
  installs alongside with `--ignore-installed --no-deps` (commit b3d533e,
  `docs/deployment/apache.md`).
- mod_wsgi embeds the system Python, so libraries in that Python are what
  Apache actually imports (7.2.0 adds a check for a mismatch).

**Rejected, for now:** the dedicated venv at `/opt/planetgen/venv`, reached
by mod_wsgi through `python-home` (commit 38a67c5, "the standard deployment
pattern Boss shared"). It was reverted because "the venv is still an open
question" to Boss, and the commit is kept in history to restore if Boss
chooses it (commit 498df96). On macOS there is no system
package manager to lean on, so a venv is used there; no further reason is
recorded.

---

## 8. Population as a separate pass, optional and off by default

**Chosen:** species, civilizations, polities and territories are built by
their own pass over what is already stored (`generate.py population`,
`planetgen/population/model.py`), not during system generation. The pass
is opt-in: `generate.py sector` and `galaxy` run it only with
`--population`, the admin Generate page never passes that flag, and
`install.sh` asks y/N with a
30-second timeout that defaults to No and skip the question with no
terminal; `POPULATION=1` runs it without asking. The pages
and the Galaxy Map's Territories button hide themselves when there is no
population data (`GET /api/population`, a polity count). Details:
`population-and-politics.md`.

**When:**
- 7.49.0 (2026-10-01, PR #169, schema v44, POP.1 to POP.4): the pass ships
  and runs after every `sector` and `galaxy` fill unless
  `--no-population` is given (Boss's decision 5, "territories recompute
  automatically after each fill").
- 7.58.1 (PR #181): the Territories button appears only once a polity
  exists.
- 7.58.2 (PR #180): optional and off by default everywhere; `--population`
  replaces `--no-population`; `/api/population` added.

**Why:**
- A separate pass reads only stored rows (each planet's evolution text and
  each system's position), so a galaxy generated before v44 can be
  populated without a regenerate, and nothing in system generation had
  to change (design doc, "When it runs").
- A watermark on `planets.id` (`population_state`) makes a rerun cheap: it
  scans only planets added since, then recomputes eras and territories,
  which are derived and can be rebuilt from scratch.
- Draws are seeded by the planet id, so traits, civilization odds and
  ages come out the same on every run.
- Off by default: the design doc records it only as Boss's choice of
  2026-10-01, superseding decision 5 the same day; no further reason is
  written down. In practice it means a scheduled `update.sh` never runs
  the pass, and a fill's time is spent on sectors only.

**Rejected:**
- Generating species and polities inside `StarSystem`: it would have
  needed a regenerate of every existing galaxy, and territories depend on
  neighboring systems a single system cannot see.
- Recomputing territories automatically after every fill (decision 5 as
  first accepted): replaced by the opt-in above.
- Taking the timeline's `technological_civilization` milestone as a living
  civilization: a 150-sector test fill put it on one system in seven,
  which made nearly every one a polity (commit b62d0aa). A 1 in 1,000
  chance (`CIVILIZATION_CHANCE`) decides instead.

---

## 9. An in-memory cache for the pages' API answers

**Chosen:** each web process keeps the API's public GET answers in memory
(`planetgen/web/lib/pagecache.py`, asked by `apiclient._request` before each call).
Any successful API write in the process clears it; at most once every 15
seconds per database it reads the galaxy content stamp (`GET
/api/galaxy/changes`) and drops that database's entries when the stamp
moved; nothing is kept past 5 minutes; it is capped at 2,000 entries and
64 MB. Every failure fails open. `page_cache` in `config.json` tunes it,
and `PLANETGEN_PAGE_CACHE=off` turns it off.

**When:** 7.56.0 (2026-10-01, PR #178, PERF.2).

**Why:**
- Every page called the API and every API route queried MySQL fresh, even
  for answers that rarely change (PERF.2). The Galaxy Map's tiles
  already had a cache (`tilecache.py`), and its content stamp could be
  reused, so the page cache needs no new database columns.
- Only cookie-less GETs are cached, so every visitor gets the same answer
  and nothing that depends on who is asking (`/auth`, `/admin`) is
  shared. Tile and stage data stay in the disk cache.
- Clearing on any write in the process makes an edit show on the very
  next page; the stamp catches writes from other processes (a
  `generate.py` job, another worker); the age limit catches what the
  stamp cannot see, such as a rename from the command line.

**Rejected:** the original plan for PERF.2, invalidating each page from its own
rows' `modified_at`, in the style of the tile cache. The commit does not
say why; the chosen scheme is simpler and needs only the one stamp call.

---

## 10. A drill-down for the Galaxy Map instead of a free camera

**Chosen:** the Galaxy Map opens on discrete stages: blocks 243, 27 and 3
sectors a side, then single sectors, alternating a 3D view of a block's
children and a top-down view of one slab, eight clicks from the galaxy to
a sector. Each stage has its own URL. The old free camera stays behind a
Free look button for now. Details: `galaxy-drilldown-navigation.md`.

**When:**
- 7.41.2 and 7.41.3 (PR #160): the block ladder and the stage API.
- 7.44.0 (2026-10-01, PR #171, MAP.16): the stages on the map.
- 7.50.0, 7.52.0, 7.58.0: the address bar, the NAV course on the map, and
  the Sector Map's pick mode, all built on the stage URLs.

**Why:** Boss asked to "select in sections" (design doc, section 1), and
the research he shared supports discrete scale tiers, a breadcrumb and
animated flights. Each stage draws at most about 900 blocks, so the
phone-performance limits of the free camera stop mattering. Boss chose
243, 27, 3, 1 ("the bigger targets") over 81, 9, 1.

**Still open:** MAP.17 (2026-10-01) asks to remove the free camera and
every drag-rotate and to start top-down, picking a wedge, then a slice,
then a block. That would settle the design doc's decision 2 against Free
look and change how the first stages work.

## 11. Third-party libraries instead of zero dependencies (planned)

**Chosen:** Boss (2026-10-03): "we have an explicit directive to move the
code base from a 0-dependency model into using 3rd party open source
libraries to simplify our own code deployment." Hand-written TOTP, QR,
Kepler solvers, SQL building, the process-pool queue, rate limiting,
the page cache, Markdown and validation give way to pyotp, segno,
scipy, astropy, SQLAlchemy with Alembic, RQ on Redis, Flask-Limiter,
cachetools, markdown and Pydantic (the JSON tile cache stays); the pages gain Shoelace, TanStack,
Xterm.js and three-mesh-bvh. The code is first reorganized into
importable packages. Details: `library-migration.md`.

**When:** planned 2026-10-03 and 2026-10-07 as phase 0 groundwork; not
built yet.

**Why:** less code to maintain, audited implementations of security code,
and several open bugs (job logs, test flakes, menus, slow map zoom) are
fixed by the libraries rather than patched.

**Turned down:** the job guide's broker-less process pool (Boss chose
Redis); Huey and Celery (section 3 of `library-migration.md`).

## 12. Names from IDs (planned)

**Chosen:** Boss (2026-10-03): "Do away with name generation for any
object, do away with our word-salad code entirely", narrowed on
2026-10-08: "keep word salad method for stars and sectors, everything else
gets a name derived from it's unique ID." Every object gets a unique ID.
Stars and sectors keep their word-salad names and, by Boss's decision
of 2026-10-08 04:00Z, planets, moons and belts keep the "<star> I"
pattern; the objects with no star-derived name (rogue planets, black
holes, neutron stars, nebulae, remnants, quasars, interstellar comets,
asteroid fields, constellations) get Boss's phoneme codec
(`planetgen.names.gated_phoneme_codec`) applied to the ID under a naming key stored in
the control database, drawn when the galaxy is created and changeable by
an admin. Details: `object-ids.md`.

**Why:** names become unique by construction (no collision rules, no
registries, no word lists to keep reproducible), and changing the key
renames the whole galaxy without touching a row.

## 13. Every sector can be generated (planned)

**Chosen:** Boss (2026-10-03 and 2026-10-07): "No sector should generate
and not fill with stuff even if there are no star systems" and "The
projected density is just to be used in generation, not to block
generation anywhere." Every address inside the galaxy's bounds can be
generated; the density only sets the odds (with a small floor
everywhere), and `sector_stats.bright_level_sol` = 0 is the one
"generated" flag. The console warns but does what it is asked.

**Why:** the "doesn't qualify" rule left sparse and edge sectors
impossible to fill, broke neighborhood runs, and kept rogue planets,
comets and fields out of empty sectors.
