# Changelog

## [5.46.29] - 2026-09-24

### Changed
- **Galaxy Map (3D): a placed (already-generated) sector's own marker is
  now colored by its real stellar density** (`system_count / edge_ly **
  3`, relative to `physical_constants.LOCAL_STELLAR_DENSITY_LY3` -- the
  real local-neighborhood average this whole generator already
  calibrates against), not just sized by raw system count. Marker size
  still scales with `system_count` as before; only the color (dim bronze
  at low density, bright gold at high) is new. The info panel gained a
  "Density" field showing the same ratio (e.g. "1.8x local average").
  `queryDb.galaxy_sectors_in_view`'s own returned shape gained `edge_ly`
  per sector to make this possible.

## [5.46.28] - 2026-09-24

### Fixed
- **Galaxy Map (3D): rapidly clicking/double-clicking or mashing the
  zoom buttons could fire off a live API request per click, faster than
  the server can service them.** `doFetch`'s own `activeAbort.abort()`
  only stops the *browser* from waiting on a superseded response -- it
  doesn't reliably stop the server from finishing a query it already
  started (Flask/WSGI doesn't check for a disconnected client mid-query
  unless specifically coded to), so rapid clicking still burned a real
  WSGI thread/DB-connection-pool slot per click even when every earlier
  response got thrown away client-side the instant the next one fired --
  a real contributor to the production connection-exhaustion pattern
  already fixed elsewhere in this and recent releases. Every interaction
  that requests an immediate fetch (click, double-click, the +/-/reset
  buttons) now funnels through one shared cap: at most 4 accepted
  immediate fetches per second. A click faster than that still moves the
  camera/selection instantly (never throttled), it just falls back to
  the existing debounced delay for its own data fetch instead of firing
  right away, so a rapid burst still settles on exactly one fetch
  shortly after it stops rather than either hammering the server once
  per click or never syncing to the final camera position at all.

## [5.46.27] - 2026-09-24

### Added
- **`browse.py`'s Sectors and Standalone Systems tables are now really
  paginated** (100 rows/page, independent Prev/Next controls per table)
  instead of a single page capped at 500 rows with a "try Search
  instead" hint and no way to ever reach anything past that cap. Each
  table's own `sector_offset`/`standalone_offset` page position is
  independent, so paginating one never resets the other back to page 1.

## [5.46.26] - 2026-09-24

### Fixed
- **Rogue planets and interstellar comets were completely invisible
  everywhere** -- generated at a non-trivial rate
  (`program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM`'s own
  `"rogue-planet": 0.1` and `"comet": 0.05`, roughly one rogue planet per
  ten star systems, far more common than a nebula) and saved to the
  `rogue_planets`/`interstellar_comets` tables the whole time, but no
  query function anywhere (`list_phenomena`, `count_phenomena`,
  `phenomenon_detail`) ever read either table, so they never appeared in
  the Phenomena listing or had a detail page of their own, despite
  existing in the database. Both are now wired up the same way
  `supernova_remnant` already was (no galaxy-frame placement columns of
  their own, so still absent from the Galaxy Map / Sector Map / NAV --
  only the flat listing and detail page gain them). `html/phenomenon.py`
  gained field specs for each type's own real columns (a rogue planet's
  `planet_type`/`mass_kg`/`composition`/etc., a comet's
  `nucleus_diameter_km`/`velocity_kms`/`is_active`/etc.).

  Nebulae, by contrast, were already fully wired up -- their apparent
  rarity is real, deliberately-researched astronomical calibration
  (`"nebula"` rate is a cited-literature `5e4 / 2e11` per star system,
  several orders of magnitude below a rogue planet's own rate), not a
  bug; expect one to actually appear only in a very large generated
  galaxy.

## [5.46.25] - 2026-09-24

### Fixed
- **`migrate_database` never actually applied 5.46.19's v26 migration**
  (the spatial indexes on `nebulae`/`asteroid_fields`/`black_holes`/
  `neutron_stars` that fix `sector.py`'s production timeout) -- the
  `_migrate_v25_to_v26` function existed but was never called from
  `migrate_database`'s own version cascade, so running `migrateDb.py`
  against an existing (pre-v26) database left it silently stuck at v25
  forever, `GET /api/health`'s `schema_current` reporting `false`
  indefinitely with no way to clear it short of applying the index by
  hand. Caught by running the full test suite against a real MySQL
  server rather than relying on syntax/logic review alone -- every
  `test_migrate_v*` regression test failed with `schema_version == 25`,
  not `SCHEMA_VERSION` (26). A fresh database (`_ensure_schema`, which
  reads the indexes straight from `schema.sql`'s own `CREATE TABLE`) was
  never affected -- only a database migrated from an earlier version.

## [5.46.24] - 2026-09-24

### Fixed
- **System Map: an orbit line drawn under a body's own live sphere
  visibly cut across it** instead of being occluded -- the sphere is
  drawn on a separate `<canvas>` layered by CSS z-index relative to the
  SVG, so within a single `<svg>` a body's own opaque sphere had no way
  to occlude a sibling orbit-line element painted in that same stacking
  context. Each scene is now built as two sibling `<svg>`s (an
  aria-hidden orbits-only layer, and the existing body-marker layer,
  both toggled together by `static/systemmap.js`'s `showScene`) stacked
  either side of the sphere canvas, so a sphere now actually covers the
  orbit line drawn under it. Belt rings (interactive markers, not
  decorative lines) stay in the body-marker layer.

### Changed
- **System Map: stars now render with a mottled granulation texture**
  (layered sine "turbulence" in both UV directions, tinted to the star's
  own spectral color) instead of a flat single-color sphere.
- **System Map: a star's own glow shell is now bigger and brighter than
  a planet's subtle atmosphere rim** (wider falloff, higher intensity,
  larger radius), reading clearly as a light source rather than the same
  faint haze a planet's atmosphere gets.

## [5.46.23] - 2026-09-24

### Changed
- **Sector Map's compass arrow (pointing toward the galactic center) now
  labels itself plain "N"**, matching a real map's compass-rose
  convention, instead of the more verbose "Galactic Center →" text.

## [5.46.22] - 2026-09-24

### Fixed
- **Galaxy Map (3D): a selected placed/planned sector could disappear
  entirely once zoomed in close to it.** The live re-fetch's own bounding
  box shrinks as the camera's orbit radius shrinks; a click/double-click
  that landed even slightly off a sector's own exact stored position
  (easy from far out, where its marker is only a handful of screen
  pixels) meant a later, smaller-radius re-fetch could legitimately no
  longer include it, and the client dropped anything missing from a
  fresh fetch. The selected entry is now pinned client-side and
  re-inserted into each fetch's own tier if the live query didn't happen
  to return it, so it stays in the scene for as long as it's selected.

## [5.46.21] - 2026-09-23

### Changed
- **Galaxy Map (3D) interaction model reworked:** left-click now only
  centers the view on the clicked dot/empty space and selects it
  (previously it also zoomed in, which punished an imprecise click by
  zooming into empty space nowhere near the intended target -- the
  likely cause of "I can't zoom into known space" once a dot was too
  small/far to click precisely from the full-galaxy starting view).
  Double-click now does what a single click used to (center, select,
  AND zoom in by one `clickZoomFactor` step). Right-click no longer does
  anything (previously zoomed out) -- the browser's own default context
  menu is left alone instead of being suppressed for nothing. The panel's
  own hint text/`aria-label` (`lib/galaxymap3d.py`) updated to match.
- **The illustrative density cloud (the "shows the spiral arms" layer)
  now renders as soft, translucent, additively-blended spheres
  (`THREE.InstancedMesh`) instead of tiny flat `THREE.Points` dots** --
  it read as a sparse scatter-plot rather than shaded spiral structure.
  Each sphere's size varies with its own `relative_density` (denser
  regions read as visibly bigger/brighter blobs) and with the camera's
  current orbit radius (so the cloud keeps a sensible relative size
  across zoom levels); additive blending lets overlapping spheres
  brighten rather than simply occlude, the cheap way many soft blobs
  merge into continuous-looking shading along a spiral arm.

## [5.46.20] - 2026-09-23

### Fixed
- **A `:80`-to-HTTPS redirect vhost with no exclusion for `/api` silently
  turns every CGI page's own internal API call into a real public
  round trip back into the same server, doubling Apache's connection
  load per page view.** Confirmed against a real production vhost:
  `html/lib/apiclient.py` defaults to `PLANETGEN_API_BASE_URL=http://
  127.0.0.1/api`, a plain-HTTP loopback call every page makes at least
  twice (its own data, plus `/api/auth/me`); a blanket `Redirect
  permanent /` (or an equivalent `RewriteRule`, including certbot's own
  default `--apache` rewrite) on the `:80` vhost catches that loopback
  call too, and `urllib` follows the redirect out through DNS/TLS/
  anything in front of the box and back in -- exactly the connection-
  exhaustion/timeout pattern (with refused connections that never reach
  the error log, since Apache never accepted them) reported alongside
  the `sector_detail`/`galaxy/view` timeouts this and recent releases
  already fixed. `examples/apache/planetgen.conf.example`'s own "HTTPS"
  section previously instructed copying the whole `:80` block (daemon
  process declaration included, which would conflict once duplicated)
  into `:443` and redirecting `:80` unconditionally -- it now excludes
  `/api` from that redirect and mounts `/api` locally on `:80` too,
  reusing one server-wide `WSGIDaemonProcess` declaration instead of two
  conflicting ones. `docs/api.md`'s "Deploying behind Apache" section
  cross-references the same warning.

## [5.46.19] - 2026-09-23

### Fixed
- **`GET /api/sectors/<id>` (`html/sector.py`'s page, and its Sector Map)
  and, under load, unrelated pages sharing the same single-process
  `planetgen-api` WSGIDaemonProcess (`html/system.py` included) were
  timing out / failing outright in production** ("Truncated or oversized
  response headers received from daemon process" and read-timeout errors
  in `planetgen.error.log`, clearing only after an Apache restart -- the
  same failure signature [5.46.16]'s `idx_sectors_center` fix addressed
  for the Galaxy Map). Root cause: `queryDb.phenomena_near_sector` called
  `_placed_phenomenon_rows` with no bounding box at all, so *every*
  `sector_detail` call did a genuine full-table scan across all four
  placed-phenomenon tables (`nebulae`/`asteroid_fields`/`black_holes`/
  `neutron_stars`), pulling every galaxy-placed phenomenon in the entire
  database into Python on every single sector page view. Once a database
  had a non-trivial number of placed phenomena, a handful of concurrent
  sector-page views were enough to hold every one of the API's 5 worker
  threads (and, in turn, its MySQL connection pool, which has no
  checkout timeout) in slow queries at once, starving every other
  request behind them until Apache was restarted.
  `_placed_phenomenon_rows` now takes an optional SQL bounding-box filter
  (the same `BETWEEN`-range-scan technique [5.46.16]'s fix used for
  `sectors`), and `phenomena_near_sector` uses it, padded by the widest
  currently-placed phenomenon radius (a cheap `MAX(radius_ly)` query, not
  a fixed assumption, so it stays exactly as correct for an arbitrarily
  large placed phenomenon as the old unconditional scan). New schema v26
  adds the matching spatial indexes
  (`idx_{nebulae,asteroid_fields,black_holes,neutron_stars}_center`) --
  **run `migrateDb.py` (or `update.sh`/`install.sh`) against any existing
  deployment's database for this fix to actually take effect**; `GET
  /api/health` (see [5.46.17]) will report `schema_current: false` in the
  meantime.

## [5.46.18] - 2026-09-23

### Fixed
- **`GET /api/health` was returning a bare 503 "error" for a database
  that's reachable but has never had `schema.sql`/`migrateDb.py` applied
  to it at all** (no `schema_migrations` table yet -- one step further
  back than "some migrations pending", which it already handled). Caught
  by CI: `test_health_ok` exercised exactly this case by accident (its
  `client` fixture's database starts completely empty) and failed after
  [5.46.17]'s health-reporting change merged. `/api/health` now reports
  this the same way it already reports a stale-but-present schema
  (`200`, `schema_current: false`, a `detail` naming the fix) rather than
  folding it into the "unreachable" 503 case, which is meant for an
  actually-unreachable server. `test_health_ok` now lays the schema down
  first (matching how a real deployment's database always already has
  one by the time its API is queried), and a new test covers the
  never-migrated-at-all case directly.

## [5.46.17] - 2026-09-23

### Fixed
- **[5.46.16]'s own changelog entry claimed a plain app restart applies a
  pending schema migration (e.g. that fix's own `idx_sectors_center`)
  "automatically" -- it doesn't.** `GET /api/health` now reports
  `schema_version`/`schema_current` (and a `detail` message naming the
  fix) by comparing the database's own `schema_migrations` table against
  the code's `SCHEMA_VERSION`, so a deployment that pulled in a
  schema-fixing code change but never actually ran `migrateDb.py` (or
  `update.sh`/`install.sh`) against its database is visible at `/api/health`
  instead of continuing to silently run the old, unmigrated schema -- the
  likely explanation if the same full-table-scan timeouts (and the site
  going down under load) kept happening after [5.46.16]'s code shipped
  but its database was never separately migrated.

## [5.46.16] - 2026-09-23

### Fixed
- **`GET /api/galaxy/view` (the interactive 3D Galaxy Map's live-viewport
  query, added in [5.46.13]) was a genuine full-table scan on `sectors`
  every single call -- confirmed in production as the site going down
  under load (`TimeoutError`/"Truncated or oversized response headers"
  from the WSGI daemon), the exact same failure mode already documented
  and fixed once before for the pre-v22 `/api/search` ([5.35.7], schema
  v22). `sectors.center_x/y/z_pc` had no index
  (`queryDb.galaxy_sectors_in_view`'s own bounding-box `WHERE` clause said
  so explicitly), and this endpoint is hit hard: once server-side on
  every galaxy-map page load (the zoomed-all-the-way-out starting view,
  spanning the whole galaxy) and repeatedly (debounced) as the 3D camera
  moves. New schema v25 adds a composite `idx_sectors_center` index
  (`_db._migrate_v24_to_v25`) so the query can range-scan instead of
  examining every row. **Run `migrateDb.py` (or `update.sh`/`install.sh`,
  which call it) against the database to pick this up** -- restarting the
  app alone does *not* apply it: the API's own connections are read-only
  (`ensure_schema=False`) and never run schema DDL at all, and even a
  read-write connection's `_ensure_schema` only runs `CREATE TABLE IF NOT
  EXISTS`, a no-op against a table that already exists (see [5.46.17]).

## [5.46.15] - 2026-09-23

### Changed
- **`galaxy.py` (the "Galaxy Map" page) now renders the 3D map directly**,
  in place of the flat SVG projection -- rather than living alongside it
  as a separate `galaxy3d.py` page (5.46.13's original approach). The
  flat SVG rendering code (`lib/galaxymap.py`'s former
  `render_galaxy_map_panel` and `static/galaxymap.js`) is removed
  entirely; `lib/galaxymap.py` keeps only the plain Quadrant/Ring
  classification math `galaxy.py`'s own data tables and
  `sector.py`/`browse.py`'s "Quadrant N" links still need.

### Fixed
- **The 3D map's illustrative density cloud was invisible from the
  starting full-galaxy view.** Its points were sized in world-space
  parsecs with perspective attenuation on -- correct for something meant
  to represent real physical size, but a 1-2 pc point shrinks to
  sub-pixel from thousands of parsecs away. Switched to a constant
  on-screen pixel size (`sizeAttenuation: false`) so the cloud stays
  visible at any zoom level.
- **The density cloud didn't read as a recognizable galaxy shape.**
  Points were drawn uniformly across the whole query volume -- at a wide
  view, almost all of that volume is near-empty halo, so only a sparse,
  shapeless scatter of points ever landed somewhere bright. Switched to
  importance sampling from a bulge+disk mixture shaped like the galaxy's
  own real mass distribution, so the cloud now visibly reads as a bright
  core plus a disk (spiral-arm contrast still comes through via each
  point's own real predicted density driving its color).
- **Placed/not-yet-generated sector markers were effectively invisible
  from a wide view** (the same world-space-sizing problem as the density
  cloud above) -- "I don't see the sectors we've generated anywhere."
  Both tiers now use the same constant-screen-size billboard technique
  (recomputed every frame from each marker's own live distance to the
  camera), so a generated sector stays a visible, clickable dot
  regardless of how far the camera currently is, without ever growing to
  dominate the view up close either.

### Added
- **`tests/galaxy_shape_visualizer_cli.py`** -- a diagnostic tool that
  renders the galaxy's real density model as an actual face-on/edge-on
  image (matplotlib), for directly eyeballing whether a set of shape
  parameters produces a recognizable spiral rather than only judging it
  through `relative_density` numbers. Optionally overlays every real,
  already-generated sector's own position when given `--mysql-*`
  connection args.

## [5.46.14] - 2026-09-23

### Changed
- **Every star/planet/moon marker on the System Map now renders its own
  live 3D sphere in place, not just a floating preview beside a click.**
  Previously, only the one planet/moon last clicked got a rotating 3D
  preview -- floated in a small box beside its flat marker (`#sysmap-
  preview`) rather than replacing it, and the whole-system view's star and
  every other unclicked body stayed flat 2D circles regardless. Every
  visible marker (star included) now gets its own sphere, sized and
  positioned to exactly cover -- and read as replacing -- its own flat
  circle, all drawn each frame through one shared WebGL canvas
  (`#sysmap-spheres-canvas`, `lib/systemmap.py` + `static/systemmap.js`)
  via a scissored sub-viewport per marker, rather than one `<canvas>`/
  context per body (which would risk exceeding a browser's cap on
  concurrent WebGL contexts in a crowded system). A star gets its own
  unlit, spectral-color-tinted sphere (`_star_color`, newly exposed as
  `data-color` same as a planet/moon's own class color) plus a matching
  glow shell; falls back to the plain flat marker for a browser that
  can't create a WebGL context at all.

## [5.46.13] - 2026-09-23

### Added
- **Interactive 3D Galaxy Map.** New "Galaxy Map (3D)" page
  (`galaxy3d.py`, linked from the existing flat Galaxy Map) -- a real
  perspective-camera WebGL scene (three.js, `lib/galaxymap3d.py` +
  `static/galaxymap3d.js`) a visitor can rotate, dolly, and click through,
  instead of only ever viewing the galaxy from directly above the disk.
  Because a real 3D camera scales sprite size with distance for free, this
  also fixes the flat map's own "star icon doesn't shrink as you zoom in"
  scaling problem, without any special-case code.
  - **Live viewport queries, not one whole-galaxy payload.** New
    `GET /api/galaxy/view` (`queryDb.galaxy_view`, backed by a new pure
    `stellarObjects.galaxyViewport` module) returns, for whatever the
    camera's current view actually covers: real, already-generated
    sectors nearby; real, not-yet-generated sector addresses this
    galaxy's own density model predicts would qualify (exact, out to a
    200 pc cap); and, for the rest of a wider view, a coarse illustrative
    density point cloud. Fetched (debounced) by the page's own
    client-side JS directly from a new browser-facing proxy,
    `galaxy_view.py`, as the camera moves -- never baked into one page
    load the way the flat map's own dataset is.
  - **Logarithmic click-to-zoom.** Left-click zooms in on whatever's
    under the cursor (a sector, a real not-yet-generated address, or
    empty space), right-click zooms out -- both by a step size that
    shrinks the closer the camera already is (big multiplicative jumps
    while zoomed out over the whole galaxy, fine ones once close to a
    single sector), rather than a flat factor that's either too slow to
    cross the galaxy or too coarse to land on one sector.
  - **Sector designation/address, surfaced and copyable.** Clicking a
    real, not-yet-generated address now shows its provisional designation
    and a "Copy CLI command" button with the exact
    `generate.py galaxy --shell K --slot N` invocation to generate it.
  - **`generate.py galaxy --shell K --slot N`.** New single-address
    generation mode (on top of existing `--shell` batch and
    `--center-sector` neighborhood modes) -- generates exactly the one
    sector slot at that address via the existing lazy-generation entry
    point (`ensure_sector_generated`), the direct path from a designation
    copied out of the new 3D map into this script.

## [5.46.12] - 2026-09-23

### Added
- **Zoomable, real-scale diagram on every stellar phenomenon's page.**
  `phenomenon.py` gains a "Diagram" panel (new `lib/phenomenonmap.py` +
  `static/phenomenonmap.js`, reusing `static/mapzoom.js`'s shared zoom/pan)
  drawn directly to astronomical-unit scale: a nebula/asteroid field's real
  `radius_ly` becomes an actual to-scale circle, zoomable from about 1 AU
  up to 1 ly across. A black hole/neutron star (whose real size is
  negligible at this scale) instead shows a small fixed illustrative dot.
- **Supernova remnants are now a full first-class phenomenon type.**
  Previously `supernova_remnants` had no web page at all. `phenomena.py`'s
  listing and `phenomenon.py`'s detail page (morphology, progenitor type,
  age, radius, any compact remnant left behind, galactic orbit) now cover
  it, plus its own real-scale Diagram panel. It has no galaxy-frame
  placement columns of its own, though (unlike the other four phenomenon
  types), so it never appears on the Galaxy Map and can't be a NAV
  endpoint -- `phenomenon.py` shows a short note explaining this instead
  of offering "Navigate from/to here" buttons that would only fail.

### Fixed
- `queryDb.nav_between` would 500 with a raw "unknown column center_x_pc"
  SQL error if ever asked to resolve a supernova remnant as a NAV
  endpoint (its table genuinely has no such column). Now raises a clean
  `ValueError`, same as any other invalid NAV request.

## [5.46.11] - 2026-09-23

### Added
- **Real interactive zoom/pan on the Galaxy Map.** With enough placed
  sectors, the core cluster used to squash into what looked like a single
  dot no matter how many sectors actually existed -- the map was one fixed,
  non-interactive SVG scaled to fit the single farthest-placed sector, and
  `?quadrant=` only cropped that same squashed drawing. Scroll/wheel to
  zoom (centered on the cursor), drag to pan, and `+`/`-`/`Reset view`
  buttons, all the way in to about 100 ly across -- implemented as `viewBox`
  mutations on the existing server-drawn SVG (new shared
  `static/mapzoom.js`, reused as-is by a future phenomenon-diagram zoom;
  `static/galaxymap.js` wires it to the galaxy map's own scale readout).
  The `?quadrant=` crop is unchanged as the map's *starting* view; zoom/pan
  layers on top of it.

### Fixed
- **A marker/label click anywhere on the Galaxy Map stopped navigating**
  partway through implementing the above: capturing the pointer on
  `pointerdown` (needed so a drag that leaves the SVG mid-gesture keeps
  panning) retargeted the resulting `click` event to the `<svg>` itself
  per the Pointer Events spec, so `navform.js`'s delegated
  `closest("[data-nav-target]")` lookup never found the actual marker.
  Deferred `setPointerCapture` until a real drag is detected instead of
  calling it unconditionally on every pointerdown.

## [5.46.10] - 2026-09-23

### Added
- **"Navigate to here" (symmetric with the existing "Navigate from here"),
  and NAV support for standalone phenomena.** `system.py` now offers both
  directions; `phenomenon.py` gains both buttons too (nebulae, asteroid
  fields, black holes, and neutron stars can now be NAV origins/
  destinations, including a full optimal route via adjacent systems, not
  just a direct course). `nav.py`'s origin picker now accepts an
  already-known destination (from a "Navigate to here" link) and carries
  it through to the course instead of re-prompting. `GET /api/nav` gained
  `from_kind`/`to_kind`/`from_type`/`to_type` query parameters for this
  (see `docs/api.md`'s NAV section) -- existing system-to-system callers
  are unaffected.

### Fixed
- **`GET /api/nav` 500'd for any route involving a phenomenon endpoint.**
  `route.positions` could end up with both an int key (a system hop) and
  a string key (a phenomenon endpoint), and Flask's default JSON
  serialization sorts dict keys, which raises `TypeError` comparing an
  int to a string. Fixed by stringifying every key at the JSON boundary
  (`api/routes.py`), confirmed against a running instance and covered by
  a regression test.

## [5.46.9] - 2026-09-23

### Changed
- **Compacted the spread-out page header on System/Sector/Galaxy/Phenomenon
  pages into one line.** The breadcrumb, Octant/binary badges, "Navigate
  from here" button, and nearest-location text used to each be a
  separately stacked, full-width block -- on `system.py` alone that was 4
  lines of near-empty vertical space before the actual content started.
  Added a shared `.page-subhead` flex row (`static/style.css`) and wired
  it into `system.py`, `sector.py`, `galaxy.py`, and `phenomenon.py`.
  Verified in a browser: the header on a binary system's page shrank from
  roughly 650px of vertical space to about 150px.

## [5.46.8] - 2026-09-23

### Changed
- **System Map: better label collision avoidance, and the 3D body preview
  now lives inside the map itself.** A planet/moon label that collided
  with a neighbor used to try only "above"/"below" before giving up and
  hiding the label entirely; it now also tries "right"/"left", then a
  further-out "above"/"below" tier (connected back to its marker with a
  short leader line) before giving up -- in a stress test, this cut the
  hidden-label rate from ~40% to ~16% for a tightly packed cluster, with
  zero label-to-label overlaps either way. Separately, the rotating 3D
  sphere preview (`#sysmap-preview`) used to sit in a fixed sidebar box
  next to the map; it now floats inside the map viewport itself, next to
  whichever marker was just clicked.

## [5.46.7] - 2026-09-23

### Fixed
- **A planet's evolutionary/civilization narrative could contradict its own
  class description.** A planet's fixed per-class flavor text (e.g. Class
  G: "a rocky, barren world with simple life") and its evolutionary/tech
  milestone (`evolution.get_evolutionary_timeline`, up to "Technological
  Civilization") were generated by two completely independent code paths
  -- the planet's own class was never passed into the evolutionary-timeline
  calculation at all, so an old/fast-evolving star, or a forced
  `INTELLIGENT_LIFE=True`, could still report a full technological
  civilization for a planet whose own class text says life there tops out
  at "simple" or "bacterial." Added `program_constants.
  PLANET_CLASS_MAX_LIFE_STAGE`, a per-class ceiling on the highest
  milestone that class's description is consistent with (Class E capped at
  the most minimal stage, F/G at simple/bacterial, L at vegetation/
  multicellular; classes with no explicit life-complexity wording in their
  description stay uncapped), and threads the planet's class through so
  both the natural age-based roll and a forced `INTELLIGENT_LIFE=True` are
  capped by it. Includes a hard invariant assertion and dedicated
  regression tests.

## [5.46.6] - 2026-09-23

### Fixed
- **A binary system's own per-star property tables didn't render.** Each
  star's `###`/`===` section header was joined to its property table by a
  single newline instead of a blank line, so `html/lib/mdconvert.py`'s
  blank-line block splitter lumped the heading and table into one block --
  neither a valid single-line heading nor a valid table -- and rendered it
  as one escaped paragraph of literal `#`/`|` characters instead of a real
  heading plus table. Happened once per star, so every binary system showed
  two broken blocks. Fixed in `systemData.py`'s close- and wide-binary
  rendering paths.
- **Generated systems could place two asteroid belts overlapping each
  other, or a planet's orbit inside an asteroid belt.** Two independent
  causes: (1) a wide (S-type) binary's cross-star clearance check
  unconditionally skipped a trailing asteroid belt when finding a star's
  "outermost" object, so belt-vs-belt (or belt-vs-planet) overlap between
  the two stars' disks was never checked at all; (2) a forced
  habitable-world/explicit-class placement drew its distance uniformly
  within the target zone with no awareness of already-placed belts, and
  could land inside one, or otherwise leave the object list no longer
  sorted by distance -- which the existing overlap correction only ever
  checks between immediate list neighbors, so a resulting overlap with a
  non-adjacent belt went uncorrected. Fixed by making cross-star clearance
  belt-aware (using a belt's own outer edge and a fixed minimum-separation
  threshold, since a belt has no mass for the existing Hill-radius
  criterion to apply to) and by making zone-forced distance selection
  avoid already-placed belt spans. Added an independent, all-pairs overlap
  invariant check (not a re-derivation of the existing correction's own
  formula) to catch any future regression of this kind.
## [5.46.5] - 2026-09-20

### Fixed
- **`src/html/search.py` and `GET /api/search` crashed with a 500 error on any
  text search.** In `queryDb.py`, the SQL LIKE clauses across all five text search
  helpers (`_search_result_sectors`, `_search_result_systems`,
  `_search_result_stars`, `_search_result_planets`, `_search_result_moons`) used
  Python string literals `"ESCAPE '\\'"`. In Python string literals, `"\\"`
  resolves to a single backslash (`\`), passing `ESCAPE '\'` to MySQL. In
  MySQL/MariaDB, `\` is an escape character in string literals, so `\'` was
  parsed as an escaped single quote that left the string literal unclosed,
  causing MySQL syntax error 1064 and triggering a 500 Internal Server Error in
  the Flask API client (`planetGen API error (500): internal server error`).
  Updated all five queries to `"ESCAPE '\\\\'"` so MySQL receives `ESCAPE '\\'`
  and correctly evaluates the escape character as a single literal backslash.

## [5.46.4] - 2026-09-19

### Fixed
- **`generate.py sector`/`galaxy` crashed and discarded an entire sector's
  worth of already-generated systems once the sector ran out of physical
  room.** `SpaceSector.add_system`'s Hill-sphere-based random placement
  raises `ValueError` once a sector's cube has no space left for another
  system without overlapping an existing one's Hill sphere -- a real,
  physically expected outcome once `--num-systems`/`--density` (compounded
  by `--min-habitable` forcing extra large, larger-Hill-sphere stars) asks
  for more systems than a sector can hold at realistic stellar spacing
  (reproduced with `--num-systems 60` in a default 11.5 ly sector, which
  only fits ~35-44). `generate_sector` let that exception propagate,
  aborting the whole run and throwing away every system already
  generated -- including the expensive planet/moon generation behind each
  one. It now stops as soon as one system can't be placed, logs how many
  of the requested systems it actually placed, and returns the sector with
  whatever fit instead of crashing or generating further systems that were
  never going to fit either.

## [5.46.3] - 2026-09-19

### Fixed
- **`generate.py sector --debug` gave no way to see where a sector's generation
  time actually went.** Added `stellarObjects.log.timed_phase`, a debug-only
  context manager that logs `"<label>: <elapsed>ms"` (timestamped, like every
  other `--debug` line) around `generate_sector`'s and `StarSystem.__init__`'s
  major phases -- config building, each system's own generation and placement,
  phenomena generation, and the star/planet/comet/life-data passes inside each
  system -- so a `--debug` run's own output doubles as a per-phase profile with
  no separate profiling flag needed. Also added an opt-in benchmark
  (`PLANETGEN_RUN_PERF_BENCHMARK=1 pytest src/tests/test_sector_generation_perf.py`)
  that captures and aggregates those phase records across several generated
  sectors, cross-checked with a `cProfile` run -- this is what surfaced the
  sector-capacity crash fixed in the next release.

## [5.46.2] - 2026-09-19

### Fixed
- **`generate.py galaxy`'s random-start mode no longer shows a bogus
  single-sector progress bar.** The previous release ([5.46.1]) restored
  the "Sectors" progress bar for `galaxy`'s three modes, but random-start
  mode's own seed sector was wrapped in its own `"Sectors (random start)"`
  task with `total=1` -- a bar that goes straight from 0 to 100% in one
  `advance` call before it can render a meaningful rate or ETA, so in
  practice it was just a single completed bar sitting on screen, not a
  genuine progress display. That was a mistake in how [5.46.1] restored
  the bar: the whole point of the restore was the *galaxy*-level bar
  (sectors being filled across a shell/neighborhood/random-start's
  surrounding radius), not a bar for generating one sector. Random-start
  mode now generates its seed sector with no progress task of its own and
  falls straight through to `run_local_neighborhood`, whose own
  `"Sectors (local neighborhood)"` task (covering every sector generated
  around that seed) is the only bar shown -- matching `--shell`/
  `--center-sector` mode, which never had this problem.

## [5.46.1] - 2026-09-19

### Changed
- **Restored the progress bar for `generate.py galaxy` (`--shell`/
  `--center-sector`/random-start), but not for `generate.py sector`.** A
  prior release removed the progress bar from both commands entirely
  while chasing a flicker/scrolling bug; the bar itself (routed through
  `progress.console.print` so it stays pinned at the bottom with no
  flicker -- see `_generation_progress`'s own docstring) was fine and is
  genuinely useful for a `galaxy` run, which can mean thousands of
  sectors. `run_sector`'s own `--num-sectors` loop stays bar-free, since
  that run is normally short enough that a bar added more noise than it
  was worth.

## [5.46.0] - 2026-09-19

### Changed
- **A qualifying sector can no longer come out completely empty.** A
  sector's own system count and each exotic-phenomenon type's own count
  are independent Poisson draws, so all of them landing on zero at once
  is a real, expected outcome -- more likely the closer local density
  sits to the 1-star-per-sector qualification threshold (e.g. ~13% at a
  mean of 2 systems, ~37% right at the threshold itself, mean 1) -- but a
  sector with nothing in it at all isn't useful to anyone visiting it.
  `generate_sector` now force-adds exactly one system when a
  `--density`-driven sector's own draws (system count, every phenomenon
  type) all came back empty. Only applies when the count came from
  `--density` (explicit or `_BatchDensity`-resolved, as in `galaxy`
  mode); an explicit `--num-systems 0` (including on the plain `sector`
  subcommand) is a deliberate request this never second-guesses.

## [5.45.1] - 2026-09-19

### Changed
- **Cut database round trips for name-uniqueness bookkeeping** -- profiling
  a real save (`SHOW GLOBAL STATUS` query counters before/after) found
  `reserve_body_name`/`confirm_body_name` (run once per planet and once
  per moon -- hundreds of times in a single sector) issuing 4 separate
  SELECT/INSERT/UPDATE round trips each, the large majority of a saved
  sector's total query count. `reserve_body_name` now combines its two
  stateless cross-level collision checks (`sector_name_registry`,
  `system_name_registry`) into one query via `UNION ALL` instead of two
  sequential ones (its `body_name_registry` lookup keeps its own
  unchanged `FOR UPDATE` lock, still queried separately). `confirm_body_name`/
  `confirm_system_name`/`confirm_sector_name` now upsert
  (`INSERT ... ON DUPLICATE KEY UPDATE`, relying on each registry
  table's own `base_name UNIQUE` constraint) instead of a SELECT to
  decide between an INSERT and an UPDATE. Measured ~30-35% fewer total
  queries per planet/moon across repeated benchmark runs (same
  generation code, a throwaway database, and `SHOW GLOBAL STATUS`
  counters compared before/after this change) -- since DB writes
  dominate real generation time (a separate finding: ~70 systems/sec of
  pure in-memory generation vs. ~11-16 systems/sec once real database
  writes are included), this directly speeds up `sector`/`galaxy`
  generation. Every existing name-uniqueness test (collision handling,
  occurrence-count tracking, suffix/diminutive progression) still passes
  unchanged -- the upsert's `ON DUPLICATE KEY UPDATE` branch updates the
  exact same columns (`occurrence_count`, `suffix_index`/
  `diminutive_index`) the old SELECT-then-UPDATE branch did, and never
  touches `first_body_id`/`first_star_system_id`/`first_sector_id`,
  matching the old code exactly.

## [5.45.0] - 2026-09-19

### Added
- **`generate.py galaxy`'s random-start mode (no `--shell`/`--center-sector`)
  gained `--min-start-density`** -- requires the randomly chosen starting
  sector's own real `relative_density` (the same "expected" figure printed
  alongside each saved sector) to be at least the given value before
  accepting it, retried the same way an already-occupied or otherwise
  non-qualifying address already was. Lets an operator skip past the
  galaxy's own vast, sparse outskirts (a plain random start lands there
  most of the time, since a volume-weighted draw favors them) and start
  somewhere with real content to look at -- e.g. `--min-start-density 1.0`
  for at least as dense as the galaxy's own real local density. Only
  applies to random-start mode, and can't be combined with
  `--density`/`--num-systems` (those override every position's density
  uniformly, leaving nothing per-position to compare against).

## [5.44.1] - 2026-09-19

### Changed
- **Removed the "Sectors"/"Sectors (shell N)"/"Sectors (local
  neighborhood)"/"Sectors (random start)" progress bar entirely** from
  `generate.py sector`/`galaxy` -- the previous release only removed the
  nested per-sector bar and tried to fix the outer one's flicker by
  routing prints through it, but the outer bar itself was still visible
  and still wasn't what was wanted. `run_sector`/`run_shell_batch`/
  `run_local_neighborhood`/`run_random_start`/`run_galaxy` no longer take
  or build a `rich.progress.Progress` at all -- every status line is a
  plain `print` again, and `generate.py` no longer imports `rich.progress`.

## [5.44.0] - 2026-09-19

### Changed
- **Removed the per-sector "systems in this sector" progress bar** that
  `generate.py sector`/`galaxy` nested under the outer "Sectors" bar --
  most sectors, especially since the density-gating fix above, hold
  anywhere from zero to a handful of systems, and system generation
  itself is fast, so a bar that flashed on and off again within a single
  frame for nearly every sector added visual noise without conveying
  anything a viewer could actually track. `generate_sector`/
  `generate_and_save_sector_at` no longer take a `progress` parameter at
  all.
- **Fixed the remaining "Sectors" bar fighting with the status text
  printed alongside it, which is what actually caused the flicker/
  scrolling** -- `run_sector`/`run_shell_batch`/`run_local_neighborhood`/
  `run_random_start` now print every "Saved sector ..." status line (and
  the per-sector system/phenomena/density summary) via
  `progress.console.print(...)` instead of the builtin `print`, the
  correct way to write to the console alongside a live `rich.progress.
  Progress` display. Printing directly to stdout while `Progress`'s own
  `Live` region is active fights with its redraws -- each raw `print`
  forced the bar to erase itself, scroll up with the new text, and get
  redrawn at the bottom again, which is what showed up as flicker/
  scrolling on a real terminal. Routed through the shared console
  instead, rich prints each status line safely above the live region and
  leaves the bar itself pinned at the bottom, redrawn in place with no
  flicker -- verified against a real pseudo-terminal (`script`), and
  unchanged (a single plain line at the end, no ANSI live redraw) when
  stdout isn't a real terminal at all (piped to a file, a CI log, etc.),
  which `rich.Console` already detects and handles on its own.

## [5.43.1] - 2026-09-19

### Added
- **`test_galaxy_gen.py` now has two end-to-end tests that run against a
  *real* `generate.py plan` skeleton** (the actual `find_shell_bands`
  scan, not the file's existing `_seed_skeleton` shortcut) rather than an
  explicit `--num-systems`/`--density` that bypasses `_BatchDensity`'s own
  density/gating logic entirely -- every density-related test before this
  did one or the other, so none of them actually exercised the "run
  `plan`, then `galaxy` with neither flag given" workflow the previous
  release's empty-sectors regression slipped through.
  `test_random_start_neighborhood_matches_the_real_skeleton_plan` runs
  `galaxy`'s own default random-start mode -- pick a location, generate
  the nearest sectors out to `--radius-pc` (trimmed to 25 ly here, from
  the real default of 100 ly, to keep the test fast; `-planets` forced so
  each system skips its own planet/moon tree, since system *count* is
  what's under test) -- then independently recomputes, against the real
  stored skeleton, whether every candidate slot in that neighborhood
  should have been saved, and checks the aggregate system count generated
  is within a statistical band of what the plan's own density predicted.
  `test_shell_batch_generates_nothing_beyond_the_real_skeletons_outer_edge`
  is its deterministic companion: a shell chosen well past the real
  skeleton's own discovered edge must generate exactly zero sectors.
  Both fail against the pre-fix code (confirmed by hand, reverting
  `generate.py` locally and re-running).

## [5.43.0] - 2026-09-19

### Fixed
- **`generate.py galaxy` (`--shell`, `--center-sector`, and the default
  random-start mode) generated and saved a real, empty (0 systems, 0
  phenomena) sector row for every not-yet-occupied slot it visited, once
  `generate.py plan`'s skeleton existed to drive per-sector density.**
  Only the lazy, visit-triggered `ensure_sector_generated` was actually
  gating generation on `galaxy_shell_band`'s stored candidate bands and
  the exact `predicted_star_count >= 1.0` threshold; the three batch/
  neighborhood modes in `run_shell_batch`/`run_local_neighborhood`/
  `run_random_start` never consulted either, so every slot below that
  threshold -- most of a realistic galaxy's volume, off the spiral arms/
  disk plane -- still got a sector saved, almost always with 0 systems
  once its (correctly tiny) relative density was fed through the Poisson
  draw. A galaxy shell/neighborhood run could come back "full of empty
  sectors" as a result, especially the default random-start mode, whose
  volume-weighted starting-shell pick favors the sparse outer galaxy.
  `_BatchDensity.resolve` (`generate.py`) now applies the same
  band-then-exact-density gate `ensure_sector_generated` already used,
  returning `None` for a non-qualifying slot so every caller skips it
  instead of generating and saving it; the admin web UI's "generate more
  sectors around this one" action (`generate_sector_neighborhood`) picked
  up the same skeleton-driven density and gating, having previously used
  a flat `--num-systems 10` for every sector regardless of position at
  all.
- Every sector-generating command (`sector`, and `galaxy`'s `--shell`/
  `--center-sector`/random-start modes) now prints, right after saving a
  sector, how many star systems of each spectral class and phenomena of
  each type it actually holds, plus that sector's actual vs. expected
  star density (`1.0` = real local stellar density) -- enough to
  sanity-check a generation run from its own console output, without a
  separate database query.

## [5.42.0] - 2026-09-19

### Added
- **A 3D body-preview sphere on the System Map.** Clicking a planet or
  moon in `system.py`'s System Map now also redraws `#sysmap-preview`: a
  small rotating three.js sphere (reusing the same vendored build the
  Sector Map uses) shaded by the body's own class color, banded with a
  tilted ring for a gas giant, and wrapped in a fresnel-glow atmosphere
  shell -- tinted by surface temperature -- whenever the body actually has
  one. The info panel also gains "Atmosphere," "Surface composition," and
  "Surface temperature" fields (`planets`/`moons.atmosphere`/
  `composition`/`surface_temperature_k`, already generated and stored,
  just not previously surfaced here). The true-position SVG diagram itself
  is unchanged -- this is an appearance preview alongside it, not a
  replacement.

## [5.41.0] - 2026-09-19

### Changed
- **The Sector Map is now a real WebGL scene instead of a CSS 3D
  illusion.** `html/lib/starmap.py` no longer builds one `<div>` per star/
  outline edge/compass arrow positioned via CSS `transform-style:
  preserve-3d` -- it now serializes the same position/size/color/label
  data it always computed into a `<script type="application/json">`
  block, and a rewritten `html/static/sectormap.js` renders it with
  three.js (a real perspective camera, GPU-billboarded sprites for stars/
  nebulae/asteroid fields/black holes/neutron stars, and a wireframe
  outline for the sector's wedge or fallback cube) -- proper perspective/
  occlusion, a scale bar that now accounts for the panel's own responsive
  size instead of assuming a fixed 320px scene, and a `<noscript>` link
  list plus a hidden screen-reader-accessible button list (a canvas has
  no focusable children of its own the way the old per-star `<div
  role="button">`s were) so the sector's systems/phenomena stay reachable
  without JavaScript or with a keyboard/screen reader alike. three.js is
  vendored at `html/static/vendor/` (bundled and minified from the `three`
  npm package, not loaded from a CDN) so `html/lib/page.py`'s existing
  `Content-Security-Policy: default-src 'self'` needs no exception for it.
- **Every star/phenomenon marker on the Sector Map now navigates via
  `data-nav-target`/`data-nav-params` (`static/navform.js`) instead of a
  plain `href`**, catching the WebGL rewrite above up to the
  no-address-bar-params convention `5.40.0` (below) introduced for the
  rest of `html/` after this branch had already diverged from it.

## [5.40.0] - 2026-09-19

### Changed
- **Every navigational link in `html/` now posts its parameters as
  hidden form fields instead of putting them in a `<a href="page.py?
  db=...&id=...">`'s query string** -- `db`, a sector/system/phenomenon
  id, a search filter, a wiki-upload/admin-action field, and the like no
  longer show up in the browser's own address bar. `lib/fmt.py`'s new
  `post_link` builds a same-effect, no-JS-required `<form method="post">`
  submit button, styled (`static/style.css`'s `.link-btn`) to be visually
  indistinguishable from the plain link it replaces; `lib/page.py`'s new
  `nav_params`/`nav_multi_params` are what a page reads a followed link's
  params back with (a POST body when present, else the GET query string,
  so a bare `QUERY_STRING`-only smoke test still works). The one
  exception is a Galaxy Map/Sector Map/NAV Map marker plotted inside an
  `<svg>` (a `<form>` can't nest inside one) -- those still navigate via
  a real, focusable `<a>`, now carrying `data-nav-target`/
  `data-nav-params` (`fmt.data_nav_params`) instead of an `href` query
  string, intercepted by the new `static/navform.js` (loaded on every
  page) to post the same hidden form a click on any other link would.
  Every such marker still has a plain, no-JS-required row in a table
  below its own map, so a marker click is never the only way to reach
  something. `index.py` no longer redirects to `browse.py?db=...` for
  this deployment's one database (a redirect's `Location` URL would
  itself show `db` in the address bar) -- it calls `browse.handler`
  in-process instead and renders the result directly, guarded by
  `browse.py`'s own `if __name__ == "__main__":` so `browse.py` reached
  directly is unaffected. This makes every page un-bookmarkable/
  un-shareable by URL and, for a map marker specifically,
  JavaScript-dependent -- a deliberate trade-off (see `lib/page.py`'s
  module docstring) for keeping database names, record ids, and search
  terms out of browser history, address bars, and referrer headers.

## [5.39.1] - 2026-09-19

### Fixed
- **CI was red on every run.** `test_db_persistence.py` still asserted
  `_db.SCHEMA_VERSION == 22` (and `migrate_database(...) == 22`) in
  sixteen places, left over from before the v23 wiki-publishing and v24
  name-uniqueness migrations bumped `SCHEMA_VERSION` to 24; each
  assertion now just compares against `_db.SCHEMA_VERSION` instead of a
  stale literal. Fixing that surfaced a second, previously-masked bug in
  the same file: eight of those tests reset `schema_migrations` to an
  older version to replay later migrations, but never dropped the v22
  search-index migration's indexes first, so replaying it against a
  freshly-bootstrapped (already-v24) test database hit a "Duplicate key
  name" error -- a new `_drop_v22_search_indexes` helper (alongside the
  existing `_drop_v20_trajectory_columns`/`_drop_v21_phenomenon_columns`)
  fixes that.
- **Two `test_galaxy_gen.py` shell-batch tests raised `TypeError`.** Their
  `_fake_generate_sector` test doubles didn't accept the `progress`
  keyword the nested-progress-bar feature added to the real
  `generate_sector`, so `monkeypatch`ing it in broke as soon as
  `generate.py galaxy --shell` started passing one.
- **A MySQL connection-pool leak was exhausting the test server's
  `max_connections` partway through a full test run.** `_db.py` caches one
  `PooledDB` per distinct connection config in a module-level dict that's
  never evicted -- fine for the handful of long-lived databases a real
  deployment ever points at, but the test suite's own `mysql_config`
  fixture hands every single test a uniquely-named throwaway database, so
  each test left one more pool (and its `mincached` real connection)
  behind for the rest of the process's life. A new `_db.close_pool()`,
  called from that fixture's teardown once its database is dropped for
  good, closes and discards the pool immediately instead.

## [5.39.0] - 2026-09-18

### Added
- **Galaxy-wide name uniqueness.** `stellarObjects/nameUniqueness.py`
  tracks every sector/system/planet-or-moon base name ever generated
  (`sector_name_registry`/`system_name_registry`/`body_name_registry`,
  schema v24) and decorates a colliding name instead of letting two rows
  anywhere in the database share a display name -- Greek/Roman letters
  for a sector or system colliding with its own kind, a diminutive prefix
  for a system colliding with its sector, and a companion suffix for a
  planet/moon colliding with anything. `_db.py`'s `insert_sector`/
  `insert_star_system`/`insert_planet`/`insert_moon` all consult and
  update these registries now. `src/dedupeNames.py` is a new one-off
  script to decorate any duplicate names an existing database already
  has from before this feature existed.
- **A "generate more sectors around this one" admin action** on
  `html/sector.py`, for any already galaxy-placed sector -- fills in
  every not-yet-generated sector within a 100 ly sphere around it
  (`POST /api/sectors/<id>/generate-neighborhood`), the same
  local-neighborhood logic `generate.py galaxy --center-sector` already
  used from the CLI, now reachable from the web interface.
- **Live progress bars for `generate.py sector`/`galaxy`.** Both now show
  nested `rich.progress` bars (elapsed and estimated-remaining time) --
  an outer "Sectors" task for galaxy's shell-batch/local-neighborhood/
  random-start modes (and sector's own `--num-sectors` loop), and an
  inner "systems in this sector" task nested under it.

### Changed
- **Removed the separate write-capable MySQL config.** `config.json`'s
  `mysql_write` section (and the `PLANETGEN_MYSQL_WRITE_USER`/
  `_PASSWORD` env vars layered over it) is gone -- the Flask API's
  `WRITE_MYSQL_CONFIG` now simply reuses `MYSQL_CONFIG`, so there's a
  single account (`config.json`'s `mysql` section /
  `PLANETGEN_MYSQL_*`) for the generation CLIs, `install.sh`/
  `migrateDb.py`, and the API's reads and writes alike -- give that one
  account whatever grants the most demanding caller needs.

## [5.38.0] - 2026-09-18

### Fixed
- **A binary system's stored `markdown_content` mixed wikitext template
  blocks into it (and vice versa for `wikitext_content`)** -- only ever
  for the secondary star's own "Star Data" table/age sentence, and only
  for a binary system (`+binary_system`), never a single star.
  `StarSystem.__init__` gave the secondary star its own `copy.deepcopy`d
  `SystemConfig` (to force `LARGE_STAR` off on it without affecting the
  primary), which also forked `MARKDOWN` into its own disconnected copy.
  `_db.py.insert_star_system` renders *both* formats from one generated
  system by toggling `system_config.MARKDOWN` and calling `str(star_system)`
  again -- a toggle that only ever reached the primary/shared config, never
  the secondary's own deep copy, leaving its data table stuck rendering in
  whichever format was current at generation time regardless of which one
  was actually being requested afterward. The secondary star now shares
  the system's own `SystemConfig` object throughout (`LARGE_STAR` is still
  forced off for it, just transiently, restored right after).
- **A binary system's own name/title had "Binary System" literally baked
  into it** (e.g. "Sol Binary System" instead of "Sol") for a 'close'
  (P-type) pair -- `doubleStar.BinaryStarProxy.name` (the proxy's `.name`
  stands in for the whole system's own identity, stored as
  `star_systems.name`) now takes the primary star's own bare name, the
  same convention a 'wide' (S-type) pair's system name already used.
- **The Sector Map and System Map labeled a binary's two stars "Primary"/
  "Secondary"** instead of a real name -- both now read "&lt;name&gt; A"/
  "&lt;name&gt; B", matching the generator's own existing convention for
  the secondary star's *stored* name (`"<primary name> B"`).
- **A wide (S-type) binary's own System Map diagram never showed the
  secondary star's planets at all, and text/markers routinely rendered
  cut off or overlapping** -- both stars' planets, and the pair's own
  real separation (routinely tens to thousands of AU, per
  `wideBinary.py`), used to share one log-scaled radial pixel budget.
  Since a wide pair's real separation is so much larger than either
  star's own planetary system, that shared scale either crushed both
  stars' planets down near the frame's center to make room for the real
  separation, or pushed the stars themselves (and their planets' label
  text) toward the frame's outer edge and straight off the visible
  canvas. `lib/systemmap.py` now gives the primary its own full-budget
  "system" scene (its own planets only) with a companion marker for the
  secondary that swaps to the secondary's *own* full-budget scene when
  clicked -- the same "drill into it" pattern a planet with moons already
  used, one level up.
- **The Sector Map's default zoom could shrink every star dot into an
  illegible, barely-visible smear** for a galaxy-placed sector -- a
  sector's on-shell wedge wireframe is deliberately allowed to draw well
  past the fixed-size scene (so it can be dragged/zoomed into fully), but
  the *default* zoom used to fit that wedge's own extent alongside every
  star/cloud's, in one combined list. A wedge can be many times wider
  than the scene regardless of how tightly clustered the sector's own
  stars actually are, so a single oversized wedge could drag the whole
  default zoom down to its own floor, shrinking every star dot along with
  it -- confirmed by rendering a realistic sector this way and finding
  its dots reduced to a handful of barely-visible pixels, easily read as
  "almost nothing rendered". The default zoom now fits the real content
  (every star/cloud) on its own whenever there is any, falling back to
  fitting the wedge/cube outline only for a genuinely empty sector (the
  original problem that outline-fitting behavior was written to fix).
- **`GET /api/nav`'s cross-sector ("galaxy" scope) route could time out**
  once the generated galaxy grew large. `navGraph.build_knn_adjacency`
  built its routing graph with an O(n&sup2;) "compare every point to every
  other point" pass -- unnoticeable for one sector's own handful of
  systems, but this same function also runs over *every* system in *every*
  galaxy-placed sector generated so far for a cross-sector route, a set
  that only ever grows as more of the galaxy gets visited/generated. Now
  builds an in-memory 3D k-d tree and queries each point's true k nearest
  neighbors through it instead (O(n log n)) -- the exact same resulting
  graph, computed roughly 10x+ faster already at a couple thousand
  systems, comfortably under a second even at 50,000.
- The Sector page now also lists every nearby exotic phenomenon
  (`queryDb.phenomena_near_sector`, the same set the Sector Map's own
  clouds/points are drawn from) in its own table below the systems one,
  each row linking to that phenomenon's `phenomenon.py` detail page --
  previously only visible on the map itself, with no plain listing.

## [5.37.0] - 2026-09-18

### Added
- **Wiki publishing is wired up.** The standalone `wikiClient` library
  (`src/wikiClient/`, unified in [5.34.0]) is now actually called: an
  "Upload to Wiki" form on both `html/system.py` and `html/sector.py`
  (admin sessions only) publishes to whichever of Wiki.js/MediaWiki is
  configured deployment-wide, backed by new `POST /api/systems/<id>/wiki`/
  `POST /api/sectors/<id>/wiki` write routes. A system publishes its
  already-generated `markdown_content`/`wikitext_content`; a sector (which
  has no persisted page of its own) gets one built fresh at upload time
  from its own current detail. `config.json` gains a `wiki` section
  (`wikijs.base_url`/`.api_token`, `mediawiki.base_url`/`.username`/
  `.password`, each independently optional -- either, both, or neither
  backend may be configured at once, letting an upload choose "the wiki
  of their choice" when both are), read by the new `GET /api/wiki-config`
  endpoint the two forms use to know which backend(s) to offer.
  `star_systems.wikijs_url`/`mediawiki_url` (present in the schema since
  [5.34.0] but never populated) and a new `sectors.wiki_url` column
  (schema v23, `stellarObjects._db._migrate_v22_to_v23`) record where each
  page ends up; once set, `html/system.py`'s Description section is
  replaced by a link to the wiki page (opening in a new tab) instead of
  the locally rendered/source view, and `html/sector.py` shows the same
  kind of link. `html/admin.py` also gains a small form to manually set or
  clear a sector's `wiki_url` directly (`PATCH /api/sectors/<id>`), for a
  sector with a hand-written page from outside this app.

## [5.36.0] - 2026-09-18

### Added
- **A sector-then-system picker on the Nav page, reachable with no
  starting system already known.** Previously `nav.py` only ever worked
  when arriving via a specific system's own "Navigate from here" button
  (`?from=<id>` required); the sidenav's new "Nav" link now reaches it
  with nothing chosen yet, and a two-step `<select>` picker (every sector,
  then every system in the chosen one -- `GET /api/sectors` then
  `GET /api/sectors/<id>`) sets `from=` the same way arriving via
  `system.py` already did. The cross-sector half of the destination picker
  (choosing `to=` once an origin is known) gets the identical two-step
  sector-then-system cascade in place of its old plain numeric
  destination-system-id field.
- **A list and detail page for exotic phenomena** (nebula/asteroid
  field/black hole/neutron star) -- this project's first per-phenomenon
  pages; previously a phenomenon had no page of its own at all, only a
  hover tooltip on the Sector Map/Galaxy Map. `phenomena.py` lists every
  phenomenon across every sector, regardless of galaxy placement (`GET
  /api/phenomena`, paginated); each row links to `phenomenon.py`'s full
  detail view (`GET /api/phenomena/<type>/<id>`, each type's own real
  columns -- a nebula's `composition`/`formation_cause`, a black hole's
  `mass_solar`/`spin`/`has_accretion_disk`, etc. -- via new `queryDb.
  list_phenomena`/`phenomenon_detail`). Both pages are linked from the
  sidenav; the Sector Map's and Galaxy Map's own phenomenon markers
  (`lib/starmap.py`/`lib/galaxymap.py`) now click through to the same
  detail page instead of only showing a tooltip.

## [5.35.7] - 2026-09-18

### Fixed
- **The Search page timed out** (`TimeoutError`/`urllib.error.URLError`
  surfaced through `html/lib/apiclient.py`, rendered as an unexpected
  error page) once the database grew past a trivial size. `GET
  /api/search` always runs its full facet-count and autocomplete query set
  up front, on every visit, regardless of whether any filter is active
  (`queryDb.search`) -- ten `GROUP BY`/`SELECT DISTINCT ... ORDER BY`
  queries, none of them backed by an index on the column they group,
  filter, or sort by (`stars.yerkes_class`, `planets`/
  `moons`.`planet_class`/`body_type`/`life_chemical`,
  `asteroid_belts.density`, and every table's own `name`), so each one was
  a genuine full-table scan/sort. New schema v22
  (`_migrate_v21_to_v22`/`schema.sql`) adds the missing indexes; a name
  *term* search (`LIKE '%text%'`, a leading wildcard) isn't sped up by any
  of them -- that would need a FULLTEXT index, out of this fix's scope --
  but the facet counts and autocomplete lists that run unconditionally on
  every visit are. `apiclient.py`'s own request timeout also widened
  15s -> 30s as a second line of defense, not a replacement for the real
  fix.

## [5.35.6] - 2026-09-18

### Fixed
- **Every galaxy-generated sector's system count was flatly stuck at 10**,
  regardless of where it actually sits in the spiral galaxy -- a bulge
  sector and a sparse outer-disk sector generated the same way.
  `ensure_sector_generated` (the visit-triggered lazy-generation path)
  already correctly drove system count from the galaxy skeleton's real
  position-based `relative_density`, but `generate.py galaxy`'s own
  batch/local-neighborhood/random-start generation -- how every sector in
  a real deployment actually gets made -- never consulted it at all: every
  sector in a run shared one flat CLI value, defaulting to `num_systems =
  10` when neither `--density` nor `--num-systems` was given.
  `validate_shared_generation_args` now leaves both unset for `galaxy`
  mode specifically in that case (`sector` mode, which has no galaxy
  position to compute a density from, is unaffected); a new `_BatchDensity`
  helper resolves each sector's own `relative_density` from the stored
  skeleton (fetched once, reused for the whole run) and feeds it through
  exactly the way `ensure_sector_generated` already does, in
  `run_shell_batch`/`run_local_neighborhood`/`run_random_start` alike. An
  explicit `--density`/`--num-systems` still applies uniformly for the
  whole run, unchanged.
- Audited the actual system-placement code path (`SpaceSector.add_system`/
  `_random_position`, via `generate_sector`'s `for system, cfg in
  zip(...): sector.add_system(...)` loop) for whether it could silently
  place fewer systems than the (now real, skeleton-driven) requested
  count -- it can't: a sector too crowded to fit the next system's minimum
  Hill-sphere separation raises `ValueError` after
  `SECTOR_MAX_PLACEMENT_ATTEMPTS` tries rather than skipping it, so an
  under-delivered density would already be a loud failure, not a silent
  one. No code change needed for this part; noted here since it was the
  other half of what was reported.

## [5.35.5] - 2026-09-18

### Fixed
- **A galaxy-placed sector's Sector Map opened looking almost empty/broken**
  -- a couple of giant wireframe edges crossing the visible crop instead of
  a wedge shape, with any star near the outline's own edge invisible
  outside the fixed, non-panning viewport. The wedge wireframe (`lib/
  starmap.py`'s `_wedge_edges_px`) is deliberately allowed to extend well
  past the fixed 320px scene (the wedge's angular patch doesn't coincide
  with a cube's flat sides), but the map always *started* at `zoom = 1`
  regardless -- confirmed by rendering the real output in a browser and
  comparing that default against manually zooming all the way out, which
  showed the exact same content correctly. `render_map_panel` now computes
  a `_default_zoom` from the actual extent of everything being drawn
  (wedge/cube vertices, every star/cloud) and starts (and "Reset view"
  returns to) that fitted zoom instead of a flat default; `sectormap.js`'s
  own `MIN_ZOOM` floor widened to match.
- **The Galaxy Map rendered as a dense, unreadable smear of overlapping
  ring labels for any sector placed far from the core**, with its own dot
  sitting right at the visible circle's edge -- reproduced directly with a
  sector at `shell_index` ~1400, which implied 157 fixed-shell-width Rings,
  each drawn as its own guide circle + label, all crammed into the same
  480px panel. `lib/galaxymap.py`'s `_rings_to_show` (which picks the
  map's *scale*, so a far sector still fits) is now decoupled from how
  many ring guides `_ring_elements` actually *draws*: past
  `_MAX_RINGS_DRAWN` (10), it switches from one guide per literal
  fixed-shell-width Ring to 10 evenly-spaced distance markers spanning the
  same range -- still real, accurate distance labels, just no longer
  cluttering the map once there would be too many to read. The common
  near-core case (few real Rings) is unaffected -- confirmed with a
  regression render.

## [5.35.3] - 2026-09-18

### Fixed
- **Site `<title>`/browser-tab text always said "planetGen"**, ignoring
  `config.json`'s own `site_name` (e.g. "Molten Aether Starmap") that
  every other page-title path already honored. `lib/page.py`'s `render()`
  hardcoded the literal string instead of calling `load_config()` the way
  `api_base_url`/`base_url` already do; `index.py`'s own "no databases"
  title had the identical hardcoded string. Both now read `site_name` from
  config.
- **The "pick a database" landing page is gone.** This project deploys as
  one branded starmap per vhost now (`config.json`'s `site_name`/
  `api_base_url`), so a picker whose choice is realistically always length
  1 just added an extra click/page load in front of every visit.
  `index.py` now redirects straight to `browse.py` for the first database
  `GET /api/databases` returns, regardless of how many exist; every other
  page's breadcrumb (`browse.py`/`galaxy.py`/`sector.py`/`search.py`/
  `system.py`) drops its now-pointless leading "Databases" link, and the
  sidenav's own "Databases" item is removed (there is no longer a picker
  page for it to reach).

## [5.35.1] - 2026-09-18

### Fixed
- **Planet/moon/star infobox fields showed literal `<sup>7</sup>` markup
  instead of a superscript 7.** `tabledisplay.py`'s scientific-notation
  formatters (`format_body_distance`/`format_star_mass`/`format_star_radius`/
  `format_star_luminosity`) emit real HTML (`"5.3 × 10<sup>7</sup> km"`),
  correct for `system.py`'s static table cells (inserted unescaped on
  purpose) but wrong for `lib/systemmap.py`'s interactive System Map: it
  carries the same strings through `data-*` attributes that
  `static/systemmap.js` reads back with `.textContent` (deliberately never
  `innerHTML`, so database-derived values can never execute as markup) --
  which shows a `<sup>` tag as literal text instead of rendering it. Added
  `tabledisplay.to_plain_text`, converting the one `<sup>N</sup>` pattern
  into real Unicode superscript digits (`10⁷`), and applied it at every
  `data-*`-building call site in `systemmap.py` (distance, mass, radius,
  luminosity); `system.py`'s own raw-HTML table cells are untouched.

## [5.35.0] - 2026-09-17

### Changed
- **Unified `systemGen.py`/`sectorGen.py`/`galaxyGen.py`/`galaxyPlan.py`/
  `phenomenonGen.py` into a single `generate.py` script.** Those five
  root-level scripts are removed; every generator in this project is now
  reached through one program and one subcommand: `generate.py
  system|sector|galaxy|plan|phenomenon [options]`. Each subcommand
  accepts exactly the option surface its old standalone script offered
  and saves to the same database -- this is a pure consolidation, not a
  behavior change. The five scripts used to import each other
  (`sectorGen.py` called into `systemGen.py`, `galaxyGen.py` called into
  `sectorGen.py`, and so on); that logic now lives together in
  `generate.py`'s own sections (system -> sector -> galaxy -> galaxy
  skeleton -> exotic phenomena -> the unified CLI itself), calling each
  other directly instead of through cross-module imports.
  `setup.py`'s `py_modules`/console-script entry point were updated to
  match (`planetgen=generate:main`, replacing the old `systemgen`/
  `sectorgen` scripts), and `src/tests/test_examples.py`/
  `test_sector_gen.py`/`test_galaxy_gen.py` (the tests that imported the
  removed modules directly) now import `generate` instead.

## [5.34.0] - 2026-09-17

### Changed
- **Unified `wikijs`/MediaWiki publishing into a single `wikiClient`
  library.** `src/wikijs/` (the standalone Wiki.js GraphQL client) is
  replaced by `src/wikiClient/`, which exposes one `WikiClient` object
  (`backend="wikijs"` or `backend="mediawiki"`) whose `create_page` works
  the same way regardless of target -- both backends are create-only and
  share one exception hierarchy (`WikiClientAuthError`/
  `WikiClientPageExistsError`/`WikiClientRequestError`). The former
  `WikiJsClient` logic moves in unchanged as `wikijs.WikiJsBackend`; new
  alongside it is `mediawiki.MediaWikiBackend`, a from-scratch, stdlib-only
  MediaWiki Action API client (Bot Password login, CSRF token, `action=edit`
  with `createonly=1`) -- this project previously had no MediaWiki API
  client at all, only a wikitext *text format* option. Nothing in the app
  calls either backend yet (still an open item, see `docs/TODO.md`); this
  is purely the shared library those still-`# TODO` call sites
  (`routes.py`, `config.py`, `appconfig.py`, `system.py`) will build on.
  Tests renamed/moved to match (`test_wikiclient_wikijs(_integration).py`)
  and a `test_wikiclient_mediawiki(_integration).py` pair added, plus
  `test_wikiclient_client.py` for the new dispatch facade.

## [5.33.0] - 2026-09-17

### Added
- **Galaxy random-start mode.** `galaxyGen.py` run with neither `--shell`
  nor `--center-sector` now picks a uniformly random (by volume, not by
  shell index -- see new `_pick_random_shell_index`) not-yet-occupied
  sector address within a real Milky-Way-scale galaxy (`--max-shell`,
  defaulting to the shell nearest new `program_constants.GALAXY_RADIUS_PC`,
  15,000 pc), generates it, then falls straight through to
  local-neighborhood mode's own logic to generate every sector within
  `--radius-pc` of it too (defaulting to new `program_constants.
  RANDOM_START_NEIGHBORHOOD_RADIUS_LY`, 100 ly, in every direction) -- so
  a bare `galaxyGen.py` with no arguments at all creates a whole small
  starmap around a fresh, randomly chosen starting point in one run.
- **Science-based exotic phenomena as part of ordinary sector
  generation.** Every generated sector (`sectorGen.py` directly, or via
  `galaxyGen.py`) now also seeds a realistically sparse population of
  `phenomenonGen.py`'s own seven phenomenon types -- black holes, neutron
  stars, nebulae, supernova remnants, rogue planets, interstellar comets,
  standalone asteroid fields -- sampled independently per type via a
  Poisson draw whose mean is a cited real (or, where flagged, a
  deliberately conservative order-of-magnitude) astrophysical rate per
  star system (new `program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM`:
  Lamberts et al. 2018 for black holes, Sartore et al. 2010 for neutron
  stars, Frew & Parker 2010 for nebulae, Diehl et al. 2006 for supernova
  remnants, Sumi et al. 2011/Mroz et al. 2017 for rogue planets), scaled
  by however many star systems the sector actually ended up with -- a
  denser sector (`--density`) gets proportionally more, and most sectors,
  realistically, get none at all, exactly as real space this size usually
  is. New `sectorGen.generate_sector_phenomena`.
- **Hill-sphere-safe placement generalized to exotic phenomena.** A
  standalone black hole/neutron star is a real, stellar-mass gravitating
  body, so `spaceSector.py`'s existing star-to-star Hill-sphere placement
  logic (`hill_radius_ly`/`required_separation_ly`) now applies to it
  too, via new `SpaceSector.add_phenomenon`/`SectorPhenomenonEntry`:
  neither can ever land within a neighboring star system's or another
  compact remnant's own Hill sphere -- and, since the same neighbor list
  now includes already-placed phenomena, a star system placed afterward
  can't land inside a black hole's/neutron star's Hill sphere either.
  Every other phenomenon type has no comparable gravitational footprint
  at this generator's own scale and is placed at a random point in the
  sector's cube instead. Persisted via schema v21: `black_holes`/
  `neutron_stars` gain the same `sector_id`/`center_x/y/z_pc`/
  `galactic_radius_pc` galaxy-frame placement shape v18 gave `nebulae`/
  `asteroid_fields` (`_migrate_v20_to_v21`), with the new position
  converted directly from each phenomenon's own real in-sector offset
  rather than independently re-randomized; the other five phenomenon
  types' own `sector_id` column -- present since v16 but never actually
  populated -- finally gets wired through `_db.save_phenomenon`. Surfaced
  in the web interface everywhere nebulae/asteroid fields already were:
  `queryDb.py`'s phenomena queries, the Sector Map (a small glowing point
  instead of a translucent cloud, since a compact remnant's own physical
  radius is negligible at this scale), and the Galaxy Map.

### Fixed
- **`SectorPhenomenonEntry.distance_to` crashed** with `TypeError:
  'SectorPhenomenonEntry' object is not iterable`. `spaceSector.
  distance_between` only recognized `SectorSystemEntry` via `isinstance`,
  so calling `.distance_to()` on the new phenomenon-entry class (added
  alongside the Hill-sphere generalization above) tried to iterate the
  entry object itself instead of reading its `.position`. Generalized to
  duck-type on "has a `.position`" instead, covering both entry classes
  (and any future one) without an `isinstance` check needing to know
  about each concretely. Caught by this release's own new test suite,
  run against a real MariaDB server rather than skipped for lack of one
  (9,265 tests, 0 failed, 0 skipped) -- not by any pre-existing test.

## [5.32.0] - 2026-09-16

### Added
- **Proper two-body (barycentric) trajectories for binary stars, planets,
  and moons.** Every orbital pair previously modeled only "the lighter
  body orbits a fixed primary" -- a poor approximation for a binary
  secondary (sampled at 0.1-0.8x the primary's mass) and for a moon near
  this generator's own mass cap (up to 1/10 its parent planet's mass,
  approaching real "double planet" ratios like Pluto/Charon). Existing
  "relative position" columns (`planets`/`moons.position_x/y/z_km`,
  `star_systems.binary_mutual_position_x/y/z_km`) are unchanged -- still
  the true separation a large amount of existing physics (insolation,
  Hill sphere, tidal locking) depends on. New columns instead add the
  ORBITED body's own small "reflex offset"/"wobble" away from its
  nominal fixed point (new shared `utils.calculate_reflex_offset`
  helper): both binary members now visibly orbit their common barycenter
  (`star_systems.binary_primary_position_*_km`/
  `binary_secondary_position_*_km`, plus the constant
  `binary_secondary_mass_fraction` and, for a 'close' pair, the
  additional `binary_planetary_wobble_*_km` from its own circumbinary
  planets), a planet-hosting star gets its own wobble from the combined
  pull of its planets (`stars.reflex_offset_*_km`), and a moon-hosting
  planet gets its own wobble from its moons
  (`planets.reflex_offset_*_km`). `_db.advance_orbital_phases` recomputes
  all of these fresh on every run (a cheap derived value, no independent
  update-guard interval of its own), introducing a new correlated-
  subquery `UPDATE` technique to sum a parent's pull from multiple
  children in one set-based statement. Persisted via schema v20 (new
  columns on the three already-existing `star_systems`/`stars`/`planets`
  tables) and a `_migrate_v19_to_v20` migration step that backfills real
  values for every pre-existing row (unlike v17's migration, every value
  here is fully derivable from data already stored).
- **Standalone Wiki.js GraphQL client for page publishing.** New
  `wikijs` package (`src/wikijs/`): `client.py`'s `WikiJsClient` creates
  pages via Wiki.js's GraphQL API using Personal API Tokens, returning a
  `WikiPage` dataclass (`id`/`path`/`title`/a constructed `url`);
  `exceptions.py`'s hierarchy (`WikiJsError`/`WikiJsAuthError`/
  `WikiJsPageExistsError`/`WikiJsRequestError`) distinguishes auth
  failures (401/403, or a GraphQL-level auth error) from a duplicate-path
  rejection (message-text matching, since Wiki.js has no stable error
  code for it across versions) from any other failure. Stdlib-only
  (`urllib.request`/`urllib.error`/`json`, no external dependency) and
  create-only -- `create_page` never checks for an existing page first,
  letting Wiki.js reject a duplicate naturally rather than racing a
  check against it. Not yet wired into the app -- TODO comments mark
  where `apiclient.py`/`routes.py`/`system.py`/`appconfig.py`/
  `config.py` will eventually call it. (Superseded in `5.34.0` by the
  unified `wikiClient` library, which folds this client in unchanged as
  one of two backends.)

## [5.31.0] - 2026-09-12

### Added
- **Galaxy-frame placement for nebulae/asteroid fields, and web-interface
  map overhaul.** Nebulae/asteroid fields previously had no location at
  all -- `sector_id` existed but was explicitly "reserved for a future
  sector-context encounter" and always NULL. `phenomenonGen.py
  --sector-id` now actually places one: a galaxy-frame sphere
  (`center_x/y/z_pc`/`galactic_radius_pc`, schema v18) centered near an
  already galaxy-placed sector (`stellarObjects._db.
  compute_phenomenon_placement`), not a sector-relative offset -- a
  nebula can span up to 200 ly, far larger than a single ~11.5 ly sector,
  so it's a real sphere that may overlap several sectors' cubes, or none.
  `sector_id` itself becomes a real (if non-authoritative) "nearest
  sector" convenience link (`ON DELETE SET NULL`, not `CASCADE`). New
  `queryDb.phenomena_near_sector` (bounding-sphere-vs-sector-cube overlap
  test) and `galaxy_placed_phenomena`, exposed via `/api/sectors/<id>`'s
  new `phenomena` key and a new `GET /api/galaxy/phenomena` endpoint.
  `html/lib/starmap.py`'s Sector Map now draws a translucent spheroid
  cloud for every nebula/asteroid field near that sector (nebula color by
  type; a mottled tan/dark texture for an asteroid field); `html/
  lib/galaxymap.py`'s Galaxy Map plots every placed one as a small fixed-
  size dot instead -- a phenomenon's own real size is shown where it
  actually fits (the Sector Map), not misleadingly as a galaxy-scale dot.
- **System Map: real body positions, not a schematic.** `html/
  lib/systemmap.py`'s System Map previously always drew every planet due
  east of its star on one fixed concentric-orbit diagram, ignoring each
  body's own real orbital angle entirely. It's now a true top-down plot:
  every star/planet/moon/belt sits at its actual angle (from
  `planets`/`moons.position_x/y_km`) and a distance from its own anchor
  that's log-scaled into one shared pixel budget spanning the whole
  scene -- a log scale is what lets a close binary's ~0.05 AU separation
  and an outer planet's 30+ AU orbit coexist on the same diagram without
  either collapsing to a point or blowing out the frame. A binary pair's
  two stars are placed at their real mass-weighted offsets from the
  system's own barycenter (new `binary_mutual_position_x/y/z_km` in
  `queryDb.system_detail`'s response, split by each star's `mass_kg`) --
  for a 'wide' (S-type) pair, this also removes the old map's "primary
  only, see the tables below for the secondary" limitation: both stars
  and each one's own independent planets now draw in the same
  true-position scene. Asteroid belts are now full rings around their own
  anchor (the natural true-position shape) rather than a one-directional
  shaded band. A once-only pairwise-repulsion pass (`_relax_markers`)
  nudges apart any two markers real placement happened to put too close
  together -- real position first, decluttering only where legibility
  actually needs it.
- **Star-bound elliptical/parabolic comets.** Adds the two comet orbit
  types the existing `InterstellarComet` (always unbound/hyperbolic)
  couldn't represent: a periodic elliptical comet and a single-apparition
  parabolic one, both bound to their star. New `keplerMotion.py` solves
  the two-body Kepler/Barker equation so an eccentric or parabolic orbit
  moves with the correct non-uniform angular speed (Kepler's second law)
  instead of a linear phase advance; new `cometData.Comet` generates
  realistic orbital elements (Jupiter-family/Halley-type/long-period
  subtypes) and scales coma/tail activity chance by perihelion distance
  rather than a flat roll. Wired all the way through: `StarSystem` gets
  an independent `comets`/`secondary_comets` list (kept separate from
  planets/asteroid belts, since a comet's distance is a continuously
  varying current position, not a fixed orbital-slot distance) and a new
  `SystemConfig.COMETS` tri-state flag; `advance_comet_orbits`
  (`updateOrbits.py`) advances each comet's orbital anomaly and
  recomputes its position on every run; `queryDb.system_detail`/
  `html/system.py` gain a parallel comets table. Persisted via schema v19
  (`comets`/`comet_composition` tables, `_migrate_v18_to_v19`).
- **Five project-wide bugs, found in a general bug-check pass:**
  `StarSystem.validate_system`'s orbital-overlap correction
  double-counted `last_planet.distance`, roughly doubling the corrected
  distance for any overlap involving an asteroid belt; `utils.
  years_to_time_string` built its total from a 365.25-day year but
  decomposed it back with a plain 365-day divisor, leaking the
  discrepancy into every displayed orbital period's "days"/"hours" (e.g.
  1.0 year rendered as "1 year and 6 hours"); the admin session cookie
  (`html/api/auth.py`) was scoped to `Path=/api`, so it was never sent
  back to the CGI admin pages served from the deployment's own document
  root, locking every admin page immediately after a successful login;
  three latent bugs in `test_query_db_planets_moons.py` itself (a
  nonexistent `StarSystem.name`, `Planet.radius_km` vs. the real
  `radius`, a missing `body_type` guard before touching
  `AsteroidBelt.moons`); and `advance_comet_orbits` now batches its
  per-row updates into one `executemany` call instead of one `execute`
  per row.

### Fixed
- **CI: schema v18's null-together CHECK on `nebulae`/`asteroid_fields`
  broke `test_migrate_v*` on real MySQL 8.0** (`pymysql.err.
  OperationalError: (3959, "Check constraint 'nebulae_chk_2' uses column
  'center_x_pc', hence column cannot be dropped or renamed.")`) --
  MySQL 8.0 refuses `DROP COLUMN` on a column an anonymous CHECK still
  references, and the migration-rollback test helper had no name to drop
  it by. Both CHECKs are now named explicitly
  (`chk_nebulae_placement`/`chk_asteroid_fields_placement`) rather than
  left anonymous; `_migrate_v17_to_v18` now also adds the same named
  constraint (via the portable `ADD CONSTRAINT ... CHECK`/
  `DROP CONSTRAINT`, not MySQL-only `DROP CHECK`), which it had
  previously omitted -- a database migrated (not freshly created) would
  otherwise silently lack this invariant's enforcement. Caught locally by
  running the full suite against a real MariaDB server (which tolerated
  the original anonymous CHECK's `DROP COLUMN` fine, masking this) rather
  than against MySQL 8.0 as CI does. The workflow's own test matrix also
  gained `fail-fast: false`, since the default had been silently
  cancelling the 3.9 job the moment 3.12 failed, hiding whether a failure
  was version-specific.

## [5.30.0] - 2026-09-12

### Added
- **Standalone asteroid fields, the seventh exotic phenomenon.**
  `phenomenonGen.py --type asteroid-field` generates a field of asteroid
  debris drifting in open interstellar space (new `asteroidFieldData.
  AsteroidField`) -- physically the same object as an in-system
  `AsteroidBelt` (density + mineral composition), just without a host
  star/orbit, reusing `AsteroidBelt`'s own composition-generation logic
  (now extracted into shared `asteroidData.generate_asteroid_composition`/
  `format_composition_summary` functions rather than duplicated).
- **Galactic-orbital motion for every standalone exotic phenomenon.**
  A black hole/neutron star with no owning system, a nebula, a supernova
  remnant, a rogue planet, an interstellar comet, and a standalone
  asteroid field are all still gravitationally part of the galaxy even
  though none is bound to any specific star -- each now gets the same
  `galactic_orbital_speed_kms`/`_period_gy`/`_phase_deg`/
  `_min_update_interval_years` quartet a lone `Star` has (new shared
  `utils.generate_galactic_orbit_fields`/`format_galactic_orbit` helpers,
  also adopted by `Star`/`BinaryStarProxy`/`BlackHole`/`NeutronStar` for
  consistency), advanced over real elapsed time by `updateOrbits.py`/
  `_db.advance_orbital_phases` the identical way a star's already is.
  `advance_orbital_phases` now returns a name-keyed dict rather than a
  positional tuple, since the set of tables it advances keeps growing.
  Persisted via schema v17 (new columns on six pre-existing tables plus
  the new `asteroid_fields`/`asteroid_field_composition` tables) and a
  `_migrate_v16_to_v17` migration step.
- **Fixed: an anchored black hole/neutron star silently lost its identity
  on reload.** `_db.load_star_system`'s single-star branch always called
  `Star.from_dict`, with no dispatch on the owning `stars.yerkes_class`
  marker (`'BH'`/`'NS'`) and no query against the `black_holes`/
  `neutron_stars` satellite tables -- a system saved via `phenomenonGen.py
  --anchor-system` reloaded (via `queryDb.py` or the Flask API) as a
  generic `Star` carrying a nonsensical Yerkes class, missing every
  remnant-specific field, and rendered with `Star`'s own paragraph text
  instead of the remnant's. New `_db._load_single_star` now dispatches
  correctly, reusing the same satellite-row-plus-base-row combination
  `_binary_proxy_row_to_dict` already does for a close binary's merged
  proxy.
- **`stellarObjects/phenomenaPlausibility.py`, a statistical anomaly
  finder for all seven exotic phenomena** -- the same two-tier design as
  the existing planet-focused `plausibility.py` (analytically-derived hard
  invariants, e.g. an event horizon radius must match the Schwarzschild
  formula for its own mass, gated by `test_phenomena_plausibility.py`;
  Tukey's-fences statistical outliers plus category-frequency comparisons
  against each phenomenon's own configured chance, reported for human
  review via the new `src/tests/phenomena_plausibility_cli.py`, never
  asserted exactly).

## [5.29.0] - 2026-09-12

### Added
- **Exotic stellar phenomena, via a new, separate `phenomenonGen.py` CLI.**
  Six phenomena -- black holes, neutron stars, nebulae, supernova remnants,
  rogue planets, and interstellar comets -- can now be generated on demand,
  each grounded in real astrophysics (Schwarzschild radius for black
  holes; NICER-measured neutron star mass/radius ranges and ATNF-catalog
  pulsar spin/field populations; the four standard ISM nebula classes;
  Sedov-Taylor blast-wave expansion for supernova remnant age/size;
  'Oumuamua/Borisov-informed interstellar comet speed/composition).
  Deliberately **not** wired into `systemGen.py`/`sectorGen.py`'s normal
  per-slot generation odds -- `StarSystem._generate_planets` never
  produces one; they're reachable only through `phenomenonGen.py`'s own
  `--type` choice (uniformly random among all six when omitted). A black
  hole or neutron star (new `compactRemnant.py`, subclassing `Star` the
  same way `doubleStar.BinaryStarProxy` does) can optionally anchor a full
  `StarSystem` via `--anchor-system` -- real pulsar planets exist (PSR
  B1257+12) -- reusing all of `StarSystem`'s existing orbit-placement/
  rendering/serialization logic unchanged; its zero-or-near-zero
  luminosity naturally collapses the habitable zone to (0, 0) AU and the
  disk-physics planet-count ceiling to zero, matching the real rarity of
  confirmed planets around compact remnants, without any special-casing.
  Persisted via six new tables (schema v16: `black_holes`/`neutron_stars`
  as satellite tables extending a `stars` row when anchored,
  `nebulae`/`supernova_remnants`/`rogue_planets`/`interstellar_comets`
  always standalone) and a `_migrate_v15_to_v16` bookkeeping-only
  migration step (the six tables are brand new, so no existing table
  needed an `ALTER TABLE`).

## [5.28.0] - 2026-09-12

### Added
- **`queryDb.py`: `planets`/`moons` CLI subcommands.** Closes the gap
  `docs/TODO.md`'s "Open items" > "Search" tracked: the `systems`
  subcommand only ever filtered by star type/sector, with no way to ask
  this CLI "every Class D planet smaller than Earth" the way the web/API
  faceted search (`GET /api/search`, `queryDb.search`) already could.
  New `list_planets`/`list_moons` (mirroring `list_systems`'s shape) take
  an exact `--class`, a `--min-radius-km`/`--max-radius-km` range, and a
  `--sector-id`/`--system-id` scope, reusing the existing
  `_append_size_clause` helper the faceted-search result panels already
  share -- new `_body_filter_clause` factors the class/size/sector/system
  `WHERE` fragment the same way `_systems_filter_clause` already does for
  `systems`, so the two CLI-side query functions can't drift out of sync
  with each other. `moons` rows also report their parent planet's name,
  since a moon's own name alone doesn't say which planet it orbits.

### Fixed
- **README.md's version badge had drifted 3 releases stale** (5.24.0 while
  `_version.py`'s `__version__` was already 5.27.0), caught by hand during
  a deploy-readiness check with no CI signal at all. New
  `src/tests/test_version_sync.py` asserts the README badge and
  `CHANGELOG.md`'s own top entry both match `__version__`, so this can't
  silently recur.
- **CI's `dependency-audit` job had no explicit `setuptools` upgrade step**,
  leaving it exposed to whatever `setuptools` version happens to ship
  preinstalled on the runner's Python image -- caught locally via
  `pip-audit` flagging `PYSEC-2026-3447` against a preinstalled 79.0.1.
  `.github/workflows/ci.yml` now upgrades `setuptools` alongside `pip`
  before installing this project's own dependencies, the same way the
  `test` job's setup already keeps `pip` itself current.

## [5.27.0] - 2026-09-12

### Added
- **S-type (wide) binary star systems.** `+binary_system` previously only
  ever generated a P-type (close/circumbinary) pair -- the two stars merged
  into one effective star (`doubleStar.BinaryStarProxy`) for planet
  placement. A new `+wide_binary`/`-wide_binary` option (random if omitted)
  selects the other real binary configuration instead: an S-type pair,
  separated by tens to thousands of AU (log-uniformly sampled, matching the
  real, roughly log-normal spread of observed wide-binary separations),
  where each star keeps its own separate identity -- its own mass,
  luminosity, habitable zone -- and hosts its own independently-generated
  planets (new `doubleStar`-sibling module `wideBinary.py`'s
  `WideBinaryPair`). Each star's maximum stable planetary orbit is capped
  by Holman & Wiegert's (1999) empirical critical-semi-major-axis formula
  for the companion's long-term perturbation, and a Gladman (1993)
  mutual-Hill-radius check additionally prunes either star's outermost
  planet if the two stars' own disks would otherwise gravitationally
  encroach on each other -- a rare safety net for tight/eccentric pairs,
  not the common case. The pair's own orbital eccentricity is sampled from
  a realistic "thermal" distribution (unlike the close pair, a wide pair
  never tidally circularizes) and feeds both the stability formula and the
  reported periapsis/apoapsis separation, though (like every other orbit
  this generator tracks) the pair's live position/phase-advance tracking
  stays circular -- a deliberate, documented simplification consistent
  with the rest of the engine. Persisted via new `star_systems`/`stars`/
  `asteroid_belts` columns (schema v15) and a `_migrate_v14_to_v15`
  migration step; also fixes a latent bug the new columns exposed in
  `_db.advance_orbital_phases`, whose old single combined `UPDATE` could
  never have advanced a wide pair's mutual-orbit phase even after this
  release, had it not been caught -- now two independently-guarded
  `UPDATE`s.

### Fixed
- **Several tests that don't pin `BINARY_SYSTEM` could intermittently fail
  once merged against [5.26.0]'s real, spectral-class-dependent binary
  chance.** A companion star landing on a system those tests otherwise
  treat as single (or as a fixed planet count) could throw off an
  unrelated assertion -- most seriously, an S-type (wide) pair's own
  `a_crit_au` stability ceiling can make a forced `HABITABLE_WORLD`/
  `ASTEROID_BELT` guarantee geometrically impossible for an especially
  luminous host (an O-type supergiant's habitable zone can sit beyond any
  sampled companion's stability limit), which `test_full_matrix.py`'s
  full-star-type sweep and several of `test_systems.py`'s tri-state-flag
  tests surfaced; `test_disk_physics.py`, `test_db_persistence.py`, and
  `test_api.py` had similar exposure via an unexpected second star's own
  planet count/DB rows. Pinned `BINARY_SYSTEM=False` in each, following
  the same precedent already established when [5.26.0] itself pinned it
  in four other tests -- binary-vs-single was always incidental to what
  each of these was actually testing.
- **CI's `test (3.9)` job failed on every PR, unrelated to whatever the PR
  actually changed.** `test_planet_physics_fixes.py` called
  `statistics.correlation`, added in Python 3.10 -- this repo's CI matrix
  still runs a `3.9` job. Replaced with a small `_pearson_correlation`
  helper (matches `statistics.correlation` exactly; verified against it
  directly) used by both `test_atmospheric_pressure_correlates_positively_
  with_gravity` and `_spearman_correlation`'s own rank-based call.

## [5.26.0] - 2026-09-11

### Changed
- **`SystemConfig.BINARY_SYSTEM` now follows the same tri-state contract
  as every other flag (`HABITABLE_WORLD`, `ASTEROID_BELT`, etc.):
  `None` (the default) is no longer treated as "always single."** It now
  rolls real chance instead, from new
  `program_constants.BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS` --
  keyed by the primary star's own spectral letter, since real
  stellar-multiplicity surveys consistently find companionship rate
  rising with primary mass rather than sitting at one flat rate: ~26% for
  M dwarfs (Duchene & Kraus 2013) up through ~44% for solar-type F/G/K
  (anchored to Raghavan et al. 2010's 46%; Duchene & Kraus's own review
  groups F/G/K together at 44+/-2%) to ~90% for O-type primaries (Moe &
  Di Stefano 2017's 94+/-14%). New `StarSystem._should_generate_binary`
  (called from `__init__`, replacing the old flat `if self.system_config.
  BINARY_SYSTEM:` check) looks this up against `self.primary_star.type[0]`
  -- run *after* the primary star already exists, specifically so its
  real, already-rolled spectral type can drive the roll. `True`/`False`
  still force the outcome exactly as before; only `None`'s meaning
  changed, from "never" to "real chance for this star." `SystemConfig.
  BINARY_SYSTEM`'s own docstring updated to match.
- Four tests that generate systems without pinning `BINARY_SYSTEM`
  (two in `test_space_sector.py`'s name/position round-trip and
  recipe-fallback-reload coverage, one more in `test_space_sector.py`'s
  file save/load round trip, one in `test_serialization.py`'s
  single-star full-object-graph round trip) were relying on the old
  "`None` always means single" behavior to keep specific expected
  names/types deterministic, or (the `test_serialization.py` one) for
  `reloaded.star is reloaded.primary_star` to hold -- true only for a
  single star, since `__init__` never repoints `primary_star` at the
  `BinaryStarProxy` it reassigns `star` to for a real binary. All four
  now pin `BINARY_SYSTEM=False` explicitly, since binary-vs-single was
  always incidental to what each was actually testing.

## [5.25.0] - 2026-09-11

### Changed
- **`StarSystem.estimate_num_objects` now derives its planet/belt ceiling
  from real protoplanetary-disk physics instead of an arbitrary curve fit
  to stellar mass.** The old formula (`BASE_MAX_SYSTEM_OBJECTS *
  (1 + log10(solar_masses))`, base 15) had no grounding in orbital
  dynamics and no relationship at all to the mutual-Hill-radius spacing
  rule `validate_system` enforces ([5.24.0]) -- two disconnected dials
  governing "how many" and "how far apart," tuned independently by feel.
  The replacement, `StarSystem._estimate_max_objects_from_disk_physics`,
  walks outward from the same inner-edge distance `_generate_planets`
  itself seeds its first slot at, and at each step: computes the local
  *isolation mass* an oligarchic-growth embryo would reach there
  (Lissauer 1993; Kokubo & Ida 2000, 2002 -- new `utils.isolation_mass_kg`,
  closed-form-solved the same way `_mutual_min_distance_au` is, since the
  embryo's own Hill radius depends on its own still-unknown mass), from
  the Minimum Mass Solar Nebula's real solid surface-density profile
  (Hayashi 1981 -- new `utils.mmsn_surface_density_gcm2`,
  `physical_constants.MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2 = 7.0 g/cm^2` at
  1 AU falling off as `distance^-1.5`, jumping `SNOW_LINE_ICE_BOOST_FACTOR
  = 30/7` beyond the snow line where ices condense), scaled for this
  star's own disk-mass budget (new `utils.disk_surface_density_scale`,
  `program_constants.DISK_MASS_STELLAR_MASS_EXPONENT = 1.8` --
  mm-continuum disk-demographics surveys, Andrews et al. 2013; Pascucci
  et al. 2016, find real disk dust mass scales roughly as
  `M_star^1.8-2.7`, not logarithmically), then advances by that same
  embryo's own mutual-Hill-radius feeding zone (the *same*
  `program_constants.MUTUAL_HILL_RADII_SEPARATION = 10` and kappa/clamp
  algebra `_mutual_min_distance_au` uses, so the count estimate and the
  spacing rule that will later constrain actual placement are provably
  consistent with each other) and counts a slot. The walk terminates at
  the disk's outer edge -- new `program_constants.
  DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER = 18` times the star's own snow
  line (new `utils.snow_line_au`, `physical_constants.
  SNOW_LINE_AU_AT_1_LSUN = 2.7`, the same `sqrt(luminosity)` shape
  `calculate_habitable_zone` already uses) -- deliberately *not*
  `star.system_perimeter` (that's the star's own galactic-tidal Hill
  sphere, tens to hundreds of thousands of AU; real disks are truncated
  far short of it by viscous spreading/photoevaporation) -- or
  `ABSOLUTE_MAX_SYSTEM_OBJECTS` isolation-mass slots, whichever comes
  first. Not every oligarch survives as a final planet: real N-body
  integrations of the subsequent giant-impact phase (Chambers 2001) show
  most merge or get ejected, so the raw slot count is scaled by new
  `program_constants.GIANT_IMPACT_SURVIVAL_FRACTION = 0.4` (tuned toward
  the middle of that literature's own range, empirically checked to keep
  a solar-mass star's typical resulting count in the same well-tested,
  playable range this generator already verified via repeated full-suite
  runs) before being returned. `estimate_num_objects`'s own override
  contract (`PLANETS`/`NUM_ORBITS`/`MAX_PLANETS`) is entirely unchanged --
  only what `max_objects` means changed. The resulting shape now tracks
  real demographics better than the old mass-scaling formula did: cool
  low-mass stars (whose smaller mutual-Hill spacing packs oligarchs more
  tightly per unit distance -- the real, observed TRAPPIST-1-style
  "compact multis favor small stars" pattern) come out *more*
  planet-rich on average than hot, luminous, high-mass stars (whose
  correspondingly larger isolation masses claim proportionally more of
  their own, larger disk per embryo, and whose short main-sequence
  lifetimes and intense UV output make real, confirmed planets around
  O/B-type stars genuinely rare) -- the reverse of the old formula's
  "bigger star, more objects" curve, and a better match to what's
  actually been observed. `BASE_MAX_SYSTEM_OBJECTS` removed (no longer
  referenced); `ABSOLUTE_MAX_SYSTEM_OBJECTS` unchanged, still the same
  hard safety cap.
- New `src/tests/test_disk_physics.py`: unit coverage for the new
  `utils` helpers directly (snow-line `sqrt(luminosity)` scaling,
  disk-density-scale monotonicity, MMSN falloff/snow-line jump,
  isolation mass matching the literature's own ~0.05-0.1 Earth-mass
  figure at 1 AU) plus `estimate_num_objects`'s override contract and the
  `ABSOLUTE_MAX_SYSTEM_OBJECTS` cap, exercised via `MAX_PLANETS=True`
  (deterministic, no `random.randint` draw) the same way test_systems.py
  already does.

## [5.24.0] - 2026-09-11

### Changed
- **Orbital spacing between adjacent planets now uses their *mutual* Hill
  radius, not either planet's own individual one.** The previous rule
  (`min_orbit_distance = 5 x this planet's own Hill radius`) was a
  reasonable approximation, but the real orbital-dynamics stability
  literature (the analytically rigorous two-planet minimum of `2*sqrt(3)`
  mutual Hill radii, Gladman 1993; a recommended ~8-10x margin for
  longer-term N-body stability, Chambers, Wetherill & Boslough 1996 and
  Smith & Lissauer 2009) expresses this in terms of the *pair's* combined
  mass and average distance instead. New `utils.mutual_hill_radius_m`
  (`((a1+a2)/2) * ((m1+m2)/(3*M_star))^(1/3)`) and
  `program_constants.MUTUAL_HILL_RADII_SEPARATION = 10` (the safer end of
  the literature's recommended range) back `StarSystem.
  _mutual_min_distance_au`, which `validate_system`'s planet-planet
  spacing check now calls instead of `max(planet.min_orbit_distance,
  last_planet.min_orbit_distance)`. Solved in closed form for the exact
  minimum distance rather than evaluated once at the pre-correction
  position and added on top: since the mutual radius depends on the
  *average* of both distances, that naive approach understates the
  requirement once the correction actually moves one of them -- confirmed
  by a real `assert_no_orbital_overlap` failure before the closed-form
  version replaced it. Reclassification (`planetPhysics.
  reconcile_zone_and_class`, [5.22.0]) can change a planet's own mass,
  which can in turn invalidate a spacing decision already made against
  its predecessor -- `validate_system` now retries the mutual-distance
  check (bounded to 3 iterations; converges in practice within 2) after
  any reclassification triggered by its own push. `Planet.
  min_orbit_distance` itself is unchanged and still single-body -- it
  remains the right tool for a *moon's* own orbital limit around its
  parent (`planetPhysics.generate_moons`), a different physical question
  from planet-to-planet spacing. `assert_no_orbital_overlap` (test_systems.py)
  updated to check the same mutual-radius formula, with a relative (not
  fixed-1e-9) tolerance -- the two independent computations of the same
  quantity can differ at the floating-point-noise level even when both
  are correct, and that noise scales with the (sometimes very large, for
  a massive star's own extreme systems) distances involved.
- **`html/search.py`/`GET /api/search` gained a min/max size filter for
  stars, planets, and moons.** `queryDb.search`'s new `sizes` argument
  (`{"star", "planet", "moon"} -> (min_km, max_km)`, either bound
  optional) filters/joins alongside the existing class/body/life-
  chemistry tags and per-entity name search, via a new shared
  `_append_size_clause` helper across `_search_result_stars`/
  `_search_result_planets`/`_search_result_moons` -- a continuous
  quantity like size has no discrete set of values to offer as a facet,
  so it's its own query-parameter pair
  (`star_min_radius_km`/`star_max_radius_km`, likewise `planet_`/
  `moon_`) rather than a tag. `GET /api/search` validates and parses
  these (`routes._parse_size_range`); `html/search.py` gained matching
  number-input fields, active-filter chips (e.g. "Planet size: 5,000–
  8,000 km"), and a Radius column on the Stars/Planets/Moons result
  tables. `queryDb.py`'s own CLI (`sectors`/`systems`/`near`
  subcommands) still has no `planets`/`moons` equivalent at all -- see
  `docs/TODO.md`'s "Open items" > "Search" for that narrower, separate,
  still-open gap.

### Fixed
- **`GET /api/health` could crash instead of returning a clean `503`.**
  `queryDb.open_readonly` raises a bare `SystemExit` (correct for its own
  CLI callers) when the configured MySQL server is unreachable --
  `SystemExit` is a `BaseException`, not an `Exception`, so left
  uncaught it would propagate straight through Flask's request dispatch
  (and the WSGI worker handling it) instead of becoming any HTTP
  response at all, from *every* route that opens a connection via
  `routes.get_db()`, not just `/health`'s own explicit check. `get_db()`
  now catches it and re-raises an `ApiError` (503 -- the same status
  `/health` already wanted to report), which the app's existing
  `ApiError` handler turns into the usual JSON error response for every
  other route, and which `/health`'s own `except Exception` catches
  directly. Regression test builds its own Flask app against a
  deliberately-unreachable config (a closed local port), so it runs
  without needing a real MySQL server the way every other API test does.
- **Sector Map star dots were plotted as if their own sector-local axes
  already ran parallel to the galaxy frame's.** The wedge outline and
  "Galactic Center" compass arrow (`_wedge_edges_px`/`_compass_html`)
  were always computed directly from galaxy-frame quantities
  (`sectors.center_x/y/z_pc`, `sector_wedge_vertices_pc`) and so were
  always correct on their own terms; star dots
  (`star_systems.position_x/y/z_mpc`) were plotted directly in the same
  scene without ever being rotated into that frame -- a design
  convention this project documents (`docs/design/
  galaxy-coordinate-system.md`'s "Cube orientation" section: radial-
  outward local `+Z`, projected-galactic-north local `+X`) but never
  actually applies at generation time. Rather than rotate stored
  positions, `starmap.py`'s new `_rotate_to_galaxy_frame` applies that
  convention at render time, computed fresh from the sector's own stored
  `center_x/y/z_pc` -- reusing `stellarObjects.sectorGeometry.
  cube_orientation`, the exact same basis that module already computes
  as the tangent-plane frame for this sector's own wedge vertices, so no
  new stored orientation column was needed. A sector with no galaxy
  placement (`center_pc=None`) keeps the previous unrotated behavior
  unchanged. New `src/tests/test_starmap.py` covers the rotation's
  length-preservation, the on-galactic-axis degeneracy case, and that
  `render_map_panel` actually renders a different on-screen position for
  a placed vs. unplaced sector.

## [5.23.0] - 2026-09-11

### Added
- **`config.json`: one unified deployment config file, replacing
  `webconfig.json`.** Every entry point in this project (generation CLIs,
  the Flask API, the `html/` CGI browser) used to read its own scattered
  `PLANETGEN_*` environment variables, each with its own hardcoded
  default -- fine per-variable, but it meant a deployment that just wants
  "one MySQL server, one account, one API base URL" still had to set half
  a dozen `SetEnv`/`EnvironmentFile` lines to get there. `webconfig.json`
  existed to solve exactly this for the web interface, but only ever
  covered `site_name`/`base_url` plus three `db_*` placeholders that
  predated the MySQL port and were never wired to anything.
  `stellarObjects.appconfig.load_config()` replaces it: a single
  `config.json` at the repo root, deep-merged onto built-in defaults, now
  covering the read-only and write-capable MySQL connections, the control
  schema name, the database-listing prefix, the API's rate limits, the
  admin cookie's `Secure` flag, the debug-page toggle, and the site's own
  name/base URL/API endpoint -- see `docs/config.md` for the full field
  list. Every `PLANETGEN_*` environment variable still works and still
  takes precedence over `config.json` (needed for, e.g.,
  `planetgen-orbits@.service`'s per-instance
  `PLANETGEN_MYSQL_DATABASE=%i`); `config.json` only adds a place to set
  the shared defaults once instead of repeating them everywhere.
  `config.json.example` (repo root) is the committed template;
  `config.json` itself is gitignored, next to the `webconfig.json` entry
  it replaces.

## [5.22.0] - 2026-09-11

### Fixed
- **`validate_system` could strand a planet outside the zone its class
  needs, and the sequential orbit-spacing loop had no outer bound.**
  Found via full `StarSystem` stress testing rather than the existing
  per-class plausibility tooling (`physical_plausibility_cli.py`), which
  only ever constructs isolated planets directly and never exercises the
  sequential placement loop or `validate_system` at all. A planet's class
  was chosen once, early, from its *initial* estimated distance -- but
  `validate_system`'s orbital-overlap correction could later push that
  distance arbitrarily far out (Hill-radius-based minimum spacing scales
  with a planet's own distance, so it compounds geometrically across a
  many-planet system) without ever re-checking whether the already-chosen
  class still made physical sense there. Measured before this fix: 10% of
  a 400-system sample had at least one ecosphere-class planet (M/K/N/etc.)
  stranded outside its own zone -- e.g. a "Class M, Earth-like world"
  ending up hundreds of thousands of AU out at ~30K -- concentrated almost
  entirely in O/B-type stars (92% of occurrences), none in G/K/M dwarfs.
  - **`planetPhysics.reconcile_zone_and_class(planet, primary_mass_kg,
    distance_override=None)`**: re-derives a body's zone from its
    *current* distance and, only if its existing class is no longer valid
    there, regenerates the class and everything derived from it
    (composition/radius/mass/density/atmosphere/period/gravity/orbital
    motion) -- the same thing real orbital migration does to a body's
    actual final conditions, not just its originally-assumed ones.
  - **`StarSystem._reconcile_moved_planet`** calls it every time
    `validate_system` moves a top-level planet's distance, for the planet
    itself and each of its moons: a moon's zone is always its parent's
    (`generate_moons`' `zone_override`), and a reclassified parent can
    come out with a different mass, which changes a moon's own
    period/position/rotation (Kepler's third law around the parent) even
    when the moon's own class didn't need to change.
  - **The sequential placement loop never bounded how far a system could
    grow.** `StarSystem._generate_planets` (split out of `__init__` so it
    can be retried -- see below) now stops adding slots once the next
    one would land beyond `star.system_perimeter` -- the star's own Hill
    sphere *relative to the galaxy*, already computed for an analogous
    purpose elsewhere (`spaceSector.py`, keeping neighboring systems'
    spheres of influence from overlapping) but never enforced during
    planet placement itself. A system that runs out of stable room this
    way ends up with fewer planets, the same outcome a real
    protoplanetary disk of finite extent would produce, rather than
    letting the geometric compounding above run unbounded.
  - **Reconciliation can occasionally reclassify away the specific body a
    requested `HABITABLE_WORLD`/`ASTEROID_BELT` guarantee was relying on**
    (previously silently masked by the same bug -- an invalid class
    sitting in the wrong zone still counted as satisfying it).
    `StarSystem.__init__` now retries the whole placement (same star,
    fresh positions and object count) up to the new
    `program_constants.MAX_SYSTEM_GENERATION_ATTEMPTS` (8) when a
    requested guarantee isn't met afterward, rather than silently
    dropping it. A smaller, complementary fix
    (`StarSystem._distance_within_zone_with_margin`) reserves a safety
    margin against the fixed `MIN_ASTEROID_BELT_SEPARATION` nudge when
    placing a forced-habitable or explicit-slot-class planet, reducing
    (not eliminating -- a large neighboring planet's own Hill-radius push
    has no fixed size to margin against, which is what the retry loop is
    for) how often the retry actually triggers.
  - Verified via a 500-system randomized stress test (mixed
    `HABITABLE_WORLD`/`ASTEROID_BELT`/`BINARY_SYSTEM` flags): zero orbital
    overlaps, zero misplaced ecosphere-class planets, zero guarantee
    failures, repeated across multiple runs.

## [5.21.0] - 2026-09-11

### Added
- **Ecosphere-zone classes now generate at a class-appropriate distance
  within the habitable zone, instead of every class sharing the same
  distance-blind draw.** Resolves `docs/TODO.md`'s "Class K (Mars analog)
  ... generated at the same zone-midpoint orbital distance as Class M"
  open item -- this generator previously picked a planet's *class* only
  after its *distance* was already fixed by unrelated orbital-spacing
  logic, so real-world position (Venus close-in, Mars farther out) had no
  influence on which class actually generated where. `program_constants.
  PLANET_CLASSES` gains a new per-class `zone_position_mode` (0.0-1.0,
  "how far through the zone's `[inner, outer]` AU range this class's real
  or reasoned analog sits") on every class with a single, fixed
  identity within the ecosphere zone: E/F/G/H/K/L/M/N/O/P/V. Class Q
  (eccentric orbit, extreme temperature swings) deliberately has none --
  it has no single fixed position by its own flavor. `planetPhysics.
  generate_planet_properties` reads it once a planet's class is settled
  and redraws `planet.distance` there via `utils.sample_bounded_bell`
  (the same bounded-bell-curve mechanism `size_mode` already uses for
  radius), for ordinary planets only -- explicitly skipped for moons,
  since a moon's own `distance` is its orbit around its *parent planet*,
  not an AU-scale position within the star's own zone `planet.
  habitable_zone` describes; redrawing it there would corrupt it, not
  correct it. `StarSystem.validate_system`'s existing orbital-overlap
  correction absorbs whatever reordering a class-biased redraw causes
  against already-placed neighbors, the same way it already absorbed
  `calculate_distance_for_class`'s explicit-slot distance nudging.
  `stellarObjects.plausibility._extract_record` now also reports
  `distance`, letting a new `test_planets.py` regression suite verify
  the bias directly (per-class mean zone-fraction close to its declared
  mode; Class N < Class M < Class K in mean orbital distance; a moon's
  distance is provably untouched).

### Changed
- **Class K and Class N retuned now that they're placed at a real,
  class-appropriate distance instead of sharing Class M's midpoint.**
  Both carried an explicit "tuned to compensate for the wrong distance"
  comment (see above) -- with `zone_position_mode` now doing the
  distance part of the work for real, their climate ranges needed
  re-deriving via `climate_tuning_cli.py` rather than staying tuned
  against the old, distance-blind placement:
  - **K (Mars analog)**: `albedo_range` raised slightly (0.34-0.42 ->
    0.36-0.44) and `atm_density_range` widened/raised (0.012-0.025 ->
    0.022-0.042) to fit the new, colder starting point. Mean
    surface_temperature ~214K (real Mars ~210K, +1.9%, was +9.9% before
    this pass) and mean atmospheric_pressure ~611Pa (real Mars ~610Pa,
    +0.2%, was -11.6%) over a 1000-sample run -- the "as close as
    achievable without a zone change" caveat the old K tuning note
    carried is resolved by the zone change.
  - **N (Venus analog)**: `greenhouse_multiplier_range` cut roughly in
    half (370-420 -> 260-295) and `atm_density_range` raised (270-320 ->
    300-350) -- the old range was deliberately oversized specifically to
    compensate for N's too-cold midpoint placement, so keeping it at the
    same size now overshoots real Venus's temperature once N is
    correctly placed near the zone's hot inner edge. Mean
    surface_temperature ~737K (real Venus 737K, +0.0%) and mean
    atmospheric_pressure ~9.17MPa (real Venus ~9.2MPa, -0.3%, was -16.8%)
    over a 1000-sample run.
  - Every other affected class (E/F/G/H/L/O/P/V) keeps its existing
    albedo/atmosphere/greenhouse ranges unchanged -- none of them chase a
    single real-world numeric target the way K/N do, and each one's
    existing relative ordering (hotter/colder than its neighbors in the
    E->F->G->M/O progression, K/L, P) still holds with the new
    distance-aware placement. Their "Verified via climate_tuning_cli.py"
    comments are refreshed with new measured means reflecting the new
    placement.

## [5.20.0] - 2026-09-11

### Changed
- **`update.sh` skips the package reinstall when there's nothing new to
  install.** Previously it unconditionally re-ran the whole of
  `install.sh` (`pip install --force-reinstall`, an NLTK re-fetch,
  re-enabling Apache modules, a full permissions pass) every single
  invocation, even when `git pull` found no new commits at all -- pure
  wasted work for a caller like `examples/maintenance/planetgen-update.timer`
  that may run this monthly for years between real updates. Now, when the
  pull is a no-op, `update.sh` runs `src/migrateDb.py` directly instead
  (a cheap, idempotent no-op once the schema is already current) and
  skips the rest; `install.sh` (migration included, as its own step 2)
  still runs in full whenever the pull actually brought new commits, same
  as before.

### Added
- **`examples/maintenance/`: `update.sh` runs on the same schedule as the
  orbit update.** New `planetgen-update.service`/`.timer` run
  `sudo ./update.sh` (git pull + `install.sh`) monthly, at a fixed time
  30 minutes ahead of `planetgen-orbits@.timer`'s own now-fixed time (both
  timers dropped `RandomizedDelaySec` in favor of this deliberate,
  guaranteed ordering) -- update.sh can `pip install --force-reinstall` a
  new version of the very `stellarObjects` code `updateOrbits.py` imports,
  so the code update needs to land first, not run independently sometime
  in the same month. `planetgen-orbits@.service` also gained an
  `After=planetgen-update.service` ordering line for the case where both
  happen to be queued together. `install-maintenance-timer.sh` installs
  and enables both by default, sharing the same `/etc/planetgen/
  maintenance.env` credentials file (`update.sh`'s `install.sh` step needs
  DB credentials for `src/migrateDb.py` too); pass `--skip-update-timer`
  to opt out of unattended code updates and keep only the orbit timer, if
  this deployment's branch should only ever be updated by a human running
  `update.sh` deliberately.

## [5.18.0] - 2026-09-11

### Added
- **`examples/maintenance/`: systemd timer for `updateOrbits.py`.** An
  Ubuntu/Debian-native alternative to the raw crontab line
  `docs/database-schema.md` already documented for running the periodic
  orbital-motion update ("once a month or so"). `planetgen-orbits@.service`/
  `.timer` are a systemd *template* unit -- the instance name (e.g.
  `planetgen-orbits@planetgen.timer`) selects which database gets
  updated, so a deployment with more than one `PLANETGEN_MYSQL_DATABASE_PREFIX`
  schema enables one timer instance per database rather than needing a
  separate script per database. `install-maintenance-timer.sh` installs
  both units, writes `/etc/planetgen/maintenance.env` (mode 600) from
  `maintenance.env.example` for the shared read-write MySQL credentials
  (skipped if that file already exists, so it never clobbers credentials
  already set up), and enables the timer for each database name given on
  its command line (defaulting to `$PLANETGEN_MYSQL_DATABASE`, or
  "planetgen"). Output is captured by journald automatically, so there's
  no logfile/logrotate entry to maintain the way the crontab example
  needs.

## [5.17.0] - 2026-09-11

### Added
- **Binary mutual-orbit position.** `updateOrbits.py`/`_db.advance_orbital_phases`
  now recomputes `star_systems.binary_mutual_position_x_km`/`_y_km`/`_z_km`
  (the secondary star's Cartesian position relative to the primary) every
  time it advances `binary_mutual_orbital_phase_deg` -- the same "position
  has no independent update of its own, it just has to move whenever
  phase does" treatment `position_x/y/z_km` already gets for planets/
  moons. Derived via `utils.orbital_position_au` from
  `binary_separation_km` and the mutual orbit's own (non-near-ecliptic,
  full `[0, 180)`/`[0, 360)`-range) inclination/ascending-node/phase --
  the pair's true "direction of orbit" was already fully captured by
  those v13 orbital elements; this just keeps the derived position
  correctly in sync with them as time passes, rather than only at
  generation. Schema v14, with a migration that backfills real position
  values for existing binary rows from their own already-stored
  separation/orbital-element columns.

## [5.16.0] - 2026-09-11

### Added
- **Star motion: galactic orbit phase, plus binary mutual orbit.** Stars
  now get the same floating-point update guard planets/moons already
  have, and binary pairs get a real orbit around each other, both
  actively advanced over time.
  - Every star gains `galactic_orbital_phase_deg` (its current angular
    position around the galactic center, randomly rolled at generation --
    the same role `orbital_phase_deg` plays for a planet/moon) and
    `galactic_min_update_interval_years` (the same guard formula,
    `utils.minimum_update_interval_years` applied to
    `galactic_orbital_period_gy * 1e9` years). `_db.advance_orbital_phases`
    now advances this phase too, guarded by its own interval, reversing
    v12's "stars have no periodic update mechanism" scoping note -- they
    do now. Both stars of a binary pair, and `star_systems`' own mirrored
    `binary_galactic_orbital_phase_deg`/`binary_galactic_min_update_interval_years`,
    always carry the identical value: a binary's AU-scale separation is
    negligible next to its light-year-scale galactic orbit, so the pair
    moves around the galaxy together, not independently
    (`StarSystem.__init__` rolls the phase once and threads it to both
    stars and the proxy).
  - Binary pairs also get their own **mutual orbit** -- the two stars
    circling their common barycenter, entirely separate from (and vastly
    faster than) the galactic orbit above:
    `star_systems.binary_mutual_orbital_period_years`/`_speed_kms`
    (Kepler's third law / circular-orbit speed --
    `planetPhysics.calculate_orbital_period_years`/
    `utils.circular_orbital_speed_kms`, the same formulas a planet's own
    orbit already uses, applied to the pair's separation/combined mass),
    `_inclination_deg`/`_ascending_node_deg`/`_phase_deg` (the same
    `utils.orbital_position_au` orbital-element convention planets/moons
    use for "direction", but drawn from the full `[0, 180)`/`[0, 360)`
    range -- a binary's mutual orbital plane has no protoplanetary-disk
    reason to prefer any alignment), and its own
    `binary_mutual_min_update_interval_years` guard. `advance_orbital_phases`
    advances `binary_mutual_orbital_phase_deg` the same way.
  - Schema v13, with a migration that backfills every real-derivable
    value (both `*_min_update_interval_years` guards, plus the mutual
    orbit's period/speed) from existing rows' own already-stored data;
    the phase/orientation columns with no derivable "correct" value get
    an arbitrary `0` placeholder, the same treatment v9 already gives
    pre-existing planets'/moons' orbital orientation.

## [5.15.0] - 2026-09-10

### Added
- **Floating-point update guard.** Every planet and moon now stores
  `min_update_interval_years` -- the shortest `elapsed_years` worth
  calling `_db.advance_orbital_phases` for. Below this threshold, the
  phase delta `elapsed_years` would add is smaller than
  `orbital_phase_deg`'s own IEEE 754 double-precision resolution, so
  `MOD(orbital_phase_deg + delta, 360)` is guaranteed to round right back
  to the exact value already stored -- a wasted write that changes
  nothing. Derived purely from `period_years`
  (`utils.minimum_update_interval_years`: `period_years *
  math.ulp(360.0) / 360`, the coarsest representable step anywhere in
  `orbital_phase_deg`'s `[0, 360)` range). Not a narrative/display stat --
  purely a guard value: `advance_orbital_phases` now skips a row's
  `UPDATE` entirely (not just a no-op write, no attempt at all) when a
  call's `elapsed_years` is below it. In practice this floor sits many
  orders of magnitude below any realistic elapsed time
  (`updateOrbits.py` runs "once a month or so"), so it exists for
  correctness against a caller advancing time in much smaller steps, not
  because today's actual usage pattern comes close to tripping it.
  Scoped to `planets`/`moons` only -- `stars`' galactic-orbit values are
  fixed forever at generation time, with no periodic update mechanism to
  guard. Persisted as `planets`/`moons.min_update_interval_years` --
  schema v12, with a migration that backfills real derived values (not a
  placeholder) for existing rows.

## [5.14.0] - 2026-09-10

### Added
- **Planet/moon position.** Every generated planet and moon now gets an
  actual 3D Cartesian position (`Planet.position_x/y/z`, in AU), relative
  to its orbital anchor -- the star (or, for a binary system, the
  `BinaryStarProxy` standing in for the system's combined center) for a
  planet, the parent planet for a moon -- continuing the same "each body
  positioned relative to its immediate primary" hierarchy
  `docs/design/galaxy-coordinate-system.md` already uses one level up for
  sectors/systems relative to the galactic center. Derived from the
  existing orbital-motion elements (`orbital_inclination_deg`/
  `orbital_ascending_node_deg`/`orbital_phase_deg`, schema v9) via the
  standard circular-orbit-to-Cartesian transform
  (`utils.orbital_position_au`). Also adds `Planet.orbital_speed_kms`, a
  circular orbit's constant tangential speed (`utils.
  circular_orbital_speed_kms`, `v = 2*pi*r/T`) -- exact here, unlike the
  galactic orbit's rotation-curve model, since a planet's/moon's period is
  already known exactly from Kepler's third law. Surfaced as a new
  "Speed" row in every planet's/moon's rendered data table. Persisted as
  `planets`/`moons`.`position_x_km`/`_y_km`/`_z_km`/`orbital_speed_kms` --
  schema v11, with a migration that backfills real derived values (not a
  placeholder) for existing rows.
- `updateOrbits.py`/`_db.advance_orbital_phases` now recomputes position
  in lockstep with `orbital_phase_deg` as real time passes (a single
  set-based SQL `UPDATE` per table, matching phase advancement's own
  performance characteristics); `StarSystem.validate_system` recomputes
  position/speed too whenever it corrects a planet's `distance` post-hoc
  to resolve an orbital overlap, closing the same kind of staleness gap
  its own docstring already flagged for `period` before this.

## [5.13.0] - 2026-09-10

### Added
- **Galactic orbit.** Every generated star system now gets a circular
  orbital speed and period around the galactic center
  (`Star.galactic_orbital_speed_kms`/`galactic_orbital_period_gy`,
  `utils.calculate_galactic_orbit`), derived from the system's actual
  distance from the galactic center where known (a sector-placed system)
  or the same fixed `physical_constants.GALACTIC_CENTER_DISTANCE_LY`
  fallback `system_perimeter`/`heliosphere_radius` already use otherwise.
  Uses a simple rotation-curve model
  (`GALACTIC_ROTATION_FLAT_VELOCITY_KMS`/`GALACTIC_ROTATION_CORE_RADIUS_PC`)
  rather than a Keplerian point-mass orbit around the galaxy's total mass,
  which would overshoot Sol's real ~220-240 km/s orbital speed by roughly
  4x -- calibrated against Sol's own distance (~206 km/s, ~236 million
  years per orbit, both close to the real Sun's measured values). Surfaced
  as a new "Galactic Orbit" row in every star's/binary pair's rendered data
  table, and persisted as `stars.galactic_orbital_speed_kms`/
  `galactic_orbital_period_gy` (plus the `star_systems.binary_galactic_orbital_*`
  equivalents for binaries) -- schema v10, see `docs/database-schema.md`.

## [5.12.0] - 2026-09-10

_Originally developed and released as `5.9.0` on a separate branch; renumbered
on merge since `5.9.0` was independently already used below for the same-day
Class W removal. No functional difference from the original release._

### Changed
- **`../src/html/` is now a thin frontend for the Flask API instead of a
  direct MySQL client.** Every CGI page (`index.py`, `browse.py`,
  `sector.py`, `system.py`, `nav.py`, `search.py`, `galaxy.py`) fetches
  its data from `GET /api/...` via a new stdlib-only HTTP client
  (`html/lib/apiclient.py`, `PLANETGEN_API_BASE_URL`) rather than opening
  its own read-only MySQL connection -- the database-querying logic those
  pages used to duplicate now lives once in `queryDb.py`, shared with the
  API. `html/lib/dbutil.py` is gone; its non-DB formatting helpers
  (`esc`/`linkify_location`/`format_density`) moved to a new
  `html/lib/fmt.py`. The API itself gained what this needed: `?db=` on
  every route (multi-schema, matching the browser's own picker --
  `stellarObjects._db.list_databases`/`resolve_database`, moved out of
  `html/lib/dbutil.py` into the package itself), `GET /api/databases`,
  `GET /api/galaxy/sectors`, `GET /api/search` (the full faceted-search
  query layer, ported from `html/search.py`), `sector_id=none` on
  `GET /api/systems` (standalone systems), and richer `GET /api/sectors`/
  `GET /api/sectors/<id>`/`GET /api/systems/<id>` responses (display-
  ready fields -- ids, quadrant/location, star summaries, galaxy
  placement -- distinct from `stellarObjects._db.load_sector`/
  `load_star_system`'s *generation* object graph, still reachable the
  same way). Apache's example vhost
  (`examples/apache/planetgen.conf.example`) mounts the API at `/api/`
  (`WSGIScriptAlias`, in its own `WSGIDaemonProcess`) on the same vhost
  that serves `../src/html/` -- landed alongside `html/api/`'s own move
  into the same `DocumentRoot` (`[5.10.0]`), resolving the open question
  `TODO.md` left from the 5.3.2/5.3.3 file-system cleanup from two
  independent, orthogonal angles at once (where the API's source lives,
  and how `html/` gets its data) that happen to combine cleanly. See
  `docs/api.md` and `docs/html-interface.md`.
- System Map: a body with `life_chemical` set (habitable) now gets a
  small green badge on its marker, surfaced in the info panel's "Life
  Chemistry" field too (`html/lib/systemmap.py`'s `has_life`,
  `static/systemmap.js`) -- the one piece of state the existing
  circle+letter marker couldn't show at a glance. `docs/TODO.md`'s long-
  stale "sprite-based graphical system view" item is resolved: the System
  Map (an interactive scaled-orbit SVG diagram, not bitmap sprite art)
  already satisfied it, per that module's own docstring.

## [5.11.1] - 2026-09-10

### Changed
- **Moon tidal locking is now physics-based, not a flat probability.**
  [5.11.0] gave a moon a flat 75% (`MOON_TIDAL_LOCK_PROBABILITY`) chance
  of being tidally locked; replaced with an actual tidal-despinning
  timescale estimate (new `planetPhysics._tidal_locking_timescale_seconds`,
  the standard simplified Murray & Dermott formula: `t_lock = (2Q/15k2) *
  omega0 * a^6 * m_moon / (G * M_primary^2 * R_moon^3)`, with fixed
  representative `Q`/`k2` values for a rocky/icy body -- new
  `physical_constants.MOON_TIDAL_DISSIPATION_Q`/`_LOVE_NUMBER_K2`).
  Verified directly against three real systems spanning many orders of
  magnitude: ~47 million years for the real Earth-Moon system (real
  estimates: tens of millions of years), ~90 years for Mars/Deimos
  (consistent with Deimos being locked given its tiny size), and ~1
  billion years for Saturn/Iapetus (consistent with real estimates of a
  billion-year-plus despinning time for that unusually slow case). A
  candidate (pre-locking) rotation period is drawn first, same as any
  planet; a moon actually ends up locked only if that timescale is
  shorter than the system's age (`star.age`, the best available proxy --
  planets/moons don't carry an independent age of their own).
- **Moon orbital distance is now drawn log-uniformly, not
  linearly-uniformly.** Surfaced by the physics-based tidal locking
  above: `generate_moons`' distance range can span many orders of
  magnitude (its outer bound reaches 1/5 of the parent planet's own Hill
  radius -- tens to hundreds of millions of km for a large planet, far
  beyond where any real large moon actually orbits, e.g. our Moon at
  384,400 km), and a plain `random.uniform` over that range spends almost
  all its density in the single largest order of magnitude -- nearly
  every generated moon landed implausibly far out, so almost none had
  time to tidally lock (measured: ~1.3% of moons locked). Real moon
  systems are much closer to log-spaced (e.g. the Galilean moons run
  421,700 / 671,100 / 1,070,400 / 1,882,700 km, each roughly 1.5-1.6x the
  last), which log-uniform sampling matches far better while still
  allowing occasional genuinely distant/irregular moons. Measured effect:
  ~1.3% -> ~14% of moons locked over the same sample size -- still a
  minority overall (this generator's moons span a much wider population
  than just the handful of large, close, well-known real moons that
  dominate popular intuition about "most moons are locked"), but an order
  of magnitude more of them landing close enough to plausibly have
  locked.

## [5.11.0] - 2026-09-10

### Added
- **Orbital and rotational motion.** Every planet/moon now has a real 3D
  orbital orientation and a rotation period, plus a live position that a
  new, separately-run script advances over time -- resolves
  `docs/TODO.md`'s "Introducing realistic orbital paths and speeds..."
  entry.
  - `Planet` gains four new attributes (`planetPhysics.
    generate_orbital_motion_properties`): `orbital_inclination_deg`/
    `orbital_ascending_node_deg` (fixed at generation time -- together
    they orient this generator's circular-orbit model in 3D; planets draw
    from a tighter, real-solar-system-like range than moons, which can be
    tilted much further), `orbital_phase_deg` (this body's current
    position angle around its orbit -- the one field that changes after
    generation), and `rotation_period_hours` (axial "day length" -- a
    static descriptive stat; no rotational phase is tracked, by design,
    since nothing consumes "which side currently faces the primary"). A
    moon has a real chance (75%, `MOON_TIDAL_LOCK_PROBABILITY`) of being
    tidally locked (rotation period equal to orbital period) -- the norm
    for real large moons, not a rare special case; otherwise rotation
    period is drawn from a `body_type`-appropriate range (`physical_constants.
    ROTATION_PERIOD_RANGE_HOURS`).
  - **New `src/updateOrbits.py`.** Advances every planet's/moon's
    `orbital_phase_deg` in the configured database based on real elapsed
    time since the last run (`stellarObjects._db.advance_orbital_phases`,
    a single set-based `UPDATE` per table, not a per-row Python loop) --
    meant to be run periodically (e.g. cron, "once a month or so"), not
    on every generation run. A new `orbit_simulation_state` singleton row
    tracks when it last ran, measured server-side via `TIMESTAMPDIFF`
    rather than trusting the calling process' own clock to agree with the
    database server's (`get_orbit_update_elapsed_years`). The first run
    against a database just establishes that reference point (nothing to
    advance yet) rather than guessing a start time.
  - **Orbital period fix, needed for this feature to mean anything.**
    `Planet.period` used to be `sqrt(distance_au^3)` unconditionally --
    correct Kepler's-third-law shorthand only for a 1-solar-mass primary.
    For an ordinary planet orbiting a star of very different mass, and
    *especially* for a moon (whose real primary is its parent planet, not
    the grandparent star `self.star` still points at for other purposes
    like life chemistry) this was wrong by orders of magnitude -- a moon's
    period came out as if it orbited its host star directly at that same
    tiny distance. New `planetPhysics.calculate_orbital_period_years(distance_au,
    primary_mass_kg)` takes the actual primary's mass
    (`Planet.__init__`'s new `primary_mass_kg` parameter, threaded in
    from `generate_moons` as the parent planet's own `mass` for a moon,
    defaulting to `star.mass` for an ordinary planet -- transparently
    correct for a binary system too via `BinaryStarProxy.mass`'s existing
    effective-mass property). Also fixed a related latent staleness bug
    found while testing this: `StarSystem.validate_system` adjusts a
    planet's `distance` after generation to resolve orbital overlap and
    already recomputed atmospheric conditions to match, but never
    recomputed `period` -- now it does.
  - Schema v8 -> v9: `orbital_inclination_deg`/`orbital_ascending_node_deg`/
    `orbital_phase_deg`/`rotation_period_hours` on `planets`/`moons`, and
    the new `orbit_simulation_state` table. `stellarObjects._db.
    _migrate_v8_to_v9` is the first real per-version migration step of the
    MySQL era (every database before this one started fresh, already at
    the then-current schema) -- existing rows default to `0`/`24` (an
    arbitrary but harmless placeholder), since this generator never ran
    its actual random generation for them; every body generated from this
    point on gets real values. See `docs/database-schema.md`.

## [5.10.1] - 2026-09-10

### Fixed
- **Subdwarf (Yerkes VI) age modeling could produce a pre-Big-Bang star.**
  `Star._calculate_initial_star_age_and_lifespan` routed every non-main-
  sequence, non-white-dwarf Yerkes class (giants, subgiants, bright
  giants, supergiants, hypergiants, *and* subdwarfs) through the same
  model: derive main-sequence lifespan from the star's own already-
  generated mass, then draw age from the post-main-sequence window. That
  model's mass-sampling side already got a reject-and-resample guard
  against pre-Big-Bang progenitor masses ([5.3.1]) for every class it
  applies to -- except Yerkes VI, deliberately excluded because its
  *entire* allowed mass range (0.1-0.8 Msun) implies a main-sequence
  lifespan longer than the age of the universe, so every draw would have
  been rejected. That left subdwarfs still running through the age *model*
  itself, unguarded, which routinely produced ages of hundreds of billions
  of years -- because real subdwarfs (sdB/sdO) are thought to form via
  binary mass-stripping near the tip of a lower/intermediate-mass
  progenitor's red-giant branch, not single-star post-main-sequence
  evolution, so this star's own post-strip mass was never a valid stand-in
  for a progenitor's main-sequence lifespan in the first place -- the
  model was wrong for this class, not just missing a guard rail (flagged
  as an open design question in `docs/TODO.md`'s "Future ideas" ever
  since).
  - Yerkes VI now has its own age-generation branch, entirely independent
    of the star's own mass: age is drawn directly from a dedicated,
    old-population-biased range (new `program_constants.SUBDWARF_MIN_AGE_GY`
    /`SUBDWARF_MAX_AGE_GY`, 1.0-13.5 Gy, with the same young/old
    `SystemConfig.AGE` bias every other branch applies), reflecting that
    real subdwarf progenitors are typically old, low-mass population stars
    (a short-lived, higher-mass star wouldn't have had time to reach the
    RGB tip and get stripped). Lifespan is that age plus a short
    remaining-phase window (new `SUBDWARF_REMAINING_PHASE_MIN_GY`/
    `_MAX_GY`, 0.05-0.3 Gy) rather than a separately-derived value, since
    the real core-helium-burning subdwarf phase is short relative to the
    age itself.
  - New regression test `test_subdwarf_age_never_exceeds_universe_age`
    (`test_star_matrix.py`) locks this in across every spectral letter.

## [5.10.0] - 2026-09-10

### Changed
- **Moved the Flask API (`src/api/`) into `src/html/`.** Resolves the
  open question `docs/TODO.md` had carried since the 5.3.2/5.3.3
  file-system cleanup: the API now lives at `src/html/api/`, served from
  the same checkout/deployment tree as the interim CGI browser instead of
  a second, separately-tracked location. `src/wsgi.py` moved alongside it
  to `src/html/wsgi.py` (same "lives next to `api/`, so `import api` just
  works via `sys.path[0]`" property as before, plus an explicit
  `sys.path` entry for `src/` now that `queryDb`/`stellarObjects` are a
  directory farther away). `pytest.ini`'s `pythonpath` gained `src/html`
  so `test_api.py`'s `from api...` imports keep resolving unchanged.
  `examples/apache/planetgen.conf.example` gained a `WSGIScriptAlias /api`
  pointing at `html/wsgi.py`, plus deny-all `<Directory>`/`<Files>` blocks
  for `html/api/` and `html/wsgi.py` itself (mirroring the existing
  `html/lib/` block -- both hold source that's imported, never meant to
  be requested directly). Along the way, corrected a pre-existing
  inaccuracy in `docs/api.md`'s Apache deployment guidance: the vhost's
  `SetEnv` directives configure the CGI browser (`mod_cgi`/`mod_cgid`
  copies them into each script's real process environment) but never
  reach `os.environ` under `mod_wsgi` -- `PLANETGEN_MYSQL_*` for the API
  needs to come from the Apache service's own process environment instead
  (e.g. `/etc/apache2/envvars`).

## [5.9.1] - 2026-09-10

### Changed
- **Class P climate tuning.** Gave Class P the same
  `atm_molar_density_range`/`atm_density_range`/`greenhouse_multiplier_range`
  treatment the M/O/H/K/L/N/E/F/G/V pass ([5.3.7]) gave the other habitable
  classes; P had stopped at just its own `albedo_range` that pass. No
  single real-world analog for P, so this isn't chasing a target delta the
  way M/K/N are -- instead the new ranges make the class's own "cold,
  glaciated"/"thinning with age" flavor text physically real: molar
  density stays near Earth's real value (P's atmosphere text names
  oxygen/nitrogen/argon, not a heavier CO2-like mix), `atm_density` is set
  thin, and `greenhouse_multiplier` is set weak so the cold comes from
  genuine physics on top of the class's already-tuned high albedo, not
  albedo alone. Verified via `climate_tuning_cli.py --class P`: mean
  surface_temperature ~219K (well below freezing, clearly colder than
  Class M's ~286K), mean atmospheric_pressure ~10.4kPa (~0.1 atm) over a
  400-sample run across the full host-star grid.

## [5.9.0] - 2026-09-10

### Removed
- **Class W ("a tidally locked world with extreme temperature
  variations") removed entirely.** Its day/night-split identity can't be
  produced from a single global `surface_temperature` scalar under this
  generator's climate model -- no per-class range (albedo, molar density,
  greenhouse multiplier, or atmosphere density) reaches it without an
  actual dayside/nightside model this generator doesn't have, flagged as
  a known gap during the [5.3.7] climate-tuning pass and left open in
  `docs/TODO.md`'s "Investigate Further" section ever since. Rather than
  build a whole day/night thermal model for one class, cut it entirely --
  removed from `PLANET_CLASSES`, `PLANET_CLASS_PROBABILITIES` (its
  0.0001 weight just dropped; these are relative weights, not a
  normalized distribution, so nothing needed redistributing),
  `HABITABLE_PLANET_CLASSES`, and `MOON_BLACKLIST`
  (`program_constants.py`), plus its color entry in
  `html/lib/systemmap.py`'s `_CLASS_COLORS`.

## [5.8.3] - 2026-09-10

### Changed
- **Provisional sector designation is now a single packed hex number.**
  `provisional_sector_designation` (added in [5.8.2]) dropped its
  `"R<ring>-Q<quadrant>-<slot>"` letter/dash format for a plain
  bit-packed hex integer -- `ring` in the high bits, `quadrant - 1` in
  the next 2, the raw `shell_slot_index` in the low 32
  (`DESIGNATION_SLOT_BITS`/`DESIGNATION_QUADRANT_BITS`), e.g.
  `"1500002EE0"` (was `"R5-Q2-2EE0"`). Genuinely reversible back to
  `(ring, quadrant, shell_slot_index)` now too -- the fixed-width bit
  fields have no ambiguous boundary the way concatenating separately-
  sized hex numbers would.

## [5.8.2] - 2026-09-10

### Added
- **Provisional sector designations for un-generated/unvisited addresses.**
  `stellarObjects.galaxyGeometry.provisional_sector_designation(shell_index,
  shell_slot_index, edge_pc, edge_ly)` builds a short, human-readable label
  for a `(shell_index, shell_slot_index)` sector address --
  `R<ring>-Q<quadrant>-<slot>`, Ring and slot index in uppercase hex, Quadrant
  a plain 1-4 digit (e.g. `"R5-Q2-2EE0"`) -- the way a real astronomical
  catalog gives a not-yet-fully-characterized object a provisional name from
  its position rather than waiting for one. `sector_ring`/`sector_quadrant`
  (new, duplicating `html/lib/galaxymap.py`'s identically-named Ring/Quadrant
  concept in the dependency-free `galaxyGeometry` module so `galaxyGen.py`'s
  CLI doesn't need to import the web front-end layer) back it, and are both
  O(1) -- no scan over a shell's other slots, consistent with this codebase's
  galaxy-skeleton design principle of never doing per-sector work
  proportional to a shell's slot count. `galaxyGen.py`'s `--shell`/
  `--center-sector` generation modes now print this designation alongside
  each newly saved sector's procedural name and raw shell/slot numbers.

## [5.8.1] - 2026-09-10

### Added
- **NAV Map: a rendered plot of the NAV feature's origin/destination/route.**
  `queryDb.nav_between` now also returns `origin_position`/
  `destination_position` (the light-year positions its `direct` course was
  computed from) and, when a route was found, each hop's own position
  (`route["positions"]`) -- surfaced from `GET /api/nav`'s JSON response
  too (see `docs/api.md`'s updated "NAV" section). `src/html/lib/navmap.py`
  builds a new "NAV Map" panel on `html/nav.py`'s results page from that
  data: a flat, static, top-down SVG plot of the galactic X-Y plane --
  origin and destination as labeled, clickable points, a dashed line for
  the direct course, and a solid polyline through the optimal route's
  intermediate hops when one exists. Modeled on `galaxymap.py`'s flat 2D
  SVG rather than `starmap.py`'s rotatable 3D CSS scene: like a galaxy
  Quadrant, it's deliberately blind to altitude (the course panel's own
  Altitude figure already covers that axis), auto-scaled to whatever
  points it's given (no fixed sector size to normalize against) with one
  uniform light-years-per-pixel ratio on both axes so azimuth angles
  aren't visually distorted, plus a `+X` compass tick tying the plot's
  orientation to the course panel's own azimuth convention and a
  light-year scale-bar legend. Closes the rendered-image gap the original
  NAV feature (CHANGELOG [5.8.0]) left open in `docs/TODO.md` and
  `src/api/routes.py`'s `systems_near` TODO comment.

## [5.8.0] - 2026-09-10

### Added
- **NAV feature: course, distance, and optimal routing between two
  systems.** Three new pure/query modules plus an API endpoint and a web
  page:
  - `stellarObjects/navigation.py` -- `course_between` (Euclidean
    distance plus galactic-plane-relative azimuth/altitude: azimuth from
    +X in the galactic X-Y plane, altitude as elevation above/below that
    plane) and `warp_travel_times` (`velocity_multiple_of_c = warp_factor
    ** (10/3)`, reported at warp 1/3/6/9, formatted via the existing
    `utils.years_to_time_string`).
  - `stellarObjects/navGraph.py` -- `build_knn_adjacency` (a symmetrized
    k-nearest-neighbor adjacency graph over a `{id: (x, y, z)}` position
    set) and `shortest_path` (Dijkstra) for the "optimal route via
    adjacent systems" half of NAV.
  - `queryDb.nav_between` -- resolves NAV availability between two
    systems (unavailable if either has no sector; same-sector always
    available; cross-sector only when both sectors have a galaxy
    placement), combining a sector's galaxy-frame center (parsecs) with a
    system's sector-local offset (milliparsecs) into one absolute
    position (no such combinator existed before this), then returns the
    direct course plus the optimal route.
  - `GET /api/nav?from=<id>&to=<id>` -- JSON endpoint over `nav_between`,
    `404` for an unknown system id, `400` (`NavUnavailable`) when NAV
    doesn't apply to the pair. See `docs/api.md`'s new "NAV" section.
  - `src/html/nav.py` -- a destination picker (a `<select>` of the
    origin's own sector-mates, plus a typed destination-id field when
    cross-sector NAV is available) and a results page (direct course,
    warp travel times, and the hop-by-hop optimal route, each hop linking
    to `system.py`). `system.py` links here ("Navigate from here")
    whenever a system has a sector.

## [5.7.0] - 2026-09-09

### Changed
- **Merged the parallel Phase 4 "galaxy density skeleton" work (schema
  v6-v8: `sector_vertices`, `galaxy_shape`, `galaxy_shell_band`,
  `galaxyGen.ensure_sector_generated` lazy generation -- see [5.4.7] and
  [5.4.8] below) into the MySQL backend from [5.5.0].** That work landed
  on `main` entirely against the pre-port SQLite backend while this
  branch's MySQL port was in flight from the same base commit, so every
  new table, `_db.py` function (`get_sector_id_at`, `save_galaxy_shape`,
  `get_galaxy_shape`, `replace_galaxy_shell_bands`,
  `get_galaxy_shell_bands`, the vertex-writing/-reading in
  `insert_sector`/`get_sector_galaxy_position`), and CLI script
  (`galaxyGen.py`'s rewrite, the new `galaxyPlan.py`) needed porting to
  MySQL/`pymysql`/`MySQLConfig` as part of reconciling the two lines of
  work, same conventions as the rest of the MySQL port: `?` placeholders
  (translated by the `Connection` wrapper), `--mysql-*`/`config=` instead
  of `--db-path`/`db_path=`, and MySQL's `INSERT ... ON DUPLICATE KEY
  UPDATE` instead of SQLite's `INSERT ... ON CONFLICT DO UPDATE` for
  `save_galaxy_shape`'s singleton-row upsert.

### Fixed
- **`stellarObjects/galaxyGeometry.py` called an undefined
  `_theta_for_index` in both `sector_position_pc` and
  `sector_wedge_vertices_pc`** -- a pre-existing bug on `main` (confirmed
  present there independent of this merge), only `_phi_for_index` had
  ever been defined despite the module's own `GOLDEN_RATIO` docstring
  describing the golden-angle azimuthal step it was supposed to compute.
  Every galaxy-placement test failed with `NameError` until this merge's
  full test run surfaced it. Added the missing function
  (`theta_i = (2*pi*i / GOLDEN_RATIO) mod 2*pi`), matching the exact
  formula `test_galaxy_geometry.py`'s own worked-example test already
  documented and asserted against.
- **`stellarObjects._db.Connection`'s `execute`/`executemany` could
  return a `tuple` instead of a `list` from `.fetchall()` when zero rows
  matched** -- confirmed by testing, `pymysql`'s own cursor returns `()`
  for no rows but a `list` when rows exist, unlike `sqlite3`'s cursor
  (always a `list`). Surfaced by a real MariaDB test run as a spurious
  `assert result == []` failure (`test_ensure_sector_generated_reports_no_content_outside_every_stored_band`)
  that never appeared without a live server. Added a `_Cursor` proxy
  wrapping every returned cursor to normalize `fetchall()` to always
  return a `list`, so no other call site (present or future) can trip
  over the same inconsistency.

## [5.6.0] - 2026-09-09

### Added
- **Rate limiting (Flask-Limiter), applied to the whole API out of the
  box.** Default limits (200/day, 50/hour per client IP -- Flask-Limiter's
  own quickstart example, overridable via `PLANETGEN_RATELIMIT_DEFAULT`)
  apply to every route except `/api/health`; every write endpoint (below)
  additionally layers a stricter 10/minute limit on top
  (`routes.WRITE_RATE_LIMIT`). Exceeding a limit returns `429
  {"error": "rate limit exceeded", "detail": "..."}` (never Flask-Limiter's
  own default plain-text body) with `Retry-After`/`X-RateLimit-*`
  headers. Storage backend defaults to in-memory
  (`PLANETGEN_RATELIMIT_STORAGE_URI`, correct for a single-process
  deployment; a multi-worker `mod_wsgi`/`gunicorn` deployment needs a
  shared backend, e.g. Redis, or each worker under-enforces the
  configured limit by tracking its own separate counters).
- **Write-stub endpoints**: `POST`/`PATCH`/`DELETE` on `/api/sectors` and
  `/api/systems`. Every one validates its JSON request body (sectors get
  a fully mapped-out `{"name", "edge_ly"}` schema; systems only check
  "is a JSON object" for now -- see docs/api.md for why that one's
  field-level schema is still an open design question) and applies
  `WRITE_RATE_LIMIT`, but always responds `501 {"error": "... is not
  implemented yet"}` -- no row is ever inserted, updated, or deleted.
  Exists now so the request/response contract is settled and testable
  before the real database logic lands; see docs/api.md's "Write
  endpoints" section for what filling them in for real will also need
  (a write-capable database account, and authentication/authorization --
  neither exists yet, and both are load-bearing on this API staying
  read-only in practice today).

## [5.5.0] - 2026-09-09

### Added
- **Flask API: pagination, input validation, health check, and JSON-only
  error handling.** `/api/sectors`/`/api/systems` now return a paginated
  `{"items", "total", "limit", "offset"}` envelope instead of a bare list
  (`limit` defaults to 100, clamped to 500; `offset` defaults to 0) --
  this project's own roadmap plans galaxy-scale generation, so an
  unbounded listing endpoint would eventually return an unbounded
  response. Every query parameter (`limit`/`offset`/`sector_id`/`radius`)
  is now validated and rejected with a `400 {"error": "..."}` rather than
  silently ignored or crashing. Added `GET /api/health` for
  liveness/readiness monitoring. Every error response -- 400, 404
  (including an unmatched route), 405, and 500 -- is now JSON, never
  Flask's default HTML error page, and a 500 never leaks exception detail
  to the client. See `docs/api.md`.
- **MySQL backend (`TODO.md` Phase 5), replacing SQLite entirely.**
  `stellarObjects/schema.sql` is now MySQL/InnoDB DDL (`BIGINT UNSIGNED`
  ids, `DOUBLE`/`VARCHAR`/`TEXT`/`LONGTEXT` typing, a real
  `schema_migrations` tracking table replacing `PRAGMA user_version`,
  every index/foreign key declared inline per table for
  `CREATE TABLE IF NOT EXISTS` idempotency). `stellarObjects/_db.py` now
  talks to MySQL via `pymysql` (pure-Python driver) through a small
  `Connection` wrapper that keeps every existing call site's
  `conn.execute(sql, params)`/`row["column"]` shape unchanged, backed by
  a real connection pool (`DBUtils.PooledDB`) per TODO.md's "add real
  connection pooling" note. Every tool that touches the database
  (`sectorGen.py`, `systemGen.py`, `galaxyGen.py`, `queryDb.py`,
  `migrateDb.py`, `src/api/`, `src/html/`) now takes `--mysql-*` flags/
  `PLANETGEN_MYSQL_*` environment variables (`stellarObjects._db.MySQLConfig`)
  instead of `--db-path`/`PLANETGEN_DB_PATH`; the CGI browser's database
  picker now lists MySQL schemas on the configured server (filtered by
  `PLANETGEN_MYSQL_DATABASE_PREFIX`) instead of `.db` files in a
  directory. A new one-time `src/migrateSqliteToMysql.py` script imports
  an existing pre-port SQLite database (already at schema v5) into
  MySQL. The SQLite-specific `v1`-`v5` in-place migration machinery
  (`migrate_database`'s per-version steps, gzip file backups,
  `BACKUP_MARKER`) is removed entirely, since every MySQL database this
  project creates now starts at the current schema directly. See
  `docs/database-schema.md`.

### Changed
- `setup.py`'s `install_requires` gained `pymysql`/`DBUtils` (core
  dependencies now, not just the `api` extra) -- every database-touching
  entry point needs them, not only the Flask API.

### Fixed
- **`test_mdconvert.py`'s own `sys.path` setup pointed at a
  nonexistent top-level `html/lib/` instead of `src/html/lib/`,
  silently masked by `test_db_migration.py` (collected first,
  alphabetically) inserting the correct path first.** Deleting
  `test_db_migration.py` (see above) exposed it; fixed the path
  computation directly.

## [5.4.8] - 2026-09-09

### Added
- **Galaxy-wide density "skeleton"** (schema v8: `galaxy_shape`,
  `galaxy_shell_band`; `stellarObjects/galaxySkeleton.py`; `galaxyPlan.py`):
  precomputes and persists where the galaxy has any content at all,
  without storing a single sector's position/density/vertices -- those are
  pure deterministic functions of `(shell_index, shell_slot_index)` plus a
  handful of galaxy-wide shape parameters, so they're recomputed on demand
  instead. `galaxy_shape` holds the galaxy's shape parameters as one
  singleton row; `galaxy_shell_band` holds one row per contiguous
  *candidate* slot-index band per shell (almost always exactly one, found
  via an exact upper bound over spiral-arm azimuth, bisected to precision)
  -- a safe superset, not a per-sector list. Reduces the galaxy's
  structural storage from the ~1 TB a naive per-sector plan table would
  need down to ~350 KB at real Milky-Way scale (~4,076 rows), built in
  parallel across shells (`multiprocessing.Pool`) in under a second. See
  `docs/design/galaxy-coordinate-system.md` section 10 for the full
  analysis and measurements.
- **`galaxyGen.ensure_sector_generated(shell_index, shell_slot_index)`**:
  the lazy, visit-triggered generation entry point built on the skeleton
  above -- returns an already-generated sector if one exists; otherwise
  checks the stored band and exact density to decide whether the address
  holds anything, and if so generates and persists it on the spot, using
  that position's own `relative_density` as the actual system-count
  multiplier (so a bulge sector and a sparse outer-disk sector generate
  proportionally different counts, not a uniform default). A new
  `sectors.UNIQUE (shell_index, shell_slot_index)` constraint (schema v8)
  turns a concurrent visit race into a recoverable `IntegrityError`
  instead of a duplicate row.
- **Galaxy disk/spiral density model implemented** (`stellarObjects/galaxyDensity.py`):
  exponential disk radial falloff x sech^2 vertical scale-height x
  logarithmic spiral-arm modulation, plus a spherical bulge, normalized so
  `relative_density == 1.0` at a calibration point. `build_galaxy_shape`
  auto-calibrates the normalization constant; `predicted_star_count` scales
  a caller-supplied baseline expected system count by relative density.
  Tested in `tests/test_galaxy_density.py`, and exercised end-to-end
  alongside sector geometry in a full (unfilled) small spiral galaxy
  simulation.

### Fixed
- **`galaxyDensity.relative_density` could raise `OverflowError` far off
  the galactic plane** relative to `disk_scale_height_pc`: its vertical
  falloff term computed `1.0 / math.cosh(x) ** 2` directly, which raises
  once `|x|` exceeds ~710 -- found by a new skeleton test using a
  toy-scale shape, reachable in production for any shape with a small
  scale height relative to its own radius. Fixed with an algebraically
  equivalent, overflow-safe rewrite (`_sech_squared`) that underflows to
  the correct `0.0` limit instead of crashing.
- **Two same-shell Voronoi tessellation bugs found by that end-to-end
  simulation, both specific to small/sparse shells** (`sectorGeometry.py`):
  (1) the same-shell candidate projection used the raw chord vector
  instead of a proper gnomonic (central) projection, understating real
  separation worse the farther away a candidate was; (2) the half-plane
  test applied after that projection used the flat-plane bisector formula,
  which only approximates the true spherical bisector in the small-angle
  limit and was too permissive on sparse shells, silently under-clipping
  cells. Both fixes are covered by new regression tests parametrized to
  include a small shell (shell 1), plus an exhaustive
  every-vertex-is-shared check for a fully-populated small shell. See
  `docs/design/galaxy-coordinate-system.md` section 9 for the full
  derivation. Re-running the simulation after both fixes: 0 unexplained
  gaps across 31,255 generated outer vertices (previously 4,798).

## [5.4.7] - 2026-09-09

### Added
- **Every galaxy-placed sector now has explicit vertices with genuinely
  zero gaps against its same-shell neighbors, stored as plain relational
  rows.** New `sector_vertices` table (schema v7) -- one row per vertex,
  no JSON or other serialized blob anywhere in the schema (matching this
  project's existing convention: variable-length structured lists always
  get their own child table, the same treatment `asteroid_belt_composition`/
  `planet_reflection_spectrum` already have). Each sector's vertices are
  built from an exact local spherical Voronoi tessellation among its
  same-shell neighbors (`stellarObjects/sectorGeometry.local_lateral_cell`)
  extruded radially between the shell's inner and outer bounding spheres
  (`prism_vertices`). Lateral sharing is exact, not approximate: each
  shared corner is the 3D circumcenter of a sector and two of its
  neighbors -- a plain geometric fact independent of which of the three
  computes it, so two real neighbors land on identical floating-point
  values (~1e-16 agreement, verified directly) rather than two nudged
  approximations. Vertex/face count varies per sector (typically 5-7, mean
  6.0) because it has to: a cube tiling of a sphere can't be gap-free in
  general (the same reason a soccer ball needs pentagons mixed with
  hexagons), so a fixed 8-vertex shape cannot exactly reconcile a sector
  with more real neighbors than it has faces -- confirmed by an earlier,
  never-released attempt at exactly that (fixed-cube corner averaging),
  which only ever reduced gaps (~35% aggregate), not eliminated them.
  Radially, adjacent shells match in total area covered, not
  vertex-for-vertex (a "non-conforming mesh interface," the same technique
  used where independently-meshed regions meet in finite-element/CFD
  meshing) -- shell k's outer bound and shell k+1's inner bound are the
  same sphere, and each shell's own sectors fully and independently tile
  it. Getting this fast required exploiting that the underlying placement
  is a Fibonacci lattice: true same-shell neighbors concentrate at
  Fibonacci-number index offsets (confirmed against real data: exact
  offsets of F_20 through F_23), which collapses same-shell neighbor
  search from ~0.4-0.5s/sector (the naive radius-search cost near a large
  outer shell's equator, exactly where the disk/spiral density model
  concentrates real generation) down to ~0.3-0.4ms/sector -- validated
  against a guaranteed-correct brute-force search across 840 cases
  spanning the full polar range and shell sizes from 3 to 211 million
  slots with zero mismatches, with small shells
  (`shell_sector_count(k) <= 2000`) falling back to unconditionally-correct
  brute force since the underlying asymptotic theory isn't reliable that
  close to the galactic core. `galaxyGen.py` computes and stores this for
  every sector it generates; `sectorGen.py`'s standalone (non-galaxy) CLI
  is unaffected, same as every other galaxy-frame data. New
  `_migrate_v6_to_v7` schema migration drops any v6 database's briefly-lived
  `vertices_pc` JSON column entirely rather than converting it into rows
  (a sector's vertices are cheap to recompute from its address if ever
  actually needed) -- `_migrate_v4_to_v5`/`_migrate_v5_to_v6` need no
  special-casing at all, since a v4 or v5 source's `sectors` table already
  has the same shape as the current one.

## [5.4.6] - 2026-09-07

### Fixed
- **System page table of contents was simply unavailable below a 90rem
  window width.** It lived exclusively in a fixed right-hand margin rail
  (`.toc` in `src/html/static/style.css`) that only existed at
  `min-width: 90rem`; narrower windows got `display: none` and no
  alternative. It's now a collapsible pulldown rendered inline above the
  description at any width, expanding on click, and only switches to the
  fixed sidebar (always shown open) once the window is wide enough for
  one. Built with a hidden checkbox + `<label>` rather than native
  `<details>`/`<summary>`: a closed `<details>`'s children turn out to
  stay out of layout/paint via internal browser state that isn't fully
  reachable through CSS overrides of `display`/`content-visibility`
  (confirmed by testing -- only toggling the element's `.open` property
  itself, not any stylesheet rule, restored it), which made "always open
  at the wide breakpoint" unreliable. A plain checkbox has no such
  internal state to fight.

## [5.4.5] - 2026-09-07

### Fixed
- **Sector Map: the 3D scene drew over surrounding page content instead
  of staying confined to its panel.** Zooming in (or rotating to an angle
  where the cube's diagonal grew past its nominal footprint) had nothing
  bounding where the map was visible, so it grew and drew over the rest
  of the page. Added `.starmap-viewport`, a fixed-size (`aspect-ratio: 1
  / 1`, matching the square footprint the old static image occupied),
  `overflow: hidden` window that the 3D scene now rotates/zooms/pans
  inside of -- clipped at its edges instead of spilling out. The 3D
  content inside renders/rotates/occludes exactly the same either way;
  this only bounds where it's visible from the outside.

## [5.4.4] - 2026-09-07

### Added
- **Sector Map is now interactive 3D: drag to rotate, scroll (or +/-
  buttons) to zoom.** Replaced the static server-baked isometric SVG
  projection in `src/html/lib/starmap.py` with a real CSS 3D scene
  (`transform-style: preserve-3d`, orthographic -- no `perspective`), so
  the browser's own compositor handles rotation and occlusion instead of
  a hand-rolled JS matrix routine; `src/html/static/sectormap.js` tracks
  two rotation angles and a zoom factor and feeds them straight to the
  scene's CSS transform. Each star dot is billboarded (counter-rotated
  every frame to keep facing the camera) so it stays a circle instead of
  going edge-on as the view turns, and dot-click detection is resolved by
  geometry (`getBoundingClientRect`) rather than native hit-testing,
  since the latter turns out unreliable for elements nested this deep in
  a rotated `preserve-3d` hierarchy.

### Changed
- **Sector Map star colors now come from the star's actual named
  spectral color** (`SPECTRAL_CLASS_COLORS`: Blue/Blue-White/White/
  Yellow-White/Yellow/Orange/Red, keyed off `star_type`'s leading letter)
  instead of a raw Kelvin-to-RGB blackbody approximation, so a "White
  Giant" reads white and a "Blue Giant" reads blue regardless of its
  exact temperature. Luminosity shades each color's vividness/lightness
  (brighter = more vivid, dimmer = more muted) and temperature nudges
  lightness within the star's own spectral band.
- **Binary systems draw two dots** (a larger primary and a smaller
  secondary, capped at 65% of the primary's radius and offset to its
  lower-right, overlapping) built from each component's own `stars` row,
  instead of one dot from the system-level `binary_type`/
  `binary_temperature_k` summary -- so each half of a binary is sized and
  colored from its own actual data.

See `docs/TODO.md` ("Near-term: interim `../src/html/` browser
enhancements") for where this started as a plan.

## [5.4.3] - 2026-09-07

### Changed
- **Web interface: system-page table of contents moved to a fixed
  right-margin rail.** The collapsed `<details>` TOC added in [5.4.2]
  still lived inline next to the description as a flex sibling, narrowing
  the prose column whenever it was open. `system.py`'s `_toc_html` now
  renders a plain, always-expanded `<nav>`, and `style.css` positions
  `.toc` fixed in the right margin (mirroring the left `.sidenav`),
  appearing only once the window is wide enough (`min-width: 90rem`) to
  hold it without crowding the main content -- narrower windows simply
  don't get one, rather than it floating over the page.

## [5.4.2] - 2026-09-07

### Removed
- **Class R ("an ejected, geologically active world") removed entirely.**
  It had `h`/`e`/`c` all `False` -- zero probability weight, unreachable
  outside a manual `zone_override`. A genuinely free-floating/rogue planet
  (no host star at all) is a real exoplanet category, but doesn't fit this
  generator's star-centric zone model -- "which zone" is the wrong
  question for an object with no star to be zoned relative to. Rather than
  leave it as permanent dead weight, cut from `PLANET_CLASSES`,
  `PLANET_CLASS_PROBABILITIES` (already 0.0000), and `MOON_BLACKLIST`
  (`program_constants.py`); the now-vacuous
  `test_known_issue_class_with_no_valid_zone_is_unreachable` regression
  test removed from `src/tests/test_planets.py`.

### Fixed
- **Naming: triple-consonant validation bug.** `utils.is_name_valid`
  enforced "no more than two consecutive vowels" but had no matching check
  for consonants, so a fixed denylist of specific clusters
  (`BAD_CONSONANTS`) was the only thing standing between a generated name
  and an arbitrary run of 3+ consonants -- empirically ~40% of generated
  star names had one. Added a `consonant_count` run-length counter
  symmetric to the existing `vowel_count` one. Verified 0/5000 across
  star/planet/moon/sector name generation post-fix.
- **Web interface: navbar had no way back to the database picker on a
  single-database deployment.** `lib/page.py`'s sidenav only linked to
  `index.py` when more than one `.db` file existed, since `index.py`
  itself auto-redirects past the picker for a single database -- but that
  meant the common single-database production deployment
  (`starmap.moltenaether.com`) had no Databases link at all. Fixed by
  always linking to `index.py?all=1`, and teaching `index.py` to render
  the full picker table (bypassing its own auto-redirect) when `?all=1`
  is present.
- **Web interface: system-page table of contents wasn't collapsible.**
  `system.py`'s `_toc_html` now builds a `<details>`/`<summary>` pair
  instead of a plain `<aside>`/`<h3>`, collapsed by default -- a native,
  no-JS toggle.
- Minor cleanups found in passing: an unused `program_constants` import in
  `sectorGen.py`, and a stray `f`-string prefix on a non-interpolated
  string in `evolution.py`.

## [5.4.1] - 2026-09-07

### Added
- **`sectorGen.py --density`: a controllable density value for sector
  generation.** Sectors previously always generated a flat count of systems
  (`--num-systems`, default 10) with no connection to sector volume -- there
  was no way to make one generated sector meaningfully denser or sparser
  than another. `--density` is a float multiplier on the real local stellar
  density this codebase already models (`physical_constants.LOCAL_STELLAR_DENSITY_LY3`
  via `SpaceSector.expected_system_count()`) -- 1.0 means a realistic sector
  this size, 2.0 twice as dense, 0.5 half. The actual per-sector count is
  drawn with the existing `_sample_poisson_count` helper (previously only
  used by `SpaceSector.grow_from_seed`), resolved fresh for each sector
  rather than once at parse time, so counts vary naturally sector to sector
  under `--num-sectors` and across `galaxyGen.py` runs, matching a real
  spatial Poisson process. Mutually exclusive with `--num-systems`; neither
  flag given keeps the original flat-10-systems default unchanged.
  `galaxyGen.py` gets the flag for free since it shares
  `sectorGen.add_shared_generation_options`/`validate_shared_generation_args`
  with no changes needed on its side.

## [5.4.0] - 2026-09-07

### Removed
- **Two brown-dwarf-scale gas-giant classes cut entirely.** Correcting
  their radius ranges to real sub-stellar physics in [5.3.9] left them as
  near-duplicates of each other (overlapping radius, differentiated only by
  an invisible density number no description text reflects), each carrying
  a vanishingly small (0.01%) generation weight, representing objects that
  aren't really planets in the first place (brown dwarfs are sub-stellar
  objects; this codebase already has a dedicated mechanism for stellar/
  sub-stellar companions via `BinaryStarProxy`/`BINARY_SYSTEM`). Their
  combined generation weight folded into the ordinary Jupiter/Saturn-class
  gas giant. All remaining references to the removed classes (and to the
  earlier-removed small-rocky-class pair) scrubbed from comments/docs,
  generically describing what they explain rather than naming
  now-nonexistent class letters.

### Added
- **Per-class `size_mode`: bell-curve size distributions for planets and
  moons.** Radius generation (`planetPhysics.generate_planet_properties`'s
  three radius-draw sites, and `generate_moons`' moon-radius draw) switched
  from a flat uniform draw across each class's declared radius range to a
  bounded Gaussian ("bell curve") draw peaking at a new per-class
  `size_mode` value (0.0-1.0, "what fraction through the available range is
  the statistically most common size") via new `utils.sample_bounded_bell`.
  For a moon, "available range" is its actual Hill-sphere/mass-capped
  window, not necessarily the class's full declared range -- `size_mode` is
  read relative to whatever range is actually being drawn from. Every
  surviving class's `size_mode` is set from a real-world single-body analog
  where one exists (Earth for M, Mars for K, Venus for N, Jupiter/Saturn
  for J, Uranus/Neptune for I/T, the real rocky-to-gaseous transition
  radius for V) or a reasoned default (small-body populations, both real
  asteroids/KBOs (Class C/D) and this generator's own hot-zone rocky
  classes (A/B), skew toward their smaller end, matching real
  size-frequency distributions). Verified empirically: generated radius
  means land within ~1% of each class's real-world anchor (e.g. Class M
  mean 6,337km vs Earth's 6,371km; Class K mean 3,410km vs Mars's 3,389.5km).
  The self-adjusting spread (`sample_bounded_bell`'s `spread_divisor`)
  keeps a legible bell shape even for a mode pinned near one edge of a
  range, via rejection sampling rather than clamping (which would pile
  spillover probability mass at the boundary).
- `PLANET_CLASS_PROBABILITIES`'s `J` weight increased 0.0529 -> 0.0531 to
  absorb the removed brown-dwarf-scale classes' combined weight.
- Atmospheric-pressure sanity bound (`test_full_matrix.py`/`test_planets.py`)
  lowered back down `1e9` -> `5e7` Pa now that the brown-dwarf-scale
  gravity/pressure extremes those classes produced are gone (observed max
  across the full star-type matrix is now ~1.1e7 Pa, from Class N).
- `test_planet_physics_fixes.py`'s gravity/pressure correlation test
  threshold recalibrated (Spearman `>0.5` -> `>0.1`) for the narrower
  gravity range the remaining gas giants span without the removed classes'
  two-orders-of-magnitude density spread; its density-blend-skip test
  rewritten to inject a temporary `density_range` onto an existing class
  via `monkeypatch` rather than depend on a specific class declaring one
  (no class currently does -- it's generic, reusable override
  infrastructure, same as every other per-class override this codebase has
  built up).

## [5.3.9] - 2026-09-07

### Removed
- **Classes X and Y removed, merged into B and A/B respectively.** Both
  were small, redundant variants of existing hot-zone rocky classes:
  - Class X ("a stripped core from a gas giant", no atmosphere, radius
    500-5000km) folded into Class B ("a small, molten world") as an
    alternate origin story ("occasionally the stripped core of a former
    gas giant") rather than a separate atmosphere-less class -- a
    Mercury-analog this close to its star already has only a negligible
    exosphere in reality, so B's existing thin atmosphere covers X's "no
    atmosphere" identity closely enough.
  - Class Y ("a 'demon' class world", radius 5000-7500km) folded into
    *both* Class A and Class B (both radius ceilings raised 5000->7500km
    to absorb Y's size range), its toxic/irradiated flavor folded into
    each class's description as a variant rather than kept as a third
    class. B's composition gained "and sulfur" for the shared
    volcanic/irradiated theme.
  - `MOON_BLACKLIST` and `PLANET_CLASS_PROBABILITIES` updated to drop X/Y
    (their combined generation weight folded into A/B rather than dropped).

### Changed
- **Gas-giant zones reworked using real exoplanet science.** Previously
  every gas/ice-giant class (I/J/S/T/U) was valid in zone `c` only --
  meaning no gas giant could ever appear close to its star or in the
  habitable zone, despite "hot Jupiters" and "warm Jupiters" being a
  standard, well-documented real classification (orbital period < 10 days
  / 10-365 days respectively) and Neptune-mass planets in or near a
  temperate zone being common too.
  - Class J now valid in `h`/`e`/`c` (hot/warm/cold Jupiter) -- a warm/cold
    Jupiter placed in zone `e` can also generate ordinary moons via the
    existing moon-generation path, including habitable-class ones (the
    "habitable exomoon around a giant planet" trope), verified working.
  - Class I now valid in `e`/`c` (warm/cold Neptune), deliberately **not**
    `h`: real close-in Neptune-mass planets are rare (the observed "hot
    Neptune desert") because stellar irradiation photoevaporates a
    Neptune-mass H/He envelope down to a bare rocky/metal core well before
    it could stay Class I -- that outcome is exactly Class B's newly-merged
    "stripped core" identity (see Removed, above).
  - Classes S/T/U left `c`-only: real directly-imaged super-Jovian/brown-dwarf
    companions are predominantly found at wide separations (formation and
    detection-bias reasons), and close-in high-mass companions, while known,
    are much rarer (the "brown dwarf desert").
  - Fixing `generate_moons` (planetPhysics.py) to respect
    `HABITABLE_WORLD=False` the same way direct planet generation already
    does -- unreachable before this change (no gas giant could ever be in
    zone `e`), but a warm/cold Jupiter placed there could otherwise roll a
    habitable-class moon even in a system explicitly configured to disallow
    habitable worlds.
- **Classes S ("supergiant") and U ("ultragiant") radius ranges corrected
  to real sub-stellar physics** -- previously 250,000-50,000,000km and
  25,000,000-60,000,000km respectively, i.e. up to ~86 solar radii, larger
  than most actual stars. Real brown dwarfs stay within ~15% of Jupiter's
  own radius (69,911km) across their *entire* 13-80 Jupiter-mass range
  (electron degeneracy pressure means more mass compresses them, R ~
  M^(-1/8)); even the smallest true hydrogen-fusing red dwarf stars are
  only ~0.1 solar radii (~69,600km, essentially Jupiter-sized). Corrected to
  S: 60,000-120,000km, U: 65,000-130,000km (deliberately overlapping --
  that overlap *is* the real physics, not an oversight). What actually
  differentiates a higher-mass sub-stellar object at essentially the same
  radius is **density**: real measured brown-dwarf densities run roughly
  10-200 g/cm^3 (~10-150x an ordinary gas giant's), so both classes gained
  a new `density_range` (S: 10.0-60.0, U: 60.0-150.0 g/cm^3). Class T ("gas
  dwarf") also corrected, 250,000-25,000,000km -> 15,000-55,000km -- its
  old floor matched S's own floor, letting a "dwarf" be exactly as large as
  a "supergiant"; now spans ice-giant-to-Saturn scale, meaningfully below
  both J and S.
- **New per-class `density_range` override** (`PLANET_CLASSES`, read by
  `planetPhysics.get_planet_mass_ranges`/`generate_planet_properties` and
  mirrored in `plausibility.theoretical_gravity_bounds_g`) -- same
  `.get(..., default)` pattern as `atm_molar_density_range` etc. Used by
  Classes S/U above.
- **Fixed a gas-giant density-blend bug that silently collapsed every gas
  giant's density to a near-zero, physically meaningless value**, previously
  flagged as necessary-but-deferred follow-up work in
  `test_gas_giant_sampled_densities_are_finite_positive_and_within_theoretical_bounds`'s
  own docstring. `generate_planet_properties`' core/envelope harmonic-mean
  density blend reused `planet.atm_density` (drawn from
  `ATMOSPHERE_DENSITY["g"]`) as the envelope term -- but that value is a
  thin, surface/pressure-layer density (a different physical layer, correct
  for the separate atmospheric-pressure/scale-height calculation), ~1000x
  lighter than a real gas-giant envelope's actual bulk density. A harmonic
  mean is dominated by whichever term is smaller almost regardless of mass
  fraction, so this collapsed density to ~0.001-0.003 g/cm^3 for every gas
  giant, independent of core density -- silently defeating the new S/U
  `density_range` work above (their much denser cores had almost no effect
  on the final blended density). Fixed by introducing
  `physical_constants.GAS_ENVELOPE_BULK_DENSITY` (0.06-0.3 g/cm^3, grounded
  in real measured "puffy" gas giants -- WASP-193b's ~0.06 g/cm^3 is the
  lowest confirmed bulk density known), used only for this blend. A class
  declaring its own `density_range` (S/U) now skips the blend entirely and
  uses that density directly -- real brown dwarfs don't have a meaningfully
  separate light envelope over a denser core the way an ordinary gas giant
  does, so blending toward a light "puffy" value would just dilute the
  elevated density right back down.
- Atmospheric-pressure sanity bound (`test_full_matrix.py`/`test_planets.py`)
  raised `5e7` -> `1e9` Pa -- Classes S/U's now-correct brown-dwarf-like
  gravity legitimately pushes pressure up to ~3.3e8 Pa across the full
  star-type matrix.
- `test_planet_physics_fixes.py`'s two gas-giant-blend tests updated to
  match the new formula (queuing an `envelope_density_gcm3` draw instead of
  reusing `atm_density`), plus a new
  `test_density_range_override_skips_the_blend` covering the S/U direct-
  density path.

## [5.3.8] - 2026-09-07

### Fixed
- **Habitable/life-bearing classes restricted to the ecosphere zone.** Class
  Q (`h`/`e`/`c` all `True`) and Class W (`h`/`e` `True`) both carry a
  `life_chemical`, but weren't restricted to zone `e` like every other
  life-bearing class already was — meaning a life-bearing Q or W world could
  be generated directly in the hot or cold zone. Both are now `e`-only
  (`h`/`c` `False`); Q's "eccentric orbit" flavor still holds fully confined
  to the ecosphere zone, and W's "tidally locked" identity is arguably more
  scientifically apt restricted to the ecosphere zone (real tidally-locked
  *habitable* worlds are a real, actively studied trope specifically because
  a cool star's habitable zone sits close enough in for tidal locking to be
  near-guaranteed — e.g. TRAPPIST-1's planets, Proxima b).
  New regression test `test_life_bearing_classes_are_ecosphere_only`
  (`src/tests/test_planets.py`) locks this in for every class with a
  `life_chemical`, not just the ones on the separately-curated
  `HABITABLE_PLANET_CLASSES` list.
- **Description/atmosphere-text cleanup across `PLANET_CLASSES`.** Several
  classes' `description` field repeated a word the render template
  (`planetData.to_paragraph_list`) already supplies via its own "with an
  atmosphere of {atmosphere}" or "with a composition of {composition}"
  clause, producing genuinely broken rendered sentences — e.g. Class B
  previously rendered "...a small, molten world **with a thin atmosphere
  with an atmosphere of** a mix of helium, sodium, and oxygen...". Fixed for
  Classes B, E, K, N, Q, W, X, and Y — descriptive detail that belonged on
  the atmosphere field (e.g. "thin", "dense, reducing") was moved there
  instead of dropped, and Y's atmosphere field (previously the only one not
  ending in a named gas mixture) reworded to "a turbulent, toxic, and
  irradiated mix of gases".

## [5.3.7] - 2026-09-07

### Changed
- **Greenhouse-formula fix and per-class climate tuning.** The prior
  `greenhouse_factor` formula (`planetPhysics.calculate_atmospheric_conditions`)
  used `atm_molar_density` as its only lever, scaled by a single shared
  `CO2_MAX_GREENHOUSE_FACTOR` cap (5) — real Earth air's own molar mass
  already produced `greenhouse_factor ≈ 3.33` under that formula, driving
  Class M's mean surface temperature to 362K instead of ~288K, and Mars and
  Venus (nearly identical real molar mass, ~100x different real greenhouse
  forcing) could never be told apart by molar mass alone.
  - `CO2_MAX_GREENHOUSE_FACTOR` (now 500) is a generous safety ceiling, not
    the calibration knob.
  - New per-class `PLANET_CLASSES` keys — `albedo_range` (extended beyond
    Class P, which introduced it), `atm_molar_density_range`,
    `atm_density_range`, and `greenhouse_multiplier_range` — give every
    tuned class independent control of composition (molar density),
    quantity (mass density), and potency (greenhouse multiplier), following
    the exact override pattern `albedo_range` established for Class P.
    Class N's old hardcoded `atm_density = 65` / `atm_molar_density = max`
    special case in `planetPhysics.generate_planet_properties` is folded
    into this same general mechanism.
  - **Tuned this pass** (via the new `climate_tuning_cli.py`, see below):
    M (Earth analog, ~286K/~99kPa), O (warm/wet ocean world, ~293K/~97kPa),
    H (hot/dry desert, ~325K/~42kPa), K (Mars analog, ~231K/~0.57kPa), L (K
    + vegetation, warmer/thicker than K, ~256K/~2.2kPa), N (Venus analog,
    ~740K/~9.4MPa), E/F/G (a young, cooling progression, ~373K -> ~329K ->
    ~292K, converging near M), and V (thick, hot Super-Earth, ~365K/~286kPa).
    Class P and W are unchanged this pass (P already had a working
    `albedo_range`; W's "extreme temperature variations" identity needs a
    day/night model this generator doesn't have, not just range tuning).
  - Text refinements: Class H's "and metals" -> "and mineral dust" (a real
    desert's atmosphere lofts particulate, not metal vapor); Class O's
    atmosphere text (previously byte-identical to Class M's) now mentions
    water vapor; Class E's vague "hydrogen compounds" now names a real
    Hadean/Archean-analog mix (water vapor, ammonia, methane); Class K's
    "carbon dioxide" -> "a thin mix of carbon dioxide and nitrogen" (names
    the real Mars-analog composition, not just class); Class V's atmosphere
    text now reflects its tuned CO2-retention (not primordial H/He) identity.
  - Class M's disabled gravity clamp (`planetPhysics.calculate_surface_gravity`)
    deleted outright (it had been commented out, inert, since the
    atmospheric-pressure formula fix; no longer worth keeping around "in
    case it needs restoring").
- **New: `src/tests/climate_tuning_cli.py`.** A human-driven tuning tool
  (mirrors `physical_plausibility_cli.py`'s batch-generate-and-report
  pattern, reusing `stellarObjects.plausibility`'s engine): generates N
  bodies of one class, reports summary stats, and — for classes with a
  direct real-world analog (M/Earth, K/Mars, N/Venus) — a delta line against
  that reference. `--albedo`/`--molar-density`/`--density`/`--greenhouse`
  flags temporarily monkeypatch that class's `PLANET_CLASSES` entry for the
  run only, so candidate values can be iterated without editing source
  between runs.
- **New: `src/tests/test_climate_tuning.py`.** Regression suite locking in
  the tuning above via generously-toleranced bands (Class M/N/K within
  Earth/Venus/Mars-like ranges) and relative orderings (N hottest/
  highest-pressure of the tuned classes; K colder/thinner than L; L colder
  than M; H hotter/drier than O; E > F > G cooling progression converging
  near M; M < V < N) rather than brittle exact-value assertions, since
  generation is inherently stochastic.
- `test_full_matrix.py`/`test_planets.py`'s atmospheric-pressure sanity
  bound raised from `1e7` to `5e7` Pa — Class N now legitimately reaches
  ~9-15MPa depending on host star luminosity (previously ~2.98MPa mean, per
  the greenhouse-formula bug above), and the old bound was sized for the
  un-tuned, incorrectly-cold N. `test_planet_physics_fixes.py`'s
  `test_class_p_has_own_albedo_range_distinct_from_default` updated to
  check P's range differs from M's own new tuned range, rather than
  asserting M has no override at all.

## [5.3.6] - 2026-09-07

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
    flagged as not yet done when this was paused mid-session. Developed on
    a branch cut before Track A's completion and the TODO/FIXME-comment
    migration (5.3.5); merged into `main` after both, with no functional
    changes needed beyond a `test_db_migration.py` assertion that had to
    decompress the (now gzip-compressed, per Track B) v3->v4 migration
    backup the same way the existing v1->v2 backup test already did.

## [5.3.5] - 2026-09-07

### Fixed
- **Atmospheric pressure is no longer independent of gravity**
  (`stellarObjects/planetPhysics.py`): the barometric-formula pressure
  calculation algebraically canceled gravity out entirely (`atmospheric_pressure
  = atm_density * R * T / atm_molar_density`), so a Neptune-gravity gas giant
  and a Jupiter-gravity one produced the same pressure. A new
  `_atmosphere_retention_factor(gravity_g)` (linear, normalized to 1.0 at
  Earth gravity) now scales an *effective* atmospheric density used only in
  the pressure calculation (not `planet.atm_density` itself, which also
  feeds the gas-giant density blend), reintroducing a real, tunable
  gravity/pressure relationship (Spearman correlation on a mixed
  terrestrial/gas-giant sample now > 0.5, vs. ~-0.11 before).
- **Class P ("cold, glaciated") is colder than Class M again**: gave Class P
  its own `albedo_range` (0.5-0.7, matching real ice/snow Bond albedo)
  instead of sharing the default rocky/Earth-like range (0.12-0.35) with
  every other terrestrial class. Previously P and M were statistically
  indistinguishable in temperature once the disabled clamp was removed (see
  `docs/analysis/habitability-atmosphere-sanity-review.md`); P's cold
  identity now comes from the unclamped physics instead of a post-hoc
  override.
- Completes Track A (see [5.3.4]'s gas-giant density/greenhouse fixes) --
  8682 tests passing, including 10 new regression tests in
  `src/tests/test_planet_physics_fixes.py`.

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
  (see `docs/webconfig.md`).
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
  See [`docs/webconfig.md`](docs/webconfig.md) for full documentation.

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