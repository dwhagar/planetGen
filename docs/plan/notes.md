# Research notes and dependency notes

Kept from the first phase plan (2026-10-01) and the dependency report
(2026-10-02). [docs/TODO.md](../TODO.md) is the master file and wins
where they disagree; the phase files list the items and their order.

## Research notes

### From "Phase 1: Foundations and fixes" (first plan)

Boss's research notes (kept in the project's shared files under `todo-tasks/research/`) proposed fixes and numbers. They were checked against the code on 2026-10-01; where they were wrong about the code, the correction is given. Their numbers are starting points to tune, not requirements.

- **Backfill level (GEN.44).** The notes give the semantics as a table:
  -1 never backfilled, 0 fully generated, a positive value the dimmest
  luminosity reached. A backfill at floor L skips a sector at 0 or at a
  level at or below L, and afterwards sets each processed sector to L
  (or 0 once it is fully generated). Correction: the notes put the
  column on `sectors`, but an unfilled sector has no `sectors` row, so
  the level needs its own table keyed by address (or a reworked
  `bright_star_blocks`, which already keeps the level per 3x3x3 block,
  schema v49). The notes' file list (`generationLimits.py`,
  `spaceSector.py`) is not where backfill lives; it is
  `backfill_bright_stars_around` and `_backfill_block` in `generate.py`
  and the block table in `_db.py`.
- **Rogue planet mix (GEN.45).** The notes propose a power-law mass
  function, dN/dM proportional to M^-0.65 from 0.01 Earth masses to 13
  Jupiter masses, giving mostly terrestrial and few giant rogues, and
  suggest an 80/20 terrestrial/gas split. Today's arithmetic is about
  91/9. Correction: there is no `phenomena_config` table; the weights
  are `ROGUE_PLANET_MASS_BINS` and the gas giant threshold in
  `program_constants.py`.
- **Two-word names (GEN.46).** The notes suggest resolving collisions
  inside two words, with catalog numbers joined by a hyphen
  (`Kepler-452`) rather than a third word, and binary suffixes folded
  into the two words. The words come from the uniqueness decorations
  (Greek prefix, "Alpha ... IV", stacked diminutives), so that is where
  the fix goes.
- **Nebulae (GEN.47).** The notes say nebula generation was cut out of
  the sector loop in `generate.py`. Correction: it isn't; nebulae are
  rolled per sector, but at rates that make them almost never appear in
  a 4 pc sector. The notes' proposal still fits the real fix: a
  multi-octave 3D noise density field with Poisson-seeded centers at
  galaxy scale, with these classes as a starting point:

  | Class | Made of | Field threshold | Radius | Near |
  |---|---|---|---|---|
  | Emission (H II) | ionized hydrogen | > 0.72 | 15 to 45 pc | O and B supergiants |
  | Reflection | silicate and carbon dust | > 0.65 | 10 to 30 pc | B and A main sequence |
  | Dark | cold molecular gas | > 0.80 | 5 to 25 pc | low-density voids |
  | Planetary | ejected envelope | event | 0.1 to 2 pc | post-AGB stars, white dwarfs |
  | Supernova remnant | blast shell | event | 5 to 20 pc | pulsars, neutron stars |

  The existing classes and rates are in
  [nebula-and-asteroid-field-classes.md](../design/nebula-and-asteroid-field-classes.md).

### From "Phase 2: Maps and navigation" (first plan)

Boss's research notes (kept in the project's shared files under `todo-tasks/research/`) proposed fixes and numbers. They were checked against the code on 2026-10-01; where they were wrong about the code, the correction is given. Their numbers are starting points to tune, not requirements.

- **The arc pick (MAP.85).** The notes lay it out in four stages: (1)
  the whole galaxy as a particle field of its stars and spiral arms,
  with no grid, block or sector lines; (2) hovering projects an arc
  through the full height of the disk, highlighted, with faint
  boundaries of the neighboring arcs; (3) clicking flies the camera to
  the arc, where the user picks a height band (a slab); (4) the slab's
  sectors appear for the next pick. Note the word: "arc" was MAP.19's
  name for one cell of the old 3x3 region pick
  ([galaxy-drilldown-navigation.md](../design/galaxy-drilldown-navigation.md));
  it now means this first pick; the design doc's section 15 ("Planned:
  the arc pick", PR #362) describes it, and is rewritten as the current
  design when MAP.85 ships.
- **Sector colors (MAP.86).** The notes' starting values:

  | State | Opacity | Saturation | Hue from |
  |---|---|---|---|
  | Unfilled | 0.03 | 0 | neutral dark grey |
  | Filled, empty | 0.15 | 0.10 | the site's accent |
  | Sparse (1 to 10 stars) | 0.20 | 0.35 | mostly M dwarfs, red-amber |
  | Dense (50+ stars) | 0.45 | 0.85 | G and A stars, yellow-white |
  | Densest | 0.70 | 1.00 | O and B stars, blue |

  Saturation scales with star density, hue is the luminosity-weighted
  average of the stars' colors by temperature, lightness scales with the
  log of total luminosity, and a block's color, saturation and opacity
  are the plain average of its sectors'. Boss asked for filled sectors
  to stay translucent, "just a hair more solid" than unfilled, so the
  top of that opacity range is a ceiling to try, not a target.
- **Rogue planets on the Sector Map (MAP.82 to MAP.84).** The notes'
  values:

  | Property | Unmarked (default) | Marked |
  |---|---|---|
  | Radius | 1.5 px | 5 px |
  | Opacity | 0.20 | 1.0 |
  | Color | muted grey (#4A5568) | bright with a glow ring (#00E5FF) |
  | Hit radius | 3 px | 12 px |
  | Button | plain | highlighted, `aria-pressed="true"` |

  The button exists today (`lib/starmap.py`, `sectormap.js`) but starts
  pressed and has no pressed style.

### From "Phase 3: Interface, API and work queue" (first plan)

Boss's research notes (kept in the project's shared files under `todo-tasks/research/`) proposed fixes and numbers. They were checked against the code on 2026-10-01; where they were wrong about the code, the correction is given. Their numbers are starting points to tune, not requirements.

- **Structured planet details (UX.30).** The notes sketch a planet card:
  a header with class, sector and habitability; columns for physical
  properties, an atmosphere profile with a bar per gas, and satellites;
  and an actions row. They suggest HTML custom elements (Web
  Components). The site renders on the server with Jinja today, so
  server-built HTML with the same layout fits better unless a part
  needs to update in the browser.
- **Action menus (UX.26, UX.31).** One "Actions" button per row opening
  a floating menu that doesn't move the page, with entries such as Edit,
  Regenerate, Backfill sector and Delete; editing opens a focused dialog
  with Save and Cancel.
- **One-row header (UX.27).** Details and actions on one row; when it
  is too narrow, the navigate buttons fold into one "Navigation" menu
  with "From here" and "To here". The notes say below 768 px; the site's
  UX rules use container queries and the size classes in TODO.md's UX
  section, so use those.
- **Icons (UX.28).** Starting list: edit (pencil), delete (trash), show
  on map (crosshair), filter (funnel), navigate (route or compass),
  actions menu (ellipsis or gear), each with `title` and `aria-label`.
- **Rows (UX.32).** One "Class" chip (for example "Class M") and a moon
  count badge that opens the moons, instead of the type chip and moon
  label.
- **Sector table (UX.24).** Order: star systems, then non-stellar
  phenomena (nebulae, belts, comets, anomalies), then rogue planets; a
  detail row is `<tr><td colspan=...>` across every column.
- **Rogue planet octants (UX.25).** The notes number the octants 1 to 8
  by the signs of (dx, dy, dz) from the sector center: 1 (+,+,+),
  2 (-,+,+), 3 (-,-,+), 4 (+,-,+), 5 (+,+,-), 6 (-,+,-), 7 (-,-,-),
  8 (+,-,-). Check it against the octant the site already shows before
  changing anything.
- **Phenomena filters (UX.33).** Filter by kind, then by class: rogue
  planets by class (and terrestrial or gas giant, mass range, octant),
  comets by orbit type (short period, long period, hyperbolic or
  interstellar, sungrazer), and the other kinds by class. Correction:
  the notes' single `phenomena` table with `category` and `class`
  columns doesn't exist; each kind has its own table
  (`rogue_planets.planet_class`, and so on), so the filter is per kind in
  `queryDb`.
- **Comet links (UX.29).** Every comet resolves to a class and links to
  its reference page (period, eccentricity, perihelion). Parabolic
  comets have no period class today.

### From "Phase 4: Accounts, the view from a planet, more galaxies" (first plan)

None of Boss's 2026-10-01 notes apply to this phase.

## Files that several items touch

| Files | Items | Order |
|---|---|---|
| generate.py: parallel scatter and progress | PERF.21, PERF.23, GEN.32 (with TEST.74's tests) | One thread (phase 0, A). |
| stellarObjects/workQueue.py | PERF.21, PERF.22, TEST.73, then ADM.15, PERF.18, PERF.19/20 | PERF.21's thread first; the others start only after it merges. |
| generate.py: bright-star backfill (backfill_bright_stars_around, _backfill_block) | GEN.44, GEN.41, GEN.42, GEN.43, PERF.18 | Fixed order GEN.44, GEN.41, then GEN.42 + GEN.43 + PERF.18 in one thread. |
| generate.py: command-line options | GEN.51, GEN.52, GEN.24, API.3 | GEN.51/52 before GEN.24's new mode; API.3's remote mode mirrors the final options. |
| generate.py: sector summary | UX.34, OPS.9 | One PR. |
| systemData.py StarSystem constructor | GEN.49, GEN.50, then GEN.52 | One thread in phase 0 (GEN.53 and GEN.54 done, PR #367); GEN.52 after it. |
| planetPhysics.py (reconcile_zone_and_class, generate_moons) and PLANET_CLASSES | GEN.33/28, GEN.27, GEN.38, GEN.60, GEN.29 | Physics bugs done (PR #350); the classes thread. |
| Random draws in every generator file | GEN.39, then GEN.56 (decided yes, Boss 01:34Z) | Touches almost every file above; land it right after PERF.21 and tell the other generation threads to merge main. |
| _db.py | API.10 (id blocks), GEN.46 then API.12 (names) | DB.2 to DB.5 done (PR #342, PR #347). |
| Galaxy schema (schema.sql, v50 today) | DB.6, DB.7, GEN.44, PERF.11, MAP.86 (if it adds a column), NAV.10, API.11, GEN.46 (if names are migrated) | One writer at a time, in this order: DB.6 (phase 0, fresh galaxy), then DB.7, GEN.44, PERF.11 with MAP.86, NAV.10, API.11. DB.8 only reads it. |
| Control schema (v7 today) | OPS.13 (key history), API.9, API.15 (call log), USR.2, USR.4, USR.7, NAV.19 | One writer at a time; OPS.13 and API.9 first, accounts later. |
| lib/systemmap.py and static/systemmap.js | MAP.57, MAP.88, MAP.89, MAP.71 | One thread: MAP.57, MAP.88, MAP.89 (the System Map lane, after the Galaxy map picker lane starts); it can use mapcore.js helpers (MAP.63, PR #351). |
| sectormap.js and lib/starmap.py | MAP.65, MAP.79, NAV.29, MAP.68 | Phase 0 fixes and the extraction done (PR #351); later items in the engine thread. |
| galaxymap3d.js, galaxystageview.js, galaxystages.js, galaxyblocks.js | MAP.60, MAP.55, MAP.85, MAP.52, MAP.86, MAP.56, MAP.53, MAP.58, MAP.78, MAP.54, MAP.76, MAP.75, MAP.59, MAP.77, NAV.31 | One ordered Galaxy Map thread: the Galaxy map picker and arc lane (MAP.60, MAP.55, MAP.85, MAP.52) in phase 0, then phases 1 and 2. |
| Galaxy tiles (queryDb tile listing, lib/galaxymap3d.py, galaxyViewport.py, tile cache) | MAP.80, MAP.86 | In that order (MAP.90 done, PR #365); any payload change bumps the tile cache. |
| galaxyGeometry.py and galaxyprisms.js | GEN.24, MAP.85 | GEN.31 (PR #353) and NAV.38 (PR #357, sectors_along_segment / sectorsAlongSegment) done. |
| bookmarks.js | MAP.55, NAV.18, USR.7, NAV.19 | MAP.81 done (PR #351): plain 1 to 9 keys. |
| static/mapcore.js (shared helpers) and static/mapcontrol.js (camera and input controller; zoom policies free, range and locked, MAP.58 uses ZOOM_LOCKED) | MAP.60, MAP.55, MAP.85, MAP.52, MAP.53, MAP.58, MAP.75, MAP.65 to MAP.68, MAP.71 | New in PR #351 (MAP.63, MAP.64); later map items build on them rather than copying helpers. |
| Generate page (generate.html) | ADM.14, ADM.16, GEN.24 | ADM.14 first. |
| System page (system.html, lib/systempage.py, system_pages.py) | UX.29, NAV.8, UX.27, UX.31, UX.32, UX.30, MAP.74 | Roughly in that order; UX.32 and UX.30 in one thread. |
| Sector page (sector_page.py, sector.html, edit_controls.html) | UX.24, UX.25, UX.26 | One thread, in that order. |
| Navigation (nav_page.py, queryDb.nav_between, navigation.py) | NAV.7, NAV.10, NAV.11, NAV.16, NAV.17 | NAV.7 first; NAV.10 and NAV.16 touch different functions. |
| Admin edits (editStore.py, adminEdits.py, api/edits.py) and the galaxy settings JSON | ADM.18, GEN.59, GEN.61, OPS.18, ADM.19, DB.10, OPS.12 | ADM.18 writes the file; GEN.59 records deltas in the control database; GEN.61 merges them daily; OPS.18 keeps 18 backups. |
| update.sh and update.ps1 (with install.* and deploy-common.*) | OPS.7, OPS.8, then OPS.13, then OPS.15, then OPS.17 | OPS.7 and OPS.8 one PR; OPS.13, OPS.15 and OPS.17 after it, in that order. |
| test_gen_bright_scatter_edges.py | TEST.76, TEST.74, PERF.21 | TEST.76 first. |

## Near-cycles and how they are broken

- **GEN.39 and PERF.21**: Each said it touches the other's tests ("whichever lands second keeps both working"). Broken by putting PERF.21 first: its tests get fault injection that works across processes instead of leaning on the seeded single-process stream, then GEN.39 makes per-sector draws reproducible at any worker count.
- **MAP.88 and MAP.89**: Both said "whichever lands second keeps the other". Fixed order: MAP.57, MAP.88, then MAP.89, one thread.
- **MAP.85 and MAP.86**: MAP.85 removes the lines and relies on color to show structure; MAP.86 is that color. Boss moved MAP.85 to phase 0 (01:46Z); MAP.86 follows next in phase 1, with MAP.86's data side (stored per-sector color, tile field) able to land first.
- **MAP.61 and the Galaxy Map items**: MAP.61 is both before and after MAP.52 to MAP.60. Split as its text says: MAP.63 and MAP.64 in phase 0, MAP.65 to MAP.68 after the selection rewrite.
- **NAV.3 and MAP.61**: The picker needs the engine and the engine's panel carries the picker's buttons. Split: NAV.13/NAV.14 (no engine) in phase 1, pick mode (NAV.15) on MAP.65's panel layer in phase 2, NAV.32 at the end.
- **API.6 and USR.2**: Settled by Boss: keys belong to accounts, so API.6 moved to phase 3+ after USR.2.
- **NAV.4 and USR.1**: Already resolved by the per-browser default: NAV.4 in phase 2, NAV.19 moves courses into accounts in phase 3+.

## Judgment calls

- **Galaxy Map drill-down bugs not in phase 0**: MAP.52, MAP.53, MAP.54, MAP.56, MAP.77 and MAP.78 are bugs, but MAP.85's arc pick replaces the code they live in, so fixing them first would be thrown away. Their root is MAP.63/MAP.64 (phase 0) and MAP.85 (phase 1).
- **MAP.80 moved up to phase 1**: The old plan put it at the end of the selection chain. The thinning lives in the tile listing and the star points, not the pick code, so it can go right after MAP.90 and MAP.87.
- **Sector Map rogue fixes and NAV.30 in phase 0** (done, PR #351): MAP.82, MAP.83, MAP.84 and NAV.30 are small fixes in today's sectormap.js and galaxymap3d.js. Fixing them before MAP.63 moves the code, with TEST.70 pinning them, keeps them fixed through the engine work. Alternative: build them on the engine in phase 2.
- **MAP.79's toggles in phase 2**: The dimming part is phase 0 (MAP.82 to MAP.84); the per-kind show/hide buttons go on the shared control set (MAP.65).
- **UX.25, UX.26, UX.27, UX.31 in phase 1**: They are bugs but use UX.28's icons, which Boss approves first. They could ship in phase 0 with text buttons and get icons later.
- **GEN.45 in phase 0**: The old plan had it after GEN.28. It only changes the rogue mass bins, and GEN.38's class choice follows from the masses, so it sits at the root. Done in PR #350, reading M^-0.65 per log mass (per unit mass gave 87% gas giants).
- **GEN.32 in the parallel thread**: It shares the scatter functions with PERF.23. It could instead wait for GEN.44 and use per-sector levels to skip finished work.
- **UX.21 last**: A bug, but a final pass over finished pages. Its one known dead control (the nebula "-" button at the 1 ly limit) could be split out into phase 0.
- **OPS.8 in phase 0**: Not a bug, but the same two files as OPS.7, so it rides in the same PR.
- **Accounts, view and GEN.9 in phase 3+**: Nothing blocks them; they are late by priority. VIEW.4 and GEN.9 float and can start any time. If a second galaxy is likely, GEN.9's plan before phase 2's schema work would say whether those migrations add a galaxy id.
