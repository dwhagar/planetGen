# Research notes and dependency notes

Kept from the first phase plan (2026-10-01) and the dependency report
(2026-10-02). [docs/TODO.md](../TODO.md) is the master file and wins
where they disagree; the phase files list the items and their order.

## Phase 0 (finished)

Phase 0 held every bug and the groundwork the new architecture needed
(package layout, third-party libraries, Alembic, per-unit seeds, the
map engine, Redis queues). Its last items shipped on 2026-10-09 (DB.13,
GEN.126; later GEN.125, GEN.56, GEN.100, GEN.98 moved in and finished),
Boss asked on 2026-10-09 08:37Z for its plan file to be deleted, and the
version moved to 8.0 (`changes/phase-0-complete.major.md`). The item
list and each item's PR are in git history (`docs/plan/phase-0-roots.md`
before this change) and in CHANGELOG.md. The "1 → 0" and "→ 0" rows in
the change list below record moves into it.

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
  it now means this first pick; the design doc's section 15 ("The arc
  pick") describes it as built. Shipped in PR #369 (MAP.60,
  MAP.55, MAP.85, MAP.52): arcs are 45 degrees wide (the width nearest
  40 degrees whose edges fall on wedge lines in every block ring) by a
  third of the disk radius, 24 arcs in all; the galaxy map is a tilted
  3D view with no grid lines, a one-line scale, and the buttons Back,
  Forward, Up, Reset, Bookmarks and Menu, as text labels with a hook
  for UX.28's icons.
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

### Sector stats defaults (PERF.11, done in PR #425)

PERF.11's open questions shipped with these defaults, which Boss can
still change: the `sector_stats` table (galaxy schema v53, shared with
GEN.44's backfill level) stores both systems and stars as "actual";
the decaying average of expected against actual is one galaxy-wide
figure over a window of 1000 fills; regenerating a sector adds a new
sample; deleting a sector leaves the average as it is.

## The 2026-10-07 rebuild: change list

### New items

From Boss's list of 2026-10-03: ADM.23, ADM.24, ADM.25, ADM.27, ADM.31, GEN.67, GEN.72, GEN.73, GEN.75, GEN.76, GEN.77, GEN.78, GEN.79, GEN.83, GEN.90, GEN.91, GEN.93, GEN.96, GEN.97, GEN.98, GEN.99, GEN.100, GEN.106, GEN.107, GEN.109, GEN.110, GEN.111, GEN.112, MAP.103, MAP.104, MAP.105, MAP.106, MAP.107, MAP.108, MAP.109, MAP.110, MAP.111, MAP.112, MAP.113, MAP.115, MAP.116, MAP.117, MAP.119, MAP.120, MAP.121, NAV.46, NAV.47, OPS.20, OPS.22, PERF.26, POP.7, POP.10, UX.42, VIEW.5, with their subitems.

From Boss's list of 2026-10-07: ADM.26, ADM.28, ADM.29, ADM.30, ADM.32, ADM.33, ADM.34, ADM.35, ADM.36, API.18, API.19, DB.13, GEN.68, GEN.74, GEN.80, GEN.81, GEN.82, GEN.101, GEN.102, GEN.103, GEN.104, GEN.108, GEN.113, GEN.114, MAP.118, MAP.122, MAP.123, MAP.124, MAP.125, MAP.126, NAV.48, NAV.49, PERF.28, PERF.29, PERF.30, SEC.31, UX.43, UX.44, UX.46, UX.47, UX.48.

From PR #434's CI (2026-10-07): TEST.88. From PR #442's run (2026-10-07): TEST.89. From Boss's message of 2026-10-07 16:26Z (three Galaxy Map problems): MAP.115 folded into MAP.116, and MAP.127 (filed then) folded into MAP.121. From Boss's install error of 2026-10-07 16:57Z: OPS.26. From his decision answers of 2026-10-07 17:11Z: OPS.27. From PR #476 (2026-10-07): GEN.116, TEST.90. From Boss's bright-star sweep of 2026-10-07 18:43Z: GEN.117. From the bulge check and Boss's choice of 2026-10-07 19:04Z: GEN.118, GEN.119. From PR #481's browser-a11y run (2026-10-07): TEST.91. From PR #496's local run (2026-10-07): TEST.92. From Foundations' full run after PRs #501 and #504 (2026-10-07): TEST.93. From Boss's prevalence decision of 2026-10-08 00:12Z: ADM.37. From Bugfixes lane 1's full run (2026-10-08): TEST.94. From Bugfixes lane 1's full run for PR #522 (2026-10-08): TEST.95, TEST.96. From Bugfixes lane 1's full run for PR #525 (2026-10-08): TEST.97. From Boss's message of 2026-10-08 01:59Z (how sectors and blocks are colored): DB.14, MAP.128, MAP.129, MAP.130, MAP.131, MAP.132. From Bugfixes lane 1's full run for PR #531 (2026-10-08): TEST.98, TEST.99, TEST.100, TEST.101. From Boss's message of 2026-10-08 02:59Z (put the phoneme codec into the project): GEN.120. From the coordinator (2026-10-08, after UX.40 closed): UX.49. From Foundations lane 1 (2026-10-08, NAV.15's body endpoints need NAV.16): NAV.50. From Foundations lane 1's full runs for PR #588 (2026-10-08): TEST.102. From Boss's approval of the UX.37 audit (2026-10-08 17:43Z): UX.50 to UX.74 (R1 UX.50, R2 UX.51, R3 UX.52, R4 UX.53, R5 UX.54, R6 UX.55, R7 UX.56, M1 UX.57, M2 UX.58, M3 UX.59, M4 UX.60, M5 UX.61, M6 UX.62, P1 UX.63, P2 UX.64, P3 UX.65, P4 UX.66, P5 UX.67, P6 UX.68, L1 UX.69, N1 UX.70, N2 UX.71, A1 UX.72, A2 UX.73, A3 UX.74). From the coordinator's report of the Python 3.9 CI failures (2026-10-08): ADM.38. From Foundations lane 1's PR #612 report (2026-10-08): TEST.103, TEST.104 (the rate-limit flake is already TEST.83).

From Boss's message of 2026-10-07 12:25Z (the galaxy's own gravity, a gap in the orbital documents): GEN.115.

From PR #431's red CI (coordinator, 2026-10-03): PERF.27, DB.12, OPS.25, MAP.114; the MySQL 8.4 deadlock is TEST.81 (noted there).

### Merged into one item (duplicates)

- "Star generation should always actually take place" (10-07) went into GEN.76 and GEN.78 (10-03).
- "bright star placement is exclusively on the central galactic plane between layers -121 and 121" (10-07) went into GEN.79 with the 10-03 bulge bug, now marked major.
- "Generating any space gives you the option to see that space in the galaxy viewer" (10-07) and the 10-03 Generate page button are ADM.31.
- "Button in Galaxy display to allow selecting a star" and the "Select Mode" (10-07) are MAP.122.
- "Plotted courses should appear on the galactic map and stay until cleared" (10-07) is in NAV.49 (NAV.5 notes it).
- The two "add a star system to a sector" asks (10-07) are ADM.32.
- "Have phenomena ... generated through the entire galaxy first" (10-07) and the 10-03 scatter move are GEN.100.
- "Research different ways to generate an object ID" (10-07) is GEN.68, the first step of the names-from-IDs stream.
- MAP.89 (fitted System Map orbit spacing with a 12 px gap floor) was closed on Boss's word (2026-10-09 07:13Z: "Retire MAP.89, we have that done already"). No PR carries it: on main the 2D System Map still places orbits on the log scale (`_radial_px`, web/maps/systemmap.py) and the 3D view's compressed mode is also logarithmic (systemscale.js), so reopen it if Boss sees rings crowded together.
- GitHub issues of 2026-10-08/09: #736 (Generate screen tabs) went into ADM.28 and #677 (search near a picked point) into NAV.45; #709 and #708 are MAP.136; #758 and #715 are MAP.139; #718 and #716 are MAP.141; #614 and #535 are ADM.41; #638 and #614 are PERF.34; #661 and #750 are PERF.32 and PERF.33; #761 is PERF.31; #777 and #778 are GEN.129, GEN.130 and their design GEN.128; #515 and #743 are ADM.43, ADM.44 and their model ADM.42.
- "We need a planet that is like earth sized but never had any life" (10-03) is GEN.28's class Z (Lifeless temperate world), noted on GEN.28.

### Folded into a larger stream (closes with it)

| Item | Folded into | Why |
|---|---|---|
| TEST.72 | SEC.29 | The TOTP code it tests moves to pyotp. |
| TEST.83 | SEC.30 | The rate limits move to Flask-Limiter on Redis. |
| PERF.19 | PERF.24 | Boss chose Redis and RQ; the audit is its first step. |
| OPS.19 | PERF.24 | The job store moves with the queue. |
| UX.2, ADM.14, UX.26, UX.31, UX.27, UX.49 | UX.40 | Menus and buttons became Shoelace components (UX.40 done, PR #544); UX.26's sector Admin panel and UX.49's form fields remain. |
| ADM.24, ADM.25, ADM.26 | ADM.22 | Logs and progress moved to SSE and Xterm.js (ADM.22 done, PR #551); these three are now done (PR #560). |
| MAP.108, MAP.107, MAP.112, NAV.46 | MAP.65, NAV.15 | The shared picking layer and pick mode. |
| GEN.33, GEN.28, GEN.27, GEN.29 | GEN.90 | The class refactor around the habitability index. |

### Dropped

- GEN.63 (planet names unique within a sector): names come from unique IDs (GEN.67), so clashes can't happen.

### Moved

| Item | Move | Why |
|---|---|---|
| NAV.7 | 1 → 0 | Moved into phase 0: the engine (MAP.67, NAV.13) needs it, and the engine fixes the map bugs. |
| GEN.33 | 1 → 2 | Moved to phase 2 under GEN.90, after the habitability score. |
| GEN.28 | 1 → 2 | Under GEN.90; class Z is Boss's Earth-size world that never had life (2026-10-03). |
| GEN.27 | 1 → 2 | Under GEN.90. |
| GEN.52 | 1 → 0 | Moved into phase 0 with GEN.48 (all bugs in phase 0). |
| TEST.75 | 1 → 0 | Moved into phase 0 with GEN.48. |
| ADM.16 | 1 → 0 | Moved into phase 0 with GEN.48. |
| GEN.48 | 1 → 0 | Moved into phase 0: Boss wants every bug in phase 0, and its controls come with it. |
| MAP.95 | 1 → 0 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. |
| NAV.13 | 1 → 0 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. |
| NAV.14 | 1 → 0 | Its breadcrumb fixes MAP.106. |
| PERF.19 | 1 → 0 | Moved into phase 0 as the first step of PERF.24 (Boss chose Redis, 2026-10-03). |
| MAP.65 | 2 → 0 | Moved into phase 0: fixes MAP.108, MAP.107, MAP.112 and NAV.46. |
| MAP.79 | 2 → 0 | Moved into phase 0 (bug); its nebula toggle closes MAP.113. |
| NAV.15 | 2 → 0 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. |
| NAV.29 | 2 → 0 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. |
| NAV.33 | 2 → 0 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. |
| NAV.16 | 2 → 1 | Moved from phase 2: NAV.3 closes in phase 1 now. |
| UX.37 | 2 → 0 | Moved into phase 0: UX.21 (a bug) needs it. |
| MAP.66 | 3 → 0 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. |
| MAP.67 | 3 → 0 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. |
| MAP.68 | 3 → 0 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. |
| NAV.32 | 3 → 0 | Moved into phase 0 with the engine (all bugs in phase 0). |
| NAV.3 | 3 → 1 | Moved from phase 3: its parts are in phases 0 and 1. |
| UX.21 | 3 → 0 | Moved into phase 0 (all bugs in phase 0), last. |

### Changed in place

DB.8, DB.7, NAV.10 (Alembic), ADM.15 (RQ workers), UX.3 (progress from RQ), UX.23 (astropy.units), ADM.18 and OPS.13 (naming key instead of word lists and corpus hashes), GEN.57 (no name collisions), PERF.18 and PERF.20 (RQ and the new caches), API.12 and VIEW.4 (naming key), NAV.5, MAP.70. Each has a "Plan (2026-10-07)" line in TODO.md.

### Boss's answers

Boss (2026-10-07 17:11Z, answering the open decisions; 2026-10-09 07:48Z, confirming them all). His words are quoted. Each answer is written into its TODO item.

- **Lane order**: "When one thread is idle waiting for CI we can start another thread on something else, the idea being to always have at least 1 thread going without even going tover the 2 dev 1 todo limit." So up to two build threads and the TODO thread; a thread waiting on CI may hand off to another.
- **Redis on Windows (OPS.21)**: "Let's say Redis in WSL". OPS.27 dropped Memurai from the installer and docs (done, PR #475).
- **Habitability score structure (GEN.84)**: approved (PHI-4's domains and tiers for display, the Xenobiology doc's three tiers as the scores behind them).
- **Wide-binary names under the codec (GEN.71)**: "No, we should never have A I or such for planet names.  Adjust the algorithm to produce 2 words from the name.  A says word 1 I, word 1 II, etc...  B planets say word 2 I, word 2 II, etc..." In GEN.71. Confirmed 2026-10-08 04:00Z: planets keep the "<star> I" pattern, so this rule stays in force.
- **GEN.29 outside phase 0**: approved (stays in phase 2 with the class refactor).
- **Front-end build (UX.40, UX.41, MAP.102)**: approved (vendored ES module builds served by Flask, no bundler).
- **Hilbert fill order (GEN.101)**: approved (keep the Hilbert order, allow a logged jump where the ball cuts the curve).
- **Orbital sectors and frames (GEN.109, GEN.115)**: "Approved, but verify the algorithm will work with our sector geometry." The check is in GEN.109.
- **Scaling the galaxy's gravity (GEN.115)**: approved.
- **Orbital epoch and step (GEN.105)**: "No 1 year per orbital update or turn, rather, once set up and configured we follow orbital paths in real time.  The update script should have an option to update for more time in 1 go if specified.  Default is 1 day = 1 day." In GEN.105 and GEN.106.
- **GEN.65**: "I do not have the error message, keep an eye out for it, but put it on the back burner for something to watch out for, design a test that will test for it in a variety of situations, I think this error occurred when I was attempting to generate  a neighborhood when the center was close to the edge of the galaxy." In GEN.65.

- **Wiki uploads on the queue (PERF.24)**: Boss (2026-10-08 00:14Z) reversed the audit decision: wiki uploads go on the Redis queue too, because he may want batch uploads later.
- **Prevalence fields (ADM.37, then ADM.45)**: Boss (2026-10-08 00:12Z): the Generate page shows each feature's real default share (habitable worlds 24.2%, asteroid belts 59%), not "0% change" (done, PR #514). Boss (2026-10-09 07:48Z): the fields are not plus or minus percentages: the user types the override percent, the form forces the set to total 100% and makes clear what to do, and there is no density-dependent share. This is ADM.45.

No decisions are open. (GEN.67, Boss 2026-10-08: stars and sectors keep word salad, planets, moons and belts keep the "<star name> I" pattern, and the codec names only the objects with no star-derived name; defaults stand for constellations and the naming key.)

## Files that several items touch

| Files | Items | Order |
|---|---|---|
| Every module (the package move) | none (OPS.24 done, PR #473) | The move is finished; the library swaps build on the new layout. |
| stellarObjects/workQueue.py, jobRunner.py, web/jobs.py | PERF.19, PERF.24, OPS.19, ADM.22, ADM.15, PERF.18, PERF.20 | PERF.19's audit, then PERF.24 with OPS.19, then ADM.22. |
| _db.py and migrateDb.py | DB.11, DB.13, GEN.71, DB.16, DB.17, API.11 (DB.7 and NAV.10 merged, PRs #813 and #814) | CI red fixes first; then DB.11; every later schema change is an Alembic migration, one writer at a time. |
| Names (names.py, bodyNames.py, nameUniqueness.py, objectId.py) | GEN.68 to GEN.73, VIEW.4, API.12 | One stream, in TODO order. |
| generate.py: qualify, density and backfill | GEN.98, GEN.100, GEN.101, GEN.42, PERF.18 (GEN.41 and GEN.43 retired in PR #812) | Phase 0 bugs first, then phase 1 galaxy gen. |
| Galaxy Map (galaxymap3d.js, galaxystageview.js, galaxystages.js, galaxyblocks.js) | MAP.102, MAP.65 to MAP.68, MAP.110, MAP.111, MAP.95, MAP.103, MAP.122, MAP.123, MAP.124 | Bugfix lane items first; the engine group next; phase 1 map items after. |
| Templates and components (base.html, style.css, edit_controls.html) | UX.2, UX.26, UX.31, UX.27, UX.49, ADM.34, UX.37, UX.21, UX.42, UX.43 | Components first, then the sweep, then wording and the visual design. |
| Planet physics and classes (planetPhysics.py, planetData.py, planetLife.py) | GEN.85 to GEN.89, GEN.33, GEN.28, GEN.27, GEN.29, GEN.91, GEN.92 | Phase 0 bug, then the habitability inputs, then the refactor. |
| Positions (updateOrbits.py, keplerMotion.py, the new position object) | GEN.74, GEN.66, GEN.104, GEN.106 to GEN.110, GEN.115, MAP.70, VIEW.5 | GEN.74 first. |
| Generate page (generate.html, generate_page.py) | UX.49, ADM.28 and its subitems, GEN.96, GEN.24 | Components, then prevalence, then the rework. |
| update.sh and update.ps1 (with install.* and deploy-common.*) | OPS.8, OPS.13, OPS.15, OPS.17 | OPS.7 (done, PR #457) went first; then Redis and pins. |
| Routing and settings files (scripts/bench_nav.py, db/corridor.py, galaxy/settings_file.py) | NAV.10 and ADM.18 (merged), NAV.12, NAV.52, NAV.53, ADM.42 | The merged items set the shape; NAV.12 and NAV.52 build on corridor.py, ADM.42 on settings_file.py. |

## Near-cycles and how they are broken

- **The library swaps and the package move**: both touch every module. Broken by order: the layout plan (OPS.23, done in PR #435) and the move (OPS.24) land first, then each swap.
- **The class refactor and the habitability index**: the classes need the scores and the scores read the classes' atmospheres. Broken by order: habitability inputs and score in phase 1 on today's classes, the refactor in phase 2.
- **MAP.61 and the map bugs**: the engine items moved into phase 0 because the 2026-10-03 map bugs fold into them.
- **API.6 and USR.2**, **NAV.4 and USR.1**: unchanged (keys and saved courses follow accounts).

## Judgment calls

- **The engine in phase 0**: MAP.65 to MAP.68, NAV.7, NAV.13 to NAV.15, NAV.29, NAV.33, NAV.32 and MAP.95 moved into phase 0 because the breadcrumb, empty-slab, filter and picker bugs are fixed by them, which is Boss's rule for architecture that takes care of a bug.
- **CI red first**: the five failures from PR #431 lead the bugfix lane so every later PR can be judged on green CI.
- **The prevalence group in phase 0**: GEN.48 is a bug, so its controls (GEN.52, ADM.16, TEST.75) came with it.
- **Kept in the bugfix lane, not folded**: UX.24, UX.29, UX.25 (already built once), OPS.6 (one check, not worth waiting for Pydantic), MAP.114 (CI red, though MAP.68 deletes sectormap.js later).
- **DB.13 in phase 0**: Boss called it "a Phaser 0 priority".
- Execution plan (Boss's priority order of 2026-10-09 07:13Z: Orbital updates, Maps, Nearby search, Routing; Repeatable galaxy last, blockers first): `docs/plan/execution-plan.html`, with lanes by wave and the Alembic and file conflicts between groups.
- Boss (2026-10-09 07:54Z): "I haven't held anything, so integrate ALL phase 1 items into the immediate TODO breakdown, focusing on unblocking things as a first priority." The execution plan's lane lists now hold all 77 non-bug Phase 1 items, ordered with the items that unblock the most work first (`docs/plan/execution-plan.html`).
- GEN.100's galaxy-wide scatter writes about 1.6e8 rows at the default galaxy size (phenomenon_scatter table). Boss may lower the scatter rate; the coordinator is asking him (2026-10-09 08:15Z).

## Design notes in docs/design

Written by the build and research threads; each is the reference for its
items. New notes are listed here when the PR that adds them merges.

| Note | Covers |
|---|---|
| [interstellar-object-rates.md](../design/interstellar-object-rates.md) | Real-world rates for interstellar objects (retuned in PR #804). |
| [star-types-by-galactic-radius.md](../design/star-types-by-galactic-radius.md) | GEN.133: every star type's rate against distance from the core, checked against the observed Milky Way; four proposals (young and intermediate populations weighted by the star-formation profile, a lower B share, a young tail in the bulge, a metallicity gradient only if planets need it), none built or filed yet (PR #807). |
| [compact-remnant-regions.md](../design/compact-remnant-regions.md) | Research on where neutron stars and black holes are in the galaxy, behind the regional rates (GEN.132, PR #807); added in PR #804. |
| [habitability-index.md](../design/habitability-index.md) | The habitability score structure, equipment profiles and tiers (GEN.84, PR #803). |
| [orbital-updates.md](../design/orbital-updates.md) | The orbital update design (GEN.105 and its chain). |
| [library-migration.md](../design/library-migration.md) | The move to third-party libraries. |
| [reproducible-galaxies.md](../design/reproducible-galaxies.md) | Seeds, version keys and the same-seed rule. |
| [generation-determinism.md](../design/generation-determinism.md) | Floating-point rules, rounding, fingerprints and CI for the same-seed guarantee (GEN.55, GEN.57 to GEN.61, TEST.77, OPS.12 to OPS.15, DB.7). |
| [db-check-and-parity-repair.md](../design/db-check-and-parity-repair.md) | Damage check, parity file and repair, Alembic progress and locks (DB.8 to DB.15). |
| [ops-scheduling-and-rotation.md](../design/ops-scheduling-and-rotation.md) | The daily maintenance run, schedules on three systems and the 18-slot backup rotation (OPS.8, OPS.13, OPS.14, OPS.16 to OPS.18, ADM.19, ADM.20, OPS.21, OPS.27). |
| [sampling-backfill-and-resume.md](../design/sampling-backfill-and-resume.md) | Backfill cost, sampling equivalence and resumable runs (GEN.40 to GEN.43, GEN.96 to GEN.99, GEN.102, PERF.18, PERF.29, PERF.30). |
| [fill-order-curves-and-core.md](../design/fill-order-curves-and-core.md) | Hilbert and other fill orders, region enumeration and the galactic core (GEN.101, ADM.29, ADM.30, GEN.97, GEN.23, GEN.24). |
| [course-routing.md](../design/course-routing.md) | Routes, the route graph, hop lengths and nearby search (NAV.9 to NAV.12, NAV.36, NAV.39 to NAV.48). |
| [course-avoidance.md](../design/course-avoidance.md) | Keep-out zones for courses around hazards (NAV.22, NAV.24 to NAV.28, NAV.51). |
| [orbital-solvers-and-integrators.md](../design/orbital-solvers-and-integrators.md) | Kepler solvers, N-body integrators, thresholds and scheduling (GEN.105 to GEN.109, GEN.111, GEN.115). |
| [galactic-potential.md](../design/galactic-potential.md) | The galactic gravity model and its gradient for object updates (GEN.108, GEN.115, GEN.106, GEN.125). |
| [collisions-and-mergers.md](../design/collisions-and-mergers.md) | Detecting and resolving collisions and mergers during an orbital update (GEN.109 to GEN.112). |
| [atmospheres-retention-and-classes.md](../design/atmospheres-retention-and-classes.md) | When a planet keeps an atmosphere, the new hot- and cold-zone classes and life stage (GEN.91, GEN.92). |
| [activity-magnetism-radiation-hydrosphere.md](../design/activity-magnetism-radiation-hydrosphere.md) | Stellar activity, magnetic fields, surface dose, water and oceans (GEN.86 to GEN.89, GEN.104). |
| [anomalies.md](../design/anomalies.md) | Anomalies and asteroid fields: classes, rates and generation rules (GEN.113, GEN.114, GEN.103, GEN.128 to GEN.130). |
| [multistar-and-compact-systems.md](../design/multistar-and-compact-systems.md) | Multi-star systems, binaries and compact-remnant systems (GEN.62, GEN.71, GEN.128 to GEN.130, GEN.103). |
| [population-and-politics.md](../design/population-and-politics.md) | Population, politics, technology levels and facilities (POP.7 to POP.10). |
| [api-design-standards.md](../design/api-design-standards.md) | API conventions: errors, paging, versioning, keys, logging and limits (API.3 to API.19, ADM.13). |
| [units-and-number-formatting.md](../design/units-and-number-formatting.md) | Units, unit ladders and number formatting in code and in the browser (UX.36). |
| [map-ui-and-frontend-libraries.md](../design/map-ui-and-frontend-libraries.md) | Map interface patterns and front-end library choices (UX.3, UX.40 to UX.49). |
| [settings-seo-and-accounts.md](../design/settings-seo-and-accounts.md) | Settings pages, search-engine behaviour and account design (ADM.42 to ADM.44, USR.1 to USR.8, API.6). |
| [sky-view.md](../design/sky-view.md) | The view from a planet's surface: stars, backfill levels and cost (VIEW.1 to VIEW.4). |
| [multiple-galaxies.md](../design/multiple-galaxies.md) | Neighbour galaxies in the sky and several galaxies in one install (GEN.9, VIEW.2, VIEW.3). |
| [navigation-frames.md](../design/navigation-frames.md) | Coordinate frames used for navigation and courses. |
| [galaxy-coordinate-system.md](../design/galaxy-coordinate-system.md) | The galaxy frame, rings, layers and slots. |
| [galaxy-disk-density.md](../design/galaxy-disk-density.md) | The disk density model behind sector fill. |
| [performance-eta-queue-and-caching.md](../design/performance-eta-queue-and-caching.md) | Profiling, progress and ETA estimates, the RQ queue and cache invalidation (PERF.20, PERF.31 to PERF.34, UX.3, ADM.15). |
| [exotic-environments-planets-and-compact-binaries.md](../design/exotic-environments-planets-and-compact-binaries.md) | Planets and compact-object orbits in extreme environments: Boss's research digested for GEN.94, GEN.95, GEN.130. |
| [globular-clusters.md](../design/globular-clusters.md) | Globular clusters: density-component storage, King sampling, metallicity, planet culls and pulsar planets; Boss's research digested with computed checks (GEN.9, GEN.130, GEN.134). |

Most of these notes are the research pass Boss asked for on 2026-10-09 (08:43Z): each gives its sources, marks numbers as seen in a search [S], computed [C] or recalled and still to verify [R], and names the TODO items it informs. The research environment could read search-result text but not the papers themselves, so every [R] needs checking once paper access is allowed.

## Not yet checked in live use

- **PERF.34 (PR #811)**: the fix (tile cache wiped every minute during a fill when more than 1,000 sectors or new bright stars forced a full rebuild; now a `busy` answer plus keeping the cache for up to 10 minutes) was found by reading the code path and has not been load-tested. Boss's next real fill is the check; reopen issue #638 if the site is still slow.
