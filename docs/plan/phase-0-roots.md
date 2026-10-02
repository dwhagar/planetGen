# Phase 0: Roots: bugs and the groundwork they need

Rebuilt on 2026-10-02 from the dependency report (every open item, its
prerequisites and the files it shares), with Boss's decisions of that
night. [docs/TODO.md](../TODO.md) is the master file: it holds each
item's full text, and its "Plan: phases" section indexes every phase.
This file gives the phase's goal, its build threads with each item's
prerequisites in order, and its open questions. Where the two disagree,
TODO.md wins; when an item ships, it leaves TODO.md and its row here is
deleted in the same PR. Research notes, the files several items share and
the judgment calls behind the placement are in [notes.md](notes.md).

## Goal

Bug fixes first (Boss 04:45Z: phase 0 is primarily bug fixes and the groundwork that goes with them; 04:55Z: "Bugs first"): the binary-pair and forcing bugs, the Galaxy Map follow-ups and drill-down fixes, the System Map, the generation bugs (with class S for GEN.38), sector stats and colors (with the per-sector stats table MAP.86 needs), the routing groundwork, the sector and system page bugs, the small page bugs, and ops and test flakes last. The parallel path lane is done (GEN.39, DB.6, OPS.10: PRs #381, #387, #391).

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Binary pairs and single-system forcing

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.85 | Name collisions can count -1 existing names and fail generation (bug) |  | Real bug, not a flake: hit 3+ tests in 2 files (bright-star scatter, galaxy gen; PR #373, PR #393 runs). nameUniqueness.py:137 rejects the -1 count _db.py produced. Root cause: _db.reserve_system_names counts a name redrawn in the same pass against a later row; fix by computing name keys once per pass (also 1.7x faster on dense sectors). After GEN.51 in the naming lane; GEN.57/GEN.63 build on a right count. |

### System Map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.57 | The System Map writes NaN or infinite positions into its SVG (bug) |  | lib/systemmap.py. |
| MAP.88 | Parts of a star system run off the edge of the System Map (bug) | MAP.57 | Fit the whole scene; same file as MAP.57 and MAP.89. |
| MAP.92 | The System Map's side panel leaves out a planet's or moon's radius and mass (bug) | MAP.88 | Boss 04:29Z. showInfo in systemmap.js and the marker attributes in lib/systemmap.py; same files as MAP.57/MAP.88. |

### Galaxy Map drill-down

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.56 | Drop the 3x3 block pick: select a slab, zoom in, select a segment (bug) |  | Moved to phase 0: MAP.85 (the code that would have replaced it) is done in PR #369, so the bug can be fixed now. Arc, slab, segment ladder in galaxystages.js. |
| MAP.76 | Leader-line layout | MAP.56 | Moved to phase 0 as groundwork: the layout half of MAP.54, same PR. Layout half of MAP.54; same PR. |
| MAP.54 | Slab leader lines instead of the slab slider (bug) | MAP.56, MAP.76 | Moved to phase 0: MAP.85 (the code that would have replaced it) is done in PR #369, so the bug can be fixed now. |
| MAP.53 | Rotate a zoomed-in wedge, and zoom it to fit the window (bug) |  | Moved to phase 0: MAP.85 (the code that would have replaced it) is done in PR #369, so the bug can be fixed now. Rotation and fit on MAP.64's controller. |
| MAP.78 | Zooming into a wedge must show the whole wedge at every drill-down level (bug) | MAP.53 | Moved to phase 0: MAP.85 (the code that would have replaced it) is done in PR #369, so the bug can be fixed now. The fit itself; with MAP.53. |
| MAP.77 | Galaxy Map draws block divisions inside a picked slab before zooming to it (bug) | MAP.56 | Moved to phase 0: MAP.85 (the code that would have replaced it) is done in PR #369, so the bug can be fixed now. Which lines show at each level, once the ladder and ghost exist. Dependency on MAP.59 dropped: slab-only lines need no ghost; MAP.59 keeps this rule instead. |
| MAP.96 | The Galaxy Map can't be turned freely: the tilt stops at straight down and at 80 degrees (bug) | MAP.53 | Boss 05:12Z. clampTilt (TOP_DOWN_PHI to MAX_TILT 80 deg) in galaxystageview.js; turn any way at every level, trackball-style; picks keep working. |
| MAP.97 | The Galaxy Map camera should go top-down for the galaxy and a slab, isometric for a block, at every zoom step (bug) | MAP.56, MAP.96 | Boss 05:12Z. cameraFor presets: galaxy and slab top-down (replaces GALAXY_TILT 35 deg), block isometric, animated both ways; open questions on carrying a manual turn over and returning to the preset. |

### Generation bugs

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.60 | Rogue gas giants get a Jupiter-sized radius at every mass (bug) |  | Moved to phase 0: a bug with nothing ahead of it. From the Physics bugs thread (PR #350): use GEN.34's giant mass-radius relation in roguePlanetData.py. |
| GEN.38 | Rocky rogue planets over 10,000 km are still classed C (bug) |  | Moved to phase 0 with its groundwork: class S (rocky super-Earth, rogue-eligible) is built here as GEN.28's first class PR. Probably solved by a rogue-eligible S class. |
| GEN.47 | Nebulae almost never appear (bug) |  | Moved to phase 0: a bug whose prerequisite (GEN.39 per-region seeds) is done. Galaxy-scale nebula field spanning sectors: every worker and every later run must agree where a cloud is, so it needs deterministic per-region draws (GEN.39, or an address hash if GEN.39 is dropped). |

### Sector stats and colors

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.44 | Store each sector's backfill level so finished sectors drop out of any backfill |  | Moved to phase 0 as groundwork for the MAP.86 bug (Boss: MAP.86's color goes in the shared per-sector stats table). One shared per-sector stats table with PERF.11 (Boss 01:46Z). The level also drives the scatter bands: -1 with stars means a failed run (wipe and redo); otherwise draw only between the new floor and the stored level (Boss 03:25Z). Galaxy schema v53 (one writer at a time; v52 is DB.6's version key and run history). Backfill code shared with PERF.18 and GEN.42. |
| PERF.11 | Store each sector's expected and actual density | GEN.44 | Moved to phase 0 as groundwork for the MAP.86 bug (Boss: MAP.86's color goes in the shared per-sector stats table). Same per-sector stats table as GEN.44 (Boss 01:46Z); MAP.86's color goes there too. |
| PERF.1 | Generation at scale | PERF.11 | Moved to phase 0 as groundwork for the MAP.86 bug (Boss: MAP.86's color goes in the shared per-sector stats table). Parent; only PERF.11 is open under it. |
| MAP.80 | Sector-level zoom on the Galaxy Map should show almost every star in the sector (bug) |  | Moved to phase 0: a bug with nothing ahead of it; tile listing, before MAP.86 in the same tile files. Judgment: moved up from the selection chain; the thinning is in the tile listing (queryDb GALAXY_TILE_* floors) and galaxymap3d.js, not the pick code. |
| MAP.86 | Sector and block colors from what is in them: filled sectors translucent (bug) | PERF.11 | Moved to phase 0 with its groundwork GEN.44 and PERF.11. Next after MAP.85: with no lines, color carries the structure. Its color goes in the shared per-sector stats table (GEN.44 + PERF.11). Tile payload change: bump the tile cache. |

### Routing groundwork

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.79 | Route edge cases, written before NAV.12 |  | Cases from the hop-length study's report. |
| NAV.34 | Courses between separately generated areas find no route: the route graph splits into islands (bug) |  | navGraph.build_knn_adjacency; needed before NAV.12. |

### Sector and system pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.28 | Investigate icons instead of words on buttons |  | Survey and icon sprite; UX.25, UX.26 and UX.27 use its icons; MAP.55 (done, PR #369) left text labels with an icon hook. |
| UX.24 | Sector contents: rogue planets after systems and phenomena, and expanded rows the full table width (bug) |  | sector_page.py / sector.html; before UX.25 and UX.26 change the same table and template. |
| UX.29 | Every comet in a system shows its type as a link (bug) |  | _comet_row_html in lib/systempage.py; before UX.30 rebuilds the rows. |
| UX.25 | Rogue planets: octant and a small map symbol beside each name (bug) | UX.24, UX.28 | Moved to phase 0: bugs whose groundwork (UX.24, UX.28) is in phase 0. Same sector table as UX.24. |
| UX.26 | Edit and admin actions as a button that opens a menu (bug) | UX.28, UX.25 | Moved to phase 0: bugs whose groundwork (UX.24, UX.28) is in phase 0. Sector page admin panel and edit_controls.html. |
| UX.31 | Editing a star system: an edit button with a quick menu, not a long panel (bug) | UX.26 | Moved to phase 0: bugs whose groundwork (UX.24, UX.28) is in phase 0. system.html edit panel (_edit_rows in system_pages.py). |
| UX.27 | System page: the system and navigation buttons on one row that doesn't overlap (bug) | UX.28 | Moved to phase 0: bugs whose groundwork (UX.24, UX.28) is in phase 0. system.html subhead; shares wording with NAV.29. |

### Small page bugs

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.2 | Menus sized to what they hold (bug) |  | style.css menus; independent. |
| ADM.14 | Line up the Generate page's text boxes, not their headings (bug) |  | generate.html field layout; before ADM.16 and GEN.24 add fields to the same page. |
| NAV.41 | The NAV page's course map is too small to read (bug) |  | Boss 04:19Z. navmap.py 360-unit square at 22rem, 9px labels; widen and enlarge. |
| UX.36 | Scientific notation starts too early for whole numbers (bug) |  | Boss 04:29Z. numberformat.js and utils.py: whole numbers scientific from 7 digits, decimals from 5. |
| UX.33 | Filter phenomena by their classes and types (bug) | GEN.47 | Moved to phase 0: filters are built from the class lists, so new classes (GEN.28) appear on their own; after GEN.47 so nebula classes exist. Filters over class lists that GEN.28 and GEN.47 change. |
| UX.38 | The nebula and remnant diagrams' "-" button does nothing at the 1 ly limit (bug) |  | Split out of UX.21 (its known dead control): lib/phenomenonmap.py lines 53 and 127, the mapzoom.js clamp; strict xfail in test_web_browser_maps.py. |

### Ops and flakes

| ID | Item | Needs | Note |
|---|---|---|---|
| OPS.6 | Admin scripts accept impossible `--mysql-port` values (bug) |  | _db.add_mysql_connection_args. |
| OPS.7 | Update asks to fill a wiped database with population data (bug) |  | update.sh / update.ps1; one PR with OPS.8. |
| OPS.19 | The Generate jobs folder is /var/lib/planetgen while the checkout is /var/lib/planetGen (bug) |  | jobs.py and deploy-paths.py defaults; update.sh moves an old lowercase jobs folder. Same update.sh as OPS.7/OPS.8. |
| UX.34 | The sector summary calls white dwarfs "B-type" and "A-type" systems (bug) |  | sector_generation_summary_lines in generate.py; one PR with OPS.9. |
| OPS.9 | Multi-line messages lose their prefix in the debug log (bug) | UX.34 | Same summary record. |
| TEST.71 | Intermittent failure in the admin planet-regenerate test (bug) |  |  |
| TEST.72 | Intermittent failure in the two-step (2FA) sign-in test (bug) |  |  |
| TEST.78 | A resume test's sector query fails under ONLY_FULL_GROUP_BY on MariaDB 10.11 (bug) |  | From the Database thread (PR #342). |
| TEST.80 | Intermittent failure in the admin change-star test (bug) |  |  |
| TEST.81 | Two processes reserving id blocks of one table can deadlock (bug) |  | _db._reserve_id_block, 1213 deadlock on MariaDB 10.11. |
| TEST.82 | Intermittent failure in the orbit-ceiling trim test (bug) |  | test_validation.py, physics area. |
| TEST.83 | The sign-in rate-limit test fails under parallel load (bug) |  | test_web_admin.py; passes alone, fails under -n auto. |
| TEST.84 | The every-column round-trip test depends on whether a quasar got placed (bug) |  | NULL_IN_THIS_GALAXY depends on the draw; about 1 in 5 on MySQL 8.0. |
| TEST.86 | Intermittent failure in the concurrent-insert recovery test (bug) |  | test_galaxy_gen.py; failed once under -n auto. |

At most two build threads run at once (Boss 02:51Z). Done lanes:
Parallel path (PRs #381, #387, #391) and Galaxy Map follow-ups (PRs
#395, #399); the binary pairs lane has TEST.85 left. Then the lanes
start in the order above as a slot frees: System Map, Galaxy Map
drill-down, Generation bugs, Sector stats and colors (MAP.80 and
MAP.86 after the drill-down merges), Routing groundwork, Sector and
system pages (text labels with an icon hook if UX.28's icon list is
not approved yet), Small page bugs, and Ops and flakes last. NAV.7
and DB.8 open phase 1.

## Open questions for Boss

- MAP.56: Drop the 3x3 block pick: select a slab, zoom in, select a segment (bug), see its entry in TODO.md.
- MAP.54: Slab leader lines instead of the slab slider (bug), see its entry in TODO.md.
- MAP.53: Rotate a zoomed-in wedge, and zoom it to fit the window (bug), see its entry in TODO.md.
- MAP.97: The Galaxy Map camera should go top-down for the galaxy and a slab, isometric for a block, at every zoom step (bug), see its entry in TODO.md.
- PERF.11: Store each sector's expected and actual density, see its entry in TODO.md.
