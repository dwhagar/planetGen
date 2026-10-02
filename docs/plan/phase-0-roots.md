# Phase 0: Roots: bugs and groundwork

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

Fix the bugs nothing else depends on and lay the groundwork everything later builds on: the parallel path first (Boss: top priority), the 128-bit galaxy seed with its stored 22-digit version key and log line, the database consistency check, the binary-pair, forcing and name bugs, the map groundwork through the arc pick (MAP.85, which Boss wants in phase 0), the routing groundwork, object references and the small page and ops fixes. The database and physics bug threads are done (PRs #342, #347, #350).

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Parallel path (top priority)

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.76 | A bright-star test breaks on Python 3.9 and 3.10 (bug) |  | Same test file as 4 of PERF.21's 14 failures; fix first so the multi-worker run is clean on the py3.9 CI leg. |
| TEST.74 | Generation tests at more than one worker | TEST.76 | Lands first in the thread: shows what PERF.21 must fix. Fresh test databases at N workers hit DB.5's schema race. |
| PERF.22 | On Python 3.12 a run hangs forever when a worker process dies (bug) |  | workQueue._dispatch; same code as PERF.21. |
| PERF.21 | Generation works with any worker count: the parallel path is built, used and tested (bug) | TEST.74, PERF.22 | TOP PRIORITY (Boss). Rewrite fault injection so it reaches spawned workers instead of relying on the seeded single-process stream. |
| TEST.73 | Intermittent failure in the parallel galaxy-run interrupt test (bug) | TEST.74 | Race in the parallel interrupt path of workQueue.py; same thread. |
| PERF.23 | The bright-star progress bar can end at 101% (bug) |  | _LayerTracker in generate.py, same parallel scatter code. |
| GEN.32 | Re-running an interrupted bright-star band draws it twice (bug) | PERF.23 | Same scatter functions (_scatter_layers, add_bright_star_band) as PERF.23. Judgment: could instead ride GEN.44's per-sector levels. |
| GEN.39 | The same seed can't reproduce the same galaxy (bug) | PERF.21 | Decided yes (Boss 01:34-01:46Z): 128-bit seed stored in the database and logged at the top of every run; seed + version (packed in hex) reproduce the galaxy; needs a fresh galaxy. It per-sector seeded draws touch the RNG calls in every generator file, so it lands right after PERF.21 and before GEN.47, GEN.42 and PERF.18, which all need 'same seed, same stars'. |
| DB.6 | Store the galaxy's 128-bit seed, the version that made it, and every generation run | GEN.39 | Sub-item of GEN.39; fresh galaxy; 22-hex-digit version and environment key (Boss 02:08Z). |
| OPS.10 | The galaxy seed and version at the top of every generation log | DB.6 | Sub-item of GEN.39. |

### Binary pairs and single-system forcing

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.53 | The two stars of a binary don't share one age (bug) |  | StarSystem constructor in systemData.py (secondary made around line 270). |
| GEN.54 | A `--star-type` secondary gets a mass that doesn't fit its type (bug) | GEN.53 | Same lines as GEN.53; one PR. |
| GEN.49 | `+habitable_world` silently fails on hot stars (bug) | GEN.54 | Same constructor (the 8-attempt loop around line 353). Soft link: GEN.37 changes which stars can satisfy it. GEN.37 (PR #350) changed planet placement: re-measure first. |
| GEN.50 | `-planets +asteroid_belt` still makes an asteroid belt (bug) | GEN.49 | Same forcing code; rejected for single systems. |

### Galaxy geometry and names

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.90 | The tile-level helper crashes on a subnormal view radius (bug) |  | galaxyViewport.py and lib/galaxymap3d.py, the tile-level code MAP.80 also changes. |
| GEN.46 | Star system names of at most two words (bug) |  | nameUniqueness.py and _db.py name reservation; API.12 downloads the name state, so settle names first. Decided: only new names follow the rule. |
| NAV.38 | Every sector a straight line passes through |  | galaxyGeometry.py with a JS twin, after GEN.31 (same file). |

### Map groundwork

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.70 | Tests for the map JavaScript |  | Pins today's map behaviour before MAP.63/64 refactor and MAP.66/68 replace the Sector Map. |
| MAP.87 | Stars on the Sector Map and Galaxy Map need to be brighter, most of all the dim ones (bug) |  | _star_light in lib/starmap.py and STAR_LOG_LUMINOSITY in galaxymap3d.js; small, land before MAP.63 moves code. |
| MAP.83 | The "Mark rogue planets" button shows when it is on (bug) |  | Button style from aria-pressed (lib/starmap.py, sectormap.js). |
| MAP.82 | Unmarked rogue planets barely visible (bug) |  | Rogue point style in sectormap.js / lib/starmap.py. |
| MAP.84 | Marked rogue planets grow and become clickable; unmarked ones stay small (bug) | MAP.82, MAP.83 | Same code; marked vs unmarked size and hit area. |
| NAV.30 | Hide "View phenomenon" and "View system" links while picking a course (bug) |  | Hide 'View phenomenon/system' in pick mode (sectormap.js line 116, galaxymap3d.js line 411). Tiny; fix now, TEST.70 keeps it. |
| MAP.81 | Ctrl+1 to Ctrl+9 bookmark keys clash with the browser's tab switching (bug) |  | Decision gate (which keys). bookmarks.js; decide before MAP.64 builds the shared key map. |
| MAP.63 | Shared map helpers in one module | TEST.70, MAP.87, MAP.84, NAV.30 | Moves helpers out of galaxymap3d.js, sectormap.js and systemmap.js; no visible change. |
| MAP.64 | One camera and input controller | MAP.63, MAP.81 | One controller with zoom policies; MAP.53, MAP.58, MAP.75, MAP.73 build on it. |
| MAP.60 | Galaxy Map scale readout: one scale line | MAP.64 | Start of the one ordered Galaxy Map thread. |
| MAP.55 | Galaxy Map buttons: a menu, with only back, forward, up, reset and bookmark showing | MAP.60, UX.28, MAP.81 | Menu button and bookmark button; icons from UX.28, keys from MAP.81. |
| MAP.85 | The galaxy pick is an arc, on a 3D galaxy with no sector lines | MAP.55, MAP.64 | Moved to phase 0 by Boss (01:46Z). Root of the new selection: every later Galaxy Map pick item is rewritten around it. Arc size default: about 40 degrees by a third of the radius. Ships with today's block shading until MAP.86 lands. |
| MAP.52 | Galaxy Map highlights the wrong area; pick a 40-degree wedge around the cursor (bug) | MAP.85 | Bug, but its code is replaced by MAP.85; its width and snapping carry into the arc. Same PR as MAP.85. |

### System Map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.57 | The System Map writes NaN or infinite positions into its SVG (bug) |  | lib/systemmap.py. |
| MAP.88 | Parts of a star system run off the edge of the System Map (bug) | MAP.57 | Fit the whole scene; same file as MAP.57 and MAP.89. |

### Object references

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.7 | One reference for every object, with its parents |  | Root of the picker, saved courses, the 3D system view and account bookmarks. |

### Ops and flakes

| ID | Item | Needs | Note |
|---|---|---|---|
| OPS.6 | Admin scripts accept impossible `--mysql-port` values (bug) |  | _db.add_mysql_connection_args. |
| OPS.7 | Update asks to fill a wiped database with population data (bug) |  | update.sh / update.ps1; one PR with OPS.8. |
| OPS.8 | Update reloads Apache itself when run as root | OPS.7 | Not a bug, but the same files as OPS.7, so it rides along. |
| TEST.71 | Intermittent failure in the admin planet-regenerate test (bug) |  |  |
| TEST.72 | Intermittent failure in the two-step (2FA) sign-in test (bug) |  |  |
| UX.34 | The sector summary calls white dwarfs "B-type" and "A-type" systems (bug) |  | sector_generation_summary_lines in generate.py; one PR with OPS.9. |
| OPS.9 | Multi-line messages lose their prefix in the debug log (bug) | UX.34 | Same summary record. |
| API.15 | Log every API call with its user, how it came in, and its HTTP response code |  | Needs no user accounts (Boss 01:31Z). |
| TEST.78 | A resume test's sector query fails under ONLY_FULL_GROUP_BY on MariaDB 10.11 (bug) |  | From the Database thread (PR #342). |

### Page groundwork and small page bugs

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.28 | Investigate icons instead of words on buttons |  | Survey and icon sprite; UX.25, UX.26, UX.27 and MAP.55 use its icons. |
| UX.23 | A shared unit-ladder module |  | Shared unit ladder; UX.22 then UX.30 build on it. |
| UX.2 | Menus sized to what they hold (bug) |  | style.css menus; independent. |
| UX.24 | Sector contents: rogue planets after systems and phenomena, and expanded rows the full table width (bug) |  | sector_page.py / sector.html; before UX.25 and UX.26 change the same table and template. |
| UX.29 | Every comet in a system shows its type as a link (bug) |  | _comet_row_html in lib/systempage.py; before UX.30 rebuilds the rows. |
| ADM.14 | Line up the Generate page's text boxes, not their headings (bug) |  | generate.html field layout; before ADM.16 and GEN.24 add fields to the same page. |

### Routing groundwork

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.79 | Route edge cases, written before NAV.12 |  | Cases from the hop-length study's report. |
| NAV.34 | Courses between separately generated areas find no route: the route graph splits into islands (bug) |  | navGraph.build_knn_adjacency; needed before NAV.12. |

### Database consistency check

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.8 | Check a galaxy database and say whether it is damaged |  | Boss 02:13Z: phase 0, its own thread. Read-only. Stats and version checks switch on once GEN.44/PERF.11 and DB.6/DB.7 land; no hard dependency. |

The parallel path thread starts first (Boss: top priority). The map
groundwork thread runs to MAP.85 (the arc pick) and MAP.52 in one PR;
MAP.60 onward can be a second thread once MAP.64 merges, since they
share files. MAP.55 uses UX.28's icons; if UX.28 isn't done, MAP.55
ships with text labels and gets icons later. Until MAP.86 (phase 1)
lands, the arc pick keeps today's block shading. Arc size: about 40
degrees of bearing by a third of the disk radius (about 27 arcs). At
most 4 build threads run at once (Boss).

## Open questions for Boss

- MAP.52: Galaxy Map highlights the wrong area; pick a 40-degree wedge around the cursor (bug), see its entry in TODO.md.
