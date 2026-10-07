# Phase 0: Roots: every bug, and the groundwork the new architecture needs

Rebuilt on 2026-10-07 from Boss's lists of 2026-10-03 and 2026-10-07
and the dependency tree (every open item, its prerequisites and the
files it shares). [docs/TODO.md](../TODO.md) is the master file: it holds each
item's full text, and its "Plan: phases" section indexes every phase.
This file gives the phase's goal, its build threads with each item's
prerequisites in order, and its open questions. Where the two disagree,
TODO.md wins; when an item ships, it leaves TODO.md and its row here is
deleted in the same PR. Research notes, the files several items share and
the judgment calls behind the placement are in [notes.md](notes.md).

## Goal

Two lanes (Boss 2026-10-03 05:38Z: "the most fundamental and needed changes first in phase 0 along with the bug fixes in 2 lanes, bugfixes and groundwork"). Bugfixes: the CI failures first, then generation, console and progress, maps and pages, prevalence (GEN.48), and ops and test flakes. Groundwork: the package layout and the move to third-party libraries with Redis (Boss's explicit directive), the data model (SQLAlchemy and Alembic, values in columns not JSON, one point-in-space object, Pydantic, scipy and astropy), names from IDs, the RQ queue with streamed logs and progress, Shoelace and TanStack components, the shared map engine, nebula shapes, and the UX sweep. Bugs that the groundwork fixes are folded into it and listed under it. GEN.65 stays held until Boss sends the error text.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Bugfixes: CI red

Done: all eight items landed in PR #442 (2026-10-07).

### Bugfixes: generation

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.65 | A generation run fails from the web UI but not from the CLI (bug) |  | Held by Boss until he gives the error text; may be the same failure as GEN.76 (empty sectors). Boss 08:08Z: high priority, top of phase 0, not started yet. Details unknown; ask Boss for the error. |
| PERF.26 | Size estimates don't match what generation stores (bug) |  |  |
| GEN.77 | Neighborhood generation fails when its first sector is below the star threshold (bug) |  |  |
| GEN.78 | Some regions have a star probability of zero (bug) |  | Also covers Boss's 2026-10-07 "Star generation should always actually take place" (merged into GEN.76 and here). |
| GEN.79 | Bright stars only land between layers -121 and 121, so the bulge never shows (bug) | GEN.78 | Major bug (Boss 2026-10-07); merges the 2026-10-03 bulge report. |

### Bugfixes: console and progress

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.81 | The console refuses runs instead of warning and doing what was asked (bug) |  |  |
| ADM.33 | The owner can override "no room" warnings and generate anyway | GEN.81 | Groundwork for GEN.81 on the web side. |

### Bugfixes: maps and pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.28 | Investigate icons instead of words on buttons |  | Icon set approved by Boss (2026-10-02); the icon sprite feeds UX.40. Survey and icon sprite; UX.25, UX.26 and UX.27 use its icons; MAP.55 (done, PR #369) left text labels with an icon hook. |
| UX.24 | Sector contents: rogue planets after systems and phenomena, and expanded rows the full table width (bug) |  | Built in the paused Sector and system pages lane, committed only in its container (may be lost). sector_page.py / sector.html; before UX.25 and UX.26 change the same table and template. |
| UX.29 | Every comet in a system shows its type as a link (bug) |  | Built in the paused lane, committed only in its container (may be lost). _comet_row_html in lib/systempage.py; before UX.30 rebuilds the rows. |
| UX.25 | Rogue planets: octant and a small map symbol beside each name (bug) | UX.24, UX.28 | Half built in the paused lane. 24, UX.28) is in phase 0. Same sector table as UX.24. |
| NAV.41 | The NAV page's course map is too small to read (bug) |  | Boss 04:19Z. navmap.py 360-unit square at 22rem, 9px labels; widen and enlarge. |
| UX.36 | Scientific notation starts too early for whole numbers (bug) |  | Boss 04:29Z. numberformat.js and utils.py: whole numbers scientific from 7 digits, decimals from 5. |
| SEC.31 | Signing in as admin works but shows a "form expired" error (bug) |  |  |
| UX.44 | Search: mutually exclusive tags should combine with OR, the rest with AND (bug) |  |  |
| MAP.117 | Surface pressure missing from the planet and moon side panel (bug) |  |  |
| MAP.118 | The Galaxy Map shows an unfilled sector above the galaxy that can't be filled (bug) |  |  |

### Bugfixes: prevalence

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.52 | Prevalence controls for sector and galaxy runs |  | Moved into phase 0 with GEN.48 (all bugs in phase 0). Probability adjustments reach systemData.py, the same constructor phase 0 thread D fixed. |
| TEST.75 | Tests for forcing and prevalence | GEN.52 | Moved into phase 0 with GEN.48. Grows with GEN.49 to GEN.52. |
| ADM.16 | Prevalence controls on the Generate page | GEN.52, ADM.14 | Moved into phase 0 with GEN.48. generate.html, after ADM.14's layout. |
| GEN.48 | Forcing options are impractical for whole sectors; replace them with prevalence controls (bug) | GEN.52, ADM.16, TEST.75 | Moved into phase 0: Boss wants every bug in phase 0, and its controls come with it. Parent; closes with its subitems. |

### Bugfixes: ops and flakes

| ID | Item | Needs | Note |
|---|---|---|---|
| OPS.6 | Admin scripts accept impossible `--mysql-port` values (bug) |  | _db.add_mysql_connection_args. |
| OPS.7 | Update asks to fill a wiped database with population data (bug) |  | update.sh / update.ps1; one PR with OPS.8. |
| TEST.78 | A resume test's sector query fails under ONLY_FULL_GROUP_BY on MariaDB 10.11 (bug) |  | From the Database thread (PR #342). |
| TEST.82 | Intermittent failure in the orbit-ceiling trim test (bug) |  | test_validation.py, physics area. |
| TEST.84 | The every-column round-trip test depends on whether a quasar got placed (bug) |  | NULL_IN_THIS_GALAXY depends on the draw; about 1 in 5 on MySQL 8.0. |
| TEST.86 | Intermittent failure in the concurrent-insert recovery test (bug) |  | test_galaxy_gen.py; failed once under -n auto. |
| TEST.88 | The facilities test fails when the drawn gas giant's sphere of influence is too small (bug) |  | test_facilities.py; a random giant's sphere of influence can be under the test's 500,000 km orbit. |
| TEST.89 | The Galaxy Map drill-down browser test fails intermittently (bug) |  | test_web_browser_maps.py drill-down by clicks; failed once under -n auto on MariaDB 10.11, passed 3 of 3 alone. |

### Groundwork: layout and libraries

| ID | Item | Needs | Note |
|---|---|---|---|
| OPS.24 | Move the code into the new package layout, one package per PR |  | Mechanical moves with no shims or wrappers (Boss 13:04Z): each PR updates every caller; every other open branch merges main after each one. |
| OPS.22 | Reorganize the code into importable Python packages with shared utility libraries | OPS.24 | Parent; closes with its subitems. Goes before the library swaps so each file moves once. |
| SEC.29 | Two-step sign-in on pyotp, QR codes on segno |  | Folds TEST.72 (the 2FA flake). |
| TEST.72 | Intermittent failure in the two-step (2FA) sign-in test (bug) | SEC.29 | Folded into SEC.29 (pyotp). |
| SEC.30 | Login and request rate limits on Flask-Limiter with Redis storage |  | Folds TEST.83 (rate-limit tests under load). |
| TEST.83 | Rate-limit tests fail under parallel load (bug) | SEC.30 | Folded into SEC.30 (Flask-Limiter). test_web_admin.py; passes alone, fails under -n auto. |
| UX.39 | Markdown rendered by the markdown library |  |  |
| OPS.20 | Move the code base from zero dependencies to third-party open-source libraries | SEC.29, UX.39, PERF.24, SEC.30, PERF.25, DB.11, ADM.21, GEN.66, UX.40, UX.41, ADM.22, MAP.102 | Parent of the library migration; closes with its subitems. Boss's explicit directive. |

### Groundwork: data model

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.11 | The database layer and migrations on SQLAlchemy and Alembic | OPS.24 | Every later schema change (DB.7, NAV.10, API.11, DB.13, GEN.106) is an Alembic migration. |
| DB.13 | Every stored value in its own column, not in JSON blocks, and indexed for search | DB.11 | Boss (2026-10-07 11:47Z): "This is to be considered a Phaser 0 priority." |
| ADM.21 | Input validation on Pydantic models | OPS.24 | Also the error objects API recipes (API.18) return. |
| GEN.66 | Physics on scipy, and astropy constants and units | OPS.24 |  |
| GEN.74 | One point-in-space object that keeps every coordinate system in step, used by every object | OPS.24 | Groundwork for the orbital work, Hill-sphere placement and light-travel positions. |

### Groundwork: names from IDs

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.68 | Research: the cheapest unique IDs for every object, unfilled sectors included |  | Asked on 2026-10-07; decides how GEN.69 builds IDs. |
| GEN.69 | A unique ID for every object, star systems and unfilled sectors included | GEN.68 |  |
| GEN.70 | A naming key in the control database, made at galaxy creation and changeable by admin |  |  |
| GEN.71 | Remove the word-salad name code, nltk and the name registries | GEN.69, GEN.70 |  |
| TEST.71 | Intermittent failure in the admin planet-regenerate test (bug) | GEN.71 | Folded into GEN.67: codec names have no apostrophes; fix the test's escaping there. |
| GEN.72 | A backfilled bright star should get a name only when its sector is generated (bug) | GEN.69 |  |
| GEN.73 | Nebulae don't get unique names (bug) | GEN.69 |  |
| GEN.67 | Names from IDs: replace word-salad name generation | GEN.68, GEN.69, GEN.70, GEN.71, GEN.72, GEN.73 | Parent; folds the two naming bugs, TEST.71, and drops GEN.63. |

### Groundwork: queue, logs and caches

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.19 | Everything the API or web site starts runs on the work queue (investigate) |  | Moved into phase 0 as the first step of PERF.24 (Boss chose Redis, 2026-10-03). Audit only; nothing moves to the queue until it runs at any worker count. |
| PERF.24 | The work queue and web jobs on Redis with RQ | PERF.19 | Replaces workQueue.py and jobRunner.py; folds OPS.19 and settles PERF.19. |
| OPS.19 | The Generate jobs folder is /var/lib/planetgen while the checkout is /var/lib/planetGen (bug) | PERF.24 | Folded into PERF.24: the job store moves with the queue. jobs.py and deploy-paths.py defaults; update.sh moves an old lowercase jobs folder. Same update.sh as OPS.7/OPS.8. |
| PERF.25 | The page cache on cachetools; the tile cache stays |  | tilecache.py stays (JSON only, Boss 13:27Z); diskcache, sqlitedict, cachelib and Flask-Caching are out. PERF.20 plans short-term API caching on top of it. |
| ADM.22 | Job logs streamed over SSE into Xterm.js, with native progress bars | PERF.24 | Folds the four log and progress bugs below. |
| ADM.23 | Log output wraps with hard line breaks (bug) | ADM.22 |  |
| ADM.24 | A failed action's log closes before it can be read (bug) | ADM.22 |  |
| ADM.25 | Error tracebacks don't reach the console and the web log window (bug) | ADM.22 |  |
| ADM.26 | The bright-star backfill shows no progress bar on the web (bug) | ADM.22 |  |

### Groundwork: web components

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.40 | Buttons, menus and dialogs from Shoelace web components | UX.28 | Folds UX.2, UX.26, UX.27, UX.31 and ADM.14; uses UX.28's approved icons. |
| UX.26 | Edit and admin actions as a button that opens a menu (bug) | UX.28, UX.25, UX.40 | Folded into UX.40; first screen of ADM.34. 24, UX.28) is in phase 0. Sector page admin panel and edit_controls.html. |
| UX.31 | Editing a star system: an edit button with a quick menu, not a long panel (bug) | UX.26 | Folded into UX.40. 24, UX.28) is in phase 0. system.html edit panel (_edit_rows in system_pages.py). |
| UX.27 | System page: the system and navigation buttons on one row that doesn't overlap (bug) | UX.28, UX.40 | Folded into UX.40. 24, UX.28) is in phase 0. system.html subhead; shares wording with NAV.29. |
| UX.2 | Menus sized to what they hold (bug) | UX.40 | Folded into UX.40. style.css menus; independent. |
| ADM.14 | Line up the Generate page's text boxes, not their headings (bug) | UX.40 | Folded into UX.40. generate.html field layout; before ADM.16 and GEN.24 add fields to the same page. |
| UX.41 | Tables on TanStack Table and TanStack Virtual |  | Folds UX.33 (phenomena filters). |
| UX.33 | Filter phenomena by their classes and types (bug) | UX.41 | Folded into UX.41. 28) appear on their own; GEN.47 is done (PR #419), so nebulae exist and nothing blocks it. Filters over class lists that GEN.28 changes. |
| ADM.34 | One admin menu per screen, holding only that screen's actions | UX.40, UX.26, UX.31 | UX.26 and UX.31 are its first two screens. |

### Groundwork: map engine

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.7 | One reference for every object, with its parents |  | Moved into phase 0: the engine (MAP.67, NAV.13) needs it, and the engine fixes the map bugs. Root of the picker, saved courses, the 3D system view and account bookmarks. |
| MAP.95 | A "Forward to current" button next to the map's Back and Forward |  | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. 94, next in the lane. Boss 04:29Z. Jumps to maxIndex of the map history (MAP.26). |
| NAV.13 | A picker module: select, step out, step in, step sideways | NAV.7 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. picker.js; needs no engine. |
| NAV.14 | One breadcrumb for every level | NAV.13 | Its breadcrumb fixes MAP.106. |
| MAP.65 | One picking, hover and info-panel layer |  | Moved into phase 0: fixes MAP.108, MAP.107, MAP.112 and NAV.46. After the selection rewrite settles the pick flow. |
| MAP.79 | Rogue planets clog the Sector Map: dim them, and a show/hide button per kind of object (bug) | MAP.65 | Moved into phase 0 (bug); its nebula toggle closes MAP.113. Judgment: the per-kind toggles go on the shared control set; the dimming already landed in phase 0 (MAP.82 to MAP.84). |
| NAV.15 | Pick mode everywhere | NAV.13, NAV.14, MAP.65 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. Pick mode in the shared panel layer. |
| NAV.29 | Replace "Nav from here" and "Nav to here" with "Start Here" and "End Here" while picking (bug) | NAV.13, NAV.15 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. Needs the step out/in of NAV.13. |
| NAV.33 | After picking one end of a course, stay at that zoom level (bug) | NAV.15 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. |
| MAP.66 | The sector as the drill-down's last stage, on the same page | MAP.65 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. Sector as the last stage. |
| MAP.67 | One URL and history scheme for every level | MAP.66, NAV.7 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. |
| MAP.68 | Remove the old Sector Map code | MAP.67, MAP.79 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. Deletes sectormap.js. |
| MAP.61 | One map engine and control set for the Galaxy Map and the Sector Map | MAP.65, MAP.66, MAP.67, MAP.68 | With the map engine in phase 0: the 2026-10-03 map bugs fold into it. Parent; closes with its subitems. |
| NAV.32 | Every Galaxy and Sector Map control works on the navigation screens (bug) | MAP.68, NAV.15, NAV.29, NAV.33, MAP.79 | Moved into phase 0 with the engine (all bugs in phase 0). Bug, but it is the 'pick mode uses the one engine' end state. |
| MAP.102 | Galaxy Map streaming with a BVH and 3D tiles, and camera-relative rendering |  | Folds the slow zoom and dense-sector crowding bugs. |
| MAP.107 | Selecting an empty slab near the core says "There is no layer x here" (bug) | MAP.65 |  |
| MAP.106 | The breadcrumb trail falls out of sync with the map (bug) | NAV.14, MAP.67, MAP.107 | Fixed by the engine's one URL and history scheme (MAP.67) and NAV.14's breadcrumb. |
| MAP.108 | Empty slabs and wedges near the core can't be selected, and the side buttons block clicks (bug) | MAP.65 |  |
| MAP.109 | Zooming in and out loads slowly (bug) | MAP.102 |  |
| MAP.110 | Slab button lines come out of numerical order (bug) |  |  |
| MAP.111 | "Generated only" should be "Charted only" and dim the stars too (bug) |  |  |
| MAP.112 | Nothing can be selected while "Generated only" is on (bug) | MAP.65, MAP.111 |  |
| MAP.113 | A nebula covering the whole sector can't be unselected, and nebulae need a show/hide toggle (bug) | MAP.79 | The toggle is one of MAP.79's per-kind buttons. |
| MAP.116 | Dense generated sectors crowd the Galaxy Map when zoomed out (bug) | MAP.102 | Also the per-zoom rule for which objects show; folds MAP.115 (Boss 2026-10-07 16:26Z). |
| MAP.127 | See into and step to the neighbouring blocks and slabs on the Galaxy Map | MAP.102, MAP.65 | Boss's third Galaxy Map problem (2026-10-07 16:26Z). |
| NAV.46 | The NAV picker can't click galaxy wedges to zoom in (bug) | NAV.15 | Fixed by pick mode on the shared layer (NAV.15). |

### Groundwork: nebula shapes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.75 | A nebula shape from metaballs and warped noise, as a mesh |  | Groundwork for the nebula map bugs; uses scikit-image for marching cubes. |
| MAP.103 | Nebulae don't show on the Galaxy Map or any other map (bug) | GEN.75 |  |
| MAP.104 | Nebula shading is missing on unfilled sectors (bug) | MAP.103 |  |
| MAP.105 | The nebula view should show its whole shape with the dimmed galaxy around it (bug) | GEN.75 | Replaces the diagram UX.38 is about. |
| UX.38 | The nebula and remnant diagrams' "-" button does nothing at the 1 ly limit (bug) | MAP.105 | Folded into MAP.105: the 3D nebula view replaces the diagram. Split out of UX.21 (its known dead control): lib/phenomenonmap.py lines 53 and 127, the mapzoom.js clamp; strict xfail in test_web_browser_maps.py. |

### Groundwork: UX sweep

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.37 | A UX sweep: remove redundant and duplicate controls so the interface gets out of the way | UX.28, MAP.95, NAV.41, UX.26, UX.27, UX.31, UX.40, ADM.34 | Moved into phase 0: UX.21 (a bug) needs it. Boss 04:42Z. Audit first (list of what to remove or merge), Boss reviews, then removals; after the controls it audits settle. |
| UX.21 | Clean up the web interface: overlapping buttons and dead controls (bug) | MAP.68, UX.26, UX.27, UX.31, NAV.32, UX.37 | Moved into phase 0 (all bugs in phase 0), last. Bug, but a final pass over the finished pages. Judgment: its one known dead control (nebula '-' at the 1 ly limit) could be split out into phase 0. The nebula "-" control is split out to phase 0 as UX.38. |

One thread at a time (Boss, 2026-10-03 and 2026-10-07). The
bugfix lane is the threads named "Bugfixes: ...", the groundwork
lane those named "Groundwork: ...". Default order: Bugfixes: CI red
first, then Groundwork: layout and libraries (so code moves once),
then the lanes alternate group by group. GEN.65 is held until Boss
sends the error text.

The Sector and system pages lane is paused: UX.24 and UX.29 were
committed only in its container (not pushed) and UX.25 is half
built, so that work may be lost. UX.28's icon set is approved.
MAP.101 is done (PR #431). OPS.23 is done (PR #435, the package layout plan). OPS.21 is done (PR #437, pins and Redis). The CI red group is done (PR #442: GEN.76, OPS.25, DB.12, PERF.27, TEST.80, TEST.81, TEST.87, MAP.114). GEN.82, UX.34, OPS.9, GEN.80 and ADM.27 are done (PR #448). PERF.28 is done (PR #454). MAP.115 is folded into MAP.116 (2026-10-07).

## Open questions for Boss

None.
