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

Two lanes (Boss 2026-10-03 05:38Z: "the most fundamental and needed changes first in phase 0 along with the bug fixes in 2 lanes, bugfixes and groundwork"). Bugfixes: the CI failures first, then generation, console and progress, maps and pages, prevalence (GEN.48), and ops and test flakes. Groundwork: the package layout and the move to third-party libraries with Redis (Boss's explicit directive), the data model (SQLAlchemy and Alembic, values in columns not JSON, one point-in-space object, Pydantic, scipy and astropy), names from IDs, the RQ queue with streamed logs and progress, Shoelace and TanStack components, the shared map engine, nebula shapes, and the UX sweep. Bugs that the groundwork fixes are folded into it and listed under it. GEN.65 is done (PR #476); GEN.116 keeps watch for Boss's error text.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Bugfixes: CI red

Done: all eight items landed in PR #442 (2026-10-07).

### Bugfixes: generation

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.116 | Error when generating a neighbourhood centred near the galaxy edge, awaiting Boss's error text (bug) |  | Back burner (Boss 17:11Z): GEN.65's 18 edge tests pass (PR #476); waits for his error text. |

### Bugfixes: console and progress

| ID | Item | Needs | Note |
|---|---|---|---|

### Bugfixes: maps and pages

| ID | Item | Needs | Note |
|---|---|---|---|

### Bugfixes: prevalence

| ID | Item | Needs | Note |
|---|---|---|---|

### Bugfixes: ops and flakes

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.92 | The web job runner's first-failure test reports the job as interrupted under full-suite load (bug) |  | test_web_generate.py; seen once in PR #496's local run, passes alone. |
| TEST.93 | Timing tests fail and MariaDB drops connections under full-suite load (bug) |  | 10 load-only failures and MariaDB connection drops in a 13,587-test run (PRs #501, #504); rate-limit cases are TEST.83. |

### Groundwork: layout and libraries

| ID | Item | Needs | Note |
|---|---|---|---|
| SEC.29 | Two-step sign-in on pyotp, QR codes on segno |  | Folds TEST.72 (the 2FA flake). |
| TEST.72 | Intermittent failure in the two-step (2FA) sign-in test (bug) | SEC.29 | Folded into SEC.29 (pyotp). |
| SEC.30 | Login and request rate limits on Flask-Limiter with Redis storage |  | Folds TEST.83 (rate-limit tests under load). |
| TEST.83 | Rate-limit tests fail under parallel load (bug) | SEC.30 | Folded into SEC.30 (Flask-Limiter). test_web_admin.py; passes alone, fails under -n auto. |
| UX.39 | Markdown rendered by the markdown library |  |  |
| OPS.20 | Move the code base from zero dependencies to third-party open-source libraries | SEC.29, UX.39, PERF.24, SEC.30, DB.11, ADM.21, GEN.66, UX.40, UX.41, ADM.22, MAP.102 | Parent of the library migration; closes with its subitems. Boss's explicit directive. |

### Groundwork: data model

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.11 | The database layer and migrations on SQLAlchemy and Alembic |  | Every later schema change (DB.7, NAV.10, API.11, DB.13, GEN.106) is an Alembic migration. |
| DB.13 | Every stored value in its own column, not in JSON blocks, and indexed for search | DB.11 | Boss (2026-10-07 11:47Z): "This is to be considered a Phaser 0 priority." |
| ADM.21 | Input validation on Pydantic models |  | Also the error objects API recipes (API.18) return. |
| GEN.66 | Physics on scipy, and astropy constants and units |  |  |
| GEN.74 | One point-in-space object that keeps every coordinate system in step, used by every object |  | Groundwork for the orbital work, Hill-sphere placement and light-travel positions. |

### Groundwork: names from IDs

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.68 | Research: the cheapest unique IDs for every object, unfilled sectors included |  | Asked on 2026-10-07; decides how GEN.69 builds IDs. |
| GEN.69 | A unique ID for every object, star systems and unfilled sectors included | GEN.68 |  |
| GEN.70 | A naming key in the control database, made at galaxy creation and changeable by admin |  |  |
| GEN.71 | Remove the word-salad name code, nltk and the name registries | GEN.69, GEN.70 |  |
| GEN.72 | A backfilled bright star should get a name only when its sector is generated (bug) | GEN.69 |  |
| GEN.73 | Nebulae don't get unique names (bug) | GEN.69 |  |
| GEN.67 | Names from IDs: replace word-salad name generation | GEN.68, GEN.69, GEN.70, GEN.71, GEN.72, GEN.73 | Parent; folds the two naming bugs and drops GEN.63. |

### Groundwork: queue, logs and caches

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.24 | The work queue and web jobs on Redis with RQ |  | Replaces workQueue.py and jobRunner.py; folds OPS.19. PERF.19's audit is done (PR #492): single-object generation is queued and waited on. Steps 1-2 merged (PR #496): more than one worker needs Redis, else one worker. |
| OPS.19 | The Generate jobs folder is /var/lib/planetgen while the checkout is /var/lib/planetGen (bug) | PERF.24 | Folded into PERF.24: the job store moves with the queue. jobs.py and deploy-paths.py defaults; update.sh moves an old lowercase jobs folder. Same update.sh as OPS.7/OPS.8. |
| ADM.22 | Job logs streamed over SSE into Xterm.js, with native progress bars | PERF.24 | Folds the four log and progress bugs below. |
| ADM.23 | Log output wraps with hard line breaks (bug) | ADM.22 |  |
| ADM.24 | A failed action's log closes before it can be read (bug) | ADM.22 |  |
| ADM.25 | Error tracebacks don't reach the console and the web log window (bug) | ADM.22 |  |
| ADM.26 | The bright-star backfill shows no progress bar on the web (bug) | ADM.22 |  |

### Groundwork: web components

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.40 | Buttons, menus and dialogs from Shoelace web components |  | Folds UX.2, UX.26, UX.27, UX.31 and ADM.14; uses UX.28's approved icons. |
| UX.26 | Edit and admin actions as a button that opens a menu (bug) | UX.40 | Folded into UX.40; first screen of ADM.34. Sector page admin panel and edit_controls.html. |
| UX.31 | Editing a star system: an edit button with a quick menu, not a long panel (bug) | UX.26 | Folded into UX.40. system.html edit panel (_edit_rows in system_pages.py). |
| UX.27 | System page: the system and navigation buttons on one row that doesn't overlap (bug) | UX.40 | Folded into UX.40. system.html subhead; shares wording with NAV.29. |
| UX.2 | Menus sized to what they hold (bug) | UX.40 | Folded into UX.40. style.css menus; independent. |
| ADM.14 | Line up the Generate page's text boxes, not their headings (bug) | UX.40 | Folded into UX.40. generate.html field layout; before GEN.24 adds fields (ADM.16 is done, PR #506) to the same page. |
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
| MAP.110 | Slab button lines come out of numerical order (bug) |  | PR #491 fixed the crossing lines only; buttons still follow on-screen order. |
| MAP.111 | "Generated only" should be "Charted only" and dim the stars too (bug) |  |  |
| MAP.112 | Nothing can be selected while "Generated only" is on (bug) | MAP.65, MAP.111 |  |
| MAP.113 | A nebula covering the whole sector can't be unselected, and nebulae need a show/hide toggle (bug) | MAP.79 | The toggle is one of MAP.79's per-kind buttons. |
| MAP.116 | Dense generated sectors crowd the Galaxy Map when zoomed out (bug) | MAP.102 | Also the per-zoom rule for which objects show; folds MAP.115 (Boss 2026-10-07 16:26Z). Must also fix GEN.117 (bright stars kept to the plane by brightest-first picking). |
| GEN.117 | The Galaxy Map shows bright stars only in a thin band on the galactic plane (bug) | MAP.116 | Map picking, not generation (lane 1, 18:47Z): 400 brightest per tile are all young plane stars. Fixed with MAP.116. |
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
| UX.37 | A UX sweep: remove redundant and duplicate controls so the interface gets out of the way | MAP.95, UX.26, UX.27, UX.31, UX.40, ADM.34 | Moved into phase 0: UX.21 (a bug) needs it. Boss 04:42Z. Audit first (list of what to remove or merge), Boss reviews, then removals; after the controls it audits settle. |
| UX.21 | Clean up the web interface: overlapping buttons and dead controls (bug) | MAP.68, UX.26, UX.27, UX.31, NAV.32, UX.37 | Moved into phase 0 (all bugs in phase 0), last. Bug, but a final pass over the finished pages. Judgment: its one known dead control (nebula '-' at the 1 ly limit) could be split out into phase 0. The nebula "-" control is split out to phase 0 as UX.38. |

One thread at a time (Boss, 2026-10-03 and 2026-10-07). The
bugfix lane is the threads named "Bugfixes: ...", the groundwork
lane those named "Groundwork: ...". Default order: Bugfixes: CI red
first, then Groundwork: layout and libraries (so code moves once),
then the lanes alternate group by group. Boss (2026-10-07 17:11Z):
"When one thread is idle waiting for CI we can start another thread on
something else, the idea being to always have at least 1 thread going
without even going tover the 2 dev 1 todo limit." GEN.65 is done (PR #476);
GEN.116 keeps watch for Boss's error text.

The Sector and system pages lane's items are all done: UX.24 and
UX.29 (PR #487), UX.28's icons and UX.25 (PR #490), rebuilt by
Bugfixes lane 1.
MAP.101 is done (PR #431). OPS.23 is done (PR #435, the package layout plan). OPS.21 is done (PR #437, pins and Redis). The CI red group is done (PR #442: GEN.76, OPS.25, DB.12, PERF.27, TEST.80, TEST.81, TEST.87, MAP.114). GEN.82, UX.34, OPS.9, GEN.80 and ADM.27 are done (PR #448). PERF.28 is done (PR #454). GEN.77, GEN.78 and GEN.79 are done (PR #461): giants now end with a short bright stretch (2% of the giant phase, 1000 to 2500 Lsun), so bright stars appear in the old disk and bulge and the bright-star count and its database space roughly double at the default threshold; sparse cells get a halo floor of a thousandth of the local density instead of being skipped. OPS.6, UX.36, MAP.117, NAV.41, OPS.7, SEC.31, UX.44 and MAP.118 are done (PR #457). GEN.81 and OPS.26 are done (PR #467). PERF.26 is done (PR #470). OPS.24 and its parent OPS.22 are done (PR #473): the package move is finished. OPS.27 is done (PR #475: Redis in WSL2 only). ADM.33 and GEN.65 are done (PR #476: owner override, and 18 edge-of-galaxy neighbourhood tests through web and CLI; GEN.116 keeps watch for Boss's error). PERF.25 is done (PR #478: the page cache on cachetools; tilecache.py stays). The ops-and-flakes bugs TEST.71, TEST.78, TEST.82, TEST.84, TEST.86, TEST.88, TEST.89 and TEST.90 are done (PR #484, which also added GEN.117's generator z-spread test). UX.24 and UX.29 are done (PR #487). UX.28 and UX.25 are done (PR #490). TEST.91 is done (PR #493). PERF.19 is done (PR #492, the work-queue audit; Boss approved queueing single-object generation: short generation goes on the queue and the page waits on it). Since PR #496 (PERF.24 steps 1-2), generation with more than one worker needs Redis and falls back to one worker without it; tests that monkeypatch generation pin PLANETGEN_WORKERS=1 or use tests/worker_patches. GEN.118 and GEN.119 are done (PR #491: the Milky Way bulge and a thick disk). PR #491 needs a galaxy reset and regenerate (update.sh, python -m planetgen.cli.reset, planetgen plan, planetgen galaxy); the plan now reaches 4.1 kpc above the plane (was 1.3 kpc) with about 21 billion qualifying sectors (was 10 billion). GEN.117 stays open with MAP.116. GEN.48, GEN.52 and TEST.75 are done (PR #502): sector and galaxy runs take --prevalence FEATURE=PERCENT instead of forcing (for example comets=+50, binary_system=-100); habitable worlds, belts, large stars and planets use shares measured over 20,000 systems; intelligent_life prevalence scales the population pass's CIVILIZATION_CHANCE; -large_star on a single system now really excludes a large star (it was a no-op); schema v54 adds system_configs.prevalence_* columns with a migration; the local full suite now needs redis-server running (16 work-queue tests fail without it). ADM.16 (the Generate page fields) is next. ADM.16 is done (PR #506: prevalence fields on the New galaxy and Generate sectors forms; they will follow the field() macro when it moves to Shoelace with UX.40). PRs #501 and #504 (PERF.24 steps 4a and part of 4b) queue neighbourhood generation and sector regenerate: the API answers 202 with a job id, and GET /api/jobs/<id> is the job-state endpoint for queued API work. Short generation and wiki uploads remain. MAP.115 is folded into MAP.116 (2026-10-07).

## Open questions for Boss

None.
