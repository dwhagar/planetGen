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

Two lanes (Boss 2026-10-03 05:38Z: "the most fundamental and needed changes first in phase 0 along with the bug fixes in 2 lanes, bugfixes and groundwork"). Bugfixes: the CI failures first, then generation, console and progress, maps and pages, prevalence (GEN.48), and ops and test flakes. Groundwork: the package layout and the move to third-party libraries with Redis (Boss's explicit directive), the data model (SQLAlchemy and Alembic, values in columns not JSON, one point-in-space object, Pydantic, scipy and astropy), names from IDs, the RQ queue with streamed logs and progress, Shoelace and TanStack components, the shared map engine, nebula shapes, and the UX sweep. Bugs that the groundwork fixes are folded into it and listed under it. GEN.65 is done (PR #476); GEN.116 is closed (Boss 2026-10-08 17:46Z).

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Bugfixes: CI red

Done: all eight items landed in PR #442 (2026-10-07).

### Bugfixes: generation

| ID | Item | Needs | Note |
|---|---|---|---|

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

### Groundwork: layout and libraries

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.39 | Markdown rendered by the markdown library |  |  |
| OPS.20 | Move the code base from zero dependencies to third-party open-source libraries | UX.39, DB.11, ADM.21, GEN.66 | Parent of the library migration; closes with its subitems. Boss's explicit directive. |

### Groundwork: data model

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.11 | The database layer and migrations on SQLAlchemy and Alembic |  | Every later schema change (DB.7, NAV.10, API.11, DB.13, GEN.106) is an Alembic migration. |
| DB.13 | Every stored value in its own column, not in JSON blocks, and indexed for search | DB.11 | Boss (2026-10-07 11:47Z): "This is to be considered a Phaser 0 priority." |
| ADM.21 | Input validation on Pydantic models |  | Also the error objects API recipes (API.18) return. |
| GEN.66 | Physics on scipy, and astropy constants and units |  |  |

### Groundwork: names from IDs

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.70 | A naming key in the control database, made at galaxy creation and changeable by admin |  |  |
| GEN.71 | Name interstellar objects, phenomena and constellations from the codec | GEN.70 | Folds in GEN.73 (nebulae unique names). |
| GEN.67 | Names from IDs for objects that have no star-derived name | GEN.70, GEN.71 | Parent; folds the two naming bugs and drops GEN.63. Boss 2026-10-08: stars, sectors, planets, moons and belts keep their names; the codec names only the objects with no star-derived name. |

### Groundwork: queue, logs and caches

| ID | Item | Needs | Note |
|---|---|---|---|

### Groundwork: web components

| ID | Item | Needs | Note |
|---|---|---|---|

### Groundwork: map engine

| ID | Item | Needs | Note |
|---|---|---|---|

### Groundwork: nebula shapes

| ID | Item | Needs | Note |
|---|---|---|---|

### Groundwork: UX sweep

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.21 | Clean up the web interface: overlapping buttons and dead controls (bug) |  | Moved into phase 0 (all bugs in phase 0), last. Bug, but a final pass over the finished pages. Judgment: its one known dead control (nebula '-' at the 1 ly limit) could be split out into phase 0. The nebula "-" control is split out to phase 0 as UX.38. |

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
MAP.101 is done (PR #431). OPS.23 is done (PR #435, the package layout plan). OPS.21 is done (PR #437, pins and Redis). The CI red group is done (PR #442: GEN.76, OPS.25, DB.12, PERF.27, TEST.80, TEST.81, TEST.87, MAP.114). GEN.82, UX.34, OPS.9, GEN.80 and ADM.27 are done (PR #448). PERF.28 is done (PR #454). GEN.77, GEN.78 and GEN.79 are done (PR #461): giants now end with a short bright stretch (2% of the giant phase, 1000 to 2500 Lsun), so bright stars appear in the old disk and bulge and the bright-star count and its database space roughly double at the default threshold; sparse cells get a halo floor of a thousandth of the local density instead of being skipped. OPS.6, UX.36, MAP.117, NAV.41, OPS.7, SEC.31, UX.44 and MAP.118 are done (PR #457). GEN.81 and OPS.26 are done (PR #467). PERF.26 is done (PR #470). OPS.24 and its parent OPS.22 are done (PR #473): the package move is finished. OPS.27 is done (PR #475: Redis in WSL2 only). ADM.33 and GEN.65 are done (PR #476: owner override, and 18 edge-of-galaxy neighbourhood tests through web and CLI; GEN.116 keeps watch for Boss's error). PERF.25 is done (PR #478: the page cache on cachetools; tilecache.py stays). The ops-and-flakes bugs TEST.71, TEST.78, TEST.82, TEST.84, TEST.86, TEST.88, TEST.89 and TEST.90 are done (PR #484, which also added GEN.117's generator z-spread test). UX.24 and UX.29 are done (PR #487). UX.28 and UX.25 are done (PR #490). TEST.91 is done (PR #493). PERF.19 is done (PR #492, the work-queue audit; Boss approved queueing single-object generation: short generation goes on the queue and the page waits on it). Since PR #496 (PERF.24 steps 1-2), generation with more than one worker needs Redis and falls back to one worker without it; tests that monkeypatch generation pin PLANETGEN_WORKERS=1 or use tests/worker_patches. GEN.118 and GEN.119 are done (PR #491: the Milky Way bulge and a thick disk). PR #491 needs a galaxy reset and regenerate (update.sh, python -m planetgen.cli.reset, planetgen plan, planetgen galaxy); the plan now reaches 4.1 kpc above the plane (was 1.3 kpc) with about 21 billion qualifying sectors (was 10 billion). GEN.117 stays open with MAP.116. GEN.48, GEN.52 and TEST.75 are done (PR #502): sector and galaxy runs take --prevalence FEATURE=PERCENT instead of forcing (for example comets=+50, binary_system=-100); habitable worlds, belts, large stars and planets use shares measured over 20,000 systems; intelligent_life prevalence scales the population pass's CIVILIZATION_CHANCE; -large_star on a single system now really excludes a large star (it was a no-op); schema v54 adds system_configs.prevalence_* columns with a migration; the local full suite now needs redis-server running (16 work-queue tests fail without it). ADM.16 (the Generate page fields) is next. ADM.16 is done (PR #506: prevalence fields on the New galaxy and Generate sectors forms; they will follow the field() macro when it moves to Shoelace with UX.40). ADM.23 is done (PR #508: log lines printed under a progress bar are no longer hard-wrapped at 80 columns, and [brackets] print as written instead of as rich markup). ADM.37 is done (PR #514: the prevalence fields start at each feature's usual share, labelled with what each is a share of). PERF.24 is done (PRs #496, #498, #501, #504, #510, #512 and #517): generation and web jobs run on RQ. System and sector wiki uploads wait up to 8 s, then answer 201, or 202 with a job id if slower; the wiki token never goes through Redis (the worker reads it from the environment and config.json); without Redis the Generate-page estimate and the one-off system page run in the web process as before. OPS.19 and ADM.22 no longer wait on it. PRs #501 and #504 (PERF.24 steps 4a and part of 4b) queue neighbourhood generation and sector regenerate: the API answers 202 with a job id, and GET /api/jobs/<id> is the job-state endpoint for queued API work. Short generation and wiki uploads remain; Boss (2026-10-08 00:14Z) put wiki uploads on the queue too, for later batch uploads. MAP.115 is folded into MAP.116 (2026-10-07). MAP.110 is done (PR #520: the Galaxy Map slab buttons are always in slab-number order, descending from above and ascending under the plane, with the lines kept uncrossed in the side layout; on a phone the buttons are in number order but the leader-line ends still use the old line-end logic, so file a bug only if Boss sees crossed lines there). TEST.92 and TEST.93 are done (PR #522: GET /api/databases counts over one connection instead of leaving a pool open per schema, which caused MariaDB's "Too many connections" under load; get_job no longer reports a just-started job as "interrupted"). OPS.19 is done (PR #525: the jobs folder defaults to /var/lib/planetGen/jobs; create-cache-dir.sh and update.sh move old jobs and remove the old lowercase folder). UX.2 is done (PR #528, the first UX.40 PR: the Shoelace component set is vendored and the header Menu and gear are sl-dropdowns). UX.40 stays open; Foundations lane 1 continues it with UX.26, UX.31, UX.27 and ADM.14. MAP.65 is done (PR #521, the shared map picking, hover and info-panel layer; it frees MAP.79, MAP.107, MAP.108, MAP.119, MAP.120 and MAP.122 to start, and MAP.112, MAP.123, MAP.61, UX.48 and NAV.45 lose this one prerequisite). MAP.111 is done (PR #531: Charted only dims uncharted stars and outlines charted blocks; MAP.112 stays open). MAP.107's cause: a click while the just-picked stage is still loading lands on the new stage from the old view (fix in galaxystageview takePick); MAP.108 part 1 (buttons covering the map) could not be reproduced and gets a regression test, part 2 (picking empty slabs and wedges) is a separate cause. UX.31 is done (PR #533, UX.40's second piece: the edit actions on the system, sector and phenomenon pages are one Shoelace menu button with dialog confirms). UX.26 stays open for the sector page's Admin panel (generate neighborhood, wiki upload), still inline forms. MAP.107 is done (PR #537: choices are taken on the displayed stage, not the one still loading). MAP.108 part 1 (buttons covering the map) is done by a regression test, as it could not be reproduced; MAP.108 stays open for part 2 (picking empty slabs and wedges). UX.27 is done (PR #539, UX.40's third piece: the system page's buttons sit on one row, with a Navigate menu where they don't fit). UX.40 stays open (ADM.14 next, then the Generate page and form remainder). MAP.112 is done (PR #541: Charted only no longer blocks picking; only choosing a NAV end still needs something generated). MAP.108 part 2 may be fixed by the same change; it stays open until the lane confirms. TEST.94 to TEST.101 are done (PR #549, the load-only test failures). MAP.108 is closed as not reproducible (both parts), covered by tests: #537 added a regression test for part 1, and #541 with MAP.65 (#521) plus a browser test that hovers and picks an empty arc cover part 2; reopen only if Boss sees it again. ADM.22 is done (PR #551: job output streams over SSE into an Xterm.js terminal with native progress bars); ADM.24, ADM.25 and ADM.26 are unblocked. GEN.120 is done (PR #554: the phoneme codec lives in src/planetgen/names/gated_phoneme_codec.py with tests; decode(phrase, domain, length=19) is exact, so GEN.67 always passes 19). TEST.101 was fixed again in #554 (another test in the same worker left Ctrl+C ignored); it stays retired. DB.14, MAP.129 and MAP.130 are done (PR #556, sector coloring: raw sector stats in galaxy schema v55, the Galaxy Map and its blocks drawn as a translucent fill by star age, density and luminosity, the average-star-color rule removed; the stage API returns `stats` instead of `look`; Boss must run sudo ./update.sh). MAP.128 stays open, narrowed to the Sector Map half. UX.26 is done (PR #558: the sector page's Admin panel is an Admin menu with Upload to Wiki and Generate neighborhood dialogs); ADM.34 is unblocked. ADM.24, ADM.25 and ADM.26 are done (PR #560: a failed job keeps its log open until Continue, full tracebacks reach the console, the job log and the debug log with a Copy log button, and the backfill and band top-up show progress bars); the report that the post-galaxy backfill shows no bar at all was not reproducible, so reopen only if Boss names a specific Generate action. ADM.34 and MAP.95 are done (PR #562: one Admin menu per page holding only that page's actions, and a Current button beside Forward on the Galaxy Map; the Galaxy Map's page-level admin stays with ADM.35); UX.37 has no prerequisites left. NAV.13 is done (PR #564: static/picker.js holds the selection with select, up, into, sideways, onChange and trail; the Galaxy Map's Up, Reset, Esc, Backspace and breadcrumb go through it; picker.trail() gives NAV.14 the parent chain). NAV.15 is split: its body-level part became NAV.50 (needs NAV.16). NAV.15 is done (PR #566: pick mode with a Cancel banner on every display, bookmarks that keep the pick, and one Use as start/destination button on the system and phenomenon pages); NAV.16 already lists NAV.7 as its prerequisite. NAV.46 is closed as already fixed (PR #568: wedge clicking in NAV pick mode has worked since the shared picker landed, so the PR adds a browser test that clicks through every stage down to a sector). MAP.66 is done (PR #548: a generated sector opens in place on the Galaxy Map as the drill-down's last stage); the sector page embedding the same engine is done in MAP.68 (PR #573). NAV.14 is done (PR #571: one-line breadcrumb at any width with a … overflow menu, shared by every page and the Galaxy Map); MAP.106 now waits on nothing else. GEN.75 moved to Foundations lane 1, and Boss approved GEN.68's unique-ID plan on 2026-10-08 at 13:33Z (comparison at /mnt/project-files/notes/gen68-object-ids.md), so GEN.69 is free and GEN.72 and GEN.73 follow it. MAP.128 is closed by Boss's decision (2026-10-08 13:35Z): sector space is not tinted at all, so the Sector Map half is dropped and only the Galaxy Map tint (PR #556) stays; MAP.131 is now the Galaxy Map only, and MAP.131 and MAP.132 no longer wait on MAP.128. The UX audit (UX.37) is started as an audit-only thread; Boss approved all 25 findings (2026-10-08 17:43Z) and they are filed as UX.50 to UX.74 in phase 0. MAP.116, MAP.115 (folded into it) and GEN.117 are done (PRs #598, #600 and #601: bright stars sampled off the plane with a v57 index, strict per-zoom star budgets shared out by sector, and comets, rogue planets and asteroid fields kept off the Galaxy Map, with a test). MAP.102 is done (PR #574 camera-relative stars, with the tiles, level of detail and caches already on main; the BVH clause was stale because stars are not pickable). MAP.109 is done (PR #603: a stated 70,000-star cap per Galaxy Map view with a guard test). UX.41 and UX.33 are done (PRs #594 to #606: every list page, admin list and search panel is a TanStack table with virtual scrolling, sorting and faceted filters, including the phenomena filters). MAP.68 is done (PR #573: the sector page's map is the Galaxy Map engine locked to its sector; sectormap.js and render_map_panel are gone and the sector scale line reads in sectors and parsecs). MAP.67 is narrowed to the system and body URL forms, which ride with MAP.62, so MAP.106 and MAP.124 no longer wait on it. MAP.79 is done (PR #577: per-kind show/hide buttons in the map Menu, kept in the URL as hide=). MAP.113 is narrowed to unselecting a nebula that covers the sector, and MAP.123 stays open for star types, highlight, the luminosity slider and the Galaxy Map. GEN.75 is done (PR #576 the shape module, #579 the stored shape and the /api/nebulae/<id>/shape mesh at schema v56, #580 containment uses the shape; remnants stay spheres, and containment already stored stays sphere-based until the affected sectors are regenerated). MAP.103 and MAP.105 are now free; MAP.104 waits on MAP.103 and UX.38 on MAP.105. MAP.113 is done (PR #584: Escape and a click on the selected object clear the selection; the toggle half was MAP.79). MAP.106 is done (PR #586: a browser test walks the Galaxy Map to a sector and back by every route and checks the breadcrumb). MAP.105 and UX.38 are done (PR #590: the nebula page shows the nebula in 3D with the galaxy's brightest stars dimmed around it, and the supernova remnant diagram's "-" button works). MAP.103 stays open for the System Map wash. MAP.103 and MAP.104 are done (PR #592: the System Map tints a system inside a nebula, and nebulae show over unfilled sectors, with a browser test); nebulae are now drawn from their shape on every map and page. UX.40 is done (PR #544, its fourth piece: the Generate and one-off system pages' text boxes line up, ADM.14, with a browser test); its remainder became UX.49 (form fields as Shoelace components), and UX.26 stays open for the sector page's Admin panel. GEN.116 is closed (Boss 2026-10-08 17:46Z: not seen since the one report; reopen only if he reports it again). TEST.102 is done (PR #610: the flaky species test planted its leftover species on a planet that was sometimes a homeworld). SEC.29 and TEST.72 are done (PR #612: two-step sign-in runs on pyotp and segno, and the hand-rolled totp.py and qrcodegen.py are deleted). NAV.29 and NAV.33 are done (PR #615: Start Here and End Here work in place on both maps, and the scene JSON no longer carries nav links); NAV.32 is unblocked. UX.63 is done (PR #618: one action bar on the system, phenomenon and sector pages). UX.64 is done (PR #620: navigation buttons are outlined, Wikitext and Markdown are outlined toggles, and bookmarks are one star toggle). UX.66 and UX.74 are done (PR #617: the sector edge wording and the Generate page's idle Current job card). GEN.68 and GEN.69 are done (PR #623: every object has a unique uid column with a unique key on 15 tables, schema v58; Boss must run sudo ./update.sh); GEN.72 and GEN.73 are free. TEST.103 is done (PR #625: the real defect was the Galaxy Menu panel opening past the right edge at 390 px, fixed in CSS; the sector 600 and 820 px cases did not reproduce, so the layout test now waits for the page to settle and names any lasting overlap). ADM.38 is done (PR #627: the worker's traceback rides on the exception's __notes__ so the Python 3.9 leg should pass; the next 3.9 CI run confirms). MAP.61 is done (PR #629: the Sector Map's wedge fills are removed, outlines only; MAP.67's system and body URL forms still ride with MAP.62). MAP.62 is six subitems (MAP.69 to MAP.74) and MAP.70 depended on GEN.74, which is now built (PRs #642, #645, #662 and #665), and Boss has given the go: the map engine lane builds MAP.62 first, then MAP.125. UX.72, UX.56 and UX.73 are done (PR #632: the gear menu holds Theme, Account, Admin and Logout, the Admin hub and a tab row are on the admin pages, and the sector wiki link form moved into the sector page's Admin menu). UX.65 is done (PR #635: the sector header chips hold facts only, the quadrant is a plain link, and the comet estimate is a Details line under the map). UX.57, UX.58, UX.60 and UX.62 are done (PR #637), and UX.50 is partly done (the Map help dialog is on the Galaxy and Sector maps; the System, NAV and diagram maps are left). UX.50, UX.59 and UX.61 are done (PR #640: Map help on the System, NAV and diagram maps, one breadcrumb trail starting at Home, and Measure distance under the System Map), which clears UX.50 and UX.57 to UX.62. NAV.32, TEST.104, SEC.30 and TEST.83 are done (PR #631: Sector Map routes cross sectors; PR #634: the TEST.104 flake; PR #643: request limits and login lockouts count in Redis, which needs `sudo ./update.sh`, and the test causes behind TEST.83). UX.68 is done (PR #646: the system Admin menu in the action bar, per-row Admin menus for planets, moons and belts). UX.67, UX.70, UX.71 and TEST.105 are done (PR #648: Wikitext and Markdown buttons, the NAV result page and landing page, and the Admin hub links' contrast). TEST.106 is done (PR #650: the layout test left the maps' Menu, Bookmarks and history popovers closed; no product change). UX.52 and UX.53 are done (PR #654: the Contents "Nearest" column with three neighbours, the system page's "Nearest" line, and a single star appearing once in the Stars table with no Role column). GEN.72 is done (PR #657: a backfilled bright star gets its word-salad name only when its sector is generated, and keeps its position ID). UX.69 is done (PR #652 and PR #659: scroll cue and phone wrapping for the sideways-scrolling tables, and secondary columns stacked under the name on a phone). UX.51 is done (PR #663: cards no longer repeat the page title). GEN.74 is done (PRs #642, #645, #662 and #665: one point-in-space object that keeps every coordinate system in step), which meets MAP.70's prerequisite. MAP.69 is done (PR #668: GET /api/systems/<id>/scene, a system scene endpoint with 3D orbits). TEST.107 was closed without a fix: test_controls_do_not_overlap passes on current main. UX.54 and UX.55 are done (PR #671: search no longer shows empty groups or echoes the query, and Home and Systems no longer repeat other pages' tables). That was the last of the 25 UX audit items, so UX.37 is done too, and UX.21 is unblocked. UX.75 is done (PR #679: the sector and phenomenon pages' Admin menus moved into the shared action bar). MAP.125 is done (PRs #678 and #683: one 3D interface from the galaxy down to a moon). MAP.62 and its subitems MAP.67 and MAP.70 to MAP.74 are done (PR #670: position modules and the clock; PR #673: the 3D system view with a free camera, scale modes and the 3D view on the system page); Boss answered the MAP.74 question (2026-10-08): system pages open in 3D by default, remembering the last view used; the map lane makes that change in a small follow-up PR. MAP.126 is done (PR #686: orbital trajectories of selected objects in their frame of reference). NAV.7 is done (PR #688: one reference for every object kind with its parent chain). TEST.108 and TEST.109 are done (PR #690: the two tests that failed on main). MAP.131 is done (PR #692: the Color by switch on the Galaxy Map, with a legend).

## Open questions for Boss

None. (GEN.67, Boss 2026-10-08: stars, sectors, planets, moons and belts keep their names; the codec names only the objects with no star-derived name.)
