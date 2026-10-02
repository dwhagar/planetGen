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

Bug fixes first (Boss 04:45Z: phase 0 is primarily bug fixes and the groundwork that goes with them; 04:55Z: "Bugs first"): the binary-pair and forcing bugs, the Galaxy Map follow-ups and drill-down fixes, the System Map, the generation bugs (with class S for GEN.38), sector stats and colors (with the per-sector stats table MAP.86 needs), the routing groundwork, the sector and system page bugs, the small page bugs, and ops and test flakes last. The parallel path lane is done (GEN.39, DB.6, OPS.10: PRs #381, #387, #391), and so are the binary pairs lane (GEN.62, GEN.51, TEST.85: PRs #393, #398, #403) the System Map lane (MAP.57, MAP.88, MAP.92: PR #405), the generation bugs lane (GEN.60, GEN.38, GEN.47: PRs #415, #419), and the Galaxy Map drill-down lane apart from the slab-button fixes MAP.98, MAP.100 and MAP.99 (PRs #408, #410, #413). GEN.65, a web-only generation error, is high priority but held until Boss says to start.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Web generation error (high priority)

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.65 | A generation run fails from the web UI but not from the CLI (bug) |  | Boss 08:08Z: high priority, top of phase 0, not started yet. Details unknown; ask Boss for the error. |

### Galaxy Map drill-down

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.98 | Slab button lines should end at the nearest edge of their slab (bug) |  | Boss 08:08Z. slabAnchor in galaxystageview.js ends the line inside the slab; end it on the prism edge nearest the button, at that edge's point nearest the button (Boss 08:35Z). |
| MAP.100 | Slab button labels on one line: "#N" and how much is charted (bug) |  | Boss 08:17Z. "#4 Unknown", "#2 < 0.01 % charted", "#6 ≈ 2.43% charted"; drop "generated" and x / total. |
| MAP.99 | Slab buttons that don't fit the window split across both sides of the map, shrink, or give way to map picking (bug) | MAP.100 | Boss 08:17Z. Two columns, one per side; smaller buttons on small screens; none at all if still too many. Lines still end per MAP.98; replaces the one column under 600 px. |

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
| UX.33 | Filter phenomena by their classes and types (bug) |  | Moved to phase 0: filters are built from the class lists, so new classes (GEN.28) appear on their own; GEN.47 is done (PR #419), so nebulae exist and nothing blocks it. Filters over class lists that GEN.28 changes. |
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
| TEST.87 | The two-process id-block test times out under full parallel load (bug) |  | test_db_id_blocks_edges.py; _queue.Empty under -n auto, passes alone 3/3. Same path as TEST.81. |

GEN.65 (web generation error) is high priority and comes first, but
is held: no thread starts on it until Boss says so (08:08Z).

At most two build threads run at once (Boss 02:51Z). Done lanes:
Parallel path (PRs #381, #387, #391), Galaxy Map follow-ups (PRs
#395, #399), Binary pairs (PRs #393, #398, #403), System Map (PR
#405) and Generation bugs (PRs #415, #419). Galaxy Map drill-down is
done apart from MAP.98, MAP.100 and MAP.99 (PRs #408, #410, #413).
Sector stats and colors is running. Then the lanes
start in the order above as a slot frees: MAP.98, MAP.100 and MAP.99 (slab buttons), Routing groundwork, Sector and
system pages (text labels with an icon hook if UX.28's icon list is
not approved yet), Small page bugs, and Ops and flakes last. NAV.7
and DB.8 open phase 1.

## Open questions for Boss

- PERF.11: Store each sector's expected and actual density, see its entry in TODO.md.
