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
| TEST.80 | Intermittent failure in the admin change-star test (bug) |  |  |
| TEST.81 | Two processes reserving id blocks of one table can deadlock (bug) |  | _db._reserve_id_block, 1213 deadlock on MariaDB 10.11. |
| TEST.82 | Intermittent failure in the orbit-ceiling trim test (bug) |  | test_validation.py, physics area. |
| TEST.83 | The sign-in rate-limit test fails under parallel load (bug) |  | test_web_admin.py; passes alone, fails under -n auto. |
| TEST.84 | The every-column round-trip test depends on whether a quasar got placed (bug) |  | NULL_IN_THIS_GALAXY depends on the draw; about 1 in 5 on MySQL 8.0. |
| TEST.85 | A bright-star layer test once hit a name collision count of -1 (bug) |  | nameUniqueness.py:137; related to GEN.44 bands and GEN.57 names. |
| TEST.86 | Intermittent failure in the concurrent-insert recovery test (bug) |  | test_galaxy_gen.py; failed once under -n auto. |
| UX.34 | The sector summary calls white dwarfs "B-type" and "A-type" systems (bug) |  | sector_generation_summary_lines in generate.py; one PR with OPS.9. |
| OPS.9 | Multi-line messages lose their prefix in the debug log (bug) | UX.34 | Same summary record. |
| API.15 | Log every API call with its user, how it came in, and its HTTP response code |  | Needs no user accounts (Boss 01:31Z). |
| TEST.78 | A resume test's sector query fails under ONLY_FULL_GROUP_BY on MariaDB 10.11 (bug) |  | From the Database thread (PR #342). |
| OPS.19 | The Generate jobs folder is /var/lib/planetgen while the checkout is /var/lib/planetGen (bug) |  | jobs.py and deploy-paths.py defaults; update.sh moves an old lowercase jobs folder. Same update.sh as OPS.7/OPS.8. |

### Page groundwork and small page bugs

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.28 | Investigate icons instead of words on buttons |  | Survey and icon sprite; UX.25, UX.26 and UX.27 use its icons; MAP.55 (done, PR #369) left text labels with an icon hook. |
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

### Parallel path (top priority)

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.6 | Store the galaxy's 128-bit seed, the version that made it, and every generation run |  | Sub-item of GEN.39; fresh galaxy; 22-hex-digit version and environment key (Boss 02:08Z). |
| OPS.10 | The galaxy seed and version at the top of every generation log | DB.6 | Sub-item of GEN.39. |

### Database consistency check

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.8 | Check a galaxy database and say whether it is damaged |  | Boss 02:13Z: phase 0, its own thread. Read-only. Stats and version checks switch on once GEN.44/PERF.11 and DB.6/DB.7 land; no hard dependency. |

### Binary pairs and single-system forcing

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.62 | Binary stars: two-word names sharing the first word, the companion's word drawn from "small" and "child" sounds, planets named for one word (bug) |  | Boss 04:03Z: wide pairs "Blue Green"/"Blue Red", companion word from small/child sounds, planets Blue I / Red I. Replaces the 03:38Z A/B rule. |

### Galaxy map picker and arc

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.91 | Hovering the map while picking a slab highlights the whole slab, not one cube |  | Follow-up to PR #369 (Boss 03:52Z): applyHover in galaxystageview.js lights the whole slab; MAP.56, MAP.54 and MAP.77 keep it. |
| NAV.40 | Bookmarks can't be used to find the start or destination once a course pick has begun (bug) |  | Boss 04:12Z. bookmarks.js menu and NAV select keep the pick; nav_page.py, galaxymap3d.py, sector.html. |

The parallel path thread starts first (Boss: top priority). The map
groundwork thread runs to MAP.85 (the arc pick) and MAP.52 in one PR;
MAP.60 onward can be a second thread once MAP.64 merges, since they
share files. MAP.55 uses UX.28's icons; if UX.28 isn't done, MAP.55
ships with text labels and gets icons later. Until MAP.86 (phase 1)
lands, the arc pick keeps today's block shading. Arc size: about 40
degrees of bearing by a third of the disk radius (about 27 arcs). At
most 4 build threads run at once (Boss).

## Open questions for Boss

None.
