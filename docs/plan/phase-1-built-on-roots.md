# Phase 1: Built on the groundwork

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

Object references and the picker finished, the database check and repair, resumable runs, the reworked Generate page (spans, radial fills, random neighborhoods, directives, show on the map), galaxy generation changes (phenomena placed galaxy-wide first, fill order, backfill from the run's edge, nebula volume backfill, star-type tiers), the habitability index, tech levels and facility types, spin and orbital-update thresholds with their limits, routing with unknown-space stops and charting along a course, the nearby search, Select mode and view filters on the maps, the fly-through Galaxy Map (MAP.146: scroll-zoom to the cursor, double-click flight, stars that fade in with distance and zoom, a see-through near field; the headline of the 8.1 release), bookmarks, admin control from every screen, in-universe wording, the API log and scopes, and reproducible galaxies up to the golden-seed test.

The headline of this phase, and of the 8.1 release it ends with, is the fly-through Galaxy Map (MAP.146, with MAP.147 to MAP.155): scroll-zoom to the cursor, double-click flight, stars that fade in smoothly with distance and zoom, and a see-through near field (Boss, 2026-10-09 23:07Z). The 8.1 version stamp happens only when the whole phase, these items included, is done.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### References

| ID | Item | Needs | Note |
|---|---|---|---|

### Picker

| ID | Item | Needs | Note |
|---|---|---|---|

### Database consistency check

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.16 | Store the generator epoch and run id on each sector instead of four version text columns | OPS.28 | Research follow-up to DB.7 (built). |
| DB.9 | Repair a damaged galaxy database from a parity file | GEN.57 | Boss 02:13Z: phase 1. Reed-Solomon parity over groups of sector exports; seed regeneration as fallback. |
| DB.17 | Repair by regenerating a damaged sector from its seed when parity cannot rebuild it | DB.9, GEN.57 | Research split of DB.9: the regenerate-from-seed fallback. |

### Generation

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.29 | Record which runs a partly filled sector still needs |  | Builds on the generated flag that GEN.76 settles (sector_stats.bright_level_sol). |
| PERF.30 | Finish an interrupted block or sector run on the next start | PERF.29 |  |

### Generate page

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.181 | Directives: refuse impossible requests up front from the compound-Poisson tables |  | Left over from GEN.96. |
| GEN.180 | Directives: a forced fill (met_forced) after K failed draws |  | Left over from GEN.96. |
| GEN.179 | Store each sector's generation directive and attempt record with the sector |  | Left over from GEN.96. |

### Galaxy gen

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.24 | Generate the galactic core on layer 0 |  | Bulk core fill runs on the parallel path; new mode on generate.html. |
| GEN.155 | A nuclear-cluster object for the Sgr A* sector (optional) |  | Research: optional. |
| GEN.99 | Nebula volume backfill with the star types the nebula needs |  |  |
| GEN.150 | H II region radius and density from the ionizing photon rate, and IMF-based nebula hosts |  | Research: sub-item of GEN.99. |
| GEN.101 | Fill order: nearest sectors first along a pruned Hilbert octree curve |  |  |
| GEN.102 | Investigate filling all near-zero-density void space at once |  |  |
| GEN.103 | Research where each star type and phenomenon belongs in the galaxy's structure |  |  |

### Habitability

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.148 | Habitability index follow-ups from the research (GEN.84 built) |  | Research follow-up to GEN.84 (built). |
| GEN.178 | Magnetic fields: the induced field of an ocean moon |  | Left over from GEN.86 (PR #908); waits on the hydrosphere model. |
| GEN.177 | Planetary magnetic fields: a stagnant-lid factor |  | Left over from GEN.86 (PR #908). |

### Tech levels and facilities

| ID | Item | Needs | Note |
|---|---|---|---|
| POP.8 | A tech-level design from the six domain indices |  |  |
| POP.9 | A tech level generated for every technological species | POP.8 |  |
| POP.7 | Tech levels for technological species | POP.8, POP.9 | Parent. |
| POP.10 | Facility types: programmable, picked from a dropdown, with affiliation and Green/Yellow/Red ratings |  |  |

### Nebula planets

| ID | Item | Needs | Note |
|---|---|---|---|

### Orbital updates

| ID | Item | Needs | Note |
|---|---|---|---|

### Routing

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.47 | Unknown-space jumps stop at scattered stars, black holes, neutron stars and quasars |  |  |

### Nearby search

| ID | Item | Needs | Note |
|---|---|---|---|

### Maps

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.119 | Expected star density editable by admins on the Galaxy Map |  |  |
| MAP.122 | A Select mode on every galaxy view: Galaxy (blocks and sectors) or Star |  | Merges the 2026-10-07 "button to select a star" ask; MAP.101 (PR #431) made stars unpickable by default. |
| MAP.146 | Fly through the galaxy: scroll-zoom, double-click flight, distance-based visibility and a see-through near field | MAP.150, MAP.151, MAP.152 | Boss 2026-10-09 22:24Z and 22:56Z. Umbrella for MAP.148 to MAP.152; replaces the fixed ladder in galaxy/drill.py and galaxyprisms.js. |
| MAP.147 | The Galaxy Map wire format: the investigation (done) and the record of what was built from it |  | Research Lane 3 report in docs/design/galaxy-map-wire-format.md; record until the build items are done. |
| MAP.150 | The free camera: wheel zoom to the cursor, double-click flight, and the observer inside, with the container named from position |  | Fly-through report item 3; needs the near field. |
| MAP.151 | The region data layer: exact-centred frame, aligned cells, per-level aggregates and slot-wrap ranges |  | Fly-through report item 4. Decide cache keys with MAP.147. |
| MAP.152 | Scale hand-offs: galaxy, sector and system cross-fade with hysteresis, and per-tile camera-relative origins | MAP.150 | Fly-through report item 5. |
| TEST.135 | Big-galaxy query-budget test for the nebula reads | TEST.133 | Foundations lane 1. |
| GEN.202 | Scattered nebulae (planetary nebulae and remnants) join the nebula index, deduplicated by object id |  | Foundations lane 1. |
| MAP.174 | Nebula cover in the stage-view cell statistics, for filled and unfilled cells |  | Foundations lane 2. |
| MAP.173 | The Galaxy Map tile carries a nebula layer: regions and singles, nested like the star lists | MAP.172 | Foundations lane 2. |
| MAP.172 | Nebula cell aggregates in the region pyramid | MAP.151 | Foundations lane 1. |
| MAP.171 | The Galaxy Map keys a massive phenomenon's visibility to its map luminosity (replaces "always lit") | MAP.152 | Boss 2026-10-10 21:38Z. Foundations lane 2 after MAP.152. |

### System Map

| ID | Item | Needs | Note |
|---|---|---|---|

### Bookmarks

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.46 | Wiping the galaxy wipes the bookmarks |  |  |
| UX.47 | A bookmark manager: list, rename, sort, group and delete | UX.46 |  |
| UX.49 | Form fields as Shoelace components (sl-input, sl-select, sl-checkbox) across the site |  | Follows UX.40 (done): the fields themselves become Shoelace components. |
| OPS.35 | A vendored-version lock file for the static libraries |  | Research: supply-chain record. |
| OPS.41 | Put the site in an "updating" state while update.sh runs long database migrations |  | Lane 1 report 2026-10-10. |

### Admin control

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.32 | Add a star system to a sector: at the emptiest spot, at given coordinates, or at random outside every Hill sphere |  | Merges two asks (one system added by the computer; placement chosen by admin). |
| MAP.120 | Bright-star backfill from the Galaxy Map's block, slab and wedge menus |  |  |
| ADM.35 | Full control from every screen: edit anything, regenerate with every input, backfill or erase what is in view | MAP.120 | Parent; MAP item MAP.120 is the map part. |

### Pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.23 | A shared unit-ladder module |  | The ladder uses astropy.units. 22 (a feature); UX.36 keeps today's formatter and needs no ladder. Shared unit ladder; UX.22 then UX.30 build on it. |
| UX.81 | Time symbols Gyr, Myr, kyr in place of Gy, My, ky; AU from 1,000,000 km; scientific text below mantissa 1e-3 |  | Research: wording and ladder fixes. |
| UX.78 | Unit preference: Automatic, Metric only or Customary | UX.23 | Research: under UX.23. |
| UX.22 | Meaningful units for every measurement | UX.23 | One quantity family per PR. |
| UX.3 | Warn every visitor while a background job changes the galaxy |  | ETA from the RQ job's published progress. ETA from progress.json, which PERF.23 caps. |
| UX.42 | In-universe wording across the interface |  | After the UX sweep removes controls, so only kept wording is changed. |

### Queue

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.15 | Change the worker count from the Queue page, with a "Ludicrous Speed" mode |  | Worker count now means RQ workers. Changing the worker count live only makes sense once any count works. |

### API

| ID | Item | Needs | Note |
|---|---|---|---|
| API.15 | Log every API call with its user, how it came in, and its HTTP response code |  | Needs no user accounts (Boss 01:31Z). |
| API.4 | API compatibility data in the docs |  | Docs and version number. |
| API.7 | Investigate and plan upload limits |  | Plan only. |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.57 | A sector's contents depend only on the seed, the version and its address |  | Name collisions go away with GEN.67; GEN.63 dropped. Also keys GEN.64's ID collision counter on address (it checks stored rows, PR #406). |
| TEST.77 | A golden-seed regression test | GEN.57, GEN.135, OPS.28 |  |
| GEN.136 | Fingerprint encoding: floats to 9 significant digits, a stored leaf digest per sector and a ring-and-layer tree | GEN.135 | Research follow-up to GEN.58 (built). |
| GEN.135 | Deterministic math helpers for every stored float, with a lint and frozen constants |  | Research: Python 3.9 and 3.10 to 3.13 draw different sectors without it. Blocks TEST.77. |
| OPS.29 | Update reload: cover the gunicorn units and non-Apache hosts, and say what a reload aborts |  | Research follow-up to OPS.8 (built). |
| OPS.28 | Generator epoch and battery digest: say whether two checkouts generate the same galaxy | GEN.135 | Research: OPS.13 is built; this is the epoch it records. OPS.14, OPS.15, OPS.12 and TEST.77 read it. |
| ADM.50 | Remove the galaxy settings files and everything around them; keep only an internal record of the current galaxy |  | Boss 2026-10-10 21:14Z; Foundations lane 1, before OPS.28. |
| OPS.37 | A Generator version number: one sequential integer, shown in admin and the API | OPS.28 | Boss 2026-10-09 20:59Z: done by the end of phase 1. Same number as OPS.28's epoch (decided). |

### Bugs from the GitHub issues

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.111 | test_ensure_sector_generated_creates_then_reuses_the_same_sector fails in a busy parallel run (bug) |  | Test flake reported 2026-10-09; Bugfixes lane 1. |
| TEST.132 | test_a_fill_after_a_layer_failed_mid_scatter_builds_no_leftover_star failed once in a full parallel run (bug) |  | Lane 2 report 2026-10-10 21:06Z; Bugfixes lane 2 after its current items. |
| TEST.137 | test_a_neighborhood_from_the_generate_page[core] fails once under load (bug) |  |  |
| TEST.138 | test_open_map_menus_hold_no_overlap[web.galaxy-1280] fails once in a full parallel run (bug) |  |  |
| TEST.147 | Two nebula-field table tests fail under load and pass alone (bug) |  |  |
| TEST.149 | Two physics tests fail under load and pass alone (bug) |  |  |
| TEST.148 | Investigate resuming a test run where it left off when the build is the same |  |  |
| DOC.44 | Add Redis to the install and deploy docs: it is required, and no guide lists it |  |  |
| DOC.43 | Fix the Python 3.12 and numba leftovers in the install and CI docs and the 3.10 check in the installer |  |  |
| DOC.42 | Describe what update.sh and install.sh really do today in INSTALL.md and the Apache and macOS guides |  |  |
| DOC.41 | Fix the activity-log file name in the fail2ban guide and the shipped jail |  |  |
| DOC.40 | Bring docs/testing.md and docs/ci-runners.md in line with ci.yml and scripts/test.py |  |  |
| DOC.39 | Fix the README facts, the documentation table and the headline rule (Boss decision) |  |  |
| DOC.38 | Refresh docs/config.md prose after PERF.64 and PERF.24: page_cache, stored_counts, PLANETGEN_WORKERS |  |  |
| DOC.37 | Rewrite the install docs for the venv install, in or right after OPS.45 |  |  |
| DOC.36 | Sync docs/cli.md with the 11 subcommands, the options and the real bright-star defaults |  |  |
| DOC.35 | Document the helper command-line tools: orbits, warm_map, worker, query, migrate, reset, lockouts, dedupe, config |  |  |
| DOC.34 | Bring docs/api.md up to date: endpoints, object references, errors and the Not done yet list |  |  |
| DOC.33 | Rewrite the database-schema.md Versioning, Persistence and Planned changes sections for Alembic and v83 |  |  |
| DOC.32 | Update database-schema.md: control schema, missing tables and columns, old module names |  |  |
| DOC.31 | Fix docs/html-interface.md: routes, header sections, search parameters, error pages and rate limits |  |  |
| DOC.30 | Fix the system-file-format.md class table, the moon blacklist and the example-systems.md module names |  |  |
| DOC.29 | A docs link check in the test suite, and a note that the web interface has no help or about page yet |  |  |
| DOC.28 | Rewrite the code map in docs/design/architecture.md for the planetgen packages |  |  |
| DOC.27 | Update the architecture.md flows and add a design-notes index at docs/design/README.md |  |  |
| DOC.26 | Give every design note a status line and refresh the stale ones |  |  |
| DOC.25 | Refresh docs/design/design-decisions.md: mark the reversed and built decisions |  |  |
| DOC.24 | Bring Python 3.12, numba and astropy into the research notes and three source comments |  |  |
| DOC.23 | Update the nebula, object ID and map-data-layer notes to what merged |  |  |
| DOC.22 | Update the orbital, slow-reads, login, ops, test-suite and db-check notes to their built state |  |  |
| DOC.21 | Rewrite the package docstrings that are wrong or missing, and the queryDb mentions |  |  |
| DOC.20 | After PERF.86: rewrite the queue and scatter notes and the queue docstrings | PERF.86 |  |
| DOC.19 | After UX.102: document the log output standard for the CLI, the Generate page and admin logs | UX.102 |  |
| DOC.18 | Last pass: re-run the documentation sweep when the other DOC items are done and Phase 1 is complete |  |  |
| TEST.144 | test_every_plain_button_does_something[web.sector] failed once in a full parallel run (bug) |  |  |
| PERF.91 | Bulk random draws for the scatter layers (4 to 12 times on batches above 1,000 draws) |  |  |
| PERF.88 | Sector-path integrator on numba: saving the sector paths takes 13 s instead of 34 s on the benchmark fill |  |  |
| TEST.141 | test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector[15-True] fails once in a full parallel run (bug) |  |  |
| TEST.140 | test_a_body_is_reached_by_its_id_and_a_row_number_names_nothing fails once in a full parallel run (bug) |  |  |
| MAP.189 | The Galaxy Map opens flown in to the created sectors | MAP.150, MAP.152 |  |
| TEST.139 | Rework, optimise and rebuild the whole test suite so it runs faster and in modules |  |  |
| TEST.130 | test_the_scene_positions_at_the_epoch_match_the_stored_ones fails when a random system has comets (bug) |  | Lane 1 report 2026-10-10 20:55Z; Bugfixes lane 1 after its current items. |
| TEST.134 | test_many_processes_starting_at_once_start_one_job fails under a parallel run (bug) |  | Flaky under load. Bugfixes lane 1. |
| TEST.129 | Two orbit-update tests fail on the MySQL 8.4 leg: the stored last_updated_at rounds up (bug) |  | Lane 1 report 2026-10-10 20:55Z; Bugfixes lane 1 after its current items. |
| TEST.128 | test_spatial_position_db and test_web_db_fields fail under a parallel run on one MySQL (bug) |  | Lane 1 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. |
| TEST.127 | test_nebula_shape_endpoint_serves_a_mesh failed once in a full parallel run (bug) |  | Lane 2 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. |
| TEST.126 | The heavy WebGL browser-map tests time out when four workers run them together (bug) |  | Lane 2 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. |
| TEST.116 | test_bughunt_end_to_end stores M2V or M6V where it expects K2V on the MySQL 8.4 and MariaDB legs once (bug) |  | Bugfixes lane 1 report 03:57Z. |

### Foundations for the issue features

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.47 | The PERF.31 benchmark records the buffer pool, table sizes and worker start-up cost |  | Generation performance study. |
| PERF.41 | Stamp the tile and page caches with a TILE_FORMAT constant instead of the version (optional) |  | Research follow-up to PERF.25 (built). |
| PERF.40 | Two shared queues, a reserved interactive worker and a real "cancel now" |  | Research follow-up to PERF.24 (built). |
| PERF.100 | Work queue step 6: delete the old machinery and the work_* tables (Alembic migration) | PERF.99 | PERF.86 step 6 |
| PERF.99 | Work queue step 5: the Generate page, Queue page, banner and job log read Redis; cancel now through the stop command | PERF.98 | PERF.86 step 5 |
| PERF.98 | Work queue step 4: the remaining stages and the sector and backfill units on the new queue | PERF.97 | PERF.86 step 4 |
| PERF.97 | Work queue step 3: the mass scatter end to end on the new queue | PERF.96 | PERF.86 step 3 |
| PERF.96 | Work queue step 2: the run record, stage counters, event stream and the ETA function, as a library | PERF.95 | PERF.86 step 2 |
| PERF.95 | Work queue step 1: the queue spine, two queues, on-demand workers with an idle timeout and the wanted-worker count |  | PERF.86 step 1 |
| ADM.54 | A Redis health check on the admin page: answering, append-only file on, memory policy, memory used, last save, worker count (optional) |  | Redis review 5.11 item 5, optional |
| OPS.48 | Redis server settings: a documented example config and an update.sh warning (append-only file, noeviction, memory limit, loopback) |  | Redis review 5.11 item 2 |
| OPS.49 | update.sh stops at the database migrate step with 'NoneType' object has no attribute 'is_auth_switch_request' (bug) |  |  |
| TEST.133 | A reusable big-galaxy query budget test: EXPLAIN every page and list query on 2,000,000 systems |  | Foundations lane 1. |
| PERF.78 | A reserved warm worker for long admin operations, and the poll pattern for them |  | Foundations lane 1, last of the PERF.71 builds. |
| PERF.75 | Keyset paging for the other data tables and the jump-by-value box (Systems name sort done) |  | Decided by default. Foundations lane 1. |
| PERF.33 | Progress bars and ETAs from measured performance |  | GitHub issue [#661](https://github.com/dwhagar/planetGen/issues/661) (the use half). Bugfixes lane 1 (Boss, 23:38Z). |

## Open questions for Boss

- UX.22: Meaningful units for every measurement, see its entry in TODO.md.
- UX.3: Warn every visitor while a background job changes the galaxy, see its entry in TODO.md.
- API.4: API compatibility data in the docs, see its entry in TODO.md.
- OPS.13: Every update records the version key, keeping the last 10, see its entry in TODO.md.
