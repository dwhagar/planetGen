# Version 8.1: Infinite zoom and the new generate scatters

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.1

## Goal

Two headlines. Infinite zoom: one camera flies from the whole galaxy to a star system (MAP.146, done: every sub-item merged and Core 2 checked each clause of its Done text against main, 3:20 AM EDT 2026-10-23 as the lane reported). The new generate scatters: accurate galactic generation and scattering out to the 262 kpc reach, including the halo, the satellites, the nebula stars and a realistic core (GEN.301). The version also holds the Phase 1 tail: the queue steps, reproducible galaxies, the description layer, tech levels, directives and every open bug and flaky test. A third strand (Boss, 2026-10-16, 11:57 PM EDT): the building blocks of the 9.0 composition headline, which go down to 8.1 because they change no stored value: the research lanes and designs (phase 1, GEN.431 to GEN.437 with GEN.409, GEN.425 and GEN.426), the data tables (phase 2, GEN.438, GEN.410, GEN.427, GEN.439), the star composition function (phase 3, GEN.440) and its sweep (phase 4, TEST.211); the master note is docs/design/composition/composition-model.md. The 8.1 stamp waits until all of it is done (REVISION_HOLD in scripts/bump_version.py).

## Where it stands

As of 10:00 PM EDT 2026-10-23: 774 of 1182 items finished (65.5 percent), 408 open, 23 dropped (PERF.162 and PERF.166 moved to the later version, Boss 12:27 PM EDT). Open by phase: phase 1 114, phase 2 41, phase 3 25, phase 4 14, phase 5 83, phase 6 35, phase 7 22, phase 8 6, phase 9 17, phase 10 4, phase 11 39, phase 12 8 (phases 6 to 11 are the composition and sky colour suite pulled in from 8.4 and 9.0 by Boss at 5:48 PM EDT; phase 12 is the headline). Only this version is being built (Boss, 7:44 PM EDT 2026-10-22); the database chain finished with DB.85 (#2582, revision 0130, no generator epoch, so epoch 51 was still free then; since, 51 to 53 are taken: GEN.307, GEN.315 slice 1 and GEN.628, and 54 is next). Closing order: the audits TEST.185, TEST.184, TEST.159 and TEST.160 (merged, re-run as TEST.283), PERF.146 (merged as #2559) with PERF.164 and PERF.165, then SEC.33 and TEST.261 in their own threads, then GEN.301. Lane queues: docs/plan/lanes.md and reports/8-1-core-queues-2026-10-22.md.

## Phase 1: Foundations and open trees

Everything with no blocker outside its own thread: the queue steps, reproducible galaxies, the description layer, tech levels, the settings page (moved up from 8.2, Boss 1:43 PM EDT), directives, housekeeping, the bugs and the backfill weeding.

### Queue, caches and performance

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.41 | Stamp the tile and page caches with a TILE_FORMAT constant instead of the version (optional) |  | Research follow-up to PERF.25 (built). |
| PERF.86 | Redo the work queue: drop the task and subtask logic and use Redis simply |  |  |
| PERF.113 | The neighbour link of a galaxy run: measure it, and split it into chunks of sectors if the time warrants |  | Audit of 2026-10-17, 4:06 PM EDT (Boss). |
| PERF.114 | The population pass of a galaxy run: measure it, and split it into chunks of sectors if the time warrants |  | Audit of 2026-10-17, 4:06 PM EDT (Boss). |
| PERF.115 | The settle stage of a galaxy run: measure it, and split it into chunks of sectors if the time warrants |  | Audit of 2026-10-17, 4:06 PM EDT (Boss). |
| PERF.116 | The habitability re-score (`plan --habitability`): measure it, and split it into chunks of systems if the time warrants |  | Audit of 2026-10-17, 4:06 PM EDT (Boss). |
| PERF.117 | The map-luminosity fill (`plan --map-luminosity`): measure it, and split it into id-range chunks if the time warrants |  | Audit of 2026-10-17, 4:06 PM EDT (Boss). |
| PERF.98 | Work queue step 4: the remaining stages and the sector and backfill units on the new queue |  | PERF.86 step 4 |
| PERF.119 | The Generate page, the API and the one-off sector and system commands start their runs through the wrapper | PERF.98 | Boss, 2026-10-17, 4:12 PM EDT: every start path uses the wrapper. |
| PERF.149 | A map visit queues the missing sector as a run instead of generating it in the web request, and the Galaxy View is closed while a galaxy generation runs | PERF.119 |  |
| PERF.120 | Every command does one thing: the umbrella for reshaping every command and entry point |  | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| PERF.123 | `galaxy` stops being a bundle: sector generation, the neighbour link, the backfill, the population pass and the settle are each their own command, and `galaxy` is a step list the wrapper builds |  | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| PERF.172 | The galaxy as one wrapper run: the population pass and the phenomenon work are units of the step list (the D2 slice of PERF.123) | PERF.124, PERF.125 | Split from PERF.123. |
| PERF.124 | `sector` and `phenomenon` only generate: the population pass leaves `sector` and the settle leaves `phenomenon` | PERF.123 | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| PERF.125 | The population pass and the territories pass are separate commands | PERF.120 | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| PERF.128 | `benchmark` is a step list handed to the wrapper, not a sequence of its own | PERF.123 | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| GEN.309 | The backfill picks the draw path per block and gives the same stars wherever it lands |  | Split from GEN.42. |
| GEN.40 | Weed out sectors by star density before the bright-star backfill | GEN.309 | Parent; closes with its subitems. |
| PERF.35 | An interval or chunk ledger for untouched sectors once block-first backfill lands | GEN.309 | Research: ledger size at full galaxy. |
| PERF.99 | Work queue step 5: the Generate page, Queue page and job log read Redis; cancel now through the stop command | PERF.98, PERF.119 | PERF.86 step 5 |
| PERF.154 | One waiting list for runs in Redis (up to 8 runs, a setting), replacing next.json, the file lock and Busy (F5) | PERF.119, PERF.123, PERF.99 |  |
| PERF.147 | A per-task progress bar with an ETA for any task expected to take, or taking, longer than 15 seconds |  | Boss, 8:46 PM EDT 2026-10-21 (relayed). |
| GEN.310 | Fill a whole block of the shell ladder: planetgen generate --block and the queue job | PERF.99 | Split from GEN.231. |
| UX.127 | The Galaxy Map block menu offers Generate again | GEN.310 | Split from GEN.231. |
| GEN.231 | Generate a whole block on the shell ladder: the generate --block mode and the block buttons come back | GEN.310, UX.127 | Maps hand-off 2026-10-15. Core, then Interface for the buttons. |
| DB.29 | update.sh and every database migration refuse to start while jobs are active, unless --yes is passed or the user confirms cancelling them | PERF.99 |  |
| PERF.100 | Work queue step 6: delete the old machinery and the work_* tables (Alembic migration) | DB.29 | PERF.86 step 6 |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.55 | Same seed, same data: a sector's contents depend only on the seed, the version and its address (internal) |  | Parent of the chain. |
| GEN.57 | A sector's contents depend only on the seed, the version and its address | GEN.633 | Name collisions go away with GEN.67; GEN.63 dropped. Also keys GEN.64's ID collision counter on address (it checks stored rows, PR #406). |
| GEN.633 | The backfill after a run puts its dimmer bright stars only in sectors not yet generated, so a sector holds different systems depending on the order the sectors were generated in (bug) |  |  |
| GEN.179 | Store each sector's generation directive and attempt record with the sector |  | Left over from GEN.96. |
| DB.45 | Storage checksum and the dirty flag: the export, the checksum and every writer |  | Split from DB.9. Core 1 owns migrations. |
| DB.84 | Every schema change on `sectors` rebuilds the whole table because of its FULLTEXT index: decide how the later DB items pay for it |  |  |
| DB.83 | `galaxy_block_stats` has no maintainer: nothing keeps the block totals current after a scatter, a fill or an edit (bug) |  |  |
| DB.46 | The parity file: groups, records, a delta update on every save and the stale-record rule | DB.45 | Split from DB.9. Core 1 owns migrations. |
| DB.47 | planetgen repair-db: rebuild damaged sectors from the parity file | DB.46 | Split from DB.9. Core 1 owns migrations. |
| DB.48 | repair-db: rebuild the table when a clustered-index page is unreadable | DB.47 | Split from DB.9. Core 1 owns migrations. |
| DB.9 | Repair a damaged galaxy database from a parity file | DB.45, DB.46, DB.47, DB.48 | Boss 02:13Z: phase 1. Reed-Solomon parity over groups of sector exports; seed regeneration as fallback. |
| DB.17 | Repair by regenerating a damaged sector from its seed when parity cannot rebuild it | GEN.57, DB.47, GEN.179 | Research split of DB.9: the regenerate-from-seed fallback. |
| OPS.15 | Each update says whether it changes generated output |  | Needs the fingerprint, so phase 2. |
| OPS.29 | Update reload: cover the gunicorn units and non-Apache hosts, and say what a reload aborts |  | Research follow-up to OPS.8 (built). |
| OPS.37 | A Generator version number: one sequential integer, shown in admin and the API |  | Boss 2026-10-09 20:59Z: done by the end of phase 1. Same number as OPS.28's epoch (decided). |
| TEST.77 | A golden-seed regression test | GEN.57 |  |

### Tech levels and species

| ID | Item | Needs | Note |
|---|---|---|---|

### Description layer

| ID | Item | Needs | Note |
|---|---|---|---|

### Admin and API housekeeping

| ID | Item | Needs | Note |
|---|---|---|---|
| API.26 | Bump API_VERSION for the shell, band and wedge sector address |  |  |
| API.4 | API compatibility data in the docs | API.26 | Docs and version number. |
| API.15 | Log every API call with its user, how it came in, and its HTTP response code |  | Needs no user accounts (Boss 01:31Z). |
| API.31 | Every outgoing web request sends a standard User-Agent with the software, its version and the API version |  | Needs nothing. Boss, 2026-10-17, 11:33 AM EDT. |
| USR.31 | The Owner account stub: one Owner row and a role check, enough for the Owner-tier settings of ADM.86 in 8.1 |  | Boss, 2026-10-18 (fork 9, stay in 8.1). |
| ADM.86 | Owner-tier settings: Owner-only options behind a recent-password prompt | USR.31 | Boss 1:43 PM EDT 2026-10-16: the settings page is 8.1; stays in 8.1 with the Owner stub (Boss, fork 9, 2026-10-18). |
| OPS.35 | A vendored-version lock file for the static libraries |  | Research: supply-chain record. |
| UX.136 | The data tables' own row fetches (scroll and sort) wait on a page job when the answer is slow, like the pages do |  |  |
| UX.121 | The logging layer and the written standard |  | Split from UX.102. |
| UX.122 | The command line, the job log files and the downloadable log use the layer | UX.121 | Split from UX.102. |
| UX.123 | The Generate page log, the admin log pages and update.sh output use the layer | UX.121 | Split from UX.102. |
| TEST.183 | A test that fails when a new output path formats its own lines | UX.122, UX.123 | Split from UX.102. |
| UX.102 | One logging standard for every output the user or an admin reads: the same timestamps, format, dividers and step markers everywhere | UX.121, UX.122, UX.123, TEST.183 |  |


### Generate page and directives

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.180 | Directives: a forced fill (met_forced) after K failed draws, and refusal of impossible requests up front |  | Left over from GEN.96. |
| GEN.229 | Repair or re-roll the wide binaries that already break the stability rule |  |  |
| UX.148 | Log output from Generate is no longer in a fixed-width font (bug) |  |  |

### Bugs and flaky tests

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.111 | test_ensure_sector_generated_creates_then_reuses_the_same_sector fails in a busy parallel run (bug) |  | Test flake reported 2026-10-09; Bugfixes lane 1. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.116 | test_bughunt_end_to_end stores M2V or M6V where it expects K2V on the MySQL 8.4 and MariaDB legs once (bug) |  | Bugfixes lane 1 report 03:57Z. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.126 | The heavy WebGL browser-map tests time out when four workers run them together (bug) |  | Lane 2 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.307 | test_systems.py::test_asteroid_belt_and_habitable_world_forced_together_both_succeed[A0V] fails about once in 5 to 8 runs (bug) |  | Flake, fix the cause |
| TEST.306 | Five physics modules are imported only by tests and no item names them |  | Owner or delete |
| TEST.305 | The sector-map browser tests wait 120 s for the star list and fail on a fresh four-core box (bug) |  | Run on a normal lane first |
| TEST.304 | Nine Playwright test files are not marked browser, so they run in the all-core phase |  | Marker fix |
| TEST.127 | test_nebula_shape_endpoint_serves_a_mesh failed once in a full parallel run (bug) |  | Lane 2 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.128 | test_spatial_position_db and test_web_db_fields fail under a parallel run on one MySQL (bug) |  | Lane 1 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.130 | test_the_scene_positions_at_the_epoch_match_the_stored_ones fails when a random system has comets (bug) |  | Lane 1 report 2026-10-10 20:55Z; Bugfixes lane 1 after its current items. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.137 | test_a_neighborhood_from_the_generate_page[core] fails once under load (bug) |  | Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.140 | test_a_body_is_reached_by_its_id_and_a_row_number_names_nothing fails once in a full parallel run (bug) |  | Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.141 | test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector[15-True] fails once in a full parallel run (bug) |  | Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.190 | A failing session fixture makes every dependent test wait its full 900 s timeout: fail the dependents at once (bug) |  | Test-infrastructure finding. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.194 | `test_system_map_selection_drill_and_measure` failed once in the browser phase of a full run under load and passed alone: find the cause and make it robust (bug) |  | Flaky test finding from the Settings lane, 5:47 PM EDT. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.237 | The wheel-anchor browser test failed alone at box load 9.7 (call 5.47 s) after the TEST.231 fix merged (bug) |  |  |
| TEST.254 | test_core_ball_batch_fill_makes_planets builds 4 systems instead of 3 for some random galaxy seeds (bug) |  |  |
| OPS.62 | The cloud environment setup script fails at its verify step because it imports `markdown`, which nothing uses (bug) |  |  |
| OPS.65 | Debug mode wrote almost nothing to stdout: planetgen system --debug printed 3 lines because the generation runs in queue workers whose stdout is /dev/null (bug) |  |  |
| OPS.66 | Remove the unused older executor in src/planetgen/queue/work.py: RQExecutor, _rq_task and the --debug FILE worker plumbing | OPS.65 |  |
| OPS.67 | Detached helpers (the background planetgen warm a run starts, an auto-started worker's own lifecycle) write the debug log file only: decide whether they relay to stdout | OPS.65 |  |
| TEST.302 | benchmark_suite.stop_workers() kills other tests' queue workers: test_benchmark_suite fails or breaks test_queue_pause, test_reset_jobs, test_workers and test_api_jobs with exit -9 when the data sets are fetched (bug) |  |  |
| TEST.303 | test_update_database_step::test_a_failed_migration_stops_the_step_and_keeps_the_version fails on untouched main on the lane machines because the test database user lacks GRANT (bug) |  |  |
| TEST.262 | test_a_double_click_flies_even_when_the_page_is_slow_to_answer_the_first_click fails once under load: the camera target stays at (3.56, -1.36, -0.29) instead of (0, 0, 0) (bug) |  |  |
| TEST.273 | test_the_index_takes_its_size_off_the_tile_cache_budget fails once in a full run: max_cache_bytes 1073577984 against 1073741824 (flake) |  |  |
| TEST.284 | test_an_uncharted_black_hole_page_names_its_class fails under parallel load in the changed gate: the page renders without the class name (bug) |  |  |
| TEST.287 | test_web_browser_fixture_maps.py::test_galaxy_breadcrumb_always_matches_the_view_its_url_names times out in the browser phase under load (bug) |  | Bug Hunt; Core 2 finding |
| TEST.289 | Two API tests fail once under the gate's four workers: test_a_neighborhood_job_runs_its_sectors_through_the_wrapper and test_a_quick_edit_answers_in_the_same_response (bug) |  |  |
| TEST.290 | test_work_units and test_galaxy_seed each fail once under the gate load: a galaxy run stores each stage once, and a new seed replaces the old only before any sector (bug) |  |  |
| TEST.292 | test_web_browser_sector_reach.py::test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star fails at 130 s under load 9.7 and takes 97 s alone (bug) |  |  |
| TEST.298 | test_admin_script_cli.py: four fixture-finalizer errors under load in a gate (bug) |  | Core 2 found it. |
| TEST.299 | test_api_species_tech_filters fails once in a full run: the first species of the tech-descending list has tech 6.55, not the 6.0 the test set (bug) |  |  |
| TEST.297 | The 13 tests of test_work_queue_failures.py and test_web_browser_sector_open.py::test_opening_a_sector_is_one_continuous_view_with_a_filled_panel fail in a gate at load 6 to 8 and pass alone (bug) |  | Core 2 found it. |
| TEST.295 | test_parallel_generation.py::test_one_worker_draws_what_two_and_four_workers_draw hit the 900 s timeout waiting for its queue workers in a saturated gate (bug) |  |  |
| TEST.276 | test_clearing_a_scatter_with_nothing_stored_says_so_and_changes_nothing fails once under parallel load with empty captured output (flake) |  |  |
| TEST.234 | test_the_review_and_error_pages_are_accessible[admin] found 606 px of width on a 390 px page once, at the site gear dropdown (bug) |  |  |
| TEST.233 | test_a_measured_feature_moves_by_the_percentage[asteroid_belt-30] failed once in a loaded gate and passed 13 of 13 alone: its seed changes every run (bug) |  |  |
| TEST.229 | The browser phase of a gate stalls for minutes at one test with a headless_shell gpu process near 200% CPU (bug) |  |  |
| GEN.542 | scripts/star_type_by_radius.py fails with a TypeError because the plan arguments now default the disk and arm options to None (bug) |  |  |
| PERF.111 | `ResponseCache.put` without a database never trims the page cache: the cache grows past `max_entries`, and the thread test fails when its last put is one of those (bug) |  | Core 2 finding in the MAP.217 gate: the db=None branch of ResponseCache.put has no eviction. |
| TEST.219 | Browser pages still have 1 to 7 requests running 8 s after the server stops: find what keeps them busy and make it faster |  | Core 1 finding at #1997: browser pages still have requests running 8 s after the server stops. |
| TEST.218 | `test_first_control_connections_at_once` leaves 3 idle connections open on `planetgen_test_*` when its fixture ends: the test should close its own connections (bug) |  | Core 1 finding at #1997: the test leaves 3 idle connections for the drop warning to end. |
| TEST.198 | `test_galaxy_map_does_not_prefetch_on_a_data_saving_connection` failed with a `/galaxy/tiles` net::ERR under load and passed alone: find the cause and make it robust (bug) |  | Flaky test finding from Core 3, 6:22 PM EDT. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| GEN.398 | The massive-star and bright-star passes place more stars than their log lines expect: 70 and 12 percent over (bug) |  | Rescatter thread finding; Core 3 proposed (scatter spine). |
| GEN.396 | Revisit the Yerkes class III luminosity range: 50 to 1000 solar luminosities sits above many dim red giants |  | Core 3 finding, 2026-10-16 10:00 PM EDT; any Core lane. |

### Nebula data foundations

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.202 | Scattered nebulae (planetary nebulae and remnants) join the nebula index, deduplicated by object id |  | Foundations lane 1. |

### Composition research and design

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.560 | Research list for the chain: ejection survival of moons, circumplanetary disks of isolated planetary-mass objects, failed-supernova and direct-collapse shares, the size law of interstellar objects |  |  |

### Admin control suite

Boss (6:01 PM EDT, 2026-10-21): the admin pages, the Generate suite and the settings are one control suite, built in this order: Settings first, then the shell, job strip and palette, then Generate. Boss answered all 13 decision cards between 6:05 and 6:09 PM EDT (design note docs/design/interface/admin-control-suite.md, section 7); his answer 1 puts all of it in 8.1 phase 1. Owner of every item: the admin redesign thread.

Shape dials (Boss, 6:41 to 6:46 PM EDT, 2026-10-21; goal issue #2355): the generator's shape inputs become five dials built from one lever table (GEN.601), so the Generate suite is not built on the retiring absolute inputs. The table, flags (GEN.602), Scale (GEN.603), Structure (GEN.604), Edge (GEN.605, both forms with a unit selector: L, pc, kpc, Mpc), Halo (GEN.606), the unfilled Shape selector (GEN.608), range sources (GEN.609), the nuclear-amplitude fix (GEN.610), the pattern dial (GEN.611), the dwarf-mass law extension (GEN.612, needs Boss's sources) and the stored pins (DB.81) sit here ahead of the Generate rows; the dials panel ADM.124 follows ADM.110; Census (GEN.607) waits on GEN.545 in the galaxy laws section. The default galaxy is also smaller (Boss, 6:57 PM EDT: still a spiral, no more than half the Milky Way's size; GEN.603 states the number): a new generator epoch, with lower system counts and cost. The Shape dropdown is filled in 9.4 by UX.134 and GEN.389. The lever-table items belong to the Galaxy dials lane (the coordinator, 6:55 PM EDT); Census has no owner until GEN.545 merges (the first release has four dials; Boss, 7:01 PM EDT).

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.125 | The Generate page and the one-off system page on the field macro |  | The Generate page conversion rides ADM.99. Admin control suite, Boss 6:01 PM EDT 2026-10-21: the admin redesign thread; research and cards first. |
| GEN.612 | Extend the halo, stellar-to-halo and size laws down to dwarf and ultra-faint masses (below 1e8 solar masses), so the shape dials' minimum mass is reachable |  | Shape dials (Boss 7:13 PM EDT 2026-10-21); needs Boss sources. |
| GEN.613 | Investigate the limits of the density algorithm: for which mass, size, flatness and bulge share the galaxy density model stays valid, and what it does outside them |  | Boss 7:33 PM EDT 2026-10-21; read-only report. |
| GEN.614 | The nuclear disc is not softened to the grid: switch it off when its length (0.057 times the bulge radius) is under 2 pc, or soften it like the nuclear cluster |  | Density limits report (GEN.613). |
| GEN.615 | Galaxies below the low-mass cut are one disc plus a bulge: no thick disc, bar, stellar halo or nuclear structures |  | Density limits report (GEN.613); cut decided 3e8 (Boss, 7:58 PM EDT). |
| GEN.618 | The stellar-to-halo inversion gives a 6e13 solar-mass halo for 2e11 solar masses of stars, the shape dials' mass ceiling: check the Moster inversion at the top end |  | Density limits report (GEN.613). |
| UX.126 | The admin pages and the edit dialogs on the field macro |  | The admin page conversion rides ADM.56 and ADM.93. Admin control suite, Boss 6:01 PM EDT 2026-10-21: the admin redesign thread; research and cards first. |
| UX.49 | Form fields as Shoelace components (sl-input, sl-select, sl-checkbox) across the site |  | Follows UX.40 (done): the fields themselves become Shoelace components. |

## Phase 2: Infinite zoom

The fly-through Galaxy Map (MAP.146) and the nebula data its tiles carry.

### Fly-through map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.75 | The mini map as a second engine view |  | Locked second camera on MAP.64. |
| MAP.140 | Double-click on a selected object goes there and opens its information |  | GitHub issue [#714](https://github.com/dwhagar/planetGen/issues/714). |
| MAP.141 | Context around the selection: faint neighbours, and the sectors above and below |  | GitHub issues [#718](https://github.com/dwhagar/planetGen/issues/718) and [#716](https://github.com/dwhagar/planetGen/issues/716), one context view. The #716 bug label was overruled. |
| MAP.249 | The Galaxy Map is laggy and slow, and changing zoom level slogs the whole page (bug) |  | Boss 6:45 PM EDT; Core 2 profiles first. |
| MAP.258 | N1: the sector-block ladder above 243 sectors and drillChildren in JavaScript |  |  |
| MAP.259 | N2: node aggregates for the pre-placed bright stars and scattered phenomena | MAP.258 |  |
| MAP.260 | N3: node read and cache: GET /galaxy/nodes, the node cache, the change feed and the warm-up | MAP.258, MAP.259 |  |
| MAP.261 | N4: the node cut on the page: galaxynodecut.js replaces the tile cut | MAP.260 |  |
| MAP.262 | N5: retire the cube tiles | MAP.261 |  |
| MAP.263 | S1: the camera is the only state: stages follow the camera, with no gesture gate, preset flights, per-stage zoom range or pan clamp |  |  |
| MAP.264 | S2: nothing is clipped when you zoom: drop setWedgeClip and markClippedStars | MAP.263 |  |
| MAP.265 | S3: cells are a tool, not a display stage: the hover and selected sector outline, an optional grid tool and a closed-form pick | MAP.263 |  |
| MAP.266 | S4: the unresolved glow: a soft light drawn from the node aggregates, switchable in the layers menu | MAP.259, MAP.261, MAP.265 |  |
| MAP.269 | Stars on the maps do not glow: a star further than about 3 sectors is a 1.5 px point with no halo and a mostly white core (bug) |  | Boss bug 8:51 PM EDT |
| MAP.267 | S5: one star buffer in relative-to-eye pairs, filled node by node (MAP.240 folded in) | MAP.261, MAP.263 |  |
| MAP.268 | S6: labels and selection on the cut: MAP.229 to MAP.233 on the projected-size query, selection by object ID across scales | MAP.263, MAP.264, MAP.265 |  |
| MAP.257 | Infinite zoom is still clunky: the zoom zones and segments other than sectors get in the way (bug) | MAP.258, MAP.259, MAP.260, MAP.261, MAP.262, MAP.263, MAP.264, MAP.265, MAP.266, MAP.267, MAP.268 |  |
| MAP.256 | The map shows the information box twice when an object is selected: once in the window and once to its right (bug) |  |  |
| MAP.121 | Every map shows and steps to its neighbouring regions, on one map engine |  | Folds MAP.127, the Galaxy Map's neighbouring blocks and slabs (Boss 2026-10-07 16:27Z: one map engine for all maps). |
| MAP.122 | A Select mode on every galaxy view: Galaxy (blocks and sectors) or Star |  | Merges the 2026-10-07 "button to select a star" ask; MAP.101 (PR #431) made stars unpickable by default. |
| MAP.171 | The Galaxy Map keys a massive phenomenon's visibility to its map luminosity (replaces "always lit") |  | Boss 2026-10-10 21:38Z. Foundations lane 2 after MAP.152. |
| MAP.195 | Galaxy Map finer tiles leave out the stars a coarser tile already sent |  |  |
| MAP.173 | The Galaxy Map tile carries a nebula layer: regions and singles, nested like the star lists |  | Foundations lane 2. |
| TEST.135 | Big-galaxy query-budget test for the nebula reads | MAP.173, GEN.202 | Foundations lane 1. |
| MAP.156 | View one layer or a range of layers top-down from the galaxy view, as a secondary option |  | Boss 2026-10-09 23:17Z; secondary option, not a main control. |
| MAP.160 | Quantise the Galaxy Map GPU buffers (deferred) |  | Deferred; MAP.147 recommendation step 3. |
| MAP.190 | Several separate regions of created sectors: which one the Galaxy Map opens on |  |  |
| MAP.189 | The Galaxy Map opens flown in to the created sectors | MAP.190 |  |
| MAP.214 | Galaxy Map default view: no density, age or luminosity shading, only stars and markers for known space |  | Split from MAP.191. |
| MAP.209 | Opening a sector from the Galaxy Map flashes one view, reloads into another, and leaves the details panel empty (bug) |  |  |
| MAP.210 | Nothing on the map can be selected, even in sector view (bug) |  |  |
| MAP.237 | The Galaxy Map's "Charted only" view is invisible again, and charted regions are shown by colour and opacity, not outlined blocks (bug) |  | Boss 7:24 PM EDT. |
| MAP.229 | Map label engine: names placed by priority, with no overlaps, shown and hidden by zoom |  | Map labels, Boss 2:30 PM EDT. |
| MAP.234 | Galaxy Map key: a Key section in the help dialog and an on-map key that names everything the map draws |  | Boss 7:15 PM EDT; the superbubble layer is its sibling item. |
| MAP.240 | Browser buffers hold camera-relative coordinates: territory dots, superbubbles, made rings, block prisms, stage views, system view dots and belts |  | GEN.372 audit; Boss card. |
| MAP.236 | Superbubbles need clearer bounds and slightly more transparency on the Galaxy Map | MAP.240 | Boss 7:24 PM EDT; look only. |
| MAP.241 | Stars on the star maps should glow brighter, not grow bigger: size stays small, glow and colour carry luminosity and temperature (bug) |  | Boss 7:43 PM EDT. |
| MAP.230 | Names of stars, systems and sectors on the Galaxy Map | MAP.229 | Map labels, Boss 2:30 PM EDT. |
| MAP.231 | The Sector Map and the System Map use the label engine | MAP.229 | Map labels, Boss 2:30 PM EDT. |


### Golden galaxy

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.182 | Put one block-first backfill into the golden galaxy | GEN.309, TEST.77 | Split from GEN.42. |
| GEN.42 | A block-first draw for the backfill (umbrella) | GEN.309, TEST.182 | Same function as PERF.18 (backfill_bright_stars_around); one thread. |

### Star-system storage

| ID | Item | Needs | Note |
|---|---|---|---|

### Composition tables

| ID | Item | Needs | Note |
|---|---|---|---|


## Phase 3: New scatters: stars, nebulae and the region query

The multi-star catalogue, the stars every nebula needs and the region query the whole sphere reads through.

### Multi-star and compact systems

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.312 | Naming for trees: path codes (Aa, Ab, Ba) with branch words at wide splits, and planet stems such as Pikkita Ba I |  | Split from GEN.129. |
| MAP.213 | Scene, maps and the JavaScript twin for trees: per-node zoom, log-radial overview and positions_at in JavaScript |  | Split from GEN.129. |
| MAP.255 | The browser twin orbitpositions.js draws a wide pair as a circle: move it onto the stored orbit tree (bug) |  |  |
| GEN.313 | Orbital update of built systems: advance_orbital_phases for nodes and the GEN.109 hooks |  | Split from GEN.129. |
| GEN.129 | Multi-star systems of up to twelve stars | GEN.312, MAP.213, GEN.313 | GitHub issue [#777](https://github.com/dwhagar/planetGen/issues/777). Build after the GEN.128 design. |
| GEN.593 | Pulsar planets: the disc range stops short of PSR B1257+12's planets, and two of the three formation paths (a captured giant, an ablated companion) leave a pulsar with no planets |  |  |
| GEN.263 | Promote every bound compact pair to the scatter, whatever its mass |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 4. Core. |
| GEN.315 | A fill builds a system from its stored members and nodes; the system pass keeps only planets, moons, belts and disks |  | Split from GEN.264. |
| GEN.631 | GEN.316 part (b), narrowed: stages D and G draw their members from draw_member_tree, trees of 3 or more stars build as member stars, and the old companion draw and the two-star probability table are deleted |  | Core 3; takes an epoch at merge |
| ADM.128 | The admin one-off system page and the API schema offer the star count and the cap (stars and max_stars) | GEN.631 |  |
| GEN.316 | Stages D and G draw a member tree for every object the fill creates, and the two-star classes are retired | GEN.315, GEN.631, GEN.332, GEN.640 | Split from GEN.264. |
| GEN.264 | The fill builds a system from its stored members; the system pass keeps only planets, moons and belts | GEN.315, GEN.316 | Multi-star systems at scatter time, Boss 2026-10-15. Step 5. Core. |
| GEN.270 | Close-binary types: contact, Algol, spectroscopic and eclipsing flags, tertiary shares, merger products and Be stars |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 14. Core. |
| GEN.271 | White-dwarf binary types: close double white dwarfs, post-common-envelope pairs, cataclysmic variables, AM CVn, symbiotic and barium stars |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 15. Core. |
| GEN.272 | Accreting compact binaries: X-ray binaries, accretion luminosity, microquasars, spider pulsars and compact triples |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 16. Core. |
| GEN.327 | The group node and the remnant draw: offset and velocity columns, a decay age, and the lightest member ejected |  | Split from GEN.273. |
| GEN.331 | Substellar companions: brown-dwarf close (the desert) and wide companions, brown-dwarf pairs and planetary-mass companions |  | Split from GEN.275. |
| GEN.332 | Circumbinary and S-type planets, circumbinary disks and rogue planet pairs at fill | GEN.315, GEN.631 | Split from GEN.275. |
| GEN.640 | GEN.316 part (c): delete BinaryStarProxy, WideBinaryPair and the writers of the binary_* columns, once circumbinary and S-type planets come from the tree | GEN.332 | Split from GEN.631. |
| GEN.275 | Substellar and circumbinary types: brown-dwarf and planetary-mass companions, circumbinary disks and planets | GEN.331, GEN.332 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 19. Core. |

### Nebula and young stars

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.357 | Generation knows the cloud it sits in: surrounding_cloud extended, containers found early, StarSystem.apply_cloud |  | Split from GEN.95. |

### Region query

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.61 | Nearest-first by a radius ladder, and a dense-core guard with an early-stop count |  | Spherical region query, stage C. Step 4. Core. |
| GEN.376 | Region test uses the haversine form for angles under 1e-3 rad in `galaxy/region.py` |  | GEN.372 audit. |

### Superbubble density design

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.391 | Superbubbles shape the stars: the design of the bubble field as a density modulation, with its bounds and touch points checked |  | Boss 9:27 PM EDT: step 1 of the superbubble density chain; Core 3 proposed. |

### Star composition

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.464 | The disk solids law reads the star's metal mass fraction Z over Z0, not 10 to the power of [Fe/H], which is off by 1.7 to 3 times in metal-poor alpha-enhanced stars |  | Body composition research lane finding F6 (A-S2); generator epoch bump. Core 4 proposed. |

## Phase 4: New scatters: components, clusters and audits

Satellites, clusters and wanderers as components, the globular cluster chain and the audits that check what the scatters make.

### Components and wanderers

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.450 | The star-formation profile's radial slope: reconcile the Schmidt index with the observed gradient, with the choice named |  | Core 4 citation re-check, 2026-10-17; follow-up of GEN.134. Core 4. |
| GEN.371 | Feasibility study: generate the satellite galaxies through the scatter, making only the region each one occupies |  | Core 1; Boss 4:02 PM EDT: feasibility study. |
| GEN.535 | Check the Magellanic gas fraction and supernova-remnant count against the observed Clouds, and the dwarfs that still come out with gas (Leo I) |  | Open inferences of GEN.528: gas fraction, low LMC remnant count, Leo I gas. |
| GEN.534 | A satellite's young stars reach the system-type tables: the young-star system types name `satellite_young` |  | Residual of GEN.528: the system-type tables do not name satellite_young. Core 3 proposed. |
| GEN.328 | Groups at fill in star-forming cores and OB regions | GEN.327 | Split from GEN.273. |
| GEN.273 | Young groups and trapezia: the `group` node, offset and velocity columns, a decay age and the remnant draw | GEN.327, GEN.328 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 17. Core. |
| GEN.365 | Hot white dwarf and planetary-nebula central star radius disagrees with the radius that its luminosity and temperature imply (bug) |  | Finding from Core 3. |

### Globular clusters

| ID | Item | Needs | Note |
|---|---|---|---|

### Whole-sphere audits

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.161 | Fill-stage and system-pass tests: stream independence, same seed same system, and the probability table |  | Phenomena catalogue, Boss 2026-10-15. Step 2 tests. Quality. |
| TEST.162 | Multiplicity audit: companion states, row bookkeeping, stability, determinism, no planets in a scatter, fill reuse, motion and the census anchors | GEN.315 | Multi-star systems at scatter time, Boss 2026-10-15. Step 9. Quality. |
| TEST.164 | Catalogue type audit: type counts, region gates, the N tail, interaction consistency and group decay | GEN.270, GEN.271, GEN.272, GEN.275, GEN.327 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 20. Quality. |
| TEST.165 | Re-measure the region query on a real generated galaxy, warm and cold |  | Spherical region query, stage C. Step 6. Quality. |

### Composition proofs

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.211 | A star composition sweep: every kind of star, remnant and brown dwarf, fractions, trends and determinism |  | Boss 11:57 PM EDT; the proof of the star composition. |
| TEST.215 | The Solar System calibration rows live in test files only: a flag fences them and a test fails when a generator table reads one | TEST.211 | Boss 1:44 AM EDT, Stand-in only: R5, calibration rows fenced into tests. |
| TEST.214 | Every registered possibility table: for every class the gated weights add up to 100 percent, no class has an empty gate, and every row has a source and a licence |  | Boss 1:40 AM EDT: the gated weights add up to 100 percent; registry walk, so later tables are covered. Core 4 proposed. |

## Phase 5: New scatters: core, halo, interface update and docs

The realistic galactic core and the halo tier, the interface update for the spherical grid, the headline item and the docs that follow.

### Galactic core and halo

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.589 | Research: the central molecular zone's upper cluster mass and its cloud-shock disruption, so the young clusters of the core match the few that are observed |  | Finding from Core 3; research. |
| GEN.377 | Scatter limits by place for the neutron star and black hole scatter below the cut: the core and halo tiers apply to it, not one global mass cut |  | Core 3 finding after GEN.333. |
| GEN.303 | Per-cell system cap as a configuration option: `generation.cell_cap`, default 5,000, lowest 500, highest 50,000 |  | Core 3 wires it; Boss 1:47 PM EDT. |

### Interface update

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.115 | Interface update: every page and control works with the spherical grid, 96-bit IDs, the new phenomena and the region queries (stage after the sphere and generation stages) | NAV.61 | Interface; queued after sphere stages (Boss, 12:52 AM EDT) |
| TEST.168 | Control sweep: a browser test drives every control on every page against a spherical galaxy with 96-bit IDs | UX.115 | Quality; queued after UX.115 |

### Habitability system

Boss, 5:38 PM EDT 2026-10-23: the habitability system is finished before 8.1 is ready, five lanes now (this overrides the 8.1-only pause for these items). The model and score are built (GEN.83 to GEN.89, GEN.148, GEN.177, GEN.178, GEN.207, GEN.208, GEN.213; GEN.146 Kopparapu zones). Lane A symbols UX.101; lane B lists UX.100; lane C map MAP.143, MAP.219, MAP.220; lane D model GEN.639; lane E help text, browser tests and the radiation link question (GEN.445 stays in 9.0 until Boss decides).

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.157 | Every sprite icon draws as a solid black shape because the site's Content-Security-Policy refuses the sprite's own style element (bug) |  | Boss bug 9:07 PM EDT |
| UX.100 | Planet and moon lists show habitability as symbols, not notes |  |  |
| MAP.143 | Color sectors by their number of habitable locations |  | GitHub issue [#717](https://github.com/dwhagar/planetGen/issues/717). |
| MAP.219 | Galaxy Map marker layer: a kind registry with per-kind toggles, the zoom budget and picking, shipping the black-hole kind |  | Split from MAP.132. |
| MAP.220 | Markers for habitable worlds and chosen classes, optional per kind | MAP.219 | Split from MAP.132. |
| GEN.639 | The planet-mass habitable zone limits and the pre-main-sequence runaway history left open by GEN.146 (design note 4.1) |  | From Blocker 3's report of GEN.146 (#2685). |
| GEN.643 | Habitable means PHI-6 overall Blue or Green: one rule, one function, read by the sector counts, the markers, the lists, the search and the descriptions |  |  |
| GEN.649 | The four life scores: Microbial, Flora, Fauna and Sentient, each 0 (none) to 9 (abundant), from the evolutionary stage and the surface conditions; human operability is removed |  |  |
| GEN.666 | The life text on a planet page contradicts itself: one paragraph says life would evolve at a normal pace, the next says hyper-accelerated (bug) |  | Boss bug 9:11 PM EDT |
| GEN.653 | Suitability's food rung reads the Flora and Fauna scores: Flora is always good food, both are better, Fauna includes all animal life | GEN.649 |  |
| GEN.654 | Aquifers: a world with a dry, frozen or too hot surface can hold liquid water underground, within reach of the 273 K level; the depth feeds the Water and Suitability rungs and the Microbial score only | GEN.649 |  |
| GEN.655 | Survival equipment: one physics/equipment.py returns the list of needs (required or recommended) for the four slots, and replaces habitability.equipment and protection.requirements |  | Survival equipment system, Boss 7:56 PM EDT 2026-10-23; design draft section 11. |
| DB.93 | Planets and moons store the equipment levels and regional temperatures and drop equipment_tier; the API serves one equipment object | GEN.655 | Survival equipment system, Boss 7:56 PM EDT 2026-10-23; design draft section 11. |
| GEN.656 | Regional temperatures (Polar, Northern, Equatorial, Southern; day side, terminator and night side for locked worlds) and the day and night swing, with the regional sentences of the detail text | GEN.655 | Survival equipment system, Boss 7:56 PM EDT 2026-10-23; design draft section 11. |
| UX.153 | Search offers a Needs tick per equipment item, the habitability levels help page and the describer sentences read the needs list, and no sentence says none | GEN.655, DB.93 | Survival equipment system, Boss 7:56 PM EDT 2026-10-23; design draft section 11. |
| UX.154 | Lists and the system page show the equipment lists as chips: solid for required, outlined for recommended | GEN.655, DB.93, UX.153 | Survival equipment system, Boss 7:56 PM EDT 2026-10-23; design draft section 11. |
| POP.15 | A per-item technology level table for the species fit: which equipment item a species can field at which tech level | GEN.655 | Survival equipment system, Boss 7:56 PM EDT 2026-10-23; design draft section 11. |
| DB.92 | Planets and moons store the depth to liquid water and the thickness of the water-bearing layer (aquifers) | GEN.654 |  |
| DB.91 | Planets and moons store the four life scores and drop phi_bio, phi_cpx and phi_tech | GEN.649 |  |
| UX.152 | Planet and moon pages, lists, the system page and the API show the four life scores in the cluster; the old Microbial, Complex life and Human operability line goes | DB.91, GEN.649 |  |
| UX.155 | The /lists planet and moon pages show the PHI-6 symbols in a Habitability column and the API fields carry them (UX.100 part 2) |  | Habitability lane B sub-item, requested by lane B 8:08 PM EDT 2026-10-23 (relayed). |
| UX.156 | The planet and moon item pages show the PHI-6 symbols (UX.100 part 3) |  | Habitability lane B sub-item, requested by lane B 8:08 PM EDT 2026-10-23 (relayed). |
| GEN.644 | Proposal: tie the PHI-6 into the planet classification (class letters, a habitability tier, or a suffix) and the rule for where a class and a PHI-6 colour must agree | GEN.643 |  |
| GEN.645 | Design: condition-based planet and moon classes: the rule table, the inputs, the legacy letters and what happens to GEN.28, GEN.90 and GEN.91 | GEN.644, GEN.643 |  |
| DOC.78 | One end-to-end draw-order table from the star to the PHI-6 for the contracts file: blocks, streams and the reproducibility rule |  |  |
| GEN.650 | Design: the planet mass law once the class no longer supplies the mass: the isolation mass at the formation distance, the snow line and the host, replacing the class mass ranges of decision 28 (design note, prerequisite of GEN.646) | GEN.645 |  |
| DOC.75 | Help and docs text for the whole habitability system: the score, the symbols, the lists, the sector colouring and the map markers | UX.100, MAP.143, MAP.220, GEN.639 |  |
| DOC.77 | Docs for the four life scores: the help text, docs/design/habitability/habitability-index.md and the database schema document | GEN.649, UX.152 |  |
| TEST.300 | End-to-end browser tests for the habitability system: symbols on lists, sector colouring and habitable-world markers | UX.100, MAP.143, MAP.220 |  |

### Radiation chain for habitability (GEN.445)

Boss, 5:42 PM EDT 2026-10-23: build the entire GEN.445 chain in 8.1; testing mode, not production, so new generator epochs, database versions and API versions are fine. These 20 items moved here from 8.4 and 9.0 in dependency order. The spine (lane I, to be started): GEN.565, GEN.557, GEN.428, GEN.429, GEN.452 (part 2, PR #2613, held until now for the 8.4 epoch; it now rides the spine's first epoch), GEN.411. Side items that feed the spine: GEN.553 (lane A, after UX.101), DB.73 (lane B, after UX.100), GEN.412 (lane D, after GEN.639). Lane F, radiation constants and ozone: GEN.453, GEN.413. Lane G, bands: GEN.417, GEN.414, GEN.443. Lane H, surface dose: GEN.459, GEN.423, GEN.461, GEN.441, GEN.444. Lane E: GEN.445 with the help text and tests.

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.553 | One star-like body draw: origin, birth, state and [Fe/H], with the initial-final relation, the remnant mass, radius and age rules in one place |  |  |
| GEN.652 | The two main-sequence tables disagree: evolve_star reads mass-only piecewise laws and star_makeup reads MIST; one becomes the source of the other (bug) | GEN.553 |  |
| GEN.662 | Cluster and satellite stars draw a halo [Fe/H] for their lifetimes (stand-in, bug) | GEN.651 | Stand-in until GEN.651 lands |
| GEN.661 | Bright-star counts still use the reference-composition fraction (stand-in, bug) | GEN.651 | Stand-in until GEN.651 lands |
| GEN.651 | Wire the star's [Fe/H] into the star evolution: birth [Fe/H] is drawn before the state and passed to the lifetime and the giant branch, and the composition function is read by more than the white-dwarf flag | GEN.553 |  |
| DB.73 | Storage for one orbited-body model: a shared parent link for every kind, free parents, the orbit node table filled, moons storable under a rogue |  |  |
| DB.94 | The object index, resolve_object and the search panels read the moons of free planets | DB.73 | Habitability lane B sub-item, requested by lane B 8:08 PM EDT 2026-10-23 (relayed). |
| DB.95 | The sector_objects view includes the moons of free planets | DB.73 | Habitability lane B sub-item, requested by lane B 8:08 PM EDT 2026-10-23 (relayed). |
| GEN.565 | The mass object core and its three family classes: StellarObject, PlanetaryObject and SmallBody, with state and parent as data, replace the nine kind classes |  |  |
| GEN.412 | Generation records the albedo and the greenhouse multiplier it drew, and the temperature is exactly what they give |  | Boss 11:24 PM EDT; the inverse feeds the migration. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.557 | Wanderers: a rogue planet becomes a planet with no parent, drawn by the one body chain with an origin table, a notional formation host, the population's age and motion | GEN.553, DB.73, GEN.565 |  |
| GEN.428 | Planets and moons draw their bulk composition fractions instead of carrying a fixed text |  | Boss 11:49 PM EDT; generator epoch bump. Core 4 proposed. |
| GEN.429 | Bulk density, mass, radius and gravity are derived from the composition fractions with the mass-radius relation | GEN.428 | Boss 11:49 PM EDT; the atmosphere chain reads this gravity. Core 4 proposed. |
| GEN.459 | One radiogenic module: heat from the isotope-wise decay of potassium, thorium, uranium and aluminium-26, drawn from the host star's abundances, replaces the bulk-silicate values that are 37 to 50 percent too high and the separate ground-dose draw (bug) | GEN.428 | Body composition research lane findings F1 and F2, cluster 4; generator epoch bump. Core 4 proposed. |
| GEN.452 | The condensation table is replaced by the checked saturation curves: the hydrogen sulphide row is a copy of the sulphur dioxide row, and carbon dioxide, sulphur dioxide and the cold gases are 2 to 20 K off (bug) |  | Atmosphere research lane finding 1 (part 01 section 7); generator epoch bump. Core 4 proposed. |
| GEN.411 | The atmosphere's molar density comes from its per-gas partial pressures, so mixture, pressure and density always agree | GEN.429, GEN.452 | Boss 11:24 PM EDT; generator epoch bump. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.461 | The hydrosphere reads the stored water mass fraction: the compressed water-layer depth, the base pressure, the high-pressure ice and the liquidus are recomputed, and the class label follows from the water share (bug) | GEN.428, GEN.429 | Body composition research lane finding F4, cluster 8 (parts C, D, E); generator epoch bump. Core 4 proposed. |
| GEN.423 | Atmosphere density, scale height and the greenhouse multiplier follow from the mixture | GEN.411, GEN.412 | Composition thread finding; Boss 11:24 PM EDT chain. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.414 | Haze, dust and cloud: derived or drawn optical depths for every planet and moon | GEN.411, GEN.412 | Boss 11:24 PM EDT; draws and derivations, nothing fixed per class. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.441 | Surface and crust composition: the minerals, ices and ocean on the ground, from the bulk fractions, the temperature and the air | GEN.428, GEN.429, GEN.459, GEN.461 | Boss 11:58 PM EDT; the ground's composition; generator epoch bump. |
| GEN.443 | Atmospheric radiation: transmission by band and the greenhouse effect from the stored mixture, the aerosols and the cloud | GEN.414 | Boss 11:58 PM EDT; band transmission and greenhouse from the mixture. |
| GEN.444 | Surface radiation: the ground-level dose from the air's shielding and the ground's own radiogenic content | GEN.443, GEN.441, GEN.459 | Boss 11:58 PM EDT; replaces the proxies of physics/radiation.py; generator epoch bump. |
| GEN.445 | Habitability, toxicity and the UV index read the derived radiation instead of the proxies | GEN.444 | Boss 11:58 PM EDT; habitability reads the derived values. |

### Docs

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.19 | After UX.102: document the log output standard for the CLI, the Generate page and admin logs | UX.102 |  |
| DOC.72 | The README lists every external file the program needs, with its URL and how to install it, and a check keeps the list complete |  |  |
| DOC.73 | INSTALL.md recommends a larger innodb_redo_log_capacity for a machine that runs the full test gate |  |  |
| DOC.20 | After PERF.86: rewrite the queue and scatter notes and the queue docstrings | PERF.86 |  |
| DOC.48 | Docs: the reach, the halo, the profile rules and the Wolf-Rayet and nebula-contents rules in the design notes |  | Whole-sphere build train, Boss 2026-10-15. Step 14. Quality. |
| DOC.50 | Docs: the schema, the multi-star note and the catalogue say where the members live | GEN.264, TEST.164 | Multi-star systems at scatter time, Boss 2026-10-15. Step 10. Quality. |
| DOC.51 | Docs: the cell key in the coordinate system and schema docs, the region note moves to built, the API states the cap and the ladder |  | Spherical region query, stage C. Step 8. Quality. |
| DOC.18 | Last pass: re-run the documentation sweep when the other DOC items are done and Phase 1 is complete | DOC.19, DOC.20, DOC.48, DOC.50, DOC.51 |  |

### Globular clusters

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.161 | Bright-first fill for cluster sectors |  | Globular-cluster chain. |
| GEN.326 | Blue stragglers in clusters: N_BSS scales as the core mass to the 0.4 | GEN.161 | Split from GEN.162. |
| GEN.162 | Planet cull and blue stragglers in clusters | GEN.326 | Globular-cluster chain. |
| GEN.163 | Type-B pulsar planets in globular clusters (GEN.130 follow-on) |  | Globular-cluster chain. |
| GEN.330 | Cluster-hosted compact binaries: the cluster rates replace the field fractions | GEN.163 | Split from GEN.274. |
| GEN.274 | Cluster-hosted compact binaries and runaway bookkeeping: `runaway_origin` and the compact companion left by a supernova | GEN.330 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 18. Core. |
| GEN.297 | Calibrate the cluster generator against the Harris file: its distributions guide the laws, and no real cluster is placed in any galaxy |  | The file arrived (Boss, 3:04 PM EDT); Core 3, cluster line. |
| GEN.366 | Check the reference satellite table against McConnachie 2012: it only tests that the derived population looks plausible |  | Finding from Core 3. |

### Galaxy types groundwork

| ID | Item | Needs | Note |
|---|---|---|---|

### Superbubble density

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.392 | A bubble field handle on the galaxy shape: seed-carried bubbles with a box cut, a certified maximum factor and a bubble-free density read | GEN.391 | Step 2 of the superbubble density chain; Core 3 proposed. |
| GEN.393 | The density terms apply the bubble factors at one call site per form, with a golden test that no bubble means identical density | GEN.392 | Step 3 of the superbubble density chain; waits for Core 3's profile work so the files are not edited twice. |
| GEN.394 | The scatter bounds and the phenomena pass respect the bubble factor: a certified bound in band_majorants and per-position thinning | GEN.393 | Step 4 of the superbubble density chain; Core 3 proposed. |
| GEN.395 | Sector build follows the bubbles: the fill count from the cell mean, star positions weighted by the factor, shown in the profile and the region query | GEN.393 | Step 5 of the superbubble density chain; Core 3 proposed. |
| TEST.206 | Superbubble density tests: rim excess, interior depletion, old stars unchanged, scatter against fill and determinism | GEN.394, GEN.395 | Step 6 of the superbubble density chain. |
| DOC.61 | Docs: how the bubbles shape the stars, the setting values and what the drawn layer now is | GEN.395 | Superbubble density chain. |

### Galaxy laws build

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.228 | test_stars_and_phenomena_fit_within_their_sectors_real_cells compares a dark cloud's centre with its sector's cell, but the cloud is filed by its centroid |  |  |
| GEN.536 | Research the velocity dispersion of molecular gas (CO) in nearby spirals to set the young stars' scale height |  | Core 1: the young scale height is 116 pc against about 50 observed, because the gas dispersion is the HI value. Research. |
| GEN.405 | Metallicity and rates by position from the galaxy's own model: the mass-metallicity relation and the local rates replace the fixed gradient and the quoted-at-the-Sun values |  | GEN.368 build, Boss approved the table 10:11 PM EDT; Core 1 proposed for GEN.400 to GEN.403, Core 3 or Core 4 for GEN.404 and GEN.405 (profile.py, activity.py), any lane for GEN.406. Core 4 proposed (coordinator 10:16 PM EDT). |
| GEN.531 | Boss's research requests from the galaxy laws research (R1 to R15): the observed stellar halo, the dwarf range, the thick disk fits and the cluster system |  | Research requests R1 to R15 of the GEN.406 note; does not gate GEN.301. |

## Phase 6: Bulk composition and the albedo records

The composition fractions that drive density, mass and gravity, stored and described, and the recording of the albedo and greenhouse values. The research, the designs and the data tables they read (the gas, materials, element and stellar tables) are in 8.1 phases 1 and 2 (Boss, 2026-10-16, 11:58 PM EDT: the building blocks of the composition headline go down to 8.1).

Moved into 8.1 from 8.4 phase 1 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Bulk composition

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.75 | Object IDs stay stable: a refill that changes a system's planet or moon count keeps every surviving body's ID, and a body that changes parent keeps its ID |  |  |
| GEN.566 | Research: the size-frequency law of small bodies, for belt members drawn on demand and the comet size cut, beside the interstellar size law of GEN.560 |  |  |
| GEN.659 | The interstellar debris density stand-in of 1e12 per pc^3 is 700 to 1,800 times too low: the measured value is 7e14 to 1.8e15 per pc^3 (bug) |  | Found by lane J, 8:08 PM EDT 2026-10-23. |
| DB.74 | One generic insert and load driven by the column list, one registry of kind tables, and the shared core table with three profile tables | DB.73, GEN.565 |  |
| GEN.570 | Field declaration: every unit-valued field of the mass object classes declares its kind and native unit, and the one field table serves the serializer | GEN.565 |  |
| GEN.575 | Sweep the mass object classes: the 35 unit-suffixed attributes and every quantity held under two units get one declared native unit | GEN.565, GEN.570 |  |
| GEN.584 | The core's kinematic record: velocity, spin axis and rate, mass and epoch in the point, and the derived columns dropped | GEN.578, GEN.565, DB.73 |  |
| TEST.244 | Tests: re-primary and reparent invariants, frame round trips, momentum conservation, a moon's velocity as the sum of its hops, a locked body's rate follows a primary change | GEN.576, GEN.578, GEN.581, GEN.583 |  |
| GEN.583 | reparent: the state-preserving change of a body's parent, shared by re-primary at an event, an ejection, a capture and a merger | GEN.578 |  |
| GEN.582 | Origin and present-state split: fields marked origin or state in the Field table, the draw record held until validation finishes, and the zone and class re-run from the seed | GEN.565, GEN.570 |  |
| GEN.581 | primary_first(members): before a system is validated the heaviest branch becomes the primary and the others are recalculated, one function for the tree, the fill, moons and free pairs | GEN.315 |  |
| DOC.71 | Docs for the System rule: database-schema.md, the system view in docs/design and the header of system_tree.py | GEN.576, DB.78, GEN.578 |  |
| TEST.242 | Tests: invariants for any system at any time: zero mass-weighted offset sum, add and remove returns exactly, a parent change keeps every absolute position, one member equals the binary formula bit for bit | GEN.576 |  |
| GEN.580 | Research: the binary and satellite fraction of small bodies by size, and the survival of sub-satellites of moons, to replace the stand-in cut of the hosting gate |  |  |
| GEN.579 | The hosting gate: can_host(body) in three parts and the gated table of the satellite share by size class, as labelled stand-in rows |  |  |
| MAP.253 | Map and scene read the parent chain: systemscene, body_positions, the browser twin orbitpositions.js and the systemmap branches drop the three origin conventions | GEN.576, DB.78, GEN.577 |  |
| GEN.578 | The point class gets the local frame and the barycentre system frame, and the system anchor is no longer the nearest star | GEN.576 |  |
| GEN.577 | Stored planet and moon states become the state of the barycentre of the body and its own members | GEN.576 |  |
| DB.78 | Central-object offsets replace the wobble: the reflex_offset columns, calculate_reflex_offset, the SQL mirror and the nine binary position columns go | GEN.576, DB.73 |  |
| GEN.576 | system(host): the one function that groups bodies by their parent link and returns the barycentre sums, the central object's offset and the valid-system check | GEN.565, GEN.315 |  |
| GEN.587 | Replace the hand-written relation sites with the relations module, by file group, bit-identical on the golden set |  |  |
| GEN.571 | Replace the ad-hoc unit conversions with the registry, by file group, bit-identical on the golden set |  |  |
| GEN.574 | Describe and API: qty() builds a Quantity and the {kind, si} shape becomes {kind, value, unit} |  |  |
| DB.86 | Six columns in the column registry are declared 'unspecified': atm_density, atm_molar_density, formed_at, parabolic_mean_anomaly and two more |  |  |
| DOC.70 | Docs for the unit registry: units-and-number-formatting.md and database-schema.md updated, the unit rules in the registry's file header |  |  |
| GEN.564 | A refill job for a galaxy under an older generator epoch: clear the fill-level rows of every sector and rebuild them without a rescatter | DB.75 |  |
| GEN.555 | The pair builder for any masses: planet pairs, brown dwarf pairs and a planet with its moons share one stability test | GEN.315, GEN.331, GEN.332, GEN.581 |  |
| GEN.554 | Remnant counts follow from the mass function and the star formation history, with a kick step, and Boss's research densities become test targets | GEN.553 |  |
| GEN.550 | Umbrella (placeholder): one recursive generation process for anything that can have something else orbiting it, from a planet with moons up to a group of star-like objects that each have a system |  |  |
| GEN.460 | Giant planets and brown dwarfs get their radius and heavy-element content from the checked laws: the heavy-element mass law, the Fortney and Baraffe tables with inflation, and the burning limits shifted by metallicity and helium (bug) | GEN.428 | Body composition research lane finding F7, cluster 6 (parts A and C); generator epoch bump. Core 4 proposed. |
| GEN.641 | Research: the interior of a giant planet, ice giant and brown dwarf: the layers (molecular, liquid and metallic hydrogen and helium, helium rain, a dilute core), whether there is a rocky or icy core, and what the core is made of |  |  |
| GEN.642 | Giants, ice giants and brown dwarfs get a full interior composition: the molecular envelope, the liquid and metallic hydrogen and helium layers, the dilute zone, and a rocky, iron or icy core with its own composition, drawn and described | GEN.641, GEN.460, GEN.428, GEN.429 |  |
| DB.89 | Giants store their interior layers: the layer boundaries, the helium share per layer, the core flag, core mass and core composition | GEN.642 |  |
| DB.69 | Planets and moons store the composition fractions | GEN.429, GEN.459, GEN.460 | Boss 11:49 PM EDT; migration revision after DB.67. Core 4 proposed. |
| GEN.430 | The planet and moon descriptions state the composition from the stored fractions | GEN.428, GEN.429 | Boss 11:49 PM EDT; the description side. Core 4 proposed. |

### Albedo and greenhouse records

| ID | Item | Needs | Note |
|---|---|---|---|

### Atmosphere tables

| ID | Item | Needs | Note |
|---|---|---|---|

## Phase 7: Surfaces, small bodies and one mixture

The ground's composition and the small bodies' composition from the bulk model, and the atmosphere's molar density derived from its partial pressures with the values that follow from it.

Moved into 8.1 from 8.4 phase 2 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Surface and small bodies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.442 | Asteroids, comets, rings and dwarf planets draw their composition from the formation distance and the snow line | GEN.428, GEN.429, GEN.459 | Boss 11:58 PM EDT; small bodies; generator epoch bump. |
| GEN.660 | The complex_organics row of the materials table has a density of 1.5 g/cm3 while Choukroun 2020, used by the comet draw, has about 2 g/cm3 (bug) |  | Found by lane K, 8:21 PM EDT 2026-10-23. |
| GEN.466 | Moons draw their origin type, ice fraction and composition from the parent, the circumplanetary ice line and the age | GEN.428, GEN.429, GEN.442 | Body composition research lane part D sections 2.8 to 2.11; generator epoch bump. Core 4 proposed. |
| GEN.559 | Unbound comets: a size law and a cut, composition from the comet table, derived activity and velocity | GEN.557, GEN.442 |  |
| GEN.558 | Moons of wanderers: ejection survivors and isolated-disk moons are table rows, and the fixed 20% boolean is gone | GEN.556, GEN.557, GEN.466 |  |
| GEN.556 | `populate(host)`: planets, belts, comets and moons come from the class tables by one recursive routine, and the system pass is absorbed | GEN.315, GEN.555, GEN.565, GEN.581 |  |
| GEN.455 | Airless bodies and comets draw an exosphere or coma from their surface composition, temperature and the star's wind: the sodium and potassium tails and the production rates | GEN.441, GEN.442 | Atmosphere research lane new work 6 (parts 04 section 5 and 05); generator epoch bump. Core 4 proposed. |

### Mixture and stored values

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.479 | The class atmosphere mixes contradict the checked chemistry: classes C, D and S have no mix, free oxygen sits in seven classes with no life rule, class L is 60 percent argon, E and F hold ammonia and methane at 340 to 400 K, N reads reducing and is oxidised, and I, J and T share one mix in every zone (bug) | GEN.411 | Atmosphere lane part 09 conflicts 3 to 8 with the code, checked on main; replaced by the recipe tables. Core 4 proposed. |
| GEN.632 | Two code paths still hold their own vapour-pressure constants: rogue_surface.nitrogen_vapor_pressure_pa and the water curve in hydrosphere.py should read the one saturation function | GEN.452 | Pool; rides the 8.4 epoch |
| DOC.74 | Design part 01, section 7.4: the CO2 row and its recommended CO2-2 sublimation fit are wrong in the cold; the table now says to use CO2-1 below 195.3 K (bug) |  | TODO thread |
| GEN.451 | The mantle redox mean follows the meteorite building blocks: a young rocky planet starts reduced, not at the modern Earth's 3.5 above the iron-wuestite buffer | GEN.428 | Atmosphere research lane finding, part 5, 2026-10-17; generator epoch bump. Core 4 proposed. |
| GEN.454 | Giant planets, ice giants and brown dwarfs draw their atmosphere mixture from the equilibrium chemistry sequence: by temperature, pressure, metallicity and quench | GEN.411, GEN.460 | Atmosphere research lane new work 5 (part 05 sections 8 and 9); generator epoch bump. Core 4 proposed. |
| GEN.422 | Stop the recalculation redraw and order the steps so pressure, mixture and temperature agree after the cold-side freeze-out | GEN.411, GEN.412, GEN.452 | Composition thread finding; Boss 11:24 PM EDT chain. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.594 | Each AD dataset (AD1 to AD8) ships its best possible fallback, used only when a run is told explicitly to disregard missing files |  |  |
| GEN.424 | The balance gas is a named gas with its own molar mass, and habitability reads and writes the mean molar mass of the stored mixture | GEN.411 | Composition thread finding; Boss 11:24 PM EDT chain. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.657 | Toxic-species rows in the equipment limits table, wired to the stored gases: required gas filter grade by IDLH, with every recalled limit checked against NIOSH and ACGIH | GEN.424, GEN.655 | Survival equipment system, Boss 7:56 PM EDT 2026-10-23; design draft section 11. |
| GEN.658 | Dust loading and dust toxicity feed the dust filter and seals level of the equipment needs | GEN.414, GEN.655 | Survival equipment system, Boss 7:56 PM EDT 2026-10-23; design draft section 11. |
| DB.67 | Planets and moons store the albedo, the greenhouse multiplier and the molar density they were built with | GEN.411, GEN.412, GEN.422, GEN.423, GEN.424 | Boss 11:24 PM EDT; migration revision after Core 1's 0102 to 0105 stack. Core 4 proposed (Boss, 11:41 PM EDT). |
| DB.71 | Planets and moons store the partial pressures of the added gases, not only the ten stored today | GEN.411, DB.67 | Atmosphere research lane finding 10 (part 05 section 13); migration revision after DB.67. Core 4 proposed. |

### The Quantity call-site sweep

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.128 | Map side panels and other browser call sites take Quantity |  | Split from GEN.300. Moved from 8.3 and rescoped onto GEN.569 and UX.140. |
| GEN.349 | API output and stored-column readers take Quantity; the grep test fails on a hand-written conversion factor; the old helpers are removed | UX.128 | Split from GEN.300. Moved from 8.3. |
| GEN.300 | Every quantity in the code moves onto the Quantity class, and hand conversions are removed | UX.128, GEN.349 | The call-site sweep, partly covered by GEN.571 and GEN.574. Moved from 8.3. |

## Phase 8: Ozone, haze, dust and cloud, and the stored shares

The ozone figures and the scatterers that the colour needs, each its own thread, and the columns for the surface and small-body shares.

Moved into 8.1 from 8.4 phase 3 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Ozone

| ID | Item | Needs | Note |
|---|---|---|---|

### Haze, dust and cloud

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.499 | The lava-world vapour species and mean molar mass from temperature and composition, replacing the heavy-oxide range of 40 to 60 g/mol |  | Body lane G-5 (part G3): lava-world vapour species and mean molar mass 25 to 33 g/mol. Core 4 proposed. |
| GEN.456 | Lava worlds and hot rocky planets draw a vapour atmosphere over the magma from the surface composition and temperature: sodium, potassium, silicon monoxide, magnesium, iron and the oxides | GEN.441, GEN.411, GEN.414, GEN.451, GEN.499 | Atmosphere research lane new work 8 (parts 05 and 06); generator epoch bump. Core 4 proposed. |

### Stored surface composition

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.70 | Planets, moons and small bodies store the surface and small-body composition shares | GEN.441, GEN.442, GEN.461, GEN.462, GEN.466 | Boss 11:58 PM EDT; migration revision, in merge order. |
| GEN.477 | The asteroid belt and field composition is drawn from the star, the distance and the belt's mass, not four names drawn uniformly from a list of 84 (bug) | GEN.442, DB.70 | Body composition research lane finding F11 (tables E1 to E4): belt and field composition drawn from the star, distance and belt mass; generator epoch bump. Core 4 proposed. |
| GEN.471 | Clathrate hydrates and eutectic brines for icy crusts and ocean moons: the dissociation pressures and the liquid range below 273 K | GEN.461 | Body composition research lane part F item F4 (later); generator epoch bump. Core 4 proposed. |
| GEN.462 | Ground albedo by band: the area-weighted colour of the ground, ice and ocean cover with space weathering, handed to the atmosphere colour once | GEN.441, GEN.442 | Body composition research lane cluster 9 (part E sections 2.6 and 2.7, part D taxon albedo); generator epoch bump. Core 4 proposed. |

## Phase 9: Stored scatterers, the text and the proofs

The columns for ozone and the aerosols, the atmosphere text written from the mixture, the consistency sweeps and the docs.

Moved into 8.1 from 8.4 phase 4 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Consistency

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.68 | Planets and moons store the ozone column and the haze, dust and cloud optical depths | GEN.414, DB.67 | Boss 11:24 PM EDT; migration only. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.646 | Class designation leaves generation: planets and moons are drawn from their zone, mass, formation distance and composition recipes with no class letter | GEN.645, GEN.428, GEN.429, GEN.411, GEN.412, GEN.442, GEN.454, GEN.461, GEN.466, DB.67, DB.69, DB.70, GEN.650 |  |
| GEN.647 | The class is assigned after the body is complete: classify(body) from the stored atmosphere, temperature, terrestrial composition and gas composition, for planets and moons | GEN.646 |  |
| GEN.648 | The generation request for a habitable world (+habitable_world) runs an ordered ladder of nine adaptations before it gives up: redraw the surface, restrict the planet, more starting water, aquifers, regenerate the planet, adapt the system, add a brighter primary, regenerate the star, then the best result | GEN.646, GEN.643 |  |
| TEST.301 | Tests for the condition-based classes: every planet and moon class from its conditions, class and PHI-6 agreement, reclassification after an edit, determinism | GEN.647, DB.90, GEN.648 |  |
| DB.90 | Planets and moons store the derived class and a reclassification step keeps it true after an edit or a rebuild | GEN.647 |  |
| UX.150 | Pages, lists, search, filters, the API and the map markers show and filter by the derived class and the PHI-6 colour; the class wording and the descriptions follow the rule table | DB.90, GEN.647 |  |
| GEN.415 | Every compound the planet descriptions mention comes from the gas table, and the atmosphere text is written from the stored mixture | GEN.411 | Boss 11:24 PM EDT; the description side of consistency. Core 4 proposed (Boss, 11:41 PM EDT). |
| TEST.209 | One consistency sweep over every planet and moon class: composition, density, temperature, ozone, aerosols and text all agree | GEN.411, GEN.412, GEN.414, GEN.415, DB.67, DB.68, GEN.422, GEN.423, GEN.424 | Boss 11:24 PM EDT; the proof. |
| DOC.76 | Docs for the condition-based classes: the rule table, the planet class pages and the help | GEN.645, DB.90 |  |
| TEST.210 | A bulk composition sweep over every class: fractions, density, mass, radius and gravity agree | GEN.428, GEN.429, GEN.430, DB.69, GEN.441, GEN.442, DB.70, GEN.459, GEN.460, GEN.461, GEN.462, GEN.466 | Boss 11:49 PM EDT; the proof. |
| TEST.240 | Round-trip and registry tests for the mass object core: every kind's stored row round-trips, and no kind table is named in a list the registry should own | GEN.565, DB.74 |  |
| PERF.136 | Measure the fill with wanderers and their moons, and decide whether the rogue cost is on by default | GEN.558, GEN.559 |  |
| TEST.239 | Tests for the unified chain: determinism by host, representation counts per kind and level, rogue and comet rates, derived remnant shares | GEN.554, GEN.558, GEN.559 |  |
| UX.138 | Pages, symbols and descriptions for the new rows: moons of wanderers, free brown dwarfs as stars, unbound comets | GEN.558, GEN.559 |  |
| DOC.63 | Docs: the atmosphere data model, the gas table and how the values agree | TEST.209 | Boss 11:24 PM EDT; docs only. |
| DOC.65 | Docs: bulk composition, the materials table and how density, mass and gravity follow | TEST.210 | Boss 11:49 PM EDT; docs only. |

## Phase 10: The colour model

The physics of the colour: the blackbody colour of the star, the Rayleigh optical depth of the mixture, the colour of a solid world's air and of a giant's cloud top, all as the absolute colour (not adapted to the star).

Moved into 8.1 from 9.0 phase 1 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Atmosphere colour

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.418 | Atmosphere colour of a solid world: the absolute colour of the air with ozone, haze, dust, cloud and the ground's bounce | GEN.414, DB.67, DB.68, GEN.441, GEN.462 | Boss 11:24 PM EDT; absolute colour, not eye-adapted. |
| GEN.668 | Atmosphere loss events (GEN.189) do not yet reduce the stored ozone column |  | Follow-up of GEN.413 |
| GEN.667 | Ozone reads the oxygen partial pressure through the mixture accessor once the child table lands | GEN.411 | Follow-up of GEN.413 |
| GEN.419 | Atmosphere colour of gas and ice giants, read at a reference level in the cloud tops | GEN.418 | Boss 11:24 PM EDT; reference level chosen in the design. |

## Phase 11: Composition and colour in the product

The colour shown in the descriptions, on the maps in place of the planet colour, the composition and radiation text of stars, surfaces and small bodies, and the documents of the colour.

Moved into 8.1 from 9.0 phase 3 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Descriptions and maps

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.420 | The descriptions show the atmosphere colour: a swatch and a word in the planet and moon facts, in every output | GEN.418, GEN.419 | Boss 11:24 PM EDT; one place serves every output. |
| MAP.244 | The maps and pages use the atmosphere colour where they showed the planet colour | GEN.420 | Boss 11:24 PM EDT; replaces the temperature gesture. |
| GEN.446 | The descriptions state the composition of stars, surfaces and small bodies, and the radiation at the surface | GEN.441, GEN.442, GEN.444, GEN.445 | Boss 11:58 PM EDT; stars, surfaces, small bodies and radiation in the text. |
| DOC.64 | Docs: the atmosphere colour, its method and its sources | GEN.420, MAP.244 | Boss 11:24 PM EDT; docs only. |

### Stand-in replacements

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.473 | Replace the asteroid taxon shares and the zone kernel, measured in the Solar System belt, by a derived or general distribution | GEN.442 | Boss 1:44 AM EDT, Stand-in only: replace the measured Solar System taxon shares by a derived or general distribution. |
| GEN.474 | Replace the Earth rock-type library, measured on Earth, by ground types derived from the bulk composition and the melting of the mantle | GEN.441 | Boss 1:44 AM EDT, Stand-in only: replace the Earth rock-type library by derived ground types. |
| GEN.475 | Meteorite groups become outputs and not inputs: regions of composition space reachable from the drawn bulk, with the measured groups kept as calibration | GEN.473, GEN.474, GEN.451 | Boss 1:44 AM EDT, Stand-in only: R3, meteorite groups as outputs of the drawn chemistry. |
| GEN.476 | Aluminium-26 and iron-60 follow from the star's formation environment, not from the Solar System's own value and range | GEN.459 | Boss 1:44 AM EDT, Stand-in only: R4, Al-26 and Fe-60 from the formation environment. |

### Weight research (atmosphere)

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.480 | Source the occurrence weights of bare rock, hydrogen-helium envelope and water-rich worlds and of the envelope metal fraction, by radius, mass and insolation (AR1) |  | Atmosphere lane research need R1 (rows b4, e3, e4, e5); replaces the seeded bracket of the DESIGN rows. |
| GEN.481 | Source the volatile delivery against formation distance: the water mass fraction by a_form and the snow-line law over time (AR2) |  | Atmosphere lane research need R2 (rows b1 to b6); replaces the seeded bracket of the DESIGN rows. |
| GEN.482 | Source the magma-ocean fO2 prior against mass and the oxidation law of the later mantle (AR3) |  | Atmosphere lane research need R3 (rows a1, a2, a3); replaces the seeded bracket of the DESIGN rows. |
| GEN.483 | Source the cosmic-shoreline constants and the pressure caps by retention regime (AR4) |  | Atmosphere lane research need R4 (rows d6); replaces the seeded bracket of the DESIGN rows. |
| GEN.484 | Source the giant planets' C/O distribution and the fraction of carbon-rich hosts (AR5) |  | Atmosphere lane research need R5 (rows f7); replaces the seeded bracket of the DESIGN rows. |
| GEN.485 | Source the hot-Jupiter thermochemistry above 2200 K, a Kzz distribution and the edges of the 1500 K and 2200 K regimes (AR6) |  | Atmosphere lane research need R6 (rows h2, h3, h4); replaces the seeded bracket of the DESIGN rows. |
| GEN.486 | Source the lava-world vapour regimes and the fraction of lava worlds that keep volatiles (AR7) |  | Atmosphere lane research need R7 (rows i6); replaces the seeded bracket of the DESIGN rows. |
| GEN.487 | Source the evolved and frost atmospheres: cold vapour pressures, the habitable-zone coefficients, the occurrence of Titan-, Pluto- and Venus-like bodies and the abiotic O2 frequency (AR8) |  | Atmosphere lane research need R8 (rows j1 to j5); replaces the seeded bracket of the DESIGN rows. |
| GEN.488 | Source the comet coma classes and their volatile ratios (AR9) |  | Atmosphere lane research need R9 (rows k1 to k5, k7); replaces the seeded bracket of the DESIGN rows. |
| GEN.489 | Source the white-dwarf and pulsar planet options (AR10) |  | Atmosphere lane research need R10 (rows a4, l1, l2, l3); replaces the seeded bracket of the DESIGN rows. |
| GEN.490 | Source the trace-gas presence weights (AR11) |  | Atmosphere lane research need R11 (rows m2 to m6, m11); replaces the seeded bracket of the DESIGN rows. |
| GEN.491 | Source the cloud and haze windows and yields (AR12) |  | Atmosphere lane research need R12 (rows n6, n7); replaces the seeded bracket of the DESIGN rows. |
| GEN.492 | Source the ground radioisotope abundances by host and the potassium depletion by body class (AR13) |  | Atmosphere lane research need R13 (rows o2, o3, o4, o6); replaces the seeded bracket of the DESIGN rows. |
| GEN.493 | Source the initial-rotation distribution and the saturation time against stellar mass for the host-star UV regime (AR14) |  | Atmosphere lane research need R14 (rows p6); replaces the seeded bracket of the DESIGN rows. |
| GEN.494 | Source the brown-dwarf type against temperature and the T/Y boundary near 500 to 700 K (AR15) |  | Atmosphere lane research need R15 (rows g1, g4); replaces the seeded bracket of the DESIGN rows. |
| GEN.501 | Source the activity factor of the NH3 photolysis rate from the Huebner and Mukherjee rate tables |  | Body lane reconciliation: NH3 lifetime 6.7e3 s quiet Sun, 5.6e3 s quiet/active mean; activity factor unknown. |

### Tuned draws and selection models

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.509 | Tuned draw builder B1: the radius valley and the envelope |  | Tunable gaps report: Tuned draw builder B1: the radius valley and the envelope. |
| GEN.510 | Tuned draw builder B2: delivery and water worlds |  | Tunable gaps report: Tuned draw builder B2: delivery and water worlds. |
| GEN.512 | Tuned draw builder B4: the iron branch |  | Tunable gaps report: Tuned draw builder B4: the iron branch. |
| GEN.515 | Tuned draw builder B7: ground and petrology |  | Tunable gaps report: Tuned draw builder B7: ground and petrology. |
| GEN.518 | Tuned draw builder B10: the chemistry module |  | Tunable gaps report: Tuned draw builder B10: the chemistry module. |
| GEN.665 | Two tuning targets are weaker than they look: the white-dwarf redox sample is 7 stars and the mass-radius cut barely filters |  | Target power |
| GEN.664 | The thin-disc offset rows disagree by 0.12 dex: Botelho solar twins (+0.08 dex) against the Hypatia thin-disc statistic (-0.044) |  | Scale choice |
| GEN.663 | A chemical-evolution birth-epoch model for thorium and uranium: row o4 draws the formation time from 3 to 13 Gyr |  | Stand-in until built |
| GEN.519 | Tuned draw builder B11: aerosol amplitude |  | Tunable gaps report: Tuned draw builder B11: aerosol amplitude. |
| GEN.520 | Tuned draw builder B12: the comet sampler |  | Tunable gaps report: Tuned draw builder B12: the comet sampler. |
| GEN.595 | Replace the Solar System comet population shape (GEN.520 D6) with a derived distribution |  | Boss, 2026-10-18 (fork 8, flagged stand-ins only). |
| GEN.521 | Tuned draw builder B13: belts |  | Tunable gaps report: Tuned draw builder B13: belts. |
| GEN.522 | Tuned draw builder B14: moons and Kuiper belt objects |  | Tunable gaps report: Tuned draw builder B14: moons and Kuiper belt objects. |
| GEN.523 | Tuned draw builder B15: host-star activity |  | Tunable gaps report: Tuned draw builder B15: host-star activity. |
| GEN.524 | Tuned draw builder B16: halogens and trace gases |  | Tunable gaps report: Tuned draw builder B16: halogens and trace gases. |

## Phase 12: The headline closes

The last audits and the headline: they test the component pass, the core and the halo, so they come after phases 4 and 5.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### The headline closes

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.175 | benchmark --suite accepts --keep but drops the scratch databases on success: make it honour --keep (bug) |  |  |
| PERF.174 | Every work unit's database connection runs the schema check again: skip it when the schema version is current |  |  |
| PERF.173 | A new-galaxy Generate run dies in stage 5 with "Execution not found in Redis": the worker warm outlives the 90 s execution record (bug) |  |  |
| PERF.155 | A unit-length ceiling: no unit of a stage runs longer than a set time (default 5 minutes) unless the stage is one unit; longer stages are split (F6) |  |  |
| TEST.261 | Release gate: a full bug hunt audit of the whole code base and its behaviour, done before every major or minor version change (8.1 first) |  | Standing release gate, Boss 2026-10-21 5:10 PM EDT; also gates every later major or minor version change, not build numbers (Boss, 8:01 PM EDT 2026-10-22). |
| SEC.33 | Release gate: a full security code audit of the whole code base, done before every major or minor version change (8.1 first) |  | Standing release gate, Boss 2026-10-21 5:10 PM EDT; also gates every later major or minor version change, not build numbers (Boss, 8:01 PM EDT 2026-10-22). |
| GEN.301 | Headline of 8.1: the new generate scatters, accurate galactic generation and scattering including the halo | GEN.450, SEC.33, TEST.261, UX.100, MAP.143, MAP.220, GEN.639, GEN.445, DOC.75, TEST.300 |  |
| GEN.637 | Re-run the GEN.625 1,000-seed all-Auto sweep on main after GEN.624 and GEN.403 (#2638) |  |  |
