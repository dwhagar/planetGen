# Version 8.1: Infinite zoom and the new generate scatters

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.1

## Goal

Two headlines. Infinite zoom: one camera flies from the whole galaxy to a star system (MAP.146). The new generate scatters: accurate galactic generation and scattering out to the 262 kpc reach, including the halo, the satellites, the nebula stars and a realistic core (GEN.301). The version also holds the Phase 1 tail: the queue steps, reproducible galaxies, the description layer, tech levels, directives and every open bug and flaky test. A third strand (Boss, 2026-10-16, 11:58 PM EDT): the building blocks of the 9.0 composition headline, which go down to 8.1 because they change no stored value: the research lanes and designs (phase 1, GEN.431 to GEN.437 with GEN.409, GEN.425 and GEN.426), the data tables (phase 2, GEN.438, GEN.410, GEN.427, GEN.439), the star composition function (phase 3, GEN.440) and its sweep (phase 4, TEST.211); the master note is docs/design/composition-model.md. The 8.1 stamp waits until all of it is done (REVISION_HOLD in scripts/bump_version.py).

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
| PERF.120 | Every command does one thing: the umbrella for reshaping every command and entry point |  | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| PERF.121 | `plan` does only the skeleton, and the scatters become one `scatter` command with exactly one kind per run (mass, luminosity, phenomena, nebula field, band) | PERF.120 | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| PERF.122 | The map-luminosity fill and the habitability re-score are their own commands, not `plan` flags | PERF.120 | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| PERF.123 | `galaxy` stops being a bundle: sector generation, the neighbour link, the backfill, the population pass and the settle are each their own command, and `galaxy` is a step list the wrapper builds | PERF.121 | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| PERF.124 | `sector` and `phenomenon` only generate: the population pass leaves `sector` and the settle leaves `phenomenon` | PERF.123 | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| PERF.125 | The population pass and the territories pass are separate commands | PERF.120 | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| PERF.128 | `benchmark` is a step list handed to the wrapper, not a sequence of its own | PERF.123 | Boss, 2026-10-17, 4:19 PM EDT (all commands). |
| GEN.309 | The backfill picks the draw path per block and gives the same stars wherever it lands |  | Split from GEN.42. |
| GEN.40 | Weed out sectors by star density before the bright-star backfill | GEN.309 | Parent; closes with its subitems. |
| PERF.35 | An interval or chunk ledger for untouched sectors once block-first backfill lands | GEN.309 | Research: ledger size at full galaxy. |
| PERF.99 | Work queue step 5: the Generate page, Queue page and job log read Redis; cancel now through the stop command | PERF.98, PERF.119 | PERF.86 step 5 |
| PERF.144 | Every job queued to RQ gets a time limit from the best available estimate of how long it takes, never no limit, so a killed worker's entry expires |  |  |
| PERF.145 | The run output repeats the stage: drop the stage-list line and the 'Stage 4 of 13: ...' line the engine logs under the wrapper's '=== Stage 4 of 13: ... ===' banner |  |  |
| GEN.600 | The updater and the downloader print the URL of every file they load |  |  |
| GEN.599 | A run with missing downloaded data files refuses unless it is told explicitly to disregard them, and then uses the best possible fallback and says so in its header |  |  |
| GEN.310 | Fill a whole block of the shell ladder: planetgen generate --block and the queue job | PERF.99 | Split from GEN.231. |
| UX.127 | The Galaxy Map block menu offers Generate again | GEN.310 | Split from GEN.231. |
| GEN.231 | Generate a whole block on the shell ladder: the generate --block mode and the block buttons come back | GEN.310, UX.127 | Maps hand-off 2026-10-15. Core, then Interface for the buttons. |
| DB.29 | update.sh and every database migration refuse to start while jobs are active, unless --yes is passed or the user confirms cancelling them | PERF.99 |  |
| PERF.100 | Work queue step 6: delete the old machinery and the work_* tables (Alembic migration) | DB.29 | PERF.86 step 6 |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.55 | Same seed, same data: a sector's contents depend only on the seed, the version and its address (internal) |  | Parent of the chain. |
| GEN.307 | The object-ID collision number keyed on the lower address |  | Split from GEN.57. |
| TEST.181 | The order-independence test: the same sectors at 1 and 4 workers and in two orders | GEN.307 | Split from GEN.57. |
| GEN.57 | A sector's contents depend only on the seed, the version and its address | GEN.307, TEST.181 | Name collisions go away with GEN.67; GEN.63 dropped. Also keys GEN.64's ID collision counter on address (it checks stored rows, PR #406). |
| GEN.179 | Store each sector's generation directive and attempt record with the sector |  | Left over from GEN.96. |
| DB.43 | Reed-Solomon codec in the standard library: planetgen/util/rscodec.py |  | Split from DB.9. Core 1 owns migrations. |
| DB.44 | Migration only: the sector storage checksum and the parity dirty flag |  | Split from DB.9. Core 1 owns migrations. |
| DB.45 | Storage checksum and the dirty flag: the export, the checksum and every writer | DB.44 | Split from DB.9. Core 1 owns migrations. |
| DB.46 | The parity file: groups, records, a delta update on every save and the stale-record rule | DB.43, DB.45 | Split from DB.9. Core 1 owns migrations. |
| DB.47 | planetgen repair-db: rebuild damaged sectors from the parity file | DB.46 | Split from DB.9. Core 1 owns migrations. |
| DB.48 | repair-db: rebuild the table when a clustered-index page is unreadable | DB.47 | Split from DB.9. Core 1 owns migrations. |
| DB.9 | Repair a damaged galaxy database from a parity file | DB.43, DB.44, DB.45, DB.46, DB.47, DB.48 | Boss 02:13Z: phase 1. Reed-Solomon parity over groups of sector exports; seed regeneration as fallback. |
| DB.17 | Repair by regenerating a damaged sector from its seed when parity cannot rebuild it | GEN.57, DB.47, GEN.179 | Research split of DB.9: the regenerate-from-seed fallback. |
| GEN.304 | The canonical fingerprint encoder: floats to 9 significant digits, SHA-256, structural and float digests |  | Split from GEN.136. |
| DB.49 | Migration only: the leaf digest column and the shell-band tree table |  | Split from GEN.136. Core 1 owns migrations. |
| GEN.305 | Store a leaf digest at save time and build the tree: roots and region digests | GEN.304, DB.49 | Split from GEN.136. |
| GEN.136 | Fingerprint tree: a stored leaf digest per sector and a shell-band tree of digests | GEN.304, DB.49, GEN.305 | Research follow-up to GEN.58 (built). |
| OPS.57 | generator_battery_digest(): a database-free digest over 6 to 8 canned cases | GEN.304 | Split from OPS.28. |
| OPS.58 | The version history row and the galaxy's internal record carry the epoch, the battery digest and the environment | OPS.57 | Split from OPS.28. |
| OPS.28 | Generator epoch and battery digest: say whether two checkouts generate the same galaxy | OPS.57, OPS.58 | Research: OPS.13 is built; this is the epoch it records. OPS.14, OPS.15, OPS.12 and TEST.77 read it. |
| OPS.15 | Each update says whether it changes generated output |  | Needs the fingerprint, so phase 2. |
| DB.16 | Store the generator epoch and run id on each sector instead of four version text columns |  | Research follow-up to DB.7 (built). |
| OPS.29 | Update reload: cover the gunicorn units and non-Apache hosts, and say what a reload aborts |  | Research follow-up to OPS.8 (built). |
| OPS.37 | A Generator version number: one sequential integer, shown in admin and the API |  | Boss 2026-10-09 20:59Z: done by the end of phase 1. Same number as OPS.28's epoch (decided). |
| TEST.77 | A golden-seed regression test | GEN.304, OPS.57, GEN.57 |  |

### Tech levels and species

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.30 | One migration for the per-species environment limits |  |  |
| POP.14 | A race is only placed where it can survive: a placement check and a consistency check over the data | DB.30 |  |

### Description layer

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.31 | Store structure, not prose, for life and flavour |  |  |

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

### Bugs and flaky tests

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.33 | MySQL 8.0.46 reflects the foreign keys fk_facilities_type and fk_facility_field_values_field as RESTRICT, so check-db calls a healthy database damaged (bug) |  | Quality run 2026-10-15, MySQL 8.0.46. Core. |
| DB.80 | `VALUES ROW(...)` key-run lists are rejected by MariaDB 10.11 although INSTALL.md supports MariaDB 10.4+ and MySQL 8.0.16+ (bug) |  |  |
| TEST.111 | test_ensure_sector_generated_creates_then_reuses_the_same_sector fails in a busy parallel run (bug) |  | Test flake reported 2026-10-09; Bugfixes lane 1. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.116 | test_bughunt_end_to_end stores M2V or M6V where it expects K2V on the MySQL 8.4 and MariaDB legs once (bug) |  | Bugfixes lane 1 report 03:57Z. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.126 | The heavy WebGL browser-map tests time out when four workers run them together (bug) |  | Lane 2 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
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
| TEST.260 | test_galaxy_absurd_radius_is_clean[200.0] hangs past its 10 s limit on main and leaves the run lock set, so every later test that starts a run fails (bug) |  |  |
| TEST.262 | test_a_double_click_flies_even_when_the_page_is_slow_to_answer_the_first_click fails once under load: the camera target stays at (3.56, -1.36, -0.29) instead of (0, 0, 0) (bug) |  |  |
| TEST.258 | Other tests that scatter a low shell of the toy shape may hold gigabytes after GEN.404: sweep them for resident memory |  |  |
| TEST.257 | test_sector_page_and_tiles_show_a_phenomenon_the_cli_added asserts an empty cloud list on an unseeded galaxy and sometimes gets a nebula (bug) |  |  |
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
| MAP.172 | Nebula cell aggregates in the region pyramid |  | Foundations lane 1. |
| MAP.174 | Nebula cover in the stage-view cell statistics, for filled and unfilled cells | MAP.172 | Foundations lane 2. |

### Composition research and design

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.465 | A second extraction pass over the sources the body and atmosphere research could not reach and that are reachable now |  | Body composition research lane master section 8; read-only. Unowned. |
| GEN.468 | Settle the water share beyond the snow line: the Marboeuf-derived 0.14 to 0.39 of the solids against the convention of about 0.5 |  | Body composition research lane master section 5 decision 2; Boss 1:21 AM EDT card Keep the switch; read-only. Unowned. |
| GEN.457 | The small bodies and moons composition design: the taxon and comet draws, compaction, the moon origin types and the ice fraction |  | Body composition research lane part D; design note, no code. Core 4 proposed. |
| GEN.560 | Research list for the chain: ejection survival of moons, circumplanetary disks of isolated planetary-mass objects, failed-supernova and direct-collapse shares, the size law of interstellar objects |  |  |

### Admin control suite

Boss (6:01 PM EDT, 2026-10-21): the admin pages, the Generate suite and the settings are one control suite, built in this order: Settings first, then the shell, job strip and palette, then Generate. Boss answered all 13 decision cards between 6:05 and 6:09 PM EDT (design note docs/design/admin-control-suite.md, section 7); his answer 1 puts all of it in 8.1 phase 1. Owner of every item: the admin redesign thread.

Shape dials (Boss, 6:41 to 6:46 PM EDT, 2026-10-21; goal issue #2355): the generator's shape inputs become five dials built from one lever table (GEN.601), so the Generate suite is not built on the retiring absolute inputs. The table, flags (GEN.602), Scale (GEN.603), Structure (GEN.604), Edge (GEN.605, form awaiting Boss), Halo (GEN.606), the unfilled Shape selector (GEN.608), range sources (GEN.609), the nuclear-amplitude fix (GEN.610), the pattern dial (GEN.611) and the stored pins (DB.81) sit here ahead of the Generate rows; the dials panel ADM.124 follows ADM.110; Census (GEN.607) waits on GEN.545 in the galaxy laws section. The default galaxy is also smaller (Boss, 6:57 PM EDT: still a spiral, no more than half the Milky Way's size; GEN.603 states the number): a new generator epoch, with lower system counts and cost. The Shape dropdown is filled in 9.4 by UX.134 and GEN.389. The lever-table items belong to the Galaxy dials lane (the coordinator, 6:55 PM EDT); Edge and Census have no owner until Boss's cards and GEN.545.

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.124 | One form-field look and error slot, piloted on the search, nearby and list-filter fields |  | Split from UX.49. |
| ADM.98 | Admin control suite research: how to pack many options into one quiet interface (admin pages, the Generate suite and the settings), a report and decision cards before anything is built |  | Admin control suite, Boss 6:01 PM EDT 2026-10-21: the admin redesign thread; research and cards first. |
| ADM.55 | Admin usability audit: walk every admin page and task, rank what is hard to use, and propose one layout |  | Audit narrowed to the pages other than Generate and Settings (key K8); UX.115 dropped. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.101 | Settings: words for every option: label, one-line summary, group, control and advanced flag in the settings model |  | Admin control suite key S1. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.102 | Settings: view classes: editable, owner, server and overridden, with the reason an option is fixed | ADM.101 | Admin control suite key S2. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| UX.125 | The Generate page and the one-off system page on the field macro | UX.124 | The Generate page conversion rides ADM.99. Admin control suite, Boss 6:01 PM EDT 2026-10-21: the admin redesign thread; research and cards first. |
| ADM.103 | Settings: the new page: grouped rows, a left rail, search, a switch that shows settings that cannot be changed (off by default), a sticky save bar, per-row reset and Details | ADM.101, ADM.102, UX.124 | Admin control suite key S3. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.104 | Settings: controls: switch, slider with a number box, duration, rate limit, chips and colour picker | ADM.103 | Admin control suite key S4. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.105 | Settings: the review page in plain language: 'Cache size: 1,024 MB to 2,048 MB', the restart list and the overridden list | ADM.103 | Admin control suite key S5. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.106 | Settings: the server-files disclosure and the overridden-option badge | ADM.102 | Admin control suite key S6. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| GEN.601 | Shape dials: the lever table, one registry row per slider (name, flag, label, unit, minimum, maximum, step, Auto default, the law it moves, range source) and the parity test |  | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| GEN.602 | Shape dials: the command-line flags built from the lever table, and the old absolute shape options retired | GEN.601 | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| GEN.603 | Shape dials: the Scale dial: mass override, compactness in sigma, density, and the any-two-of-three pin logic with the third shown | GEN.601 | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| GEN.604 | Shape dials: the Structure dial: flatness (scale length over scale height), bulge share and bulge size (scale radius over L) | GEN.601 | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| GEN.605 | Shape dials: the Edge dial: where the galaxy ends, with the resulting radius and the expected sector and system counts shown before a plan | GEN.601 | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| GEN.606 | Shape dials: the Halo dial: dark halo mass over stellar mass and the stellar halo share, with r200 shown and a size-halo warning | GEN.601 | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| GEN.608 | Shape dials: the Shape selector, added to the table and the Generate form unfilled: key `spiral`, flag `--shape`, one option for now | GEN.601 | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| GEN.609 | Shape dials: check the range source of every slider (size-law mass span, flatness, bulge share, bulge size, edge, halo ratios) and record each as verified, derived, bracket or placeholder |  | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| GEN.610 | The `--nuclear-amplitude` help text defines 1 as 'the observed' nuclear cluster and disc (2.5e7 and 1.5e9 solar masses), which are Milky Way values: replace them with a draw from the galaxy's own laws | GEN.601 | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| DB.81 | Shape dials: the plan stores which sliders were pinned and at what value, so the galaxy page can say 'mass pinned at X, the rest drawn' | GEN.601 | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| GEN.611 | Shape dials: the pattern dial: arm contrast (crest minus trough over crest plus trough, default 0.4), with arm count and pitch angle, replacing the two arm densities | GEN.601 | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| ADM.118 | Admin registry: pages, tasks and settings declare label, summary, keywords, role, danger level and URL, and the rail, hub and search read it | ADM.101 | Admin control suite key K1. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.92 | Admin shell: one layout, one navigation and the shared table, form, confirmation and job-log components | ADM.118, ADM.98 | Split from ADM.56. Key K2. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.119 | Job strip on every admin page: a thin strip that shows a running job, Cancel and the log | ADM.92 | Admin control suite key K3. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.120 | Command palette and /admin/search: server search over the registry, Ctrl+K and /, results for pages, tasks, settings, jobs and sectors | ADM.118, ADM.92 | Admin control suite key K4. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.122 | Admin keyboard set and a help overlay | ADM.92 | Admin control suite key K6. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.123 | Queue one next job: a single next slot that starts when the running job ends, with Cancel, shown in the job strip and History | ADM.119 | Admin control suite key K9. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.107 | Generate: the frame: status strip, pinned current-job panel (idle shows the last result), History tab and the Work Queue link |  | Admin control suite key G1. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.108 | Generate: the task flow engine: a Start screen with six cards, a step rail, the draft carried in the URL, a Review step with the estimate, summary and command line, and Run | ADM.107, UX.125, ADM.123 | Admin control suite key G2. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.109 | Generate: the Add sectors task: two choices, the four-tab place picker and how-big with a live count | ADM.108 | Admin control suite key G3. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.110 | Generate: the Start a new galaxy task: seed-decides shape fields, the first neighbourhood and a typed confirmation | ADM.108, GEN.601, GEN.602 | Admin control suite key G4. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.124 | Generate: the shape dials panel, built from the lever table: five dials of sliders with Auto, the pinned-versus-drawn read-out, any-two-of-three for mass, size and density, and the Shape dropdown | GEN.601, GEN.602, GEN.603, GEN.604, GEN.608, ADM.108 | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| ADM.111 | Generate: the Scatter task: rebuild, redo some, add a dimmer step, with one Advanced block of scatter limits | ADM.108 | Admin control suite key G5. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.112 | Generate: prevalence, star mix and overrides as one shared Advanced block, rendered once | ADM.108 | Admin control suite key G6. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.113 | Generate: the one-off system on the frame, with tri-state rows as segmented controls | ADM.107, UX.125 | Admin control suite key G7. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.114 | Generate: the Check and Wipe tasks and the History tab | ADM.107 | Admin control suite key G8. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.115 | Generate: Run again from History, and named recipes saved per admin in the control database | ADM.108, ADM.114 | Admin control suite key G9. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.116 | Generate: remove the old Generate page pieces and update docs/html-interface.md and the architecture table | ADM.109, ADM.110, ADM.111, ADM.112, ADM.113, ADM.114, ADM.124 | Admin control suite key G10. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.117 | Generate: browser tests for every task end to end, in both themes, without script and at phone width | ADM.109, ADM.110, ADM.111, ADM.112, ADM.113, ADM.114, ADM.124 | Admin control suite key G11. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.121 | Typed commands that open a task's Review step pre-filled and never run by themselves | ADM.120, ADM.108 | Admin control suite key K5. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.93 | Move every existing admin page onto the shell and fix the defects the audit ranks above its line | ADM.92, ADM.55 | Split from ADM.56. Key K7. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| UX.126 | The admin pages and the edit dialogs on the field macro | UX.124 | The admin page conversion rides ADM.56 and ADM.93. Admin control suite, Boss 6:01 PM EDT 2026-10-21: the admin redesign thread; research and cards first. |
| UX.49 | Form fields as Shoelace components (sl-input, sl-select, sl-checkbox) across the site | UX.124 | Follows UX.40 (done): the fields themselves become Shoelace components. |
| ADM.56 | Streamline and fix the admin pages from the audit: one admin shell, one navigation and one pattern for tables, forms, confirmations and job logs | ADM.92, ADM.93 | Admin control suite, Boss 6:01 PM EDT 2026-10-21: the admin redesign thread; research and cards first. |
| ADM.99 | Generate suite redesign: the Generate page rebuilt inside the admin control suite so a long list of options is easy to use | ADM.107, ADM.108, ADM.109, ADM.110, ADM.111, ADM.112, ADM.113, ADM.114, ADM.115, ADM.116, ADM.117, ADM.124 | Parent: closes with its children G1 to G11. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |
| ADM.100 | Settings page redesign: plain-language groups and descriptions, and settings that cannot be changed hidden by a toggle that is on by default | ADM.101, ADM.102, ADM.103, ADM.104, ADM.105, ADM.106 | Parent: closes with its children S1 to S6. Admin control suite, Boss 6:01 PM EDT 2026-10-21; his 13 cards were answered 6:05 to 6:09 PM EDT. |

## Phase 2: Infinite zoom

The fly-through Galaxy Map (MAP.146) and the nebula data its tiles carry.

### Fly-through map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.75 | The mini map as a second engine view |  | Locked second camera on MAP.64. |
| MAP.139 | The Galaxy View uses its spare space: an info box with a Details link, and menu items |  | GitHub issues [#758](https://github.com/dwhagar/planetGen/issues/758) and [#715](https://github.com/dwhagar/planetGen/issues/715), one layout change. |
| MAP.140 | Double-click on a selected object goes there and opens its information |  | GitHub issue [#714](https://github.com/dwhagar/planetGen/issues/714). |
| MAP.141 | Context around the selection: faint neighbours, and the sectors above and below |  | GitHub issues [#718](https://github.com/dwhagar/planetGen/issues/718) and [#716](https://github.com/dwhagar/planetGen/issues/716), one context view. The #716 bug label was overruled. |
| MAP.249 | The Galaxy Map is laggy and slow, and changing zoom level slogs the whole page (bug) |  | Boss 6:48 PM EDT; Core 2 profiles first. |
| MAP.251 | When a map has to drop stars or objects, drop them evenly across the viewing space and keep the most massive and brightest of the galaxy last (bug) |  |  |
| MAP.121 | Every map shows and steps to its neighbouring regions, on one map engine |  | Folds MAP.127, the Galaxy Map's neighbouring blocks and slabs (Boss 2026-10-07 16:27Z: one map engine for all maps). |
| MAP.146 | Fly through the galaxy: scroll-zoom, double-click flight, distance-based visibility and a see-through near field |  | Boss 2026-10-09 22:24Z and 22:56Z. Umbrella for MAP.148 to MAP.152; replaces the fixed ladder in galaxy/drill.py and galaxyprisms.js. |
| MAP.122 | A Select mode on every galaxy view: Galaxy (blocks and sectors) or Star |  | Merges the 2026-10-07 "button to select a star" ask; MAP.101 (PR #431) made stars unpickable by default. |
| MAP.171 | The Galaxy Map keys a massive phenomenon's visibility to its map luminosity (replaces "always lit") |  | Boss 2026-10-10 21:38Z. Foundations lane 2 after MAP.152. |
| MAP.195 | Galaxy Map finer tiles leave out the stars a coarser tile already sent |  |  |
| MAP.173 | The Galaxy Map tile carries a nebula layer: regions and singles, nested like the star lists | MAP.172 | Foundations lane 2. |
| TEST.135 | Big-galaxy query-budget test for the nebula reads | MAP.172, MAP.173, GEN.202 | Foundations lane 1. |
| MAP.156 | View one layer or a range of layers top-down from the galaxy view, as a secondary option |  | Boss 2026-10-09 23:17Z; secondary option, not a main control. |
| MAP.160 | Quantise the Galaxy Map GPU buffers (deferred) |  | Deferred; MAP.147 recommendation step 3. |
| MAP.190 | Several separate regions of created sectors: which one the Galaxy Map opens on |  |  |
| MAP.189 | The Galaxy Map opens flown in to the created sectors | MAP.190 |  |
| MAP.214 | Galaxy Map default view: no density, age or luminosity shading, only stars and markers for known space |  | Split from MAP.191. |
| MAP.209 | Opening a sector from the Galaxy Map flashes one view, reloads into another, and leaves the details panel empty (bug) |  |  |
| MAP.210 | Nothing on the map can be selected, even in sector view (bug) |  |  |
| MAP.237 | The Galaxy Map's "Charted only" view is invisible again, and charted regions are shown by colour and opacity, not outlined blocks (bug) |  | Boss 7:24 PM EDT. |
| MAP.229 | Map label engine: names placed by priority, with no overlaps, shown and hidden by zoom |  | Map labels, Boss 2:30 PM EDT. |
| MAP.234 | Galaxy Map key: a Key section in the help dialog and an on-map key that names everything the map draws |  | Boss 7:16 PM EDT; the superbubble layer is its sibling item. |
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
| GEN.410 | One table of the gases and compounds that can be in an atmosphere, with every compound the descriptions mention |  | Boss 11:24 PM EDT; the one table the rest reads. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.496 | The photolysis lifetime table: species lifetimes at 1 AU for the quiet Sun, the active-Sun factors and the species-specific stellar-spectrum ratios | GEN.410 | Body lane G-2: photolysis lifetimes as constants of the species table. Core 4 proposed. |
| GEN.427 | One table of the bulk materials with their equation-of-state constants and the source of each |  | Boss 11:49 PM EDT; the one materials table. Core 4 proposed. |
| GEN.470 | One saturation-pressure function p_sat(species, T) for the gases and ices: the Fray and Schmitt polynomials, water by Feistel and Wagner, and its tests | GEN.410, GEN.427 | Body composition research lane part F items F1 and F1t; no behaviour change until GEN.452 reads it. Core 4 proposed. |
| GEN.498 | One oxygen-fugacity scale for the generator: the reference fixed, delta-IW stored on it, and the buffers converted |  | Body lane G-4 (part G3): one IW reference scale (Hirschmann 2021), delta-IW stored on it. Core 4 proposed. |
| GEN.497 | The cold-end vapour pressures: the laboratory-corrected N2, CO, CH4, O2 and NH3 below 50 K and the Span and Wagner carbon dioxide sublimation curve | GEN.470 | Body lane G-3 (part G4): the corrected cold-end tables and Span and Wagner CO2 extend GEN.470. Core 4 proposed. |


## Phase 3: New scatters: stars, nebulae and the region query

The multi-star catalogue, the stars every nebula needs and the region query the whole sphere reads through.

### Multi-star and compact systems

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.311 | Built systems store and load their orbit tree: insert_star_system, load_star_system and serialization read and write orbit_nodes |  | Split from GEN.129. |
| GEN.312 | Naming for trees: path codes (Aa, Ab, Ba) with branch words at wide splits, and planet stems such as Pikkita Ba I |  | Split from GEN.129. |
| MAP.213 | Scene, maps and the JavaScript twin for trees: per-node zoom, log-radial overview and positions_at in JavaScript | GEN.311 | Split from GEN.129. |
| GEN.313 | Orbital update of built systems: advance_orbital_phases for nodes and the GEN.109 hooks | GEN.311 | Split from GEN.129. |
| GEN.129 | Multi-star systems of up to twelve stars | GEN.311, GEN.312, MAP.213, GEN.313 | GitHub issue [#777](https://github.com/dwhagar/planetGen/issues/777). Build after the GEN.128 design. |
| TEST.241 | Test: no absolute position is stored for anything inside a system, so every 8.1 item that stores a position stores it parent-relative |  |  |
| GEN.593 | Pulsar planets: the disc range stops short of PSR B1257+12's planets, and two of the three formation paths (a captured giant, an ablated companion) leave a pulsar with no planets |  |  |
| GEN.263 | Promote every bound compact pair to the scatter, whatever its mass |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 4. Core. |
| GEN.315 | A fill builds a system from its stored members and nodes; the system pass keeps only planets, moons, belts and disks | GEN.311 | Split from GEN.264. |
| GEN.316 | Stages D and G draw a member tree for every object the fill creates, and the two-star classes are retired | GEN.315 | Split from GEN.264. |
| GEN.264 | The fill builds a system from its stored members; the system pass keeps only planets, moons and belts | GEN.315, GEN.316 | Multi-star systems at scatter time, Boss 2026-10-15. Step 5. Core. |
| GEN.267 | A MAX_STARS setting and a --stars N directive, with a prevalence knob for five or more stars |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 11. Core. |
| GEN.269 | Interaction states in the star model: a Roche classifier, `interaction` on every pair, member state and flags, rejuvenated age and the new kinds |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 13. Core. |
| GEN.270 | Close-binary types: contact, Algol, spectroscopic and eclipsing flags, tertiary shares, merger products and Be stars | GEN.269 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 14. Core. |
| GEN.253 | White-dwarf pollution flag and cool white dwarf presets |  | Phenomena catalogue, Boss 2026-10-15. Step 6, later. Core. |
| GEN.271 | White-dwarf binary types: close double white dwarfs, post-common-envelope pairs, cataclysmic variables, AM CVn, symbiotic and barium stars | GEN.269 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 15. Core. |
| GEN.272 | Accreting compact binaries: X-ray binaries, accretion luminosity, microquasars, spider pulsars and compact triples | GEN.269 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 16. Core. |
| GEN.329 | Runaway bookkeeping: runaway_origin and the compact companion left by a supernova |  | Split from GEN.274. |
| GEN.327 | The group node and the remnant draw: offset and velocity columns, a decay age, and the lightest member ejected | GEN.329 | Split from GEN.273. |
| GEN.331 | Substellar companions: brown-dwarf close (the desert) and wide companions, brown-dwarf pairs and planetary-mass companions |  | Split from GEN.275. |
| GEN.332 | Circumbinary and S-type planets, circumbinary disks and rogue planet pairs at fill | GEN.315 | Split from GEN.275. |
| GEN.275 | Substellar and circumbinary types: brown-dwarf and planetary-mass companions, circumbinary disks and planets | GEN.331, GEN.332 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 19. Core. |

### Nebula and young stars

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.258 | Protostars embedded in giant molecular clouds |  | Phenomena catalogue, Boss 2026-10-15. Step 8, later. Core. |
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
| GEN.369 | Satellites drawn from the galaxy seed: count, distance, mass, size, profile, [Fe/H] and age are calculated from stated laws, with no Milky Way option |  | Core 1; Boss 4:01 PM EDT: Draw them. |
| GEN.535 | Check the Magellanic gas fraction and supernova-remnant count against the observed Clouds, and the dwarfs that still come out with gas (Leo I) |  | Open inferences of GEN.528: gas fraction, low LMC remnant count, Leo I gas. |
| GEN.534 | A satellite's young stars reach the system-type tables: the young-star system types name `satellite_young` |  | Residual of GEN.528: the system-type tables do not name satellite_young. Core 3 proposed. |
| GEN.530 | A pre-placed star's hosted nebula shows from its scatter row on the maps and pages |  | Residual of GEN.320: the hosted nebula of a pre-placed star from its scatter row. |
| GEN.328 | Groups at fill in star-forming cores and OB regions | GEN.327 | Split from GEN.273. |
| GEN.273 | Young groups and trapezia: the `group` node, offset and velocity columns, a decay age and the remnant draw | GEN.327, GEN.328 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 17. Core. |
| GEN.365 | Hot white dwarf and planetary-nebula central star radius disagrees with the radius that its luminosity and temperature imply (bug) |  | Finding from Core 3. |

### Globular clusters

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.325 | Planet cull in dense components by the survival formula and the metallicity factor |  | Split from GEN.162. |

### Whole-sphere audits

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.184 | Whole-sphere tests, nebula and Wolf-Rayet side: contents per class, the Wolf-Rayet count, tiling, determinism and the King solver |  | Split from TEST.159. |
| TEST.161 | Fill-stage and system-pass tests: stream independence, same seed same system, and the probability table |  | Phenomena catalogue, Boss 2026-10-15. Step 2 tests. Quality. |
| TEST.162 | Multiplicity audit: companion states, row bookkeeping, stability, determinism, no planets in a scatter, fill reuse, motion and the census anchors | GEN.315 | Multi-star systems at scatter time, Boss 2026-10-15. Step 9. Quality. |
| TEST.164 | Catalogue type audit: type counts, region gates, the N tail, interaction consistency and group decay | GEN.269, GEN.270, GEN.271, GEN.272, GEN.275, GEN.327, GEN.329 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 20. Quality. |
| TEST.165 | Re-measure the region query on a real generated galaxy, warm and cold |  | Spherical region query, stage C. Step 6. Quality. |

### Composition proofs

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.211 | A star composition sweep: every kind of star, remnant and brown dwarf, fractions, trends and determinism |  | Boss 11:58 PM EDT; the proof of the star composition. |
| TEST.215 | The Solar System calibration rows live in test files only: a flag fences them and a test fails when a generator table reads one | TEST.211 | Boss 1:44 AM EDT, Stand-in only: R5, calibration rows fenced into tests. |
| TEST.214 | Every registered possibility table: for every class the gated weights add up to 100 percent, no class has an empty gate, and every row has a source and a licence | GEN.410, GEN.427 | Boss 1:40 AM EDT: the gated weights add up to 100 percent; registry walk, so later tables are covered. Core 4 proposed. |

## Phase 5: New scatters: core, halo, interface update and docs

The realistic galactic core and the halo tier, the interface update for the spherical grid, the headline item and the docs that follow.

### Galactic core and halo

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.589 | Research: the central molecular zone's upper cluster mass and its cloud-shock disruption, so the young clusters of the core match the few that are observed |  | Finding from Core 3; research. |
| GEN.377 | Scatter limits by place for the neutron star and black hole scatter below the cut: the core and halo tiers apply to it, not one global mass cut |  | Core 3 finding after GEN.333. |
| GEN.303 | Per-cell system cap as a configuration option: `generation.cell_cap`, default 5,000, lowest 500, highest 50,000 |  | Core 3 wires it; Boss 1:47 PM EDT. |
| GEN.284 | Planet survival at the galactic core: a metallicity factor and a survival formula instead of the hard cut | GEN.325 | Core 3; step 4 of the galactic core note |

### Interface update

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.115 | Interface update: every page and control works with the spherical grid, 96-bit IDs, the new phenomena and the region queries (stage after the sphere and generation stages) | NAV.61 | Interface; queued after sphere stages (Boss, 12:52 AM EDT) |
| TEST.168 | Control sweep: a browser test drives every control on every page against a spherical galaxy with 96-bit IDs | UX.115 | Quality; queued after UX.115 |

### Docs

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.19 | After UX.102: document the log output standard for the CLI, the Generate page and admin logs | UX.102 |  |
| DOC.72 | The README lists every external file the program needs, with its URL and how to install it, and a check keeps the list complete |  |  |
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
| GEN.162 | Planet cull and blue stragglers in clusters | GEN.325, GEN.326 | Globular-cluster chain. |
| GEN.163 | Type-B pulsar planets in globular clusters (GEN.130 follow-on) |  | Globular-cluster chain. |
| GEN.330 | Cluster-hosted compact binaries: the cluster rates replace the field fractions | GEN.329, GEN.163, GEN.325 | Split from GEN.274. |
| GEN.274 | Cluster-hosted compact binaries and runaway bookkeeping: `runaway_origin` and the compact companion left by a supernova | GEN.329, GEN.330 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 18. Core. |
| GEN.297 | Calibrate the cluster generator against the Harris file: its distributions guide the laws, and no real cluster is placed in any galaxy |  | The file arrived (Boss, 3:04 PM EDT); Core 3, cluster line. |
| GEN.366 | Check the reference satellite table against McConnachie 2012: it only tests that the derived population looks plausible | GEN.369 | Finding from Core 3. |

### Galaxy types groundwork

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.378 | The galaxy type families: the checked design with its scaling relation table and the draw order |  | Boss 8:48 PM EDT, for 9.4: step 1 of the galaxy types chain; Core 1 proposed. |

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
| GEN.403 | The globular cluster count follows the halo mass: the system mass from the halo mass replaces the fixed count of 157 |  | GEN.368 build, Boss approved the table 10:11 PM EDT; Core 1 proposed for GEN.400 to GEN.403, Core 3 or Core 4 for GEN.404 and GEN.405 (profile.py, activity.py), any lane for GEN.406. |
| GEN.405 | Metallicity and rates by position from the galaxy's own model: the mass-metallicity relation and the local rates replace the fixed gradient and the quoted-at-the-Sun values |  | GEN.368 build, Boss approved the table 10:11 PM EDT; Core 1 proposed for GEN.400 to GEN.403, Core 3 or Core 4 for GEN.404 and GEN.405 (profile.py, activity.py), any lane for GEN.406. Core 4 proposed (coordinator 10:16 PM EDT). |
| GEN.545 | The mean system mass (0.4 solar, recalled within a factor of 1.5) sets every galaxy's system count once GEN.404 lands: derive it from the generator's own mass distribution |  |  |
| GEN.607 | Shape dials: the Census dial: star richness (systems per solar mass, equivalently the mean system mass), advanced and hidden by default | GEN.601, GEN.545 | Shape dials (Boss 6:45 PM EDT 2026-10-21). |
| GEN.598 | CLUSTER_MEAN_SYSTEM_MASS_SOL (0.4, recalled) sets the system count of globular clusters and satellites: derive it as GEN.545 derives the galaxy's | GEN.545 |  |
| GEN.597 | A rocky planet of a wide binary gets flux_variation above 1 (the insolation test passes only on seeds 0 to 59) (bug) |  |  |
| GEN.596 | A companion below 0.01 solar masses makes multistar.draw_member_tree raise 'no cooling track' (bug) |  |  |
| GEN.531 | Boss's research requests from the galaxy laws research (R1 to R15): the observed stellar halo, the dwarf range, the thick disk fits and the cluster system |  | Research requests R1 to R15 of the GEN.406 note; does not gate GEN.301. |

## Phase 6: The headline closes

The last audits and the headline: they test the component pass, the core and the halo, so they come after phases 4 and 5.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every point release.

### The headline closes

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.185 | Whole-sphere tests, component side: pre-placed share per population, no double count of a cluster star, nothing hosted in a halo or old cluster cell, D and E members inside their radius |  | Split from TEST.159. |
| TEST.159 | Whole-sphere tests: required stars per nebula kind, Wolf-Rayet counts, tiling at the reach, profile and empty-sector determinism, and the pre-placed share per population | TEST.184, TEST.185 | Whole-sphere build train, Boss 2026-10-15. Step 15. Quality. |
| TEST.160 | Representation audit: seeded sectors in each region class show every kind at its expected frequency, with systems and required contents checked |  | Phenomena catalogue, Boss 2026-10-15. Step 4. Quality. |
| TEST.196 | Close-binary check: generated close-binary fractions by primary mass and [Fe/H] bin against the published ranges, with out-of-range bins flagged |  | Core 4; Boss 6:15 PM EDT: expected-numbers check for GEN.103. |
| PERF.146 | A performance profile of every generation option and scatter, with a repeatable benchmark command, a report and a ranked optimization list: one of the last three 8.1 stages | TEST.159, TEST.160 |  |
| TEST.261 | Release gate: a full bug hunt audit of the whole code base and its behaviour, done before every point release (8.1 first) | TEST.159, TEST.160, PERF.146 | Standing release gate, Boss 2026-10-21 5:10 PM EDT; also gates every later point release. |
| SEC.33 | Release gate: a full security code audit of the whole code base, done before every point release (8.1 first) | TEST.159, TEST.160, PERF.146 | Standing release gate, Boss 2026-10-21 5:10 PM EDT; also gates every later point release. |
| GEN.301 | Headline of 8.1: the new generate scatters, accurate galactic generation and scattering including the halo | GEN.284, TEST.159, TEST.160, TEST.196, GEN.450, GEN.530, SEC.33, TEST.261 |  |
