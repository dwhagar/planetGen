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
| PERF.20 | Short-term caching through the work queue and API |  | Plans on RQ results and the new caches. Plan with Boss. |
| PERF.36 | Memory and request guard: list at most about 50,000 candidate cells, stream the rest, and refuse enumerations above about 2 million cells in a web request |  | Research: protects the Generate page estimates. |
| PERF.41 | Stamp the tile and page caches with a TILE_FORMAT constant instead of the version (optional) |  | Research follow-up to PERF.25 (built). |
| PERF.75 | Keyset paging for the other data tables and the jump-by-value box (Systems name sort done) |  | Decided by default. Foundations lane 1. |
| PERF.86 | Redo the work queue: drop the task and subtask logic and use Redis simply |  |  |
| PERF.91 | Bulk random draws for the scatter layers (4 to 12 times on batches above 1,000 draws) |  |  |
| PERF.106 | The plan command on the run record, with stage timing, the killed-worker handler and job timeouts |  | Split from PERF.98. |
| PERF.107 | The sector units and the stages around them: start, neighborhood, sectors, link, population and settle | PERF.106 | Split from PERF.98. |
| PERF.108 | The backfill as queue units: chunks of at least 10,000 cells, one enumeration (folds in PERF.18) | PERF.107 | Split from PERF.98. |
| PERF.98 | Work queue step 4: the remaining stages and the sector and backfill units on the new queue | PERF.106, PERF.107, PERF.108 | PERF.86 step 4 |
| GEN.308 | A block-first sampler for the backfill: a Poisson candidate count from the block's density bound, then thinning |  | Split from GEN.42. |
| GEN.309 | The backfill picks the draw path per block and gives the same stars wherever it lands | GEN.308, PERF.108 | Split from GEN.42. |
| GEN.40 | Weed out sectors by star density before the bright-star backfill | GEN.309 | Parent; closes with its subitems. |
| PERF.35 | An interval or chunk ledger for untouched sectors once block-first backfill lands | GEN.309 | Research: ledger size at full galaxy. |
| PERF.99 | Work queue step 5: the Generate page, Queue page and job log read Redis; cancel now through the stop command | PERF.98 | PERF.86 step 5 |
| GEN.310 | Fill a whole block of the shell ladder: planetgen generate --block and the queue job | PERF.99 | Split from GEN.231. |
| UX.127 | The Galaxy Map block menu offers Generate again | GEN.310 | Split from GEN.231. |
| GEN.231 | Generate a whole block on the shell ladder: the generate --block mode and the block buttons come back | GEN.310, UX.127 | Maps hand-off 2026-10-15. Core, then Interface for the buttons. |
| DB.29 | update.sh and every database migration refuse to start while jobs are active, unless --yes is passed or the user confirms cancelling them | PERF.99 |  |
| PERF.100 | Work queue step 6: delete the old machinery and the work_* tables (Alembic migration) | DB.29 | PERF.86 step 6 |
| PERF.105 | Slow page requests queue a job that builds the answer into Redis, with a status route (the mechanism) |  | Boss 2026-10-16, 12:21 PM EDT; Core 2, after MAP.152 and MAP.146. |
| UX.118 | Slow pages show a progress bar while their data loads: the search page and the phenomena list pages | PERF.105 | Boss 2026-10-16; Core 2, after PERF.105. |
| OPS.59 | `generate galaxy --shell 0` does not stop on SIGTERM: the busy `worker --burst` child keeps running and the parent is left a zombie (bug) |  | Finding from a lane. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |

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
| GEN.218 | One habitability describer, in the short wording |  |  |
| GEN.219 | The habitability levels page and tooltips are generated from the limits |  |  |
| DB.50 | Migration only: structure columns for life and flavour |  | Split from DB.31. Core 1 owns migrations. |
| DB.51 | Generation stores structure, not prose; the readers and the string searches go | DB.50, GEN.218 | Split from DB.31. Core 1 owns migrations. |
| UX.119 | An HTML renderer for the Document: facts, rating and composition tables, with unit spans, in both themes |  | Split from UX.105. |
| UX.120 | The system page and the one-off system page use it: every body's detail from the Document | UX.119 | Split from UX.105. |
| UX.105 | The system page builds every body's detail from the document model as HTML tables | UX.119, UX.120 |  |
| GEN.220 | "Habitable" means the index | GEN.218, UX.120 |  |
| UX.108 | Describe phenomena, sectors, galaxy and facilities in the same Document |  |  |
| GEN.221 | Short plain prose for the fixed text |  |  |
| UX.109 | Remove the old per-class description text | UX.120, UX.108 |  |
| DB.52 | Migration only: drop the prose columns and tables | DB.51 | Split from DB.31. Core 1 owns migrations. |
| DB.31 | Store structure, not prose, for life and flavour | DB.50, DB.51, DB.52 |  |
| TEST.155 | Wording and index agree: property test | GEN.218 |  |
| UX.107 | Units in the description go through the ladders |  |  |
| DOC.45 | Reference docs follow the description layer | UX.105, GEN.218, GEN.219 |  |
| UX.110 | The phenomena count on the overview page is wrong (bug) |  | Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |

### Admin and API housekeeping

| ID | Item | Needs | Note |
|---|---|---|---|
| API.26 | Bump API_VERSION for the shell, band and wedge sector address |  |  |
| API.4 | API compatibility data in the docs | API.26 | Docs and version number. |
| API.15 | Log every API call with its user, how it came in, and its HTTP response code |  | Needs no user accounts (Boss 01:31Z). |
| API.31 | Every outgoing web request sends a standard User-Agent with the software, its version and the API version |  | Needs nothing. Boss, 2026-10-17, 11:33 AM EDT. |
| OPS.35 | A vendored-version lock file for the static libraries |  | Research: supply-chain record. |
| UX.124 | One form-field look and error slot, piloted on the search, nearby and list-filter fields |  | Split from UX.49. |
| UX.125 | The Generate page and the one-off system page on the field macro | UX.124 | Split from UX.49. |
| UX.126 | The admin pages and the edit dialogs on the field macro | UX.124 | Split from UX.49. |
| UX.49 | Form fields as Shoelace components (sl-input, sl-select, sl-checkbox) across the site | UX.124, UX.125, UX.126 | Follows UX.40 (done): the fields themselves become Shoelace components. |
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
| UX.3 | Warn every visitor while a background job changes the galaxy |  | ETA from the RQ job's published progress. ETA from progress.json, which PERF.23 caps. |

### Bugs and flaky tests

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.33 | MySQL 8.0.46 reflects the foreign keys fk_facilities_type and fk_facility_field_values_field as RESTRICT, so check-db calls a healthy database damaged (bug) |  | Quality run 2026-10-15, MySQL 8.0.46. Core. |
| TEST.111 | test_ensure_sector_generated_creates_then_reuses_the_same_sector fails in a busy parallel run (bug) |  | Test flake reported 2026-10-09; Bugfixes lane 1. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.116 | test_bughunt_end_to_end stores M2V or M6V where it expects K2V on the MySQL 8.4 and MariaDB legs once (bug) |  | Bugfixes lane 1 report 03:57Z. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.126 | The heavy WebGL browser-map tests time out when four workers run them together (bug) |  | Lane 2 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.127 | test_nebula_shape_endpoint_serves_a_mesh failed once in a full parallel run (bug) |  | Lane 2 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.128 | test_spatial_position_db and test_web_db_fields fail under a parallel run on one MySQL (bug) |  | Lane 1 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.129 | Two orbit-update tests fail on the MySQL 8.4 leg: the stored last_updated_at rounds up (bug) |  | Lane 1 report 2026-10-10 20:55Z; Bugfixes lane 1 after its current items. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.130 | test_the_scene_positions_at_the_epoch_match_the_stored_ones fails when a random system has comets (bug) |  | Lane 1 report 2026-10-10 20:55Z; Bugfixes lane 1 after its current items. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.134 | test_many_processes_starting_at_once_start_one_job fails under a parallel run (bug) |  | Flaky under load. Bugfixes lane 1. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.137 | test_a_neighborhood_from_the_generate_page[core] fails once under load (bug) |  | Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.140 | test_a_body_is_reached_by_its_id_and_a_row_number_names_nothing fails once in a full parallel run (bug) |  | Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.141 | test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector[15-True] fails once in a full parallel run (bug) |  | Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.170 | test_math_check::test_whole_run_is_fast and the Galaxy Map breadcrumb folding test fail under full-suite load and pass alone (bug) |  | Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.171 | test_gen_resume::test_sectors_a_forced_scatter_skipped_fill_correctly_afterwards failed once under gate load with no scattered system (0 > 0) and passes alone (bug) |  | Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.172 | test_galaxy_map_forward_to_current_jumps_to_the_newest_view failed once: clicking a Galaxy Map block with 0 sectors charted did nothing, and it passes alone (bug) |  | Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.174 | test_web_a11y::test_page_is_accessible_and_fits[web.sector-desktop-light] failed once under gate load and passes alone (bug) |  | Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| UX.133 | `GET /polities` and `GET /species` answer 404 before the population pass has run: show an empty page instead (bug) |  | Finding from a lane. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.189 | A test that runs `migrate` leaves maintenance mode on: restore `state/maintenance.json` in the fixture (bug) |  | Test-infrastructure finding. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.190 | A failing session fixture makes every dependent test wait its full 900 s timeout: fail the dependents at once (bug) |  | Test-infrastructure finding. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.191 | The plan tests write real settings JSON into `/var/lib/planetGen/jobs/galaxy-settings` and log "No naming key drawn ... Unknown database": give them a temporary directory and a test database |  | Test-infrastructure finding. |
| TEST.194 | `test_system_map_selection_drill_and_measure` failed once in the browser phase of a full run under load and passed alone: find the cause and make it robust (bug) |  | Flaky test finding from the Settings lane, 5:47 PM EDT. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.213 | `test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star` failed in a loaded browser run: the wheel loop opened the system before the disc reached full share (bug) |  | Core 1 finding, 2026-10-17; any lane, Bug Hunt proposed. |
| TEST.220 | test_the_scatter_lists_each_shells_phenomena_by_kind_like_the_star_scatter finds `neutron-star` inside `magnetar neutron-star` (bug) |  | Core 2 finding, 2026-10-17; Bug Hunt proposed. |
| TEST.221 | test_a_runaway_system_moves_beyond_the_rotation_curve_and_its_stars_and_bodies_share_it raises StopIteration when the first system's only planet entry is an asteroid belt (bug) |  | Core 5 finding, 2026-10-17; Core 5 carries the fix. |
| TEST.216 | `test_systems_inside_a_nebula_point_at_it_and_the_innermost_wins` failed in the full parallel gate: the nested system's `inside_nebula_id` was None, and it passes alone (bug) |  | Core 4 finding (UX.114 gate, 2:44 AM EDT as the lane reported): fails only in the full parallel run, passes alone on main. |
| PERF.111 | `ResponseCache.put` without a database never trims the page cache: the cache grows past `max_entries`, and the thread test fails when its last put is one of those (bug) |  | Core 2 finding in the MAP.217 gate: the db=None branch of ResponseCache.put has no eviction. |
| TEST.219 | Browser pages still have 1 to 7 requests running 8 s after the server stops: find what keeps them busy and make it faster |  | Core 1 finding at #1997: browser pages still have requests running 8 s after the server stops. |
| TEST.218 | `test_first_control_connections_at_once` leaves 3 idle connections open on `planetgen_test_*` when its fixture ends: the test should close its own connections (bug) |  | Core 1 finding at #1997: the test leaves 3 idle connections for the drop warning to end. |
| TEST.198 | `test_galaxy_map_does_not_prefetch_on_a_data_saving_connection` failed with a `/galaxy/tiles` net::ERR under load and passed alone: find the cause and make it robust (bug) |  | Flaky test finding from Core 3, 6:22 PM EDT. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.203 | `test_one_worker_leaves_the_runs_own_random_stream_alone` failed once under heavy load with WorkerDied "the task's job ended canceled" and passes alone (bug) |  | Settings lane finding. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.199 | `test_skeleton_bound_is_an_upper_bound_at_every_angle` fails on main for a far-tail case: the halo subtraction loses precision, 2.4e-9 against a 1e-9 tolerance (bug) |  | Core 2 found it; fails on main; also a GEN.372 finding. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| GEN.397 | `planetgen plan` stage 2 prints the progress weight as the number of phenomena to place (bug) |  | Rescatter thread finding, 2026-10-16; any Core lane. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| GEN.398 | The massive-star and bright-star passes place more stars than their log lines expect: 70 and 12 percent over (bug) |  | Rescatter thread finding; Core 3 proposed (scatter spine). |
| DOC.62 | Docs: the generation cost figures and the redo scope are stale (the phenomena cut, rows, bytes and time; `--redo-scatters` leaves satellites and clusters alone) |  | Rescatter thread finding; docs only, any lane. |
| GEN.463 | The neutron star's radius is derived from sourced mass-radius relations instead of the fixed range of 10 to 13 km (bug) |  | Body composition research lane finding F5; Boss 1:22 AM EDT card Replace both; generator epoch bump. Core 3 proposed. |
| GEN.467 | The main-sequence lifetime is derived from sourced stellar relations instead of the fixed 10 Gyr for one solar mass and a power law of -2.5 (bug) |  | Body composition research lane finding F5; Boss 1:22 AM EDT card Replace both; generator epoch bump. Core 3 proposed. |
| GEN.396 | Revisit the Yerkes class III luminosity range: 50 to 1000 solar luminosities sits above many dim red giants |  | Core 3 finding, 2026-10-16 10:00 PM EDT; any Core lane. |
| GEN.502 | The giant branch finishes the stellar tables check: the red giant time law as a function of core mass, the tip by metallicity and the branch duration against the observed horizontal-branch ratio |  | Residual of GEN.407: law as a function, tip by metallicity, branch duration. |
| TEST.208 | `scripts/test.py --resume` reuses an earlier pass across builds when none of that test's own dependencies changed |  | Boss, 10:37 PM EDT; saves hours per cancelled gate; Bug Hunt proposed. |

### Nebula data foundations

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.202 | Scattered nebulae (planetary nebulae and remnants) join the nebula index, deduplicated by object id |  | Foundations lane 1. |
| MAP.172 | Nebula cell aggregates in the region pyramid |  | Foundations lane 1. |
| MAP.174 | Nebula cover in the stage-view cell statistics, for filled and unfilled cells | MAP.172 | Foundations lane 2. |

### Composition research and design

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.409 | The atmosphere colour design: the checked plan for the gas table, the mixture, the stored values, ozone, aerosols and the colour model |  | Boss 11:24 PM EDT; design note, no code. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.426 | The bulk composition design: the fractions, the equation of state, the draws and the derived density, mass and gravity |  | Boss 11:49 PM EDT; design note, no code. Core 4 proposed. |
| GEN.465 | A second extraction pass over the sources the body and atmosphere research could not reach and that are reachable now |  | Body composition research lane master section 8; read-only. Unowned. |
| GEN.468 | Settle the water share beyond the snow line: the Marboeuf-derived 0.14 to 0.39 of the solids against the convention of about 0.5 |  | Body composition research lane master section 5 decision 2; Boss 1:21 AM EDT card Keep the switch; read-only. Unowned. |
| GEN.457 | The small bodies and moons composition design: the taxon and comet draws, compaction, the moon origin types and the ice fraction |  | Body composition research lane part D; design note, no code. Core 4 proposed. |
| GEN.435 | The stellar composition design: the abundance relations, what is stored and what is derived, and the remnants |  | Boss 11:58 PM EDT; design note, no code. |
| GEN.436 | The radiation design: atmospheric transmission and greenhouse by band, shielding by column, and the ground's own dose |  | Boss 11:58 PM EDT; design note, no code. |
| GEN.437 | The composition model: how the element, compound, material and body tables link, and what is stored against derived | GEN.409, GEN.426, GEN.435, GEN.436 | Boss 11:58 PM EDT; design note, no code. |


## Phase 2: Infinite zoom

The fly-through Galaxy Map (MAP.146) and the nebula data its tiles carry.

### Fly-through map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.75 | The mini map as a second engine view |  | Locked second camera on MAP.64. |
| MAP.139 | The Galaxy View uses its spare space: an info box with a Details link, and menu items |  | GitHub issues [#758](https://github.com/dwhagar/planetGen/issues/758) and [#715](https://github.com/dwhagar/planetGen/issues/715), one layout change. |
| MAP.140 | Double-click on a selected object goes there and opens its information |  | GitHub issue [#714](https://github.com/dwhagar/planetGen/issues/714). |
| MAP.141 | Context around the selection: faint neighbours, and the sectors above and below |  | GitHub issues [#718](https://github.com/dwhagar/planetGen/issues/718) and [#716](https://github.com/dwhagar/planetGen/issues/716), one context view. The #716 bug label was overruled. |
| MAP.245 | The wheel zoom's remaining parts: the column-density ray fallback, the clearance from the nearest pickable thing and the half-way double click |  | Split from MAP.225 (wheel part merged, #2031). |
| MAP.121 | Every map shows and steps to its neighbouring regions, on one map engine |  | Folds MAP.127, the Galaxy Map's neighbouring blocks and slabs (Boss 2026-10-07 16:27Z: one map engine for all maps). |
| MAP.146 | Fly through the galaxy: scroll-zoom, double-click flight, distance-based visibility and a see-through near field | MAP.245 | Boss 2026-10-09 22:24Z and 22:56Z. Umbrella for MAP.148 to MAP.152; replaces the fixed ladder in galaxy/drill.py and galaxyprisms.js. |
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
| MAP.239 | Raise the Galaxy Map tile root to the 262 kpc reach: `TILE_ROOT` in `galaxy/viewport.py` and `galaxymap3d.js` |  | GEN.372 audit; Boss card 7:36-7:38 PM EDT. |
| MAP.230 | Names of stars, systems and sectors on the Galaxy Map | MAP.229 | Map labels, Boss 2:30 PM EDT. |
| MAP.231 | The Sector Map and the System Map use the label engine | MAP.229 | Map labels, Boss 2:30 PM EDT. |


### Golden galaxy

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.182 | Put one block-first backfill into the golden galaxy | GEN.309, TEST.77 | Split from GEN.42. |
| GEN.42 | A block-first draw for the backfill (umbrella) | GEN.308, GEN.309, TEST.182 | Same function as PERF.18 (backfill_bright_stars_around); one thread. |

### Star-system storage

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.53 | Built-system star storage: the orbit_nodes table, stars.parent_node_id, planets.host_node_id, the catalogue columns and the disk descriptor (one Core 1 revision) |  | Split from GEN.129. Core 1 owns migrations. |

### Composition tables

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.472 | One gated possibility-table format and loader: every composition category is a table of possibilities with a class gate and a weight, drawn by the normalised weights of the rows that pass | GEN.437 | Boss 1:40 AM EDT: gated possibility tables with weighted probabilities; the one format and loader every composition table uses. Core 4 proposed. |
| GEN.500 | The licence audit: the licence column filled for every data source, and the permission list for NIST, JPL, Leiden, Bremen and Karkoschka |  | Body lane G-6 (part G4): the licence audit; NIST SRD is copyrighted, JPL 19-5 is not public domain. Core 4 proposed. |
| GEN.438 | One table of the elements and their isotopes: atomic data and the radioactive isotopes' half-lives and decay energies | GEN.437, GEN.500 | Boss 11:58 PM EDT; atomic data and the radioactive isotopes; no abundance pattern. |
| GEN.410 | One table of the gases and compounds that can be in an atmosphere, with every compound the descriptions mention | GEN.409, GEN.437, GEN.438, GEN.472, GEN.500 | Boss 11:24 PM EDT; the one table the rest reads. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.496 | The photolysis lifetime table: species lifetimes at 1 AU for the quiet Sun, the active-Sun factors and the species-specific stellar-spectrum ratios | GEN.410 | Body lane G-2: photolysis lifetimes as constants of the species table. Core 4 proposed. |
| GEN.427 | One table of the bulk materials with their equation-of-state constants and the source of each | GEN.426, GEN.437, GEN.438, GEN.472, GEN.500 | Boss 11:49 PM EDT; the one materials table. Core 4 proposed. |
| GEN.470 | One saturation-pressure function p_sat(species, T) for the gases and ices: the Fray and Schmitt polynomials, water by Feistel and Wagner, and its tests | GEN.410, GEN.427 | Body composition research lane part F items F1 and F1t; no behaviour change until GEN.452 reads it. Core 4 proposed. |
| GEN.498 | One oxygen-fugacity scale for the generator: the reference fixed, delta-IW stored on it, and the buffers converted | GEN.437 | Body lane G-4 (part G3): one IW reference scale (Hirschmann 2021), delta-IW stored on it. Core 4 proposed. |
| GEN.497 | The cold-end vapour pressures: the laboratory-corrected N2, CO, CH4, O2 and NH3 below 50 K and the Span and Wagner carbon dioxide sublimation curve | GEN.470 | Body lane G-3 (part G4): the corrected cold-end tables and Span and Wagner CO2 extend GEN.470. Core 4 proposed. |
| GEN.439 | The stellar composition relations: abundances by population, age and [Fe/H], and the remnants' layers | GEN.435, GEN.437, GEN.438, GEN.472 | Boss 11:58 PM EDT; a constants-and-functions module, no database access. |


## Phase 3: New scatters: stars, nebulae and the region query

The multi-star catalogue, the stars every nebula needs and the region query the whole sphere reads through.

### Multi-star and compact systems

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.311 | Built systems store and load their orbit tree: insert_star_system, load_star_system and serialization read and write orbit_nodes | DB.53 | Split from GEN.129. |
| GEN.312 | Naming for trees: path codes (Aa, Ab, Ba) with branch words at wide splits, and planet stems such as Pikkita Ba I |  | Split from GEN.129. |
| MAP.213 | Scene, maps and the JavaScript twin for trees: per-node zoom, log-radial overview and positions_at in JavaScript | GEN.311 | Split from GEN.129. |
| GEN.313 | Orbital update of built systems: advance_orbital_phases for nodes and the GEN.109 hooks | GEN.311 | Split from GEN.129. |
| GEN.314 | Planet zones in a hierarchy and flux sums: S-type and P-type windows per node (Holman-Wiegert) and the flux a planet receives from every member |  | Split from GEN.129. |
| GEN.129 | Multi-star systems of up to twelve stars | DB.53, GEN.311, GEN.312, MAP.213, GEN.313, GEN.314 | GitHub issue [#777](https://github.com/dwhagar/planetGen/issues/777). Build after the GEN.128 design. |
| GEN.317 | Pulsar planets and the radiation environment of neutron stars | GEN.314 | Split from GEN.130. |
| GEN.318 | White-dwarf and planetary-nebula-star rows of the system pass: engulfment of inner planets and the planetary nebula's circumbinary disk | DB.53, GEN.314 | Split from GEN.130. |
| GEN.130 | Compact-remnant systems at sector fill: companions and pulsar planets for neutron stars and black holes | GEN.317, GEN.318 | Phenomena catalogue, Boss 2026-10-15. Step 3, rescoped. Core. |
| GEN.252 | Brown dwarfs as stars: an L, T and Y state for companions and primaries, with the brown-dwarf desert |  | Phenomena catalogue, Boss 2026-10-15. Step 5, later. Core. |
| GEN.263 | Promote every bound compact pair to the scatter, whatever its mass |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 4. Core. |
| GEN.315 | A fill builds a system from its stored members and nodes; the system pass keeps only planets, moons, belts and disks | GEN.311, GEN.314 | Split from GEN.264. |
| GEN.316 | Stages D and G draw a member tree for every object the fill creates, and the two-star classes are retired | GEN.315 | Split from GEN.264. |
| GEN.264 | The fill builds a system from its stored members; the system pass keeps only planets, moons and belts | GEN.315, GEN.316 | Multi-star systems at scatter time, Boss 2026-10-15. Step 5. Core. |
| GEN.267 | A MAX_STARS setting and a --stars N directive, with a prevalence knob for five or more stars |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 11. Core. |
| GEN.269 | Interaction states in the star model: a Roche classifier, `interaction` on every pair, member state and flags, rejuvenated age and the new kinds | DB.53 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 13. Core. |
| GEN.270 | Close-binary types: contact, Algol, spectroscopic and eclipsing flags, tertiary shares, merger products and Be stars | GEN.269 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 14. Core. |
| GEN.253 | White-dwarf pollution flag and cool white dwarf presets | GEN.318 | Phenomena catalogue, Boss 2026-10-15. Step 6, later. Core. |
| GEN.271 | White-dwarf binary types: close double white dwarfs, post-common-envelope pairs, cataclysmic variables, AM CVn, symbiotic and barium stars | GEN.269 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 15. Core. |
| GEN.272 | Accreting compact binaries: X-ray binaries, accretion luminosity, microquasars, spider pulsars and compact triples | GEN.269 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 16. Core. |
| GEN.449 | The companion table for mid and early B stars is raised to the observed frequencies |  | Core 4 citation re-check, 2026-10-17; follow-up of GEN.103. Core 4. |
| GEN.448 | Massive stars get binaries between 0.25 and 50 AU, so the metallicity trend reaches the separations where it was measured |  | Core 4 citation re-check, 2026-10-17; follow-up of GEN.103. Core 4. |
| GEN.447 | The wide-pair encounter limit uses the relative-velocity dispersion and the mass-weighted perturber mass |  | Core 4 citation re-check, 2026-10-17; follow-up of GEN.103. Core 4. |
| GEN.329 | Runaway bookkeeping: runaway_origin and the compact companion left by a supernova | DB.53 | Split from GEN.274. |
| GEN.327 | The group node and the remnant draw: offset and velocity columns, a decay age, and the lightest member ejected | DB.53, GEN.329 | Split from GEN.273. |
| GEN.331 | Substellar companions: brown-dwarf close (the desert) and wide companions, brown-dwarf pairs and planetary-mass companions | GEN.252 | Split from GEN.275. |
| GEN.332 | Circumbinary and S-type planets, circumbinary disks and rogue planet pairs at fill | DB.53, GEN.314, GEN.315 | Split from GEN.275. |
| GEN.275 | Substellar and circumbinary types: brown-dwarf and planetary-mass companions, circumbinary disks and planets | GEN.331, GEN.332 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 19. Core. |

### Nebula and young stars

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.529 | The pre-main-sequence track takes its Hayashi temperature, birthline radius and energy factor from a sourced table, and the state reaches stars still contracting after 10 Myr |  | Residual of GEN.242: three inferred constants and the 10 Myr cut; needs a sourced Baraffe 2015 or Siess 2000 table. |
| GEN.243 | Wolf-Rayet stars and ring nebulae: a W state for the last 0.3 Myr of very massive stars, and nebula class X |  | Whole-sphere build train, Boss 2026-10-15. Step 13. Core. |
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
| GEN.440 | Star composition is derived from the stored [Fe/H], age and population, for every star, remnant and brown dwarf | GEN.439, GEN.472 | Boss 11:58 PM EDT; a pure function of stored values; no column, no stored change. |
| GEN.464 | The disk solids law reads the star's metal mass fraction Z over Z0, not 10 to the power of [Fe/H], which is off by 1.7 to 3 times in metal-poor alpha-enhanced stars | GEN.440 | Body composition research lane finding F6 (A-S2); generator epoch bump. Core 4 proposed. |

## Phase 4: New scatters: components, clusters and audits

Satellites, clusters and wanderers as components, the globular cluster chain and the audits that check what the scatters make.

### Components and wanderers

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.450 | The star-formation profile's radial slope: reconcile the Schmidt index with the observed gradient, with the choice named |  | Core 4 citation re-check, 2026-10-17; follow-up of GEN.134. Core 4. |
| GEN.469 | The disk's star-formation history from the gas law: the local B-star share is 0.2 percent against the observed 0.04 percent |  | GEN.134 follow-up (proposal 2); generator epoch bump, joins the one rescatter. Core 4. |
| GEN.371 | Feasibility study: generate the satellite galaxies through the scatter, making only the region each one occupies |  | Core 1; Boss 4:02 PM EDT: feasibility study. |
| GEN.369 | Satellites drawn from the galaxy seed: count, distance, mass, size, profile, [Fe/H] and age are calculated from stated laws, with no Milky Way option | GEN.371 | Core 1; Boss 4:01 PM EDT: Draw them. |
| GEN.528 | The rest of GEN.319: a Magellanic cell shows every class, and the scatter's planetary nebulae are weighted by the old share |  | Residual of GEN.319: the Magellanic cell classes and the planetary nebula scatter weight. |
| GEN.527 | The second star of a wide pair in a saved runaway system loads back 0.4 m/s off the system's velocity: find the cause instead of loosening the test (bug) |  | Core 1 finding at #1997: wide-pair second star 0.4 m/s of 57 km/s off since GEN.336; test tolerance widened to 1e-4. |
| GEN.238 | Wanderers: hypervelocity stars reach out to the grid, and kicked neutron stars and black holes are an ejecta population |  | Whole-sphere build train, Boss 2026-10-15. Step 6. Core. |
| GEN.321 | Emission regions D and E: the ionizing members of a young cluster in a Plummer sphere, cut by sector cube |  | Split from GEN.241. |
| GEN.530 | A pre-placed star's hosted nebula shows from its scatter row on the maps and pages |  | Residual of GEN.320: the hosted nebula of a pre-placed star from its scatter row. |
| GEN.241 | Every nebula gets the stars it needs as real systems made with the nebula, inside its bounds | GEN.321 | Whole-sphere build train, Boss 2026-10-15. Step 12. Core. Folds in MAP.182. |
| GEN.328 | Groups at fill in star-forming cores and OB regions | GEN.327, GEN.321 | Split from GEN.273. |
| GEN.273 | Young groups and trapezia: the `group` node, offset and velocity columns, a decay age and the remnant draw | GEN.327, GEN.328 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 17. Core. |
| GEN.244 | Young open clusters and OB associations as components (later) | GEN.241, GEN.321 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.298 | Remnants read the superbubble: the remnant size at build and the remnant factor use the bubble a supernova sits in |  |  |
| GEN.365 | Hot white dwarf and planetary-nebula central star radius disagrees with the radius that its luminosity and temperature imply (bug) |  | Finding from Core 3. |

### Globular clusters

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.325 | Planet cull in dense components by the survival formula and the metallicity factor |  | Split from GEN.162. |

### Whole-sphere audits

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.184 | Whole-sphere tests, nebula and Wolf-Rayet side: contents per class, the Wolf-Rayet count, tiling, determinism and the King solver | GEN.243 | Split from TEST.159. |
| TEST.161 | Fill-stage and system-pass tests: stream independence, same seed same system, and the probability table |  | Phenomena catalogue, Boss 2026-10-15. Step 2 tests. Quality. |
| TEST.162 | Multiplicity audit: companion states, row bookkeeping, stability, determinism, no planets in a scatter, fill reuse, motion and the census anchors | GEN.315 | Multi-star systems at scatter time, Boss 2026-10-15. Step 9. Quality. |
| TEST.164 | Catalogue type audit: type counts, region gates, the N tail, interaction consistency and group decay | GEN.269, GEN.270, GEN.271, GEN.272, GEN.275, GEN.327, GEN.329 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 20. Quality. |
| TEST.165 | Re-measure the region query on a real generated galaxy, warm and cold |  | Spherical region query, stage C. Step 6. Quality. |

### Composition proofs

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.211 | A star composition sweep: every kind of star, remnant and brown dwarf, fractions, trends and determinism | GEN.440 | Boss 11:58 PM EDT; the proof of the star composition. |
| TEST.215 | The Solar System calibration rows live in test files only: a flag fences them and a test fails when a generator table reads one | GEN.472, TEST.211 | Boss 1:44 AM EDT, Stand-in only: R5, calibration rows fenced into tests. |
| TEST.214 | Every registered possibility table: for every class the gated weights add up to 100 percent, no class has an empty gate, and every row has a source and a licence | GEN.472, GEN.410, GEN.427, GEN.439 | Boss 1:40 AM EDT: the gated weights add up to 100 percent; registry walk, so later tables are covered. Core 4 proposed. |

## Phase 5: New scatters: core, halo, interface update and docs

The realistic galactic core and the halo tier, the interface update for the spherical grid, the headline item and the docs that follow.

### Galactic core and halo

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.259 | Intermediate-mass black hole share gives about 100,000 in the Galaxy against a handful of candidates (bug) |  | Phenomena catalogue, Boss 2026-10-15. Question for Boss. Core. |
| GEN.282 | Galactic core gas and young terms: the central molecular zone, the nebula rate and young clusters (Arches, Quintuplet, the central parsec) | GEN.244 | Core 3; step 3 of the galactic core note |
| GEN.377 | Scatter limits by place for the neutron star and black hole scatter below the cut: the core and halo tiers apply to it, not one global mass cut |  | Core 3 finding after GEN.333. |
| GEN.303 | Per-cell system cap as a configuration option: `generation.cell_cap`, default 5,000, lowest 500, highest 50,000 |  | Core 3 wires it; Boss 1:47 PM EDT. |
| GEN.284 | Planet survival at the galactic core: a metallicity factor and a survival formula instead of the hard cut | GEN.325 | Core 3; step 4 of the galactic core note |

### Interface update

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.243 | The Galaxy Map cannot select or enter anything in the halo: the halo's scatter objects are not pickable at galaxy scale and a zoom toward one stops without opening a view (bug) | MAP.239 | Boss 9:21 PM EDT; bug; Core 2 proposed. |
| UX.115 | Interface update: every page and control works with the spherical grid, 96-bit IDs, the new phenomena and the region queries (stage after the sphere and generation stages) | NAV.61, MAP.243 | Interface; queued after sphere stages (Boss, 12:52 AM EDT) |
| TEST.168 | Control sweep: a browser test drives every control on every page against a spherical galaxy with 96-bit IDs | UX.115 | Quality; queued after UX.115 |

### Docs

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.19 | After UX.102: document the log output standard for the CLI, the Generate page and admin logs | UX.102 |  |
| DOC.20 | After PERF.86: rewrite the queue and scatter notes and the queue docstrings | PERF.86 |  |
| DOC.24 | Bring Python 3.12, numba and astropy into the research notes and three source comments |  |  |
| DOC.48 | Docs: the reach, the halo, the profile rules and the Wolf-Rayet and nebula-contents rules in the design notes | GEN.241, GEN.243 | Whole-sphere build train, Boss 2026-10-15. Step 14. Quality. |
| DOC.50 | Docs: the schema, the multi-star note and the catalogue say where the members live | GEN.264, TEST.164, GEN.130 | Multi-star systems at scatter time, Boss 2026-10-15. Step 10. Quality. |
| DOC.51 | Docs: the cell key in the coordinate system and schema docs, the region note moves to built, the API states the cap and the ladder |  | Spherical region query, stage C. Step 8. Quality. |
| DOC.52 | A docs index page: links to every HTML page under docs/, each with a description of what it contains |  | Docs index. Quality. |
| DOC.67 | An architecture map Boss can drill down: project areas, modules and functions with the call flow between them, generated from the code |  | Architecture map. Core 5 (recommended). |
| DOC.53 | The tuning module header names the rule for a new signed constant: the fuzz tests accept only a dict or a positive "fall" |  | any Core lane; Core 2 finding |
| DOC.18 | Last pass: re-run the documentation sweep when the other DOC items are done and Phase 1 is complete | DOC.19, DOC.20, DOC.24, DOC.48, DOC.50, DOC.51, DOC.52, DOC.53 |  |
| DOC.60 | Settings note: read sections A2 and A6 of the settings and accounts design against the built one-file settings and fix what is stale |  | Settings lane report, 2026-10-16; docs only. Any lane. |

### Globular clusters

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.161 | Bright-first fill for cluster sectors |  | Globular-cluster chain. |
| GEN.326 | Blue stragglers in clusters: N_BSS scales as the core mass to the 0.4 | GEN.161 | Split from GEN.162. |
| GEN.162 | Planet cull and blue stragglers in clusters | GEN.325, GEN.326 | Globular-cluster chain. |
| GEN.163 | Type-B pulsar planets in globular clusters (GEN.130 follow-on) | GEN.317 | Globular-cluster chain. |
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
| GEN.401 | The shape of the stars from the stellar mass: disk length and height, thick disk, bulge, bar and arms replace the generate defaults |  | GEN.368 build, Boss approved the table 10:11 PM EDT; Core 1 proposed for GEN.400 to GEN.403, Core 3 or Core 4 for GEN.404 and GEN.405 (profile.py, activity.py), any lane for GEN.406. |
| GEN.402 | The central black hole mass and the nuclear cluster follow the bulge: the M-bulge and M-sigma relations replace the fixed ranges | GEN.401 | GEN.368 build, Boss approved the table 10:11 PM EDT; Core 1 proposed for GEN.400 to GEN.403, Core 3 or Core 4 for GEN.404 and GEN.405 (profile.py, activity.py), any lane for GEN.406. |
| GEN.403 | The globular cluster count follows the halo mass: the system mass from the halo mass replaces the fixed count of 157 |  | GEN.368 build, Boss approved the table 10:11 PM EDT; Core 1 proposed for GEN.400 to GEN.403, Core 3 or Core 4 for GEN.404 and GEN.405 (profile.py, activity.py), any lane for GEN.406. |
| GEN.404 | The Sun radii are removed: every use of a home or Sun radius becomes a calculation at the object's own position | GEN.401 | GEN.368 build, Boss approved the table 10:11 PM EDT; Core 1 proposed for GEN.400 to GEN.403, Core 3 or Core 4 for GEN.404 and GEN.405 (profile.py, activity.py), any lane for GEN.406. Core 4 proposed (coordinator 10:16 PM EDT). |
| GEN.405 | Metallicity and rates by position from the galaxy's own model: the mass-metallicity relation and the local rates replace the fixed gradient and the quoted-at-the-Sun values | GEN.404 | GEN.368 build, Boss approved the table 10:11 PM EDT; Core 1 proposed for GEN.400 to GEN.403, Core 3 or Core 4 for GEN.404 and GEN.405 (profile.py, activity.py), any lane for GEN.406. Core 4 proposed (coordinator 10:16 PM EDT). |
| GEN.531 | Boss's research requests from the galaxy laws research (R1 to R15): the observed stellar halo, the dwarf range, the thick disk fits and the cluster system |  | Research requests R1 to R15 of the GEN.406 note; does not gate GEN.301. |

## Phase 6: The headline closes

The last audits and the headline: they test the component pass, the core and the halo, so they come after phases 4 and 5.

### The headline closes

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.185 | Whole-sphere tests, component side: pre-placed share per population, no double count of a cluster star, nothing hosted in a halo or old cluster cell, D and E members inside their radius | GEN.321 | Split from TEST.159. |
| TEST.159 | Whole-sphere tests: required stars per nebula kind, Wolf-Rayet counts, tiling at the reach, profile and empty-sector determinism, and the pre-placed share per population | TEST.184, TEST.185 | Whole-sphere build train, Boss 2026-10-15. Step 15. Quality. |
| TEST.160 | Representation audit: seeded sectors in each region class show every kind at its expected frequency, with systems and required contents checked | GEN.282, GEN.321 | Phenomena catalogue, Boss 2026-10-15. Step 4. Quality. |
| TEST.196 | Close-binary check: generated close-binary fractions by primary mass and [Fe/H] bin against the published ranges, with out-of-range bins flagged |  | Core 4; Boss 6:15 PM EDT: expected-numbers check for GEN.103. |
| GEN.301 | Headline of 8.1: the new generate scatters, accurate galactic generation and scattering including the halo | GEN.238, GEN.241, GEN.243, GEN.244, GEN.259, GEN.282, GEN.284, GEN.298, TEST.159, TEST.160, TEST.196, GEN.447, GEN.448, GEN.449, GEN.450, GEN.463, GEN.467, GEN.469, GEN.502, GEN.528, GEN.529, GEN.530 |  |
