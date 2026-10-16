# Version 8.1: Infinite zoom and the new generate scatters

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.1

## Goal

Two headlines. Infinite zoom: one camera flies from the whole galaxy to a star system (MAP.146). The new generate scatters: accurate galactic generation and scattering out to the 262 kpc reach, including the halo, the satellites, the nebula stars and a realistic core (GEN.301). The version also holds the Phase 1 tail: the queue steps, reproducible galaxies, the description layer, tech levels, directives and every open bug and flaky test. The 8.1 stamp waits until all of it is done (REVISION_HOLD in scripts/bump_version.py).

## Phase 1: Foundations and open trees

Everything with no blocker outside its own thread: the queue steps, reproducible galaxies, the description layer, tech levels, the settings page (moved up from 8.2, Boss 1:43 PM EDT), directives, housekeeping, the bugs and the backfill weeding.

### Queue, caches and performance

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.18 | Migration helpers for slow DDL: online indexes, instant columns and batched updates |  | Research: slow DDL for the big tables. |
| OPS.15 | Each update says whether it changes generated output |  | Needs the fingerprint, so phase 2. |
| PERF.18 | Run the GEN.30 bright-star backfill in parallel on the work queue |  | Backfill blocks become RQ jobs. Same stars as the one-process backfill for one seed needs GEN.39. |
| PERF.20 | Short-term caching through the work queue and API (needs planning) |  | Plans on RQ results and the new caches. Plan with Boss. |
| PERF.29 | Record which runs a partly filled sector still needs |  | Builds on the generated flag that GEN.76 settles (sector_stats.bright_level_sol). |
| PERF.30 | Finish an interrupted block or sector run on the next start | PERF.29 |  |
| PERF.36 | Memory and request guard: never list more than about 50,000 candidate cells, and refuse huge enumerations in a web request |  | Research: protects the Generate page estimates. |
| PERF.41 | Stamp the tile and page caches with a TILE_FORMAT constant instead of the version (optional) |  | Research follow-up to PERF.25 (built). |
| PERF.47 | The PERF.31 benchmark records the buffer pool, table sizes and worker start-up cost |  | Generation performance study. |
| PERF.48 | Low priority: a numeric-only INSERT formatter or C driver for bright_stars and phenomenon_scatter |  | Generation performance study. |
| PERF.75 | Keyset paging for the other data tables and the jump-by-value box (Systems name sort done) |  | Decided by default. Foundations lane 1. |
| PERF.86 | Redo the work queue: drop the task and subtask logic and use Redis simply |  |  |
| PERF.88 | Sector-path integrator on numba: saving the sector paths takes 13 s instead of 34 s on the benchmark fill |  |  |
| PERF.91 | Bulk random draws for the scatter layers (4 to 12 times on batches above 1,000 draws) |  |  |
| PERF.97 | Work queue step 3: the mass scatter end to end on the new queue |  | PERF.86 step 3 |
| PERF.98 | Work queue step 4: the remaining stages and the sector and backfill units on the new queue | PERF.97 | PERF.86 step 4 |
| PERF.99 | Work queue step 5: the Generate page, Queue page and job log read Redis; cancel now through the stop command | PERF.98 | PERF.86 step 5 |
| GEN.231 | Generate a whole block on the shell ladder: the generate --block mode and the block buttons come back | PERF.99 | Maps hand-off 2026-10-15. Core, then Interface for the buttons. |
| PERF.100 | Work queue step 6: delete the old machinery and the work_* tables (Alembic migration) | PERF.99 | PERF.86 step 6 |
| PERF.105 | Slow page requests queue a job that builds the answer into Redis, with a status route (the mechanism) |  | Boss 2026-10-16, 12:21 PM EDT; Core 2, after MAP.152 and MAP.146. |
| UX.118 | Slow pages show a progress bar while their data loads: the search page and the phenomena list pages | PERF.105 | Boss 2026-10-16; Core 2, after PERF.105. |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.50 | Remove the galaxy settings files and everything around them; keep only an internal record of the current galaxy |  | Boss 2026-10-10 21:14Z; Foundations lane 1, before OPS.28. |
| GEN.55 | Same seed, same data: a sector's contents depend only on the seed, the version and its address (internal) |  | Parent of the chain. |
| GEN.57 | A sector's contents depend only on the seed, the version and its address |  | Name collisions go away with GEN.67; GEN.63 dropped. Also keys GEN.64's ID collision counter on address (it checks stored rows, PR #406). |
| DB.9 | Repair a damaged galaxy database from a parity file | GEN.57 | Boss 02:13Z: phase 1. Reed-Solomon parity over groups of sector exports; seed regeneration as fallback. |
| DB.17 | Repair by regenerating a damaged sector from its seed when parity cannot rebuild it | DB.9, GEN.57 | Research split of DB.9: the regenerate-from-seed fallback. |
| GEN.135 | Deterministic math helpers for every stored float, with a lint and frozen constants |  | Research: Python 3.9 and 3.10 to 3.13 draw different sectors without it. Blocks TEST.77. |
| GEN.136 | Fingerprint encoding: floats to 9 significant digits, a stored leaf digest per sector and a ring-and-layer tree | GEN.135 | Research follow-up to GEN.58 (built). |
| OPS.28 | Generator epoch and battery digest: say whether two checkouts generate the same galaxy | GEN.135 | Research: OPS.13 is built; this is the epoch it records. OPS.14, OPS.15, OPS.12 and TEST.77 read it. |
| DB.16 | Store the generator epoch and run id on each sector instead of four version text columns | OPS.28 | Research follow-up to DB.7 (built). |
| OPS.29 | Update reload: cover the gunicorn units and non-Apache hosts, and say what a reload aborts |  | Research follow-up to OPS.8 (built). |
| OPS.37 | A Generator version number: one sequential integer, shown in admin and the API | OPS.28 | Boss 2026-10-09 20:59Z: done by the end of phase 1. Same number as OPS.28's epoch (decided). |
| TEST.77 | A golden-seed regression test | GEN.57, GEN.135, OPS.28 |  |

### Tech levels and species

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.29 | update.sh and every database migration refuse to start while jobs are active, unless --yes is passed or the user confirms cancelling them |  |  |
| DB.30 | One migration for the per-species environment limits | DB.29 |  |
| POP.14 | A race is only placed where it can survive: a placement check and a consistency check over the data |  |  |

### Description layer

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.31 | Store structure, not prose, for life and flavour |  |  |
| GEN.218 | One habitability describer, in the short wording |  |  |
| GEN.219 | The habitability levels page and tooltips are generated from the limits |  |  |
| GEN.220 | "Habitable" means the index | GEN.218 |  |
| GEN.221 | Short plain prose for the fixed text |  |  |
| TEST.155 | Wording and index agree: property test | GEN.218 |  |
| UX.105 | The system page builds every body's detail from the document model as HTML tables |  |  |
| DOC.45 | Reference docs follow the description layer | UX.105, GEN.218, GEN.219 |  |
| UX.107 | Units in the description go through the ladders |  |  |
| UX.108 | Describe phenomena, sectors, galaxy and facilities in the same Document |  |  |
| UX.109 | Remove the old per-class description text | UX.105, UX.108 |  |
| UX.110 | The phenomena count on the overview page is wrong (bug) |  |  |

### Admin and API housekeeping

| ID | Item | Needs | Note |
|---|---|---|---|
| API.4 | API compatibility data in the docs |  | Docs and version number. |
| API.7 | Investigate and plan upload limits |  | Plan only. |
| API.15 | Log every API call with its user, how it came in, and its HTTP response code |  | Needs no user accounts (Boss 01:31Z). |
| API.26 | Bump API_VERSION for the shell, band and wedge sector address |  |  |
| OPS.35 | A vendored-version lock file for the static libraries |  | Research: supply-chain record. |
| UX.49 | Form fields as Shoelace components (sl-input, sl-select, sl-checkbox) across the site |  | Follows UX.40 (done): the fields themselves become Shoelace components. |
| UX.102 | One logging standard for every output the user or an admin reads: the same timestamps, format, dividers and step markers everywhere |  |  |

### Settings page

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.87 | Backup archive: one bz2 tarball per backup holding the settings, the control database and the parity data as parts, with retention and restore of the whole or one part |  | Boss 1:59 and 2:00 PM EDT: one bz2 tarball per backup; retention settings.backups_kept (default 20). |
| ADM.83 | Settings writer: validate, write and back up the settings overlay, with an audit row per change | ADM.87 | Boss 1:43 PM EDT: the configuration page moved up from 8.2 (GitHub issue [#515](https://github.com/dwhagar/planetGen/issues/515)). |
| ADM.84 | Settings page under Admin: every option with its help, edit, validate and save, with a diff before saving | ADM.83 |  |
| ADM.85 | Settings backups, restore, import and export on the settings page | ADM.84, ADM.87 |  |
| ADM.43 | A full configuration page under Admin | ADM.83, ADM.84, ADM.85, ADM.87 | GitHub issue [#515](https://github.com/dwhagar/planetGen/issues/515). Umbrella; the Owner tier is ADM.86 in 8.2. |
| ADM.44 | Web, Open Graph and SEO settings | ADM.43 | GitHub issue [#743](https://github.com/dwhagar/planetGen/issues/743). Moved from 8.2 with ADM.43. |

### Generate page and directives

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.179 | Store each sector's generation directive and attempt record with the sector |  | Left over from GEN.96. |
| GEN.180 | Directives: a forced fill (met_forced) after K failed draws |  | Left over from GEN.96. |
| GEN.181 | Directives: refuse impossible requests up front from the compound-Poisson tables |  | Left over from GEN.96. |
| GEN.229 | Repair or re-roll the wide binaries that already break the stability rule |  |  |
| UX.3 | Warn every visitor while a background job changes the galaxy |  | ETA from the RQ job's published progress. ETA from progress.json, which PERF.23 caps. |

### Bugs and flaky tests

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.33 | MySQL 8.0.46 reflects the foreign keys fk_facilities_type and fk_facility_field_values_field as RESTRICT, so check-db calls a healthy database damaged (bug) |  | Quality run 2026-10-15, MySQL 8.0.46. Core. |
| TEST.111 | test_ensure_sector_generated_creates_then_reuses_the_same_sector fails in a busy parallel run (bug) |  | Test flake reported 2026-10-09; Bugfixes lane 1. |
| TEST.116 | test_bughunt_end_to_end stores M2V or M6V where it expects K2V on the MySQL 8.4 and MariaDB legs once (bug) |  | Bugfixes lane 1 report 03:57Z. |
| TEST.126 | The heavy WebGL browser-map tests time out when four workers run them together (bug) |  | Lane 2 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. |
| TEST.127 | test_nebula_shape_endpoint_serves_a_mesh failed once in a full parallel run (bug) |  | Lane 2 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. |
| TEST.128 | test_spatial_position_db and test_web_db_fields fail under a parallel run on one MySQL (bug) |  | Lane 1 report 2026-10-10 20:50Z; Bugfixes lane 2 after its current items. |
| TEST.129 | Two orbit-update tests fail on the MySQL 8.4 leg: the stored last_updated_at rounds up (bug) |  | Lane 1 report 2026-10-10 20:55Z; Bugfixes lane 1 after its current items. |
| TEST.130 | test_the_scene_positions_at_the_epoch_match_the_stored_ones fails when a random system has comets (bug) |  | Lane 1 report 2026-10-10 20:55Z; Bugfixes lane 1 after its current items. |
| TEST.134 | test_many_processes_starting_at_once_start_one_job fails under a parallel run (bug) |  | Flaky under load. Bugfixes lane 1. |
| TEST.137 | test_a_neighborhood_from_the_generate_page[core] fails once under load (bug) |  |  |
| TEST.138 | test_open_map_menus_hold_no_overlap[web.galaxy-1280] fails once in a full parallel run (bug) |  |  |
| TEST.140 | test_a_body_is_reached_by_its_id_and_a_row_number_names_nothing fails once in a full parallel run (bug) |  |  |
| TEST.141 | test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector[15-True] fails once in a full parallel run (bug) |  |  |
| TEST.144 | test_every_plain_button_does_something[web.sector] failed once in a full parallel run (bug) |  |  |
| TEST.148 | Investigate resuming a test run where it left off when the build is the same |  |  |
| TEST.156 | Rewrite the tests that name ring, layer and slot |  |  |
| TEST.167 | test_random_start_mode_respects_min_start_density hit a 900 s timeout in a full run and passes alone (bug) |  | Test failure under load. Quality. |
| TEST.170 | test_math_check::test_whole_run_is_fast and the Galaxy Map breadcrumb folding test fail under full-suite load and pass alone (bug) |  |  |
| TEST.171 | test_gen_resume::test_sectors_a_forced_scatter_skipped_fill_correctly_afterwards failed once under gate load with no scattered system (0 > 0) and passes alone (bug) |  |  |
| TEST.172 | test_galaxy_map_forward_to_current_jumps_to_the_newest_view failed once: clicking a Galaxy Map block with 0 sectors charted did nothing, and it passes alone (bug) |  |  |
| TEST.173 | test_gen_bright_scatter_edges::test_a_fill_after_a_shell_failed_mid_scatter_builds_no_leftover_star failed once under gate load with an empty ring in _nearest (bug) |  |  |
| TEST.174 | test_web_a11y::test_page_is_accessible_and_fits[web.sector-desktop-light] failed once under gate load and passes alone (bug) |  |  |

### Backfill weeding

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.42 | A pass that drops sectors from a region by probability |  | Same function as PERF.18 (backfill_bright_stars_around); one thread. |
| GEN.40 | Weed out sectors by star density before the bright-star backfill | GEN.42 | Parent; closes with its subitems. |
| PERF.35 | An interval or chunk ledger for untouched sectors once block-first backfill lands | GEN.42 | Research: ledger size at full galaxy. |

## Phase 2: Infinite zoom

The fly-through Galaxy Map (MAP.146) and the nebula data its tiles carry.

### Fly-through map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.75 | The mini map as a second engine view |  | Locked second camera on MAP.64. |
| MAP.59 | Make it plain that a zoomed-in slab is a slab, not a wedge | MAP.75 |  |
| MAP.121 | Every map shows and steps to its neighbouring regions, on one map engine |  | Folds MAP.127, the Galaxy Map's neighbouring blocks and slabs (Boss 2026-10-07 16:27Z: one map engine for all maps). |
| MAP.122 | A Select mode on every galaxy view: Galaxy (blocks and sectors) or Star |  | Merges the 2026-10-07 "button to select a star" ask; MAP.101 (PR #431) made stars unpickable by default. |
| MAP.139 | The Galaxy View uses its spare space: an info box with a Details link, and menu items |  | GitHub issues [#758](https://github.com/dwhagar/planetGen/issues/758) and [#715](https://github.com/dwhagar/planetGen/issues/715), one layout change. |
| MAP.140 | Double-click on a selected object goes there and opens its information |  | GitHub issue [#714](https://github.com/dwhagar/planetGen/issues/714). |
| MAP.141 | Context around the selection: faint neighbours, and the sectors above and below |  | GitHub issues [#718](https://github.com/dwhagar/planetGen/issues/718) and [#716](https://github.com/dwhagar/planetGen/issues/716), one context view. The #716 bug label was overruled. |
| MAP.147 | The Galaxy Map wire format: the investigation (done) and the record of what was built from it |  | Research Lane 3 report in docs/design/galaxy-map-wire-format.md; record until the build items are done. |
| MAP.152 | Scale hand-offs: galaxy, sector and system cross-fade with hysteresis, and per-tile camera-relative origins |  | Fly-through report item 5. |
| MAP.146 | Fly through the galaxy: scroll-zoom, double-click flight, distance-based visibility and a see-through near field | MAP.152 | Boss 2026-10-09 22:24Z and 22:56Z. Umbrella for MAP.148 to MAP.152; replaces the fixed ladder in galaxy/drill.py and galaxyprisms.js. |
| MAP.156 | View one layer or a range of layers top-down from the galaxy view, as a secondary option |  | Boss 2026-10-09 23:17Z; secondary option, not a main control. |
| MAP.160 | Quantise the Galaxy Map GPU buffers (deferred) |  | Deferred; MAP.147 recommendation step 3. |
| MAP.171 | The Galaxy Map keys a massive phenomenon's visibility to its map luminosity (replaces "always lit") | MAP.152 | Boss 2026-10-10 21:38Z. Foundations lane 2 after MAP.152. |
| MAP.189 | The Galaxy Map opens flown in to the created sectors | MAP.152 |  |
| MAP.190 | Several separate regions of created sectors: which one the Galaxy Map opens on | MAP.189 |  |
| MAP.191 | Galaxy Map default view: no density shading, only nebulae shade the stars behind them |  |  |
| MAP.195 | Galaxy Map finer tiles leave out the stars a coarser tile already sent |  |  |
| MAP.209 | Opening a sector from the Galaxy Map flashes one view, reloads into another, and leaves the details panel empty (bug) |  |  |
| MAP.210 | Nothing on the map can be selected, even in sector view (bug) |  |  |

### Nebula data on the Galaxy Map

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.202 | Scattered nebulae (planetary nebulae and remnants) join the nebula index, deduplicated by object id |  | Foundations lane 1. |
| MAP.172 | Nebula cell aggregates in the region pyramid |  | Foundations lane 1. |
| MAP.173 | The Galaxy Map tile carries a nebula layer: regions and singles, nested like the star lists | MAP.172 | Foundations lane 2. |
| MAP.174 | Nebula cover in the stage-view cell statistics, for filled and unfilled cells |  | Foundations lane 2. |
| TEST.135 | Big-galaxy query-budget test for the nebula reads |  | Foundations lane 1. |


## Phase 3: New scatters: stars, nebulae and the region query

The multi-star catalogue, the stars every nebula needs and the region query the whole sphere reads through.

### Multi-star and compact systems

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.129 | Multi-star systems of up to twelve stars |  | GitHub issue [#777](https://github.com/dwhagar/planetGen/issues/777). Build after the GEN.128 design. |
| GEN.130 | Compact-remnant systems at sector fill: companions and pulsar planets for neutron stars and black holes |  | Phenomena catalogue, Boss 2026-10-15. Step 3, rescoped. Core. |
| GEN.252 | Brown dwarfs as stars: an L, T and Y state for companions and primaries, with the brown-dwarf desert |  | Phenomena catalogue, Boss 2026-10-15. Step 5, later. Core. |
| GEN.263 | Promote every bound compact pair to the scatter, whatever its mass |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 4. Core. |
| GEN.264 | The fill builds a system from its stored members; the system pass keeps only planets, moons and belts |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 5. Core. |
| GEN.267 | A MAX_STARS setting and a --stars N directive, with a prevalence knob for five or more stars |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 11. Core. |
| GEN.269 | Interaction states in the star model: a Roche classifier, `interaction` on every pair, member state and flags, rejuvenated age and the new kinds |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 13. Core. |
| GEN.270 | Close-binary types: contact, Algol, spectroscopic and eclipsing flags, tertiary shares, merger products and Be stars | GEN.269 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 14. Core. |
| GEN.271 | White-dwarf binary types: close double white dwarfs, post-common-envelope pairs, cataclysmic variables, AM CVn, symbiotic and barium stars | GEN.269 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 15. Core. |
| GEN.272 | Accreting compact binaries: X-ray binaries, accretion luminosity, microquasars, spider pulsars and compact triples | GEN.269 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 16. Core. |
| GEN.273 | Young groups and trapezia: the `group` node, offset and velocity columns, a decay age and the remnant draw |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 17. Core. |
| GEN.274 | Cluster-hosted compact binaries and runaway bookkeeping: `runaway_origin` and the compact companion left by a supernova |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 18. Core. |
| GEN.275 | Substellar and circumbinary types: brown-dwarf and planetary-mass companions, circumbinary disks and planets |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 19. Core. |

### Nebula and young stars

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.242 | A pre-main-sequence state for young stars, so protostars can exist as real stars |  | Whole-sphere build train, Boss 2026-10-15. Step 11. Core. |
| GEN.241 | Every nebula gets the stars it needs as real systems made with the nebula, inside its bounds | GEN.242 | Whole-sphere build train, Boss 2026-10-15. Step 12. Core. Folds in MAP.182. |
| GEN.243 | Wolf-Rayet stars and ring nebulae: a W state for the last 0.3 Myr of very massive stars, and nebula class X | GEN.241 | Whole-sphere build train, Boss 2026-10-15. Step 13. Core. |
| GEN.258 | Protostars embedded in giant molecular clouds | GEN.241, GEN.242 | Phenomena catalogue, Boss 2026-10-15. Step 8, later. Core. |

### Region query

| ID | Item | Needs | Note |
|---|---|---|---|
| API.27 | A read-only profile lookup for any address in the reach, so an uncharted sector can be previewed |  | Whole-sphere build train, Boss 2026-10-15. Step 8. Interface. |
| GEN.276 | planetgen/galaxy/cellkey.py: the cell key function and the sphere-to-key-runs cover, numpy only |  | Spherical region query, stage C. Step 2. Core. |
| NAV.45 | "What's within N pc" from the Galaxy Map and Sector Map |  | Boss 04:39Z. Map action opening NAV.44; drawn sphere optional. |
| NAV.60 | What is nearby reads by the cell key: one statement for every kind, including objects with no sector and unbuilt scatter rows | GEN.276 | Spherical region query, stage C. Step 3. Core. |
| NAV.61 | Nearest-first by a radius ladder, and a dense-core guard with an early-stop count | NAV.60 | Spherical region query, stage C. Step 4. Core. |
| UX.113 | The sector page shows the profile of any address in the reach, and heavy cells and far neighbourhoods warn before they run | API.27 | Whole-sphere build train, Boss 2026-10-15. Step 8. Interface. |
| UX.114 | Scatter-row pages list the stored members and the map shows a system mass |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 8. Interface. |

## Phase 4: New scatters: components, clusters and audits

Satellites, clusters and wanderers as components, the globular cluster chain and the audits that check what the scatters make.

### Components and wanderers

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.99 | Nebula volume backfill with the star types the nebula needs |  |  |
| GEN.102 | Investigate filling all near-zero-density void space at once |  |  |
| GEN.103 | Research where each star type and phenomenon belongs in the galaxy's structure |  |  |
| GEN.134 | Tune the star populations to the observed star-formation profile by galactic radius |  | From GEN.133's four unbuilt proposals (star-types-by-galactic-radius.md); Boss to confirm which he wants. |
| GEN.237 | Satellite and cluster components: one table for the Magellanic Clouds, Sagittarius, dwarf galaxies and globular clusters |  | Whole-sphere build train, Boss 2026-10-15. Step 5. Core. |
| GEN.238 | Wanderers: hypervelocity stars reach out to the grid, and kicked neutron stars and black holes are an ejecta population |  | Whole-sphere build train, Boss 2026-10-15. Step 6. Core. |
| GEN.240 | The three scatters fill the whole sphere with one algorithm: band groups, the phenomena pass on the object-first sampler, and a component pass | GEN.237 | Whole-sphere build train, Boss 2026-10-15. Step 6. Core. |
| GEN.244 | Young open clusters and OB associations as components (later) | GEN.237, GEN.241 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.298 | Remnants read the superbubble: the remnant size at build and the remnant factor use the bubble a supernova sits in |  |  |

### Globular clusters

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.161 | Bright-first fill for cluster sectors |  | Globular-cluster chain. |
| GEN.162 | Planet cull and blue stragglers in clusters |  | Globular-cluster chain. |
| GEN.163 | Type-B pulsar planets in globular clusters (GEN.130 follow-on) | GEN.130 | Globular-cluster chain. |

### Whole-sphere audits

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.159 | Whole-sphere tests: required stars per nebula kind, Wolf-Rayet counts, tiling at the reach, profile and empty-sector determinism, and the pre-placed share per population | GEN.241, GEN.243 | Whole-sphere build train, Boss 2026-10-15. Step 15. Quality. |
| TEST.160 | Representation audit: seeded sectors in each region class show every kind at its expected frequency, with systems and required contents checked |  | Phenomena catalogue, Boss 2026-10-15. Step 4. Quality. |
| TEST.161 | Fill-stage and system-pass tests: stream independence, same seed same system, and the probability table |  | Phenomena catalogue, Boss 2026-10-15. Step 2 tests. Quality. |
| TEST.162 | Multiplicity audit: companion states, row bookkeeping, stability, determinism, no planets in a scatter, fill reuse, motion and the census anchors | GEN.264 | Multi-star systems at scatter time, Boss 2026-10-15. Step 9. Quality. |
| TEST.164 | Catalogue type audit: type counts, region gates, the N tail, interaction consistency and group decay | GEN.269, GEN.270, GEN.271, GEN.272, GEN.273, GEN.274, GEN.275 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 20. Quality. |
| TEST.165 | Re-measure the region query on a real generated galaxy, warm and cold | GEN.276 | Spherical region query, stage C. Step 6. Quality. |
| TEST.166 | The cell key follows an object across a cube boundary, and the point columns and the key round-trip bit for bit |  | Spherical region query, stage C. Step 7. Quality. |

## Phase 5: New scatters: core, halo, interface update and docs

The realistic galactic core and the halo tier, the interface update for the spherical grid, the headline item and the docs that follow.

### Galactic core and halo

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.24 | Generate the galactic core on layer 0 |  | Bulk core fill runs on the parallel path; new mode on generate.html. |
| GEN.210 | The galactic core is one round disc sector in the middle of the grid, 4 pc across and 4 pc tall |  |  |
| GEN.259 | Intermediate-mass black hole share gives about 100,000 in the Galaxy against a handful of candidates (bug) |  | Phenomena catalogue, Boss 2026-10-15. Question for Boss. Core. |
| GEN.280 | Galactic core density: nuclear star cluster and nuclear disc terms in the density model, replacing the Gaussian core of ADM.49 |  | Core 3; step 1 of the galactic core note; needs a Core 1 migration revision |
| GEN.282 | Galactic core gas and young terms: the central molecular zone, the nebula rate and young clusters (Arches, Quintuplet, the central parsec) | GEN.237, GEN.244 | Core 3; step 3 of the galactic core note |
| GEN.283 | Scatter limits by place: a core luminosity floor and mass limit, a halo luminosity floor (default 500, lowest 50), and the baseline everywhere else |  | Core 3; step 3b of the galactic core note; shares one Core 1 migration with GEN.281 |
| GEN.303 | Per-cell system cap as a configuration option: `generation.cell_cap`, default 5,000, lowest 500, highest 50,000 | GEN.280, ADM.84 | Core 3 wires it; Boss 1:47 PM EDT. |
| GEN.284 | Planet survival at the galactic core: a metallicity factor and a survival formula instead of the hard cut | GEN.280, GEN.162 | Core 3; step 4 of the galactic core note |
| GEN.285 | Nuclear mass components in the galactic potential for the orbital update | GEN.280 | Core 3; step 5 of the galactic core note |
| GEN.286 | Horizontal-branch and AGB phases in the stellar evolution model |  | Core 3; HB and AGB phases (Boss, 1:12 AM EDT); plan in note 5.3; no migration; interim [Fe/H] defaults until GEN.158 |
| GEN.287 | The stellar evolution model spreads giant-branch stars evenly in log luminosity, which probably overstates the bright giants by about a factor of two against real halo counts (bug) |  | Core 3; folded into GEN.286 unless Core 3 says otherwise |
| UX.116 | Sliders for every scatter and fill limit: one limit_slider with a number box, in three new segments under the scatter-limits row of the Generate page | GEN.283 | Core 3; sliders for the scatter and fill limits |
| GEN.301 | Headline of 8.1: the new generate scatters, accurate galactic generation and scattering including the halo | GEN.237, GEN.238, GEN.240, GEN.241, GEN.242, GEN.243, GEN.244, GEN.259, GEN.280, GEN.282, GEN.283, GEN.284, GEN.285, GEN.286, GEN.287, GEN.298, UX.116, TEST.159, TEST.160 |  |

### Interface update

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.115 | Interface update: every page and control works with the spherical grid, 96-bit IDs, the new phenomena and the region queries (stage after the sphere and generation stages) | GEN.237, GEN.240, NAV.61, API.27, UX.113 | Interface; queued after sphere stages (Boss, 12:52 AM EDT) |
| TEST.168 | Control sweep: a browser test drives every control on every page against a spherical galaxy with 96-bit IDs | UX.115 | Quality; queued after UX.115 |

### Docs

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.4 | Correct the stale statements the research found in docs, docstrings and comments |  | Docs and comments only; no behaviour change. |
| DOC.18 | Last pass: re-run the documentation sweep when the other DOC items are done and Phase 1 is complete |  |  |
| DOC.19 | After UX.102: document the log output standard for the CLI, the Generate page and admin logs | UX.102 |  |
| DOC.20 | After PERF.86: rewrite the queue and scatter notes and the queue docstrings | PERF.86 |  |
| DOC.24 | Bring Python 3.12, numba and astropy into the research notes and three source comments |  |  |
| DOC.47 | Rebuild docs/design/architecture.html from the new architecture.md |  |  |
| DOC.48 | Docs: the reach, the halo, the profile rules and the Wolf-Rayet and nebula-contents rules in the design notes | GEN.237, GEN.241, GEN.243 | Whole-sphere build train, Boss 2026-10-15. Step 14. Quality. |
| DOC.49 | Docs: point the rate, anomaly, multi-star and nebula notes at the phenomena catalogue | GEN.130 | Phenomena catalogue, Boss 2026-10-15. Step 9. Quality. |
| DOC.50 | Docs: the schema, the multi-star note and the catalogue say where the members live | GEN.264, TEST.164 | Multi-star systems at scatter time, Boss 2026-10-15. Step 10. Quality. |
| DOC.51 | Docs: the cell key in the coordinate system and schema docs, the region note moves to built, the API states the cap and the ladder | GEN.276, NAV.60 | Spherical region query, stage C. Step 8. Quality. |
| DOC.52 | A docs index page: links to every HTML page under docs/, each with a description of what it contains |  | Docs index. Quality. |
| DOC.53 | The tuning module header names the rule for a new signed constant: the fuzz tests accept only a dict or a positive "fall" |  | any Core lane; Core 2 finding |
