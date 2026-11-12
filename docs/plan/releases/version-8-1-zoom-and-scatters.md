# Version 8.1: Infinite zoom and the new generate scatters

Core design principle (Boss, 2:58 PM EST on 2026-11-08): everything in a system affects everything else in the system, by 3D position, orbit-averaged. Items GEN.1441 to GEN.1454 build it; the rules and the register of host-only paths are in [system-interactions.md](../../design/generation/system-interactions.md).

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.1

## Goal

Two headlines. Infinite zoom: one camera flies from the whole galaxy to a star system (MAP.146, done: every sub-item merged and Lane 12 checked each clause of its Done text against main, 3:20 AM EDT 2026-10-23 as the lane reported). The new generate scatters: accurate galactic generation and scattering out to the 262 kpc reach, including the halo, the satellites, the nebula stars and a realistic core (GEN.301). The version also holds the Phase 1 tail: the queue steps, reproducible galaxies, the description layer, tech levels, directives and every open bug and flaky test. A third strand (Boss, 2026-10-16, 11:57 PM EDT): the building blocks of the 9.0 composition headline, which go down to 8.1 because they change no stored value: the research lanes and designs (phase 1, GEN.431 to GEN.437 with GEN.409, GEN.425 and GEN.426), the data tables (phase 2, GEN.438, GEN.410, GEN.427, GEN.439), the star composition function (phase 3, GEN.440) and its sweep (phase 4, TEST.211); the master note is docs/design/composition/composition-model.md. The 8.1 stamp waits until all of it is done (REVISION_HOLD in scripts/bump_version.py).

## Where it stands

As of 12:34 PM EDT 2026-10-28: 1160 of 1590 items finished (73.0 percent), 430 open, 36 dropped (TEST.374 dropped as not reproduced on main, Lane 04 10:31 PM EDT as the lane reported; GEN.564 folded into GEN.764, Lane 11 7:32 PM EDT; TEST.366 dropped as no longer reproducing after #3125, Lane 04 7:20 PM EDT as the lane reported; TEST.335 dropped as not reproduced on main, Lane 04 7:14 PM EDT as the lane reported; TEST.360 and TEST.361 dropped as not reproduced on main, Lane 04 6:35 PM EDT as the lane reported; TEST.337 and TEST.339 dropped as not reproduced on main, Bug Hunt 10:41 AM EDT as the lane reported; TEST.326 dropped as a duplicate of TEST.313, TODO 4:07 AM EDT; DB.95 dropped by TODO on lane B's finding, 2:53 AM EDT: nothing reads the view and DB.94 covers free moons; PERF.162 and PERF.166 moved to the later version, Boss 12:27 PM EDT; GEN.695 to GEN.698 dropped, Boss 1:22 AM EDT: class B never has hydrogen or helium). Open by phase: phase 1 83, phase 2 44, phase 3 24, phase 4 13, phase 5 55, phase 6 29, phase 7 17, phase 8 4, phase 9 31, phase 10 2, phase 11 131, phase 12 19 (phases 6 to 11 are the composition and sky colour suite pulled in from 8.4 and 9.0 by Boss at 5:48 PM EDT; phase 12 is the headline). Only this version is being built (Boss, 7:44 PM EDT 2026-10-22); the database chain finished with DB.85 (#2582, revision 0130, no generator epoch, so epoch 51 was still free then; since, 51 to 53 are taken: GEN.307, GEN.315 slice 1 and GEN.628, and 54 is next). Closing order: the audits TEST.185, TEST.184, TEST.159 and TEST.160 (merged, re-run as TEST.283), PERF.146 (merged as #2559) with PERF.164 and PERF.165, then SEC.33 and TEST.261 in their own threads, then GEN.301. Lane queues: docs/plan/lanes.md and reports/8-1-core-queues-2026-10-22.md.

## Phase 1: Foundations and open trees

Everything with no blocker outside its own thread: the queue steps, reproducible galaxies, the description layer, tech levels, the settings page (moved up from 8.2, Boss 1:43 PM EDT), directives, housekeeping, the bugs and the backfill weeding.

### Queue, caches and performance

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.40 | Weed out sectors by star density before the bright-star backfill |  | Parent; closes with its subitems. |
| PERF.35 | An interval or chunk ledger for untouched sectors once block-first backfill lands |  | Research: ledger size at full galaxy. |
| PERF.147 | A per-task progress bar with an ETA for any task expected to take, or taking, longer than 15 seconds |  | Boss, 8:46 PM EDT 2026-10-21 (relayed). |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.179 | Store each sector's generation directive and attempt record with the sector |  | Left over from GEN.96. |
| DB.17 | Repair by regenerating a damaged sector from its seed when parity cannot rebuild it | GEN.179 | Research split of DB.9: the regenerate-from-seed fallback. |

### Tech levels and species

| ID | Item | Needs | Note |
|---|---|---|---|

### Description layer

| ID | Item | Needs | Note |
|---|---|---|---|

### Admin and API housekeeping

| ID | Item | Needs | Note |
|---|---|---|---|


### Generate page and directives

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.283 | Update the single-system generation interface (the individual system generator page) to match everything added since it was last touched |  |  |

### Bugs and flaky tests

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.344 | Gate: browser and WebGL tests run as their own serial phase, with a worker cap by free memory |  |  |
| TEST.190 | A failing session fixture makes every dependent test wait its full 900 s timeout: fail the dependents at once (bug) |  | Test-infrastructure finding. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.237 | The wheel-anchor browser test failed alone at box load 9.7 (call 5.47 s) after the TEST.231 fix merged (bug) |  |  |
| GEN.396 | Revisit the Yerkes class III luminosity range: 50 to 1000 solar luminosities sits above many dim red giants |  | Lane 13 finding, 2026-10-16 10:00 PM EDT; any Core lane. |
| GEN.1400 | Class-letter leftovers: m_count, the unused trojan_bodies and hycean_candidate letter code, the dead pass-throughs in bodylook and habitability_world, the stale tuning docstring (audit F7 and F8) |  | Class-last audit, train 18. |
| GEN.1403 | Water-limited weathering: the carbon dioxide drawdown scales with land times the share of that land that is wet, instead of land alone (the GEN.909 option 3 lever) |  | Class M rarity report, 2026-11-08. |
| UX.289 | Any step that takes more than 5 seconds shows a progress bar from the 5-second mark: one shared mechanism for generation stages, the uncharted-sector index, the upper-level roll-up and every long step |  | Boss, 2026-11-08. |
| UX.291 | Phone width: the star-system map card is 338 px high against the half-window limit of 329 px because its three rows wrap (bug) |  | Lane 03, 2026-11-08. |
| UX.292 | At phone width an expanded star-system card with long rows scrolls itself instead of its rest region: the first view alone fills the card's maximum height (bug) |  | Lane 03, 2026-11-08. |
| GEN.1387 | Barred world redraw rate is 22.75 per system against a ceiling of 0.0183; barred_world_limit is 1 and 2 in test_barred_world[104] and [160] (bug) |  | Bug Hunt first pass, 2026-11-07. |
| TEST.695 | Stale or failing tests, group one: point_columns (sky_tier_draws), column round trip (irregular_moons), materials table, command docs (link-neighbors), TODO codes shown, JS unit and raw formatting, web_formulas, multistar registry MS4 and MS6, sql portability, db boundary values, update_database_step, literal insert, old_urls belt_ref, parent_link, container_plan, nebula_shape, object_rows, outline_reach, describe_render, property_sheet, api_link_ids (bug) |  | Bug Hunt first pass, 2026-11-07. |
| GEN.1390 | Generation sanity reds: moon history, class_sample P, daylight, prevalence habitable_world, climate sweep, atmosphere_state, biogenic oxygen, orbit_elements, life_schedule, population_streams, moon_disk_zone, sweep_bulk, forced_star (a subdwarf comes out V), invalid slot type (bug) |  | Bug Hunt first pass, 2026-11-07. |
| TEST.696 | Browser test failures: Help menu expectation (9 of 10), plane view and side strip clicks 'not visible', edit menu dialog never opens, fixture maps (9 of 227), generate flow phone-light, map card, density, hover, key, galaxy closed (bug) |  | Bug Hunt first pass, 2026-11-07. |
| PERF.254 | Memory governor for the generation run: admit and shed workers from the live free memory, learn each stage's worker size, and measure the cut of a stage as its own step |  | Boss, 2026-11-08. |
| DOC.114 | Full refresh of all documents of all types against main: READMEs, design notes, plan files, reports, change notes, generated pages and help text |  | Boss, 2026-11-08. |
| DOC.115 | One unified planning source: the plan files consolidated so that docs/TODO.md stays the item list and every other plan file points at it, with item numbers and links kept |  | Boss, 2026-11-08. |
| DOC.117 | History sweep: read every commit, item, choice, issue and pull request, alive or dead, and add past choices, rejected paths and trends to the design documents |  | Boss, 2026-11-08. |
| OPS.122 | Job model: a job record with declared dependencies, subjob spawning and a per-job memory claim, replacing the stage barrier |  | Boss, 2026-11-08 (proposed). |
| OPS.123 | Stage opening (the cut), conclude and finish hooks become single-purpose commands of their own instead of running inside the last unit's horse |  | Boss, 2026-11-08 (proposed). |
| OPS.124 | Memory broker: a per-job memory claim, measured and adjusted by load (PERF.254 step 1 becomes its first form) |  | Boss, 2026-11-08 (proposed). |
| OPS.125 | DROPPED (revive only if Boss asks): optional Dask distributed engine behind a switch for the generation stages (PERF.254 step 2) |  | Boss, 2026-11-08 (proposed). |
| OPS.126 | Sector generation as single-sector jobs fed by a spawner |  | Boss, 2026-11-08 (proposed). |
| OPS.127 | Migrate the run stage by stage to the job model, each stage compared byte for byte with today |  | Boss, 2026-11-08 (proposed). |
| OPS.128 | Run planner: every task is planned first as a stored sequence of single-purpose commands for every start path, then executed (takes over the job spec of PERF.203) |  | Boss, 2026-11-08 (proposed). |
| OPS.129 | Split the bundled commands into single-purpose commands: the phenomena finish (about five), run-population (about four) and every other prepare, expand, finish and conclude hook |  | Boss, 2026-11-08 (proposed). |
| OPS.131 | Independence harness for overlapping jobs: each command declares its reads and writes, a SQL tap test fails any command that touches an undeclared table, and a pair overlaps only after a test proves it |  | Boss, 2026-11-08 (proposed). |
| OPS.132 | Overlap rollout: per-pair overlap PRs (first the commands that need only the skeleton: sectors, mass scatter, phenomena) and narrower edges (link and settle per sector, backfill per cell) |  | Boss, 2026-11-08 (proposed). |

### Nebula data foundations

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.202 | Scattered nebulae (planetary nebulae and remnants) join the nebula index, deduplicated by object id |  | Foundations lane 1. |

### Composition research and design

| ID | Item | Needs | Note |
|---|---|---|---|

### Admin control suite

Boss (6:01 PM EDT, 2026-10-21): the admin pages, the Generate suite and the settings are one control suite, built in this order: Settings first, then the shell, job strip and palette, then Generate. Boss answered all 13 decision cards between 6:05 and 6:09 PM EDT (design note docs/design/interface/admin-control-suite.md, section 7); his answer 1 puts all of it in 8.1 phase 1. Owner of every item: the admin redesign thread.

Shape dials (Boss, 6:41 to 6:46 PM EDT, 2026-10-21; goal issue #2355): the generator's shape inputs become five dials built from one lever table (GEN.601), so the Generate suite is not built on the retiring absolute inputs. The table, flags (GEN.602), Scale (GEN.603), Structure (GEN.604), Edge (GEN.605, both forms with a unit selector: L, pc, kpc, Mpc), Halo (GEN.606), the unfilled Shape selector (GEN.608), range sources (GEN.609), the nuclear-amplitude fix (GEN.610), the pattern dial (GEN.611), the dwarf-mass law extension (GEN.612, needs Boss's sources) and the stored pins (DB.81) sit here ahead of the Generate rows; the dials panel ADM.124 follows ADM.110; Census (GEN.607) waits on GEN.545 in the galaxy laws section. The default galaxy is also smaller (Boss, 6:57 PM EDT: still a spiral, no more than half the Milky Way's size; GEN.603 states the number): a new generator epoch, with lower system counts and cost. The Shape dropdown is filled in 9.4 by UX.134 and GEN.389. The lever-table items belong to the Galaxy dials lane (the coordinator, 6:55 PM EDT); Census has no owner until GEN.545 merges (the first release has four dials; Boss, 7:01 PM EDT).

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.612 | Extend the halo, stellar-to-halo and size laws down to dwarf and ultra-faint masses (below 1e8 solar masses), so the shape dials' minimum mass is reachable |  | Shape dials (Boss 7:13 PM EDT 2026-10-21); needs Boss sources. |
| GEN.618 | The stellar-to-halo inversion gives a 6e13 solar-mass halo for 2e11 solar masses of stars, the shape dials' mass ceiling: check the Moster inversion at the top end |  | Density limits report (GEN.613). |

## Phase 2: Infinite zoom

The fly-through Galaxy Map (MAP.146) and the nebula data its tiles carry.

### Fly-through map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.75 | The mini map as a second engine view |  | Locked second camera on MAP.64. |
| MAP.140 | Double-click on a selected object goes there and opens its information |  | GitHub issue [#714](https://github.com/dwhagar/planetGen/issues/714). |
| MAP.257 | Infinite zoom is still clunky: the zoom zones and segments other than sectors get in the way (bug) |  |  |
| PERF.200 | Galaxy Map node cut bookkeeping: about half the busy time in the profile is string building, geometry recomputation and cache re-insertion |  | optional until a lag sighting; Lane 12 |
| TEST.135 | Big-galaxy query-budget test for the nebula reads | GEN.202 | Foundations lane 1. |
| MAP.160 | Quantise the Galaxy Map GPU buffers (deferred) |  | Deferred; MAP.147 recommendation step 3. |
| MAP.189 | The Galaxy Map opens flown in to the created sectors |  |  |
| DB.109 | Optional parity follow-ups: refresh the parity file after a generation run, and keep a before-image on the edit paths so an edit is a parity delta, not a group recompute |  |  |


### Golden galaxy

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.182 | Put one block-first backfill into the golden galaxy |  | Split from GEN.42. |
| GEN.42 | A block-first draw for the backfill (umbrella) | TEST.182 | Same function as PERF.18 (backfill_bright_stars_around); one thread. |

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
| GEN.313 | Orbital update of built systems: advance_orbital_phases for nodes and the GEN.109 hooks |  | Split from GEN.129. |
| GEN.129 | Multi-star systems of up to twelve stars | GEN.312, GEN.313 | GitHub issue [#777](https://github.com/dwhagar/planetGen/issues/777). Build after the GEN.128 design. |
| GEN.593 | Pulsar planets: the disc range stops short of PSR B1257+12's planets, and two of the three formation paths (a captured giant, an ablated companion) leave a pulsar with no planets |  |  |
| GEN.263 | Promote every bound compact pair to the scatter, whatever its mass |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 4. Core. |
| GEN.315 | A fill builds a system from its stored members and nodes; the system pass keeps only planets, moons, belts and disks |  | Split from GEN.264. |
| ADM.128 | The admin one-off system page and the API schema offer the star count and the cap (stars and max_stars) |  |  |
| GEN.316 | Stages D and G draw a member tree for every object the fill creates, and the two-star classes are retired | GEN.315, GEN.332, GEN.640 | Split from GEN.264. |
| GEN.264 | The fill builds a system from its stored members; the system pass keeps only planets, moons and belts | GEN.315, GEN.316 | Multi-star systems at scatter time, Boss 2026-10-15. Step 5. Core. |
| GEN.669 | Residuals of the generator laws after GEN.450 and GEN.270: the star-formation profile shares and the Algol rate against their observed windows |  | Lane 14, 10:20 PM EDT 2026-10-23; no epoch until a law changes. |
| GEN.708 | White dwarf merger products draw their 0.8 to 1.1 solar mass bump: the star model gives the white dwarf its mass before the member tree, so the tree only flags 10 percent as merger products |  |  |
| GEN.707 | Register the white dwarf binary stand-ins in the gap plan registry: kind, evidence label, knobs, validity range and falsifier test for each |  |  |
| GEN.331 | Substellar companions: brown-dwarf close (the desert) and wide companions, brown-dwarf pairs and planetary-mass companions |  | Split from GEN.275. |
| GEN.332 | Circumbinary and S-type planets, circumbinary disks and rogue planet pairs at fill | GEN.315 | Split from GEN.275. |
| GEN.640 | GEN.316 part (c): delete BinaryStarProxy, WideBinaryPair and the writers of the binary_* columns, once circumbinary and S-type planets come from the tree | GEN.332 | Split from GEN.631. |

### Nebula and young stars

| ID | Item | Needs | Note |
|---|---|---|---|

### Region query

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.376 | Region test uses the haversine form for angles under 1e-3 rad in `galaxy/region.py` |  | GEN.372 audit. |

### Superbubble density design

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.391 | Superbubbles shape the stars: the design of the bubble field as a density modulation, with its bounds and touch points checked |  | Boss 9:27 PM EDT: step 1 of the superbubble density chain; Lane 13 proposed. |

### Star composition

| ID | Item | Needs | Note |
|---|---|---|---|

## Phase 4: New scatters: components, clusters and audits

Satellites, clusters and wanderers as components, the globular cluster chain and the audits that check what the scatters make.

### Components and wanderers

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.534 | A satellite's young stars reach the system-type tables: the young-star system types name `satellite_young` |  | Residual of GEN.528: the system-type tables do not name satellite_young. Lane 13 proposed. |

### Globular clusters

| ID | Item | Needs | Note |
|---|---|---|---|

### Whole-sphere audits

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.162 | Multiplicity audit: companion states, row bookkeeping, stability, determinism, no planets in a scatter, fill reuse, motion and the census anchors | GEN.315 | Multi-star systems at scatter time, Boss 2026-10-15. Step 9. Quality. |
| TEST.164 | Catalogue type audit: type counts, region gates, the N tail, interaction consistency and group decay |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 20. Quality. |

### Composition proofs

| ID | Item | Needs | Note |
|---|---|---|---|

## Phase 5: New scatters: core, halo, interface update and docs

The realistic galactic core and the halo tier, the interface update for the spherical grid, the headline item and the docs that follow.

### Galactic core and halo

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.377 | Scatter limits by place for the neutron star and black hole scatter below the cut: the core and halo tiers apply to it, not one global mass cut |  | Lane 13 finding after GEN.333. |

### Interface update

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.115 | Interface update: every page and control works with the spherical grid, 96-bit IDs, the new phenomena and the region queries (stage after the sphere and generation stages) |  | Interface; queued after sphere stages (Boss, 12:52 AM EDT) |
| DB.112 | Drop the stored sector quadrant column and the unused quadrant filter code in db/query.py (migration) |  |  |
| TEST.168 | Control sweep: a browser test drives every control on every page against a spherical galaxy with 96-bit IDs | UX.115 | Quality; queued after UX.115 |

### Habitability system

Boss, 5:38 PM EDT 2026-10-23: the habitability system is finished before 8.1 is ready, five lanes now (this overrides the 8.1-only pause for these items). The model and score are built (GEN.83 to GEN.89, GEN.148, GEN.177, GEN.178, GEN.207, GEN.208, GEN.213; GEN.146 Kopparapu zones). Lane A symbols UX.101; lane B lists UX.100; lane C map MAP.143, MAP.219, MAP.220; lane D model GEN.639; lane E help text, browser tests and the radiation link question (GEN.445 stays in 9.0 until Boss decides).

| ID | Item | Needs | Note |
|---|---|---|---|

### Radiation chain for habitability (GEN.445)

Boss, 5:42 PM EDT 2026-10-23: build the entire GEN.445 chain in 8.1; testing mode, not production, so new generator epochs, database versions and API versions are fine. These 20 items moved here from 8.4 and 9.0 in dependency order. The spine (lane I, to be started): GEN.565, GEN.557, GEN.428, GEN.429, GEN.452 (part 2, PR #2613, held until now for the 8.4 epoch; it now rides the spine's first epoch), GEN.411. Side items that feed the spine: GEN.553 (lane A, after UX.101), DB.73 (lane B, after UX.100), GEN.412 (lane D, after GEN.639). Lane F, radiation constants and ozone: GEN.453, GEN.413. Lane G, bands: GEN.417, GEN.414, GEN.443. Lane H, surface dose: GEN.459, GEN.423, GEN.461, GEN.441, GEN.444. Lane E: GEN.445 with the help text and tests.

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.553 | One star-like body draw: origin, birth, state and [Fe/H], with the initial-final relation, the remnant mass, radius and age rules in one place |  |  |
| GEN.565 | The mass object core and its three family classes: StellarObject, PlanetaryObject and SmallBody, with state and parent as data, replace the nine kind classes |  |  |
| GEN.444 | Surface radiation: the ground-level dose from the air's shielding and the ground's own radiogenic content |  | Boss 11:58 PM EDT; replaces the proxies of physics/radiation.py; generator epoch bump. |

### Docs

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.50 | Docs: the schema, the multi-star note and the catalogue say where the members live | GEN.264, TEST.164 | Multi-star systems at scatter time, Boss 2026-10-15. Step 10. Quality. |
| DOC.18 | Last pass: re-run the documentation sweep when the other DOC items are done and Phase 1 is complete | DOC.50 |  |

### Globular clusters

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.161 | Bright-first fill for cluster sectors |  | Globular-cluster chain. |
| GEN.326 | Blue stragglers in clusters: N_BSS scales as the core mass to the 0.4 | GEN.161 | Split from GEN.162. |
| GEN.163 | Type-B pulsar planets in globular clusters (GEN.130 follow-on) |  | Globular-cluster chain. |
| GEN.330 | Cluster-hosted compact binaries: the cluster rates replace the field fractions | GEN.163 | Split from GEN.274. |
| GEN.274 | Cluster-hosted compact binaries and runaway bookkeeping: `runaway_origin` and the compact companion left by a supernova | GEN.330 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 18. Core. |
| GEN.297 | Calibrate the cluster generator against the Harris file: its distributions guide the laws, and no real cluster is placed in any galaxy |  | The file arrived (Boss, 3:04 PM EDT); Lane 13, cluster line. |

### Galaxy types groundwork

| ID | Item | Needs | Note |
|---|---|---|---|

### Superbubble density

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.392 | A bubble field handle on the galaxy shape: seed-carried bubbles with a box cut, a certified maximum factor and a bubble-free density read | GEN.391 | Step 2 of the superbubble density chain; Lane 13 proposed. |
| GEN.393 | The density terms apply the bubble factors at one call site per form, with a golden test that no bubble means identical density | GEN.392 | Step 3 of the superbubble density chain; waits for Lane 13's profile work so the files are not edited twice. |
| GEN.394 | The scatter bounds and the phenomena pass respect the bubble factor: a certified bound in band_majorants and per-position thinning | GEN.393 | Step 4 of the superbubble density chain; Lane 13 proposed. |
| GEN.395 | Sector build follows the bubbles: the fill count from the cell mean, star positions weighted by the factor, shown in the profile and the region query | GEN.393 | Step 5 of the superbubble density chain; Lane 13 proposed. |
| TEST.206 | Superbubble density tests: rim excess, interior depletion, old stars unchanged, scatter against fill and determinism | GEN.394, GEN.395 | Step 6 of the superbubble density chain. |

### Galaxy laws build

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.228 | test_stars_and_phenomena_fit_within_their_sectors_real_cells compares a dark cloud's centre with its sector's cell, but the cloud is filed by its centroid |  |  |
| GEN.405 | Metallicity and rates by position from the galaxy's own model: the mass-metallicity relation and the local rates replace the fixed gradient and the quoted-at-the-Sun values |  | GEN.368 build, Boss approved the table 10:11 PM EDT; Lane 11 proposed for GEN.400 to GEN.403, Lane 13 or Lane 14 for GEN.404 and GEN.405 (profile.py, activity.py), any lane for GEN.406. Lane 14 proposed (coordinator 10:16 PM EDT). |
| GEN.531 | Boss's research requests from the galaxy laws research (R1 to R15): the observed stellar halo, the dwarf range, the thick disk fits and the cluster system |  | Research requests R1 to R15 of the GEN.406 note; does not gate GEN.301. |

## Phase 6: Bulk composition and the albedo records

The composition fractions that drive density, mass and gravity, stored and described, and the recording of the albedo and greenhouse values. The research, the designs and the data tables they read (the gas, materials, element and stellar tables) are in 8.1 phases 1 and 2 (Boss, 2026-10-16, 11:58 PM EDT: the building blocks of the composition headline go down to 8.1).

Moved into 8.1 from 8.4 phase 1 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Bulk composition

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.75 | Object IDs stay stable: a refill that changes a system's planet or moon count keeps every surviving body's ID, and a body that changes parent keeps its ID |  |  |
| DB.74 | One generic insert and load driven by the column list, one registry of kind tables, and the shared core table with three profile tables | GEN.565 |  |
| GEN.570 | Field declaration: every unit-valued field of the mass object classes declares its kind and native unit, and the one field table serves the serializer | GEN.565 |  |
| GEN.575 | Sweep the mass object classes: the 35 unit-suffixed attributes and every quantity held under two units get one declared native unit | GEN.565, GEN.570 |  |
| GEN.584 | The core's kinematic record: velocity, spin axis and rate, mass and epoch in the point, and the derived columns dropped | GEN.578, GEN.565 |  |
| GEN.583 | reparent: the state-preserving change of a body's parent, shared by re-primary at an event, an ejection, a capture and a merger | GEN.578 |  |
| GEN.582 | Origin and present-state split: fields marked origin or state in the Field table, the draw record held until validation finishes, and the zone and class re-run from the seed | GEN.565, GEN.570 |  |
| GEN.581 | primary_first(members): before a system is validated the heaviest branch becomes the primary and the others are recalculated, one function for the tree, the fill, moons and free pairs | GEN.315 |  |
| DOC.71 | Docs for the System rule: database-schema.md, the system view in docs/design and the header of system_tree.py | GEN.576, DB.78, GEN.578 |  |
| MAP.253 | Map and scene read the parent chain: systemscene, body_positions, the browser twin orbitpositions.js and the systemmap branches drop the three origin conventions | GEN.576, DB.78, GEN.577 |  |
| GEN.578 | The point class gets the local frame and the barycentre system frame, and the system anchor is no longer the nearest star | GEN.576 |  |
| GEN.577 | Stored planet and moon states become the state of the barycentre of the body and its own members | GEN.576 |  |
| DB.78 | Central-object offsets replace the wobble: the reflex_offset columns, calculate_reflex_offset, the SQL mirror and the nine binary position columns go | GEN.576 |  |
| GEN.576 | system(host): the one function that groups bodies by their parent link and returns the barycentre sums, the central object's offset and the valid-system check | GEN.565, GEN.315 |  |
| GEN.587 | Replace the hand-written relation sites with the relations module, by file group, bit-identical on the golden set |  |  |
| GEN.571 | Replace the ad-hoc unit conversions with the registry, by file group, bit-identical on the golden set |  |  |
| GEN.555 | The pair builder for any masses: planet pairs, brown dwarf pairs and a planet with its moons share one stability test | GEN.315, GEN.331, GEN.332, GEN.581 |  |
| GEN.554 | Remnant counts follow from the mass function and the star formation history, with a kick step, and Boss's research densities become test targets | GEN.553 |  |
| GEN.550 | Umbrella (placeholder): one recursive generation process for anything that can have something else orbiting it, from a planet with moons up to a group of star-like objects that each have a system |  |  |

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
| GEN.559 | Unbound comets: a size law and a cut, composition from the comet table, derived activity and velocity |  |  |
| GEN.558 | Moons of wanderers: ejection survivors and isolated-disk moons are table rows, and the fixed 20% boolean is gone | GEN.556 |  |
| GEN.556 | `populate(host)`: planets, belts, comets and moons come from the class tables by one recursive routine, and the system pass is absorbed | GEN.315, GEN.555, GEN.565, GEN.581 |  |
| GEN.455 | Airless bodies and comets draw an exosphere or coma from their surface composition, temperature and the star's wind: the sodium and potassium tails and the production rates |  | Atmosphere research lane new work 6 (parts 04 section 5 and 05); generator epoch bump. Lane 14 proposed. |

### Mixture and stored values

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.594 | Each AD dataset (AD1 to AD8) ships its best possible fallback, used only when a run is told explicitly to disregard missing files |  |  |

### The Quantity call-site sweep

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.349 | API output and stored-column readers take Quantity; the grep test fails on a hand-written conversion factor; the old helpers are removed |  | Split from GEN.300. Moved from 8.3. |

## Phase 8: Ozone, haze, dust and cloud, and the stored shares

The ozone figures and the scatterers that the colour needs, each its own thread, and the columns for the surface and small-body shares.

Moved into 8.1 from 8.4 phase 3 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Ozone

| ID | Item | Needs | Note |
|---|---|---|---|

### Haze, dust and cloud

| ID | Item | Needs | Note |
|---|---|---|---|

### Stored surface composition

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.477 | The asteroid belt and field composition is drawn from the star, the distance and the belt's mass, not four names drawn uniformly from a list of 84 (bug) |  | Body composition research lane finding F11 (tables E1 to E4): belt and field composition drawn from the star, distance and belt mass; generator epoch bump. Lane 14 proposed. |
| GEN.471 | Clathrate hydrates and eutectic brines for icy crusts and ocean moons: the dissociation pressures and the liquid range below 273 K |  | Body composition research lane part F item F4 (later); generator epoch bump. Lane 14 proposed. |

## Phase 9: Stored scatterers, the text and the proofs

The columns for ozone and the aerosols, the atmosphere text written from the mixture, the consistency sweeps and the docs.

Moved into 8.1 from 8.4 phase 4 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Consistency

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.688 | Generation requests name conditions, not class letters: slots and the admin class change take conditions, +habitable_world steers on the PHI-6, the directive counts is_habitable | GEN.684 | Class design step C13. |
| GEN.684 | Albedo and greenhouse come from the composition; delete the two class ranges |  | Class design step C9. |
| TEST.314 | Only classify writes the class column, and the class is recomputed when an input is edited |  | Class design step C5. |
| UX.150 | Pages, lists, search, filters, the API and the map markers show and filter by the derived class and the PHI-6 colour; the class wording and the descriptions follow the rule table |  |  |
| DOC.76 | Docs for the condition-based classes: the rule table, the planet class pages and the help |  |  |
| TEST.240 | Round-trip and registry tests for the mass object core: every kind's stored row round-trips, and no kind table is named in a list the registry should own | GEN.565, DB.74 |  |
| PERF.136 | Measure the fill with wanderers and their moons, and decide whether the rogue cost is on by default | GEN.558, GEN.559 |  |
| TEST.239 | Tests for the unified chain: determinism by host, representation counts per kind and level, rogue and comet rates, derived remnant shares | GEN.554, GEN.558, GEN.559 |  |
| UX.138 | Pages, symbols and descriptions for the new rows: moons of wanderers, free brown dwarfs as stars, unbound comets | GEN.558, GEN.559 |  |

## Phase 10: The colour model

The physics of the colour: the blackbody colour of the star, the Rayleigh optical depth of the mixture, the colour of a solid world's air and of a giant's cloud top, all as the absolute colour (not adapted to the star).

Moved into 8.1 from 9.0 phase 1 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Atmosphere colour

| ID | Item | Needs | Note |
|---|---|---|---|

## Phase 11: Composition and colour in the product

The colour shown in the descriptions, on the maps in place of the planet colour, the composition and radiation text of stars, surfaces and small bodies, and the documents of the colour.

Moved into 8.1 from 9.0 phase 3 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Descriptions and maps

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.446 | The descriptions state the composition of stars, surfaces and small bodies, and the radiation at the surface | GEN.444 | Boss 11:58 PM EDT; stars, surfaces, small bodies and radiation in the text. |

### Stand-in replacements

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.473 | Replace the asteroid taxon shares and the zone kernel, measured in the Solar System belt, by a derived or general distribution |  | Boss 1:43 AM EDT, Stand-in only: replace the measured Solar System taxon shares by a derived or general distribution. |
| GEN.474 | Replace the Earth rock-type library, measured on Earth, by ground types derived from the bulk composition and the melting of the mantle |  | Boss 1:43 AM EDT, Stand-in only: replace the Earth rock-type library by derived ground types. |
| GEN.1021 | Replace the 67P calibration of the comet species factors with a second comet or a laboratory trapping table |  |  |
| GEN.1022 | One latent ice formation temperature per comet drives its CO/H2O, N2 and Ar together (option B of the comet species model) | GEN.1021 |  |
| GEN.1227 | Store the coma class and the dominant volatiles of an interstellar comet, so the map can colour its coma | GEN.1021 |  |
| GEN.475 | Meteorite groups become outputs and not inputs: regions of composition space reachable from the drawn bulk, with the measured groups kept as calibration | GEN.473, GEN.474 | Boss 1:43 AM EDT, Stand-in only: R3, meteorite groups as outputs of the drawn chemistry. |
| GEN.1364 | Core composition split: Fe-Ni, S, Si, O, C and H with a redox and pressure tie, eutectic melting, inner-core growth and the compositional dynamo input |  | |
| GEN.1369 | Vocabulary reconciliation: each of the 84 asteroid component names is drawable (a speciation step names its host mineral) or tagged vocabulary-only, and a test enforces it |  | |
| GEN.1370 | White dwarf core and pollutant: C/O ratio and 22Ne by progenitor, and the polluted white dwarf's element ratios from the belt body that fell in |  | |
| GEN.1371 | Heavy elements: the s- and r-process split from [Ba/Fe] and [Eu/Fe], REE, PGE, Au and Pb, and the highly siderophile partitioning between core and mantle | GEN.1364 | |
| GEN.1373 | Ice giant materials: rock-heavy mixed interiors, diamond precipitation, MgO and SiO2 and water miscibility, ammonia in the superionic phase |  | |
| GEN.1425 | Venus chemistry set V1 to V10 as registered stand-ins and the chain gates that add them to any planet when the physics supports it |  | Boss, 2026-11-08. |
| UX.295 | The class page Conditions line prints raw pascals and kelvins and ignores the visitor's unit setting |  | Lane 11, 2026-11-08. |
| UX.297 | A groundless body (sub-Neptune or gas dwarf) shows its pressure as a reference level, not a surface pressure, on the sheet, the page and the docs | none | Lane 14, 2026-11-11. |
| UX.298 | The single-system generator page cannot ask for a rogue planet or a starless result | UX.283 | Lane 06 first pass, 2026-11-11. |
| UX.299 | The single-system generator page has no weirdness dial or oddities control | UX.283 | Lane 06 first pass, 2026-11-11. |
| UX.300 | The single-system generator page has no debug-areas choice for the debug log | UX.283 | Lane 06 first pass, 2026-11-11. |
| UX.301 | Phone-width and accessibility check of the single-system generator page and its result document | UX.283 | Lane 06 first pass, 2026-11-11. |
| UX.302 | The single-system generator page has no life-class, class-letter or fudge override on its form | UX.283 | Lane 06 first pass, 2026-11-11. |
| UX.303 | Dose rates on the sheets use a unit ladder (nSv/yr to MSv/yr) instead of scientific notation, so "1.32 x 10^4 mSv/yr" reads "13.2 Sv/yr" | none | Boss live report, 2026-11-11. |
| UX.304 | Other rate rows still go scientific or lack a ladder: magnetic field in gauss (uT/mT), flare rates (1/yr), accretion rate (Msun/yr), 1/Gy rates (.4g) | UX.303 | UX.303 builder, 2026-11-11. |
| GEN.1443 | Every bright source in the sky and in the day and night text: second suns, white dwarfs, neutron stars and accretion disks as seen from a body |  | Boss, 2026-11-08. |
| GEN.1447 | Kozai and secular pumping between all bodies: a companion star, sibling planets and the host drive eccentricity and inclination of every planet and moon over the age of the system | none | Boss, 2026-11-08. |
| GEN.1449 | Moon systems under every perturber: the Hill limit, stability, resonance and spin of moons count the star, the sibling planets, a companion and the other moons | none | Boss, 2026-11-08. |
| GEN.1454 | Spin and obliquity under every body: the torques of sibling planets, a companion star and the moons on a planet's spin axis over time | none | Boss, 2026-11-08. |
| GEN.1464 | Belts, asteroids and comets as tiny gravity sources: counted in the secular and tide terms of the planets, and a large asteroid or dwarf planet raising tides on its own moons |  | Boss, 2026-11-08, approved for Train 21. |
| GEN.1466 | Galactic tide and passing stars on the outer edge of a system: the local galactic field sets the Oort-type cloud, the outer limit and the stability of very wide companions |  | Boss, 2026-11-08, approved for Train 21. |
| GEN.1467 | Mutual tides, spin and apsidal precession between the stars of a binary or triple |  | Boss, 2026-11-08, approved for Train 21. |
| GEN.1468 | An N-body spot-check of the stability screens: a short integration of a sample of generated systems against the Hill spacing, Kozai and Laplace-Lagrange screens (report and test only) |  | Boss, 2026-11-08, approved for Train 21. |
| GEN.1485 | Cost the sparse mass rings by block in sky_tiers.build, so a ring that holds few stars is not walked cell by cell | none | Lane 13 read for GEN.1481, 2026-11-09. |
| GEN.1486 | Finer ring edges for the sky backfill: measure, then try narrower rings so fewer stars are drawn than shown | GEN.1485 | Lane 13 read for GEN.1481, 2026-11-09. |
| GEN.1492 | A companion star's secular excitation of a belt: forced eccentricity and inclination of belt members from the other stars of the system |  | Lane 22 finding in GEN.1451, 2026-11-09. |
| GEN.1493 | The air reads the dust and comet influx: late accretion from the belts and giants of the system feeds the planet's atmosphere, which needs the planet and belt generation order changed |  | Lane 22 finding in GEN.1453, 2026-11-09. |
| GEN.1498 | The refill-effect-ledger job: fill ledger rows for systems drawn before the ledger existed, changing no result | GEN.1497 | Split from GEN.1465 at Lane 12's request, 2026-11-09. |
| GEN.1499 | Planets forced by a companion star: the companion's stellar term in the planets' secular exchange | GEN.1492 | Split from GEN.1492 at Lane 22's request, 2026-11-09. |
| GEN.1500 | Mutual inclination of the orbits of the stars in a binary or triple | none | Split from GEN.1492 at Lane 22's request, 2026-11-09. |
| GEN.1501 | Dark members of the orbit tree give light: neutron star cooling, black hole accretion disk, brown dwarf cooling | GEN.1443 | Split from GEN.1443 at Lane 11's request, 2026-11-09. |
| GEN.1514 | Pebble flux decays smoothly after the first million years instead of a step at the pebble end (GEN.1473 stand-in DE3) | GEN.1473 | Lane 14 follow-up, 2026-11-11. |
| GEN.1515 | Pebble lag window scales with the disk's drift speed instead of a fixed 0.05 Myr per AU (GEN.1473 stand-in DE5) | GEN.1473 | Lane 14 follow-up, 2026-11-11. |
| GEN.1516 | Nitrogen, sulfur and noble gases follow the drifted ice weights in a giant's envelope enrichments (N2, NH3, H2S, Ar, Kr, Xe volatility lines) | GEN.1473 | Lane 14 follow-up, 2026-11-11. |
| GEN.1517 | (bug) The water drift model starts the gas inside each ice line as pure evaporated ice | GEN.1473 | Lane 14, 2026-11-11. |
| GEN.1518 | Giant gas-gathering age from a pebble accretion efficiency (stand-in DE9) | GEN.1517 | Lane 14, 2026-11-11. |
| GEN.1519 | Provenance of the planetesimals a giant accretes (where they formed), so noble gases reach the envelope through solids that froze them | GEN.1514 | Lane 14, 2026-11-11. |
| GEN.1520 | The disk outer edge: research the observed gas and dust disk sizes by stellar mass and replace the 18 snow line multiplier (tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER) | GEN.1514 | Lane 14, 2026-11-11. |
| GEN.1521 | Chain-aware spacing test for tight chains of small planets (N-body screen follow-up) | GEN.1468 | Lane 04 report, 2026-11-11. |
| GEN.1523 | Wire the four built-but-unread physics modules: atmosphere_retention, helium_fractionation, hycean_candidate and small_body_composition.draw_ring | none | Wiring guard after retiring GEN.1047 and GEN.1055, 2026-11-11. |
| GEN.1525 | Replace the shortlist constants that held data can check or derive: MIST ages, neutron star EOS maxima, Harris clusters, the mass function integral, the occurrence fit, pre-main-sequence tracks | GEN.1524 | Lane 21 GEN.945 audit, 2026-11-11. |
| GEN.1526 | One definition for repeated physical constants: ICE_DENSITY_KG_M3, BAND_COUNT, AVOGADRO, GAS_CONSTANT, LIQUID_WATER_DENSITY_KG_M3, WATER_MELTING_K | none | Lane 21 GEN.945 audit, 2026-11-11. |
| GEN.1527 | Wire helium fractionation (W16): the envelope-loss history reaches the envelope mixture's helium budget at the capture step | none | Lane 21 split from GEN.1523, 2026-11-11. |
| GEN.1528 | Effect-ledger storage of dark sources: store the member number and kind of a dark source (brown dwarf, neutron star, black hole) so a planet's ledger share for light includes them | GEN.1501 | Lane 11 follow-up, 2026-11-11. |
| GEN.1529 | A standalone neutron star reads its surface temperature from dark_light.surface_temperature_k at the age since collapse, so a neutron star in a tree and a lone one agree | GEN.1501 | Lane 11 follow-up, 2026-11-11. |
| GEN.1530 | Stripped helium star light: a luminosity from its mass (stellar-state tables) and an XUV share, in dark_light.source_of (stand-in row DL7 notes the gap) | GEN.1501 | Lane 11 follow-up, 2026-11-11. |
| GEN.1531 | Roche-lobe overflow feeds a black hole: a companion that fills its lobe makes the hole a bright X-ray source (transfer rate from the companion's evolution state) | GEN.1501 | Lane 11 follow-up, 2026-11-11. |
| GEN.1533 | Give the GEN.1497 body history table a reader: nothing in production code reads the stored history yet | GEN.1497 | Lane 14 sweep, 2026-11-11. |
| GEN.1534 | A Gaia nearby-stars sample fetched at setup replaces the stand-in separation distribution SEP1 (HYG is incomplete and mixes wide pairs) | none | Lane 21 GEN.946 triage, 2026-11-11. |
| GEN.1535 | Obliquity variation, part 3: wire the companion star's mode, store the range (two columns) for the seasonal swing, the Cassini state from the real precession constant, an oddity name for chaotic obliquity | GEN.1454 | Lane 21, 2026-11-11. |
| GEN.1536 | The Cassini state of a planet's tilt from its real precession constant (needs the stored J2 or C22 figure), replacing the log-uniform stand-in | GEN.1535 | Lane 21, 2026-11-11. |
| GEN.1537 | Class and life census re-tag after trains 22 and 23 (epochs 144 to 146): re-measure the census shares and re-baseline the class and life ranges the trains moved | GEN.1501 | Lane 11, 2026-11-11. |
| GEN.1538 | Design check: settle the orbit before making the air (birth orbit, settle, air once) instead of making the air and remaking it after the settle | none | Boss, 2026-11-11. |
| GEN.1539 | Scale the outer-disk belt start (tuning.PLANET_OUTER_GROWTH_AU, 11.36 au for every star) with the star's mass or disk size | GEN.1520 | Lane 14, 2026-11-11. |
| GEN.1540 | (bug) A rare captured moon's position radius is outside its allowed range (a wrong or stale frame), so test_orbital_motion fails by draw luck | none | Lane 23, 2026-11-11. |
| GEN.1541 | The tide a large asteroid or dwarf planet raises on its own moons (GEN.1464 part C) | GEN.1464 | Lane 23, 2026-11-11. |
| GEN.1542 | Per-planet dust supply for the Venus sulfur cloud set (the V8 species), so a planet's dust decides whether the V8 layer exists | GEN.1425 | Lane 14, 2026-11-11. |
| GEN.1543 | (bug) The dust influx of physics/influx.py treats a system's belt and comet dust as an isotropic stream, so the Earth gets about 45 kg a year against 1e7 measured: fix the geometry and the scale, then rewire the V8 sulfur cloud dust to it | GEN.1542 | Lane 14, 2026-11-11. |
| GEN.1544 | (bug) A planet's cosmogenic isotope rows (Al-26, Be-10) are first made before its magnetic field exists, so they use a field factor of 1.0 until the final-age pass remakes them | none | Train Manager (filed 2026-11-11). |
| GEN.1545 | (bug) A star whose age is raised for its planets is not evolved again, so its light, temperature, radius, type, X-rays, flares and spin stay at the first-guess age | none | Boss, 2026-11-11. |
| GEN.1546 | Forward build, item 1: remove the age floor (age_needs, the age-raising code, the trim and belt rebuild after the bump, its retry reason) so the star is made once at its final age | none | Boss, 2026-11-11. |
| GEN.1547 | Forward build, item 2: key the dice of every end stage to the body (air, albedo, surface water, spin, giant interior), so a remake or a move never re-rolls | GEN.1546 | Boss, 2026-11-11. |
| GEN.1548 | Forward build, item 3: reorder the system build so each thing is made once (star, bare bodies, spacing, trims and belts, orbits settled, air, water and surface on the final orbit, class and life last) | GEN.1546, GEN.1547 | Boss, 2026-11-11. |
| GEN.1549 | Census re-tag after the forward build (items 1 to 3): re-measure the class, life, star and age census shares and re-baseline the ranges the new order moved | GEN.1548 | Boss, 2026-11-11. |
| GEN.1550 | (bug) StarSystem._set_debris_disk clears a young star's class 0, I or II disk before its class check, so a young star that gets a belt loses its disk | none | Boss, 2026-11-11, young-stars research. |
| GEN.1551 | Each star gets its own disk dispersal time, drawn once from its own stream: a star holds a class 0, I or II disk only while younger than that time, and a close companion shortens it (GEN.1475 reads this time) | GEN.1546, GEN.1550 | Boss, 2026-11-11, young-stars research. |
| GEN.1552 | Bodies appear only once the star's own disk is gone: the fixed 10 Myr planet cut becomes the star's own dispersal time, with belts as the debris that follows | GEN.1551, GEN.1546 | Boss, 2026-11-11, young-stars research. |
| GEN.1553 | Census (add to GEN.1549): share of stars with a disk, with bodies, with neither, by age bin 0 to 1, 1 to 3, 3 to 6, 6 to 10 and 10 to 30 Myr, against the survey curve | GEN.1549, GEN.1551, GEN.1552 | Boss, 2026-11-11, young-stars research. |
| GEN.1554 | K and M dwarf zone edges against the greenhouse spectrum: K5V worlds are cold in the inner third of the zone | none | Boss, 2026-11-11, Class M fixes. |
| GEN.1555 | (bug) The surface air solver has two stable answers near the moist limit and water's boiling point, and the fixed point flips or never settles (a hot CO2 air against a cool one, or a steam air against a frozen airless body) | none | Lane 14 and Lane 21 reports, 2026-11-11. |
| GEN.1556 | The Cassini state of a moon's tilt from its real nodal regression (the planet's bulge and the star), part 2 of GEN.1536 | GEN.1536 | Lane 21 report, 2026-11-11. |
| GEN.1557 | (bug) The class T and class I sentences describe a small core and mostly ices, but the stored bodies are rock-dominated (class T rock 50 to 70 percent, envelope 1 to 18 percent), and the day-counter line says "No solid surface" without saying the rock lies deep under the envelope | none | Boss live report, 2026-11-11. |
| GEN.1558 | Retire the stepped pressure caps (recipe row d6: d3 = 1 kPa and d4 = 10 kPa) so supply, sinks and escape set the air pressure, with the built energy-limited loss added to the escape subtraction in the same single pass | GEN.1546, GEN.1547 | Boss, 2026-11-11, class census. |
| GEN.1559 | Rewire burn-down: every physics module, function, registry row, data file, dataset and column that the wiring guard lists as unread and not by design gets wired, retired or marked by design, class-moving ones first | none | Boss, 2026-11-11, wiring audit. |
| GEN.1560 | Earth-mass planets in the temperate zone of FGK stars: the isolation-mass law gives a median 0.3 Earth masses at 0.7 to 2.5 AU, so add the merge of oligarchs (giant-impact stage) or fit the mass function to [R0572] | none | Census thread, 2026-11-11. |
| GEN.1561 | Belt mass scale: belts hold a median 0.21 Earth masses, 25 percent hold over 1, and the largest member exceeds Earth's radius in 57 belts; decide whether planet-sized belts are intended and, if not, merge the tail into planets | GEN.1560 | Census thread, 2026-11-11. |
| GEN.1562 | Sulfur dioxide airs are never made: classes A and Y are unreachable (1,391 of 1,392 airs hold under 3 percent SO2); decide whether volcanic bodies (Io-like, tidally heated, close super-Earths) get SO2 airs and wire the supply, otherwise mark A and Y reserved | none | Census thread, 2026-11-11. |
| GEN.1563 | Ammonia-water fraction never reaches 0.1, so class X is unreachable (157 bodies pass every other condition); find where ammonia_water_mass_fraction is set and why it stays under 0.1 | none | Census thread, 2026-11-11. |
| GEN.1564 | Freeze-out for every gas, sulfur (S2) included: add S2 to the saturation table and run the freeze-out for all gases, including hycean and ice worlds at low temperature (93 gas cases on 61 bodies are supersaturated, 36 bodies hold hot S2 at cold surfaces) | none | Census thread, 2026-11-11. |
| GEN.1565 | Oxygen worlds keep too much carbon dioxide: all 9 temperate oxygen-rich worlds have CO2 over 1 kPa or a sea over 95 percent; reopen the GEN.1309 carbon draw-down question with this measurement | GEN.1403 | Census thread, 2026-11-11. |
| GEN.1566 | Ocean cover is bimodal: 47 of 77 temperate sea worlds are 95 to 100 percent ocean and only 11 are 20 to 80 percent; research the water law (mantle cap, relief, delivered water) in the papers first, then fix the cover distribution (the biggest remaining physics lever after the pressure caps) | none | Census thread, 2026-11-11. |
| GEN.1567 | Split ocean_fraction into liquid cover and a new ice_cover_fraction (land ice and sea ice apart), so the H, M and O gates read liquid cover only and ice worlds stop storing 1.0; land fraction for weathering is 1 minus liquid cover | none | Census thread, 2026-11-11. |
| GEN.1568 | Temperature guards pile bodies up: the 2,000 K cap holds 131 bodies and the 375 K guard 32; list the classes they feed and check class N against the tier T5 note | none | Census thread, 2026-11-11. |
| GEN.1569 | One-off check of the class mix against Kepler: gas dwarfs are 30 percent of bodies and 37 percent of planets, T tripled while V and E fell; compare with the observed planet-radius mix | GEN.1549 | Census thread, 2026-11-11. |
| GEN.1570 | Class edges for H, M and O on liquid cover: O ocean_fraction 0.95 or more (was 0.8), M under 0.95 (was under 0.8), H unchanged (under 20 percent); no body changes letter today | GEN.1567 | Boss, 2026-11-11, water-land research. |
| GEN.1571 | Re-run the ocean cover histogram after the water-law decision and check that O is rarer than M | GEN.1566 | Boss, 2026-11-11, water-land research. |
| GEN.1572 | Optional ice-capped flag when sea ice plus land ice exceeds 25 percent of the surface (not counted in liquid cover) | GEN.1567 | Boss, 2026-11-11, water-land research. |
| GEN.1573 | Water-load basin deepening in hydrosphere.ocean_fraction: ocean-crust depth grows with the water load up to d_o^max = 11.4 km / (g / g_Earth), so a wet planet keeps some land ([R0914: Eq. 10]), with the Earth anchor kept | GEN.1567 | Boss, 2026-11-11, water-law research. |
| GEN.1574 | Mantle water partition as a seafloor-pressure feedback instead of the log-uniform 1 to 8 ocean mantle buffer (MANTLE_CAP_OCEANS), so the mantle takes up less when the seafloor pressure is high | GEN.1566 | Boss, 2026-11-11, water-law research. |
| GEN.1575 | Check the delivered-water high tail (GEN.849) against the formation papers: read [R2381] in full plus [R3193, R3683] for the spread of water delivered to terrestrial planets | none | Boss, 2026-11-11, water-law research. |
| GEN.1576 | (bug) A planet whose thermal and gravitational tides balance away from synchronism (body_spin kind asynchronous) can spin past the breakup limit, because the balance rate is not capped at spin.BREAKUP_FRACTION (0.85) of the breakup rate | none | Lane 12 finding, 2026-11-11. |
| GEN.1577 | Upper-atmosphere temperature-pressure profile of a hot giant (1e-6 to 1e3 bar), so the cloud decks and trace-gas ceilings of GEN.994 (hot_condensates) have a profile to read | GEN.994 | Lane 14 triage finding, 2026-11-11. |
| GEN.1578 | Research the moist limit as a function of the air (replacing the fixed 375 K MOIST_LIMIT_K guard): the greenhouse state at which a wet or dense air runs away, fitted from the greenhouse table | none | Lane 14 guard report, 2026-11-11. |
| GEN.1579 | Replace MOIST_LIMIT_K in surface_air with the moist limit as a function of the air, so 375 K stops pinning dense nitrogen and CO2 airs | GEN.1578 | Lane 14 guard report, 2026-11-11. |
| GEN.1580 | Raise the 2,000 K upper temperature cap to a sourced, very high physical ceiling (the ionization or plasma temperature of the air's main gas at its pressure, tens of thousands of K), after checking nothing assumes 2,000 K | none | Boss, 2026-11-11, upper cap. |
| TEST.714 | Tests that pin a seed for a rare case find the case by searching seeds instead | none | Train Manager finding, Train 22, 2026-11-09. |
| TEST.715 | The constants completeness gate: turn the GEN.945 census into a test with a committed baseline of unregistered constant names that may only shrink | none | Lane 21 GEN.945 audit, 2026-11-11. |
| OPS.133 | scripts/train.py stamps each train epoch with the union of the stages its items declared, plus a fifth stage 'derived' for refills, and takes the epoch out of the map node-cache key (audit C7) |  | Generation rescan audit, 2026-11-09. |
| OPS.139 | Heartbeat part 6: the admin Queue page shows each live unit's last beat age, counter, phase, rate and state badge | OPS.138 | Boss approved the heartbeat, 5:16 PM EST on 2026-11-09. |
| OPS.143 | Heartbeat part 10: tests for every rule with a fake clock, and a full-galaxy trial on the lane machine | OPS.137, OPS.140, OPS.141 | Boss approved the heartbeat, 5:16 PM EST on 2026-11-09. |
| OPS.146 | Heartbeat part 12: a restarted job runs in debug mode with its own debug log (config, on by default) | OPS.140 | Boss, 5:29 PM EST on 2026-11-09. |
| OPS.147 | Per-stage golden digests so a train's declared stages can be checked against what its golden moved | OPS.133 | Lane 05 audit, 2026-11-11. |
| GEN.1459 | Settled eccentricity under every body: the pumping of siblings, moons and a companion against the tidal wearing-down | none | Boss, 2026-11-08. |
| MAP.409 | Page error: wheel into a system throws "Cannot read properties of null (reading 'update')" about 60 times and the scene never fades in |  | Lane 03, 2026-11-08. |
| MAP.410 | Page error: editing the waypoints throws "Cannot read properties of null (reading 'textContent')" |  | Lane 03, 2026-11-08. |
| GEN.1431 | Flavor text life subtree: lines for all nine life classes at rungs 2 to 5, gated on the life's solvent medium, energy source, rung and animal grade (replaces the multicellular-only trigger) | GEN.1430 | Boss, 2026-11-08. |
| GEN.1433 | Flavor text for rogue planets, belts, comets, stars and phenomena: one root of the tree for each scope | GEN.1430 | Boss, 2026-11-08. |
| GEN.1434 | Derived sky and clock fields from physics: gravitational time dilation, the host's sky size, the eclipse calendar and the solar-lens focus distance |  | Boss, 2026-11-08. |
| GEN.1437 | Design notes for the chaotic triple-star zone and the dimming star (debris stream or disintegrating planet with irregular transits) |  | Boss, 2026-11-08. |
| GEN.1438 | Flavor tree tests: classless gates, quantities exist, every line reachable, coverage per life class and scope, order independence, plain words, settings pruning | GEN.1430 | Boss, 2026-11-08. |
| GEN.1382 | Ore processes build, weathering and pegmatite deposits (bauxite and laterite, REE and Li pegmatites) | GEN.1371 | |
| GEN.1383 | Ore processes build, carbonatite, kimberlite and unconformity uranium deposits | GEN.1382 | |
| TEST.686 | Real-material closing gate: a composition census over a generated sample (each new element and mineral appears where physics says and nowhere it should not) and a check that every stand-in the set added has its GEN.725 row | GEN.1383, GEN.1370, GEN.1369 | |
| GEN.1384 | Reduced-body recipes read the sulfide, nitride and carbide rows: the metal ground and the reduced taxa name pyrite, pyrrhotite, pentlandite, chalcopyrite, oldhamite, niningerite, alabandite, daubreelite, osbornite and sinoite instead of troilite and kamacite stand-ins |  |  |
| GEN.1391 | Minor core elements: Cr, Mn, P and Cu partition between the core and the mantle, read from the star element set (GEN.1360) and the core split (GEN.1364) | GEN.1364 | Lane 23 follow-on to GEN.1364. |
| GEN.1392 | Reduced cores and the dynamo energy scale: let rock FeO reach below IW-3.3 so Mercury-like Si-rich cores occur, and re-tune the compositional dynamo energy that stopped most small-body dynamos | GEN.1364 | Lane 23 follow-on to GEN.1364. |
| GEN.1395 | Tune a gas giant's carbon-to-oxygen ratio to observed Sun-like giants and to what physics expects |  | Boss, 9:43 PM EST on 2026-11-07. |

### Weight research (atmosphere)

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.483 | Source the cosmic-shoreline constants and the pressure caps by retention regime (AR4) |  | Atmosphere lane research need R4 (rows d6); replaces the seeded bracket of the DESIGN rows. |
| GEN.492 | Source the ground radioisotope abundances by host and the potassium depletion by body class (AR13) |  | Atmosphere lane research need R13 (rows o2, o3, o4, o6); replaces the seeded bracket of the DESIGN rows. |

### Tuned draws and selection models

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.510 | Tuned draw builder B2: delivery and water worlds |  | Tunable gaps report: Tuned draw builder B2: delivery and water worlds. |
| GEN.515 | Tuned draw builder B7: ground and petrology |  | Tunable gaps report: Tuned draw builder B7: ground and petrology. |
| GEN.519 | Tuned draw builder B11: aerosol amplitude |  | Tunable gaps report: Tuned draw builder B11: aerosol amplitude. |
| GEN.595 | Replace the Solar System comet population shape (GEN.520 D6) with a derived distribution |  | Boss, 2026-10-18 (fork 8, flagged stand-ins only). |
| GEN.522 | Tuned draw builder B14: moons and Kuiper belt objects |  | Tunable gaps report: Tuned draw builder B14: moons and Kuiper belt objects. |
| GEN.523 | Tuned draw builder B15: host-star activity |  | Tunable gaps report: Tuned draw builder B15: host-star activity. |
| GEN.524 | Tuned draw builder B16: halogens and trace gases |  | Tunable gaps report: Tuned draw builder B16: halogens and trace gases. |

### Physics rectification (Boss, 5:30 PM EDT 2026-10-24): belts, fields, irregular moons, overlap-free systems and the oddity table

Boss accepted the physics fix plan (reports/physics-fix-plan-2026-10-24.md) and put the physics fixes into this release: "I don't want to go forward without correcting the physics." These items moved here from 8.4. They are prerequisites of the headline GEN.301.

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1006 | Stars and bodies in one system never overlap by design: every member knows the others when any is rebuilt, then fit and adjust afterwards |  |  |
| GEN.1017 | Overlap-free stars P0: give the star tree the primary's [Fe/H] (and the build's age and mass) so the tree and the build evolve each companion the same way (replaces the dropped GEN.1007) |  |  |
| GEN.1016 | Overlap-free stars P10: per-level hierarchical stability margin and the 2+2 quadruple fit, if the verifier shows it is needed |  |  |
| GEN.1015 | Overlap-free stars P9: delete the GEN.1004 quick fix (built_tree_clears, the clearance branch of _companions_fit, _built_stars_clear), no shim | GEN.1014 |  |
| GEN.1014 | Overlap-free stars P8: a clearance ledger and planet, moon and belt sampling from it |  |  |
| GEN.1105 | Field members (individual bodies by the belt cut rule) with `<field ID>.<rank hex>` designations, status flags, and the sector and field pages | GEN.1103, GEN.1070 |  |
| GEN.1103 | AsteroidField columns, radius derived from age, the field density leaving 0, registry text changed |  |  |
| GEN.1092 | Wire the scatter phase and the belt origin: scatter tag on stored members, the outer disc fed from the scatter step, origin and formation shown on the belt page | GEN.1070 |  |
| GEN.1086 | System map and pages for irregular moons (eccentric, inclined, retrograde orbits on the moon-orbit view, a body page per irregular) and the clearance ledger taking the irregular zone as occupied space |  |  |
| TEST.660 | The moon class ranges are re-baselined with the real body counts after the last 8.1 train |  |  |
| GEN.1073 | Bring-to-epoch for belts: grinding recomputed from the stored age, members above the grinding scale kept | GEN.1069, GEN.1006 |  |
| GEN.1072 | White dwarf and giant-branch belt events wired into remnant_orbits and the engulfment step, replacing the uniform expansion factor of belts | GEN.1068 |  |
| GEN.1070 | Store the individual belt members (rows with parent link and rank), the designation display, the system map drawing them as dots at true orbits, and a body page per member; API version bump | GEN.1069 |  |
| GEN.1509 | Design: belt bodies at or above the smallest planet class floor get full planet status | GEN.1070 | Boss, 6:08 PM EST on 2026-11-09. |
| GEN.1510 | Build: belt bodies at or above the class floor run through the planet generator | GEN.1509 | Boss, 6:08 PM EST on 2026-11-09. |
| GEN.1511 | Promoted belt bodies appear as planets in names, pages, lists, maps and the census | GEN.1510 | Boss, 6:08 PM EST on 2026-11-09. |
| GEN.1069 | Store the belt mass, Dc, age history text and formation snow line (DB.70 column, migration) and wire belt_mass into AsteroidBelt, replacing the three density words | GEN.477, GEN.1068 |  |
| GEN.1068 | Belt birth from the architecture: a belt is a gap or an unfinished tail in the planet chain, not a 10 percent slot roll, and sees the giants (GEN.728 slot loop, system.py) | GEN.1014 |  |
| GEN.1231 | White dwarf pollution flag reads the belt's scattered rate (is_polluted, pollution_rate_g_s) and the census band 25 to 50 percent | GEN.1068 |  |
| GEN.1150 | Lagrange trojans T2: the generator draws a trojan swarm at the L4 and L5 of each giant planet or brown dwarf from the scatter trojan share; new columns co_orbital_host_id, lagrange_point, libration amplitude and phase (db/column_units.py); orbit updater libration overlay; the object and sector pages and the designation rule "<host> L4 trojan region" | GEN.1006 | Trojans T2, epoch lock, train 8 |
| GEN.1151 | Lagrange trojans T3: whole planets and moons (Class M and Class J allowed by the mass rule) at L4 and L5 from the disc formation model only; no forced rate, the observed 5 percent upper bound (WR10) is a plausibility test only; folds in the weird-worlds W7 row | GEN.1150 | Trojans T3, epoch lock, train 8 or 9 |
| MAP.354 | Lagrange trojans T4: System Map icon at the centre of each Lagrange point that has objects, labelled with the host and point (for example Jupiter-L4); shares the geometry of the MAP.169 gravity layer | GEN.1150 | Trojans T4, with T2 |
| TEST.492 | Census re-baseline for every oddity share and the prevalence suite rows (one pass after train 9) | GEN.1047 |  |

### Physics rectification: test health, constant audits and the star mass function

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.395 | A planet distribution test with observed-rate bands, run after the class cut-over on the final classes and on galaxy-placed stars: giants and hot Jupiters per host mass, planets per mass bin, 1 to 20 Earth-mass planets per star, habitable-class bodies per star, and a rule for classes O, P and Q |  |  |
| TEST.417 | Whole-sector and whole-galaxy tests take 20 to 45 minutes: make them smaller or cached so a gate run finishes (bug) |  | found by Aux 8 |
| TEST.455 | A completeness gate: every numeric constant in the generation and physics modules has a stand-in registry row or is read from a table; the list of unregistered constants may only shrink | GEN.945 |  |
| GEN.922 | The star mass function is built right: initial mass function and present-day populations match the observed curve segment by segment, for every mass range and every star-like object type (parent, replaces the one-number GEN.914 tweak) |  |  |

### Physics rectification: life classes, chemistries and their pages

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.878 | Chemistry-specific energy: the pigment's absorption band times the star's spectrum is the multiplier on LIGHT_SHARE, instead of one bolometric curve for every metabolism |  | Life chemistry review; moved to 8.1 by Boss, 8:00 PM EDT 2026-10-28. |
| GEN.1307 | Class 5a life (sulfuric-acid cloud layers) has no surface-temperature window, so a 2595 K vapor world stores animals aloft (bug) |  |  |

### Physics rectification: infrared opacity, tides and the class re-sweep

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1019 | Venus clouds and the CO2 far-wing line shape in the infrared model (the IR3 far wing is a labelled stand-in) |  |  |
| GEN.1144 | Post-train-6 class re-sweep on the new main, then decide the E and F merge, the T versus N boundary and the H and L reach |  |  |

### Physics rectification: names, designations and Object IDs

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1141 | Designation history part e: bodies, belt members and captured bodies follow a parent ID change, and tombstones on delete |  |  |
| API.38 | API and search for local names: name is the local name, a new sector field carries the place, a typed name returns every match with its sector, an ID typed in full still opens its object |  |  |
| UX.259 | Uncharted star and phenomenon pages, lists, infobox and map labels show the local name with its sector, not 'Uncharted G2V star <row id>' |  |  |

### Physics rectification: strength, rings, trojans and small-body follow-ups

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1179 | Optional: Hoek-Brown geological strength index as a body parameter for rock |  |  |
| GEN.1175 | Porosity weakening of body strength from a sourced law instead of the placeholder exponential |  |  |
| GEN.1167 | The stellar-to-halo relation [R3181] has 0.3 to 0.4 dex more stars than [R0319] in halos of 1e10.75 to 1e11.5 and 0.06 to 0.26 dex fewer in halos of 1e12.5 to 1e14.5: compare the published fits and pick the best-supported relation |  |  |
| GEN.1204 | Dark members of an orbit tree (brown dwarf, neutron star, black hole in dark_members) store density, binding energy and the shape values the stars, rogues and comets store |  |  |
| GEN.1159 | Store ice fraction, bulk density, roundness ratio, yield strength and binding energy on planets, moons, asteroids and belt or field members (columns, migration, API version) |  |  |
| GEN.1201 | A per-body shape-setting temperature: the roundness engine uses the temperature a body took its shape at, but no body stores one, so strength is the cold strength today | GEN.1159 |  |
| GEN.1157 | Earth-type moons around gas giants: how often a giant has a moon the size of a planet, and how big, from the physics of moon formation and capture |  |  |
| GEN.1153 | Day and night state for stars and rogue planets (follow-up to GEN.1097, which covers planets and moons): a star stores the constant "no night"; a rogue planet gets the same state, night hours and solar days as an in-system planet from its own spin, with no star the answer is no solar day |  | Aux 2 report 10-30 |
| GEN.1152 | Mass literals that differ from constants.py (changing them changes output): moon_composition JUPITER_MASS_KG 1.89813e27, ring_origin star_mass_kg 1.989e30, giant_interior 317.8 and EARTH_MASS_JUPITER 317.83, stellar_abundances JUPITER_MASS_SOL, tuning 317.8; use the constants.py values with an epoch |  | Sources 1 report 10-30, train-bound |
| GEN.1129 | Fingerprint split: a sector's digest (GEN.58) becomes "made here" and "arrived here" so rogues landing from other sectors do not change the sector's own digest |  |  |

### Physics rectification: greenhouse refits, region runs and tests red on main

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1202 | The generator's disc gas surface density at 8 kpc is 4.3 solar masses per pc2 against 13 measured, which holds the disc cluster maximum mass too low |  |  |
| GEN.1128 | Region run command: a sphere of a given radius from a start point becomes whole sectors, cut into units of a block of sectors by the galaxy's own (shell, band) arcs, nearest bands first, with `--grow-shells K`, all writing the chart ledger |  |  |

### Physics rectification: the class-by-domain matrix and stored liquids and air

| ID | Item | Needs | Note |
|---|---|---|---|

### Physics rectification: citations and the gate's test health

| ID | Item | Needs | Note |
|---|---|---|---|

### Physics rectification: best-fit classes, orbit spacing and rogue wiring

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1314 | Fudge table part 4: class_fit reads fudge_table.json for the edge shoulders, the hand-picked CF and DM rows are deleted, the miss budget and matches for former holes are stored | GEN.1313 |  |
| GEN.1315 | Fudge table part 5: outcomes read the match (habitable count by threshold, the barred-world loop, +habitable_world, search filters) | GEN.1314 |  |
| GEN.1316 | Fudge table part 6: the stand-in ranking report gets a share-of-fudge column; condition-based-planet-classes.md and the /classes pages explain the fudge | GEN.1313 |  |
| GEN.1317 | Fudge table part 7: update checks the data against the epoch and regenerates the fudge table automatically when the epoch has changed | GEN.1313 |  |
| OPS.118 | Find what sent SIGTERM to the work horse on unit 11 of Boss's warm-up (status 15): no code path in the program sends it, so it is outside the program or a manual stop (bug) |  |  |
| DB.115 | A generation job stage fails with InterfaceError (0, '') after a dropped MySQL connection: find out why MySQL dropped the connection (the masking rollback and the sky-tier retry are fixed) (bug) |  |  |

### Physics rectification: winds and ring darkening

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1110 | Average wind speeds by day and by night in each of the seven latitude bands: the day-to-night thermal wind from the pressure difference, zonal jets and trade winds from the equator-to-pole temperature difference and the Rossby number, faded out below about 1 kPa |  |  |

### Physics rectification: stars, spectra and climate follow-ups

| ID | Item | Needs | Note |
|---|---|---|---|

### Physics rectification: older physics follow-ups, units and the storage rows

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.820 | Residuals of the bulk iron branch table (GEN.512): the coreless branch, a slope other than 1, the 0.05 dex scatter stand-in and a larger tuning sample |  |  |
| GEN.835 | The tidal heat flux (tidal_flux_w_m2) is wired into the ground draw's sulfur rule |  |  |
| DB.101 | Storage and units for every physics rectification output: a column, a column_units entry and a migration shipped with its item |  |  |
| GEN.805 | The surface dose and the habitability radiation component read the stored flare irradiation index and magnetopause (G2, step 2) | GEN.804 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
| GEN.944 | Later: an age-metallicity relation inside each population and a 10 to 13 Gy thick-disk window; variable-star flags (Cepheids, Mira and long-period variables, carbon and S stars) |  |  |
| GEN.949 | Wire the held methane absorption tables (AD1 Karkoschka k-tables, PDS, NMSU) into the visible colour of giant planets |  |  |
| GEN.993 | Widen the core light-element total in bulk_draw to 0.02 to 0.25 (0.25 dex) instead of the uniform 0.03 to 0.12 |  |  |
| GEN.994 | Wire the hot-giant regime edges (equilibrium_edges) into the atmosphere row for hot giants |  |  |
| GEN.91 | Classes like S and V in the hot and cold zones |  |  |

### Temperature coverage, nebula composition and belt rules (rows filed in the 8.2 file by mistake, moved here)

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1267 | A generated nebula's own shape follows the same model-density contour the map now draws: bounds, containment and the stars placed inside it |  |  |
| MAP.423 | (bug) Hot-star nebulae (emission, reflection, planetary, remnants) near the Sun and in the outer disk are effectively not drawn on the map | GEN.1267 | Lane 12 census, 2026-11-11. |
| MAP.424 | (bug) The Galaxy Map's sector list query (GET /api/galaxy/sectors) hits the 10 s database statement timeout on Boss's server, so /galaxy returns 504 and the cache warm-up fails with QUERY_TIMEOUT | none | Boss live report, 2026-11-11. |
| MAP.425 | (bug) On the 3D system view the highlighted orbit ring of a moon does not pass through the highlighted moon (the ring is drawn from inclination and node measured in the planet's equator, the moon from a state vector in the star's frame) | none | Boss live report, 2026-11-11. |
| TEST.717 | (bug) test_math_gate and test_datasets_gate fail together on main: the status-command subprocess tests see /var/cache/planetgen/datasets from the machine instead of the test's own data directory | none | Lane 23 report, 2026-11-11. |
| TEST.718 | Release gate: zero unread items that are not marked by design, beside SEC.33 and TEST.261, so a version cannot be released with built physics that nothing reads | GEN.1559 | Boss, 2026-11-11. |
| TEST.719 | Make the class census a repeatable script (scripts/class_census.py from reports/class-census-2026-11-11/collect.py and report_tables.py), with forced-young and caps-off options, for GEN.1549 and every later train | none | Census thread, 2026-11-11. |
| TEST.720 | Pre-gate quick lane (3 to 5 minutes) that the author runs before pushing and the gate runs before the long run: release note, pyflakes, spelling, own-log-format, citations, settings, wiring guard, constants, and merges cleanly with main | none | Gate lane report, 2026-11-11. |
| TEST.721 | Stop committing generated files that conflict, or give them a merge step that regenerates them (docs/generated, stand-in-values.md, overview.json, stand-ins.json, the references files) | none | Gate lane report, 2026-11-11. |
| TEST.722 | Re-run only the failed files after a fix, and pick tests from the changed modules with a stated time limit per PR (the rest left to the end-of-version run) | none | Gate lane report, 2026-11-11. |
| TEST.723 | Fail loudly when database tests are skipped for no MySQL login, and run a short smoke test of main after each merge | none | Gate lane report, 2026-11-11. |
| TEST.724 | Cut gate flakes: each shard machine gets its own Redis database and MySQL prefix, at most 3 test workers with nothing else running, and the rare unseeded failures go on scripts/data/known_reds.json | none | Boss, 2026-11-11, train management. |
| TEST.725 | Output-changing fixes found in a train gate ride the next train (unless they crash, lose data or change a visible page), with a known-issue line in the train PR and in the merge message, a cap near 1 percent of draws, and the fix first in the next train | none | Boss, 2026-11-11, train management. |
| GEN.1272 | Nebula composition part d: the nebula page, the describe sheet and the Properties sheet show the composition with units and derived marks |  |  |
| PERF.240 | A system that bars habitable worlds costs about 12 s to draw against 3 s, because the late remakes and planet redraws that remove life-class bodies run again and again (bug) |  |  |
| TEST.582 | Temperature coverage sweep: set the Blue cut and the Green and Yellow edges on a sweep that includes cold-mean worlds with warm regions and hot-mean worlds with cool regions, then re-baseline the TEST.430 prevalence census |  |  |

## Phase 12: The headline closes

The last audits and the headline: they test the component pass, the core and the halo, so they come after phases 4 and 5.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### The headline closes

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.261 | Release gate: a full bug hunt audit of the whole code base and its behaviour, done before every major or minor version change (8.1 first) |  | Standing release gate, Boss 2026-10-21 5:10 PM EDT; also gates every later major or minor version change, not build numbers (Boss, 8:01 PM EDT 2026-10-22). |
| TEST.380 | First live test pass of species, civilizations, polities and territories on a generated galaxy, before the 8.1 stamp |  |  |
| SEC.33 | Release gate: a full security code audit of the whole code base, done before every major or minor version change (8.1 first) |  | Standing release gate, Boss 2026-10-21 5:10 PM EDT; also gates every later major or minor version change, not build numbers (Boss, 8:01 PM EDT 2026-10-22). |
| GEN.301 | Headline of 8.1: the new generate scatters, accurate galactic generation and scattering including the halo | SEC.33, TEST.261 |  |
