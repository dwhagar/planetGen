# Version 8.1: Infinite zoom and the new generate scatters

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
| PERF.226 | Every generation stage has a measured worker memory, and the worker cap uses it, so no stage past the scatter can run the host out of memory |  | Boss run OOM, hal. |
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
| GEN.1402 | A rocky planet whose orbit eccentricity settled after its air was made keeps the air, temperature, dose and UV of the older light: remake its end stages when the eccentricity moved past a tolerance, even if the star's age did not (bug) |  | Found by Lane 05, 2026-11-08. |
| GEN.1403 | Water-limited weathering: the carbon dioxide drawdown scales with land times the share of that land that is wet, instead of land alone (the GEN.909 option 3 lever) |  | Class M rarity report, 2026-11-08. |
| GEN.1404 | Land fraction against the water inventory: read the Cowan and Abbot mantle-water land fraction paper and the mantle water cap work, and report whether the model's 95 to 100 percent ocean cover for most sea worlds is what the physics gives |  | Class M rarity report, 2026-11-08. |
| GEN.1405 | Oxygen store delay: source the time an oxygen world takes to fill its stores before oxygen rises (EARTH_STORE_MOL, T9, L3, tuned to a 0.7 Gyr delay) from the geochemistry |  | Class M rarity report, 2026-11-08. |
| GEN.1406 | Inner-planet Earth-size share against the Kepler occurrence rate: the generator gives 3.3 percent of stars an Earth-size inner planet (0.8 to 1.5 Earth radii) against the 16.5 percent the survey implies |  | Class M rarity report, 2026-11-08. |
| GEN.1407 | K-dwarf zone edges against the greenhouse spectrum: the K5V ecosphere is cold in its inner third; split the probe by flux to see whether the edges or the greenhouse move with the star's light |  | Class M rarity report, 2026-11-08. |
| TEST.689 | src/tests/test_property_sheet.py fails on main: the oddities column is not on the sheet (bug) |  |  |
| TEST.691 | The ozone-loss radiation test fails on main (bug) |  |  |
| TEST.692 | src/tests/test_planet_classes.py::test_the_class_sweep_finds_every_stored_letter_is_the_derived_one fails on main (bug) |  |  |
| PERF.249 | A hypergiant system takes 5 to 11 minutes and 8.4 GB and ends with 0 planets (bug) |  | Bug Hunt first pass, 2026-11-07. |
| UX.289 | Any step that takes more than 5 seconds shows a progress bar from the 5-second mark: one shared mechanism for generation stages, the uncharted-sector index, the upper-level roll-up and every long step |  | Boss, 2026-11-08. |
| UX.291 | Phone width: the star-system map card is 338 px high against the half-window limit of 329 px because its three rows wrap (bug) |  | Lane 03, 2026-11-08. |
| UX.292 | At phone width an expanded star-system card with long rows scrolls itself instead of its rest region: the first view alone fills the card's maximum height (bug) |  | Lane 03, 2026-11-08. |
| GEN.1387 | Barred world redraw rate is 22.75 per system against a ceiling of 0.0183; barred_world_limit is 1 and 2 in test_barred_world[104] and [160] (bug) |  | Bug Hunt first pass, 2026-11-07. |
| TEST.695 | Stale or failing tests, group one: point_columns (sky_tier_draws), column round trip (irregular_moons), materials table, command docs (link-neighbors), TODO codes shown, JS unit and raw formatting, web_formulas, multistar registry MS4 and MS6, sql portability, db boundary values, update_database_step, literal insert, old_urls belt_ref, parent_link, container_plan, nebula_shape, object_rows, outline_reach, describe_render, property_sheet, api_link_ids (bug) |  | Bug Hunt first pass, 2026-11-07. |
| GEN.1390 | Generation sanity reds: moon history, class_sample P, daylight, prevalence habitable_world, climate sweep, atmosphere_state, biogenic oxygen, orbit_elements, life_schedule, population_streams, moon_disk_zone, sweep_bulk, forced_star (a subdwarf comes out V), invalid slot type (bug) |  | Bug Hunt first pass, 2026-11-07. |
| TEST.696 | Browser test failures: Help menu expectation (9 of 10), plane view and side strip clicks 'not visible', edit menu dialog never opens, fixture maps (9 of 227), generate flow phone-light, map card, density, hover, key, galaxy closed (bug) |  | Bug Hunt first pass, 2026-11-07. |
| TEST.697 | test_sweep_life.test_an_index_that_reads_the_letter_is_reported fails on main: the class-blind score check no longer reports an index that reads the letter (bug) |  | Lane 02, 2026-11-08. |
| TEST.700 | test_species_limits needs 66 minutes alone because its fixture hunts a civilized system (about 23 s a system, about 1 in 12 qualifies), so the train gate's 25 minute cap always cuts it; pin a system or speed the fixture (bug) |  | Lane 02, 2026-11-08. |
| PERF.251 | The block light totals bar sits at 100% while the upper-level roll-up and the uncharted-sector index run with no bar or live log; check whether the backfill leaves stale glow sums (bug) |  | Boss live observation, 2026-11-08. |
| PERF.252 | A full-galaxy run was killed by SIGKILL (out of memory) 45 minutes after the block light totals, with no log line and a stage label that does not match (bug) |  | Boss's full-galaxy run 20261108-020807-50fea318. |
| PERF.253 | test_each_tristate_flag_forced[B2V-True-MOONS] runs past the 900 s timeout in habitability_world.world_from (temperature_coverage.series_of_values) (bug) |  | Lane 05, 2026-11-08. |

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
| GEN.1377 | Per-body isotope record, non-fictional part: D/H, 15N/14N, 13C/12C, 18O/16O, 36Ar/38Ar, xenon isotopes and the short-lived Fe-60, Mn-53, I-129 and Pu-244 |  | |
| GEN.1420 | Major isotopes exist on every body: the per-body isotope record carries the natural isotope mix of every element the body holds (deuterium included), from the isotope table, with the measured anomalies of GEN.1377 on top |  | Boss, 2026-11-08. |
| GEN.1421 | Tritium and the cosmogenic trace isotopes exist: a steady-state inventory from production and decay for tritium, carbon-14, beryllium-10 and aluminium-26 on bodies and in air and water | GEN.1420 | Boss, 2026-11-08. |
| GEN.1422 | Superionic water and ammonia boundary table: the deep-layer boundary of ice giants and water worlds from the observed superionic phases, replacing the weird-worlds G4 text line |  | Boss, 2026-11-08. |
| GEN.1423 | Silicon carbide (moissanite) and cohenite (Fe3C) mineral rows for carbon-rich worlds |  | Boss, 2026-11-08. |
| GEN.1424 | The H2SO4-H2O vapour-pressure curve and the sulfuric acid cloud condensation window (the missing prerequisite of the Venus chemistry set) |  | Boss, 2026-11-08. |
| GEN.1425 | Venus chemistry set V1 to V10 as registered stand-ins and the chain gates that add them to any planet when the physics supports it | GEN.1424 | Boss, 2026-11-08. |
| GEN.1426 | Super-Earth mantle density splits: the high-pressure mantle phase changes (post-perovskite and the deeper transitions) as density steps in the interior of massive rocky planets |  | Boss, 2026-11-08. |
| GEN.1427 | Neutron-star nuclear pasta: the inner crust layer between the outer crust and the core, in the neutron star interior and radius model |  | Boss, 2026-11-08. |
| GEN.1428 | GEN.725 evidence-status column: every stand-in row says whether the thing is observed in a planet, observed in a lab, predicted only, or disputed |  | Boss, 2026-11-08. |
| GEN.1429 | GEN.725 WATCH lines for the predicted-only and disputed materials of the unseen compounds report | GEN.1428 | Boss, 2026-11-08. |
| GEN.1382 | Ore processes build, weathering and pegmatite deposits (bauxite and laterite, REE and Li pegmatites) | GEN.1371 | |
| GEN.1383 | Ore processes build, carbonatite, kimberlite and unconformity uranium deposits | GEN.1382 | |
| TEST.686 | Real-material closing gate: a composition census over a generated sample (each new element and mineral appears where physics says and nowhere it should not) and a check that every stand-in the set added has its GEN.725 row | GEN.1383, GEN.1377, GEN.1370, GEN.1369 | |
| GEN.1384 | Reduced-body recipes read the sulfide, nitride and carbide rows: the metal ground and the reduced taxa name pyrite, pyrrhotite, pentlandite, chalcopyrite, oldhamite, niningerite, alabandite, daubreelite, osbornite and sinoite instead of troilite and kamacite stand-ins |  |  |
| GEN.1391 | Minor core elements: Cr, Mn, P and Cu partition between the core and the mantle, read from the star element set (GEN.1360) and the core split (GEN.1364) | GEN.1364 | Lane 23 follow-on to GEN.1364. |
| GEN.1392 | Reduced cores and the dynamo energy scale: let rock FeO reach below IW-3.3 so Mercury-like Si-rich cores occur, and re-tune the compositional dynamo energy that stopped most small-body dynamos | GEN.1364 | Lane 23 follow-on to GEN.1364. |
| GEN.1395 | Tune a gas giant's carbon-to-oxygen ratio to observed Sun-like giants and to what physics expects |  | Boss, 9:43 PM EST on 2026-11-07. |
| GEN.1396 | The snow-line jump in solid density (2.17): move it from the old table to the disc's solids if that is the sounder physics |  | Boss, 9:43 PM EST on 2026-11-07. |

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
| GEN.1047 | Wire the oddity table and the weirdness dial: epoch, stored oddities text column and numeric columns (db/column_units.py), the "Odd: ..." line with its cause sentence on the object page, the API field, an admin force-oddity like FORCE_CLASS, the dial in settings (default 1, 0 turns all off) |  |  |
| GEN.1015 | Overlap-free stars P9: delete the GEN.1004 quick fix (built_tree_clears, the clearance branch of _companions_fit, _built_stars_clear), no shim | GEN.1014 |  |
| GEN.1014 | Overlap-free stars P8: a clearance ledger and planet, moon and belt sampling from it |  |  |
| GEN.1105 | Field members (individual bodies by the belt cut rule) with `<field ID>.<rank hex>` designations, status flags, and the sector and field pages | GEN.1103, GEN.1070 |  |
| GEN.1103 | AsteroidField columns, radius derived from age, the field density leaving 0, registry text changed |  |  |
| GEN.1092 | Wire the scatter phase and the belt origin: scatter tag on stored members, the outer disc fed from the scatter step, origin and formation shown on the belt page | GEN.1070 |  |
| GEN.1086 | System map and pages for irregular moons (eccentric, inclined, retrograde orbits on the moon-orbit view, a body page per irregular) and the clearance ledger taking the irregular zone as occupied space |  |  |
| TEST.660 | The class census reads irregular (captured) moons, and the moon class ranges are re-baselined after the last 8.1 train |  |  |
| GEN.1073 | Bring-to-epoch for belts: grinding recomputed from the stored age, members above the grinding scale kept | GEN.1069, GEN.1006 |  |
| GEN.1072 | White dwarf and giant-branch belt events wired into remnant_orbits and the engulfment step, replacing the uniform expansion factor of belts | GEN.1068 |  |
| GEN.1070 | Store the individual belt members (rows with parent link and rank), the designation display, the system map drawing them as dots at true orbits, and a body page per member; API version bump | GEN.1069 |  |
| GEN.1069 | Store the belt mass, Dc, age history text and formation snow line (DB.70 column, migration) and wire belt_mass into AsteroidBelt, replacing the three density words | GEN.477, GEN.1068 |  |
| GEN.1068 | Belt birth from the architecture: a belt is a gap or an unfinished tail in the planet chain, not a 10 percent slot roll, and sees the giants (GEN.728 slot loop, system.py) | GEN.1014 |  |
| GEN.1231 | White dwarf pollution flag reads the belt's scattered rate (is_polluted, pollution_rate_g_s) and the census band 25 to 50 percent | GEN.1068 |  |
| GEN.1150 | Lagrange trojans T2: the generator draws a trojan swarm at the L4 and L5 of each giant planet or brown dwarf from the scatter trojan share; new columns co_orbital_host_id, lagrange_point, libration amplitude and phase (db/column_units.py); orbit updater libration overlay; the object and sector pages and the designation rule "<host> L4 trojan region" | GEN.1006 | Trojans T2, epoch lock, train 8 |
| GEN.1151 | Lagrange trojans T3: whole planets and moons (Class M and Class J allowed by the mass rule) at L4 and L5 from the disc formation model only; no forced rate, the observed 5 percent upper bound (WR10) is a plausibility test only; folds in the weird-worlds W7 row | GEN.1150 | Trojans T3, epoch lock, train 8 or 9 |
| MAP.354 | Lagrange trojans T4: System Map icon at the centre of each Lagrange point that has objects, labelled with the host and point (for example Jupiter-L4); shares the geometry of the MAP.169 gravity layer | GEN.1150 | Trojans T4, with T2 |
| GEN.1058 | Wire the G1 giant-impact tilt and the ring plane following the equator | GEN.1047 |  |
| TEST.492 | Census re-baseline for every oddity share and the prevalence suite rows (one pass after train 9) | GEN.1047 |  |

### Physics rectification: test health, constant audits and the star mass function

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.395 | A planet distribution test with observed-rate bands, run after the class cut-over on the final classes and on galaxy-placed stars: giants and hot Jupiters per host mass, planets per mass bin, 1 to 20 Earth-mass planets per star, habitable-class bodies per star, and a rule for classes O, P and Q |  |  |
| TEST.417 | Whole-sector and whole-galaxy tests take 20 to 45 minutes: make them smaller or cached so a gate run finishes (bug) |  | found by Aux 8 |
| TEST.455 | A completeness gate: every numeric constant in the generation and physics modules has a stand-in registry row or is read from a table; the list of unregistered constants may only shrink | GEN.945 |  |
| GEN.945 | Audit every numeric constant the generator uses against the data we already download or hold: which are registered stand-ins, which come from a table, and which are unregistered hand-tuned numbers; ranked findings with a data source named for each |  |  |
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
| GEN.1309 | Free oxygen is rare: only 24 of 1446 worlds with air reach the 16 percent oxygen that classes M, H, L, O, P and R need, and those 24 sit outside the classes' pressure, gravity and temperature ranges (bug) |  |  |
| GEN.1311 | Fudge table part 1: a sensitivity harness that moves one stand-in at a time across its registered range and records how far each class quantity moves |  |  |
| GEN.1312 | Fudge table part 2: run the sensitivity harness on a lane machine (150 systems from each of 3 seeds) to build the stand-in by end-result matrix | GEN.1311 |  |
| GEN.1313 | Fudge table part 3: physics/data/fudge_table.json, its build script, its test and the registry rows for the widths | GEN.1312 |  |
| GEN.1314 | Fudge table part 4: class_fit reads fudge_table.json for the edge shoulders, the hand-picked CF and DM rows are deleted, the miss budget and matches for former holes are stored | GEN.1313 |  |
| GEN.1315 | Fudge table part 5: outcomes read the match (habitable count by threshold, the barred-world loop, +habitable_world, search filters) | GEN.1314 |  |
| GEN.1316 | Fudge table part 6: the stand-in ranking report gets a share-of-fudge column; condition-based-planet-classes.md and the /classes pages explain the fudge | GEN.1313 |  |
| GEN.1317 | Fudge table part 7: update checks the data against the epoch and regenerates the fudge table automatically when the epoch has changed | GEN.1313 |  |
| GEN.1318 | Fudge table revisit: rerun the sensitivity and minimal-fudge sweep after the oxygen fixes (GEN.1309, GEN.1306 parts 2 and 3) and re-derive the widths and cap tiers | GEN.1309, GEN.1312, GEN.1313 |  |
| GEN.1333 | A forced class G or P request is reachable in about 1 of 24 systems, and fails every time after the GEN.1323 envelopes: make the forced-class draw reach every class it offers |  |  |
| GEN.1343 | Super-puff planets (W15): an envelope-bearing sub-Neptune branch through the air solver, class and radius, so the super_puff module has a body to inflate |  |  |
| GEN.1347 | Tidal inflation of super-puff radius (the tidal-inflation equation in the GEN.1343 report, which needs the planet's drawn obliquity) so super-puffs stay inflated past 1 Gyr |  |  |
| GEN.1349 | Sealed-ocean calcium: the model leaches all of the crust's 8.5 weight percent calcium to solution as a free cation, which sets lime-buffered pH 12 to 14 in 60 percent of sealed oceans; find what share a leached basalt gives up as free calcium against what stays in carbonate and silicate minerals |  |  |
| GEN.1353 | Class and life census re-tag after train 15 (epoch 137): re-measure the census shares and re-baseline the class and life ranges the train moved |  |  |
| OPS.118 | Find what sent SIGTERM to the work horse on unit 11 of Boss's warm-up (status 15): no code path in the program sends it, so it is outside the program or a manual stop (bug) |  |  |
| OPS.115 | update.sh and set-permissions.sh warn about files that already exist on every update (config.json in both the checkout and the state folder): existing files are normal, so they are an info note, not a warning (bug) |  |  |
| DB.115 | A generation job stage fails with InterfaceError (0, '') after a dropped MySQL connection: find out why MySQL dropped the connection (the masking rollback and the sky-tier retry are fixed) (bug) |  |  |
| GEN.1283 | The habitability re-score gives a different first-band day length than generation stored (14.21 against 13.60 hours on a 'wide' system) (bug) |  |  |

### Physics rectification: winds and ring darkening

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1110 | Average wind speeds by day and by night in each of the seven latitude bands: the day-to-night thermal wind from the pressure difference, zonal jets and trade winds from the equator-to-pole temperature difference and the Rossby number, faded out below about 1 kPa |  |  |
| GEN.1055 | Ring age and pollution darkening: a ring stores its age and its ice share and darkness follow from pollution over time (R6); ring shadow climate hook (R8) |  |  |

### Physics rectification: stars, spectra and climate follow-ups

| ID | Item | Needs | Note |
|---|---|---|---|

### Physics rectification: older physics follow-ups, units and the storage rows

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.820 | Residuals of the bulk iron branch table (GEN.512): the coreless branch, a slope other than 1, the 0.05 dex scatter stand-in and a larger tuning sample |  |  |
| PERF.198 | Speed up util/detmath.matmul and cache equal gas_equilibrium.solve calls (finding 6) |  |  |
| GEN.835 | The tidal heat flux (tidal_flux_w_m2) is wired into the ground draw's sulfur rule |  |  |
| DB.101 | Storage and units for every physics rectification output: a column, a column_units entry and a migration shipped with its item |  |  |
| GEN.804 | Atmosphere escape and retention read the stored XUV exposure, flare irradiation and magnetopause (G2, step 1) |  | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
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
| GEN.1268 | Nebulae carry a real composition: gas and dust fractions by kind, metals from the place and age, ejecta from the progenitor |  |  |
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
| TEST.665 | test_class_sample.py::test_a_forced_rare_class_is_in_the_system[P] is red on main: no system with a forced class P world in 24 tries (bug) |  |  |
| GEN.301 | Headline of 8.1: the new generate scatters, accurate galactic generation and scattering including the halo | SEC.33, TEST.261 |  |
