# Version 8.1: Infinite zoom and the new generate scatters

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.1

## Goal

Two headlines. Infinite zoom: one camera flies from the whole galaxy to a star system (MAP.146, done: every sub-item merged and Core 2 checked each clause of its Done text against main, 3:20 AM EDT 2026-10-23 as the lane reported). The new generate scatters: accurate galactic generation and scattering out to the 262 kpc reach, including the halo, the satellites, the nebula stars and a realistic core (GEN.301). The version also holds the Phase 1 tail: the queue steps, reproducible galaxies, the description layer, tech levels, directives and every open bug and flaky test. A third strand (Boss, 2026-10-16, 11:57 PM EDT): the building blocks of the 9.0 composition headline, which go down to 8.1 because they change no stored value: the research lanes and designs (phase 1, GEN.431 to GEN.437 with GEN.409, GEN.425 and GEN.426), the data tables (phase 2, GEN.438, GEN.410, GEN.427, GEN.439), the star composition function (phase 3, GEN.440) and its sweep (phase 4, TEST.211); the master note is docs/design/composition/composition-model.md. The 8.1 stamp waits until all of it is done (REVISION_HOLD in scripts/bump_version.py).

## Where it stands

As of 12:34 PM EDT 2026-10-28: 1160 of 1590 items finished (73.0 percent), 430 open, 36 dropped (TEST.374 dropped as not reproduced on main, Bug 1 10:31 PM EDT as the lane reported; GEN.564 folded into GEN.764, Core 1 7:32 PM EDT; TEST.366 dropped as no longer reproducing after #3125, Bug 1 7:20 PM EDT as the lane reported; TEST.335 dropped as not reproduced on main, Bug 1 7:14 PM EDT as the lane reported; TEST.360 and TEST.361 dropped as not reproduced on main, Bug 1 6:35 PM EDT as the lane reported; TEST.337 and TEST.339 dropped as not reproduced on main, Bug Hunt 10:41 AM EDT as the lane reported; TEST.326 dropped as a duplicate of TEST.313, TODO 4:07 AM EDT; DB.95 dropped by TODO on lane B's finding, 2:53 AM EDT: nothing reads the view and DB.94 covers free moons; PERF.162 and PERF.166 moved to the later version, Boss 12:27 PM EDT; GEN.695 to GEN.698 dropped, Boss 1:22 AM EDT: class B never has hydrogen or helium). Open by phase: phase 1 83, phase 2 44, phase 3 24, phase 4 13, phase 5 55, phase 6 29, phase 7 17, phase 8 4, phase 9 31, phase 10 2, phase 11 131, phase 12 19 (phases 6 to 11 are the composition and sky colour suite pulled in from 8.4 and 9.0 by Boss at 5:48 PM EDT; phase 12 is the headline). Only this version is being built (Boss, 7:44 PM EDT 2026-10-22); the database chain finished with DB.85 (#2582, revision 0130, no generator epoch, so epoch 51 was still free then; since, 51 to 53 are taken: GEN.307, GEN.315 slice 1 and GEN.628, and 54 is next). Closing order: the audits TEST.185, TEST.184, TEST.159 and TEST.160 (merged, re-run as TEST.283), PERF.146 (merged as #2559) with PERF.164 and PERF.165, then SEC.33 and TEST.261 in their own threads, then GEN.301. Lane queues: docs/plan/lanes.md and reports/8-1-core-queues-2026-10-22.md.

## Phase 1: Foundations and open trees

Everything with no blocker outside its own thread: the queue steps, reproducible galaxies, the description layer, tech levels, the settings page (moved up from 8.2, Boss 1:43 PM EDT), directives, housekeeping, the bugs and the backfill weeding.

### Queue, caches and performance

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.40 | Weed out sectors by star density before the bright-star backfill |  | Parent; closes with its subitems. |
| PERF.35 | An interval or chunk ledger for untouched sectors once block-first backfill lands |  | Research: ledger size at full galaxy. |
| PERF.227 | Starting a single scatter part from the Generate console on a fresh download kills the web daemon: five 'Truncated or oversized response headers' errors, an error 500, and the job reads 'stopped without finishing' (bug) |  | Boss bug, scatter start kills the daemon. |
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
| GEN.229 | Repair or re-roll the wide binaries that already break the stability rule |  |  |

### Bugs and flaky tests

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.344 | Gate: browser and WebGL tests run as their own serial phase, with a worker cap by free memory |  |  |
| TEST.190 | A failing session fixture makes every dependent test wait its full 900 s timeout: fail the dependents at once (bug) |  | Test-infrastructure finding. Bug Hunt proposed (Boss started the lane, 10:20 PM EDT). |
| TEST.237 | The wheel-anchor browser test failed alone at box load 9.7 (call 5.47 s) after the TEST.231 fix merged (bug) |  |  |
| TEST.290 | test_work_units and test_galaxy_seed each fail once under the gate load: a galaxy run stores each stage once, and a new seed replaces the old only before any sector (bug) |  |  |
| TEST.229 | The browser phase of a gate stalls for minutes at one test with a headless_shell gpu process near 200% CPU (bug) |  |  |
| GEN.396 | Revisit the Yerkes class III luminosity range: 50 to 1000 solar luminosities sits above many dim red giants |  | Core 3 finding, 2026-10-16 10:00 PM EDT; any Core lane. |

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
| GEN.613 | Investigate the limits of the density algorithm: for which mass, size, flatness and bulge share the galaxy density model stays valid, and what it does outside them |  | Boss 7:33 PM EDT 2026-10-21; read-only report. |
| GEN.618 | The stellar-to-halo inversion gives a 6e13 solar-mass halo for 2e11 solar masses of stars, the shape dials' mass ceiling: check the Moster inversion at the top end |  | Density limits report (GEN.613). |

## Phase 2: Infinite zoom

The fly-through Galaxy Map (MAP.146) and the nebula data its tiles carry.

### Fly-through map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.75 | The mini map as a second engine view |  | Locked second camera on MAP.64. |
| MAP.140 | Double-click on a selected object goes there and opens its information |  | GitHub issue [#714](https://github.com/dwhagar/planetGen/issues/714). |
| MAP.141 | Context around the selection: faint neighbours, and the sectors above and below |  | GitHub issues [#718](https://github.com/dwhagar/planetGen/issues/718) and [#716](https://github.com/dwhagar/planetGen/issues/716), one context view. The #716 bug label was overruled. |
| MAP.269 | Stars on the maps do not glow: a star further than about 3 sectors is a 1.5 px point with no halo and a mostly white core (bug) |  | Boss bug 8:51 PM EDT |
| MAP.257 | Infinite zoom is still clunky: the zoom zones and segments other than sectors get in the way (bug) |  |  |
| PERF.200 | Galaxy Map node cut bookkeeping: about half the busy time in the profile is string building, geometry recomputation and cache re-insertion |  | optional until a lag sighting; Core 2 |
| MAP.171 | The Galaxy Map keys a massive phenomenon's visibility to its map luminosity (replaces "always lit") |  | Boss 2026-10-10 21:38Z. Foundations lane 2 after MAP.152. |
| TEST.135 | Big-galaxy query-budget test for the nebula reads | GEN.202 | Foundations lane 1. |
| MAP.156 | View one layer or a range of layers top-down from the galaxy view, as a secondary option |  | Boss 2026-10-09 23:17Z; secondary option, not a main control. |
| MAP.160 | Quantise the Galaxy Map GPU buffers (deferred) |  | Deferred; MAP.147 recommendation step 3. |
| MAP.189 | The Galaxy Map opens flown in to the created sectors |  |  |
| MAP.397 | Some View menu options show only the word "Forward" and cannot be clicked (bug) |  |  |
| DB.109 | Optional parity follow-ups: refresh the parity file after a generation run, and keep a before-image on the edit paths so an edit is a parity delta, not a group recompute |  |  |
| API.36 | Optional: a region_names search table and a "Go to" box, refilled when the naming key changes (regions chain, step 6) |  |  |


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
| GEN.669 | Residuals of the generator laws after GEN.450 and GEN.270: the star-formation profile shares and the Algol rate against their observed windows |  | Core 4, 10:20 PM EDT 2026-10-23; no epoch until a law changes. |
| GEN.708 | White dwarf merger products draw their 0.8 to 1.1 solar mass bump: the star model gives the white dwarf its mass before the member tree, so the tree only flags 10 percent as merger products |  |  |
| GEN.707 | Register the white dwarf binary stand-ins in the gap plan registry: kind, evidence label, knobs, validity range and falsifier test for each |  |  |
| GEN.331 | Substellar companions: brown-dwarf close (the desert) and wide companions, brown-dwarf pairs and planetary-mass companions |  | Split from GEN.275. |
| GEN.332 | Circumbinary and S-type planets, circumbinary disks and rogue planet pairs at fill | GEN.315 | Split from GEN.275. |
| GEN.640 | GEN.316 part (c): delete BinaryStarProxy, WideBinaryPair and the writers of the binary_* columns, once circumbinary and S-type planets come from the tree | GEN.332 | Split from GEN.631. |
| GEN.275 | Substellar and circumbinary types: brown-dwarf and planetary-mass companions, circumbinary disks and planets | GEN.331, GEN.332 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 19. Core. |

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
| GEN.391 | Superbubbles shape the stars: the design of the bubble field as a density modulation, with its bounds and touch points checked |  | Boss 9:27 PM EDT: step 1 of the superbubble density chain; Core 3 proposed. |

### Star composition

| ID | Item | Needs | Note |
|---|---|---|---|

## Phase 4: New scatters: components, clusters and audits

Satellites, clusters and wanderers as components, the globular cluster chain and the audits that check what the scatters make.

### Components and wanderers

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.371 | Feasibility study: generate the satellite galaxies through the scatter, making only the region each one occupies |  | Core 1; Boss 4:02 PM EDT: feasibility study. |
| GEN.534 | A satellite's young stars reach the system-type tables: the young-star system types name `satellite_young` |  | Residual of GEN.528: the system-type tables do not name satellite_young. Core 3 proposed. |
| GEN.273 | Young groups and trapezia: the `group` node, offset and velocity columns, a decay age and the remnant draw |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 17. Core. |

### Globular clusters

| ID | Item | Needs | Note |
|---|---|---|---|

### Whole-sphere audits

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.162 | Multiplicity audit: companion states, row bookkeeping, stability, determinism, no planets in a scatter, fill reuse, motion and the census anchors | GEN.315 | Multi-star systems at scatter time, Boss 2026-10-15. Step 9. Quality. |
| TEST.164 | Catalogue type audit: type counts, region gates, the N tail, interaction consistency and group decay | GEN.275 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 20. Quality. |

### Composition proofs

| ID | Item | Needs | Note |
|---|---|---|---|

## Phase 5: New scatters: core, halo, interface update and docs

The realistic galactic core and the halo tier, the interface update for the spherical grid, the headline item and the docs that follow.

### Galactic core and halo

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.377 | Scatter limits by place for the neutron star and black hole scatter below the cut: the core and halo tiers apply to it, not one global mass cut |  | Core 3 finding after GEN.333. |
| GEN.303 | Per-cell system cap as a configuration option: `generation.cell_cap`, default 5,000, lowest 500, highest 50,000 |  | Core 3 wires it; Boss 1:47 PM EDT. |

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
| GEN.644 | Proposal: tie the PHI-6 into the planet classification (class letters, a habitability tier, or a suffix) and the rule for where a class and a PHI-6 colour must agree |  |  |
| TEST.300 | End-to-end browser tests for the habitability system: symbols on lists, sector colouring and habitable-world markers |  |  |

### Radiation chain for habitability (GEN.445)

Boss, 5:42 PM EDT 2026-10-23: build the entire GEN.445 chain in 8.1; testing mode, not production, so new generator epochs, database versions and API versions are fine. These 20 items moved here from 8.4 and 9.0 in dependency order. The spine (lane I, to be started): GEN.565, GEN.557, GEN.428, GEN.429, GEN.452 (part 2, PR #2613, held until now for the 8.4 epoch; it now rides the spine's first epoch), GEN.411. Side items that feed the spine: GEN.553 (lane A, after UX.101), DB.73 (lane B, after UX.100), GEN.412 (lane D, after GEN.639). Lane F, radiation constants and ozone: GEN.453, GEN.413. Lane G, bands: GEN.417, GEN.414, GEN.443. Lane H, surface dose: GEN.459, GEN.423, GEN.461, GEN.441, GEN.444. Lane E: GEN.445 with the help text and tests.

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.553 | One star-like body draw: origin, birth, state and [Fe/H], with the initial-final relation, the remnant mass, radius and age rules in one place |  |  |
| GEN.565 | The mass object core and its three family classes: StellarObject, PlanetaryObject and SmallBody, with state and parent as data, replace the nine kind classes |  |  |
| GEN.557 | Wanderers: a rogue planet becomes a planet with no parent, drawn by the one body chain with an origin table, a notional formation host, the population's age and motion | GEN.553, GEN.565 |  |
| GEN.444 | Surface radiation: the ground-level dose from the air's shielding and the ground's own radiogenic content |  | Boss 11:58 PM EDT; replaces the proxies of physics/radiation.py; generator epoch bump. |

### Docs

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.50 | Docs: the schema, the multi-star note and the catalogue say where the members live | GEN.264, TEST.164 | Multi-star systems at scatter time, Boss 2026-10-15. Step 10. Quality. |
| DOC.51 | Docs: the cell key in the coordinate system and schema docs, the region note moves to built, the API states the cap and the ladder |  | Spherical region query, stage C. Step 8. Quality. |
| DOC.18 | Last pass: re-run the documentation sweep when the other DOC items are done and Phase 1 is complete | DOC.50, DOC.51 |  |

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
| GEN.405 | Metallicity and rates by position from the galaxy's own model: the mass-metallicity relation and the local rates replace the fixed gradient and the quoted-at-the-Sun values |  | GEN.368 build, Boss approved the table 10:11 PM EDT; Core 1 proposed for GEN.400 to GEN.403, Core 3 or Core 4 for GEN.404 and GEN.405 (profile.py, activity.py), any lane for GEN.406. Core 4 proposed (coordinator 10:16 PM EDT). |
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
| TEST.244 | Tests: re-primary and reparent invariants, frame round trips, momentum conservation, a moon's velocity as the sum of its hops, a locked body's rate follows a primary change | GEN.576, GEN.578, GEN.581, GEN.583 |  |
| GEN.583 | reparent: the state-preserving change of a body's parent, shared by re-primary at an event, an ejection, a capture and a merger | GEN.578 |  |
| GEN.582 | Origin and present-state split: fields marked origin or state in the Field table, the draw record held until validation finishes, and the zone and class re-run from the seed | GEN.565, GEN.570 |  |
| GEN.581 | primary_first(members): before a system is validated the heaviest branch becomes the primary and the others are recalculated, one function for the tree, the fill, moons and free pairs | GEN.315 |  |
| DOC.71 | Docs for the System rule: database-schema.md, the system view in docs/design and the header of system_tree.py | GEN.576, DB.78, GEN.578 |  |
| TEST.242 | Tests: invariants for any system at any time: zero mass-weighted offset sum, add and remove returns exactly, a parent change keeps every absolute position, one member equals the binary formula bit for bit | GEN.576 |  |
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
| GEN.559 | Unbound comets: a size law and a cut, composition from the comet table, derived activity and velocity | GEN.557 |  |
| GEN.558 | Moons of wanderers: ejection survivors and isolated-disk moons are table rows, and the fixed 20% boolean is gone | GEN.556, GEN.557 |  |
| GEN.556 | `populate(host)`: planets, belts, comets and moons come from the class tables by one recursive routine, and the system pass is absorbed | GEN.315, GEN.555, GEN.565, GEN.581 |  |
| GEN.455 | Airless bodies and comets draw an exosphere or coma from their surface composition, temperature and the star's wind: the sodium and potassium tails and the production rates |  | Atmosphere research lane new work 6 (parts 04 section 5 and 05); generator epoch bump. Core 4 proposed. |

### Mixture and stored values

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.594 | Each AD dataset (AD1 to AD8) ships its best possible fallback, used only when a run is told explicitly to disregard missing files |  |  |

### The Quantity call-site sweep

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.349 | API output and stored-column readers take Quantity; the grep test fails on a hand-written conversion factor; the old helpers are removed |  | Split from GEN.300. Moved from 8.3. |
| GEN.300 | Every quantity in the code moves onto the Quantity class, and hand conversions are removed | GEN.349 | The call-site sweep, partly covered by GEN.571 and GEN.574. Moved from 8.3. |

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
| GEN.477 | The asteroid belt and field composition is drawn from the star, the distance and the belt's mass, not four names drawn uniformly from a list of 84 (bug) |  | Body composition research lane finding F11 (tables E1 to E4): belt and field composition drawn from the star, distance and belt mass; generator epoch bump. Core 4 proposed. |
| GEN.471 | Clathrate hydrates and eutectic brines for icy crusts and ocean moons: the dissociation pressures and the liquid range below 273 K |  | Body composition research lane part F item F4 (later); generator epoch bump. Core 4 proposed. |

## Phase 9: Stored scatterers, the text and the proofs

The columns for ozone and the aerosols, the atmosphere text written from the mixture, the consistency sweeps and the docs.

Moved into 8.1 from 8.4 phase 4 by Boss, 5:48 PM EDT 2026-10-23: build the entire composition and sky colour suite now, testing mode (new generator epochs, database versions and API versions are fine). Lanes J to M take it (J bulk composition, K surfaces and small bodies, L the atmosphere mixture, haze and proofs, M atmosphere and sky colour).

### Consistency

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.647 | The class is assigned after the body is complete: classify(body) from the stored atmosphere, temperature, terrestrial composition and gas composition, for planets and moons |  |  |
| GEN.690 | Delete PLANET_CLASSES generation fields, PLANET_CLASS_PROBABILITIES, LIFE_PLANET_CLASSES, PLANET_CLASS_MAX_LIFE_STAGE and MOON_BLACKLIST; update the golden outputs | GEN.689 | Class design step C15. |
| GEN.689 | The class cut-over: classify at the end of generation, validation becomes total and derived, one generator epoch, a migration backfills the column | GEN.681, GEN.682, GEN.684, GEN.685, GEN.687 | Class design step C14. |
| GEN.688 | Generation requests name conditions, not class letters: slots and the admin class change take conditions, +habitable_world steers on the PHI-6, the directive counts is_habitable | GEN.681, GEN.682, GEN.684, GEN.685 | Class design step C13. |
| GEN.687 | Rogue planets get their class from classify; delete choose_rogue_planet_class | GEN.558, GEN.682 | Class design step C12. |
| GEN.685 | Water comes from the composition and the origin; delete WATER_MASS_FRACTION[class] |  | Class design step C10. |
| GEN.684 | Albedo and greenhouse come from the composition; delete the two class ranges |  | Class design step C9. |
| GEN.682 | Radius, density and gravity come from the bulk composition; delete radius_range, size_mode and density_range | GEN.681 | Class design step C7. |
| GEN.681 | The kind of a planet or moon comes from formation (mass against the snow line and the host) and replaces the class draw | GEN.556 | Class design step C6. |
| TEST.314 | Only classify writes the class column, and the class is recomputed when an input is edited |  | Class design step C5. |
| GEN.648 | The generation request for a habitable world (+habitable_world) runs an ordered ladder of nine adaptations before it gives up: redraw the surface, restrict the planet, more starting water, aquifers, regenerate the planet, adapt the system, add a brighter primary, regenerate the star, then the best result |  |  |
| DB.90 | Planets and moons store the derived class and a reclassification step keeps it true after an edit or a rebuild | GEN.647 |  |
| UX.150 | Pages, lists, search, filters, the API and the map markers show and filter by the derived class and the PHI-6 colour; the class wording and the descriptions follow the rule table | DB.90, GEN.647 |  |
| DOC.76 | Docs for the condition-based classes: the rule table, the planet class pages and the help | DB.90 |  |
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
| GEN.476 | Aluminium-26 and iron-60 follow from the star's formation environment, not from the Solar System's own value and range |  | Boss 1:43 AM EDT, Stand-in only: R4, Al-26 and Fe-60 from the formation environment. |

### Weight research (atmosphere)

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.483 | Source the cosmic-shoreline constants and the pressure caps by retention regime (AR4) |  | Atmosphere lane research need R4 (rows d6); replaces the seeded bracket of the DESIGN rows. |
| GEN.492 | Source the ground radioisotope abundances by host and the potassium depletion by body class (AR13) |  | Atmosphere lane research need R13 (rows o2, o3, o4, o6); replaces the seeded bracket of the DESIGN rows. |
| GEN.493 | Source the initial-rotation distribution and the saturation time against stellar mass for the host-star UV regime (AR14) |  | Atmosphere lane research need R14 (rows p6); replaces the seeded bracket of the DESIGN rows. |
| GEN.494 | Source the brown-dwarf type against temperature and the T/Y boundary near 500 to 700 K (AR15) |  | Atmosphere lane research need R15 (rows g1, g4); replaces the seeded bracket of the DESIGN rows. |

### Tuned draws and selection models

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.510 | Tuned draw builder B2: delivery and water worlds |  | Tunable gaps report: Tuned draw builder B2: delivery and water worlds. |
| GEN.515 | Tuned draw builder B7: ground and petrology |  | Tunable gaps report: Tuned draw builder B7: ground and petrology. |
| GEN.519 | Tuned draw builder B11: aerosol amplitude |  | Tunable gaps report: Tuned draw builder B11: aerosol amplitude. |
| GEN.595 | Replace the Solar System comet population shape (GEN.520 D6) with a derived distribution |  | Boss, 2026-10-18 (fork 8, flagged stand-ins only). |
| GEN.521 | Tuned draw builder B13: belts |  | Tunable gaps report: Tuned draw builder B13: belts. |
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
| GEN.1049 | Wire W5 resonance heating and the W3 high-eccentricity and W4 misaligned-orbit tails in planet_architecture | GEN.1047, GEN.1016 |  |
| GEN.1015 | Overlap-free stars P9: delete the GEN.1004 quick fix (built_tree_clears, the clearance branch of _companions_fit, _built_stars_clear), no shim | GEN.1014 |  |
| GEN.1014 | Overlap-free stars P8: a clearance ledger and planet, moon and belt sampling from it |  |  |
| GEN.1107 | Bring-to-epoch for fields: age advances, radius follows, retirement past 1 ly marks or deletes the row (Boss to pick; default retire) | GEN.1103, GEN.1006 |  |
| GEN.1106 | Field render: a 3D view for the field page (today _NO_VIEW), the sector map drawing the cloud or annulus plus the individual points at v x age, the render seed from the uid not the database id | GEN.1105 |  |
| GEN.1105 | Field members (individual bodies by the belt cut rule) with `<field ID>.<rank hex>` designations, status flags, and the sector and field pages | GEN.1103, GEN.1070 |  |
| GEN.1104 | Supernova host rule: a young supernova remnant and its remnant core draw the unbound belt as a field, like NEBULA_HOST_RULES | GEN.1103 |  |
| GEN.1103 | AsteroidField columns, radius derived from age, the field density leaving 0, registry text changed |  |  |
| GEN.1092 | Wire the scatter phase and the belt origin: scatter tag on stored members, the outer disc fed from the scatter step, origin and formation shown on the belt page | GEN.1070, GEN.1071 |  |
| GEN.1086 | System map and pages for irregular moons (eccentric, inclined, retrograde orbits on the moon-orbit view, a body page per irregular) and the clearance ledger taking the irregular zone as occupied space | GEN.1084 |  |
| GEN.1085 | Belt ledger link: a captured member leaves its belt's mass total; until then captured bodies come from the outer-disc reservoir only | GEN.1084, GEN.1069 |  |
| GEN.1084 | Store the individual irregular moons and the new moon columns (migration, API version), and call the stability test for irregulars in drop_unstable_moons and in bring-to-epoch when a planet moves | GEN.1083, GEN.1070 |  |
| GEN.1083 | Wire irregular moons into generate_moons as a second tier above the regular moons (MOON_MIN_MASS_KG left alone for the regular list) with the captured-origin weight | GEN.1082 |  |
| GEN.1082 | The architecture records its capture windows (the instability event, the number of giant-giant close encounters, the gas-disc lifetime already kept in DiskSystem) as a field on the system |  |  |
| GEN.1073 | Bring-to-epoch for belts: grinding recomputed from the stored age, members above the grinding scale kept | GEN.1069, GEN.1006 |  |
| GEN.1072 | White dwarf and giant-branch belt events wired into remnant_orbits and the engulfment step, replacing the uniform expansion factor of belts | GEN.1068 |  |
| GEN.1203 | Polluted white dwarf census: 25 to 50 percent of white dwarf systems polluted once belt birth sees the giants | GEN.1068, GEN.1072 |  |
| GEN.1071 | Outer (Kuiper-type) disc built as a belt from the same chain, linked to the debris descriptor so the descriptor and the belt are one thing | GEN.1068, GEN.1069 |  |
| GEN.1278 | Outer-disc mass and width from a planetesimal-formation model, so A and F stars hold the debris disks surveys find (A stars have none against 24 percent, Sun-like 8 against 14) | GEN.1071 |  |
| GEN.1070 | Store the individual belt members (rows with parent link and rank), the designation display, the system map drawing them as dots at true orbits, and a body page per member; API version bump | GEN.1069 |  |
| GEN.1069 | Store the belt mass, Dc, age history text and formation snow line (DB.70 column, migration) and wire belt_mass into AsteroidBelt, replacing the three density words | GEN.477, GEN.1068 |  |
| GEN.1241 | Wire the system body view: the ledger and barycentre function read the typed tables for the system page |  |  |
| GEN.1068 | Belt birth from the architecture: a belt is a gap or an unfinished tail in the planet chain, not a 10 percent slot roll, and sees the giants (GEN.728 slot loop, system.py) | GEN.1014 |  |
| GEN.1231 | White dwarf pollution flag reads the belt's scattered rate (is_polluted, pollution_rate_g_s) and the census band 25 to 50 percent | GEN.1068 |  |
| GEN.1023 | Two planets of a close pair in triple seed 208 sit 9.93 mutual Hill radii apart where the generator promises 10 (bug) |  |  |
| GEN.1150 | Lagrange trojans T2: the generator draws a trojan swarm at the L4 and L5 of each giant planet or brown dwarf from the scatter trojan share; new columns co_orbital_host_id, lagrange_point, libration amplitude and phase (db/column_units.py); orbit updater libration overlay; the object and sector pages and the designation rule "<host> L4 trojan region" | GEN.1006 | Trojans T2, epoch lock, train 8 |
| GEN.1151 | Lagrange trojans T3: whole planets and moons (Class M and Class J allowed by the mass rule) at L4 and L5 from the disc formation model only; no forced rate, the observed 5 percent upper bound (WR10) is a plausibility test only; folds in the weird-worlds W7 row | GEN.1150 | Trojans T3, epoch lock, train 8 or 9 |
| MAP.354 | Lagrange trojans T4: System Map icon at the centre of each Lagrange point that has objects, labelled with the host and point (for example Jupiter-L4); shares the geometry of the MAP.169 gravity layer | GEN.1150 | Trojans T4, with T2 |
| GEN.1048 | Wire W1 (tilt-driven seasonal climate), W2 (Venus-type spin), W10 (iron-rich cores) and W15 (super-puffs) into the bulk, spin and climate inputs | GEN.1047 |  |
| GEN.1050 | Wire W14 (the dying-planet flag) and W9 (rings from decayed moons) | GEN.1047 |  |
| GEN.1051 | Wire the W19 nightside collapse and the W20 bistable branches (the history draw between branches) | GEN.1047 |  |
| GEN.1058 | Wire the G1 giant-impact tilt and the ring plane following the equator | GEN.1047 |  |
| TEST.492 | Census re-baseline for every oddity share and the prevalence suite rows (one pass after train 9) | GEN.1047 |  |

### Physics rectification: test health, constant audits and the star mass function

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.395 | A planet distribution test with observed-rate bands, run after the class cut-over on the final classes and on galaxy-placed stars: giants and hot Jupiters per host mass, planets per mass bin, 1 to 20 Earth-mass planets per star, habitable-class bodies per star, and a rule for classes O, P and Q | GEN.689 |  |
| TEST.417 | Whole-sector and whole-galaxy tests take 20 to 45 minutes: make them smaller or cached so a gate run finishes (bug) |  | found by Aux 8 |
| TEST.455 | A completeness gate: every numeric constant in the generation and physics modules has a stand-in registry row or is read from a table; the list of unregistered constants may only shrink | GEN.945 |  |
| TEST.461 | test_db_column_round_trip fails on planets.porosity: the column is never written (bug) |  |  |
| GEN.945 | Audit every numeric constant the generator uses against the data we already download or hold: which are registered stand-ins, which come from a table, and which are unregistered hand-tuned numbers; ranked findings with a data source named for each |  |  |
| GEN.951 | Companion mass ratio depends on the primary mass and the orbital period, with the twin excess ([R3108: Table 13]) instead of a flat 0.1 to 1 (low priority) |  |  |
| GEN.932 | The star mass function epoch: golden, design docs, census targets and tags, retire the old tags; items 3 to 7 land as one merge in one lock turn |  |  |
| GEN.922 | The star mass function is built right: initial mass function and present-day populations match the observed curve segment by segment, for every mass range and every star-like object type (parent, replaces the one-number GEN.914 tweak) |  |  |
| GEN.918 | The galaxy's star count follows its own star-formation history: the density model counts systems per parsec cubed while the history counts stars formed, so the galaxy holds 0.46 of the stars its history implies (root of GEN.900, GEN.901, GEN.902) |  |  |

### Physics rectification: life classes, chemistries and their pages

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.915 | Life depends only on the conditions, never on the planet class letter: every combination of conditions that can host life can evolve it on any class |  |  |
| TEST.430 | Class census: a standing seeded report of every planet and moon class share and overall colour per star type and zone, with stored expected ranges and anomaly flags |  |  |
| GEN.862 | Design: the life class system, merged as docs/design/habitability/life-classes.md (nine classes, the stage ladder, the ledger, scoring and storage) |  | Boss, 9:10 PM EDT 2026-10-28: carbon life classes move into 8.1. |
| GEN.878 | Chemistry-specific energy: the pigment's absorption band times the star's spectrum is the multiplier on LIGHT_SHARE, instead of one bolometric curve for every metabolism |  | Life chemistry review; moved to 8.1 by Boss, 8:00 PM EDT 2026-10-28. |
| GEN.879 | A far-red oxygenic chemistry (chlorophyll d and f) that K and M stars can draw, so the most common stars can have oxygenic life and animals | GEN.878 | Life chemistry review; moved to 8.1 by Boss, 8:00 PM EDT 2026-10-28. |
| GEN.880 | Chemolithoautotroph producers (H2 methanogens and acetogens, sulfur and iron oxidisers, radiolytic H2) replace Melanin as the producer under ice and in aquifers, with an oxidant budget that caps aphotic animals |  | Life chemistry review; moved to 8.1 by Boss, 8:00 PM EDT 2026-10-28. |
| GEN.876 | Adaptation traits built: drawn in the life pass, stored in a text column, listed on the object page (adapted to: radiation, sulfide) and read by the habitability scores |  | Life chemistry review; moved to 8.1 by Boss, 8:00 PM EDT 2026-10-28. |
| UX.262 | Planet and moon page Life section: the adaptation traits of the living class (the last piece of UX.237) | GEN.876 |  |
| GEN.881 | Retinal becomes a photoheterotroph supplement that needs a carbon fixer, and Melanin becomes the radiotolerance adaptation with dose in its draw, not a radiation food | GEN.876 | Life chemistry review; moved to 8.1 by Boss, 8:00 PM EDT 2026-10-28. |
| TEST.436 | Prevalence tests for the chemistries, the adaptation traits and the animal routes: shares per chemistry per star letter, per route and per trait inside researched ranges | TEST.430, GEN.879, GEN.880, GEN.876 | Life chemistry review; moved to 8.1 by Boss, 8:00 PM EDT 2026-10-28. |
| GEN.888 | The acid medium model: a sulfuric acid cloud layer and acid pools from the temperature profile and the SO2 and water shares, filling the dead solvent() sulfuric branch |  | Life class system (order final, Boss 9:20 PM EDT). |
| GEN.863 | Life class scores: every world stores its habitability per life class, shown on the object page, the API and the map filters | GEN.862 | A migration; moved to 8.1 by Boss, 9:10 PM EDT 2026-10-28. |
| GEN.854 | Most class M planets do not look habitable: only a tiny share are, and most need a hardened pressure suit (bug) |  |  |
| GEN.1209 | Several life types coexist on one planet: a per-body list of living classes, each with its stage, animal grade, kinds (microbes, flora, fauna) and prevalence | GEN.862 |  |
| GEN.1213 | Atmospheric composition accounts for the whole mix of living classes on a planet: each class changes the air in its own way, summed in timeline order | GEN.862 |  |

### Physics rectification: infrared opacity, tides and the class re-sweep

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.852 | Infrared opacity leftovers: Venus clouds and wings (the IR3 far wing is a labelled stand-in), more gas tables, tropopause quantisation, and species without data in the cp over R path (OCS and others) |  |  |
| GEN.1059 | The GEN.412 inversion of a stored giant row ignores the giant's own heat, so its stored albedo and greenhouse factor no longer reproduce its temperature from starlight alone (bug) |  |  |
| GEN.1019 | Venus clouds and the CO2 far-wing line shape in the infrared model (the IR3 far wing is a labelled stand-in) | GEN.852 |  |
| GEN.842 | Follow-ups to elliptical orbits, tidal heating and planet spin: the approximations GEN.727, GEN.737 and GEN.731 left |  | Math 2 follow-up |
| TEST.499 | RP9: rogue yield calibration against microlensing and the terrestrial scale: rogues per star against [R4254, R4255] (central 21, lower 1 sigma 8 over the probed range), the [R3191] Jupiter limit and [R0258]'s 2.5 terrestrial bodies per star |  |  |
| GEN.1144 | Post-train-6 class re-sweep on the new main, then decide the E and F merge, the T versus N boundary and the H and L reach |  |  |

### Physics rectification: names, designations and Object IDs

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1141 | Designation history part e: bodies, belt members and captured bodies follow a parent ID change, and tombstones on delete |  |  |
| API.38 | API and search for local names: name is the local name, a new sector field carries the place, a typed name returns every match with its sector, an ID typed in full still opens its object | GEN.1254 |  |
| UX.259 | Uncharted star and phenomenon pages, lists, infobox and map labels show the local name with its sector, not 'Uncharted G2V star <row id>' | GEN.1254 |  |
| GEN.1257 | Charting keeps the old ID-derived name as the first Formerly line when a word-salad name replaces it | GEN.1254 |  |
| GEN.1256 | One uid-derived name scheme for every ID-named object: rogues, nebulae, remnants, quasars, comets, asteroid fields and black holes stop using the position ID as a name | GEN.1254 |  |
| GEN.1254 | Every scatter star and scattered phenomenon gets an Object ID (uid) when the scatter is made: uid column on bright_stars and phenomenon_scatter, serial = rank in the sector |  |  |

### Physics rectification: strength, rings, trojans and small-body follow-ups

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1186 | Small evidence updates to registered stand-ins from the strength papers review |  |  |
| GEN.1185 | Interior activation energy range 260 to 390 kJ/mol instead of 200 to 300, with the reference viscosity retuned |  |  |
| GEN.1179 | Optional: Hoek-Brown geological strength index as a body parameter for rock |  |  |
| GEN.1178 | Strength of a body sampled from its range so the round boundary is fuzzy, as in the real Solar System | GEN.1159 |  |
| GEN.1175 | Porosity weakening of body strength from a sourced law instead of the placeholder exponential |  |  |
| GEN.1174 | Decide the pre-main-sequence time definition: the generator's time to the main sequence is 3 to 4 times BHAC15 for 0.1 to 0.3 solar masses and 1.5 to 2 times from 0.4 up |  |  |
| GEN.1171 | Ring lifetime by ring mass and planet instead of one fixed decay, because fallen-moon rings now persist on 13 percent of rocky planets |  |  |
| GEN.1167 | The stellar-to-halo relation [R3181] has 0.3 to 0.4 dex more stars than [R0319] in halos of 1e10.75 to 1e11.5 and 0.06 to 0.26 dex fewer in halos of 1e12.5 to 1e14.5: compare the published fits and pick the best-supported relation |  |  |
| GEN.1204 | Dark members of an orbit tree (brown dwarf, neutron star, black hole in dark_members) store density, binding energy and the shape values the stars, rogues and comets store |  |  |
| GEN.1161 | Asteroid field mass rule and a member ledger like the belts' (individual bodies above 1e19 kg stored, the rest drawn on demand) |  |  |
| GEN.1159 | Store ice fraction, bulk density, roundness ratio, yield strength and binding energy on planets, moons, asteroids and belt or field members (columns, migration, API version) |  |  |
| GEN.1201 | A per-body shape-setting temperature: the roundness engine uses the temperature a body took its shape at, but no body stores one, so strength is the cold strength today | GEN.1159 |  |
| GEN.1157 | Earth-type moons around gas giants: how often a giant has a moon the size of a planet, and how big, from the physics of moon formation and capture |  |  |
| GEN.1153 | Day and night state for stars and rogue planets (follow-up to GEN.1097, which covers planets and moons): a star stores the constant "no night"; a rogue planet gets the same state, night hours and solar days as an in-system planet from its own spin, with no star the answer is no solar day |  | Aux 2 report 10-30 |
| GEN.1152 | Mass literals that differ from constants.py (changing them changes output): moon_composition JUPITER_MASS_KG 1.89813e27, ring_origin star_mass_kg 1.989e30, giant_interior 317.8 and EARTH_MASS_JUPITER 317.83, stellar_abundances JUPITER_MASS_SOL, tuning 317.8; use the constants.py values with an epoch |  | Sources 1 report 10-30, train-bound |
| GEN.1130 | The +intelligent_life and +habitable_world request options become probability tweaks only: they raise the odds of such a world in the generated system and never force a placement, a ladder of adaptations or a civilization |  |  |
| GEN.1129 | Fingerprint split: a sector's digest (GEN.58) becomes "made here" and "arrived here" so rogues landing from other sectors do not change the sector's own digest |  |  |

### Physics rectification: greenhouse refits, region runs and tests red on main

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1194 | Runaway greenhouse: the radiative-convective model has no plateau in outgoing flux, so water window opacity, continuum temperature dependence and wings must be reworked until it plateaus |  |  |
| GEN.1202 | The generator's disc gas surface density at 8 kpc is 4.3 solar masses per pc2 against 13 measured, which holds the disc cluster maximum mass too low |  |  |
| GEN.1190 | Replace unsourced rheology constants with the sourced ones: partial-melt factor, Andrade exponent and ratio, ice grain-boundary sliding and the ice melting curves |  |  |
| GEN.1189 | Greenhouse refit method: the refit table is blind to the methane share of a hydrogen air and a dense carbon dioxide air is 2 to 8 percent warm, so find a better fit than one L-BFGS network |  |  |
| GEN.1128 | Region run command: a sphere of a given radius from a start point becomes whole sectors, cut into units of a block of sectors by the galaxy's own (shell, band) arcs, nearest bands first, with `--grow-shells K`, all writing the chart ledger |  |  |
| GEN.1127 | chart_ledger table (tile, run, epoch, time, topped_up_at) with the charted state of a sector derived from it, written by the sector fill units and backfilled from sector_stats.filled_at; the ON DELETE CASCADE removed from position-owned objects (rogue_planets first) so deleting or regenerating a sector removes only objects whose birth tile it is |  |  |
| TEST.498 | A rogue planet is functionally and programmatically identical to a bound planet: the same columns, the same moon and composition generation, the same classes, differing only in the rogue note and having no host star |  |  |
| TEST.569 | Red on main: the reproducibility guard (five call sites), the mass-conversion guard (two tests), the pinned galaxy and the sector panel test (bug) |  |  |
| TEST.580 | Four test groups are red on main: test_gen_cloud_at_generation (3), test_gen_kind_factors (halo versus disk planets), test_burn_down_register (6.170 rows for MAP.275), test_gen_magellanic_classes (supernova factor 0.92 not above 1.0) (bug) |  |  |
| TEST.585 | Three more tests red on main: the accessibility sample_params error for /references/method/<area>, test_health_reports_the_reason_for_maintenance, and test_sweep_atmosphere (bug) |  |  |
| TEST.579 | Remeasure tuning.PREVALENCE_BASE_SHARES: asteroid_belt reads 0.662 on a 600-system draw against the constant 0.039, so test_measured_base_shares_still_hold and every belt prevalence case rest on a stale or wrong base (bug) |  |  |

### Physics rectification: the class-by-domain matrix and stored liquids and air

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1221 | Class-by-domain matrix as data: six domains plus regional, day and night temperature and wet bulb, stored relative humidity and a groundwater flag, and a phase gate per solvent |  |  |
| GEN.1258 | Store an ammonia-water fraction for icy bodies so class X (ammonia-water) can join the class-by-domain fit table | GEN.1221 |  |
| GEN.1253 | toxin_factor has no SO2 or H2S term, so worlds with an acid-eater air read as clean for animals (bug) |  |  |
| GEN.1252 | Methane fluxes of fuel-limited life classes are multiplied by their share of net primary production: sum fuel-limited rows once per fuel (bug) |  |  |
| GEN.1250 | With HABITABLE_WORLD off, about 3 percent of sampled dwarf systems still hold a planet that meets the request's conditions (bug) |  |  |
| GEN.1249 | A sea labelled methane because ice supplies it can store a mix that is mostly ethane: make the label and the mix agree |  |  |
| GEN.1217 | Design: the surface liquid composition model: what is stored (fractions of water, salts, acid, ammonia, CO2, hydrocarbons), how it is derived, and how life changes it |  |  |
| GEN.1211 | A new class for hot thin-air rocks (about 1.6 percent of bodies: 1.4 to 2 g, 340 to 1100 K, 1 to 10 kPa, no ocean) is Class Z of the class code system (Boss confirmed) |  |  |
| GEN.1207 | Best-fit classes, step 3: the classifier replaces the first-match tree, stores class_match and the runner-up letter and match, and the class pages show the percent match |  |  |
| GEN.1281 | Heavy worlds with a plain nitrogen air and no oxygen or methane (about 1.3 percent of bodies) fit no class: widen a class or add one | GEN.1207 |  |

### Physics rectification: citations and the gate's test health

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.99 | Test lock: a citation written outside the registry and references.md fails the gate | DOC.98 | Boss, citations only in the registry. |
| DOC.98 | Citations live only in the central registry and references.md: move them out of every other page, design document and code comment |  | Boss, citations only in the registry. |
| DOC.97 | Citation links in Help > References: every citation on a References page is a link to its entry on the reference list and to its document information |  | Boss, citation links. |
| PERF.241 | A system with every orbit flag forced on takes 15 to 96 seconds on main, so 19 tests time out at 900 s (bug) |  |  |
| TEST.616 | tests/test_sweep_atmosphere.py::test_every_body_of_the_sweep_obeys_every_rule: the mean-molar-mass and air-density rules miss by about 1e-8 relative on 4 bodies (bug) |  |  |
| TEST.614 | tests/test_prevalence.py::test_resolve_keeps_forced_flags_and_the_forcing_rules fails: HABITABLE_WORLD is True but the test expects None (bug) |  |  |
| TEST.606 | The habitability re-score of a wide binary gives a different band1_day_hours than generation stored (18.487 against 18.297, planets row 6), so two re-score tests fail (bug) |  |  |
| TEST.602 | test_fingerprint_tree.py::test_every_saved_sector_has_a_leaf_and_the_stored_nodes_are_the_tree is red on plain main: 24 sectors have a NULL leaf_digest (bug) |  |  |
| VIEW.42 | The System sky route reports built_to_pc 0.0 after a full chart because it does not read the VIEW.25 ledger (bug) |  |  |
| TEST.600 | Census re-tag and prevalence re-baseline after GEN.1275 and TEST.579 land (train 12) |  |  |
| TEST.598 | test_fingerprint (its table list lacks designations and chart_ledger) and test_bring_to_epoch (its epoch counts) are red on plain main (bug) |  |  |
| TEST.591 | test_number_guard (Oxford comma wording) and test_habitability_fill (life_microbial of a rescore differs from generation) are red on plain main (bug) |  |  |

### Physics rectification: best-fit classes, orbit spacing and rogue wiring

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1208 | Best-fit classes, wording: drop the age and stage words from the class tier labels and the "Checked after" and "leftover" wording from the class pages |  |  |
| PERF.242 | Lower or short-circuit REQUEST_PLANET_DRAW_BUDGET for hot O and B stars: a forced habitable world makes the placement loop draw the whole system 8 times and spend the 60-planet request budget | GEN.1282 |  |
| GEN.1282 | Settle orbit spacing before any atmosphere is solved, so a planet pushed outward by space_orbits is built once, not twice (about 60 percent of the cost of a hot O star system) |  |  |
| GEN.1297 | values_from_body does not hand day_share_mean to the latitude bands: generation computes every rocky band for a point star while the re-score and API use the star's real angular size (bug) |  |  |
| GEN.1298 | Life class 8 lives on a solvent it does not allow: test_life_standing::test_no_class_lives_on_a_solvent_it_does_not_allow fails on main (F5V / Ec11: class 8 needs ammonia_water_mass_fraction of 0.1 or more, the body has 0.023) (bug) |  |  |
| GEN.1296 | The volatile-bearing branch (i6) of the lava vapor stage raises "no convergence" from gas_equilibrium.solve for carbon at held fO2 when hot and at low pressure (4000 K at 1 bar, 6000 K at 100 bar) (bug) |  |  |
| GEN.1284 | Generated systems hold less mass than the density profile expects: 78 percent in a solar-circle sample (bug) |  |  |
| GEN.1283 | The habitability re-score gives a different first-band day length than generation stored (14.21 against 13.60 hours on a 'wide' system) (bug) |  |  |
| GEN.1287 | Wire rogue_mass_function into rogue_topup.draw_mass_earth in place of tuning.ROGUE_PLANET_MASS_BINS (RM1, RM2) |  |  |
| GEN.1288 | cluster_unbinding feeds the ejection record: the share of planets a birth cluster has unbound is added to the rogues' ejections |  |  |
| GEN.1289 | Decide rogue_density (RSM1 to RSM4): wire the smoothed star density into the virtual-donor sector draw, or delete it |  |  |
| GEN.1295 | Rogue ages and the rogue stand-ins read the ejection records: the delay (RG2), the closest approach (RK7, MO2, TU3) and the per-star rate (TU1) come from each loss in the architecture, not from drawn distributions |  |  |

### Physics rectification: winds and ring darkening

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1110 | Average wind speeds by day and by night in each of the seven latitude bands: the day-to-night thermal wind from the pressure difference, zonal jets and trade winds from the equator-to-pole temperature difference and the Rossby number, faded out below about 1 kPa |  |  |
| GEN.1055 | Ring age and pollution darkening: a ring stores its age and its ice share and darkness follow from pollution over time (R6); ring shadow climate hook (R8) |  |  |

### Physics rectification: stars, spectra and climate follow-ups

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.998 | The default galaxy disc gas depletion time is 0.85 Gy against the 4.2 Gy the Galaxy measures: check the value and its source (bug) |  |  |
| GEN.995 | Draw the binary-stripped Wolf-Rayet channel so the Galaxy Wolf-Rayet count reaches the observed range (bug) |  |  |
| GEN.962 | Rocky planets hold 1.5 to 2.0 moons of 1e20 kg and over each (the Solar System has 1 for 4 rocky planets): ground the moon rate in giant-impact and capture frequency |  |  |
| GEN.743 | Neutron star equation of state: integrate the tables for radius and structure, and add helium-core and oxygen-neon white dwarfs | GEN.968 |  |
| GEN.968 | Neutron star structure (a): a TOV integration over the equation-of-state tables held on disk, tested against NICER radii |  |  |
| GEN.744 | Alpha-enhanced stellar tracks and an age-metallicity relation from the galaxy's own history |  |  |
| GEN.745 | Stellar spectra by temperature, gravity and activity (model atmospheres) instead of a blackbody with a seeded UV share | GEN.970 |  |
| GEN.970 | Stellar spectra (a): a spectral library entry (BT-Settl, MUSCLES for M dwarfs) and a lookup by temperature, gravity and activity, with tests on the Sun and an M dwarf |  |  |
| GEN.752 | Climate feedback on the surface temperature fixed point: ice-albedo, water vapour and cloud albedo | GEN.974 |  |
| GEN.974 | Climate feedback 3: cloud albedo on the fixed point (cloud cover changes the reflected share) |  |  |
| GEN.753 | Validation sweep for each physics area: compare generated populations with published statistics | GEN.503, GEN.504 |  |

### Physics rectification: older physics follow-ups, units and the storage rows

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.820 | Residuals of the bulk iron branch table (GEN.512): the coreless branch, a slope other than 1, the 0.05 dex scatter stand-in and a larger tuning sample |  |  |
| PERF.198 | Speed up util/detmath.matmul and cache equal gas_equilibrium.solve calls (finding 6) |  |  |
| API.33 | API field descriptions, the help pages and the CLI output: every number states its unit (scan) (bug) |  |  |
| TEST.389 | test_golden_seed one-versus-four-workers reports a sector-structure difference after epoch 101 (GEN.512); unconfirmed, the clean-main run timed out at 900 s (bug) |  |  |
| GEN.835 | The tidal heat flux (tidal_flux_w_m2) is wired into the ground draw's sulfur rule |  |  |
| TEST.550 | Five census tests stay red on main until the census re-tag for epochs 129 and 130 (bug) |  |  |
| DB.101 | Storage and units for every physics rectification output: a column, a column_units entry and a migration shipped with its item |  |  |
| GEN.804 | Atmosphere escape and retention read the stored XUV exposure, flare irradiation and magnetopause (G2, step 1) |  | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
| GEN.805 | The surface dose and the habitability radiation component read the stored flare irradiation index and magnetopause (G2, step 2) | GEN.804 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
| GEN.806 | Delete the duplicate XUV recomputation in physics/mixture_subject.py and read xuv_exposure_index (G2, step 3) | GEN.804 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
| GEN.807 | Stellar activity (flare rate, XUV saturation) follows the star's rotation: rotation_period_hours, flare_alpha and xuv_saturated get readers (G2, step 4) |  | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
| GEN.808 | The ozone-loss flag feeds the UV limit on life and the surface look (G2, step 5) | GEN.805 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
| GEN.801 | Wire the giant-planet interior outputs into the field, the heat and the description (G3) |  | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each. |
| GEN.944 | Later: an age-metallicity relation inside each population and a 10 to 13 Gy thick-disk window; variable-star flags (Cepheids, Mira and long-period variables, carbon and S stars) |  |  |
| GEN.949 | Wire the held methane absorption tables (AD1 Karkoschka k-tables, PDS, NMSU) into the visible colour of giant planets |  |  |
| GEN.993 | Widen the core light-element total in bulk_draw to 0.02 to 0.25 (0.25 dex) instead of the uniform 0.03 to 0.12 |  |  |
| GEN.994 | Wire the hot-giant regime edges (equilibrium_edges) into the atmosphere row for hot giants |  |  |
| GEN.91 | Classes like S and V in the hot and cold zones |  |  |

### Temperature coverage, nebula composition and belt rules (rows filed in the 8.2 file by mistake, moved here)

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1236 | Two class N planets sit at the 2000 K surface_air cap but keep the uncapped greenhouse factor, so stored temperature and greenhouse disagree (bug) |  |  |
| GEN.1238 | A Venus twin (Venus mass at 0.72 AU, G2V) reaches 3,100 bar and 1,279 K against Venus's 92 bar and 737 K, with a tenth of the CO2 mass the pressure implies (bug) |  |  |
| GEN.1240 | The greenhouse table is held at 1000 bar, so steam airs of 1e4 to 1e5 bar show 700 to 1300 K, too cool (bug) |  |  |
| GEN.1247 | Outer asteroid belts store about 4.5e25 kg (7 Earth masses) across 11 to 74 AU at 6 to 9 Gy, far above a Kuiper-belt-like 0.03 Earth masses, while inner belts hold 1e20 to 1e24 kg: check the belt mass law at large distance and old age (bug) |  |  |
| GEN.1259 | A planetless stretch becomes a belt only if enough leftover rocky mass could survive there: check the belt rule of GEN.1068 and GEN.1156 against the physics and the observed belt rates (about two thirds of systems have a belt now, from 4 percent before train 7) |  |  |
| GEN.1260 | The latitude band table's Day/Night and Longest day/night rows are empty or wrong for slow spinners: rewrite them as a time-stepped spell calculation (bug) |  |  |
| GEN.1261 | Cold belts are 3 to 4 times too common: add a physical loss step for solids that never formed planetesimals, and cap a gap belt by the isolation mass a planet would have taken | GEN.1259 |  |
| GEN.1262 | Temperature coverage part 1: the shared illumination series (physics/surface_cycle.py): cosine of the zenith angle per latitude-longitude cell and time step, with GEN.1260's spell search moved onto it | GEN.1260 |  |
| GEN.1263 | Temperature coverage part 2: the temperature series: one energy balance per cell, anchored so the area-and-time mean equals the stored mean | GEN.1262, GEN.1039, GEN.972 |  |
| GEN.1264 | Temperature coverage part 3: what is read off the series: Blue, Green, cold and hot shares, livable time, coldest and hottest hour, best site, longest spell out of the livable band | GEN.1263 |  |
| GEN.1265 | Temperature coverage part 4 (lock): the temperature score from the coverage share, the new stored columns, and the migration | GEN.1264 |  |
| GEN.1266 | Temperature coverage part 5 (lock): gear from the best site with a worst-conditions note, hours outside from the best site with a worst-to-best range, the text and the design doc | GEN.1265 |  |
| GEN.1267 | A generated nebula's own shape follows the same model-density contour the map now draws: bounds, containment and the stars placed inside it |  |  |
| GEN.1268 | Nebulae carry a real composition: gas and dust fractions by kind, metals from the place and age, ejecta from the progenitor |  |  |
| GEN.1271 | Nebula composition part c: supernova remnants as swept-up gas plus ejecta, checked against Cassiopeia A |  |  |
| GEN.1272 | Nebula composition part d: the nebula page, the describe sheet and the Properties sheet show the composition with units and derived marks | GEN.1276, GEN.1277 |  |
| GEN.367 | The Magellanic Clouds' gas is not in `kind_factors` and `nebula_field`: the Clouds hold stars but no nebulae |  | Finding from Core 3. |
| GEN.1273 | Nebula composition part e: the Magellanic Clouds' gas metallicity as a place value | GEN.367 |  |
| GEN.1275 | The habitable_world option (and its prevalence setting) does not mean what the prevalence count measures: it bars only Earth-like ecosphere worlds, the count includes every body of a life-carrying class (bug) |  |  |
| GEN.1276 | Nebula composition part g: store the composition of every nebula and supernova remnant (new columns, progenitor mass, epoch, migration, backfill, units and wiring register) | GEN.1271 |  |
| GEN.1277 | Nebula composition part h: store the composition of the field clouds (nebula_field), including the progenitor mass where a cloud has one | GEN.1276 |  |
| PERF.236 | The prevalence draws with a forced habitable_world are too slow to finish: 600 systems did not finish in 5,000 s (about 0.9 s a system at baseline, several times that forced), and asteroid_belt +30 did not finish in 6,600 s (bug) |  |  |
| PERF.238 | Temperature coverage cost: the illumination and temperature series must stay about 20 ms a rocky body in bulk generation and in the bring-to-epoch recompute | GEN.1263 |  |
| PERF.240 | A system that bars habitable worlds costs about 12 s to draw against 3 s, because the late remakes and planet redraws that remove life-class bodies run again and again (bug) | GEN.1275 |  |
| TEST.582 | Temperature coverage sweep: set the Blue cut and the Green and Yellow edges on a sweep that includes cold-mean worlds with warm regions and hot-mean worlds with cool regions, then re-baseline the TEST.430 prevalence census | GEN.1265 |  |

## Phase 12: The headline closes

The last audits and the headline: they test the component pass, the core and the halo, so they come after phases 4 and 5.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### The headline closes

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.261 | Release gate: a full bug hunt audit of the whole code base and its behaviour, done before every major or minor version change (8.1 first) |  | Standing release gate, Boss 2026-10-21 5:10 PM EDT; also gates every later major or minor version change, not build numbers (Boss, 8:01 PM EDT 2026-10-22). |
| TEST.380 | First live test pass of species, civilizations, polities and territories on a generated galaxy, before the 8.1 stamp |  |  |
| SEC.33 | Release gate: a full security code audit of the whole code base, done before every major or minor version change (8.1 first) |  | Standing release gate, Boss 2026-10-21 5:10 PM EDT; also gates every later major or minor version change, not build numbers (Boss, 8:01 PM EDT 2026-10-22). |
| GEN.301 | Headline of 8.1: the new generate scatters, accurate galactic generation and scattering including the halo | SEC.33, TEST.261, TEST.300, GEN.743, GEN.744, GEN.745, GEN.752, GEN.753 |  |
