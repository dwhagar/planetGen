# Version 8.4: New worlds

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.4

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

The composition suite (the bulk composition of every body, surfaces and small bodies, the one atmosphere mixture, ozone, haze, dust and cloud, and the consistency proofs) moved into 8.1 phases 6 to 9 on Boss's word (5:48 PM EDT 2026-10-23: build it now, testing mode). The habitability pages moved into 8.1 phase 5 (5:38 PM EDT). What stays here: the planet class refactor around the habitability index (classes R to Z), nebula conditions on planets, asteroid fields as object systems and Proposed by the TODO thread; not yet approved by Boss.

## Phase 1: New worlds

Planet classes, nebula conditions and asteroid fields (the habitability pages moved into 8.1 by Boss, 5:38 PM EDT 2026-10-23).

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### Planet classes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.33 | One class per PR, each with its tests |  | Moved to phase 2 under GEN.90, after the habitability score. One class per PR (R and S first). Built on fixed physics so new classes aren't tuned to wrong masses, moons or zones. |
| GEN.355 | Atmosphere retention module with its test table | GEN.33 | Split from GEN.91. |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | GEN.33 | Under GEN.90. Same reconcile/zone code as phase 0's physics fixes. |
| GEN.350 | Class R, sub-Neptune, with the envelope-survival rule and the Hycean flag (class T stays: it is the sub-Neptune row of the classifier) | GEN.33, GEN.355 | Split from GEN.28. |
| GEN.351 | Class U, icy world | GEN.33 | Split from GEN.28. |
| GEN.352 | Classes W and X, subsurface ocean bodies | GEN.33 | Split from GEN.28. |
| GEN.353 | Class Y, Titan-like world | GEN.355 | Split from GEN.28. |
| GEN.354 | Class Z, lifeless temperate world |  | Split from GEN.28. |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | GEN.350, GEN.351, GEN.352, GEN.353, GEN.354 | Under GEN.90; class Z is Boss's Earth-size world that never had life (2026-10-03). Class S landed with GEN.38 (PR #415); the other six classes here. Closes with GEN.33's PRs; PLANET_CLASSES in program_constants.py. |
| GEN.356 | Class N re-scoped as the hot, lifeless world; its life chemical and hidden timeline removed | GEN.355 | Split from GEN.91. |
| GEN.91 | Classes like S and V in the hot and cold zones | GEN.355, GEN.356 |  |
| GEN.145 | Class S atmosphere rule: S keeps air unless the shoreline ratio is over 30 | GEN.355 | Research follow-up to class S (built). |
| GEN.29 | Sweep every planet class for sense once the new ones are in (bug) | GEN.28, GEN.27, GEN.91 | Under GEN.90: the refactor is the sweep. Bug, but by definition a sweep after the new classes; it can't go earlier. Includes rocky rogues of 10-16 Earth masses (up to 17,600 km) that get S as nearest fit (PR #415). |
| GEN.90 | Refactor the planet classes around the habitability index | GEN.33, GEN.28, GEN.27, GEN.29, GEN.91 | Parent; takes GEN.33, GEN.28, GEN.27 and GEN.29 as its subitems. |

### Nebula conditions on planets

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.358 | GEN.94's rules and the four exotic-environment additions reach planet formation |  | Split from GEN.95. |
| GEN.95 | Nebula conditions applied when planets and surfaces are generated | GEN.358 |  |
| GEN.93 | Nebula conditions in planet generation | GEN.95 | Parent. |
| GEN.775 | Replace the planet-minus-star Fe/Mg scatter and super-Mercury tail (registry row BC8) with a larger sample of rocky planets with precise masses and stellar Fe/Mg |  |  |
| GEN.776 | Replace the rogue-planet volatile loss at ejection (registry row BC11) with a measured atmosphere or free-floating-planet volatile survey |  | Aux 14, GEN.734: replacing item of a registry row; pool, after 8.1. |
| GEN.777 | Replace the core light-element split (registry row BC18) with a second body that has seismic or core-density data |  | Aux 14, GEN.734: replacing item of a registry row; pool, after 8.1. |
| GEN.794 | Replace the NH3 photolysis activity factor 1.078 (GEN.501, registry stand-in 1.05 to 1.10) with the Huebner and Mukherjee 2015 rates or a computation from the Cheng 2006 cross section |  | Citations lane, GEN.793: pool, after 8.1. |

### Asteroid fields

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.360 | Asteroid fields and belts render from their uid: one seeded generator in a shared asteroidfield.js for dust, rocks and heroes |  | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
| GEN.361 | The asteroid field page has a 3D view | GEN.360 | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
| GEN.362 | Rocks are instanced meshes in chunks with three levels of detail, within a frame budget | GEN.361 | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
| GEN.363 | Belt bodies move on the GPU: the vertex shader solves Kepler from per-instance elements | GEN.362 | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
| MAP.274 | Land colour from the bulk composition: oxidation and rock type |  | Body rendering R5; Boss's plan of 10:26 PM EDT 2026-10-23. |
| MAP.275 | Gas and ice giant painter: bands, storms and thermal glow from rotation, temperature and the atmosphere colour, rings with resonance gaps from the moons | GEN.670, GEN.419 | Body rendering R6a; Boss's plan of 10:26 PM EDT 2026-10-23. |
| MAP.279 | Small bodies: lumpy meshes for moons and asteroids under about 400 km, from the stored shape |  | Body rendering R13; Boss's plan of 10:26 PM EDT 2026-10-23. |
| DB.99 | Store the ground colour and the ice fraction on planets and moons | MAP.274 | Body rendering R2a; Boss's plan of 10:26 PM EDT 2026-10-23. |

### Wiring stored outputs

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.804 | Atmosphere escape and retention read the stored XUV exposure, flare irradiation and magnetopause (G2, step 1) |  | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
| GEN.805 | The surface dose and the habitability radiation component read the stored flare irradiation index and magnetopause (G2, step 2) | GEN.804 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
| GEN.806 | Delete the duplicate XUV recomputation in physics/mixture_subject.py and read xuv_exposure_index (G2, step 3) | GEN.804 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
| GEN.807 | Stellar activity (flare rate, XUV saturation) follows the star's rotation: rotation_period_hours, flare_alpha and xuv_saturated get readers (G2, step 4) |  | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
| GEN.808 | The ozone-loss flag feeds the UV limit on life and the surface look (G2, step 5) | GEN.805 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.800, Boss 10:35 PM EDT. |
| GEN.800 | Wire the stored magnetosphere and stellar activity outputs into escape, dose and habitability (G2) | GEN.804, GEN.805, GEN.806, GEN.807, GEN.808 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each. |
| GEN.809 | A giant's dipole moment and class come from the metallic-hydrogen shell and the dilute core (G3, step 1) | GEN.733 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.801, Boss 10:35 PM EDT. |
| GEN.810 | A giant's internal heat reads the interior layers and the helium shell (G3, step 2) | GEN.733 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.801, Boss 10:35 PM EDT. |
| GEN.811 | The description of a giant says what its interior is (G3, step 3) | GEN.809, GEN.810 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain GEN.801, Boss 10:35 PM EDT. |
| GEN.801 | Wire the giant-planet interior outputs into the field, the heat and the description (G3) | GEN.733, GEN.809, GEN.810, GEN.811 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each. |

## Phase 2: Settlements and nation-states

Nation-states, the capital and ten largest cities, and the station and space port scatter for every planet with a technological civilization at a sufficient tech level, all fixed until 9.2 (Boss, 9:27 PM EDT 2026-10-24).

### Settlements

| ID | Item | Needs | Note |
|---|---|---|---|
| POP.17 | Which worlds get settlements: the tech gate and the planet inputs (habitable land area, carrying capacity) for a species whose civilization is dated | DB.103, GEN.774, POP.30 | Boss, 9:27 PM EDT 2026-10-24: settlements tree; fixed until 9.2. |
| POP.18 | Nation-states on a planet: how many, their land shares, names, governments, power and which are dominant | POP.17 |  |
| POP.19 | Cities: the capital and the ten largest cities of a planet, their populations, names and positions on land | POP.17 |  |
| POP.20 | Cities belong to nation-states: a probability draw weighted by city size and nation dominance, with the capital in the dominant nation | POP.18, POP.19 |  |
| POP.21 | Station and space port scatter: orbital stations, space ports and starbases placed for each nation and polity from the tech gate, as stored facilities | POP.20 | Carves the station part out of POP.11. |
| POP.22 | The settlements stage: its own command, its place in the galaxy run and the rescan rules | POP.21 |  |
| POP.23 | Pages and API for nations, cities and ports: the planet page sections, nation and city pages, species and polity links | POP.22 |  |
| TEST.381 | First live test pass of nations, cities and ports on a generated galaxy, before the 8.4 stamp | POP.23, POP.27 |  |
| POP.40 | Species traits take radiation tolerance, pressure and air tolerance, thermal range and day length from the homeworld (G1, step 1) | POP.17 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain POP.38, Boss 10:35 PM EDT. |
| POP.41 | Facility ratings (health, housing, maintenance, resources) read the host's dose, pressure, air and temperature (G1, step 2) |  | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain POP.38, Boss 10:35 PM EDT. |
| POP.42 | Facility need and upkeep read the equipment kit the host requires (eq_body_*, eq_rad_suit*, eq_gas_filter, eq_co2_scrub, eq_dust, eq_supply*, eq_shelter*, outside_hours) (G1, step 3) | POP.41 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain POP.38, Boss 10:35 PM EDT. |
| POP.43 | Settlement and facility pages say why: the species traits and the host's dose, air and kit behind each rating (G1, step 4) | POP.40, POP.42 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain POP.38, Boss 10:35 PM EDT. |
| POP.38 | Species traits and facility needs read the homeworld and host: dose, pressure, air, thermal range, day length and the equipment kit (G1) | POP.17, POP.40, POP.41, POP.42, POP.43 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each. |
| POP.44 | Register the facility constants (_WEIGHTS, FACILITY_HOST_RICHNESS, FACILITY_KIND_SIZE, FACILITY_RATING_*) in the GEN.725 registry with tier and tweak level (G6, step 1) |  | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain POP.39, Boss 10:35 PM EDT. |
| POP.45 | A facility's richness reads the host's stored composition and resource scores instead of a per-type constant (G6, step 2) | POP.44, GEN.799 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain POP.39, Boss 10:35 PM EDT. |
| POP.39 | Facility constants are registered with a tier and tweak level, and a facility's richness reads the host's stored composition (G6) | GEN.799, POP.44, POP.45 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each. |

### Resource pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.223 | Resource pages: belts, fields, giants and icy bodies show their resources, and lists, search and the API filter and sort by resource and score | GEN.799 | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |
| TEST.386 | First live test pass of the resource locations on a generated galaxy, before the 8.4 stamp | MAP.312, UX.223 | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |
