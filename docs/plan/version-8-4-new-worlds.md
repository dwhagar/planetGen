# Version 8.4: New worlds

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.4

## Goal

Phases 1 to 4 hold the generation and storage half of the composition of every body and of the atmosphere colour: the bulk composition that drives mass, density and gravity, the surface and small-body composition, the one atmosphere mixture, the stored albedo and greenhouse values, ozone, haze, dust and cloud, and the consistency proofs (Boss, 2026-10-16, 11:25 PM EDT, 11:49 PM EDT and 11:58 PM EDT). Their research, designs and data tables are in 8.1; the headline that completes them is 9.0. Phase 5: the planet class refactor around the habitability index (classes R to Z), nebula conditions on planets, asteroid fields as object systems and the habitability symbols on the pages. Proposed by the TODO thread; not yet approved by Boss.

## Phase 1: Bulk composition and the albedo records

The composition fractions that drive density, mass and gravity, stored and described, and the recording of the albedo and greenhouse values. The research, the designs and the data tables they read (the gas, materials, element and stellar tables) are in 8.1 phases 1 and 2 (Boss, 2026-10-16, 11:58 PM EDT: the building blocks of the composition headline go down to 8.1).

### Bulk composition

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.458 | Each system draws its snow line, 0.2 dex around the 170 K blackbody position, and stores it once; the fixed 2.7 AU constant goes |  | Body composition research lane finding F8, cluster 3; generator epoch bump. Core 4 proposed. |
| GEN.428 | Planets and moons draw their bulk composition fractions instead of carrying a fixed text | GEN.427, GEN.440, GEN.458 | Boss 11:49 PM EDT; generator epoch bump. Core 4 proposed. |
| GEN.429 | Bulk density, mass, radius and gravity are derived from the composition fractions with the mass-radius relation | GEN.428 | Boss 11:49 PM EDT; the atmosphere chain reads this gravity. Core 4 proposed. |
| GEN.459 | One radiogenic module: heat from the isotope-wise decay of potassium, thorium, uranium and aluminium-26, drawn from the host star's abundances, replaces the bulk-silicate values that are 37 to 50 percent too high and the separate ground-dose draw (bug) | GEN.428, GEN.438 | Body composition research lane findings F1 and F2, cluster 4; generator epoch bump. Core 4 proposed. |
| GEN.460 | Giant planets and brown dwarfs get their radius and heavy-element content from the checked laws: the heavy-element mass law, the Fortney and Baraffe tables with inflation, and the burning limits shifted by metallicity and helium (bug) | GEN.427, GEN.428 | Body composition research lane finding F7, cluster 6 (parts A and C); generator epoch bump. Core 4 proposed. |
| DB.69 | Planets and moons store the composition fractions | GEN.429, GEN.458, GEN.459, GEN.460, GEN.498 | Boss 11:49 PM EDT; migration revision after DB.67. Core 4 proposed. |
| GEN.430 | The planet and moon descriptions state the composition from the stored fractions | GEN.428, GEN.429, GEN.221 | Boss 11:49 PM EDT; the description side. Core 4 proposed. |

### Albedo and greenhouse records

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.412 | Generation records the albedo and the greenhouse multiplier it drew, and the temperature is exactly what they give |  | Boss 11:24 PM EDT; the inverse feeds the migration. Core 4 proposed (Boss, 11:41 PM EDT). |

### Atmosphere tables

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.478 | The atmosphere possibility tables: the 16 categories and 104 rows of part 09, loaded through the gated table loader, with every row that has no sourced weight flagged and drawn from its bracket | GEN.410 | Atmosphere lane part 09: the 16 weighted atmosphere tables (104 rows) as the data of GEN.472. Core 4 proposed. |

## Phase 2: Surfaces, small bodies and one mixture

The ground's composition and the small bodies' composition from the bulk model, and the atmosphere's molar density derived from its partial pressures with the values that follow from it.

### Surface and small bodies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.461 | The hydrosphere reads the stored water mass fraction: the compressed water-layer depth, the base pressure, the high-pressure ice and the liquidus are recomputed, and the class label follows from the water share (bug) | GEN.428, GEN.429, GEN.458 | Body composition research lane finding F4, cluster 8 (parts C, D, E); generator epoch bump. Core 4 proposed. |
| GEN.441 | Surface and crust composition: the minerals, ices and ocean on the ground, from the bulk fractions, the temperature and the air | GEN.428, GEN.429, GEN.459, GEN.461 | Boss 11:58 PM EDT; the ground's composition; generator epoch bump. |
| GEN.495 | The comet population table: species ranges and log-normal scatters replace the draw centred on one comet, with a carbon-chain-depleted flag | GEN.457 | Body lane G-1 (2:43 AM EDT as the lane reported): the comet population table replaces the 67P-centred draw. Core 4 proposed. |
| GEN.442 | Asteroids, comets, rings and dwarf planets draw their composition from the formation distance and the snow line | GEN.428, GEN.429, GEN.458, GEN.459, GEN.470, GEN.495 | Boss 11:58 PM EDT; small bodies; generator epoch bump. |
| GEN.466 | Moons draw their origin type, ice fraction and composition from the parent, the circumplanetary ice line and the age | GEN.428, GEN.429, GEN.442, GEN.458 | Body composition research lane part D sections 2.8 to 2.11; generator epoch bump. Core 4 proposed. |
| GEN.455 | Airless bodies and comets draw an exosphere or coma from their surface composition, temperature and the star's wind: the sodium and potassium tails and the production rates | GEN.441, GEN.442, GEN.410, GEN.478 | Atmosphere research lane new work 6 (parts 04 section 5 and 05); generator epoch bump. Core 4 proposed. |

### Mixture and stored values

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.411 | The atmosphere's molar density comes from its per-gas partial pressures, so mixture, pressure and density always agree | GEN.410, GEN.429, GEN.478 | Boss 11:24 PM EDT; generator epoch bump. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.479 | The class atmosphere mixes contradict the checked chemistry: classes C, D and S have no mix, free oxygen sits in seven classes with no life rule, class L is 60 percent argon, E and F hold ammonia and methane at 340 to 400 K, N reads reducing and is oxidised, and I, J and T share one mix in every zone (bug) | GEN.478, GEN.411 | Atmosphere lane part 09 conflicts 3 to 8 with the code, checked on main; replaced by the recipe tables. Core 4 proposed. |
| GEN.452 | The condensation table is replaced by the checked saturation curves: the hydrogen sulphide row is a copy of the sulphur dioxide row, and carbon dioxide, sulphur dioxide and the cold gases are 2 to 20 K off (bug) | GEN.410, GEN.470 | Atmosphere research lane finding 1 (part 01 section 7); generator epoch bump. Core 4 proposed. |
| GEN.451 | The mantle redox mean follows the meteorite building blocks: a young rocky planet starts reduced, not at the modern Earth's 3.5 above the iron-wuestite buffer | GEN.410, GEN.428, GEN.478, GEN.498 | Atmosphere research lane finding, part 5, 2026-10-17; generator epoch bump. Core 4 proposed. |
| GEN.454 | Giant planets, ice giants and brown dwarfs draw their atmosphere mixture from the equilibrium chemistry sequence: by temperature, pressure, metallicity and quench | GEN.410, GEN.411, GEN.460, GEN.478 | Atmosphere research lane new work 5 (part 05 sections 8 and 9); generator epoch bump. Core 4 proposed. |
| GEN.422 | Stop the recalculation redraw and order the steps so pressure, mixture and temperature agree after the cold-side freeze-out | GEN.411, GEN.412, GEN.452 | Composition thread finding; Boss 11:24 PM EDT chain. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.423 | Atmosphere density, scale height and the greenhouse multiplier follow from the mixture | GEN.410, GEN.411, GEN.412 | Composition thread finding; Boss 11:24 PM EDT chain. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.424 | The balance gas is a named gas with its own molar mass, and habitability reads and writes the mean molar mass of the stored mixture | GEN.410, GEN.411 | Composition thread finding; Boss 11:24 PM EDT chain. Core 4 proposed (Boss, 11:41 PM EDT). |
| DB.67 | Planets and moons store the albedo, the greenhouse multiplier and the molar density they were built with | GEN.411, GEN.412, GEN.422, GEN.423, GEN.424 | Boss 11:24 PM EDT; migration revision after Core 1's 0102 to 0105 stack. Core 4 proposed (Boss, 11:41 PM EDT). |
| DB.71 | Planets and moons store the partial pressures of the added gases, not only the ten stored today | GEN.411, GEN.410, DB.67 | Atmosphere research lane finding 10 (part 05 section 13); migration revision after DB.67. Core 4 proposed. |

## Phase 3: Ozone, haze, dust and cloud, and the stored shares

The ozone figures and the scatterers that the colour needs, each its own thread, and the columns for the surface and small-body shares.

### Ozone

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.413 | Ozone in the mixture: the column and partial pressure from the oxygen, the star's UV and the Chapman balance | GEN.411, GEN.478 | Boss 11:24 PM EDT; real research, generated where missing. Core 4 proposed (Boss, 11:41 PM EDT). |

### Haze, dust and cloud

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.414 | Haze, dust and cloud: derived or drawn optical depths for every planet and moon | GEN.411, GEN.412, GEN.478 | Boss 11:24 PM EDT; draws and derivations, nothing fixed per class. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.499 | The lava-world vapour species and mean molar mass from temperature and composition, replacing the heavy-oxide range of 40 to 60 g/mol | GEN.410, GEN.427 | Body lane G-5 (part G3): lava-world vapour species and mean molar mass 25 to 33 g/mol. Core 4 proposed. |
| GEN.456 | Lava worlds and hot rocky planets draw a vapour atmosphere over the magma from the surface composition and temperature: sodium, potassium, silicon monoxide, magnesium, iron and the oxides | GEN.441, GEN.410, GEN.411, GEN.414, GEN.451, GEN.478, GEN.499 | Atmosphere research lane new work 8 (parts 05 and 06); generator epoch bump. Core 4 proposed. |

### Stored surface composition

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.70 | Planets, moons and small bodies store the surface and small-body composition shares | GEN.441, GEN.442, GEN.461, GEN.462, GEN.466 | Boss 11:58 PM EDT; migration revision, in merge order. |
| GEN.477 | The asteroid belt and field composition is drawn from the star, the distance and the belt's mass, not four names drawn uniformly from a list of 84 (bug) | GEN.442, GEN.457, GEN.458, DB.70 | Body composition research lane finding F11 (tables E1 to E4): belt and field composition drawn from the star, distance and belt mass; generator epoch bump. Core 4 proposed. |
| GEN.471 | Clathrate hydrates and eutectic brines for icy crusts and ocean moons: the dissociation pressures and the liquid range below 273 K | GEN.461, GEN.427 | Body composition research lane part F item F4 (later); generator epoch bump. Core 4 proposed. |
| GEN.462 | Ground albedo by band: the area-weighted colour of the ground, ice and ocean cover with space weathering, handed to the atmosphere colour once | GEN.441, GEN.442 | Body composition research lane cluster 9 (part E sections 2.6 and 2.7, part D taxon albedo); generator epoch bump. Core 4 proposed. |

## Phase 4: Stored scatterers, the text and the proofs

The columns for ozone and the aerosols, the atmosphere text written from the mixture, the consistency sweeps and the docs.

### Consistency

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.68 | Planets and moons store the ozone column and the haze, dust and cloud optical depths | GEN.413, GEN.414, DB.67 | Boss 11:24 PM EDT; migration only. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.415 | Every compound the planet descriptions mention comes from the gas table, and the atmosphere text is written from the stored mixture | GEN.410, GEN.411, GEN.413, GEN.221 | Boss 11:24 PM EDT; the description side of consistency. Core 4 proposed (Boss, 11:41 PM EDT). |
| TEST.209 | One consistency sweep over every planet and moon class: composition, density, temperature, ozone, aerosols and text all agree | GEN.411, GEN.412, GEN.413, GEN.414, GEN.415, DB.67, DB.68, GEN.422, GEN.423, GEN.424 | Boss 11:24 PM EDT; the proof. |
| TEST.210 | A bulk composition sweep over every class: fractions, density, mass, radius and gravity agree | GEN.428, GEN.429, GEN.430, DB.69, GEN.441, GEN.442, DB.70, GEN.458, GEN.459, GEN.460, GEN.461, GEN.462, GEN.466 | Boss 11:49 PM EDT; the proof. |
| DOC.63 | Docs: the atmosphere data model, the gas table and how the values agree | TEST.209 | Boss 11:24 PM EDT; docs only. |
| DOC.65 | Docs: bulk composition, the materials table and how density, mass and gravity follow | TEST.210 | Boss 11:49 PM EDT; docs only. |

## Phase 5: New worlds

Planet classes, nebula conditions, asteroid fields and habitability on the pages.

### Planet classes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.33 | One class per PR, each with its tests |  | Moved to phase 2 under GEN.90, after the habitability score. One class per PR (R and S first). Built on fixed physics so new classes aren't tuned to wrong masses, moons or zones. |
| GEN.146 | Teff-dependent habitable zone from the Kopparapu table |  | Research: in GEN.91's dependency chain. |
| GEN.355 | Atmosphere retention module with its test table | GEN.33, GEN.410 | Split from GEN.91. |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | GEN.33, GEN.146 | Under GEN.90. Same reconcile/zone code as phase 0's physics fixes. |
| GEN.350 | Class R, sub-Neptune, with the envelope-survival rule, the Hycean flag and the retirement of class T | GEN.33, GEN.355, GEN.146, GEN.428, GEN.411 | Split from GEN.28. |
| GEN.351 | Class U, icy world | GEN.33, GEN.428, GEN.411 | Split from GEN.28. |
| GEN.352 | Classes W and X, subsurface ocean bodies | GEN.33, GEN.428, GEN.411 | Split from GEN.28. |
| GEN.353 | Class Y, Titan-like world | GEN.355, GEN.146, GEN.428, GEN.411 | Split from GEN.28. |
| GEN.92 | Life and its highest stage follow the habitability score | GEN.350, GEN.351, GEN.352, GEN.353 |  |
| GEN.354 | Class Z, lifeless temperate world | GEN.92, GEN.146, GEN.428, GEN.411 | Split from GEN.28. |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | GEN.350, GEN.351, GEN.352, GEN.353, GEN.354 | Under GEN.90; class Z is Boss's Earth-size world that never had life (2026-10-03). Class S landed with GEN.38 (PR #415); the other six classes here. Closes with GEN.33's PRs; PLANET_CLASSES in program_constants.py. |
| GEN.356 | Class N re-scoped as the hot, lifeless world; its life chemical and hidden timeline removed | GEN.355, GEN.92, GEN.428, GEN.411 | Split from GEN.91. |
| GEN.91 | Classes like S and V in the hot and cold zones | GEN.355, GEN.356 |  |
| GEN.145 | Class S atmosphere rule: S keeps air unless the shoreline ratio is over 30 | GEN.355, GEN.92 | Research follow-up to class S (built). |
| GEN.29 | Sweep every planet class for sense once the new ones are in (bug) | GEN.28, GEN.27, GEN.91, GEN.92, GEN.146, TEST.209, TEST.210 | Under GEN.90: the refactor is the sweep. Bug, but by definition a sweep after the new classes; it can't go earlier. Includes rocky rogues of 10-16 Earth masses (up to 17,600 km) that get S as nearest fit (PR #415). |
| GEN.90 | Refactor the planet classes around the habitability index | GEN.33, GEN.28, GEN.27, GEN.29, GEN.91, GEN.92 | Parent; takes GEN.33, GEN.28, GEN.27 and GEN.29 as its subitems. |

### Nebula conditions on planets

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.358 | GEN.94's rules and the four exotic-environment additions reach planet formation | GEN.357 | Split from GEN.95. |
| GEN.95 | Nebula conditions applied when planets and surfaces are generated | GEN.357, GEN.358 |  |
| GEN.93 | Nebula conditions in planet generation | GEN.95 | Parent. |

### Habitability on the pages

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.143 | Color sectors by their number of habitable locations |  | GitHub issue [#717](https://github.com/dwhagar/planetGen/issues/717). |
| UX.101 | The PHI-6 symbols as vector shapes in the icon library, both themes, coloured by rating |  |  |
| UX.100 | Planet and moon lists show habitability as symbols, not notes | UX.101 |  |

### Asteroid fields

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.360 | Asteroid fields and belts render from their uid: one seeded generator in a shared asteroidfield.js for dust, rocks and heroes |  | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
| GEN.361 | The asteroid field page has a 3D view | GEN.360 | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
| GEN.362 | Rocks are instanced meshes in chunks with three levels of detail, within a frame budget | GEN.361 | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
| GEN.363 | Belt bodies move on the GPU: the vertex shader solves Kepler from per-instance elements | GEN.362 | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
