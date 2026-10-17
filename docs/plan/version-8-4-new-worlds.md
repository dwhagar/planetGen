# Version 8.4: New worlds

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.4

## Goal

Phases 1 to 4 hold the foundations of the 9.0 atmosphere colour (the sources and design, the one gas table and mixture, the stored albedo and greenhouse values, ozone, haze, dust and cloud, and the consistency proof; Boss, 2026-10-16, 11:25 PM EDT). Phase 5: the planet class refactor around the habitability index (classes R to Z), nebula conditions on planets, asteroid fields as object systems and the habitability symbols on the pages. Proposed by the TODO thread; not yet approved by Boss.

## Phase 1: Research and design for the atmosphere

The base of the 9.0 atmosphere colour (Boss, 2026-10-16, 11:25 PM EDT: foundations back from 9.0): the verified sources and the checked design.

### Research and design

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.408 | Research: verify the atmosphere sources and say what is missing |  | Boss 11:24 PM EDT; a read-only research thread, any lane; no code. |
| GEN.409 | The atmosphere colour design: the checked plan for the gas table, the mixture, the stored values, ozone, aerosols and the colour model | GEN.408 | Boss 11:24 PM EDT; design note, no code. Core 4 proposed (Boss, 11:41 PM EDT). |

## Phase 2: One table, one mixture, stored values

The one gas-and-compound table, the molar density derived from the partial pressures, and the albedo and greenhouse values recorded and stored.

### Table, mixture and stored values

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.410 | One table of the gases and compounds that can be in an atmosphere, with every compound the descriptions mention | GEN.409 | Boss 11:24 PM EDT; the one table the rest reads. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.411 | The atmosphere's molar density comes from its per-gas partial pressures, so mixture, pressure and density always agree | GEN.410 | Boss 11:24 PM EDT; generator epoch bump. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.412 | Generation records the albedo and the greenhouse multiplier it drew, and the temperature is exactly what they give | GEN.409 | Boss 11:24 PM EDT; the inverse feeds the migration. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.422 | Stop the recalculation redraw and order the steps so pressure, mixture and temperature agree after the cold-side freeze-out | GEN.411, GEN.412 | Composition thread finding; Boss 11:24 PM EDT chain. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.423 | Atmosphere density, scale height and the greenhouse multiplier follow from the mixture | GEN.410, GEN.411, GEN.412 | Composition thread finding; Boss 11:24 PM EDT chain. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.424 | The balance gas is a named gas with its own molar mass, and habitability reads and writes the mean molar mass of the stored mixture | GEN.410, GEN.411 | Composition thread finding; Boss 11:24 PM EDT chain. Core 4 proposed (Boss, 11:41 PM EDT). |
| DB.67 | Planets and moons store the albedo, the greenhouse multiplier and the molar density they were built with | GEN.411, GEN.412, GEN.422, GEN.423, GEN.424 | Boss 11:24 PM EDT; migration revision after Core 1's 0102 to 0105 stack. Core 4 proposed (Boss, 11:41 PM EDT). |

## Phase 3: Ozone, haze, dust and cloud

The ozone figures and the scatterers that the colour needs, each its own thread.

### Ozone

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.413 | Ozone in the mixture: the column and partial pressure from the oxygen, the star's UV and the Chapman balance | GEN.408, GEN.411 | Boss 11:24 PM EDT; real research, generated where missing. Core 4 proposed (Boss, 11:41 PM EDT). |

### Haze, dust and cloud

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.414 | Haze, dust and cloud: derived or drawn optical depths for every planet and moon | GEN.411, GEN.412 | Boss 11:24 PM EDT; draws and derivations, nothing fixed per class. Core 4 proposed (Boss, 11:41 PM EDT). |

## Phase 4: Stored scatterers, the text and the proof

The columns for ozone and the aerosols, the description text written from the mixture, the consistency sweep and the docs.

### Consistency

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.68 | Planets and moons store the ozone column and the haze, dust and cloud optical depths | GEN.413, GEN.414, DB.67 | Boss 11:24 PM EDT; migration only. Core 4 proposed (Boss, 11:41 PM EDT). |
| GEN.415 | Every compound the planet descriptions mention comes from the gas table, and the atmosphere text is written from the stored mixture | GEN.410, GEN.411, GEN.413 | Boss 11:24 PM EDT; the description side of consistency. Core 4 proposed (Boss, 11:41 PM EDT). |
| TEST.209 | One consistency sweep over every planet and moon class: composition, density, temperature, ozone, aerosols and text all agree | GEN.411, GEN.412, GEN.413, GEN.414, GEN.415, DB.67, DB.68, GEN.422, GEN.423, GEN.424 | Boss 11:24 PM EDT; the proof. |
| DOC.63 | Docs: the atmosphere data model, the gas table and how the values agree | TEST.209 | Boss 11:24 PM EDT; docs only. |

## Phase 5: New worlds

Planet classes, nebula conditions, asteroid fields and habitability on the pages.

### Planet classes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.33 | One class per PR, each with its tests |  | Moved to phase 2 under GEN.90, after the habitability score. One class per PR (R and S first). Built on fixed physics so new classes aren't tuned to wrong masses, moons or zones. |
| GEN.146 | Teff-dependent habitable zone from the Kopparapu table |  | Research: in GEN.91's dependency chain. |
| GEN.355 | Atmosphere retention module with its test table | GEN.33 | Split from GEN.91. |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | GEN.33, GEN.146 | Under GEN.90. Same reconcile/zone code as phase 0's physics fixes. |
| GEN.350 | Class R, sub-Neptune, with the envelope-survival rule, the Hycean flag and the retirement of class T | GEN.33, GEN.355, GEN.146 | Split from GEN.28. |
| GEN.351 | Class U, icy world | GEN.33 | Split from GEN.28. |
| GEN.352 | Classes W and X, subsurface ocean bodies | GEN.33 | Split from GEN.28. |
| GEN.353 | Class Y, Titan-like world | GEN.355, GEN.146 | Split from GEN.28. |
| GEN.92 | Life and its highest stage follow the habitability score | GEN.350, GEN.351, GEN.352, GEN.353 |  |
| GEN.354 | Class Z, lifeless temperate world | GEN.92, GEN.146 | Split from GEN.28. |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | GEN.350, GEN.351, GEN.352, GEN.353, GEN.354 | Under GEN.90; class Z is Boss's Earth-size world that never had life (2026-10-03). Class S landed with GEN.38 (PR #415); the other six classes here. Closes with GEN.33's PRs; PLANET_CLASSES in program_constants.py. |
| GEN.356 | Class N re-scoped as the hot, lifeless world; its life chemical and hidden timeline removed | GEN.355, GEN.92 | Split from GEN.91. |
| GEN.91 | Classes like S and V in the hot and cold zones | GEN.355, GEN.356 |  |
| GEN.145 | Class S atmosphere rule: S keeps air unless the shoreline ratio is over 30 | GEN.355, GEN.92 | Research follow-up to class S (built). |
| GEN.29 | Sweep every planet class for sense once the new ones are in (bug) | GEN.28, GEN.27, GEN.91, GEN.92, GEN.146 | Under GEN.90: the refactor is the sweep. Bug, but by definition a sweep after the new classes; it can't go earlier. Includes rocky rogues of 10-16 Earth masses (up to 17,600 km) that get S as nearest fit (PR #415). |
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
| UX.32 | Planet rows show the class only, without the type and moon labels |  | Same rows as UX.105 (absorbed UX.30); one thread. |
| UX.101 | The PHI-6 symbols as vector shapes in the icon library, both themes, coloured by rating |  |  |
| UX.100 | Planet and moon lists show habitability as symbols, not notes | UX.101 |  |

### Asteroid fields

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.360 | Asteroid fields and belts render from their uid: one seeded generator in a shared asteroidfield.js for dust, rocks and heroes |  | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
| GEN.361 | The asteroid field page has a 3D view | GEN.360 | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
| GEN.362 | Rocks are instanced meshes in chunks with three levels of detail, within a frame budget | GEN.361 | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
| GEN.363 | Belt bodies move on the GPU: the vertex shader solves Kepler from per-instance elements | GEN.362 | Split from GEN.112, plan approved by Boss 3:00 PM EDT. |
