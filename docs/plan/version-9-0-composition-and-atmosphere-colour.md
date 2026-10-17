# Version 9.0: Composition and atmosphere colour

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 9.0

## Goal

What every body, star and atmosphere is made of, worked out from real science, and what follows from it: the colour of the air and the radiation in the air and on the ground. Boss (2026-10-16, 11:25 PM EDT and 11:58 PM EDT): atmosphere composition and the composition of planets, moons, stars and everything is the joint headline of 9.0 with the atmosphere colour; every body has figures on its composition based on actual science, not made up where it can be avoided, and drawing the density first and filling in a fitting composition is as acceptable as deriving it. The atmosphere colour (sky colour and atmosphere colour are the same thing) is the absolute colour, not adapted to the star, and covers gas giants; it shows in the descriptions and replaces the planet colour on the maps. The atmospheric radiation (band transmission and the greenhouse effect from the mixture) and the surface radiation (the dose from the air's shielding and the ground's own composition) are derived the same way. Everything else is built before: the research, the designs and the data tables (the elements, the gases, the materials, the stars) are in 8.1 phases 1 and 2, the star composition function in 8.1 phase 3, and the generation and storage in 8.4 phases 1 to 4, so 9.0 holds the colour, the radiation and what shows them. A major version needs Boss's word. The master note is [composition-model.md](../design/composition-model.md); the atmosphere plan is [atmosphere-colour.md](../design/atmosphere-colour.md), based on [sky-colour-analysis.md](../design/sky-colour-analysis.md). The view from a planet and the gravity map are 9.1, and the backups and repair chain is 8.9 (Boss, 11:26 PM EDT).

## Phase 1: The colour model

The physics of the colour: the blackbody colour of the star, the Rayleigh optical depth of the mixture, the colour of a solid world's air and of a giant's cloud top, all as the absolute colour (not adapted to the star).

### Atmosphere colour

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.416 | Blackbody colour of a star: Planck to the CIE curves to sRGB, as the first part of the photometry library |  | Boss 11:24 PM EDT; VIEW.11 extends it. |
| GEN.417 | Rayleigh optical depth of the mixture: per-gas cross sections from the table and the air column from pressure and gravity | GEN.410, GEN.411 | Boss 11:24 PM EDT; pure function of the stored mixture. |
| GEN.418 | Atmosphere colour of a solid world: the absolute colour of the air with ozone, haze, dust, cloud and the ground's bounce | GEN.416, GEN.417, GEN.413, GEN.414, DB.67, DB.68, GEN.441 | Boss 11:24 PM EDT; absolute colour, not eye-adapted. |
| GEN.419 | Atmosphere colour of gas and ice giants, read at a reference level in the cloud tops | GEN.418 | Boss 11:24 PM EDT; reference level chosen in the design. |

## Phase 2: Atmospheric and surface radiation

The radiation derived from the stored composition: the band transmission and greenhouse of the atmosphere, the dose at the ground from the air's shielding and the ground's own content, and the habitability values that read them.

### Atmospheric and surface radiation

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.443 | Atmospheric radiation: transmission by band and the greenhouse effect from the stored mixture, the aerosols and the cloud | GEN.417, GEN.413, GEN.414, GEN.436 | Boss 11:58 PM EDT; band transmission and greenhouse from the mixture. |
| GEN.444 | Surface radiation: the ground-level dose from the air's shielding and the ground's own radiogenic content | GEN.443, GEN.441, GEN.438 | Boss 11:58 PM EDT; replaces the proxies of physics/radiation.py; generator epoch bump. |
| GEN.445 | Habitability, toxicity and the UV index read the derived radiation instead of the proxies | GEN.444 | Boss 11:58 PM EDT; habitability reads the derived values. |

## Phase 3: Composition and colour in the product

The colour shown in the descriptions, on the maps in place of the planet colour, the composition and radiation text of stars, surfaces and small bodies, and the documents of the colour.

### Descriptions and maps

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.420 | The descriptions show the atmosphere colour: a swatch and a word in the planet and moon facts, in every output | GEN.418, GEN.419 | Boss 11:24 PM EDT; one place serves every output. |
| MAP.244 | The maps and pages use the atmosphere colour where they showed the planet colour | GEN.420 | Boss 11:24 PM EDT; replaces the temperature gesture. |
| GEN.446 | The descriptions state the composition of stars, surfaces and small bodies, and the radiation at the surface | GEN.440, GEN.441, GEN.442, GEN.444, GEN.445 | Boss 11:58 PM EDT; stars, surfaces, small bodies and radiation in the text. |
| DOC.64 | Docs: the atmosphere colour, its method and its sources | GEN.420, MAP.244 | Boss 11:24 PM EDT; docs only. |

## Phase 4: The proofs and the headline

The radiation sweep, the documents of the whole composition model and the headline.

### The headline

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.212 | A radiation sweep over every class: transmission, greenhouse, shielding and dose agree with the composition | GEN.443, GEN.444, GEN.445, GEN.446 | Boss 11:58 PM EDT; the radiation proof. |
| DOC.66 | Docs: the composition model of stars, bodies and atmospheres, the radiation it derives and the sources of every number | TEST.212, TEST.211, DOC.63, DOC.65 | Boss 11:58 PM EDT; docs only. |
| GEN.421 | Headline of 9.0: what every body and atmosphere is made of, and the atmosphere colour that follows from it | GEN.420, MAP.244, DOC.64, GEN.446, TEST.209, TEST.210, TEST.211, TEST.212, DOC.66 | Headline of 9.0 (Boss, 11:25 PM EDT). |
