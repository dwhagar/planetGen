# Version 9.0: Atmosphere colour

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 9.0

## Goal

The colour of the air of any world, worked out from the atmosphere it really has: the one gas-and-compound table, the mixture of partial pressures that sets the molar density, ozone, haze, dust and cloud, and the albedo and greenhouse values stored with every planet and moon. The atmosphere colour (Boss: sky colour and atmosphere colour are the same thing) is the absolute colour, not adapted to the star, and it covers gas giants too; it shows in the descriptions and replaces the planet colour on the maps where that is what it meant. Boss (2026-10-16, 11:24 PM EDT, 11:25 PM EDT): this is the seminal feature of 9.0, and every value must agree across the board. The foundations (the table, the mixture, the stored albedo and greenhouse values, ozone, haze, dust and cloud, and their tests) are built in 8.4 phases 1 to 4, so 9.0 holds the colour itself and what shows it. A major version needs Boss's word. The plan note is [atmosphere-colour.md](../design/atmosphere-colour.md), based on [sky-colour-analysis.md](../design/sky-colour-analysis.md). The view from a planet and the gravity map are now 9.1, and the backups and repair chain is 8.9 (Boss, 11:26 PM EDT: "Everything for 9.0 before splits between 9.1 and 8.9").

## Phase 1: The colour model

The physics of the colour: the blackbody colour of the star, the Rayleigh optical depth of the mixture, the colour of a solid world's air and of a giant's cloud top, all as the absolute colour (not adapted to the star).

### Atmosphere colour

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.416 | Blackbody colour of a star: Planck to the CIE curves to sRGB, as the first part of the photometry library |  | Boss 11:24 PM EDT; VIEW.11 extends it. |
| GEN.417 | Rayleigh optical depth of the mixture: per-gas cross sections from the table and the air column from pressure and gravity | GEN.410, GEN.411 | Boss 11:24 PM EDT; pure function of the stored mixture. |
| GEN.418 | Atmosphere colour of a solid world: the absolute colour of the air with ozone, haze, dust, cloud and the ground's bounce | GEN.416, GEN.417, GEN.413, GEN.414, DB.67, DB.68 | Boss 11:24 PM EDT; absolute colour, not eye-adapted. |
| GEN.419 | Atmosphere colour of gas and ice giants, read at a reference level in the cloud tops | GEN.418 | Boss 11:24 PM EDT; reference level chosen in the design. |



## Phase 2: Atmosphere colour in the product

The colour shown in the descriptions, on the maps in place of the planet colour, and in the documents.

### Descriptions and maps

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.420 | The descriptions show the atmosphere colour: a swatch and a word in the planet and moon facts, in every output | GEN.418, GEN.419 | Boss 11:24 PM EDT; one place serves every output. |
| MAP.244 | The maps and pages use the atmosphere colour where they showed the planet colour | GEN.420 | Boss 11:24 PM EDT; replaces the temperature gesture. |
| DOC.64 | Docs: the atmosphere colour, its method and its sources | GEN.420, MAP.244 | Boss 11:24 PM EDT; docs only. |
| GEN.421 | Headline of 9.0: the atmosphere colour of every world, from its real atmosphere, with every stored value consistent | GEN.420, MAP.244, DOC.64, TEST.209 | Headline of 9.0 (Boss, 11:25 PM EDT). |
