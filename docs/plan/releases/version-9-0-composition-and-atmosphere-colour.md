# Version 9.0: Composition and atmosphere colour

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 9.0

Status: paused until 8.1 is complete (its composition, radiation and colour work moved into 8.1 phases 5 to 11 on Boss's word, 5:48 PM EDT 2026-10-23; only the proofs and the headline stay here). Earlier note (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

What every body, star and atmosphere is made of, worked out from real science, and what follows from it: the colour of the air and the radiation in the air and on the ground. Boss (2026-10-16, 11:25 PM EDT and 11:58 PM EDT): atmosphere composition and the composition of planets, moons, stars and everything is the joint headline of 9.0 with the atmosphere colour; every body has figures on its composition based on actual science, not made up where it can be avoided, and drawing the density first and filling in a fitting composition is as acceptable as deriving it. The atmosphere colour (sky colour and atmosphere colour are the same thing) is the absolute colour, not adapted to the star, and covers gas giants; it shows in the descriptions and replaces the planet colour on the maps. The atmospheric radiation (band transmission and the greenhouse effect from the mixture) and the surface radiation (the dose from the air's shielding and the ground's own composition) are derived the same way. Everything else is built before: the research, the designs and the data tables (the elements, the gases, the materials, the stars) are in 8.1 phases 1 and 2, the star composition function in 8.1 phase 3, and the generation and storage in 8.4 phases 1 to 4, so 9.0 holds the colour, the radiation and what shows them. A major version needs Boss's word. The master note is [composition-model.md](../../design/composition/composition-model.md); the atmosphere plan is [atmosphere-colour.md](../../design/composition/atmosphere-colour.md), based on [sky-colour-analysis.md](../../design/composition/sky-colour-analysis.md). The view from a planet and the gravity map are 9.1, and the backups and repair chain is 8.9 (Boss, 11:26 PM EDT).

## Phase 1: The proofs and the headline

The radiation sweep, the documents of the whole composition model and the headline.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### The headline

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.212 | A radiation sweep over every class: transmission, greenhouse, shielding and dose agree with the composition | GEN.443, GEN.444, GEN.445, GEN.446 | Boss 11:58 PM EDT; the radiation proof. |
| DOC.66 | Docs: the composition model of stars, bodies and atmospheres, the radiation it derives and the sources of every number | TEST.212, TEST.211, DOC.63, DOC.65 | Boss 11:58 PM EDT; docs only. |
| GEN.421 | Headline of 9.0: what every body and atmosphere is made of, and the atmosphere colour that follows from it | GEN.420, MAP.244, DOC.64, GEN.446, TEST.209, TEST.210, TEST.211, TEST.212, DOC.66 | Headline of 9.0 (Boss, 11:25 PM EDT). |
| UX.158 | Info panel line: how this colour was worked out, which inputs and whether partial or full quality | MAP.271 | Body rendering R12; proposed 9.0. |
| TEST.308 | Body rendering tests: look record stability, painter determinism, tier switching and the memory ceiling | MAP.272 | Body rendering R11; proposed 9.0. |
| GEN.670 | Research for body rendering: band, storm and Rossby scale laws, vegetation tint by host spectrum, ring colour by composition |  | Body rendering R10; proposed 9.0. |
| MAP.278 | Phenomena painters: rogue planet, comet, nebula, remnant, asteroid field and black hole disc | MAP.271, MAP.272, GEN.488 | Body rendering R9; proposed 9.0. |
| MAP.277 | Atmosphere glow from pressure and composition, and the day and night terminator | MAP.271, DB.96 | Body rendering R8; proposed 9.0. |
| MAP.276 | Star painter: blackbody colour shared with star_spectrum, class texture, spots and glow | MAP.271, MAP.272 | Body rendering R7; proposed 9.0. |
| MAP.275 | Gas and ice giant painter: bands, storms, chemistry tint, thermal glow and rings with gaps | MAP.271, MAP.272, GEN.422 | Body rendering R6; proposed 9.0. |
| MAP.274 | Land colour from the bulk composition: oxidation and rock type | GEN.411, GEN.428, GEN.429, MAP.273 | Body rendering R5; proposed 9.0. |
| MAP.273 | Solid world painter: ocean, land, ice and cloud from the stored fractions and the atmosphere colour | MAP.271, MAP.272, GEN.441, GEN.461, GEN.462 | Body rendering R4; proposed 9.0. |
| MAP.272 | Shared body painter module: seeded noise, texture cache, paint queue and size tiers |  | Body rendering R3; proposed 9.0. |
| DB.96 | Store the atmosphere colour and the cloud cover on planets and moons | GEN.414 | Body rendering R2; proposed 9.0. |
| MAP.271 | The bodylook record: the server builds one data-look record per body from stored columns |  | Body rendering R1; proposed 9.0. |
