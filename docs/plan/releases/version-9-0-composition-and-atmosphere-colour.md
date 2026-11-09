# Version 9.0: Composition, atmosphere colour and the view from a planet

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 9.0

Status: paused until 8.1 is complete (its composition, radiation and colour work moved into 8.1 phases 5 to 11 on Boss's word, 5:48 PM EDT 2026-10-23; only the proofs and the headline stay here). Earlier note (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

Boss (5:31 PM EST on 2026-11-05): 9.1 is folded into 9.0, so this version holds three things: the composition and atmosphere colour, the view from a planet (phase 2, with its sky data, research and dust model built in 8.2 as foundations) and the gravity map (finished in phase 3 from its foundations in 8.2, 8.3 and 8.5). The dynamic population and trade moved to 9.2 (Boss, 6:00 PM EST on 2026-11-05). The numbers 9.2 onward stay as they are.

What every body, star and atmosphere is made of, worked out from real science, and what follows from it: the colour of the air and the radiation in the air and on the ground. Boss (2026-10-16, 11:25 PM EDT and 11:58 PM EDT): atmosphere composition and the composition of planets, moons, stars and everything is the joint headline of 9.0 with the atmosphere colour; every body has figures on its composition based on actual science, not made up where it can be avoided, and drawing the density first and filling in a fitting composition is as acceptable as deriving it. The atmosphere colour (sky colour and atmosphere colour are the same thing) is the absolute colour, not adapted to the star, and covers gas giants; it shows in the descriptions and replaces the planet colour on the maps. The atmospheric radiation (band transmission and the greenhouse effect from the mixture) and the surface radiation (the dose from the air's shielding and the ground's own composition) are derived the same way. Everything else is built before: the research, the designs and the data tables (the elements, the gases, the materials, the stars) are in 8.1 phases 1 and 2, the star composition function in 8.1 phase 3, and the generation and storage in 8.4 phases 1 to 4, so 9.0 holds the colour, the radiation and what shows them. A major version needs Boss's word. The master note is [composition-model.md](../../design/composition/composition-model.md); the atmosphere plan is [atmosphere-colour.md](../../design/composition/atmosphere-colour.md), based on [sky-colour-analysis.md](../../design/composition/sky-colour-analysis.md). The view from a planet and the gravity map are 9.1, and the backups and repair chain is 8.9 (Boss, 11:26 PM EDT).

## Phase 1: The proofs and the headline

The radiation sweep, the documents of the whole composition model and the headline.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### The headline

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.212 | A radiation sweep over every class: transmission, greenhouse, shielding and dose agree with the composition | GEN.444, GEN.446 | Boss 11:58 PM EDT; the radiation proof. |
| DOC.66 | Docs: the composition model of stars, bodies and atmospheres, the radiation it derives and the sources of every number | TEST.212 | Boss 11:58 PM EDT; docs only. |
| GEN.421 | Headline of 9.0: what every body and atmosphere is made of, and the atmosphere colour that follows from it | GEN.446, TEST.212, DOC.66 | Headline of 9.0 (Boss, 11:25 PM EDT). |
| MAP.284 | Night-side thermal glow on hot worlds, with no city lights until 9.2 | MAP.280 | Body rendering R15; Boss's plan of 10:26 PM EDT 2026-10-23. |
| MAP.283 | Comet colour from the coma species |  | Body rendering R9b; Boss's plan of 10:26 PM EDT 2026-10-23. |
| MAP.282 | Life tint: the pigment colour of the life chemistry, faint, on worlds rated Yellow or better; life that is not green mixes in the colour of its own evolutionary chemistry | MAP.280 | Body rendering R14; Boss's plan of 10:26 PM EDT 2026-10-23. |
| MAP.281 | Giant chemistry bands once the band strengths exist |  | Body rendering R6b; Boss's plan of 10:26 PM EDT 2026-10-23. |
| MAP.280 | Cloud layer, terminator and cover on solid worlds; cloud-deck worlds paint the deck as the surface |  | Body rendering R4b; Boss's plan of 10:26 PM EDT 2026-10-23. |
| MAP.416 | Census of cloud and haze data by class | none | Boss, 2026-11-09 (render report); haze and clouds need not be in 8.1. |
| MAP.418 | (bug) Cloud and haze colour come from the dominant condensate or haze material | MAP.416 | Boss, 2026-11-09 (render report); haze and clouds need not be in 8.1. |
| MAP.420 | Coloured haze shell and limb, sized from tau_haze | MAP.418 | Boss, 2026-11-09 (render report); haze and clouds need not be in 8.1. |
| MAP.421 | Giants get a haze and cloud-deck tint | MAP.418 | Boss, 2026-11-09 (render report); haze and clouds need not be in 8.1. |
| MAP.422 | The info panel line says which of cloud and haze were drawn, and from what | MAP.418 | Boss, 2026-11-09 (render report); haze and clouds need not be in 8.1. |

### Resource layer

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.313 | Galaxy Map resource richness layer: a heat map of the resource score by sector, by resource kind, with the richest sectors listed | MAP.312, DB.106 | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |

## Phase 2: Sky rendering (the view from a planet, folded in from 9.1)

The view from a planet, moved here from 9.0 and back with it (Boss, 10:27 PM EDT; 9.1 folded into 9.0 at 5:35 PM EST on 2026-11-05): the photometry library, the sky catalogue, the star map and the PNG. It builds on the atmosphere colour of 9.0 and on the sky data and research of 8.9.

### Sky rendering

| ID | Item | Needs | Note |
|---|---|---|---|
| VIEW.14 | Orientation and horizon: the planet's own sky, time of day and twilight | VIEW.12, VIEW.13, VIEW.6 | Split from VIEW.2. |
| VIEW.2 | A starmap seen from a planet. RESEARCH WITH BOSS FIRST | VIEW.11, VIEW.12, VIEW.13, VIEW.14 |  |
| VIEW.15 | Whole-sky PNG: Hammer-Aitoff projection, written with numpy and zlib, cached and served by the API | VIEW.11, VIEW.12, VIEW.13 | Split from VIEW.3. |
| VIEW.16 | Constellations: the seeded partition over the brightest stars, stored per system, named by the codec, with Bayer labels | VIEW.4, VIEW.12 | Split from VIEW.3. |
| VIEW.17 | Constellation overlay: SVG lines and names over the sky chart, and the planet-sky card | VIEW.15, VIEW.16 | Split from VIEW.3. |
| MAP.285 | Close-up tier: fragment-shader detail over about 500 px, shared with the sky view | VIEW.14 | Body rendering R16; Boss's plan of 10:26 PM EDT 2026-10-23. |
| VIEW.3 | Render the view as a PNG, with constellations | VIEW.15, VIEW.16, VIEW.17 |  |
| VIEW.18 | Headline of 9.0: the view from a planet, with its starmap, its picture and its constellations | VIEW.2, VIEW.3 | Headline of 9.0 (Boss, 11:42 PM EDT). |

## Phase 3: The gravity map closes

The Sector Map and System Map layers, the wording and the umbrella.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### Gravity layers

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.170 | The Galaxy Map gravity layer: the galaxy potential and region aggregates, coarse and optional |  | Later; only if PERF.72 says it is cheap. |
| MAP.167 | Gravity map: a heat map of the gravitational field inside a sector (the seminal feature of 9.0) | GEN.199 | Seminal feature of 9.0 (Boss 2026-10-10). Parent. |
