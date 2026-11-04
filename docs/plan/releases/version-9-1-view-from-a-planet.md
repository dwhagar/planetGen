# Version 9.1: The view from a planet

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 9.1

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

A heat map of the gravitational field inside a sector and across a star system (MAP.167), moved from 9.0 because Boss does not think it is the best 9.0 headline. The view from a planet (phase 1: the photometry library, the sky catalogue, the star map and the PNG) joins it here from 9.0 (Boss, 2026-10-16, 11:26 PM EDT: "Everything for 9.0 before splits between 9.1 and 8.9"), because it builds on the atmosphere colour of 9.0; its sky data and research are in 8.9. Boss (2026-10-16, 11:42 PM EDT, card "View from a planet"): the view from a planet is the headline of 9.1 and the gravity map is the second feature.

## Phase 1: Sky rendering

The view from a planet, moved here from 9.0 (Boss, 10:27 PM EDT): the photometry library, the sky catalogue, the star map and the PNG. It builds on the atmosphere colour of 9.0 and on the sky data and research of 8.9.

### Sky rendering

| ID | Item | Needs | Note |
|---|---|---|---|
| VIEW.11 | Photometry library: apparent magnitude, bolometric correction, blackbody colour and the extinction law | VIEW.1 | Split from VIEW.2. Extends GEN.416 (the blackbody colour, 9.0). |
| VIEW.12 | Sky catalogue: the stars an observer sees, with direction, brightness and light-time | VIEW.11, VIEW.7 | Split from VIEW.2. |
| VIEW.13 | Galactic band: the glow of the unresolved stars and dust along each line of sight | VIEW.11, VIEW.8 | Split from VIEW.2. |
| VIEW.14 | Orientation and horizon: the planet's own sky, time of day and twilight | VIEW.12, VIEW.13, VIEW.6 | Split from VIEW.2. |
| VIEW.2 | A starmap seen from a planet. RESEARCH WITH BOSS FIRST | VIEW.11, VIEW.12, VIEW.13, VIEW.14 |  |
| VIEW.15 | Whole-sky PNG: Hammer-Aitoff projection, written with numpy and zlib, cached and served by the API | VIEW.11, VIEW.12, VIEW.13 | Split from VIEW.3. |
| VIEW.16 | Constellations: the seeded partition over the brightest stars, stored per system, named by the codec, with Bayer labels | VIEW.4, VIEW.12 | Split from VIEW.3. |
| VIEW.17 | Constellation overlay: SVG lines and names over the sky chart, and the planet-sky card | VIEW.15, VIEW.16 | Split from VIEW.3. |
| MAP.285 | Close-up tier: fragment-shader detail over about 500 px, shared with the sky view | VIEW.14 | Body rendering R16; Boss's plan of 10:26 PM EDT 2026-10-23. |
| VIEW.3 | Render the view as a PNG, with constellations | VIEW.15, VIEW.16, VIEW.17 |  |
| VIEW.18 | Headline of 9.1: the view from a planet, with its starmap, its picture and its constellations | VIEW.2, VIEW.3 | Headline of 9.1 (Boss, 11:42 PM EDT). |

## Phase 2: The gravity grid

The per-zone grid, its cost and the API route.

### Gravity grid

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.72 | Research the cost of a sector's gravity grid and where to cut between exact and aggregated sources | GEN.198 | Research. |
| GEN.199 | Per-zone gravity grid for a sector: stored in a cache and rebuilt when the orbits update | GEN.198, PERF.72 |  |
| API.24 | GET /api/sectors/<id>/gravity: the per-zone gravity grid as JSON | GEN.199 |  |

## Phase 3: The gravity layers

The Sector Map and System Map layers, the wording and the umbrella.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### Gravity layers

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.168 | The Sector Map gravity layer: coloured zones for pull, well depth and tidal strength | API.24 |  |
| MAP.169 | The System Map gravity layer: the orbital plane as a heat map with Lagrange points and Hill spheres | MAP.168 |  |
| MAP.170 | The Galaxy Map gravity layer: the galaxy potential and region aggregates, coarse and optional | MAP.168 | Later; only if PERF.72 says it is cheap. |
| UX.94 | Gravity map wording: the legend, the mode names and the Map help text | MAP.168 |  |
| MAP.167 | Gravity map: a heat map of the gravitational field inside a sector (the seminal feature of 9.0) | GEN.199, MAP.168, MAP.169, UX.94, API.24 | Seminal feature of 9.0 (Boss 2026-10-10). Parent. |

## Phase 4: Living population and trade

Populations and trade routes become dynamic (Boss, 9:38 PM EDT 2026-10-24).

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### Dynamic population and trade

| ID | Item | Needs | Note |
|---|---|---|---|
| POP.35 | Dynamic world populations: populations step forward with the nightly update from the stored growth state | POP.30 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| POP.36 | Dynamic trade routes: routes grow, shift and close as populations, ports and territories change, using the same hop pathfinding | POP.33, POP.35, MAP.309 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
