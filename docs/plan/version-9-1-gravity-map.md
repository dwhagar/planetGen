# Version 9.1: The gravity map

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 9.1

## Goal

A heat map of the gravitational field inside a sector and across a star system (MAP.167), moved from 9.0 because Boss does not think it is the best 9.0 headline.

## Phase 1: The gravity grid

The per-zone grid, its cost and the API route.

### Gravity grid

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.72 | Research the cost of a sector's gravity grid and where to cut between exact and aggregated sources | GEN.198, PERF.88 | Research. |
| GEN.199 | Per-zone gravity grid for a sector: stored in a cache and rebuilt when the orbits update | GEN.198, PERF.72 |  |
| API.24 | GET /api/sectors/<id>/gravity: the per-zone gravity grid as JSON | GEN.199 |  |

## Phase 2: The gravity layers

The Sector Map and System Map layers, the wording and the umbrella.

### Gravity layers

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.168 | The Sector Map gravity layer: coloured zones for pull, well depth and tidal strength | API.24, MAP.217 |  |
| MAP.169 | The System Map gravity layer: the orbital plane as a heat map with Lagrange points and Hill spheres | MAP.168, MAP.217, GEN.289 |  |
| MAP.170 | The Galaxy Map gravity layer: the galaxy potential and region aggregates, coarse and optional | MAP.168, MAP.217 | Later; only if PERF.72 says it is cheap. |
| UX.94 | Gravity map wording: the legend, the mode names and the Map help text | MAP.168 |  |
| MAP.167 | Gravity map: a heat map of the gravitational field inside a sector (the seminal feature of 9.0) | GEN.199, MAP.168, MAP.169, UX.94, API.24 | Seminal feature of 9.0 (Boss 2026-10-10). Parent. |
