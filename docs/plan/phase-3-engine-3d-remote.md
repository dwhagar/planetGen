# Phase 3: Infinite zoom, 3D, remote generation and anomalies

Rebuilt on 2026-10-07 from Boss's lists of 2026-10-03 and 2026-10-07
and the dependency tree (every open item, its prerequisites and the
files it shares). [docs/TODO.md](../TODO.md) is the master file: it holds each
item's full text, and its "Plan: phases" section indexes every phase.
This file gives the phase's goal, its build threads with each item's
prerequisites in order, and its open questions. Where the two disagree,
TODO.md wins; when an item ships, it leaves TODO.md and its row here is
deleted in the same PR. Research notes, the files several items share and
the judgment calls behind the placement are in [notes.md](notes.md).

## Goal

The 3D system view and infinite zoom from galaxy to moon on one interface, orbital trajectories in their own frame, courses that bend around gravity wells, remote generation through the API (with galaxy-scale recipes) reproducing what the server would make and the anomalies chosen in phase 2.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### 3D system

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.167 | Gravity map: a heat map of the gravitational field inside a sector (the seminal feature of 9.0) | GEN.199, MAP.168, MAP.169, UX.94, API.24 | Seminal feature of 9.0 (Boss 2026-10-10). Parent. |
| API.24 | GET /api/sectors/<id>/gravity: the per-zone gravity grid as JSON | GEN.199 |  |
| MAP.168 | The Sector Map gravity layer: coloured zones for pull, well depth and tidal strength | API.24, MAP.125 |  |
| MAP.169 | The System Map gravity layer: the orbital plane as a heat map with Lagrange points and Hill spheres | MAP.168 |  |
| UX.94 | Gravity map wording: the legend, the mode names and the Map help text | MAP.168 |  |

### Courses

| ID | Item | Needs | Note |
|---|---|---|---|

### API

| ID | Item | Needs | Note |
|---|---|---|---|
| API.13 | Generation without a database | API.12 | Local generation with all workers. |
| API.14 | Upload routes, compressed, in batches | API.7, API.11 |  |
| API.8 | Verify uploaded data before it is finalized | API.11, API.14 |  |
| ADM.13 | Incomplete uploads page | API.10, API.8 |  |
| API.3 | Remote generate: generate on a local machine, upload through the API | API.10, API.11, API.12, API.13, API.14, API.8, ADM.13 | Parent; remote mode mirrors the CLI options GEN.51/52 settle. |
| API.17 | Remote generation reproduces what the server would make | API.12, API.13, GEN.57 |  |

### Database consistency check

| ID | Item | Needs | Note |
|---|---|---|---|

### Daily maintenance

| ID | Item | Needs | Note |
|---|---|---|---|

### Anomalies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.114 | Add the chosen anomalies to the starmap | GEN.113 | Placeholder until the analysis picks them. |
| GEN.163 | Type-B pulsar planets in globular clusters (GEN.130 follow-on) | GEN.130, GEN.158, GEN.159 | Globular-cluster chain. |
| GEN.162 | Planet cull and blue stragglers in clusters | GEN.158, GEN.160 | Globular-cluster chain. |
| GEN.161 | Bright-first fill for cluster sectors | GEN.160 | Globular-cluster chain. |
| GEN.160 | Cluster density in the sector gate, with a "cluster" population | GEN.159 | Globular-cluster chain. |
| GEN.159 | Globular clusters: cluster table, King tables and the Milky Way catalogue | GEN.158 | Globular-cluster chain. |
| GEN.158 | Add a metallicity value to stars |  | Globular-cluster chain. |

### Infinite zoom

| ID | Item | Needs | Note |
|---|---|---|---|

### Recipes

| ID | Item | Needs | Note |
|---|---|---|---|
| API.19 | Galaxy-scale recipes: build a whole galaxy, piece by piece, from JSON | API.18 | The far end; framework and capability, not a full Milky Way. |

## Open questions for Boss

- NAV.6: Courses that steer clear of gravity wells, see its entry in TODO.md.
- API.8: Verify uploaded data before it is finalized, see its entry in TODO.md.
- ADM.13: Incomplete uploads page, see its entry in TODO.md.
