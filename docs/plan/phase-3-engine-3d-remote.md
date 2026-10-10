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

### Courses

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.22 | Courses inside a sector and a system |  |  |
| NAV.23 | Open a saved course on the map | NAV.4, NAV.21 |  |
| NAV.5 | Show a course on the Galaxy Map | NAV.21, NAV.22, NAV.23 | Courses stay drawn until cleared (NAV.49). Parent; most of it exists (MAP.27). |
| NAV.25 | Find the obstacles along a path |  | Corridor query from NAV.10; sectors along the line from NAV.38. |
| NAV.26 | Bend the path around keep-out spheres | NAV.25 |  |
| NAV.56 | Census of overlapping keep-out spheres in a generated galaxy | NAV.51 | Research: decides the overlap policy. |
| NAV.55 | A tuning block for the keep-out knobs |  | Research: keep-out knobs. |
| NAV.27 | Moving bodies inside a system | NAV.26 |  |
| NAV.28 | Show and save the adjusted course | NAV.26, NAV.4 |  |
| NAV.51 | Courses route around asteroid fields | NAV.25, NAV.26 | Boss 2026-10-09 01:02Z. Parked with the NAV chain (2026-10-09). |
| NAV.54 | Keep-out radii for asteroid fields, supermassive holes, moons and nebulae (NAV.24 built) | NAV.51 | Research: NAV.24 (built) follow-ups. |
| NAV.6 | Courses that steer clear of gravity wells | NAV.25, NAV.26, NAV.27, NAV.28, NAV.51 | Parent; closes with its subitems. |

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
