# Phase 3: One engine, 3D and remote generation

Rebuilt on 2026-10-02 from the dependency report (every open item, its
prerequisites and the files it shares), with Boss's decisions of that
night. [docs/TODO.md](../TODO.md) is the master file: it holds each
item's full text, and its "Plan: phases" section indexes every phase.
This file gives the phase's goal, its build threads with each item's
prerequisites in order, and its open questions. Where the two disagree,
TODO.md wins; when an item ships, it leaves TODO.md and its row here is
deleted in the same PR. Research notes, the files several items share and
the judgment calls behind the placement are in [notes.md](notes.md).

## Goal

The three maps on one engine, the 3D system view, courses that bend around gravity wells, and remote generation through the API, reproducing what the server would make.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Engine

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.66 | The sector as the drill-down's last stage, on the same page | MAP.65, TEST.70, MAP.80 | Sector as the last stage. |
| MAP.67 | One URL and history scheme for every level | MAP.66, NAV.7 |  |
| MAP.68 | Remove the old Sector Map code | MAP.67, MAP.79 | Deletes sectormap.js. |
| MAP.61 | One map engine and control set for the Galaxy Map and the Sector Map | MAP.63, MAP.64, MAP.65, MAP.66, MAP.67, MAP.68 | Parent; closes with its subitems. |

### 3D system

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.71 | Scale modes that keep everything visible | MAP.69, MAP.89 | Its compressed scale should reuse MAP.89's fitted knots. |
| MAP.72 | Rendering at system scale | MAP.69 |  |
| MAP.73 | Free camera on the shared engine | MAP.72, MAP.64, NAV.13 |  |
| MAP.74 | The 3D view on the system page, the flat diagram kept | MAP.73, MAP.71, UX.27, UX.31 | system.html, after the page's button and edit rework. |
| MAP.62 | A full 3D star system view with a free camera | MAP.69, MAP.70, MAP.71, MAP.72, MAP.73, MAP.74 | Parent; closes with its subitems. |

### Picker

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.32 | Every Galaxy and Sector Map control works on the navigation screens (bug) | MAP.68, NAV.15, NAV.29, NAV.30, NAV.31, NAV.33, MAP.79 | Bug, but it is the 'pick mode uses the one engine' end state. |
| NAV.3 | One shared picker for the Galaxy, Sector and System displays | NAV.13, NAV.14, NAV.15, NAV.16, NAV.32 | Parent; closes with its subitems. |

### Courses

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.22 | Courses inside a sector and a system | NAV.20, NAV.16, MAP.66, MAP.73 |  |
| NAV.23 | Open a saved course on the map | NAV.4, NAV.20, NAV.21 |  |
| NAV.5 | Show a course on the Galaxy Map | NAV.20, NAV.21, NAV.22, NAV.23 | Parent; most of it exists (MAP.27). |
| NAV.25 | Find the obstacles along a path | NAV.10, NAV.24, NAV.38 | Corridor query from NAV.10; sectors along the line from NAV.38. |
| NAV.26 | Bend the path around keep-out spheres | NAV.25 |  |
| NAV.27 | Moving bodies inside a system | NAV.26, MAP.70, NAV.16 |  |
| NAV.28 | Show and save the adjusted course | NAV.26, NAV.4, NAV.20 |  |
| NAV.6 | Courses that steer clear of gravity wells | NAV.24, NAV.25, NAV.26, NAV.27, NAV.28 | Parent; closes with its subitems. |

### API

| ID | Item | Needs | Note |
|---|---|---|---|
| API.13 | Generation without a database | API.12, PERF.21 | Local generation with all workers. |
| API.14 | Upload routes, compressed, in batches | API.7, API.11 |  |
| API.8 | Verify uploaded data before it is finalized | API.11, API.14 |  |
| ADM.13 | Incomplete uploads page | API.10, API.8 |  |
| API.3 | Remote generate: generate locally, upload through the API | GEN.51, GEN.52, API.9, API.10, API.11, API.12, API.13, API.14, API.8, ADM.13 | Parent; remote mode mirrors the CLI options GEN.51/52 settle. |
| API.17 | Remote generation reproduces what the server would make | API.12, API.13, GEN.57, GEN.58 |  |

### Pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.21 | Clean up the web interface: overlapping buttons and dead controls (bug) | MAP.55, MAP.60, MAP.68, UX.26, UX.27, UX.31, NAV.32 | Bug, but a final pass over the finished pages. Judgment: its one known dead control (nebula '-' at the 1 ly limit) could be split out into phase 0. |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.59 | Edits and time evolution recorded as layers on top of the seed | GEN.56, GEN.58, ADM.18 | Boss 02:20Z: net admin changes per object and per-thing regeneration seeds in ADM.18's JSON. |

## Open questions for Boss

- MAP.74: The 3D view on the system page, the flat diagram kept, see its entry in TODO.md.
- NAV.6: Courses that steer clear of gravity wells, see its entry in TODO.md.
- API.8: Verify uploaded data before it is finalized, see its entry in TODO.md.
- ADM.13: Incomplete uploads page, see its entry in TODO.md.
