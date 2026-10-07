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

The 3D system view and infinite zoom from galaxy to moon on one interface, orbital trajectories in their own frame, courses that bend around gravity wells, remote generation through the API (with galaxy-scale recipes) reproducing what the server would make, repair that reads the newest settings JSON, and the anomalies chosen in phase 2.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### 3D system

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.71 | Scale modes that keep everything visible | MAP.69, MAP.89 | Its compressed scale should reuse MAP.89's fitted knots. |
| MAP.72 | Rendering at system scale | MAP.69 |  |
| MAP.73 | Free camera on the shared engine | MAP.72, NAV.13 |  |
| MAP.74 | The 3D view on the system page, the flat diagram kept | MAP.73, MAP.71, UX.27, UX.31 | system.html, after the page's button and edit rework. |
| MAP.62 | A full 3D star system view with a free camera | MAP.69, MAP.70, MAP.71, MAP.72, MAP.73, MAP.74 | Parent; closes with its subitems. |
| MAP.126 | Show the orbital trajectories of selected objects in their frame of reference | MAP.70, MAP.73 |  |

### Courses

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.22 | Courses inside a sector and a system | NAV.20, NAV.16, MAP.66, MAP.73 |  |
| NAV.23 | Open a saved course on the map | NAV.4, NAV.20, NAV.21 |  |
| NAV.5 | Show a course on the Galaxy Map | NAV.20, NAV.21, NAV.22, NAV.23 | Courses stay drawn until cleared (NAV.49). Parent; most of it exists (MAP.27). |
| NAV.25 | Find the obstacles along a path | NAV.10, NAV.24 | Corridor query from NAV.10; sectors along the line from NAV.38. |
| NAV.26 | Bend the path around keep-out spheres | NAV.25 |  |
| NAV.27 | Moving bodies inside a system | NAV.26, MAP.70, NAV.16 |  |
| NAV.28 | Show and save the adjusted course | NAV.26, NAV.4, NAV.20 |  |
| NAV.6 | Courses that steer clear of gravity wells | NAV.24, NAV.25, NAV.26, NAV.27, NAV.28 | Parent; closes with its subitems. |

### API

| ID | Item | Needs | Note |
|---|---|---|---|
| API.13 | Generation without a database | API.12 | Local generation with all workers. |
| API.14 | Upload routes, compressed, in batches | API.7, API.11 |  |
| API.8 | Verify uploaded data before it is finalized | API.11, API.14 |  |
| ADM.13 | Incomplete uploads page | API.10, API.8 |  |
| API.3 | Remote generate: generate on a local machine, upload through the API | API.9, API.10, API.11, API.12, API.13, API.14, API.8, ADM.13 | Parent; remote mode mirrors the CLI options GEN.51/52 settle. |
| API.17 | Remote generation reproduces what the server would make | API.12, API.13, GEN.57, GEN.58 |  |

### Database consistency check

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.10 | Repair reads the newest settings JSON and the pending deltas | DB.9, GEN.61, OPS.18 | Falls back to the next backup if the newest JSON is damaged. |

### Daily maintenance

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.20 | A "merge now" button on the Admin dashboard (low priority) | OPS.16, GEN.61, OPS.18 | Boss 02:31Z: phase 3, low priority. Same lock and rules as the daily run; counts toward the day's slot. |

### Anomalies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.114 | Add the chosen anomalies to the starmap | GEN.113 | Placeholder until the analysis picks them. |

### Infinite zoom

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.125 | Infinite zoom: one 3D interface from the galaxy down to a moon | MAP.61, MAP.62 | The end state of the one-engine and 3D-system work. |

### Recipes

| ID | Item | Needs | Note |
|---|---|---|---|
| API.19 | Galaxy-scale recipes: build a whole galaxy, piece by piece, from JSON | API.18 | The far end; framework and capability, not a full Milky Way. |

## Open questions for Boss

- MAP.74: The 3D view on the system page, the flat diagram kept, see its entry in TODO.md.
- NAV.6: Courses that steer clear of gravity wells, see its entry in TODO.md.
- API.8: Verify uploaded data before it is finalized, see its entry in TODO.md.
- ADM.13: Incomplete uploads page, see its entry in TODO.md.
