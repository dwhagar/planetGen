# Phase 2: Maps, picker and backfill

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

The Galaxy Map built out around the arc pick, the shared picker and courses (with unknown-space jumps marked), the parallel backfill and density pass, the update's check for changed output, the daily maintenance run (positional update, merge of the day's admin changes into a new settings JSON, 18 backups), and the API pieces remote generation needs first.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Galaxy Map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.56 | Drop the 3x3 block pick: select a slab, zoom in, select a segment (bug) |  | Arc, slab, segment ladder in galaxystages.js. |
| MAP.53 | Rotate a zoomed-in wedge, and zoom it to fit the window (bug) |  | Rotation and fit on MAP.64's controller. |
| MAP.58 | Galaxy Map zoom limits: a short manual range on the galaxy wedge, locked below it | MAP.53 | A zoom policy of MAP.64. |
| MAP.78 | Zooming into a wedge must show the whole wedge at every drill-down level (bug) | MAP.53 | The fit itself; with MAP.53. |
| MAP.76 | Leader-line layout | MAP.56 | Layout half of MAP.54; same PR. |
| MAP.54 | Slab leader lines instead of the slab slider (bug) | MAP.56, MAP.76 |  |
| MAP.75 | The mini map as a second engine view | MAP.54 | Locked second camera on MAP.64. |
| MAP.59 | Make it plain that a zoomed-in slab is a slab, not a wedge | MAP.54, MAP.53, MAP.75 |  |
| MAP.77 | Galaxy Map draws block divisions inside a picked slab before zooming to it (bug) | MAP.56, MAP.59 | Which lines show at each level, once the ladder and ghost exist. |

### Engine

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.65 | One picking, hover and info-panel layer | MAP.56 | After the selection rewrite settles the pick flow. |
| MAP.79 | Rogue planets clog the Sector Map: dim them, and a show/hide button per kind of object (bug) | MAP.65 | Judgment: the per-kind toggles go on the shared control set; the dimming already landed in phase 0 (MAP.82 to MAP.84). |

### Picker

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.15 | Pick mode everywhere | NAV.13, NAV.14, MAP.65 | Pick mode in the shared panel layer. |
| NAV.29 | Replace "Nav from here" and "Nav to here" with "Start Here" and "End Here" while picking (bug) | NAV.13, NAV.15 | Needs the step out/in of NAV.13. |
| NAV.33 | After picking one end of a course, stay at that zoom level (bug) | NAV.15 |  |
| NAV.31 | Galaxy wedges don't highlight on the navigation screens (bug) | NAV.15 | The highlight is MAP.85's arc highlight. |
| NAV.16 | NAV endpoints can be any object | NAV.7 | navigation.py legs, nav_page.py endpoints. |

### Courses

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.20 | Draw the direct line and the route apart | NAV.13 |  |
| NAV.21 | Fit the view to the whole course | MAP.53, MAP.58 | Same fit as MAP.53, exempt from MAP.58's lock. |
| NAV.17 | A saved course record with both forms | NAV.7 |  |
| NAV.18 | Save, list, open, rename and delete, per browser | NAV.17 | Sibling of bookmarks.js. |
| NAV.4 | Save a course | NAV.17, NAV.18 | Per browser now; NAV.19 moves it into accounts in phase 3+. |
| NAV.24 | A keep-out radius for every kind of object | GEN.47 | Nebula keep-out question needs nebulae that actually exist. |
| NAV.36 | Unknown-space jumps drawn red and glowing | NAV.12, UX.35, NAV.20 | Route strip, navmap.py and the Galaxy Map course; legend; reduced motion; both themes. |
| NAV.39 | Saved courses remember their unknown-space jumps and check them again | NAV.12, NAV.17 | With NAV.17. |

### Classes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.29 | Sweep every planet class for sense once the new ones are in (bug) | GEN.28, GEN.27, GEN.38 | Bug, but by definition a sweep after the new classes; it can't go earlier. |

### Pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.32 | Planet rows show the class only, without the type and moon labels | UX.29 | Same rows as UX.30; one thread. |
| UX.30 | Planet information without the Markdown render | UX.22, UX.32, UX.29 | Uses the unit ladders; shows the composition rows DB.2 now reads (PR #347). |
| UX.33 | Filter phenomena by their classes and types (bug) | GEN.28, GEN.47 | Filters over class lists that GEN.28 and GEN.47 change. |

### Backfill

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.42 | A pass that drops sectors from a region by probability | GEN.41, PERF.11 | Same function as PERF.18 (backfill_bright_stars_around); one thread. |
| GEN.43 | Don't over-filter: keep bright stars in odd places | GEN.42 |  |
| PERF.18 | Run the GEN.30 bright-star backfill in parallel on the work queue | GEN.44, PERF.19 | Same stars as the one-process backfill for one seed needs GEN.39. |
| GEN.40 | Weed out sectors by star density before the bright-star backfill | GEN.41, GEN.42, GEN.43, GEN.44 | Parent; closes with its subitems. |

### Queue

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.20 | Short-term caching through the work queue and API (needs planning) | PERF.19 | Plan with Boss. |

### API

| ID | Item | Needs | Note |
|---|---|---|---|
| API.5 | API version and compatibility checking | API.4 |  |
| API.10 | Reservations: claimed sectors and id blocks per run | API.9 | Reserved id blocks build on DB.3's id-block fix (PR #347). |
| API.11 | Staging tables | API.10 | Galaxy schema migration (staging); after NAV.10 in the writer queue. |
| API.12 | The download: seed, skeleton and name state | API.5 | Downloads the seed (what a seed means is GEN.39) and the name state (rules from GEN.46, done in PR #370). |
| API.16 | The API reports the galaxy's seed, version and run history | API.5 |  |

### 3D system

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.69 | A system scene endpoint with 3D orbits | NAV.7, MAP.57 | Scene endpoint with references. |
| MAP.70 | Positions at any time | MAP.69 | Python twin feeds NAV.27. |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.17 | The Generate page shows the galaxy's seed and version |  |  |
| OPS.15 | Each update says whether it changes generated output | OPS.13, GEN.58 | Needs the fingerprint, so phase 2. |

### Daily maintenance

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.61 | The daily merge folds pending admin changes into a new JSON file | GEN.59, ADM.18 | Boss 02:28Z: JSON changes only with the day's deltas. |
| OPS.18 | Settings JSON backups kept in 18 slots: 7 daily, 4 weekly, 6 monthly, 1 yearly | GEN.61 | Grandfather-father-son rotation; unit test with simulated dates. |
| OPS.16 | A daily maintenance script for Linux, macOS and Windows | GEN.61, OPS.18 | scripts/maintenance.sh and .ps1: positional update, delta merge, rotation; lock; optional OPS.15 check. |
| OPS.17 | Install and update set up the daily maintenance schedule | OPS.16, OPS.7, OPS.8, OPS.13 | Same scripts as OPS.7/OPS.8/OPS.13 (install/update, deploy-common), after them. |
| ADM.19 | The Admin dashboard lists the 18 settings backups for download | ADM.18, OPS.18 |  |

## Open questions for Boss

- MAP.56: Drop the 3x3 block pick: select a slab, zoom in, select a segment (bug), see its entry in TODO.md.
- MAP.53: Rotate a zoomed-in wedge, and zoom it to fit the window (bug), see its entry in TODO.md.
- MAP.58: Galaxy Map zoom limits: a short manual range on the galaxy wedge, locked below it, see its entry in TODO.md.
- MAP.54: Slab leader lines instead of the slab slider (bug), see its entry in TODO.md.
- NAV.24: A keep-out radius for every kind of object, see its entry in TODO.md.
- UX.32: Planet rows show the class only, without the type and moon labels, see its entry in TODO.md.
- UX.30: Planet information without the Markdown render, see its entry in TODO.md.
