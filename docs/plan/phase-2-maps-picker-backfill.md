# Phase 2: Classes, orbits, courses and recipes

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

The planet class refactor around the habitability index (with GEN.33, GEN.28, GEN.27 and GEN.29), nebula conditions on planets, n-body orbital updates and rogue collisions, editable trajectories, courses and waypoints, generate-by-recipe in the API, the backfill density pass, daily maintenance, the pilot-style visual design, light-travel positions, the asteroid-field and anomaly plans, and the API pieces remote generation needs first.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Classes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.33 | One class per PR, each with its tests | GEN.89 | Moved to phase 2 under GEN.90, after the habitability score. One class per PR (R and S first). Built on fixed physics so new classes aren't tuned to wrong masses, moons or zones. |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | GEN.33 | Under GEN.90; class Z is Boss's Earth-size world that never had life (2026-10-03). Class S landed with GEN.38 (PR #415); the other six classes here. Closes with GEN.33's PRs; PLANET_CLASSES in program_constants.py. |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | GEN.33 | Under GEN.90. Same reconcile/zone code as phase 0's physics fixes. |
| GEN.91 | Classes like S and V in the hot and cold zones | GEN.33, GEN.85 |  |
| GEN.92 | Life and its highest stage follow the habitability score | GEN.89, GEN.28 |  |
| GEN.29 | Sweep every planet class for sense once the new ones are in (bug) | GEN.28, GEN.27, GEN.91, GEN.92 | Under GEN.90: the refactor is the sweep. Bug, but by definition a sweep after the new classes; it can't go earlier. Includes rocky rogues of 10-16 Earth masses (up to 17,600 km) that get S as nearest fit (PR #415). |
| GEN.90 | Refactor the planet classes around the habitability index | GEN.33, GEN.28, GEN.27, GEN.29, GEN.91, GEN.92 | Parent; takes GEN.33, GEN.28, GEN.27 and GEN.29 as its subitems. |

### Galaxy Map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.58 | Galaxy Map zoom limits: a short manual range on the galaxy wedge, locked below it |  | A zoom policy of MAP.64. |
| MAP.75 | The mini map as a second engine view |  | Locked second camera on MAP.64. |
| MAP.59 | Make it plain that a zoomed-in slab is a slab, not a wedge | MAP.75 |  |

### Courses

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.20 | Draw the direct line and the route apart | NAV.13 |  |
| NAV.21 | Fit the view to the whole course | MAP.58 | Same fit as MAP.53, exempt from MAP.58's lock. |
| NAV.17 | A saved course record with both forms | NAV.7 |  |
| NAV.18 | Save, list, open, rename and delete, per browser | NAV.17 | Sibling of bookmarks.js. |
| NAV.4 | Save a course | NAV.17, NAV.18 | Per browser now; NAV.19 moves it into accounts in phase 3+. |
| NAV.24 | A keep-out radius for every kind of object |  | Unblocked: nebulae exist since GEN.47 (PR #419), so the nebula keep-out question can be settled. |
| NAV.36 | Unknown-space jumps drawn red and glowing | NAV.12, UX.35, NAV.20 | Route strip, navmap.py and the Galaxy Map course; legend; reduced motion; both themes. |
| NAV.39 | Saved courses remember their unknown-space jumps and check them again | NAV.12, NAV.17 | With NAV.17. |
| NAV.49 | Waypoints: pick objects in Star select mode and plot a course through them, kept on the map until cleared | MAP.122, NAV.16, NAV.17 | Merges "Plotted courses should appear on the galactic map and stay until cleared". |

### Pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.32 | Planet rows show the class only, without the type and moon labels |  | Same rows as UX.30; one thread. |
| UX.30 | Planet information without the Markdown render | UX.22, UX.32 | Uses the unit ladders; shows the composition rows DB.2 now reads (PR #347). |
| UX.43 | A visual design built like a pilot's starmap and navigation console | UX.42, UX.37 |  |

### Backfill

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.42 | A pass that drops sectors from a region by probability | GEN.41 | Same function as PERF.18 (backfill_bright_stars_around); one thread. |
| GEN.43 | Don't over-filter: keep bright stars in odd places | GEN.42 |  |
| PERF.18 | Run the GEN.30 bright-star backfill in parallel on the work queue |  | Backfill blocks become RQ jobs. Same stars as the one-process backfill for one seed needs GEN.39. |
| GEN.40 | Weed out sectors by star density before the bright-star backfill | GEN.41, GEN.42, GEN.43 | Parent; closes with its subitems. |

### Queue

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.20 | Short-term caching through the work queue and API (needs planning) |  | Plans on RQ results and the new caches. Plan with Boss. |

### API

| ID | Item | Needs | Note |
|---|---|---|---|
| API.5 | API version and compatibility checking | API.4 |  |
| API.10 | Reservations: claimed sectors and id blocks per run | API.9 | Reserved id blocks build on DB.3's id-block fix (PR #347). |
| API.11 | Staging tables | API.10 | Galaxy schema migration (staging); after NAV.10 in the writer queue. |
| API.12 | The download: seed, skeleton and name state | API.5 | Downloads the naming key, not name registries. Downloads the seed (what a seed means is GEN.39) and the name state (rules from GEN.46, done in PR #370). |
| API.16 | The API reports the galaxy's seed, version and run history | API.5 |  |

### 3D system

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.69 | A system scene endpoint with 3D orbits | NAV.7 | Scene endpoint with references. |
| MAP.70 | Positions at any time | MAP.69 | Positions through the point-in-space object (GEN.74). Python twin feeds NAV.27. |

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
| OPS.17 | Install and update set up the daily maintenance schedule | OPS.16, OPS.8, OPS.13 | Same scripts as OPS.7/OPS.8/OPS.13 (install/update, deploy-common), after them. |
| ADM.19 | The Admin dashboard lists the 18 settings backups for download | ADM.18, OPS.18 |  |

### Picker

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.45 | "What's within N pc" from the Galaxy Map and Sector Map | NAV.44, MAP.65, NAV.15 | Boss 04:39Z. Map action opening NAV.44; drawn sphere optional. |

### Bookmarks

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.48 | Charted regions as bookmarks that frame and outline the region | UX.47, MAP.65 |  |
| UX.45 | Bookmark management | UX.46, UX.47, UX.48 | Parent; per browser until accounts (USR.7) move them. |

### Orbital updates

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.115 | The galaxy's own gravity: a smooth disk, bulge and halo potential | GEN.106, GEN.108 | Model from Boss's "Computational Astrodynamics.md" (2026-10-07): bulge, disk and halo potential. |
| GEN.109 | N-body influence from the nearest 10 bodies of equal or larger mass, with a Hill-radius warning | GEN.106, GEN.108, GEN.115 |  |
| ADM.36 | Change an object's trajectory vector | GEN.74, GEN.109 |  |
| GEN.110 | Rogue planet collisions: asteroid fields, merged giants and new stars | GEN.109 |  |
| GEN.105 | Orbital updates | GEN.106, GEN.107, GEN.108, GEN.115, GEN.109, GEN.110 | Parent of the orbital update work. |

### Nebula planets

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.95 | Nebula conditions applied when planets and surfaces are generated | GEN.94, GEN.89, GEN.75 |  |
| GEN.93 | Nebula conditions in planet generation | GEN.94, GEN.95 | Parent. |

### Asteroid fields

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.112 | Plan asteroid fields and belts as object systems for rendering |  | Plan only. |

### Anomalies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.113 | Analyze the anomaly docs: which anomalies to add and how |  | Boss: "probably Phase 2 and 3 just to analyze those". |

### Maps

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.121 | Every map shows and steps to its neighbouring regions, on one map engine | MAP.66, MAP.102 | Folds MAP.127, the Galaxy Map's neighbouring blocks and slabs (Boss 2026-10-07 16:27Z: one map engine for all maps). |

### Recipes

| ID | Item | Needs | Note |
|---|---|---|---|
| API.18 | Generate by recipe: JSON for sectors, systems, planets, moons and phenomena | ADM.21, GEN.96, API.9 | Validation through the Pydantic models; 400 for nonsense, 422 for validation failures. |

### View

| ID | Item | Needs | Note |
|---|---|---|---|
| VIEW.5 | Light-travel positions: where an object appears to a distant observer | GEN.74, MAP.70 | Groundwork for VIEW.2. |

## Open questions for Boss

- MAP.58: Galaxy Map zoom limits: a short manual range on the galaxy wedge, locked below it, see its entry in TODO.md.
- NAV.24: A keep-out radius for every kind of object, see its entry in TODO.md.
- UX.32: Planet rows show the class only, without the type and moon labels, see its entry in TODO.md.
- UX.30: Planet information without the Markdown render, see its entry in TODO.md.
