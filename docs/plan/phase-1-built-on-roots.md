# Phase 1: Built on the roots

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

Work that needs phase 0 in place: galaxy generation on the parallel path with the per-sector stats table (GEN.44 and PERF.11), prevalence controls, planet classes, reproducible galaxies up to the golden-seed test, sector colors, routing, the first picker pieces and the queue.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Classes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.33 | One class per PR, each with its tests | GEN.34, GEN.35, GEN.36, GEN.25, GEN.37 | One class per PR (R and S first). Built on fixed physics so new classes aren't tuned to wrong masses, moons or zones. |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | GEN.33 | Closes with GEN.33's PRs; PLANET_CLASSES in program_constants.py. |
| GEN.38 | Rocky rogue planets over 10,000 km are still classed C (bug) | GEN.33, GEN.45 | Probably solved by a rogue-eligible S class. |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | GEN.33, GEN.37 | Same reconcile/zone code as phase 0's physics fixes. |

### Prevalence

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.51 | Forcing options only for single-system generation | GEN.50 | generate.py sector/galaxy argument parsing. |
| GEN.52 | Prevalence controls for sector and galaxy runs | GEN.51 | Probability adjustments reach systemData.py, the same constructor phase 0 thread D fixed. |
| TEST.75 | Tests for forcing and prevalence | GEN.52 | Grows with GEN.49 to GEN.52. |
| ADM.16 | Prevalence controls on the Generate page | GEN.52, ADM.14 | generate.html, after ADM.14's layout. |
| GEN.48 | Forcing options are impractical for whole sectors; replace them with prevalence controls (bug) | GEN.49, GEN.50, GEN.51, GEN.52, ADM.16, TEST.75 | Parent; closes with its subitems. |

### Galaxy gen

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.24 | Generate the galactic core on layer 0 | GEN.31, PERF.21, ADM.14 | Bulk core fill runs on the parallel path; new mode on generate.html. |
| GEN.44 | Store each sector's backfill level so finished sectors drop out of any backfill | PERF.21, GEN.32 | One shared per-sector stats table with PERF.11 (Boss 01:46Z). Galaxy schema v51 (one writer at a time). Backfill code shared with PERF.18 and GEN.42. |
| PERF.11 | Store each sector's expected and actual density | GEN.44 | Same per-sector stats table as GEN.44 (Boss 01:46Z); MAP.86's color goes there too. |
| PERF.1 | Generation at scale (new parent) | PERF.11 | Parent; only PERF.11 is open under it. |
| GEN.41 | Investigate: how much backfill work a density pre-pass would save | GEN.44 | Investigation; go/no-go for GEN.42. |
| GEN.47 | Nebulae almost never appear (bug) | PERF.21, GEN.39 | Galaxy-scale nebula field spanning sectors: every worker and every later run must agree where a cloud is, so it needs deterministic per-region draws (GEN.39, or an address hash if GEN.39 is dropped). |

### System Map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.89 | System Map: space orbits with a fitted scale and a minimum ring gap instead of plain log | MAP.88 | Same file as MAP.88; also changes systemmap.js kmToPx (MAP.63 touched it too). |

### Galaxy Map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.86 | Sector and block colors from what is in them: filled sectors translucent (bug) | MAP.85, PERF.11 | Next after MAP.85: with no lines, color carries the structure. Its color goes in the shared per-sector stats table (GEN.44 + PERF.11). Tile payload change: bump the tile cache. |

### Galaxy tiles

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.80 | Sector-level zoom on the Galaxy Map should show almost every star in the sector (bug) | MAP.90, MAP.87 | Judgment: moved up from the selection chain; the thinning is in the tile listing (queryDb GALAXY_TILE_* floors) and galaxymap3d.js, not the pick code. |

### References

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.8 | Pages and anchors for stars, planets, moons and belts | NAV.7 | System page anchors (system.html). |
| NAV.9 | Search and locate return references for every kind | NAV.7 | queryDb search and galaxy_locate. |

### Picker

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.13 | A picker module: select, step out, step in, step sideways | NAV.7 | picker.js; needs no engine. |
| NAV.14 | One breadcrumb for every level | NAV.13 |  |

### Routing

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.10 | Routing that scales past a few thousand systems | NAV.12 | Galaxy schema migration for position indexes; queue behind PERF.11. Must keep NAV.12's guarantee (a route always exists). |
| NAV.11 | Travel times for the system-to-system route too | NAV.10, NAV.35 | Times per hop, including unknown-space jumps. |

### Pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.22 | Meaningful units for every measurement | UX.23 | One quantity family per PR. |
| UX.25 | Rogue planets: octant and a small map symbol beside each name (bug) | UX.24, UX.28 | Same sector table as UX.24. |
| UX.26 | Edit and admin actions as a button that opens a menu (bug) | UX.28, UX.25 | Sector page admin panel and edit_controls.html. |
| UX.31 | Editing a star system: an edit button with a quick menu, not a long panel (bug) | UX.26 | system.html edit panel (_edit_rows in system_pages.py). |
| UX.27 | System page: the system and navigation buttons on one row that doesn't overlap (bug) | UX.28 | system.html subhead; shares wording with NAV.29. |
| UX.3 | Warn visitors while a background job changes the galaxy | PERF.21, PERF.23 | ETA from progress.json, which PERF.23 caps. |

### Queue

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.19 | Everything the API or web site starts runs on the work queue (investigate) | PERF.21 | Audit only; nothing moves to the queue until it runs at any worker count. |
| ADM.15 | Change the worker count from the Queue page, with a "Ludicrous Speed" mode | PERF.21, TEST.73 | Changing the worker count live only makes sense once any count works. |

### API

| ID | Item | Needs | Note |
|---|---|---|---|
| API.4 | API compatibility data in the docs |  | Docs and version number. |
| API.7 | Investigate and plan upload limits |  | Plan only. |
| API.9 | Key scopes |  | Control schema migration (v8). Judgment: if user keys belong to accounts (API.6's open question), USR.2's table design comes first. |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| OPS.11 | Define "the same galaxy" and which versions stay reproducible |  | Design note. |
| GEN.56 | Every random draw in generation comes from the derived seeds | GEN.39 | Touches every generator module. |
| GEN.57 | A sector's contents depend only on the seed, the version and its address | PERF.21, GEN.39, GEN.56, GEN.46, GEN.44 |  |
| DB.7 | The version that generated each sector, and a warning for mixed-version galaxies | DB.6 |  |
| GEN.58 | A fingerprint of a galaxy's generated content | OPS.11 | Judgment: phase 1 so the golden test guards later changes. |
| TEST.77 | A golden-seed regression test | GEN.57, GEN.58 |  |

### Route display

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.36 | Unknown-space jumps drawn in glowing red | NAV.35 | NAV map now; course drawings later. |
| UX.35 | The route shown horizontally, wrapping onto several lines on narrow screens | NAV.12 | NAV page layout. |

## Open questions for Boss

- PERF.11: Store each sector's expected and actual density, see its entry in TODO.md.
- NAV.8: Pages and anchors for stars, planets, moons and belts, see its entry in TODO.md.
- UX.22: Meaningful units for every measurement, see its entry in TODO.md.
- UX.3: Warn visitors while a background job changes the galaxy, see its entry in TODO.md.
- API.4: API compatibility data in the docs, see its entry in TODO.md.
