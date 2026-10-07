# Phase 1: Built on the groundwork

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

Object references and the picker finished, the database check and repair, resumable runs, the reworked Generate page (spans, radial fills, random neighborhoods, directives, show on the map), galaxy generation changes (phenomena placed galaxy-wide first, fill order, backfill from the run's edge, nebula volume backfill, star-type tiers), the habitability index, tech levels and facility types, spin and orbital-update thresholds with their limits, routing with unknown-space stops and charting along a course, the nearby search, Select mode and view filters on the maps, bookmarks, admin control from every screen, in-universe wording, the API log and scopes, and reproducible galaxies up to the golden-seed test.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### References

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.8 | Pages and anchors for stars, planets, moons and belts | NAV.7 | System page anchors (system.html). |
| NAV.9 | Search and locate return references for every kind | NAV.7 | queryDb search and galaxy_locate. |
| NAV.16 | NAV endpoints can be any object | NAV.7 | Moved from phase 2: NAV.3 closes in phase 1 now. navigation.py legs, nav_page.py endpoints. |

### Picker

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.3 | One shared picker for the Galaxy, Sector and System displays | NAV.13, NAV.14, NAV.15, NAV.16, NAV.32 | Moved from phase 3: its parts are in phases 0 and 1. Parent; closes with its subitems. |

### Database consistency check

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.8 | Check a galaxy database and say whether it is damaged | DB.11 | Schema check reads Alembic's revision; the name-registry check goes with GEN.71. 9 after it. Boss 02:13Z: phase 0, its own thread. Read-only. Stats and version checks switch on once GEN.44/PERF.11 and DB.6/DB.7 land; no hard dependency. |
| DB.9 | Repair a damaged galaxy database from a parity file | DB.8, GEN.57, GEN.58, OPS.14 | Boss 02:13Z: phase 1. Reed-Solomon parity over groups of sector exports; seed regeneration as fallback. |

### Generation

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.29 | Record which runs a partly filled sector still needs | GEN.76 | Builds on the generated flag that GEN.76 settles (sector_stats.bright_level_sol). |
| PERF.30 | Finish an interrupted block or sector run on the next start | PERF.29 |  |

### Generate page

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.29 | Fill a span of layers, rings or columns |  |  |
| ADM.30 | Radial generation: a cylinder of N sectors around a point |  |  |
| ADM.31 | Every generate action offers to show what it made on the Galaxy Map |  | Merges two asks: the 2026-10-03 button and the 2026-10-07 "see that space". |
| GEN.96 | Generation directives for a sector (an override button) | GEN.52 | A subset of what API recipes (API.18) later take. |
| GEN.97 | Generate N random neighborhoods |  |  |
| ADM.28 | A simpler Generate page: layer specs, a Customize window and plain controls | UX.40, ADM.29, ADM.30, ADM.31, GEN.97 | Parent of the Generate page items; after the web components land. |

### Galaxy gen

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.24 | Generate the galactic core on layer 0 | ADM.14 | Bulk core fill runs on the parallel path; new mode on generate.html. |
| GEN.41 | Investigate: how much backfill work a density pre-pass would save |  | Investigation; go/no-go for GEN.42. |
| GEN.98 | Bright-star backfill from the farthest generated boundary outward |  |  |
| GEN.100 | Phenomena placed galaxy-wide first and kept when sectors fill |  | Merges the 2026-10-03 scatter-order item and the 2026-10-07 "generated through the entire galaxy first". |
| GEN.99 | Nebula volume backfill with the star types the nebula needs | GEN.75, GEN.100 |  |
| GEN.101 | Fill order: nearest sectors first along a pruned Hilbert octree curve |  |  |
| GEN.102 | Investigate filling all near-zero-density void space at once | GEN.78 |  |
| GEN.103 | Research where each star type and phenomenon belongs in the galaxy's structure |  |  |

### Habitability

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.84 | Habitability design: one score structure and reconciled thresholds |  | Decision for Boss inside it: which score structure. |
| GEN.85 | Atmosphere species, partial pressures and mantle redox for every planet | GEN.84 |  |
| GEN.86 | Stellar activity (XUV, flares) and planetary magnetic fields | GEN.84 |  |
| GEN.87 | Surface radiation dose | GEN.85, GEN.86 |  |
| GEN.88 | Hydrosphere and ocean chemistry | GEN.85 | Reuses rogueSurface's ice-shell and ocean functions. |
| GEN.89 | The habitability score for every planet and moon | GEN.85, GEN.86, GEN.87, GEN.88 |  |
| GEN.83 | A planetary habitability index (PHI) | GEN.84, GEN.85, GEN.86, GEN.87, GEN.88, GEN.89 | Parent; the class refactor (GEN item GEN.90) follows it. |

### Tech levels and facilities

| ID | Item | Needs | Note |
|---|---|---|---|
| POP.8 | A tech-level design from the six domain indices |  |  |
| POP.9 | A tech level generated for every technological species | POP.8, GEN.80 |  |
| POP.7 | Tech levels for technological species | POP.8, POP.9 | Parent. |
| POP.10 | Facility types: programmable, picked from a dropdown, with affiliation and Green/Yellow/Red ratings |  |  |

### Nebula planets

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.94 | Feasibility study: can planets form in each nebula class, and what changes |  | None of the habitability docs cover nebulae; this is new research. |

### Orbital updates

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.104 | A spin vector and a realistic axial tilt for every rotating object | GEN.74 | Rules from "Observational Kinetics for Rotational Vectors.md". |
| GEN.106 | Movement thresholds and a next-update-due column | GEN.74, DB.11 |  |
| GEN.107 | The update reports how many objects moved, changed sector, or entered or left a nebula | GEN.106, GEN.75 |  |
| GEN.108 | Orbital math limits: where each method breaks down and what happens there | GEN.106 | Boss (2026-10-07 11:47Z): "we need to ensure that reasonable limitations for when the math breaks down at the edge cases." |

### Routing

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.10 | Routing that scales past a few thousand systems | DB.11 | Its position indexes are an Alembic migration. Galaxy schema migration for position indexes (v53 taken by PR #425). Built with NAV.12 (no hop limit, a route always exists). |
| NAV.12 | No maximum hop length: a route always reaches the nearest star it can, across any number of sectors | NAV.10 | Anchor (Boss 01:53Z game mechanic). Built with NAV.10. Per-hop unknown-space flag in /api/nav. |
| UX.35 | The NAV page route shown horizontally, wrapping onto several lines on narrow screens |  | In parallel with NAV.12; replaces ol.nav-route in nav.html. |
| NAV.11 | Travel times for the system-to-system route too | NAV.10, NAV.12 | Times per hop, including unknown-space jumps; total assumes a stop at every system (Boss 04:19Z); open question on a stay per stop. |
| NAV.42 | Each route stop shows the course and distance to the next stop | UX.35 | Boss 04:19Z. format_course per hop, frame per pair. |
| NAV.47 | Unknown-space jumps stop at scattered stars, black holes, neutron stars and quasars | NAV.12, GEN.100 |  |
| NAV.48 | Offer to generate the uncharted sectors that block a course | NAV.12 |  |

### Nearby search

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.43 | Find everything within a distance of a place: the query and the API | NAV.7 | Boss 04:39Z. Replaces systems_within_radius (one sector, systems only); enumerate_sectors_within_radius then per-sector reads; open questions: max distance, ungenerated sectors, generated only. |
| NAV.44 | A "What's nearby" page: pick a place, enter a distance in parsecs, list what is there | NAV.43 | Boss 04:39Z. Page with place picker, distance in pc, kind filters, 50-row pages. |

### Maps

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.119 | Expected star density editable by admins on the Galaxy Map | MAP.65 |  |
| MAP.122 | A Select mode on every galaxy view: Galaxy (blocks and sectors) or Star | MAP.65 | Merges the 2026-10-07 "button to select a star" ask; MAP.101 (PR #431) made stars unpickable by default. |
| MAP.123 | Show or hide star types and phenomena, and set the luminosity floor, on the Galaxy and Sector Maps | MAP.65, MAP.79 | Extends MAP.79's per-kind buttons to star types and the Galaxy Map. |
| MAP.124 | The Galaxy Map opens zoomed to fit all charted space | MAP.67 |  |

### System Map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.89 | System Map: space orbits with a fitted scale and a minimum ring gap instead of plain log |  | Same file as MAP.88 (done, PR #405); also changes systemmap.js kmToPx (MAP.63 touched it too). |

### Bookmarks

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.46 | Wiping the galaxy wipes the bookmarks |  |  |
| UX.47 | A bookmark manager: list, rename, sort, group and delete | UX.46 |  |

### Admin control

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.32 | Add a star system to a sector: at the emptiest spot, at given coordinates, or at random outside every Hill sphere | GEN.74 | Merges two asks (one system added by the computer; placement chosen by admin). |
| MAP.120 | Bright-star backfill from the Galaxy Map's block, slab and wedge menus | MAP.65 |  |
| ADM.35 | Full control from every screen: edit anything, regenerate with every input, backfill or erase what is in view | ADM.34, MAP.120 | Parent; MAP item MAP.120 is the map part. |

### Pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.23 | A shared unit-ladder module | GEN.66 | The ladder uses astropy.units. 22 (a feature); UX.36 keeps today's formatter and needs no ladder. Shared unit ladder; UX.22 then UX.30 build on it. |
| UX.22 | Meaningful units for every measurement | UX.23 | One quantity family per PR. |
| UX.3 | Warn every visitor while a background job changes the galaxy | ADM.22 | ETA from the RQ job's published progress. ETA from progress.json, which PERF.23 caps. |
| UX.42 | In-universe wording across the interface | UX.37 | After the UX sweep removes controls, so only kept wording is changed. |

### Queue

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.15 | Change the worker count from the Queue page, with a "Ludicrous Speed" mode | PERF.24 | Worker count now means RQ workers. Changing the worker count live only makes sense once any count works. |

### API

| ID | Item | Needs | Note |
|---|---|---|---|
| API.15 | Log every API call with its user, how it came in, and its HTTP response code |  | Needs no user accounts (Boss 01:31Z). |
| API.4 | API compatibility data in the docs |  | Docs and version number. |
| API.7 | Investigate and plan upload limits |  | Plan only. |
| API.9 | Key scopes |  | Control schema migration (v8). Decided: user keys belong to accounts, so API.6 waits for USR.2 (phase 3+); API.9's scopes don't. |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.56 | Every random draw in generation comes from the derived seeds |  | Touches every generator module. |
| GEN.57 | A sector's contents depend only on the seed, the version and its address | GEN.56 | Name collisions go away with GEN.67; GEN.63 dropped. Also keys GEN.64's ID collision counter on address (it checks stored rows, PR #406). |
| DB.7 | The version that generated each sector, and a warning for mixed-version galaxies | DB.11 | An Alembic migration. |
| GEN.58 | A fingerprint of a galaxy's generated content |  | Judgment: phase 1 so the golden test guards later changes. |
| TEST.77 | A golden-seed regression test | GEN.57, GEN.58 |  |
| OPS.8 | Update reloads Apache itself when run as root | OPS.7 | 13 in the same update scripts. Not a bug, but the same files as OPS.7, so it rides along. |
| OPS.13 | Every update records the version key, keeping the last 10 | OPS.7, OPS.8 | No corpus or name-list hashes once GEN.71 lands; the lock hashes stay. update.sh / update.ps1 after OPS.7 and OPS.8; control-database history table. Open question on "recalculate the seed value". |
| OPS.14 | A warning when the running version key differs from the galaxy's | DB.7, OPS.13 | Feeds GEN.58's output and OPS.12. |
| ADM.18 | The galaxy's creation settings saved as a JSON file, downloadable from the Admin dashboard | DB.7, OPS.13, GEN.70 | Stores the naming key instead of the word list. Boss 02:13Z: phase 1. Includes the key history; dated backup on every change. |
| GEN.59 | Admin changes stored as a net difference from the generated galaxy | GEN.56, GEN.58, ADM.18 | Boss 02:20Z: net diff by stable address path, regenerate seeds, in ADM.18's JSON. |

## Open questions for Boss

- NAV.8: Pages and anchors for stars, planets, moons and belts, see its entry in TODO.md.
- NAV.11: Travel times for the system-to-system route too, see its entry in TODO.md.
- NAV.43: Find everything within a distance of a place: the query and the API, see its entry in TODO.md.
- UX.22: Meaningful units for every measurement, see its entry in TODO.md.
- UX.3: Warn every visitor while a background job changes the galaxy, see its entry in TODO.md.
- API.4: API compatibility data in the docs, see its entry in TODO.md.
- GEN.57: A sector's contents depend only on the seed, the version and its address, see its entry in TODO.md.
- OPS.13: Every update records the version key, keeping the last 10, see its entry in TODO.md.
