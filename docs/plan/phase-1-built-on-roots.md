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
| NAV.8 | Pages and anchors for stars, planets, moons and belts |  | System page anchors (system.html). |
| NAV.9 | Search and locate return references for every kind |  | queryDb search and galaxy_locate. |

### Picker

| ID | Item | Needs | Note |
|---|---|---|---|

### Database consistency check

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.8 | Check a galaxy database and say whether it is damaged |  | Schema check reads Alembic's revision; the name-registry check goes with GEN.71. 9 after it. Boss 02:13Z: phase 0, its own thread. Read-only. Stats and version checks switch on once GEN.44/PERF.11 and DB.6/DB.7 land; no hard dependency. |
| DB.9 | Repair a damaged galaxy database from a parity file | DB.8, GEN.57, OPS.14 | Boss 02:13Z: phase 1. Reed-Solomon parity over groups of sector exports; seed regeneration as fallback. |

### Generation

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.29 | Record which runs a partly filled sector still needs |  | Builds on the generated flag that GEN.76 settles (sector_stats.bright_level_sol). |
| PERF.30 | Finish an interrupted block or sector run on the next start | PERF.29 |  |

### Generate page

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.29 | Fill a span of layers, rings or columns |  |  |
| ADM.30 | Radial generation: a cylinder of N sectors around a point |  |  |
| ADM.31 | Every generate action offers to show what it made on the Galaxy Map |  | Merges two asks: the 2026-10-03 button and the 2026-10-07 "see that space". |
| GEN.96 | Generation directives for a sector (an override button) |  | A subset of what API recipes (API.18) later take. |
| GEN.97 | Generate N random neighborhoods |  |  |
| ADM.45 | Prevalence fields take the override share directly and must total 100% |  | Boss 2026-10-09 07:48Z; follows ADM.37. |
| ADM.28 | A simpler Generate page: layer specs, a Customize window and plain controls | ADM.29, ADM.30, ADM.31, GEN.97 | Parent of the Generate page items; after the web components land. |

### Galaxy gen

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.24 | Generate the galactic core on layer 0 |  | Bulk core fill runs on the parallel path; new mode on generate.html. |
| GEN.41 | Investigate: how much backfill work a density pre-pass would save |  | Investigation; go/no-go for GEN.42. |
| GEN.99 | Nebula volume backfill with the star types the nebula needs |  |  |
| GEN.101 | Fill order: nearest sectors first along a pruned Hilbert octree curve |  |  |
| GEN.102 | Investigate filling all near-zero-density void space at once |  |  |
| GEN.103 | Research where each star type and phenomenon belongs in the galaxy's structure |  |  |

### Habitability

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.85 | Atmosphere species, partial pressures and mantle redox for every planet |  |  |
| GEN.86 | Stellar activity (XUV, flares) and planetary magnetic fields |  |  |
| GEN.87 | Surface radiation dose | GEN.85, GEN.86 |  |
| GEN.88 | Hydrosphere and ocean chemistry | GEN.85 | Reuses rogueSurface's ice-shell and ocean functions. |
| GEN.89 | The habitability score for every planet and moon | GEN.85, GEN.86, GEN.87, GEN.88 |  |
| GEN.83 | A planetary habitability index (PHI) | GEN.85, GEN.86, GEN.87, GEN.88, GEN.89 | Parent; the class refactor (GEN item GEN.90) follows it. |

### Tech levels and facilities

| ID | Item | Needs | Note |
|---|---|---|---|
| POP.8 | A tech-level design from the six domain indices |  |  |
| POP.9 | A tech level generated for every technological species | POP.8 |  |
| POP.7 | Tech levels for technological species | POP.8, POP.9 | Parent. |
| POP.10 | Facility types: programmable, picked from a dropdown, with affiliation and Green/Yellow/Red ratings |  |  |

### Nebula planets

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.94 | Feasibility study: can planets form in each nebula class, and what changes |  | None of the habitability docs cover nebulae; this is new research. |

### Orbital updates

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.104 | A spin vector and a realistic axial tilt for every rotating object |  | Rules from "Observational Kinetics for Rotational Vectors.md". |
| GEN.107 | The update reports how many objects moved, changed sector, or entered or left a nebula |  |  |
| GEN.108 | Orbital math limits: where each method breaks down and what happens there |  | Boss (2026-10-07 11:47Z): "we need to ensure that reasonable limitations for when the math breaks down at the edge cases." |

### Routing

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.10 | Routing that scales past a few thousand systems |  | Its position indexes are an Alembic migration. Galaxy schema migration for position indexes (v53 taken by PR #425). Built with NAV.12 (no hop limit, a route always exists). |
| NAV.12 | No maximum hop length: a route always reaches the nearest star it can, across any number of sectors | NAV.10 | Anchor (Boss 01:53Z game mechanic). Built with NAV.10. Per-hop unknown-space flag in /api/nav. |
| UX.35 | The NAV page route shown horizontally, wrapping onto several lines on narrow screens |  | In parallel with NAV.12; replaces ol.nav-route in nav.html. |
| NAV.11 | Travel times for the system-to-system route too | NAV.10, NAV.12 | Times per hop, including unknown-space jumps; total assumes a stop at every system (Boss 04:19Z); open question on a stay per stop. |
| NAV.42 | Each route stop shows the course and distance to the next stop | UX.35 | Boss 04:19Z. format_course per hop, frame per pair. |
| NAV.47 | Unknown-space jumps stop at scattered stars, black holes, neutron stars and quasars | NAV.12 |  |
| NAV.48 | Offer to generate the uncharted sectors that block a course | NAV.12 |  |

### Nearby search

| ID | Item | Needs | Note |
|---|---|---|---|

### Maps

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.119 | Expected star density editable by admins on the Galaxy Map |  |  |
| MAP.122 | A Select mode on every galaxy view: Galaxy (blocks and sectors) or Star |  | Merges the 2026-10-07 "button to select a star" ask; MAP.101 (PR #431) made stars unpickable by default. |

### System Map

| ID | Item | Needs | Note |
|---|---|---|---|

### Bookmarks

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.46 | Wiping the galaxy wipes the bookmarks |  |  |
| UX.47 | A bookmark manager: list, rename, sort, group and delete | UX.46 |  |
| UX.49 | Form fields as Shoelace components (sl-input, sl-select, sl-checkbox) across the site |  | Follows UX.40 (done): the fields themselves become Shoelace components. |

### Admin control

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.32 | Add a star system to a sector: at the emptiest spot, at given coordinates, or at random outside every Hill sphere |  | Merges two asks (one system added by the computer; placement chosen by admin). |
| MAP.120 | Bright-star backfill from the Galaxy Map's block, slab and wedge menus |  |  |
| ADM.35 | Full control from every screen: edit anything, regenerate with every input, backfill or erase what is in view | MAP.120 | Parent; MAP item MAP.120 is the map part. |

### Pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.23 | A shared unit-ladder module |  | The ladder uses astropy.units. 22 (a feature); UX.36 keeps today's formatter and needs no ladder. Shared unit ladder; UX.22 then UX.30 build on it. |
| UX.22 | Meaningful units for every measurement | UX.23 | One quantity family per PR. |
| UX.3 | Warn every visitor while a background job changes the galaxy |  | ETA from the RQ job's published progress. ETA from progress.json, which PERF.23 caps. |
| UX.42 | In-universe wording across the interface |  | After the UX sweep removes controls, so only kept wording is changed. |

### Queue

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.15 | Change the worker count from the Queue page, with a "Ludicrous Speed" mode |  | Worker count now means RQ workers. Changing the worker count live only makes sense once any count works. |

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
| GEN.57 | A sector's contents depend only on the seed, the version and its address |  | Name collisions go away with GEN.67; GEN.63 dropped. Also keys GEN.64's ID collision counter on address (it checks stored rows, PR #406). |
| DB.7 | The version that generated each sector, and a warning for mixed-version galaxies |  | An Alembic migration. |
| TEST.77 | A golden-seed regression test | GEN.57 |  |
| OPS.14 | A warning when the running version key differs from the galaxy's | DB.7 | Feeds GEN.58's output and OPS.12. |
| ADM.18 | The galaxy's creation settings saved as a JSON file, downloadable from the Admin dashboard | DB.7 | Stores the naming key instead of the word list. Boss 02:13Z: phase 1. Includes the key history; dated backup on every change. |
| GEN.59 | Admin changes stored as a net difference from the generated galaxy | ADM.18 | Boss 02:20Z: net diff by stable address path, regenerate seeds, in ADM.18's JSON. |

### Bugs from the GitHub issues

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.34 | The site stays responsive during heavy generation jobs (bug) | PERF.31 | GitHub issues [#638](https://github.com/dwhagar/planetGen/issues/638) and [#614](https://github.com/dwhagar/planetGen/issues/614) (second half). |

### Foundations for the issue features

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.31 | Investigate: where generation spends its time, from the plan to a finished galaxy |  | GitHub issues [#761](https://github.com/dwhagar/planetGen/issues/761) and [#750](https://github.com/dwhagar/planetGen/issues/750) (benchmark half). |
| PERF.32 | Generation performance stats: rates recorded per run, deleted on every new version |  | GitHub issues [#661](https://github.com/dwhagar/planetGen/issues/661) (store half) and #750. |
| GEN.128 | Design: multi-star hierarchies and compact-object primaries |  | GitHub issues [#777](https://github.com/dwhagar/planetGen/issues/777) and [#778](https://github.com/dwhagar/planetGen/issues/778): the research both builds wait on. |
| ADM.42 | One settings model describes every config.json option |  | Foundation for GitHub issues [#515](https://github.com/dwhagar/planetGen/issues/515) and [#743](https://github.com/dwhagar/planetGen/issues/743). |

## Open questions for Boss

- NAV.8: Pages and anchors for stars, planets, moons and belts, see its entry in TODO.md.
- NAV.11: Travel times for the system-to-system route too, see its entry in TODO.md.
- UX.22: Meaningful units for every measurement, see its entry in TODO.md.
- UX.3: Warn every visitor while a background job changes the galaxy, see its entry in TODO.md.
- API.4: API compatibility data in the docs, see its entry in TODO.md.
- OPS.13: Every update records the version key, keeping the last 10, see its entry in TODO.md.
