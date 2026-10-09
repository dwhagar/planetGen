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
| DB.16 | Store the generator epoch and run id on each sector instead of four version text columns | OPS.28 | Research follow-up to DB.7 (built). |
| DOC.5 | Rewrite the object ID docs: object-ids.md, database-schema.md and api.md | GEN.170 | Object-ID research. |
| TEST.110 | Object ID tests: identical IDs on 1 and 4 workers, none reused, none missing | GEN.171, GEN.172, GEN.176 | Object-ID research. |
| GEN.176 | A nebula or remnant is born in the sector holding the centre of the space it occupies | DB.20 | Object-ID research. |
| GEN.175 | Regenerating a phenomenon sets its uid to NULL (bug) |  | Object-ID research. The one-line keep can go first. |
| GEN.174 | Bodies an admin adds are saved with a NULL uid (bug) | GEN.172 | Object-ID research. Closed by the run-time birth item. |
| GEN.173 | Deleting a body and then adding one fails with IntegrityError 1062 on uq_planets_uid (bug) | GEN.171, GEN.172 | Object-ID research. Closed by the fill and run-time birth items. |
| GEN.172 | Run-time births get object IDs from the counters | DB.20 | Object-ID research. |
| GEN.171 | The sector fill gives object IDs by generation rank | DB.20 | Object-ID research. |
| DB.20 | Object IDs in the schema: uid becomes BINARY(10), unique on its own, plus an id_counters table | GEN.170 | Object-ID research. |
| GEN.170 | Object ID layout: an 80-bit ID of birth sector, serial and body number, with pack, unpack, format and parse functions |  | Object-ID research. First of the object-ID items; nothing built until Boss asks. |
| DB.9 | Repair a damaged galaxy database from a parity file | DB.8, GEN.57, OPS.14 | Boss 02:13Z: phase 1. Reed-Solomon parity over groups of sector exports; seed regeneration as fallback. |
| DB.17 | Repair by regenerating a damaged sector from its seed when parity cannot rebuild it | DB.9, GEN.57, OPS.14 | Research split of DB.9: the regenerate-from-seed fallback. |

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
| GEN.155 | A nuclear-cluster object for the Sgr A* sector (optional) |  | Research: optional. |
| GEN.99 | Nebula volume backfill with the star types the nebula needs |  |  |
| GEN.150 | H II region radius and density from the ionizing photon rate, and IMF-based nebula hosts |  | Research: sub-item of GEN.99. |
| GEN.101 | Fill order: nearest sectors first along a pruned Hilbert octree curve |  |  |
| GEN.102 | Investigate filling all near-zero-density void space at once |  |  |
| GEN.103 | Research where each star type and phenomenon belongs in the galaxy's structure |  |  |

### Habitability

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.148 | Habitability index follow-ups from the research (GEN.84 built) |  | Research follow-up to GEN.84 (built). |
| GEN.86 | Stellar activity (XUV, flares) and planetary magnetic fields |  |  |
| GEN.87 | Surface radiation dose | GEN.86 |  |
| GEN.88 | Hydrosphere and ocean chemistry |  | Reuses rogueSurface's ice-shell and ocean functions. |
| GEN.89 | The habitability score for every planet and moon | GEN.86, GEN.87, GEN.88 |  |
| GEN.83 | A planetary habitability index (PHI) | GEN.86, GEN.87, GEN.88, GEN.89 | Parent; the class refactor (GEN item GEN.90) follows it. |

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

### Orbital updates

| ID | Item | Needs | Note |
|---|---|---|---|

### Routing

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.52 | Port `join_islands` and the k-d tree to cKDTree |  | Research: performance cliff in the built router. |
| UX.35 | The NAV page route shown horizontally, wrapping onto several lines on narrow screens |  | In parallel with NAV.12; replaces ol.nav-route in nav.html. |
| NAV.11 | Travel times for the system-to-system route too |  | Times per hop, including unknown-space jumps; total assumes a stop at every system (Boss 04:19Z); open question on a stay per stop. |
| NAV.42 | Each route stop shows the course and distance to the next stop | UX.35 | Boss 04:19Z. format_course per hop, frame per pair. |
| NAV.47 | Unknown-space jumps stop at scattered stars, black holes, neutron stars and quasars |  |  |
| NAV.48 | Offer to generate the uncharted sectors that block a course |  |  |

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
| OPS.35 | A vendored-version lock file for the static libraries |  | Research: supply-chain record. |

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
| UX.81 | Time symbols Gyr, Myr, kyr in place of Gy, My, ky; AU from 1,000,000 km; scientific text below mantissa 1e-3 |  | Research: wording and ladder fixes. |
| UX.78 | Unit preference: Automatic, Metric only or Customary | UX.23 | Research: under UX.23. |
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
| API.22 | An API version number: one sequential integer, shown in admin and in the status response |  | Boss 2026-10-09 20:59Z: done by the end of phase 1; do it early so later API changes bump it. |
| API.23 | The object ID as the public reference: pages, URLs, the API, wiki links and objectref use it in place of row ids | API.22, GEN.171, GEN.172 | Object-ID research. Breaking: bumps the API version (API.22). |
| API.7 | Investigate and plan upload limits |  | Plan only. |
| API.9 | Key scopes |  | Control schema migration (v8). Decided: user keys belong to accounts, so API.6 waits for USR.2 (phase 3+); API.9's scopes don't. |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.57 | A sector's contents depend only on the seed, the version and its address |  | Name collisions go away with GEN.67; GEN.63 dropped. Also keys GEN.64's ID collision counter on address (it checks stored rows, PR #406). |
| TEST.77 | A golden-seed regression test | GEN.57, GEN.135, OPS.28 |  |
| GEN.156 | Pin astropy to CODATA 2018 and IAU 2015 constants before its first import (GEN.66 follow-up) |  | Research follow-up to GEN.66 (built). |
| GEN.136 | Fingerprint encoding: floats to 9 significant digits, a stored leaf digest per sector and a ring-and-layer tree | GEN.135 | Research follow-up to GEN.58 (built). |
| GEN.135 | Deterministic math helpers for every stored float, with a lint and frozen constants |  | Research: Python 3.9 and 3.10 to 3.13 draw different sectors without it. Blocks TEST.77. |
| OPS.14 | A warning when the running version key differs from the galaxy's |  | Feeds GEN.58's output. |
| OPS.29 | Update reload: cover the gunicorn units and non-Apache hosts, and say what a reload aborts |  | Research follow-up to OPS.8 (built). |
| OPS.28 | Generator epoch and battery digest: say whether two checkouts generate the same galaxy | GEN.135 | Research: OPS.13 is built; this is the epoch it records. OPS.14, OPS.15, OPS.12 and TEST.77 read it. |
| OPS.37 | A Generator version number: one sequential integer, shown in admin and the API | OPS.28 | Boss 2026-10-09 20:59Z: done by the end of phase 1. Same number as OPS.28's epoch (decided). |

### Bugs from the GitHub issues

| ID | Item | Needs | Note |
|---|---|---|---|

### Foundations for the issue features

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.31 | Investigate: where generation spends its time, from the plan to a finished galaxy |  | GitHub issues [#761](https://github.com/dwhagar/planetGen/issues/761) and [#750](https://github.com/dwhagar/planetGen/issues/750) (benchmark half). |
| PERF.47 | The PERF.31 benchmark records the buffer pool, table sizes and worker start-up cost | PERF.31 | Generation performance study. |
| PERF.46 | Planets and moons: set the position once per body |  | Generation performance study. |
| DB.19 | Compact or derive phenomenon rows (needed only if the mass cut is lowered to 10 solar masses or less) |  | Decided: Boss accepted the 20 solar mass cut. |
| GEN.169 | Decide the phenomenon scatter rates: regional factors and the 0.1% intermediate-mass black holes |  | Needs Boss to decide. |
| PERF.41 | Stamp the tile and page caches with a TILE_FORMAT constant instead of the version (optional) | PERF.31 | Research follow-up to PERF.25 (built). |
| PERF.40 | Two shared queues, a reserved interactive worker and a real "cancel now" | PERF.39 | Research follow-up to PERF.24 (built). |
| PERF.39 | Every API job costs 2.4 s and 195 MB: import lazily and cap the burst workers |  | Research follow-up to PERF.24 (built); before batch wiki uploads. |
| PERF.38 | Cache fixes for the Galaxy Map under a fill: single-flight tile builds, a busy rule for the page cache, a deletion epoch in place of COUNT(*) |  | Research follow-up to PERF.34 (built). |
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
