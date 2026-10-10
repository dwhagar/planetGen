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

Object references and the picker finished, the database check and repair, resumable runs, the reworked Generate page (spans, radial fills, random neighborhoods, directives, show on the map), galaxy generation changes (phenomena placed galaxy-wide first, fill order, backfill from the run's edge, nebula volume backfill, star-type tiers), the habitability index, tech levels and facility types, spin and orbital-update thresholds with their limits, routing with unknown-space stops and charting along a course, the nearby search, Select mode and view filters on the maps, the fly-through Galaxy Map (MAP.146: scroll-zoom to the cursor, double-click flight, stars that fade in with distance and zoom, a see-through near field; the headline of the 8.1 release), bookmarks, admin control from every screen, in-universe wording, the API log and scopes, and reproducible galaxies up to the golden-seed test.

The headline of this phase, and of the 8.1 release it ends with, is the fly-through Galaxy Map (MAP.146, with MAP.147 to MAP.155): scroll-zoom to the cursor, double-click flight, stars that fade in smoothly with distance and zoom, and a see-through near field (Boss, 2026-10-09 23:07Z). The 8.1 version stamp happens only when the whole phase, these items included, is done.

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
| DB.16 | Store the generator epoch and run id on each sector instead of four version text columns | OPS.28 | Research follow-up to DB.7 (built). |
| DOC.5 | Rewrite the object ID docs: object-ids.md, database-schema.md and api.md |  | Object-ID research. |
| TEST.110 | Object ID tests: identical IDs on 1 and 4 workers, none reused, none missing | GEN.172, GEN.176 | Object-ID research. |
| GEN.176 | A nebula or remnant is born in the sector holding the centre of the space it occupies |  | Object-ID research. |
| GEN.172 | Run-time births get object IDs from the counters |  | Object-ID research. |
| DB.9 | Repair a damaged galaxy database from a parity file | GEN.57 | Boss 02:13Z: phase 1. Reed-Solomon parity over groups of sector exports; seed regeneration as fallback. |
| DB.21 | A deep pass for the database check: validate every star system, with the estimated time shown first | PERF.33 | Boss 23:32Z; Foundations lane 1, small, near the end of its list. |
| DB.17 | Repair by regenerating a damaged sector from its seed when parity cannot rebuild it | DB.9, GEN.57 | Research split of DB.9: the regenerate-from-seed fallback. |

### Generation

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.29 | Record which runs a partly filled sector still needs |  | Builds on the generated flag that GEN.76 settles (sector_stats.bright_level_sol). |
| PERF.30 | Finish an interrupted block or sector run on the next start | PERF.29 |  |

### Generate page

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.181 | Directives: refuse impossible requests up front from the compound-Poisson tables |  | Left over from GEN.96. |
| GEN.180 | Directives: a forced fill (met_forced) after K failed draws |  | Left over from GEN.96. |
| GEN.179 | Store each sector's generation directive and attempt record with the sector |  | Left over from GEN.96. |
| GEN.186 | Random neighborhoods: an option to keep away from filled space |  | Left over from GEN.97 (merged, PR #950). |
| TEST.122 | Browser map tests fail on plain main in Bugfixes lane 2's container (fixture maps, controls, system-page maps) (bug) |  |  |
| ADM.49 | Galaxy shape density settings: the user changes the density range of the spiral arms, the inter-arm space, the core and the bulge |  | Boss 06:16Z via coordinator; unassigned. |

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
| GEN.178 | Magnetic fields: the induced field of an ocean moon |  | Left over from GEN.86 (PR #908); waits on the hydrosphere model. |
| GEN.177 | Planetary magnetic fields: a stagnant-lid factor |  | Left over from GEN.86 (PR #908). |

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
| NAV.11 | Travel times for the system-to-system route too |  | Times per hop, including unknown-space jumps; total assumes a stop at every system (Boss 04:19Z); open question on a stay per stop. |
| NAV.42 | Each route stop shows the course and distance to the next stop |  | Boss 04:19Z. format_course per hop, frame per pair. |
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
| UX.87 | The system list shows uncharted systems: every scattered star, with its location and a way to generate it |  | Boss issue #929, 01:51Z. |
| MAP.162 | A sector holding scattered objects but never generated can still be opened, marked uncharted |  | Boss issue #928, 01:44Z. |
| MAP.146 | Fly through the galaxy: scroll-zoom, double-click flight, distance-based visibility and a see-through near field | MAP.148, MAP.149, MAP.150, MAP.151, MAP.152, MAP.154, MAP.155, MAP.157, MAP.158, MAP.159 | Boss 2026-10-09 22:24Z and 22:56Z. Umbrella for MAP.148 to MAP.152; replaces the fixed ladder in galaxy/drill.py and galaxyprisms.js. |
| MAP.147 | The Galaxy Map wire format: the investigation (done) and the record of what was built from it | MAP.157, MAP.158, MAP.159 | Research Lane 3 report in docs/design/galaxy-map-wire-format.md; record until the build items are done. |
| MAP.157 | Trim the Galaxy Map tile JSON and serve it from prebuilt, precompressed bytes |  | MAP.147 recommendation step 1; no client change. |
| MAP.158 | A gentler tile prefetch and an IndexedDB tile cache instead of localStorage |  | MAP.147 recommendation step 1; client only. |
| MAP.161 | Load the Galaxy Map faster on a first visit: bundle or preload its scripts |  | Boss 2026-10-09 23:29Z: moved to Phase 1. Wire format report, finding 7; client side, no prerequisite. |
| MAP.148 | The star visibility law: apparent-magnitude opacity, flux-based brightness and an on-screen limit from a histogram |  | Fly-through report item 1; builds after MAP.153 (its first stage). |
| MAP.149 | The near field: depth fade, a see-through focus tube, drawing from inside a container, and picking that matches what is drawn |  | Fly-through report item 2; can start now. Folds MAP.121 blocker fade and MAP.141 context. |
| MAP.150 | The free camera: wheel zoom to the cursor, double-click flight, and the observer inside, with the container named from position | MAP.149 | Fly-through report item 3; needs the near field. |
| MAP.151 | The region data layer: exact-centred frame, aligned cells, per-level aggregates and slot-wrap ranges |  | Fly-through report item 4. Decide cache keys with MAP.147. |
| MAP.152 | Scale hand-offs: galaxy, sector and system cross-fade with hysteresis, and per-tile camera-relative origins | MAP.148, MAP.150, MAP.154 | Fly-through report item 5. |
| MAP.154 | Nested bright-star lists on the server, so every parent list is a subset of its child's |  | Zoom visibility note stage 2; shares tile keys with MAP.147 and MAP.151. |
| MAP.159 | Packed binary Galaxy Map tiles (quantised planes) with the nested tile lists, on one cache stamp bump | MAP.154, MAP.158 | MAP.147 recommendation step 2; one stamp bump with MAP.154 and MAP.151. |
| MAP.155 | Other objects fade in too: point objects from level 8, a size ramp for cloud sprites, and stars that grow from a faint dot |  | Zoom visibility note stage 3. |
| UX.91 | Planet and moon description carries a full PHI-4 explanation, each colour factor and why |  |  |
| PERF.57 | Skip empty stretches in a galactic scatter by combining layers into growing groups |  |  |
| GEN.196 | One "Redo scatters" box on Generate: choose which scatters to redo, with new settings for each | GEN.195 |  |
| MAP.166 | Galaxy Map "Dimmest star shown" says "every star" only when the view is complete |  |  |
| GEN.195 | A separate mass limit for neutron stars and black holes, and the central black hole or quasar always created |  |  |
| UX.90 | Explain the habitability chips in the web interface: a visible legend with the equipment labels, colours, thresholds and scores |  |  |
| MAP.165 | Scattered phenomena store a mass so the Galaxy Map sizes them exactly |  |  |
| PERF.56 | Record how long each stage of a staged job takes, with the settings it ran with |  |  |
| UX.89 | Staged jobs show wrong stage counts and numbers, and skipped stages are not listed (bug) |  |  |

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
| API.23 | The object ID as the public reference: pages, URLs, the API, wiki links and objectref use it in place of row ids | API.22, GEN.172 | Object-ID research. Breaking: bumps the API version (API.22). |
| API.7 | Investigate and plan upload limits |  | Plan only. |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.57 | A sector's contents depend only on the seed, the version and its address |  | Name collisions go away with GEN.67; GEN.63 dropped. Also keys GEN.64's ID collision counter on address (it checks stored rows, PR #406). |
| TEST.77 | A golden-seed regression test | GEN.57, GEN.135, OPS.28 |  |
| GEN.156 | Pin astropy to CODATA 2018 and IAU 2015 constants before its first import (GEN.66 follow-up) |  | Research follow-up to GEN.66 (built). |
| GEN.136 | Fingerprint encoding: floats to 9 significant digits, a stored leaf digest per sector and a ring-and-layer tree | GEN.135 | Research follow-up to GEN.58 (built). |
| GEN.135 | Deterministic math helpers for every stored float, with a lint and frozen constants |  | Research: Python 3.9 and 3.10 to 3.13 draw different sectors without it. Blocks TEST.77. |
| OPS.29 | Update reload: cover the gunicorn units and non-Apache hosts, and say what a reload aborts |  | Research follow-up to OPS.8 (built). |
| OPS.28 | Generator epoch and battery digest: say whether two checkouts generate the same galaxy | GEN.135 | Research: OPS.13 is built; this is the epoch it records. OPS.14, OPS.15, OPS.12 and TEST.77 read it. |
| OPS.37 | A Generator version number: one sequential integer, shown in admin and the API | OPS.28 | Boss 2026-10-09 20:59Z: done by the end of phase 1. Same number as OPS.28's epoch (decided). |

### Bugs from the GitHub issues

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.111 | test_ensure_sector_generated_creates_then_reuses_the_same_sector fails in a busy parallel run (bug) |  | Test flake reported 2026-10-09; Bugfixes lane 1. |
| TEST.123 | test_a_loaded_sector_knows_every_objects_cell_and_velocity fails once under full-suite load (bug) |  |  |
| TEST.116 | test_bughunt_end_to_end stores M2V or M6V where it expects K2V on the MySQL 8.4 and MariaDB legs once (bug) |  | Bugfixes lane 1 report 03:57Z. |

### Foundations for the issue features

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.31 | Investigate: where generation spends its time, from the plan to a finished galaxy |  | GitHub issues [#761](https://github.com/dwhagar/planetGen/issues/761) and [#750](https://github.com/dwhagar/planetGen/issues/750) (benchmark half). |
| PERF.47 | The PERF.31 benchmark records the buffer pool, table sizes and worker start-up cost | PERF.31 | Generation performance study. |
| PERF.46 | Planets and moons: set the position once per body |  | Generation performance study. |
| GEN.169 | Decide the phenomenon scatter rates: regional factors and the 0.1% intermediate-mass black holes |  | Needs Boss to decide. |
| GEN.187 | Bright-star back scatter by mass: rings of 1, 2, 5 and 8 solar masses around filled space |  | Boss issue #952, 03:43Z. |
| PERF.41 | Stamp the tile and page caches with a TILE_FORMAT constant instead of the version (optional) | PERF.31 | Research follow-up to PERF.25 (built). |
| PERF.40 | Two shared queues, a reserved interactive worker and a real "cancel now" | PERF.39 | Research follow-up to PERF.24 (built). |
| PERF.39 | Every API job costs 2.4 s and 195 MB: import lazily and cap the burst workers |  | Research follow-up to PERF.24 (built); before batch wiki uploads. |
| PERF.38 | Cache fixes for the Galaxy Map under a fill: single-flight tile builds, a busy rule for the page cache, a deletion epoch in place of COUNT(*) |  | Research follow-up to PERF.34 (built). |
| PERF.33 | Progress bars and ETAs from measured performance |  | GitHub issue [#661](https://github.com/dwhagar/planetGen/issues/661) (the use half). Bugfixes lane 1 (Boss, 23:38Z). |
| PERF.55 | One global progress bar for generation jobs that run in phases, with an ETA across all phases | PERF.33 | Boss 06:44Z via coordinator; unassigned. |
| GEN.128 | Design: multi-star hierarchies and compact-object primaries |  | GitHub issues [#777](https://github.com/dwhagar/planetGen/issues/777) and [#778](https://github.com/dwhagar/planetGen/issues/778): the research both builds wait on. |

## Open questions for Boss

- NAV.8: Pages and anchors for stars, planets, moons and belts, see its entry in TODO.md.
- NAV.11: Travel times for the system-to-system route too, see its entry in TODO.md.
- UX.22: Meaningful units for every measurement, see its entry in TODO.md.
- UX.3: Warn every visitor while a background job changes the galaxy, see its entry in TODO.md.
- API.4: API compatibility data in the docs, see its entry in TODO.md.
- OPS.13: Every update records the version key, keeping the last 10, see its entry in TODO.md.
