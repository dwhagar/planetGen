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

Opens with object references (NAV.7) and the database consistency check (DB.8, then DB.9 repair). Then the work that needs phase 0 in place: prevalence controls, the other new planet classes, reproducible galaxies up to the golden-seed test (update key history, creation settings JSON, admin changes as a net diff), routing with no hop limit and the nearby search, the first picker pieces, the unit ladder, the API call log and the queue.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### References

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.7 | One reference for every object, with its parents |  | Moved from phase 0: groundwork for phase 1 features (picker, saved courses, nearby search); no phase 0 bug needs it. Root of the picker, saved courses, the 3D system view and account bookmarks. |
| NAV.8 | Pages and anchors for stars, planets, moons and belts | NAV.7 | System page anchors (system.html). |
| NAV.9 | Search and locate return references for every kind | NAV.7 | queryDb search and galaxy_locate. |

### Database consistency check

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.8 | Check a galaxy database and say whether it is damaged |  | Moved from phase 0: a new tool, not a bug fix; first item of phase 1 with DB.9 after it. Boss 02:13Z: phase 0, its own thread. Read-only. Stats and version checks switch on once GEN.44/PERF.11 and DB.6/DB.7 land; no hard dependency. |
| DB.9 | Repair a damaged galaxy database from a parity file | DB.8, GEN.57, GEN.44, GEN.58, OPS.14 | Boss 02:13Z: phase 1. Reed-Solomon parity over groups of sector exports; seed regeneration as fallback. |

### Classes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.33 | One class per PR, each with its tests |  | One class per PR (R and S first). Built on fixed physics so new classes aren't tuned to wrong masses, moons or zones. |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | GEN.33 | Class S landed with GEN.38 (PR #415); the other six classes here. Closes with GEN.33's PRs; PLANET_CLASSES in program_constants.py. |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | GEN.33 | Same reconcile/zone code as phase 0's physics fixes. |

### Prevalence

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.52 | Prevalence controls for sector and galaxy runs |  | Probability adjustments reach systemData.py, the same constructor phase 0 thread D fixed. |
| TEST.75 | Tests for forcing and prevalence | GEN.52 | Grows with GEN.49 to GEN.52. |
| ADM.16 | Prevalence controls on the Generate page | GEN.52, ADM.14 | generate.html, after ADM.14's layout. |
| GEN.48 | Forcing options are impractical for whole sectors; replace them with prevalence controls (bug) | GEN.52, ADM.16, TEST.75 | Parent; closes with its subitems. |

### Galaxy gen

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.24 | Generate the galactic core on layer 0 | ADM.14 | Bulk core fill runs on the parallel path; new mode on generate.html. |
| GEN.41 | Investigate: how much backfill work a density pre-pass would save | GEN.44 | Investigation; go/no-go for GEN.42. |

### Galaxy Map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.95 | A "Forward to current" button next to the map's Back and Forward |  | Moved from phase 0: a new button, not a bug; same files as MAP.94, next in the lane. Boss 04:29Z. Jumps to maxIndex of the map history (MAP.26). |

### System Map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.89 | System Map: space orbits with a fitted scale and a minimum ring gap instead of plain log |  | Same file as MAP.88 (done, PR #405); also changes systemmap.js kmToPx (MAP.63 touched it too). |

### Picker

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.13 | A picker module: select, step out, step in, step sideways | NAV.7 | picker.js; needs no engine. |
| NAV.14 | One breadcrumb for every level | NAV.13 |  |

### Routing

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.10 | Routing that scales past a few thousand systems |  | Galaxy schema migration for position indexes; queue behind PERF.11. Built with NAV.12 (no hop limit, a route always exists). |
| NAV.12 | No maximum hop length: a route always reaches the nearest star it can, across any number of sectors | NAV.34, TEST.79, NAV.10 | Anchor (Boss 01:53Z game mechanic). Built with NAV.10. Per-hop unknown-space flag in /api/nav. |
| UX.35 | The NAV page route shown horizontally, wrapping onto several lines on narrow screens |  | In parallel with NAV.12; replaces ol.nav-route in nav.html. |
| NAV.11 | Travel times for the system-to-system route too | NAV.10, NAV.12 | Times per hop, including unknown-space jumps; total assumes a stop at every system (Boss 04:19Z); open question on a stay per stop. |
| NAV.42 | Each route stop shows the course and distance to the next stop | UX.35 | Boss 04:19Z. format_course per hop, frame per pair. |

### Pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.23 | A shared unit-ladder module |  | Moved from phase 0: groundwork for UX.22 (a feature); UX.36 keeps today's formatter and needs no ladder. Shared unit ladder; UX.22 then UX.30 build on it. |
| UX.22 | Meaningful units for every measurement | UX.23 | One quantity family per PR. |
| UX.3 | Warn every visitor while a background job changes the galaxy |  | ETA from progress.json, which PERF.23 caps. |

### Queue

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.19 | Everything the API or web site starts runs on the work queue (investigate) |  | Audit only; nothing moves to the queue until it runs at any worker count. |
| ADM.15 | Change the worker count from the Queue page, with a "Ludicrous Speed" mode |  | Changing the worker count live only makes sense once any count works. |

### API

| ID | Item | Needs | Note |
|---|---|---|---|
| API.15 | Log every API call with its user, how it came in, and its HTTP response code |  | Moved from phase 0: a feature; no bug needs it. Needs no user accounts (Boss 01:31Z). |
| API.4 | API compatibility data in the docs |  | Docs and version number. |
| API.7 | Investigate and plan upload limits |  | Plan only. |
| API.9 | Key scopes |  | Control schema migration (v8). Decided: user keys belong to accounts, so API.6 waits for USR.2 (phase 3+); API.9's scopes don't. |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.56 | Every random draw in generation comes from the derived seeds |  | Touches every generator module. |
| GEN.57 | A sector's contents depend only on the seed, the version and its address | GEN.56, GEN.44 | Also keys GEN.64's ID collision counter on address (it checks stored rows, PR #406). |
| DB.7 | The version that generated each sector, and a warning for mixed-version galaxies |  |  |
| GEN.58 | A fingerprint of a galaxy's generated content |  | Judgment: phase 1 so the golden test guards later changes. |
| TEST.77 | A golden-seed regression test | GEN.57, GEN.58 |  |
| GEN.59 | Admin changes stored as a net difference from the generated galaxy | GEN.56, GEN.58, ADM.18 | Boss 02:20Z: net diff by stable address path, regenerate seeds, in ADM.18's JSON. Moved from phase 3. |
| OPS.8 | Update reloads Apache itself when run as root | OPS.7 | Moved from phase 0: not a bug; rides with OPS.13 in the same update scripts. Not a bug, but the same files as OPS.7, so it rides along. |
| OPS.13 | Every update records the version key, keeping the last 10 | OPS.7, OPS.8 | update.sh / update.ps1 after OPS.7 and OPS.8; control-database history table. Open question on "recalculate the seed value". |
| OPS.14 | A warning when the running version key differs from the galaxy's | DB.7, OPS.13 | Feeds GEN.58's output and OPS.12. |
| ADM.18 | The galaxy's creation settings saved as a JSON file, downloadable from the Admin dashboard | DB.7, OPS.13 | Boss 02:13Z: phase 1. Includes the key history; dated backup on every change. |
| GEN.63 | Planet names are unique within a sector | GEN.57 | Boss 04:03Z. Address-keyed clash rule from GEN.57; TEST.85 done (PR #403); rogue planets named by GEN.64's ID (PR #406). |

### Nearby search

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.43 | Find everything within a distance of a place: the query and the API | NAV.7 | Boss 04:39Z. Replaces systems_within_radius (one sector, systems only); enumerate_sectors_within_radius then per-sector reads; open questions: max distance, ungenerated sectors, generated only. |
| NAV.44 | A "What's nearby" page: pick a place, enter a distance in parsecs, list what is there | NAV.43 | Boss 04:39Z. Page with place picker, distance in pc, kind filters, 50-row pages. |

## Open questions for Boss

- NAV.8: Pages and anchors for stars, planets, moons and belts, see its entry in TODO.md.
- NAV.11: Travel times for the system-to-system route too, see its entry in TODO.md.
- UX.22: Meaningful units for every measurement, see its entry in TODO.md.
- UX.3: Warn every visitor while a background job changes the galaxy, see its entry in TODO.md.
- API.4: API compatibility data in the docs, see its entry in TODO.md.
- OPS.13: Every update records the version key, keeping the last 10, see its entry in TODO.md.
- NAV.43: Find everything within a distance of a place: the query and the API, see its entry in TODO.md.
